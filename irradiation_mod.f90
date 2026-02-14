!==============================================================
! irradiation_mod.f90  (Replacement, Scheme C with variable delay)
!
! Responsibilities:
!   - Measure instantaneous inner accretion rate mdot_in(t) from Sigma/nu
!   - Maintain a time history buffer of mdot_in(t) in PHYSICAL time [s]
!   - Compute delayed mdot_in(t - delay(t)) by linear interpolation (Scheme C)
!   - Set L_irr and LOH24 Eq.(16) parameters (rin_cgs, A1, L1, beta1, beta2, Q12)
!
! Notes:
!   - delay(t) can change each timestep; the buffer interpolation remains valid.
!   - Sigma evolution is in dimensionless time, but buffer uses t_phys = t_nd * t0.
!==============================================================
module irradiation_mod
  use kind_params, only : dp, i4b
  use constants,   only : pi, gg, cc
  use mod_global,  only : r, nu, sigmat, use_irradiation, use_irradiation_delay, use_be_decretion, &
                          use_finite_irradiation_source, M_star, kappa_es, L_star, R_star, R0, &
                          q, eta_acc, f_edd_cap, nu0_dim, nu0_nd, t0, tau_irr_lag_mode, tau_irr_lag_nd
  use units_disk_mod, only : mdot_edd_unit, r_dim
  use disk_flux_mod, only : measure_mdot_inner_from_arrays
  implicit none

  real(dp), parameter :: cos_inc_min = 1.0e-3_dp
  real(dp), parameter :: albedo      = 0.9_dp
  real(dp), parameter :: Cgap        = 2.0_dp
  real(dp), parameter :: Lirr_floor  = 1.0e-30_dp

  ! Public irradiation state used elsewhere
  real(dp) :: L_irr   = 0.0_dp
  real(dp) :: rin_cgs = 0.0_dp
  real(dp) :: A1      = 0.0_dp
  real(dp) :: L1      = 0.0_dp
  real(dp) :: beta1   = 0.0_dp
  real(dp) :: beta2   = 0.0_dp
  real(dp) :: Q12     = 0.0_dp

  ! Diagnostics
  real(dp) :: delay_sec_last  = 0.0_dp
  real(dp) :: Ledd_last       = 0.0_dp
  real(dp) :: Lacc_inst_last  = 0.0_dp
  real(dp) :: Lacc_lag_last   = 0.0_dp

  ! Buffer for delayed irradiation (use_irradiation_delay from mod_global)
  integer(i4b) :: ibuf = 0
  real(dp), allocatable :: t_hist(:)     ! [s]
  real(dp), allocatable :: mdot_hist(:)  ! [g/s]
  real(dp), allocatable :: t_arrive_hist(:) ! [s]
  real(dp), allocatable :: Lacc_hist(:)     ! [erg/s]
  integer(i4b) :: nbuf
  logical  :: buffer_ready = .false.
  integer(i4b) :: ibuf_next = 1
  integer(i4b) :: head_arrive = 1
  integer(i4b) :: nvalid_arrive = 0

  ! ---- hard cap to avoid overflow in (t_phys + delay_sec) etc. ----
  real(dp), parameter :: delay_cap = 1.0e60_dp   ! [s] huge but safe in dp arithmetic

contains

!---------------------------------------------------------------
! Initialize buffers (call once at start, and also when you reset)
!---------------------------------------------------------------
subroutine init_Larrive_buffer()
  use kind_params, only : dp, i4b
  implicit none
  integer(i4b) :: i
  do i = 1, nbuf
     t_arrive_hist(i) = -1.0_dp
     Lacc_hist(i)     =  0.0_dp
  end do
  ibuf_next    = 1
  !buffer_ready = .true.
end subroutine init_Larrive_buffer

!---------------------------------------------------------------
! Push one (arrival time, luminosity) sample into ring buffer
!---------------------------------------------------------------
subroutine push_Larrive(t_arrive, Lacc)
  real(dp), intent(in) :: t_arrive, Lacc

  ! write at current head
  t_arrive_hist(head_arrive) = t_arrive
  Lacc_hist(head_arrive)     = Lacc

  ! update valid counter
  if (nvalid_arrive < nbuf) nvalid_arrive = nvalid_arrive + 1
  !buffer_ready = (nvalid_arrive >= 2)

  ! advance ring pointer
  head_arrive = head_arrive + 1
  if (head_arrive > nbuf) head_arrive = 1
end subroutine push_Larrive

!---------------------------------------------------------------
! Sample luminosity that ARRIVES at time t_target.
!
! Policy:
! - If t_target is earlier than the earliest arrival time, return 0
!   (nothing has arrived yet).
! - If t_target is later than the latest arrival time, return the
!   latest value (or 0; choose what you want). I set it to latest
!   value to avoid sudden dropouts, but you can switch to 0.
! - Linear interpolation between bracketing arrival times.
!
! This is intentionally similar to your existing function, but it
! fixes the O(n^2) bracketing and uses sorting of valid indices.
!---------------------------------------------------------------
real(dp) function sample_Lacc_at(t_target)
  ! Sample luminosity that ARRIVES at time t_target (arrival-time model).
  ! - t_target < earliest arrival -> 0 (nothing arrived yet)
  ! - t_target > latest arrival   -> latest value
  ! - else: linear interpolation between bracketing arrival times
  use kind_params, only : dp, i4b
  implicit none
  real(dp), intent(in) :: t_target

  integer(i4b) :: i, nvalid, k
  integer(i4b), allocatable, save :: idx_work(:)
  integer(i4b), save :: idx_work_size = 0
  real(dp) :: tmin, tmax
  real(dp), parameter :: eps_t = 1.0e-12_dp

  if (.not. use_irradiation_delay) then
     sample_Lacc_at = Lacc_inst_last
     return
  end if

  sample_Lacc_at = 0.0_dp
  if (.not. buffer_ready) return

  if (.not. allocated(idx_work) .or. idx_work_size < nbuf) then
    if (allocated(idx_work)) deallocate(idx_work)
    allocate(idx_work(nbuf))
    idx_work_size = nbuf
  end if

  nvalid = 0
  do i = 1, nbuf
    if (t_arrive_hist(i) >= 0.0_dp) then
      nvalid = nvalid + 1
      idx_work(nvalid) = i
    end if
  end do

  if (nvalid == 0) return
  if (nvalid == 1) then
    if (t_target >= t_arrive_hist(idx_work(1)) - eps_t) sample_Lacc_at = Lacc_hist(idx_work(1))
    return
  end if

  call sort_idx_by_key_real_safe(t_arrive_hist, idx_work, nvalid)
  tmin = t_arrive_hist(idx_work(1))
  tmax = t_arrive_hist(idx_work(nvalid))

  if (t_target < tmin - eps_t) return
  if (t_target >= tmax - eps_t) then
    sample_Lacc_at = Lacc_hist(idx_work(nvalid))
    return
  end if

  k = bracket_rightmost_leq(t_arrive_hist, idx_work, nvalid, t_target)
  k = max(1, min(k, nvalid - 1))
  sample_Lacc_at = lerp_by_idx(t_arrive_hist, Lacc_hist, idx_work(k), idx_work(k+1), t_target)
end function sample_Lacc_at


!============================================================
! Helper: linear interpolation using two original indices i0,i1
!============================================================
real(dp) function lerp_by_idx(tarr, yarr, i0, i1, t)
  use kind_params, only : dp, i4b
  implicit none
  real(dp), intent(in) :: tarr(:), yarr(:)
  integer(i4b), intent(in) :: i0, i1
  real(dp), intent(in) :: t
  real(dp) :: t0, t1, y0, y1, denom
  real(dp), parameter :: tiny = 1.0e-99_dp

  t0 = tarr(i0); t1 = tarr(i1)
  y0 = yarr(i0); y1 = yarr(i1)
  denom = max(t1 - t0, tiny)
  lerp_by_idx = y0 + (y1 - y0) * (t - t0) / denom
end function lerp_by_idx


!============================================================
! Helper: binary search on sorted idx(:) by key=tarr(idx)
! returns k in [1, n-1] such that t(k) <= x < t(k+1)
!============================================================
integer(i4b) function bracket_rightmost_leq(tarr, idx, n, x) result(k)
  use kind_params, only : dp, i4b
  implicit none
  real(dp), intent(in) :: tarr(:)
  integer(i4b), intent(in) :: idx(:)
  integer(i4b), intent(in) :: n
  real(dp), intent(in) :: x
  integer(i4b) :: lo, hi, mid

  lo = 1
  hi = n - 1   ! ensure mid+1 exists

  do while (lo < hi)
    mid = (lo + hi + 1) / 2
    if (tarr(idx(mid)) <= x) then
      lo = mid
    else
      hi = mid - 1
    end if
  end do
  k = lo
end function bracket_rightmost_leq


!============================================================
! Safe fast sort: idx(1:n) sorted by key(idx(*))
! - iterative quicksort with bounds guards
! - small segments -> safe insertion sort
!============================================================
subroutine sort_idx_by_key_real_safe(key, idx, n)
  use kind_params, only : dp, i4b
  implicit none
  real(dp), intent(in) :: key(:)
  integer(i4b), intent(inout) :: idx(:)
  integer(i4b), intent(in) :: n

  integer(i4b), parameter :: SMALL = 24
  integer(i4b) :: lo_stack(128), hi_stack(128)
  integer(i4b) :: top, lo, hi
  integer(i4b) :: i, j, p, tmp
  real(dp) :: pivot

  if (n <= 1) return

  top = 1
  lo_stack(top) = 1
  hi_stack(top) = n

  do while (top > 0)
    lo = lo_stack(top)
    hi = hi_stack(top)
    top = top - 1

    if (hi - lo + 1 <= SMALL) then
      call insertion_sort_idx_by_key_real_safe(key, idx, lo, hi)
      cycle
    end if

    p = (lo + hi) / 2
    pivot = key(idx(p))

    i = lo
    j = hi
    do
      ! Move i right, but never beyond hi
      do while (i <= hi)
        if (key(idx(i)) >= pivot) exit
        i = i + 1
      end do

      ! Move j left, but never below lo
      do while (j >= lo)
        if (key(idx(j)) <= pivot) exit
        j = j - 1
      end do

      if (i >= j) exit

      tmp = idx(i); idx(i) = idx(j); idx(j) = tmp
      i = i + 1
      j = j - 1
    end do

    ! Partitions: [lo..j], [j+1..hi]
    if (j - lo > hi - (j+1)) then
      if (lo < j) then
        top = top + 1; lo_stack(top) = lo;   hi_stack(top) = j
      end if
      if (j+1 < hi) then
        top = top + 1; lo_stack(top) = j+1; hi_stack(top) = hi
      end if
    else
      if (j+1 < hi) then
        top = top + 1; lo_stack(top) = j+1; hi_stack(top) = hi
      end if
      if (lo < j) then
        top = top + 1; lo_stack(top) = lo;   hi_stack(top) = j
      end if
    end if

    if (top > size(lo_stack)) then
      write(*,'(a,i0)') 'FATAL: quicksort stack overflow. Increase stack size. top=', top
      stop
    end if
  end do
end subroutine sort_idx_by_key_real_safe


!============================================================
! Safe insertion sort:
! IMPORTANT: do not rely on short-circuit .and.
!============================================================
subroutine insertion_sort_idx_by_key_real_safe(key, idx, lo, hi)
  use kind_params, only : dp, i4b
  implicit none
  real(dp), intent(in) :: key(:)
  integer(i4b), intent(inout) :: idx(:)
  integer(i4b), intent(in) :: lo, hi

  integer(i4b) :: i, j, v
  real(dp) :: kv

  do i = lo + 1, hi
    v  = idx(i)
    kv = key(v)
    j = i - 1

    do while (j >= lo)
      if (key(idx(j)) > kv) then
        idx(j+1) = idx(j)
        j = j - 1
      else
        exit
      end if
    end do

    idx(j+1) = v
  end do
end subroutine insertion_sort_idx_by_key_real_safe


subroutine init_irradiation_buffer(nbuf_in)
  integer(i4b), intent(in) :: nbuf_in

  nbuf = nbuf_in

  if (allocated(t_arrive_hist)) deallocate(t_arrive_hist)
  if (allocated(Lacc_hist)) deallocate(Lacc_hist)
  if (allocated(t_hist)) deallocate(t_hist)
  if (allocated(mdot_hist)) deallocate(mdot_hist)
  allocate(t_arrive_hist(nbuf), Lacc_hist(nbuf))
  allocate(t_hist(nbuf), mdot_hist(nbuf))

  ! --- L_arrive buffer (future irradiation arrival) ---
  t_arrive_hist(:) = -1.0_dp   ! invalid
  Lacc_hist(:)     = 0.0_dp
  head_arrive      = 1
  nvalid_arrive    = 0

  ! --- mdot history buffer (past sampling) ---
  t_hist(:)    = -1.0_dp
  mdot_hist(:) = 0.0_dp
  ibuf = 0

  buffer_ready = .true.
end subroutine init_irradiation_buffer

!---------------------------------------------------------------
! Save/load arrival-time buffer for checkpoint (use_irradiation_delay only)
!---------------------------------------------------------------
subroutine save_irradiation_buffer(iu)
  integer(i4b), intent(in) :: iu
  logical :: written
  written = use_irradiation_delay .and. allocated(t_arrive_hist)
  write(iu) written
  if (.not. written) return
  write(iu) nbuf, head_arrive, nvalid_arrive
  write(iu) t_arrive_hist(:)
  write(iu) Lacc_hist(:)
end subroutine save_irradiation_buffer

subroutine load_irradiation_buffer(iu)
  integer(i4b), intent(in) :: iu
  logical :: written
  integer(i4b) :: nbuf_in, head_in, nvalid_in
  real(dp), allocatable :: t_tmp(:), L_tmp(:)
  read(iu) written
  if (.not. written) return
  read(iu) nbuf_in, head_in, nvalid_in
  if (use_irradiation_delay .and. allocated(t_arrive_hist) .and. nbuf_in == nbuf) then
    read(iu) t_arrive_hist(:)
    read(iu) Lacc_hist(:)
    head_arrive   = head_in
    nvalid_arrive = nvalid_in
  else
    allocate(t_tmp(nbuf_in), L_tmp(nbuf_in))
    read(iu) t_tmp(:)
    read(iu) L_tmp(:)
    deallocate(t_tmp, L_tmp)
  end if
end subroutine load_irradiation_buffer

subroutine set_irradiation_luminosity_from_arrays(t_nd_now, sigma_nd, nu_nd, do_push)
  real(dp), intent(in) :: t_nd_now
  real(dp), intent(in) :: sigma_nd(:), nu_nd(:)
  logical,  intent(in), optional :: do_push

  real(dp) :: t_phys, mdot_inst
  real(dp) :: Ledd, Lacc_inst
  logical :: push

  push = .true.
  if (present(do_push)) push = do_push

  if (.not. use_irradiation) then
    call set_zero_irradiation()
    return
  end if

  if (use_be_decretion) then
    L_irr = max(0.0_dp, L_star)
    call set_loh24_params(L_irr)
    return
  end if

  t_phys = t_nd_now * t0

  call measure_mdot_inner_from_arrays(sigma_nd, nu_nd, mdot_inst)

  Ledd = 4.0_dp*pi*gg*M_star*cc / max(kappa_es, 1.0e-99_dp)
  Lacc_inst = max(0.0_dp, eta_acc) * max(0.0_dp, mdot_inst) * cc*cc

  ! Arrival-time model: push (t_arrive, Lacc) where t_arrive = t_phys + delay
  ! Then sample L_irr = luminosity arriving at current time t_phys
  if (use_irradiation_delay) then
    block
      real(dp) :: delay_sec, t_arrive, Lacc_capped
      delay_sec = compute_delay_sec(nu_inner_nd=nu_nd(1))
      t_arrive = t_phys + delay_sec
      Lacc_capped = Lacc_inst
      if (f_edd_cap > 0.0_dp) Lacc_capped = min(Lacc_inst, f_edd_cap * Ledd)
      if (push .and. buffer_ready) call push_Larrive(t_arrive, Lacc_capped)
      L_irr = sample_Lacc_at(t_phys)
      delay_sec_last = delay_sec
      Lacc_lag_last  = L_irr
    end block
  else
    if (f_edd_cap > 0.0_dp) then
      L_irr = min(Lacc_inst, f_edd_cap * Ledd)
    else
      L_irr = Lacc_inst
    end if
    delay_sec_last  = 0.0_dp
    Lacc_lag_last   = Lacc_inst
  end if

  Ledd_last       = Ledd
  Lacc_inst_last  = Lacc_inst

  call set_loh24_params(L_irr)
end subroutine set_irradiation_luminosity_from_arrays

  subroutine set_zero_irradiation()
    L_irr   = 0.0_dp
    rin_cgs = 0.0_dp
    A1      = 0.0_dp
    L1      = 0.0_dp
    beta1   = 0.0_dp
    beta2   = 0.0_dp
    Q12     = 0.0_dp
  end subroutine set_zero_irradiation

  subroutine set_loh24_params(Lirr_in)
    real(dp), intent(in) :: Lirr_in
    rin_cgs = R0
    A1      = 1.0_dp - albedo
    L1      = 1.0_dp / (1.0_dp+q) * max(0.0_dp, Lirr_in)
    beta1   = 1.5_dp * q**2 / (1.0_dp + q)**2 / Cgap**2
    beta2   = 1.5_dp / (1.0_dp + q)**2 / Cgap**2
    Q12     = q
  end subroutine set_loh24_params

  subroutine push_mdot_history(t_phys, mdot_phys)
    real(dp), intent(in) :: t_phys, mdot_phys
    ibuf = ibuf + 1
    if (ibuf > nbuf) ibuf = 1
    t_hist(ibuf)    = t_phys
    mdot_hist(ibuf) = mdot_phys
  end subroutine push_mdot_history

  real(dp) function compute_delay_sec(nu_inner_nd)
    ! Delay time [s].
    !   tau_irr_lag_mode='explicit': tau_irr_lag_nd * t0
    !   tau_irr_lag_mode='viscous':  tau_irr_lag_nd * (r_in^2 / nu_in) at inner edge
    real(dp), intent(in), optional :: nu_inner_nd

    real(dp) :: r_in_cgs, nu_dim_in
    real(dp), parameter :: tiny = 1.0e-99_dp

    compute_delay_sec = 0.0_dp
    if (.not. use_irradiation_delay) return

    if (trim(adjustl(tau_irr_lag_mode)) == 'viscous') then
       ! Viscous time at inner edge: t_visc = r^2 / nu
       if (.not. present(nu_inner_nd)) then
          compute_delay_sec = tau_irr_lag_nd * t0
          compute_delay_sec = min(compute_delay_sec, delay_cap)
          return
       end if
       r_in_cgs = r_dim(r(1))
       nu_dim_in = 0.0_dp
       if (nu0_dim > 0.0_dp .and. nu0_nd > 0.0_dp) then
          nu_dim_in = (nu_inner_nd / max(nu0_nd, tiny)) * nu0_dim
       end if
       if (r_in_cgs <= 0.0_dp .or. nu_dim_in <= 0.0_dp) then
          compute_delay_sec = delay_cap
       else
          compute_delay_sec = tau_irr_lag_nd * (r_in_cgs * r_in_cgs) / nu_dim_in
          compute_delay_sec = min(compute_delay_sec, delay_cap)
       end if
    else
       ! Explicit: fixed delay in dimensionless units
       compute_delay_sec = tau_irr_lag_nd * t0
       compute_delay_sec = min(compute_delay_sec, delay_cap)
    end if
  end function compute_delay_sec

real(dp) function sample_mdot_at(t_target)
  ! Sample mdot at arbitrary time using (t_hist, mdot_hist) buffer.
  ! - Works even if delay changes each step
  ! - Safe for t_target < 0 (clamps to earliest sample)
  ! - Fast: sort indices O(n log n), then binary search O(log n)
  !
  use kind_params, only : dp, i4b
  implicit none
  real(dp), intent(in) :: t_target

  integer(i4b) :: i, nvalid, k
  real(dp) :: tmin, tmax

  ! Work arrays (reuse)
  integer(i4b), allocatable, save :: idx_work(:)
  integer(i4b), save :: idx_work_size = 0

  sample_mdot_at = 0.0_dp
  if (.not. buffer_ready) return

  ! Ensure workspace allocated
  if (.not. allocated(idx_work) .or. idx_work_size < nbuf) then
    if (allocated(idx_work)) deallocate(idx_work)
    allocate(idx_work(nbuf))
    idx_work_size = nbuf
  end if

  ! Gather valid indices: you already mark invalid with t_hist < 0
  nvalid = 0
  do i = 1, nbuf
    if (t_hist(i) >= 0.0_dp) then
      nvalid = nvalid + 1
      idx_work(nvalid) = i
    end if
  end do
  if (nvalid < 2) return

  ! Sort by time
  call sort_idx_by_key_real_safe(t_hist, idx_work, nvalid)

  tmin = t_hist(idx_work(1))
  tmax = t_hist(idx_work(nvalid))

  ! Clamp outside range
  if (t_target <= tmin) then
    sample_mdot_at = mdot_hist(idx_work(1))
    return
  end if
  if (t_target >= tmax) then
    sample_mdot_at = mdot_hist(idx_work(nvalid))
    return
  end if

  ! Find bracket k such that t(k) <= t_target < t(k+1)
  k = bracket_rightmost_leq(t_hist, idx_work, nvalid, t_target)

  ! Linear interpolation
  sample_mdot_at = lerp_by_idx(t_hist, mdot_hist, idx_work(k), idx_work(k+1), t_target)

end function sample_mdot_at


  subroutine compute_shadow_from_HoverR(nr, r_cgs, H_cgs, shadow)
    use kind_params, only : dp, i4b
    implicit none
    integer(i4b), intent(in)  :: nr
    real(dp),     intent(in)  :: r_cgs(nr), H_cgs(nr)
    logical,      intent(out) :: shadow(nr)

    integer(i4b) :: i
    real(dp) :: hoverr, hoverr_max

    shadow(:) = .false.
    if (nr <= 0) return

    shadow(1) = .false.
    hoverr_max = H_cgs(1) / max(r_cgs(1), 1.0e-99_dp)

    do i = 2, nr
       hoverr = H_cgs(i) / max(r_cgs(i), 1.0e-99_dp)
       if (hoverr <= hoverr_max) then
          shadow(i) = .true.
       else
          shadow(i) = .false.
          hoverr_max = hoverr
       end if
    end do
  end subroutine compute_shadow_from_HoverR

subroutine compute_shadow_loggrid_hyst(nr, r_cgs, H_cgs, halfwin, eps_on, eps_off, shadow)
  ! Shadowing with a true ON/OFF hysteresis around an envelope (running max).
  !
  ! - Build h = H/r
  ! - Smooth h on a log-r grid (moving average in index space)
  ! - Maintain a running maximum h_max (envelope)
  ! - Shadow turns ON if h <= (1-eps_on)*h_max
  ! - Shadow turns OFF only if h >= (1-eps_off)*h_max   (eps_off < eps_on)
  !
  ! This suppresses flickering near the threshold while preserving the
  ! "one tall bump can shadow the outside" geometry, but makes it less jumpy.

  use kind_params, only : dp, i4b
  implicit none
  integer(i4b), intent(in)  :: nr, halfwin
  real(dp),     intent(in)  :: r_cgs(nr), H_cgs(nr)
  real(dp),     intent(in)  :: eps_on, eps_off
  logical,      intent(out) :: shadow(nr)

  integer(i4b) :: i, j, j1, j2, cnt
  real(dp) :: hoverr(nr), hoverr_s(nr), hoverr_max
  real(dp), parameter :: tiny = 1.0e-99_dp

  shadow(:) = .false.
  if (nr <= 0) return

  ! Build H/r
!$omp parallel do default(shared) private(i)
  do i = 1, nr
     hoverr(i) = H_cgs(i) / max(r_cgs(i), tiny)
  end do
!$omp end parallel do

  ! Smooth on log-r grid (index window)
!$omp parallel do default(shared) private(i, j, j1, j2, cnt)
  do i = 1, nr
     j1  = max(1, i-halfwin)
     j2  = min(nr, i+halfwin)
     cnt = j2 - j1 + 1
     hoverr_s(i) = 0.0_dp
     do j = j1, j2
        hoverr_s(i) = hoverr_s(i) + hoverr(j)
     end do
     hoverr_s(i) = hoverr_s(i) / max(real(cnt,dp), 1.0_dp)
  end do
!$omp end parallel do

  ! Envelope + true hysteresis
  shadow(1)  = .false.
  hoverr_max = hoverr_s(1)

  do i = 2, nr

     ! Always update the envelope with the current smoothed value
     hoverr_max = max(hoverr_max, hoverr_s(i))

     if (.not. shadow(i-1)) then
        ! Currently illuminated: decide whether to ENTER shadow
        if (hoverr_s(i) <= (1.0_dp - eps_on) * hoverr_max) then
           shadow(i) = .true.
        else
           shadow(i) = .false.
        end if
     else
        ! Currently shadowed: decide whether to EXIT shadow
        if (hoverr_s(i) >= (1.0_dp - eps_off) * hoverr_max) then
           shadow(i) = .false.
        else
           shadow(i) = .true.
        end if
     end if

  end do

end subroutine compute_shadow_loggrid_hyst


subroutine compute_shadow_loggrid_smooth(nr, r_cgs, H_cgs, halfwin, eps_shadow, shadow)
  use kind_params, only : dp, i4b
  implicit none
  integer(i4b), intent(in)  :: nr, halfwin
  real(dp),     intent(in)  :: r_cgs(nr), H_cgs(nr), eps_shadow
  logical,      intent(out) :: shadow(nr)

  integer(i4b) :: i, j, j1, j2, cnt
  real(dp) :: hoverr(nr), hoverr_s(nr), hoverr_max

  shadow(:) = .false.
  if (nr <= 0) return

  !-----------------------------------------
  ! Build H/r (dimensionless); safe divide
  !-----------------------------------------
!$omp parallel do default(shared) private(i)
  do i = 1, nr
     hoverr(i) = H_cgs(i) / max(r_cgs(i), 1.0e-99_dp)
  end do
!$omp end parallel do

  !-----------------------------------------
  ! Smooth H/r on a log-r grid using a simple
  ! moving average in index space.
  ! (Log grid => equal spacing in ln r, so this is OK.)
  !-----------------------------------------
!$omp parallel do default(shared) private(i, j, j1, j2, cnt)
  do i = 1, nr
     j1  = max(1, i-halfwin)
     j2  = min(nr, i+halfwin)
     cnt = j2 - j1 + 1
     hoverr_s(i) = 0.0_dp
     do j = j1, j2
        hoverr_s(i) = hoverr_s(i) + hoverr(j)
     end do
     hoverr_s(i) = hoverr_s(i) / max(real(cnt,dp), 1.0_dp)
  end do
!$omp end parallel do

  !-----------------------------------------
  ! Shadow decision with hysteresis:
  ! illuminated only if H/r exceeds previous max
  ! by a margin eps_shadow.
  !
  ! shadow=.true. means blocked (no irradiation).
  !-----------------------------------------
  shadow(1) = .false.
  hoverr_max = hoverr_s(1)

  do i = 2, nr
     if (hoverr_s(i) <= (1.0_dp - eps_shadow) * hoverr_max) then
        shadow(i) = .true.
     else
        shadow(i) = .false.
        hoverr_max = max(hoverr_max, hoverr_s(i))
     end if
  end do

end subroutine compute_shadow_loggrid_smooth


  !-----------------------------------------------------------------------
  subroutine compute_Y_dYdXi(n, xi, r_cgs, H_cgs, Y, dYdXi)
    ! Compute Y = H/r and dY/dxi on a nonuniform xi grid using 2nd-order formulas.
    integer(i4b), intent(in) :: n
    real(dp), intent(in) :: xi(n), r_cgs(n), H_cgs(n)
    real(dp), intent(out) :: Y(n), dYdXi(n)

    integer(i4b) :: i
    real(dp) :: x0, x1, x2, f0, f1, f2

!$omp parallel do default(shared) private(i)
    do i = 1, n
      if (r_cgs(i) > 0.0_dp) then
        Y(i) = H_cgs(i) / r_cgs(i)
      else
        Y(i) = 0.0_dp
      end if
    end do
!$omp end parallel do

    if (n < 3) then
      dYdXi(:) = 0.0_dp
      return
    end if

    ! One-sided 2nd-order at i=1 (quadratic through points 1,2,3)
    x0 = xi(1); x1 = xi(2); x2 = xi(3)
    f0 = Y(1);  f1 = Y(2);  f2 = Y(3)
    dYdXi(1) = f0*(2.0_dp*x0-x1-x2)/((x0-x1)*(x0-x2)) + &
               f1*(x0-x2)/((x1-x0)*(x1-x2)) + &
               f2*(x0-x1)/((x2-x0)*(x2-x1))

    ! Central 2nd-order for interior points (quadratic through i-1,i,i+1)
!$omp parallel do default(shared) private(i, x0, x1, x2, f0, f1, f2)
    do i = 2, n-1
      x0 = xi(i-1); x1 = xi(i); x2 = xi(i+1)
      f0 = Y(i-1);  f1 = Y(i);  f2 = Y(i+1)
      dYdXi(i) = f0*(x1-x2)/((x0-x1)*(x0-x2)) + &
                 f1*(2.0_dp*x1-x0-x2)/((x1-x0)*(x1-x2)) + &
                 f2*(x1-x0)/((x2-x0)*(x2-x1))
    end do
!$omp end parallel do

    ! One-sided 2nd-order at i=n (quadratic through n,n-1,n-2)
    x0 = xi(n); x1 = xi(n-1); x2 = xi(n-2)
    f0 = Y(n);  f1 = Y(n-1);  f2 = Y(n-2)
    dYdXi(n) = f0*(2.0_dp*x0-x1-x2)/((x0-x1)*(x0-x2)) + &
               f1*(x0-x2)/((x1-x0)*(x1-x2)) + &
               f2*(x0-x1)/((x2-x0)*(x2-x1))
  end subroutine compute_Y_dYdXi

  !-----------------------------------------------------------------------
  subroutine compute_Y_dYdXi_sg(n, xi, r_cgs, H_cgs, halfwin, poly_order, Y, dYdXi)
    use kind_params, only: dp, i4b
    implicit none
    integer(i4b), intent(in) :: n, halfwin, poly_order
    real(dp), intent(in)  :: xi(n), r_cgs(n), H_cgs(n)
    real(dp), intent(out) :: Y(n), dYdXi(n)

    integer(i4b) :: i
    real(dp) :: yraw(n), ypre(n)

    ! Build raw Y = H/r
!$omp parallel do default(shared) private(i)
    do i = 1, n
       if (r_cgs(i) > 0.0_dp) then
          yraw(i) = H_cgs(i) / r_cgs(i)
       else
          yraw(i) = 0.0_dp
       end if
    end do
!$omp end parallel do

    ! One-pass 3-point median prefilter (robust against spikes)
    ypre = yraw
    if (n >= 3) then
       do i = 2, n-1
          ypre(i) = median3(yraw(i-1), yraw(i), yraw(i+1))
       end do
       ypre(1) = yraw(1)
       ypre(n) = yraw(n)
    end if

    ! Savitzky–Golay-like smoothing + derivative on a nonuniform grid
    call sg_smooth_and_deriv_nonuniform(n, xi, ypre, halfwin, poly_order, Y, dYdXi)

  contains

    pure real(dp) function median3(a, b, c)
      real(dp), intent(in) :: a, b, c
      if ((a <= b .and. b <= c) .or. (c <= b .and. b <= a)) then
         median3 = b
      else if ((b <= a .and. a <= c) .or. (c <= a .and. a <= b)) then
         median3 = a
      else
         median3 = c
      end if
    end function median3

  end subroutine compute_Y_dYdXi_sg

  !-----------------------------------------------------------------------
  subroutine sg_smooth_and_deriv_nonuniform(n, x, y, halfwin, p, y_sm, dy_dx)
    use kind_params, only: dp, i4b
    implicit none
    integer(i4b), intent(in) :: n, halfwin, p
    real(dp), intent(in) :: x(n), y(n)
    real(dp), intent(out) :: y_sm(n), dy_dx(n)

    integer(i4b) :: i, j, k, m, i1, i2, np, nw, w
    real(dp) :: dx, wgt, dxmax, lam, diag_scale
    real(dp), allocatable :: ATA(:,:), ATb(:), a(:)
    integer(i4b), allocatable :: ipiv(:)
    integer(i4b) :: info

    ! Basis size guard (p<=15 assumed in your original code)
    if (p > 15) stop "p too large for basis guard"

    np = p + 1

    allocate(ATA(np,np), ATb(np), a(np), ipiv(np))

    do i = 1, n

       ! Window indices (clamped at boundaries)
       w = 2*halfwin + 1
       if (w > n) then
          i1 = 1; i2 = n
       else if (i <= halfwin + 1) then
          i1 = 1; i2 = w
       else if (i >= n - halfwin) then
          i2 = n; i1 = n - w + 1
       else
          i1 = i - halfwin
          i2 = i + halfwin
       end if
       nw = i2 - i1 + 1

       ! If too few points, fall back
       if (nw < np) then
          y_sm(i) = y(i)
          if (i == 1) then
             dy_dx(i) = (y(2) - y(1)) / max(x(2) - x(1), 1.0e-99_dp)
          else if (i == n) then
             dy_dx(i) = (y(n) - y(n-1)) / max(x(n) - x(n-1), 1.0e-99_dp)
          else
             dy_dx(i) = (y(i+1) - y(i-1)) / max(x(i+1) - x(i-1), 1.0e-99_dp)
          end if
          cycle
       end if

       ATA(:,:) = 0.0_dp
       ATb(:)   = 0.0_dp

       ! dxmax for distance-based weights (simple triangular window)
       dxmax = max( abs(x(i1) - x(i)), abs(x(i2) - x(i)) )
       dxmax = max(dxmax, 1.0e-99_dp)

       do j = i1, i2
          dx = x(j) - x(i)

          ! Weight: triangular in |dx| (robust near boundaries)
          wgt = 1.0_dp - abs(dx)/dxmax
          if (wgt < 0.0_dp) wgt = 0.0_dp

          ! Build basis for polynomial in (x - x(i)) to improve conditioning
          a(1) = 1.0_dp
          do k = 2, np
             a(k) = a(k-1) * dx
          end do

          do k = 1, np
             ATb(k) = ATb(k) + wgt * a(k) * y(j)
             do m = 1, np
                ATA(k,m) = ATA(k,m) + wgt * a(k) * a(m)
             end do
          end do
       end do

       ! Light ridge regularization to avoid near-singular windows:
       ! scale lambda by a typical diagonal magnitude
       diag_scale = 0.0_dp
       do k = 1, np
          diag_scale = max(diag_scale, abs(ATA(k,k)))
       end do
       diag_scale = max(diag_scale, 1.0e-99_dp)
       lam = 1.0e-12_dp * diag_scale   ! start small; tune if needed

       do k = 1, np
          ATA(k,k) = ATA(k,k) + lam
       end do

       ! Solve ATA * coeff = ATb using LAPACK DGESV (LU with pivoting)
       ! Note: DGESV overwrites ATA and ATb; solution returned in ATb
       call dgesv(np, 1, ATA, np, ipiv, ATb, np, info)

       if (info /= 0) then
          ! Fallback: do not update; keep raw or finite diff
          y_sm(i)  = y(i)
          if (i == 1) then
             dy_dx(i) = (y(2) - y(1)) / max(x(2) - x(1), 1.0e-99_dp)
          else if (i == n) then
             dy_dx(i) = (y(n) - y(n-1)) / max(x(n) - x(n-1), 1.0e-99_dp)
          else
             dy_dx(i) = (y(i+1) - y(i-1)) / max(x(i+1) - x(i-1), 1.0e-99_dp)
          end if
       else
          ! Coefficients: a0 = ATb(1), a1 = ATb(2)
          y_sm(i)  = ATb(1)
          if (np >= 2) then
             dy_dx(i) = ATb(2)
          else
             dy_dx(i) = 0.0_dp
          end if
       end if

    end do

    deallocate(ATA, ATb, a, ipiv)
  end subroutine sg_smooth_and_deriv_nonuniform
  
  !-----------------------------------------------------------------------
  subroutine solve_small_linear(n, A_in, b_in, x, info)
    use kind_params, only: dp, i4b
    implicit none
    integer(i4b), intent(in)  :: n
    real(dp),     intent(in)  :: A_in(n,n)
    real(dp),     intent(in)  :: b_in(n)
    real(dp),     intent(out) :: x(n)
    integer(i4b), intent(out) :: info

    real(dp) :: A(n,n), b(n)
    integer(i4b) :: i, j, k, piv
    real(dp) :: maxv, tmp, factor, eps, colmax

    A = A_in
    b = b_in
    x = 0.0_dp
    info = 0

    eps = epsilon(1.0_dp)

    ! Forward elimination with partial pivoting
    do k = 1, n-1

       ! Compute a scale for the k-th column (for robust singularity check)
       colmax = 0.0_dp
       do i = k, n
          colmax = max(colmax, abs(A(i,k)))
       end do
       if (colmax <= 0.0_dp) then
          info = k   ! column is exactly zero -> singular
          return
       end if

       piv = k
       maxv = abs(A(k,k))
       do i = k+1, n
          if (abs(A(i,k)) > maxv) then
             maxv = abs(A(i,k))
             piv  = i
          end if
       end do

       ! Relative singularity check
       if (maxv <= eps * colmax) then
          info = k
          return
       end if

       ! Row swap if needed
       if (piv /= k) then
          do j = k, n
             tmp = A(k,j); A(k,j) = A(piv,j); A(piv,j) = tmp
          end do
          tmp = b(k); b(k) = b(piv); b(piv) = tmp
       end if

       ! Eliminate below pivot
       do i = k+1, n
          factor = A(i,k) / A(k,k)
          A(i,k) = 0.0_dp
          do j = k+1, n
             A(i,j) = A(i,j) - factor * A(k,j)
          end do
          b(i) = b(i) - factor * b(k)
       end do

    end do

    ! Final pivot check
    colmax = abs(A(n,n))
    if (colmax <= eps * max(1.0_dp, colmax)) then
       info = n
       return
    end if

    ! Back substitution
    x(n) = b(n) / A(n,n)
    do i = n-1, 1, -1
       tmp = b(i)
       do j = i+1, n
          tmp = tmp - A(i,j) * x(j)
       end do
       if (abs(A(i,i)) <= eps * max(1.0_dp, abs(A(i,i)))) then
          info = i
          return
       end if
       x(i) = tmp / A(i,i)
    end do

  end subroutine solve_small_linear


  subroutine compute_Qirr_eq16_raw(n, xi, Y, dYdXi, rin_cgs, A1, L1, Q12, beta1, beta2, Qirr_raw_out)
    use kind_params, only : dp, i4b
    implicit none
    integer(i4b), intent(in) :: n
    real(dp), intent(in)  :: xi(n), Y(n), dYdXi(n)
    real(dp), intent(in)  :: rin_cgs, A1, L1, Q12, beta1, beta2
    real(dp), intent(out) :: Qirr_raw_out(n)

    integer(i4b) :: i
    real(dp) :: pref, xi2, term1, term2, bracket
    real(dp) :: qraw

    if (rin_cgs <= 0.0_dp) then
      Qirr_raw_out(:) = 0.0_dp
      return
    end if

    pref = (A1 * L1) / (2.0_dp * pi * rin_cgs * rin_cgs)

!$omp parallel do default(shared) private(i, xi2, term1, term2, bracket)
    do i = 1, n
       ! --- Put the Eq.16 expression here, WITHOUT max(.,0) ---
       ! qraw = ... (your current compute_Qirr_eq16 formula) ...
       if (xi(i) <= 0.0_dp) then
          Qirr_raw_out(i) = 0.0_dp
       else
          xi2 = xi(i)*xi(i)
          term1 = (1.0_dp + Q12) * dYdXi(i)
          term2 = (beta1 + Q12*beta2) / xi2 * ( Y(i)/xi(i) - 0.5_dp*dYdXi(i) )
          bracket = term1 - term2
          qraw = pref * (1.0_dp/xi(i)) * bracket
       end if
       Qirr_raw_out(i) = qraw
       !write (*, '("i =", i3, ": Y =", 1pe12.4, ", dYdXi =", 1pe12.4, &
       !       " -> Qirr =", 1pe12.4)') i, Y(i), dYdXi(i), Qirr_raw_out(i)
    end do
!$omp end parallel do
  end subroutine compute_Qirr_eq16_raw


   subroutine compute_Qirr_eq16(n, xi, Y, dYdXi, rin_cgs, A1, L1, Q12, beta1, beta2, Qirr_out)
    use kind_params, only : dp, i4b
    implicit none
    integer(i4b), intent(in) :: n
    real(dp), intent(in)  :: xi(n), Y(n), dYdXi(n)
    real(dp), intent(in)  :: rin_cgs, A1, L1, Q12, beta1, beta2
    real(dp), intent(out) :: Qirr_out(n)

    real(dp) :: Qraw(n)

    call compute_Qirr_eq16_raw(n, xi, Y, dYdXi, rin_cgs, A1, L1, Q12, beta1, beta2, Qraw)
    Qirr_out(:) = max(Qraw(:), 0.0_dp)
  end subroutine compute_Qirr_eq16


  subroutine build_Qirr_profile_eq16_with_raw(n, r_cgs, H_cgs_in, Qirr_raw_out, Qirr_out, dYdXi_out, Y_out)
    use kind_params, only : dp, i4b
    use mod_global, only : use_finite_irradiation_source, R_star
    use units_disk_mod, only : r_dim
    implicit none
    integer(i4b), intent(in) :: n
    real(dp), intent(in)  :: r_cgs(n)
    real(dp), intent(in)  :: H_cgs_in(n)     ! CGS
    real(dp), intent(out) :: Qirr_raw_out(n) ! raw Eq.16 output (can be negative)
    real(dp), intent(out) :: Qirr_out(n)     ! after non-negativity enforcement
    real(dp), intent(out) :: dYdXi_out(n)
    real(dp), intent(out) :: Y_out(n)

    real(dp) :: xi(n), Y(n), dYdXi(n)
    real(dp) :: alpha_star, alpha_flare, factor, pref_star, Qirr_flat_star
    integer(i4b) :: i

    if (rin_cgs <= 0.0_dp .or. L1 <= 0.0_dp .or. A1 <= 0.0_dp) then
       Qirr_raw_out(:) = 0.0_dp
       Qirr_out(:)     = 0.0_dp
       dYdXi_out(:)    = 0.0_dp
       Y_out(:)        = 0.0_dp
       return
    end if

    do i = 1, n
       xi(i)    = r_cgs(i) / rin_cgs
    end do

    call compute_Y_dYdXi_sg(n=n, xi=xi, r_cgs=r_cgs, H_cgs=H_cgs_in, &
                          halfwin=7, poly_order=2, Y=Y, dYdXi=dYdXi)

    ! Raw Eq.16 (can be negative)
    call compute_Qirr_eq16_raw(n=n, xi=xi, Y=Y, dYdXi=dYdXi, rin_cgs=rin_cgs, &
                          A1=A1, L1=L1, Q12=Q12, beta1=beta1, beta2=beta2, &
                          Qirr_raw_out=Qirr_raw_out)

    ! Enforce non-negative irradiation for the solver
    Qirr_out(:) = max(Qirr_raw_out(:), 0.0_dp)

    ! Finite-size stellar disk: additive grazing angle (Chiang & Goldreich; LOH24)
    ! alpha_total = alpha_star + alpha_flare
    !   alpha_star = 0.4*R_star/r (finite stellar disk)
    !   alpha_flare = xi*dY/dxi = r*d(H/r)/dr (disk flaring, already in LOH24)
    ! alpha_flare >= alpha_star: Qirr = Qirr_LOH24 * (alpha_star+alpha_flare)/alpha_flare
    ! alpha_flare <  alpha_star: Qirr = (alpha_star+alpha_flare)/alpha_star * Qirr_flat_star
    !   (smooth, no blow-up; Qirr_flat_star = pref*0.2*(R_star/rin)/xi^3)
    if (use_finite_irradiation_source .and. R_star > 0.0_dp) then
       pref_star = (A1 * L1) / (2.0_dp * pi * rin_cgs * rin_cgs)
       do i = 1, n
          if (r_cgs(i) > 0.0_dp) then
             alpha_star = 0.4_dp * R_star / r_cgs(i)
             alpha_flare = max(xi(i) * dYdXi(i), 0.0_dp)
             Qirr_flat_star = pref_star * 0.2_dp * (R_star / rin_cgs) / (xi(i)**3)
             if (alpha_flare >= alpha_star) then
                factor = (alpha_star + alpha_flare) / alpha_flare
                Qirr_out(i) = Qirr_out(i) * factor
             else
                factor = (alpha_star + alpha_flare) / alpha_star
                Qirr_out(i) = factor * Qirr_flat_star
             end if
          end if
       end do
    end if

    dYdXi_out(:) = dYdXi(:)
    Y_out(:)     = Y(:)
  end subroutine build_Qirr_profile_eq16_with_raw


subroutine build_Qirr_profile_eq16(n, r_cgs, H_cgs_in, Qirr_out, dYdXi_out)
  use kind_params, only : dp, i4b
  use mod_global, only : use_finite_irradiation_source, R_star
  use units_disk_mod, only : r_dim
  implicit none
  integer(i4b), intent(in) :: n
  real(dp), intent(in)  :: r_cgs(n)
  real(dp), intent(in)  :: H_cgs_in(n)     ! CGS
  real(dp), intent(out) :: Qirr_out(n)
  real(dp), intent(out) :: dYdXi_out(n)

  real(dp) :: xi(n), Y(n), dYdXi(n)
  real(dp) :: alpha_star, alpha_flare, factor, pref_star, Qirr_flat_star
  integer(i4b) :: i

  if (rin_cgs <= 0.0_dp .or. L1 <= 0.0_dp .or. A1 <= 0.0_dp) then
     Qirr_out(:)  = 0.0_dp
     dYdXi_out(:) = 0.0_dp
     return
  end if

  do i = 1, n
     xi(i)    = r_cgs(i) / rin_cgs
  end do

  call compute_Y_dYdXi_sg(n=n, xi=xi, r_cgs=r_cgs, H_cgs=H_cgs_in, &
                        halfwin=7, poly_order=2, Y=Y, dYdXi=dYdXi)
  call compute_Qirr_eq16(n=n, xi=xi, Y=Y, dYdXi=dYdXi, rin_cgs=rin_cgs, &
                        A1=A1, L1=L1, Q12=Q12, beta1=beta1, beta2=beta2, &
                        Qirr_out=Qirr_out)

  ! Finite-size stellar disk: additive grazing angle (same as _with_raw)
  if (use_finite_irradiation_source .and. R_star > 0.0_dp) then
     pref_star = (A1 * L1) / (2.0_dp * pi * rin_cgs * rin_cgs)
     do i = 1, n
        if (r_cgs(i) > 0.0_dp) then
           alpha_star = 0.4_dp * R_star / r_cgs(i)
           alpha_flare = max(xi(i) * dYdXi(i), 0.0_dp)
           Qirr_flat_star = pref_star * 0.2_dp * (R_star / rin_cgs) / (xi(i)**3)
           if (alpha_flare >= alpha_star) then
              factor = (alpha_star + alpha_flare) / alpha_flare
              Qirr_out(i) = Qirr_out(i) * factor
           else
              factor = (alpha_star + alpha_flare) / alpha_star
              Qirr_out(i) = factor * Qirr_flat_star
           end if
        end if
     end do
  end if

  dYdXi_out(:) = dYdXi(:)
end subroutine build_Qirr_profile_eq16


subroutine fix_Qirr_negative_spikes(n, r_cgs, shadow, Qirr_raw, Qirr_prof, &
                                   max_run, Qirr_floor_abs, use_log_r)
  ! Post-process irradiation profile:
  ! - Copy raw -> prof
  ! - For short negative runs (length <= max_run) in non-shadowed zones,
  !   replace by interpolation between the nearest positive neighbors.
  ! - For other negatives, apply an absolute floor (or keep as-is and let caller zero).
  !
  ! Notes:
  ! - This routine does NOT enforce shadow=0. Do it outside (caller), to keep
  !   geometry and masking logic in one place.
  ! - The interpolation is done in r (or log r), but only between valid neighbors
  !   that are non-shadowed and have positive Qirr_raw.

  use kind_params, only: dp, i4b
  implicit none
  integer(i4b), intent(in) :: n
  real(dp),     intent(in) :: r_cgs(n)
  logical,      intent(in) :: shadow(n)
  real(dp),     intent(in) :: Qirr_raw(n)
  real(dp),     intent(out):: Qirr_prof(n)
  integer(i4b), intent(in) :: max_run
  real(dp),     intent(in) :: Qirr_floor_abs
  logical,      intent(in) :: use_log_r

  integer(i4b) :: i, i1, i2, L, jL, jR, j
  real(dp) :: xL, xR, x, qL, qR, t
  real(dp), parameter :: tiny = 1.0e-99_dp
  logical :: near_shadow

  ! Start from raw
  Qirr_prof(:) = Qirr_raw(:)

  i = 1
  do while (i <= n)

     ! Only consider negative points in NON-shadow zones
     if ( (.not. shadow(i)) .and. (Qirr_raw(i) < 0.0_dp) ) then

        ! Identify the contiguous negative run [i1, i2]
        i1 = i
        i2 = i
        do while (i2 < n)
           if ( (.not. shadow(i2+1)) .and. (Qirr_raw(i2+1) < 0.0_dp) ) then
              i2 = i2 + 1
           else
              exit
           end if
        end do
        L = i2 - i1 + 1

        ! Find left neighbor jL: must be in the SAME non-shadow contiguous segment.
        jL = i1 - 1
        do while (jL >= 1)
           if (shadow(jL)) then
              jL = 0              ! hit a shadow barrier: do NOT look further left
              exit
           end if
           if (Qirr_raw(jL) > 0.0_dp) exit
           jL = jL - 1
        end do

        ! Find right neighbor jR: must be in the SAME non-shadow contiguous segment.
        jR = i2 + 1
        do while (jR <= n)
           if (shadow(jR)) then
              jR = n + 1          ! hit a shadow barrier: do NOT look further right
              exit
           end if
           if (Qirr_raw(jR) > 0.0_dp) exit
           jR = jR + 1
        end do

        if ( (L <= max_run) .and. (jL >= 1) .and. (jR <= n) ) then
           ! We have a short negative run [i1,i2] and two valid positive anchors (jL,jR).

           qL = Qirr_raw(jL)
           qR = Qirr_raw(jR)

           ! --- Do NOT "heal" negatives if this run touches a shadow boundary ---
           ! Rationale: shadow transitions are intrinsically discontinuous; filling there
           ! tends to create artificial ramps/spikes.
           !
           ! Here, "touches a shadow boundary" means either anchor is shadowed OR the run
           ! is adjacent to a shadowed cell (i1-1 or i2+1).
           near_shadow = .false.

           if (shadow(jL) .or. shadow(jR)) near_shadow = .true.
           if (i1 > 1) then
              if (shadow(i1-1)) near_shadow = .true.
           end if
           if (i2 < n) then
              if (shadow(i2+1)) near_shadow = .true.
           end if

           if (.not. near_shadow) then
              ! Conservative fill: avoid creating steep artificial ramps.
              ! Using the smaller positive neighbor is robust when one side is a sharp dip.
              do j = i1, i2
                 Qirr_prof(j) = min(qL, qR)
              end do
              ! Alternative (smoother but still bounded): geometric mean
              ! do j = i1, i2
              !    Qirr_prof(j) = sqrt(max(qL,tiny)*max(qR,tiny))
              ! end do
           else
              ! Near shadow boundary: do not fill; just apply a floor (typically 0).
              do j = i1, i2
                 Qirr_prof(j) = max(Qirr_prof(j), Qirr_floor_abs)
              end do
           end if

        else
           ! One-sided or no-sided: do NOT "repair" (avoid artificial steps near shadow boundary).
           do j = i1, i2
              Qirr_prof(j) = Qirr_raw(j)     ! keep raw
           end do
        end if

        i = i2 + 1
     else
        i = i + 1
     end if

  end do

  ! Optional: apply floor everywhere in non-shadow zones to avoid tiny/negative remnants
  !do i = 1, n
  !   if (.not. shadow(i)) then
  !      Qirr_prof(i) = max(Qirr_prof(i), Qirr_floor_abs)
  !   end if
  !end do

end subroutine fix_Qirr_negative_spikes

subroutine compute_illum_factor_loggrid_hyst(nr, r_cgs, H_cgs, halfwin, eps_on, eps_off, &
                                            w_relax, illum)
  ! Continuous illumination factor (0..1) using the same envelope idea.
  !
  ! f = (H/r)_smooth / envelope
  ! If f <= 1-eps_on  -> illum = 0   (shadow)
  ! If f >= 1-eps_off -> illum = 1   (lit)
  ! Between them      -> linear ramp (or keep previous if you prefer)
  !
  ! w_relax: under-relaxation for stability in coupled T<->H iterations
  !          (0<w_relax<=1). Use ~0.1-0.3 to suppress oscillation.

  use kind_params, only : dp, i4b
  implicit none
  integer(i4b), intent(in)  :: nr, halfwin
  real(dp),     intent(in)  :: r_cgs(nr), H_cgs(nr)
  real(dp),     intent(in)  :: eps_on, eps_off, w_relax
  real(dp),     intent(out) :: illum(nr)

  integer(i4b) :: i, j, j1, j2, cnt
  real(dp) :: hoverr(nr), hoverr_s(nr)
  real(dp) :: env, env_old, f, fon, foff, illum_new
  real(dp), parameter :: tiny = 1.0e-99_dp

  if (nr <= 0) return

  ! Build H/r
  do i = 1, nr
     hoverr(i) = H_cgs(i) / max(r_cgs(i), tiny)
  end do

  ! Smooth on log-r grid (index window)
  do i = 1, nr
     j1  = max(1, i-halfwin)
     j2  = min(nr, i+halfwin)
     cnt = j2 - j1 + 1
     hoverr_s(i) = 0.0_dp
     do j = j1, j2
        hoverr_s(i) = hoverr_s(i) + hoverr(j)
     end do
     hoverr_s(i) = hoverr_s(i) / max(real(cnt,dp), 1.0_dp)
  end do

  fon  = 1.0_dp - eps_on
  foff = 1.0_dp - eps_off   ! foff > fon

  ! Initialize
  illum(:) = 1.0_dp
  env = max(hoverr_s(1), tiny)
  illum(1) = 1.0_dp

  do i = 2, nr
     env_old = env
     env     = max(env, hoverr_s(i))  ! running max envelope

     f = hoverr_s(i) / max(env, tiny)

     ! Base continuous mapping
     if (f <= fon) then
        illum_new = 0.0_dp
     else if (f >= foff) then
        illum_new = 1.0_dp
     else
        ! linear ramp across the hysteresis band
        illum_new = (f - fon) / max(foff - fon, tiny)
     end if

     ! Under-relaxation to suppress oscillation in coupled iterations
     illum(i) = (1.0_dp - w_relax) * illum(i-1) + w_relax * illum_new
  end do

end subroutine compute_illum_factor_loggrid_hyst

end module irradiation_mod

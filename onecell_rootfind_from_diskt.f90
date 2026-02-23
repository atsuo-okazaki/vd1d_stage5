program onecell_rootfind_from_diskt
  ! Standalone one-cell thermal balance root-finding using disk_t*.dat output.
  ! Reads:
  !   # M_star = ... [g]  from header
  ! and data columns:
  !   xi r sigma H Y dY/dxi rho nu Qvis Qirr Qrad Tmid Tsurf kappaR tauR
  !
  ! Solves f(T) = Qplus_visc(T) + Qplus_irr(T) - Qminus(T) = 0
  ! by bracketing + bisection, calling heating_cooling_cell() each time.

  use kind_params, only: dp, i4b
  use mod_global, only : nt
  use setup_mod, only : setup
  use opacity_table_mod, only : init_opacity_tables
  use disk_thermal_mod, only: heating_cooling_cell
  implicit none

  real(dp), parameter :: G_cgs = 6.67430e-8_dp

  character(len=256) :: fname, arg
  real(dp) :: r_target
  real(dp) :: M_star_g
  integer(i4b) :: unit, ios

  ! Selected row variables
  real(dp) :: xi, r_cgs, Sigma_cgs, H_loc, Y, dYdxi
  real(dp) :: rho_loc, nu_dim_loc
  real(dp) :: Qvis_file, Qirr_in, Qrad_file
  real(dp) :: Tmid_init, Tsurf, kappaR_file, tauR_file

  real(dp) :: OmegaK_cgs
  logical :: shadow

  ! Root-find work variables
  real(dp) :: T0, Tlo, Thi, flo, fhi, Tmid, fmid
  real(dp) :: kappa_loc, tau_loc, Qplus_visc, Qplus_irr, Qminus
  real(dp) :: Qirr_work
  integer(i4b) :: it
  logical :: bracket_ok

  call get_command_argument(1, fname)
  call get_command_argument(2, arg)
  if (len_trim(fname) == 0 .or. len_trim(arg) == 0) then
    write(*,*) 'Usage: ./onecell_rootfind <disk_t00020001.dat> <r_target_cgs>'
    stop 2
  endif
  read(arg,*) r_target

  ! --- minimal init (NO time evolution) ---
  call setup()
  call init_opacity_tables()

  ! Read file: parse M_star from header, select closest row to r_target
  call read_mstar_and_pick_row(trim(fname), r_target, M_star_g, &
                              xi, r_cgs, Sigma_cgs, H_loc, Y, dYdxi, rho_loc, nu_dim_loc, &
                              Qvis_file, Qirr_in, Qrad_file, Tmid_init, Tsurf, kappaR_file, tauR_file)

  OmegaK_cgs = sqrt(G_cgs * M_star_g / (r_cgs**3))

  ! Approximate shadow flag (adjust if you have a better definition)
  shadow = (abs(Qirr_in) <= 0.0_dp)

  write(*,'(a)') '=== Picked row (closest to r_target) ==='
  write(*,'(a,1p,e14.6)') 'M_star[g]  = ', M_star_g
  write(*,'(a,1p,e14.6)') 'r_target   = ', r_target
  write(*,'(a,1p,e14.6)') 'r_cgs      = ', r_cgs
  write(*,'(a,1p,e14.6)') 'Sigma_cgs  = ', Sigma_cgs
  write(*,'(a,1p,e14.6)') 'H_loc      = ', H_loc
  write(*,'(a,1p,e14.6)') 'rho_loc    = ', rho_loc
  write(*,'(a,1p,e14.6)') 'nu_dim_loc = ', nu_dim_loc
  write(*,'(a,1p,e14.6)') 'OmegaK     = ', OmegaK_cgs
  write(*,'(a,l1)')       'shadow     = ', shadow
  write(*,'(a,1p,e14.6)') 'Tmid_init  = ', Tmid_init
  write(*,'(a,1p,e14.6)') 'Qvis(file) = ', Qvis_file
  write(*,'(a,1p,e14.6)') 'Qirr(file) = ', Qirr_in
  write(*,'(a,1p,e14.6)') 'Qrad(file) = ', Qrad_file
  write(*,'(a,1p,e14.6)') 'kappaR(file)=', kappaR_file
  write(*,'(a,1p,e14.6)') 'tauR(file)  =', tauR_file

  ! Initialize inout variables for heating_cooling_cell
  kappa_loc = kappaR_file   ! reasonable initial guess
  tau_loc   = tauR_file
  Qirr_work = Qirr_in

  T0 = max(1.0_dp, Tmid_init)

  call eval_f(T0, fmid, kappa_loc, tau_loc, Qplus_visc, Qplus_irr, Qminus, Qirr_work)
  write(*,'(a)') '--- init (evaluated by heating_cooling_cell) ---'
  call print_state(T0, fmid, kappa_loc, tau_loc, Qplus_visc, Qplus_irr, Qminus, Qirr_work)

  ! Bracket search around T0
  Tlo = max(1.0_dp, T0/2.0_dp)
  Thi = min(1.0e8_dp, T0*2.0_dp)
  bracket_ok = .false.

  do it = 1, 50
    call eval_f(Tlo, flo, kappa_loc, tau_loc, Qplus_visc, Qplus_irr, Qminus, Qirr_work)
    call eval_f(Thi, fhi, kappa_loc, tau_loc, Qplus_visc, Qplus_irr, Qminus, Qirr_work)
    if (flo*fhi < 0.0_dp) then
      bracket_ok = .true.
      exit
    endif
    Tlo = max(1.0_dp, Tlo/2.0_dp)
    Thi = min(1.0e8_dp, Thi*2.0_dp)
  enddo

  if (.not. bracket_ok) then
    write(*,'(a)') 'ERROR: failed to bracket root for f(T)=0.'
    write(*,'(a,1p,3e14.6)') 'T0,Tlo,Thi=', T0, Tlo, Thi
    write(*,'(a,1p,2e14.6)') 'flo,fhi   =', flo, fhi
    stop 1
  endif

  ! Bisection
  do it = 1, 250
    Tmid = 0.5_dp*(Tlo+Thi)
    call eval_f(Tmid, fmid, kappa_loc, tau_loc, Qplus_visc, Qplus_irr, Qminus, Qirr_work)

    if (abs(fmid) <= 1.0e-8_dp*max(1.0_dp, abs(Qminus))) exit
    if (abs(Thi-Tlo) <= 1.0e-10_dp*max(1.0_dp, abs(Tmid))) exit

    call eval_f(Tlo, flo, kappa_loc, tau_loc, Qplus_visc, Qplus_irr, Qminus, Qirr_work)
    if (flo*fmid < 0.0_dp) then
      Thi = Tmid
    else
      Tlo = Tmid
    endif
  enddo

  write(*,'(a)') '--- root (evaluated by heating_cooling_cell) ---'
  call print_state(Tmid, fmid, kappa_loc, tau_loc, Qplus_visc, Qplus_irr, Qminus, Qirr_work)

contains

  subroutine eval_f(Tmid, f, kappa_loc, tau_loc, Qplus_visc, Qplus_irr, Qminus, Qirr_work)
    real(dp), intent(in) :: Tmid
    real(dp), intent(out) :: f
    real(dp), intent(inout) :: kappa_loc, tau_loc, Qirr_work
    real(dp), intent(out) :: Qplus_visc, Qplus_irr, Qminus
    real(dp) :: kappa_planck_loc

    call heating_cooling_cell(r_cgs, Sigma_cgs, OmegaK_cgs, shadow, Tmid, &
                              H_loc, rho_loc, nu_dim_loc, kappa_loc, kappa_planck_loc, &
                              tau_loc, Qplus_visc, Qplus_irr, Qminus,     &
                              Qirr_in=Qirr_work)
    f = (Qplus_visc + Qplus_irr) - Qminus
  end subroutine eval_f

  subroutine print_state(Tmid, f, kappa_loc, tau_loc, Qplus_visc, Qplus_irr, Qminus, Qirr_work)
    real(dp), intent(in) :: Tmid, f, kappa_loc, tau_loc, Qplus_visc, Qplus_irr, Qminus, Qirr_work
    write(*,'(a,1p,e14.6)') 'Tmid      = ', Tmid
    write(*,'(a,1p,e14.6)') 'kappa_loc = ', kappa_loc
    write(*,'(a,1p,e14.6)') 'tau_loc   = ', tau_loc
    write(*,'(a,1p,e14.6)') 'Qvis      = ', Qplus_visc
    write(*,'(a,1p,e14.6)') 'Qirr      = ', Qplus_irr
    write(*,'(a,1p,e14.6)') 'Qrad      = ', Qminus
    write(*,'(a,1p,e14.6)') 'Qirr_in   = ', Qirr_work
    write(*,'(a,1p,e14.6)') 'f=Q+-Q-   = ', f
    write(*,'(a,1p,e14.6)') '|R|       = ', abs(f)/max(1.0_dp, abs(Qminus))
  end subroutine print_state

  subroutine read_mstar_and_pick_row(fname, r_target, M_star_g, &
                                    xi, r_cgs, Sigma_cgs, H_loc, Y, dYdxi, rho_loc, nu_dim_loc, &
                                    Qvis, Qirr, Qrad, Tmid, Tsurf, kappaR, tauR)
    character(len=*), intent(in) :: fname
    real(dp), intent(in) :: r_target
    real(dp), intent(out) :: M_star_g
    real(dp), intent(out) :: xi, r_cgs, Sigma_cgs, H_loc, Y, dYdxi, rho_loc, nu_dim_loc
    real(dp), intent(out) :: Qvis, Qirr, Qrad, Tmid, Tsurf, kappaR, tauR

    integer(i4b) :: unit, ios
    character(len=1024) :: line
    real(dp) :: xi_t, r_t, sig_t, H_t, Y_t, dY_t, rho_t, nu_t
    real(dp) :: Qv_t, Qi_t, Qr_t, Tc_t, Ts_t, kap_t, tau_t
    real(dp) :: best_dr
    logical :: have_mstar

    M_star_g = -1.0_dp
    have_mstar = .false.
    best_dr = huge(1.0_dp)

    open(newunit=unit, file=trim(fname), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: cannot open ', trim(fname)
      stop 2
    endif

    do
      read(unit,'(A)',iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle

      if (line(1:1) == '#') then
        if (index(line, 'M_star') > 0) then
          call extract_mstar_g(line, M_star_g, have_mstar)
        endif
        cycle
      endif

      ! Data row: 15 columns as documented in the header
      read(line,*,iostat=ios) xi_t, r_t, sig_t, H_t, Y_t, dY_t, rho_t, nu_t, &
                               Qv_t, Qi_t, Qr_t, Tc_t, Ts_t, kap_t, tau_t
      if (ios /= 0) cycle

      if (abs(r_t - r_target) < best_dr) then
        best_dr = abs(r_t - r_target)
        xi = xi_t; r_cgs = r_t; Sigma_cgs = sig_t
        H_loc = H_t; Y = Y_t; dYdxi = dY_t
        rho_loc = rho_t; nu_dim_loc = nu_t
        Qvis = Qv_t; Qirr = Qi_t; Qrad = Qr_t
        Tmid = Tc_t; Tsurf = Ts_t
        kappaR = kap_t; tauR = tau_t
      endif
    enddo

    close(unit)

    if (.not. have_mstar) then
      write(*,*) 'ERROR: could not find M_star in header.'
      stop 2
    endif
    if (best_dr == huge(1.0_dp)) then
      write(*,*) 'ERROR: no data rows read.'
      stop 2
    endif
  end subroutine read_mstar_and_pick_row

  subroutine extract_mstar_g(line, M_star_g, ok)
    character(len=*), intent(in) :: line
    real(dp), intent(out) :: M_star_g
    logical, intent(out) :: ok
    integer(i4b) :: p_eq, p_br, p_com, ierr
    character(len=256) :: sub

    ok = .false.
    M_star_g = -1.0_dp

    p_eq  = index(line, '=')
    if (p_eq <= 0) return

    ! Take substring after '='
    sub = adjustl(line(p_eq+1:))

    ! Cut at '[' or ',' if present
    p_br = index(sub, '[')
    p_com = index(sub, ',')
    if (p_br > 0) sub = sub(:p_br-1)
    if (p_com > 0) sub = sub(:p_com-1)

    read(sub,*,iostat=ierr) M_star_g
    if (ierr == 0 .and. M_star_g > 0.0_dp) ok = .true.
  end subroutine extract_mstar_g

end program onecell_rootfind_from_diskt

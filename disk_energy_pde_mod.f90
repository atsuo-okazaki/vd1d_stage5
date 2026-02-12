!==============================================================
! disk_energy_pde_mod.f90  (FINAL replacement module)
!
! Purpose:
!   - Evolve midplane temperature Tmid with local (no radial transport) energy ODE:
!       Sigma * cV * dT/dt = Qvis(T,Sigma) + Qirr(r,t) - Qrad(T,Sigma)
!
!   - Rebuild disk vertical structure and viscosity consistently from an
!     existing Tmid slice (for initialization/restart):
!       rebuild_structure_from_current_T(it, Sigma_new, nu_new)
!
! Key design choices (MVP-stable):
!   - Q terms are in CGS; thus dt_phys = dt_nd * t0 is required.
!   - Qirr(r,t) is computed from LOH24 Eq.(16) using irradiation_mod parameters.
!   - Geometry used for Qirr is "lagged" by one time slice (H from previous step)
!     in the temperature advance routine; in rebuild it uses current H slice.
!
! Dependencies:
!   - disk_thermal_mod::heating_cooling_cell must accept Qirr_in=...
!   - irradiation_mod must provide: rin_cgs, A1, L1, Q12, beta1, beta2,
!       compute_Y_dYdXi, compute_Qirr_eq16
!==============================================================
module disk_energy_pde_mod
  use kind_params, only : dp, i4b
  use constants,   only : kb, mp, mu
  use mod_global,  only : nr, r, t0, nu0_dim, nu0_nd, nu_conv, &
                          Tmid, H, rho, kappaR, tauR, Qvis, Qrad, Qirr, dYdXi, &
                          use_irradiation, is_shadow
  use units_disk_mod, only : r_dim, omegaK_dim, sigma_dim
  use radiation_params_mod, only : T_floor, T_ceiling
  use disk_thermal_mod, only : heating_cooling_cell
  use irradiation_mod, only : rin_cgs, A1, L1, Q12, beta1, beta2, &
                              compute_Y_dYdXi, compute_Qirr_eq16
  implicit none

contains

  !------------------------------------------------------------
  subroutine advance_temperature_and_structure(itp1, dt_nd, sigma_new_nd, nu_new_nd)
    ! Advance T from (itp1-1) to itp1, then recompute structure and viscosity.
    integer(i4b), intent(in)    :: itp1
    real(dp),     intent(in)    :: dt_nd
    real(dp),     intent(in)    :: sigma_new_nd(nr)
    real(dp),     intent(inout) :: nu_new_nd(nr)

    real(dp) :: dt_phys
    real(dp) :: r_cgs(nr), Sigma_cgs(nr), OmegaK_cgs(nr)
    real(dp) :: H_geom(nr), Y(nr), dY(nr), Qirr_prof(nr)
    real(dp) :: T_old, T_new
    real(dp) :: H_loc, rho_loc, kappa_loc, tau_loc, nu_dim_loc
    real(dp) :: Qv, Qi_tmp, Qc
    integer(i4b) :: i
    logical :: shadow_cell(nr)

    real(dp), parameter :: tiny_sigma = 1.0e-8_dp

    dt_phys = dt_nd * t0

    call build_shadow_flags_safe(nr, is_shadow)

    ! Units for this step
!$omp parallel do default(shared) private(i)
    do i = 1, nr
      r_cgs(i)     = r_dim(r(i))
      Sigma_cgs(i) = sigma_dim(sigma_new_nd(i))
      if (r_cgs(i) > 0.0_dp) then
        OmegaK_cgs(i) = omegaK_dim(r(i))
      else
        OmegaK_cgs(i) = 0.0_dp
      end if
    end do
!$omp end parallel do

    ! Geometry for irradiation: use previous H slice (lagged geometry)
    if (use_irradiation .and. itp1 > 1) then
      H_geom(:) = H(itp1-1, :)
    else
      H_geom(:) = 0.0_dp
    end if

    ! LOH24 Eq.(16) irradiation profile
    call build_Qirr_profile(nr, r, r_cgs, H_geom, &
                            Y, dY, Qirr_prof, dYdXi(itp1, :))

    ! Temperature update per cell (implicit Euler + Newton)
    do i = 1, nr

      if (Sigma_cgs(i) <= tiny_sigma .or. OmegaK_cgs(i) <= 0.0_dp) then
        Tmid(itp1,i)   = 0.0_dp
        H(itp1,i)      = 0.0_dp
        rho(itp1,i)    = 0.0_dp
        kappaR(itp1,i) = 0.0_dp
        tauR(itp1,i)   = 0.0_dp
        Qvis(itp1,i)   = 0.0_dp
        Qrad(itp1,i)   = 0.0_dp
        Qirr(itp1,i)   = 0.0_dp
        if (itp1 > 1) then
          nu_new_nd(i) = nu_conv(itp1-1,i)
        end if
        cycle
      end if

      T_old = 0.0_dp
      if (itp1 > 1) T_old = Tmid(itp1-1, i)
      if (T_old <= 0.0_dp) T_old = max(T_floor, 0.5_dp*(T_floor + T_ceiling))

      call solve_T_implicit_cell(r_cgs(i), Sigma_cgs(i), OmegaK_cgs(i), shadow_cell(i), &
                                 dt_phys, T_old, Qirr_prof(i), T_new)
      Tmid(itp1,i) = T_new

      ! Recompute structure + Q terms consistently at T_new using provided Qirr_in
      call heating_cooling_cell(r_cgs(i), Sigma_cgs(i), OmegaK_cgs(i), shadow_cell(i), T_new, &
                                H_loc, rho_loc, nu_dim_loc, kappa_loc, tau_loc, Qv, Qi_tmp, Qc, &
                                Qirr_in=Qirr_prof(i))
      H(itp1,i)      = H_loc
      rho(itp1,i)    = rho_loc
      kappaR(itp1,i) = kappa_loc
      tauR(itp1,i)   = tau_loc
      Qvis(itp1,i)   = Qv
      Qrad(itp1,i)   = Qc
      Qirr(itp1,i)   = Qirr_prof(i)

      ! Map nu_dim [cm^2/s] to nu_nd
      if (nu0_dim > 0.0_dp) then
        nu_new_nd(i) = (nu_dim_loc / nu0_dim) * nu0_nd
      else
        if (itp1 > 1) nu_new_nd(i) = nu_conv(itp1-1,i)
      end if
      nu_conv(itp1,i) = nu_new_nd(i)

    end do

    ! if you keep history shadow array in mod_global:
    is_shadow(itp1, :) = shadow_cell(:)

  end subroutine advance_temperature_and_structure

  !------------------------------------------------------------
  subroutine rebuild_structure_from_current_T(it, sigma_new_nd, nu_new_nd)
    ! Rebuild vertical structure + Q terms + nu for a given time slice "it"
    ! using the EXISTING Tmid(it,:) values (no temperature evolution here).
    !
    ! Intended uses:
    !   - After init_initial_conditions() sets Sigma and Tmid at it=1
    !   - After reading a checkpoint if you want to refresh derived arrays
    !
    integer(i4b), intent(in)    :: it
    real(dp),     intent(in)    :: sigma_new_nd(nr)
    real(dp),     intent(inout) :: nu_new_nd(nr)

    real(dp) :: r_cgs(nr), Sigma_cgs(nr), OmegaK_cgs(nr)
    real(dp) :: H_geom(nr), Y(nr), dY(nr), Qirr_prof(nr)
    real(dp) :: T_use
    real(dp) :: H_loc, rho_loc, kappa_loc, tau_loc, nu_dim_loc
    real(dp) :: Qv, Qi_tmp, Qc
    integer(i4b) :: i
    logical :: shadow_cell(nr)

    real(dp), parameter :: tiny_sigma = 1.0e-8_dp

    call build_shadow_flags_safe(nr, is_shadow)

    ! Units for this slice
!$omp parallel do default(shared) private(i)
    do i = 1, nr
      r_cgs(i)     = r_dim(r(i))
      Sigma_cgs(i) = sigma_dim(sigma_new_nd(i))
      if (r_cgs(i) > 0.0_dp) then
        OmegaK_cgs(i) = omegaK_dim(r(i))
      else
        OmegaK_cgs(i) = 0.0_dp
      end if
    end do
!$omp end parallel do

    ! For rebuild, use the current H slice if already present, else use zeros.
    ! (If H(it,:) is not yet initialized, the first call yields Qirr~0 and
    !  you can call rebuild again after you have a consistent H.)
    H_geom(:) = 0.0_dp
    if (it >= 1) then
      H_geom(:) = H(it, :)
    end if

    call build_Qirr_profile(nr, r, r_cgs, H_geom, &
                            Y, dY, Qirr_prof, dYdXi(it, :))

    do i = 1, nr

      if (Sigma_cgs(i) <= tiny_sigma .or. OmegaK_cgs(i) <= 0.0_dp) then
        H(it,i)      = 0.0_dp
        rho(it,i)    = 0.0_dp
        kappaR(it,i) = 0.0_dp
        tauR(it,i)   = 0.0_dp
        Qvis(it,i)   = 0.0_dp
        Qrad(it,i)   = 0.0_dp
        Qirr(it,i)   = 0.0_dp
        if (it > 1) nu_new_nd(i) = nu_conv(it-1,i)
        cycle
      end if

      T_use = Tmid(it, i)
      if (T_use <= 0.0_dp) T_use = max(T_floor, 0.5_dp*(T_floor + T_ceiling))
      T_use = max(T_floor, min(T_use, T_ceiling))

      call heating_cooling_cell(r_cgs(i), Sigma_cgs(i), OmegaK_cgs(i), shadow_cell(i), T_use, &
                                H_loc, rho_loc, nu_dim_loc, kappa_loc, tau_loc, Qv, Qi_tmp, Qc, &
                                Qirr_in=Qirr_prof(i))

      H(it,i)      = H_loc
      rho(it,i)    = rho_loc
      kappaR(it,i) = kappa_loc
      tauR(it,i)   = tau_loc
      Qvis(it,i)   = Qv
      Qrad(it,i)   = Qc
      Qirr(it,i)   = Qirr_prof(i)

      if (nu0_dim > 0.0_dp) then
        nu_new_nd(i) = (nu_dim_loc / nu0_dim) * nu0_nd
      else
        if (it > 1) nu_new_nd(i) = nu_conv(it-1,i)
      end if
      nu_conv(it,i) = nu_new_nd(i)

    end do

  end subroutine rebuild_structure_from_current_T

  !============================================================
  ! Internal helpers
  !============================================================

  subroutine build_Qirr_profile(n, xi_nd, r_cgs, H_cgs, Y, dY, Qirr_prof, dY_store)
    integer(i4b), intent(in)  :: n
    real(dp),     intent(in)  :: xi_nd(n)
    real(dp),     intent(in)  :: r_cgs(n)
    real(dp),     intent(in)  :: H_cgs(n)
    real(dp),     intent(out) :: Y(n), dY(n), Qirr_prof(n)
    real(dp),     intent(out) :: dY_store(n)

    if (use_irradiation .and. rin_cgs > 0.0_dp .and. L1 > 0.0_dp) then
      call compute_Y_dYdXi(n, xi_nd, r_cgs, H_cgs, Y, dY)
      dY_store(:) = dY(:)
      call compute_Qirr_eq16(n, xi_nd, Y, dY, rin_cgs, A1, L1, Q12, beta1, beta2, Qirr_prof)
    else
      Y(:)         = 0.0_dp
      dY(:)        = 0.0_dp
      dY_store(:)  = 0.0_dp
      Qirr_prof(:) = 0.0_dp
    end if
  end subroutine build_Qirr_profile

  subroutine solve_T_implicit_cell(r_cgs, Sigma_cgs, OmegaK_cgs, shadow, dt_phys, T_old, Qirr_in, T_new)
    ! Implicit Euler:
    !   R(T) = (T - T_old)/dt_phys - (Qvis(T)+Qirr - Qrad(T)) / (Sigma*cV) = 0
    ! Numerical derivative dR/dT by finite difference.
    use, intrinsic :: ieee_arithmetic
    real(dp), intent(in)  :: r_cgs, Sigma_cgs, OmegaK_cgs, dt_phys, T_old, Qirr_in
    logical,  intent(in)  :: shadow
    real(dp), intent(out) :: T_new

    integer(i4b), parameter :: itmax = 30
    real(dp),    parameter :: tolrel = 1.0e-3_dp
    real(dp),    parameter :: tiny  = 1.0e-99_dp

    real(dp) :: T, R, Rp, dRdT, dT, eps
    real(dp) :: cv
    integer(i4b) :: it

    ! Constant cV per unit mass (ideal monoatomic, gamma=5/3): cV = (3/2) kB / (mu mp)
    cv = 1.5_dp * kb / (mu * mp)

    T = max(T_floor, min(T_old, T_ceiling))

    do it = 1, itmax

      call eval_residual(r_cgs, Sigma_cgs, OmegaK_cgs, shadow, dt_phys, T_old, Qirr_in, cv, T, R)

      ! Convergence in relative temperature update sense will be handled below
      eps = max(1.0e-2_dp*T, 1.0_dp)
      call eval_residual(r_cgs, Sigma_cgs, OmegaK_cgs, shadow, dt_phys, T_old, Qirr_in, cv, T+eps, Rp)
      dRdT = (Rp - R) / max(eps, 1.0e-99_dp)

      if (abs(dRdT) < tiny) exit

      dT = -R / dRdT

      ! Damping for robustness
      dT = max(min(dT, 0.5_dp*T), -0.5_dp*T)

      T = T + dT
      T = max(T_floor, min(T, T_ceiling))

      if (.not. ieee_is_finite(T)) then
        T = max(T_floor, min(T_old, T_ceiling))
        exit
      end if

      if (abs(dT)/max(T,1.0_dp) < tolrel) exit
    end do

    T_new = max(T_floor, min(T, T_ceiling))
  end subroutine solve_T_implicit_cell

  subroutine eval_residual(r_cgs, Sigma_cgs, OmegaK_cgs, shadow, dt_phys, T_old, Qirr_in, cv, T, R)
    real(dp), intent(in)  :: r_cgs, Sigma_cgs, OmegaK_cgs, dt_phys, T_old, Qirr_in, cv, T
    logical,  intent(in)  :: shadow
    real(dp), intent(out) :: R

    real(dp) :: H_loc, rho_loc, nu_dim, kappa_loc, tau_loc, Qv, Qi_tmp, Qc
    real(dp) :: rhs

    call heating_cooling_cell(r_cgs, Sigma_cgs, OmegaK_cgs, shadow, T, &
                              H_loc, rho_loc, nu_dim, kappa_loc, tau_loc, Qv, Qi_tmp, Qc, &
                              Qirr_in=Qirr_in)

    rhs = (Qv + Qirr_in - Qc) / max(Sigma_cgs*cv, 1.0e-99_dp)
    R   = (T - T_old)/max(dt_phys,1.0e-99_dp) - rhs
  end subroutine eval_residual

  subroutine build_shadow_flags_safe(nr_local, is_shadow_out)
    use mod_global, only : use_irradiation
    integer(i4b), intent(in)  :: nr_local
    logical,      intent(out) :: is_shadow_out(nr_local)

    if (.not. use_irradiation) then
      is_shadow_out(:) = .true.
    else
      ! Placeholder: no shadow anywhere
      is_shadow_out(:) = .false.
    end if
  end subroutine build_shadow_flags_safe

  subroutine apply_radial_Tdiff_implicit(n, dt, r_cgs, Sigma_cgs, H_cgs, kappa, cv, T_inout)
  ! Implicit step for radial radiative diffusion of midplane temperature.
  !
  ! Solves:
  !   (T^{n+1} - T^*)/dt = (1/(Sigma*cv)) * (1/r) d/dr [ r K dT/dr ]
  !
  ! where (vertically integrated conductivity)
  !   K = (64*sbc/3) * H^2 * T^3 / (kappa * Sigma)
  !
  ! Boundary condition (default): zero radial heat flux at both ends
  !   dT/dr = 0  => no conductive energy flow across boundaries.
  !
  ! Notes:
  ! - This step should be applied AFTER you update T by local source terms,
  !   i.e., T_inout holds T^* on entry, and is overwritten by T^{n+1}.
  ! - Use CGS units consistently.

  use kind_params, only : dp, i4b
  use constants, only : sbc
  implicit none

  integer(i4b), intent(in)    :: n
  real(dp),     intent(in)    :: dt
  real(dp),     intent(in)    :: r_cgs(n), Sigma_cgs(n), H_cgs(n), kappa(n)
  real(dp),     intent(in)    :: cv
  real(dp),     intent(inout) :: T_inout(n)

  integer(i4b) :: i
  real(dp) :: rE(n+1), drC(n), drE(n-1)
  real(dp) :: Kc(n), Ke(n-1)
  real(dp) :: a(n), b(n), c(n), rhs(n)
  real(dp) :: tiny, Ti, sigi, Hi, kappai
  real(dp) :: facL, facR, denom, w
  real(dp) :: alphaL, alphaR
  real(dp), parameter :: kappa_floor = 1.0e-3_dp 
  real(dp), parameter :: Sigma_floor = 1.0e-8_dp  
  real(dp), parameter :: Kcap = 1.0e27_dp

  tiny = 1.0e-99_dp

  !---------------------------------------------
  ! Build cell edges rE and spacings for a general grid
  ! rE(i)   : edge between cell i-1 and i (i=2..n)
  ! rE(1)   : inner boundary edge
  ! rE(n+1) : outer boundary edge
  !
  ! Here we set interior edges as geometric means (robust on log grids).
  ! Boundaries are extrapolated geometrically.
  !---------------------------------------------
  do i = 2, n
     rE(i) = sqrt(max(r_cgs(i-1),tiny) * max(r_cgs(i),tiny))
  end do
  ! Extrapolate boundary edges assuming local log-spacing
  rE(1)   = max(r_cgs(1),tiny)**2 / max(rE(2),tiny)
  rE(n+1) = max(r_cgs(n),tiny)**2 / max(rE(n),tiny)

  do i = 1, n
     drC(i) = rE(i+1) - rE(i)
  end do
  do i = 1, n-1
     drE(i) = r_cgs(i+1) - r_cgs(i)
  end do

  !---------------------------------------------
  ! Build K at cell centers using T^* (current T_inout)
  ! Kc = (64*sbc/3) * H^2 * T^3 / (kappa*Sigma)
  !---------------------------------------------
  do i = 1, n
     Ti     = max(T_inout(i), 1.0_dp)   ! prevent T^3 blow-ups near 0
     sigi   = max(Sigma_cgs(i), tiny)
     if (sigi <= Sigma_floor) then
        Kc(i) = 0.0_dp
     else
        Hi     = max(H_cgs(i), tiny)
        kappai = max(kappa(i), kappa_floor)
        Kc(i)  = (64.0_dp*sbc/3.0_dp) * (Hi*Hi) * (Ti**3) / (kappai*sigi)
        Kc(i)  = min(Kc(i), Kcap)
     end if
  end do

  ! Edge K: harmonic mean is more stable across sharp jumps (kappa/ Sigma transitions)
  do i = 1, n-1
     denom = (1.0_dp/max(Kc(i),tiny)) + (1.0_dp/max(Kc(i+1),tiny))
     Ke(i) = 2.0_dp / max(denom, tiny)
  end do

  !---------------------------------------------
  ! Assemble tridiagonal system for T^{n+1}
  !
  ! Discretization of: (1/r) d/dr ( r K dT/dr )
  ! Using finite volume:
  !   d/dr(...) ~ (F_{i+1/2} - F_{i-1/2}) / (r_i * dV_i)
  ! with F_{i+1/2} = rE(i+1) * Ke(i) * (T_{i+1}-T_i) / drE(i)
  !
  ! Overall diffusion operator:
  !   RHS_diff(i) = (1/(Sigma_i*cv)) * (1/r_i) * (F_{i+1/2} - F_{i-1/2}) / drC(i)
  !
  ! Implicit Euler:
  !   T^{n+1} - dt*RHS_diff(T^{n+1}) = T^*
  !---------------------------------------------
  rhs(:) = T_inout(:)
  a(:) = 0.0_dp
  b(:) = 1.0_dp
  c(:) = 0.0_dp

  do i = 2, n-1
     sigi = max(Sigma_cgs(i), tiny)

     facL = (dt/(sigi*cv)) * (1.0_dp/max(r_cgs(i),tiny)) * (1.0_dp/max(drC(i),tiny))
     facR = facL

     ! Contributions from left and right edge fluxes
     alphaL = facL * ( rE(i)   * Ke(i-1) / max(drE(i-1),tiny) )
     alphaR = facR * ( rE(i+1) * Ke(i)   / max(drE(i),tiny)   )

     a(i) = -alphaL
     c(i) = -alphaR
     b(i) = 1.0_dp + alphaL + alphaR
  end do

  !---------------------------------------------
  ! Neumann BC: dT/dr = 0 at inner/outer boundaries
  !
  ! Easiest robust implementation in FV:
  !   Set boundary fluxes to zero:
  !     F_{1-1/2} = 0 and F_{n+1/2} = 0
  ! which removes coupling beyond domain.
  !
  ! This yields one-sided operator at i=1 and i=n.
  !---------------------------------------------
  ! i = 1
  i = 1
  sigi = max(Sigma_cgs(i), tiny)
  facR = (dt/(sigi*cv)) * (1.0_dp/max(r_cgs(i),tiny)) * (1.0_dp/max(drC(i),tiny))
  alphaR = facR * ( rE(i+1) * Ke(i) / max(drE(i),tiny) )
  a(i) = 0.0_dp
  c(i) = -alphaR
  b(i) = 1.0_dp + alphaR

  ! i = n
  i = n
  sigi = max(Sigma_cgs(i), tiny)
  facL = (dt/(sigi*cv)) * (1.0_dp/max(r_cgs(i),tiny)) * (1.0_dp/max(drC(i),tiny))
  alphaL = facL * ( rE(i) * Ke(i-1) / max(drE(i-1),tiny) )
  a(i) = -alphaL
  c(i) = 0.0_dp
  b(i) = 1.0_dp + alphaL

  !---------------------------------------------
  ! Solve tridiagonal system: a(i) T_{i-1} + b(i) T_i + c(i) T_{i+1} = rhs(i)
  ! Thomas algorithm (in-place).
  !---------------------------------------------
  call solve_tridiag_thomas(n, a, b, c, rhs, T_inout)

contains

  subroutine solve_tridiag_thomas(n, a, b, c, d, x)
    use kind_params, only : dp, i4b
    implicit none
    integer(i4b), intent(in) :: n
    real(dp), intent(in)    :: a(n), b(n), c(n), d(n)
    real(dp), intent(out)   :: x(n)

    integer(i4b) :: i
    real(dp) :: cp(n), dpv(n), denom, tiny2
    tiny2 = 1.0e-99_dp

    denom = b(1)
    cp(1)  = c(1) / max(denom, tiny2)
    dpv(1) = d(1) / max(denom, tiny2)

    do i = 2, n
       denom = b(i) - a(i)*cp(i-1)
       cp(i)  = c(i) / max(denom, tiny2)
       dpv(i) = (d(i) - a(i)*dpv(i-1)) / max(denom, tiny2)
    end do

    x(n) = dpv(n)
    do i = n-1, 1, -1
       x(i) = dpv(i) - cp(i)*x(i+1)
    end do
  end subroutine solve_tridiag_thomas

end subroutine apply_radial_Tdiff_implicit

end module disk_energy_pde_mod

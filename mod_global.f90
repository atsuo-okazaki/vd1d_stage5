!==============================================================
! mod_global.f90
!   Global parameters and shared simulation state
!
!   - The diffusion solver uses dimensionless variables:
!       r'     = r / R0
!       t'     = t / t0
!       dt'    = dt / t0
!       Sigma' = Sigma / sigma_init
!
!   - Physical scales (R0, t0, sigma_init, Temp0, M_star, etc.)
!     are defined here and used through units_disk_mod.
!
!   - Physical diagnostic fields (Tmid, rho, H, kappaR)
!     are stored in CGS units.
!==============================================================
module mod_global
  use kind_params, only : dp, i4b
  use constants,  only : gg, msun
  implicit none

  !------------------------------------------------------------
  ! Numerical sizes (set during initialization)
  !------------------------------------------------------------
  integer(i4b) :: nr = 0      ! number of radial zones
  integer(i4b) :: nt = 0      ! number of time steps

  !------------------------------------------------------------
  ! Dimensionless simulation time window
  !   t_nd is advanced in main loop; these define the intended run range.
  !------------------------------------------------------------
  real(dp) :: t_sim_start = 0.0_dp
  real(dp) :: t_sim_end   = 0.0_dp
  real(dp) :: t_sim_end_eff = 0.0_dp  ! actual end time after rounding by dt and nt

  !------------------------------------------------------------
  ! Dimensionless radial domain
  !------------------------------------------------------------
  real(dp) :: rin_nd  = 0.0_dp
  real(dp) :: rout_nd = 0.0_dp

  !------------------------------------------------------------
  ! Dimensionless grid and solution arrays
  !------------------------------------------------------------
  real(dp), allocatable :: r(:)
  real(dp), allocatable :: sigmat(:,:)

  !------------------------------------------------------------
  ! Physical disk quantities (diagnostic)
  !------------------------------------------------------------
  real(dp), allocatable :: nu(:)
  real(dp), allocatable :: nu_conv(:,:)
  real(dp), allocatable :: Tmid(:,:)
  real(dp), allocatable :: rho(:,:)
  real(dp), allocatable :: H(:,:)
  real(dp), allocatable :: kappaR(:,:)
  real(dp), allocatable :: kappa_planck(:,:)
  real(dp), allocatable :: tauR(:,:)
  integer(i4b), allocatable :: k_iter(:)
  integer(i4b), allocatable :: m_iter(:)
  logical, allocatable :: is_shadow(:,:)   ! (nt,nr) if you store; else remove

  !---------------------------------------------------------
  ! current state vectors advanced by substeps
  !---------------------------------------------------------
  real(dp), allocatable :: sigma_cur(:), nu_cur(:)
  real(dp), allocatable :: Tmid_cur(:), H_cur(:), rho_cur(:), kappa_cur(:), kappa_planck_cur(:), tau_cur(:)
  real(dp), allocatable :: Qvis_cur(:), Qrad_cur(:), Qirr_cur(:), dYdXi_cur(:)
  logical, allocatable  :: shadow_cur(:)
  integer(i4b) :: k_iter_cur, m_iter_cur

  !------------------------------------------------------------
  ! Dimensionless time-stepping control
  !------------------------------------------------------------
  real(dp) :: dt    = 0.0_dp
  real(dp) :: alpha = 0.5_dp
  real(dp) :: dt_out         ! output interval in nd units
  real(dp) :: dt_seed        ! suggested dt for next substep (nd)
  real(dp) :: t_nd  = 0.0_dp ! current output time (nd)

  !------------------------------------------------------------
  ! Physical scaling parameters
  !------------------------------------------------------------
  real(dp) :: R0                    ! radial scale [cm]
  real(dp) :: t0     = 0.0_dp       ! time scale [s], computed in init_units
  real(dp) :: sigma_init = 0.0_dp   ! surface-density scale [g/cm^2]
  real(dp) :: Temp0  = 0.0_dp       ! temperature scale [K]

  !--------------------------------------------------------
  ! Stellar parameters for the Be star
  !--------------------------------------------------------
  real(dp) :: M_star      ! stellar mass [g]
  real(dp) :: R_star      ! stellar radius [cm]
  real(dp) :: Teff_star   ! stellar effective temperature [K]
  real(dp) :: L_star      ! stellar luminosity [erg/s]
  real(dp) :: q           ! binary mass ratio (0<= q <= 1)

  !--------------------------------------------------------
  ! Irradiation parameters
  !--------------------------------------------------------
  real(dp) :: eta_acc = 1.0_dp/6.0_dp ! Efficiency of accretion luminosity
  real(dp) :: f_edd_cap = 0.0_dp      ! Upper limit of mdot_inc_edd for accretion

  !--------------------------------------------------------
  ! Wind mass-loss parameters
  !--------------------------------------------------------
  real(dp) :: mdot_w_msunyr ! wind mass-loss rate [g/s]
  real(dp) :: vinf_w        ! wind terminal velocity [cm/s]
  real(dp) :: beta_w        ! beta-law index for wind velocity
  real(dp) :: f_rho_wind    ! parameter

  !-----------------------------------------------------------
  ! Density threshold for applying S03 opacity table [g/cm^3]
  !-----------------------------------------------------------
  real(dp), parameter :: rho_cut = 1.0e-17_dp

  !------------------------------------------------------------
  ! Viscosity (physical) parameters
  !   alphaSS : Shakura-Sunyaev alpha
  !   hr0     : aspect ratio H/R at r' = 1
  !   delta   : power-law index for viscosity profile
  !   nu0_dim : viscosity at r' = 1 in CGS units [cm^2/s]
  !------------------------------------------------------------
  real(dp) :: alphaSS = 0.1_dp
  real(dp) :: hr0     = 0.01_dp
  real(dp) :: delta   = 0.0_dp
  real(dp) :: nu0_nd  = 0.0_dp
  real(dp) :: nu0_dim = 0.0_dp   ! set in init_units()

  real(dp) :: fwhm    = 0.0_dp

  !------------------------------------------------------------
  ! Opacity parameters
  !------------------------------------------------------------
  real(dp), parameter :: kappa_es = 0.34_dp      ! electron scattering opacity [cm^2/g]
  real(dp), parameter :: kappa0   = 5.24e24_dp   ! Kramers coefficient [cm^5 g^-2 K^(7/2)]

  !---------------------------------------------------------
  ! Physics / mode switches
  !---------------------------------------------------------
  logical :: use_energy_balance  = .true.  ! .true.: Sigma-T-nu iteration
                                           ! .false.: isothermal nu(r)
  logical :: use_inflow          = .true.  ! .true.: mass injection
  logical :: use_be_decretion    = .false. ! .true.: Be decretion disk BC
  logical :: use_wind_truncation = .false. ! .true.: wind-driven ablation
  logical :: use_irradiation     = .false. ! .true.: irradiation heating
  logical :: use_irradiation_delay = .false. ! .true.: delay L_irr by tau_irr_lag (accretion history)
  logical :: use_finite_irradiation_source = .false. ! .true.: multiply Qirr by (2/pi)*arctan(R_star/r)

  !---------------------------------------------------------
  ! Energy PDE switch (MVP): evolve Tmid by BE ODE per radius
  !---------------------------------------------------------
  logical :: use_energy_pde = .false.   ! .true.: BE temperature evolution (no Sigma-structure iteration)

  !---------------------------------------------------------
  ! Irradiation delay (for inner-edge Mdot-proportional component)
  !   tau_irr_lag_mode : 'explicit' = fixed delay, 'viscous' = t_visc at inner edge
  !   tau_irr_lag_nd   : explicit: delay in t/t0; viscous: prefactor for t_visc
  !---------------------------------------------------------
  character(len=16) :: tau_irr_lag_mode = 'explicit'
  real(dp) :: tau_irr_lag_nd = 0.0_dp

  !---------------------------------------------------------
  ! Outer inflow (mass supply at outer radii)
  ! All variables are dimensionless unless noted.
  !---------------------------------------------------------
  real(dp) :: t_in_start = 0.0_dp
  real(dp) :: t_in_end   = 0.0_dp
  real(dp) :: rinj_min   = 0.0_dp
  real(dp) :: rinj_max   = 0.0_dp

  real(dp) :: mdot_inj_edd  = 0.0_dp ! input: Edd-normalized accretion supply
  real(dp) :: mdot_inj_nd   = 0.0_dp ! Dimensionless accretion rate 
                                     ! normalized by sigma_init*R0^2/t0
  real(dp) :: mdot_inj_phys = 0.0_dp ! accretion rate [g/s]
  real(dp) :: mdot_inj_msunyr = 0.0_dp ! accretion rate [msun/yr]

  !---------------------------------------------------------
  ! Inner boundary BC for mass flux:
  !   0 : open       (mass can flow out freely)
  !   1 : no-outflow (forbid outward flux)
  !   2 : zero-flux  (no flux across outer boundary)
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Inner boundary condition API (mode-independent)
  !---------------------------------------------------------
  integer(i4b), parameter :: INNER_TORQUE_FREE = 0
  integer(i4b), parameter :: INNER_FIXED_MDOT  = 1 ! decretion only (for now)
  integer(i4b), parameter :: INNER_ZERO_FLUX   = 2 ! accretion onto a solid surfae

  integer(i4b) :: inner_bc_type = INNER_TORQUE_FREE

  ! Inner boundary mass flux (PDE-normalized, inward positive)
  real(dp) :: mdot_inner_nd   = 0.0_dp
  real(dp) :: mdot_inner_phys = 0.0_dp
  real(dp) :: mdot_inner_edd  = 0.0_dp   ! diagnostic only

  !---------------------------------------------------------
  ! Outer boundary BC for mass flux:
  !   0 : open       (mass can flow out freely)
  !   1 : no-outflow (forbid outward flux)
  !   2 : zero-flux  (no flux across outer boundary)
  ! For Be decretion disk we want 0.
  !---------------------------------------------------------
  integer(i4b) :: outer_bc_type

  !---------------------------------------------------------
  ! Isothermal viscosity profile:
  !   nu(r) = nu0_nd * r**p_nu_isothermal
  !   For Keplerian disk with constant sound speed, p = 3/2.
  !---------------------------------------------------------
  real(dp) :: p_nu_isothermal = 1.5_dp

  !---------------------------------------------------------
  ! Disk edge location
  !---------------------------------------------------------
  real(dp), allocatable :: r_edge(:)      ! disk edge radius [cm]
  integer(i4b), allocatable :: i_edge(:)  ! index of edge cell

  !---------------------------------------------------------
  ! Heating and cooling rates
  !---------------------------------------------------------
  real(dp), allocatable :: Qvis(:,:), Qrad(:,:), Qirr(:,:), dYdXi(:,:)

contains

  !----------------------------------------------------------
  ! Initialize physical scales (t0, nu0_dim)
  !
  !   Omega0  = sqrt(G M / R0^3)
  !   t0      = 1 / Omega0
  !   H0      = hr0 * R0
  !   nu0_dim = 2.0_dp / 3.0_dp * alphaSS * (H0 * Omega0)^2 / Omega0
  !
  ! This must be called after reading M_star, R0, hr0, alphaSS.
  !----------------------------------------------------------
  subroutine init_units()
    use constants, only : gg, pi, cc
    implicit none
    real(dp) :: Omega0, H0

    Omega0 = sqrt( gg * M_star / R0**3 )
    t0     = 1.0_dp / Omega0
    sigma_init = 4.0_dp * pi * gg * M_star / (kappa_es * cc) / (R0*R0) &
         * t0

    H0  = hr0 * R0
    nu0_dim = 2.0_dp / 3.0_dp * alphaSS * (H0 * Omega0)**2 / Omega0
    ! => nu0 = 2/3 * alphaSS * H0^2 * Omega0
  end subroutine init_units

  !----------------------------------------------------------
  ! Allocate global arrays for given nr, nt
  !----------------------------------------------------------
  subroutine allocate_global(nr_in, nt_in)
    implicit none
    integer(i4b), intent(in) :: nr_in, nt_in
    integer(i4b) :: ierr

    nr = nr_in
    nt = nt_in

    if (allocated(r))      deallocate(r)
    if (allocated(nu))     deallocate(nu)
    if (allocated(nu_conv))     deallocate(nu_conv)
    if (allocated(sigmat)) deallocate(sigmat)
    if (allocated(Tmid))   deallocate(Tmid)
    if (allocated(H))      deallocate(H)
    if (allocated(rho))    deallocate(rho)
    if (allocated(kappaR)) deallocate(kappaR)
    if (allocated(kappa_planck)) deallocate(kappa_planck)
    if (allocated(tauR))   deallocate(tauR)
    if (allocated(Qvis)) deallocate(Qvis)
    if (allocated(Qrad)) deallocate(Qrad)
    if (allocated(Qirr)) deallocate(Qirr)
    if (allocated(dYdXi)) deallocate(dYdXi)
    if (allocated(is_shadow)) deallocate(is_shadow)
    if (allocated(k_iter)) deallocate(k_iter)
    if (allocated(m_iter)) deallocate(m_iter)
    if (allocated(r_edge)) deallocate(r_edge)
    if (allocated(i_edge)) deallocate(i_edge)

    allocate(r(nr), nu(nr), nu_conv(nt,nr), sigmat(nt,nr), &
         Tmid(nt,nr), H(nt,nr), rho(nt,nr), kappaR(nt,nr), kappa_planck(nt,nr), &
         tauR(nt,nr), Qvis(nt,nr), Qrad(nt,nr),            &
         Qirr(nt,nr), dYdXi(nt,nr), is_shadow(nt,nr),         &
         k_iter(nt), m_iter(nt), r_edge(nt), i_edge(nt), stat=ierr)

    if (ierr /= 0) then
       write(*,'(a,i6)') 'FATAL: allocation failed, stat = ', ierr
       stop
    end if
  end subroutine allocate_global

  subroutine alloc_cur_state(nr_in)
    integer(i4b), intent(in) :: nr_in
    allocate(sigma_cur(nr_in), nu_cur(nr_in))
    allocate(Tmid_cur(nr_in), H_cur(nr_in), rho_cur(nr_in), kappa_cur(nr_in), kappa_planck_cur(nr_in), tau_cur(nr_in))
    allocate(Qvis_cur(nr_in), Qrad_cur(nr_in), Qirr_cur(nr_in), dYdXi_cur(nr_in))
    allocate(shadow_cur(nr_in))
  end subroutine alloc_cur_state

  !----------------------------------------------------------
  ! Deallocate (optional)
  !----------------------------------------------------------
  subroutine finalize_global()
    implicit none
    if (allocated(r))      deallocate(r)
    if (allocated(nu))     deallocate(nu)
    if (allocated(sigmat)) deallocate(sigmat)
    if (allocated(Tmid))   deallocate(Tmid)
    if (allocated(H))      deallocate(H)
    if (allocated(rho))    deallocate(rho)
    if (allocated(kappaR)) deallocate(kappaR)
    if (allocated(kappa_planck)) deallocate(kappa_planck)
    if (allocated(tauR)) deallocate(tauR)
    if (allocated(Qvis)) deallocate(Qvis)
    if (allocated(Qrad)) deallocate(Qrad)
    if (allocated(Qirr)) deallocate(Qirr)
    if (allocated(dYdXi)) deallocate(dYdXi)
    if (allocated(is_shadow)) deallocate(is_shadow)
    if (allocated(k_iter)) deallocate(k_iter)
    if (allocated(m_iter)) deallocate(m_iter)
    if (allocated(r_edge)) deallocate(r_edge)
    if (allocated(i_edge)) deallocate(i_edge)
  end subroutine finalize_global

  !----------------------------------------------------------
  ! Print global info (debugging)
  !----------------------------------------------------------
  subroutine print_global_info()
    implicit none
    write(*,'(a,i6,a,i9)') 'nr = ', nr, ', nt = ', nt
    write(*,'(a,1pe12.4)') 'dt (dimensionless) = ', dt
    write(*,'(a,1pe12.4,a,1pe12.4)') 'rin_nd = ', rin_nd, ', rout_nd = ', rout_nd
    write(*,'(a,1pe12.4)') 'R0 [cm]   = ', R0
    write(*,'(a,1pe12.4)') 't0 [s]    = ', t0
    write(*,'(a,1pe12.4)') 'sigma_init    = ', sigma_init
    write(*,'(a,1pe12.4)') 'Temp0 [K] = ', Temp0
  end subroutine print_global_info

end module mod_global

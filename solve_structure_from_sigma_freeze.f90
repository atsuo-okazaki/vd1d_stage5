subroutine solve_structure_from_sigma_freeze( nr, r, sigma_in, T_seed,         &
     nu_seed, H_seed, rho_seed, kappa_seed, tau_seed,                         &
     Qvis_seed, Qrad_seed, Qirr_seed, dYdXi_seed, shadow_seed,                &
     tau_freeze,                                                              &
     nu_out, T_out, H_out, rho_out, kappa_out, tau_out,                       &
     Qvis_out, Qrad_out, Qirr_out, dYdXi_out, shadow_out )

  use kind_params, only : dp, i4b
  use disk_energy_mod, only : solve_structure_from_sigma
  implicit none

  integer(i4b), intent(in) :: nr
  real(dp),     intent(in) :: r(nr), sigma_in(nr), T_seed(nr)
  real(dp),     intent(in) :: nu_seed(nr), H_seed(nr), rho_seed(nr)
  real(dp),     intent(in) :: kappa_seed(nr), tau_seed(nr)
  real(dp),     intent(in) :: Qvis_seed(nr), Qrad_seed(nr), Qirr_seed(nr)
  real(dp),     intent(in) :: dYdXi_seed(nr)
  logical,      intent(in) :: shadow_seed(nr)

  real(dp),     intent(in) :: tau_freeze

  real(dp),     intent(out) :: nu_out(nr), T_out(nr), H_out(nr), rho_out(nr)
  real(dp),     intent(out) :: kappa_out(nr), tau_out(nr)
  real(dp),     intent(out) :: Qvis_out(nr), Qrad_out(nr), Qirr_out(nr)
  real(dp),     intent(out) :: dYdXi_out(nr)
  logical,      intent(out) :: shadow_out(nr)

  integer(i4b) :: i
  ! 1-point buffers for the underlying solver
  real(dp) :: r1(1), sigma1(1), Tseed1(1)
  real(dp) :: nu1(1), T1(1), H1(1), rho1(1)
  real(dp) :: kappa1(1), tau1(1)
  real(dp) :: Qvis1(1), Qrad1(1), Qirr1(1)
  real(dp) :: dYdXi1(1)
  logical  :: shadow1(1)

  !--------------------------------------------
  ! Default: copy seeds (freeze behavior)
  !--------------------------------------------
  nu_out(:)     = nu_seed(:)
  T_out(:)      = T_seed(:)
  H_out(:)      = H_seed(:)
  rho_out(:)    = rho_seed(:)
  kappa_out(:)  = kappa_seed(:)
  tau_out(:)    = tau_seed(:)
  Qvis_out(:)   = Qvis_seed(:)
  Qrad_out(:)   = Qrad_seed(:)
  Qirr_out(:)   = Qirr_seed(:)
  dYdXi_out(:)  = dYdXi_seed(:)
  shadow_out(:) = shadow_seed(:)

  !--------------------------------------------
  ! Update only cells with tau_seed < tau_freeze
  !--------------------------------------------
  do i = 1, nr
     if (tau_seed(i) < tau_freeze) then

        r1(1)     = r(i)
        sigma1(1) = sigma_in(i)
        Tseed1(1) = T_seed(i)

        call solve_structure_from_sigma( 1, r1, sigma1, Tseed1,                &
             nu1, T1, H1, rho1,                                               &
             kappa1, tau1,                                                    &
             Qvis1, Qrad1, Qirr1,                                              &
             dYdXi1, shadow1 )

        ! Commit outputs for this cell
        nu_out(i)     = nu1(1)
        T_out(i)      = T1(1)
        H_out(i)      = H1(1)
        rho_out(i)    = rho1(1)
        kappa_out(i)  = kappa1(1)
        tau_out(i)    = tau1(1)
        Qvis_out(i)   = Qvis1(1)
        Qrad_out(i)   = Qrad1(1)
        Qirr_out(i)   = Qirr1(1)
        dYdXi_out(i)  = dYdXi1(1)
        shadow_out(i) = shadow1(1)

     end if
  end do

end subroutine solve_structure_from_sigma_freeze

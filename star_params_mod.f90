module star_params_mod
  use kind_params, only : dp
  implicit none

  ! Model type label: 'be', 'bh_disk', etc.
  character(len=16) :: model_type

  ! Stellar / central object parameters in convenient units
  real(dp) :: M_star_msun  ! mass [Msun]
  real(dp) :: R_star_rsun  ! radius [Rsun]
  real(dp) :: Teff_star_in ! effective temperature [K]

  ! Scale radius options
  character(len=16) :: r0_mode       ! 'star', 'rg', 'explicit'
  real(dp)          :: r0_over_Rstar ! if r0_mode='star'
  real(dp)          :: r0_over_Rg    ! if r0_mode='rg'
  real(dp)          :: R0_input      ! [cm] if r0_mode='explicit'

  ! Namelists
  namelist /star_params/  model_type, M_star_msun, R_star_rsun, Teff_star_in
  namelist /scale_params/ r0_mode, r0_over_Rstar, r0_over_Rg, R0_input

contains

  subroutine set_default_star_params()
    ! Set default values for stellar and scale parameters.
    model_type    = 'be'
    M_star_msun   = 19.0_dp
    R_star_rsun   = 8.0_dp
    Teff_star_in  = 26000.0_dp

    r0_mode       = 'star'
    r0_over_Rstar = 1.0_dp
    r0_over_Rg    = 10.0_dp
    R0_input      = 0.0_dp
  end subroutine set_default_star_params

end module star_params_mod

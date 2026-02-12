module radiation_params_mod
  !! Shared physical and numerical parameters for the Be disk
  !! thermal / radiative modules.

  use kind_params, only : dp
  implicit none

  !--------------------------------------------------------------
  ! Optical–depth / irradiation parameters
  !--------------------------------------------------------------
  real(dp), parameter :: tau_min     = 1.0e-6_dp
  !! Minimum Rosseland optical depth used in regularization
  !! to avoid division by zero or extremely thin limits.

  !--------------------------------------------------------------
  ! Thermal–solver thresholds
  !--------------------------------------------------------------
  real(dp), parameter :: Sigma_min_energy = 1.0e-8_dp
  !! Below this surface density [g/cm^2] we skip the full
  !! energy balance and fall back to a simple prescription.

  real(dp), parameter :: T_floor    = 1.0e1_dp
  !! Lower bound for the midplane temperature [K] (safety).

  real(dp), parameter :: T_ceiling = 1.0e8_dp
  !! Upper bound for the midplane temperature [K] (safety).

end module radiation_params_mod

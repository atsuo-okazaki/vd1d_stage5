!==============================================================
! constants.f90
!   Basic mathematical / physical / astronomical constants
!   (cgs units where applicable)
!==============================================================
module constants
  use kind_params, only : dp
  implicit none

  ! Mathematical constants
  real(dp), parameter :: pi   = acos(-1.0_dp)

  ! Physical constants (cgs)
  real(dp), parameter :: gg   = 6.67430e-8_dp   ! gravitational constant [cm^3 g^-1 s^-2]
  real(dp), parameter :: cc   = 2.99792458e10_dp ! speed of light [cm s^-1]
  real(dp), parameter :: Rg   = 8.314462618e7_dp     ! gas constant [erg mol^-1 K^-1]
  real(dp), parameter :: kb   = 1.380649e-16_dp ! Boltzmann constant [erg K^-1]
  real(dp), parameter :: sbc  = 5.6704e-5_dp    ! Stefan-Boltzmann constant
  real(dp), parameter :: mH   = 1.6735575e-24_dp  ! hydrogen mass [g]
  real(dp), parameter :: mp   = 1.6726e-24_dp     ! proton mass [g]

  ! Astronomical constants (cgs)
  real(dp), parameter :: Msun = 1.98847e33_dp     ! solar mass [g]
  real(dp), parameter :: Rsun = 6.957e10_dp       ! solar radius [cm]
  real(dp), parameter :: Lsun = 3.828e33_dp       ! solar luminosity [erg/s]
  real(dp), parameter :: AU   = 1.495978707e13_dp ! astronomical unit [cm]
  real(dp), parameter :: year = 3.15576e7_dp      ! year [s]
  real(dp), parameter :: day  = 8.64e4_dp         ! day [s]
  real(dp), parameter :: mu   = 0.62_dp           ! mean molecular weight

end module constants

module analytic_lbp_mod
  use kind_params, only : dp
  implicit none
contains

  ! Return the analytic Lynden-Bell & Pringle ring solution
  ! for constant kinematic viscosity nu0_nd.
  !
  ! Input:
  !   r      : radius
  !   t      : time
  !   r0     : initial ring radius
  !   nu0_nd : constant viscosity
  !   M0     : total initial disk mass
  !
  ! Output:
  !   sigma_analytic : surface density at (r,t)
  !
  function sigma_lbp_ring(r, t, r0, nu0_nd, M0) result(sigma)
    use bessel_i_mod, only : bessik
    real(dp), intent(in) :: r, t, r0, nu0_nd, M0
    real(dp) :: sigma
    real(dp) :: x, tau, arg, prefac, bessel_val
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    real(dp) :: ri, rk, rip, rkp

    if (t <= 0.0_dp) then
       ! At t=0 the solution is a delta-function ring; in practice
       ! you should compare to numeric solution only for t > 0.
       sigma = 0.0_dp
       return
    end if

    x   = r / r0
    tau = 12.0_dp * nu0_nd * t / (r0*r0)

    ! Avoid division by zero
    if (tau <= 0.0_dp) then
       sigma = 0.0_dp
       return
    end if

    arg   = 2.0_dp * x / tau
    ! You need an implementation of I_{1/4}(z) here.
    ! For example from a special-function library or Numerical Recipes.
    call bessik(arg, 0.25_dp, ri, rk, rip, rkp)
    bessel_val = ri

    Prefac = (M0 / (pi * r0*r0)) * tau**(-1.0_dp) * x**(-0.25_dp)
    sigma  = prefac * exp( -(1.0_dp + x*x)/tau ) * bessel_val
  end function sigma_lbp_ring

end module analytic_lbp_mod

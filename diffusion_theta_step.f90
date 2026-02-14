subroutine diffusion_theta_step(nr, r, nu, Sigma_old, dt, theta, Sigma_new)
  use kind_params, only : dp, i4b
  use constants,   only : pi
  use mod_global, only  : inner_bc_type, mdot_inner_nd, INNER_FIXED_MDOT
  implicit none

  integer(i4b), intent(in)  :: nr
  real(dp),    intent(in)   :: r(nr), nu(nr)
  real(dp),    intent(in)   :: Sigma_old(nr)
  real(dp),    intent(in)   :: dt, theta
  real(dp),    intent(out)  :: Sigma_new(nr)

  real(dp) :: L_old(nr)
  real(dp) :: a(nr), b(nr), c(nr)
  real(dp) :: rhs(nr)
  integer(i4b) :: i

  real(dp) :: r_edge(nr+1), dr_cell1, A_cell1, S_in

  call compute_LSigma(nr, r, nu, Sigma_old, L_old)

!$omp parallel do default(shared) private(i)
  do i = 1, nr
     rhs(i) = Sigma_old(i) + dt*(1.0_dp - theta)*L_old(i)
  end do
!$omp end parallel do

  ! Add theta contribution of constant inner flux if applicable
  if (inner_bc_type == INNER_FIXED_MDOT) then
     r_edge(1) = r(1) - 0.5_dp*(r(2) - r(1))
     r_edge(2) = 0.5_dp*(r(1) + r(2))
     dr_cell1  = r_edge(2) - r_edge(1)

     A_cell1 = 2.0_dp*pi*r(1)*dr_cell1
     S_in = mdot_inner_nd / A_cell1
     !rhs(1) = rhs(1) + theta*dt*S_in
     rhs(1) = rhs(1) + dt*S_in
  end if

  call build_tridiag_coeff(nr, r, nu, dt, theta, Sigma_old, a, b, c)
  call solve_tridiag(nr, a, b, c, rhs, Sigma_new)

end subroutine diffusion_theta_step

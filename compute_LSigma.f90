subroutine compute_LSigma(nr, r, nu, Sigma, LSigma)
  use kind_params, only: dp, i4b
  use constants,   only: pi
  use mod_global,  only: outer_bc_type, inner_bc_type, &
                         mdot_inner_nd, mdot_inner_phys, &
                         INNER_TORQUE_FREE, INNER_FIXED_MDOT
  use units_disk_mod, only : mdot_edd_unit
  implicit none
  integer(i4b), intent(in)  :: nr
  real(dp),    intent(in)   :: r(nr), nu(nr), Sigma(nr)
  real(dp),    intent(out)  :: LSigma(nr)

  real(dp) :: r_edge(nr+1), dr_cell(nr)
  real(dp) :: mdot_edge(nr+1)
  real(dp) :: g_left, g_right, dGdr, r_face
  integer(i4b) :: i

  ! -- 1. cell edges (same as your code)
  r_edge(1) = r(1) - 0.5_dp * ( r(2) - r(1) )
  do i = 2, nr
     r_edge(i) = 0.5_dp * ( r(i-1) + r(i) )
  end do
  r_edge(nr+1) = r(nr) + 0.5_dp * ( r(nr) - r(nr-1) )

  do i = 1, nr
     dr_cell(i) = r_edge(i+1) - r_edge(i)
  end do

  ! -- 2. fluxes mdot_edge(j)

  ! inner boundary
  select case (inner_bc_type)
  case (INNER_FIXED_MDOT)
     ! prescribed mass flux at inner edge (PDE-normalized, inward positive)
     mdot_edge(1) = mdot_inner_nd

  case default  ! INNER_TORQUE_FREE
     g_left  = 0.0_dp
     g_right = nu(1) * Sigma(1) * sqrt(r(1))
     dGdr    = ( g_right - g_left ) / ( r(1) - r_edge(1) )
     r_face  = r_edge(1)
     mdot_edge(1) = 6.0_dp * pi * sqrt(r_face) * dGdr
  end select

  ! -- Mdot at the inner radius
  mdot_inner_nd   = mdot_edge(1)
  mdot_inner_phys = mdot_inner_nd * mdot_edd_unit()

  ! internal faces (same)
  do i = 2, nr
     g_left  = nu(i-1) * Sigma(i-1) * sqrt(r(i-1))
     g_right = nu(i  ) * Sigma(i  ) * sqrt(r(i  ))
     dGdr    = ( g_right - g_left ) / ( r(i) - r(i-1) )
     r_face  = 0.5_dp * ( r(i-1) + r(i) )
     mdot_edge(i) = 6.0_dp * pi * sqrt(r_face) * dGdr
  end do

  ! outer boundary + outer BC (same)
  g_left  = nu(nr) * Sigma(nr) * sqrt(r(nr))
  g_right = 0.0_dp
  dGdr    = ( g_right - g_left ) / ( r_edge(nr+1) - r(nr) )
  r_face  = r_edge(nr+1)
  mdot_edge(nr+1) = 6.0_dp * pi * sqrt(r_face) * dGdr

  select case (outer_bc_type)
  case (0)
  case (1)
     mdot_edge(nr+1) = max(0.0_dp, mdot_edge(nr+1))
  case (2)
     mdot_edge(nr+1) = 0.0_dp
  case default
     stop 'unknown outer_bc_type'
  end select

  ! -- 3. LSigma (same)
  do i = 1, nr
     LSigma(i) = ( mdot_edge(i) - mdot_edge(i+1) ) &
                 / ( 2.0_dp * pi * r(i) * dr_cell(i) )
  end do
end subroutine compute_LSigma

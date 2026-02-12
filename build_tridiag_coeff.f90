subroutine build_tridiag_coeff(nr, r, nu, dt, theta, Sigma, a, b, c)
  use kind_params, only : dp, i4b
  use constants,   only : pi
  use mod_global, only  : outer_bc_type, inner_bc_type, INNER_FIXED_MDOT
  implicit none
  integer(i4b), intent(in)  :: nr
  real(dp),    intent(in)   :: r(nr), nu(nr), Sigma(nr)
  real(dp),    intent(in)   :: dt, theta
  real(dp),    intent(out)  :: a(nr), b(nr), c(nr)

  integer(i4b) :: i
  real(dp) :: r_edge(nr+1), dr_cell(nr)
  real(dp) :: A_L, A_C, B_C, B_R
  real(dp) :: aL_i, bL_i, cL_i
  real(dp) :: A_cell, dr_im1, dr_i, r_face, coef

  !----------------------------------------------------
  ! 1. Geometry: cell edges r_edge and widths dr_cell
  !    (Option B: compute locally)
  !----------------------------------------------------
  r_edge(1) = r(1) - 0.5_dp*(r(2) - r(1))
  do i = 2, nr
     r_edge(i) = 0.5_dp*(r(i-1) + r(i))
  end do
  r_edge(nr+1) = r(nr) + 0.5_dp*(r(nr) - r(nr-1))

  do i = 1, nr
     dr_cell(i) = r_edge(i+1) - r_edge(i)
  end do

  !----------------------------------------------------
  ! 2. Loop over cells and build aL_i, bL_i, cL_i
  !
  !   For each cell i, we use:
  !      dSigma_i/dt = ( mdot_i - mdot_{i+1} ) / A_cell
  !   with
  !      mdot_i   = A_L * Sigma_{i-1} + A_C * Sigma_i
  !      mdot_i+1 = B_C * Sigma_i     + B_R * Sigma_{i+1}
  !
  !   => L[Σ]_i = aL_i Σ_{i-1} + bL_i Σ_i + cL_i Σ_{i+1}
  !      aL_i = A_L / A_cell
  !      bL_i = (A_C - B_C) / A_cell
  !      cL_i = - B_R / A_cell
  !
  !   where A_cell = 2π r_i Δr_i
  !----------------------------------------------------
  do i = 1, nr

     if (i == 1) then

        dr_im1 = r(1) - r_edge(1)
        r_face = r_edge(1)
        ! coef for inner face j=1 (used only in non-decretion mode)
        coef   = -6.0_dp*pi*sqrt(r_face)/dr_im1

        if (inner_bc_type == INNER_FIXED_MDOT) then
           A_L = 0.0_dp
           A_C = 0.0_dp
        else
           A_L = 0.0_dp
           A_C = coef * ( + nu(1)*sqrt(r(1)) )
        end if

        ! mdot_2 between r1 and r2 (same for both modes)
        dr_i   = r(2) - r(1)
        r_face = 0.5_dp*(r(1) + r(2))
        coef   = -6.0_dp*pi*sqrt(r_face)/dr_i

        B_C = coef * ( - nu(1)*sqrt(r(1)) )   ! Sigma(1)
        B_R = coef * ( + nu(2)*sqrt(r(2)) )   ! Sigma(2)

        A_cell = 2.0_dp*pi*r(1)*dr_cell(1)

        aL_i = 0.0_dp
        bL_i = (A_C - B_C) / A_cell
        cL_i = - B_R        / A_cell

     else if (i == nr) then

        ! Face j = nr (between r_{nr-1} & r_{nr})
        dr_im1 = r(nr) - r(nr-1)
        r_face = 0.5_dp*(r(nr-1) + r(nr))
        coef   = -6.0_dp*pi*sqrt(r_face)/dr_im1

        A_L = coef * (-nu(nr-1)*sqrt(r(nr-1)))   ! Sigma(nr-1)
        A_C = coef * ( nu(nr  )*sqrt(r(nr  )))   ! Sigma(nr)

        ! Outer face j = nr+1
        select case (outer_bc_type)
        case (0)
           ! open: same as in compute_LSigma:
           dr_i   = r_edge(nr+1) - r(nr)
           r_face = r_edge(nr+1)
           coef   = -6.0_dp*pi*sqrt(r_face)/dr_i

           ! mdot_{nr+1} = coef * (0 - nu(nr)*sqrt(r(nr))*Sigma(nr))
           B_C = coef * ( - nu(nr)*sqrt(r(nr)) )
        case default
           ! zero-flux or no-outflow: mdot_{nr+1} does not depend
           ! linearly on Sigma in our simple treatment here
           B_C = 0.0_dp
        end select

        B_R = 0.0_dp

        A_cell = 2.0_dp*pi*r(nr)*dr_cell(nr)

        aL_i = A_L / A_cell
        bL_i = (A_C - B_C) / A_cell
        cL_i = 0.0_dp

     else
        !---------------------------------------------
        ! Internal cell 2 <= i <= nr-1
        !   Uses mdot_edge(i)     (between r_{i-1}, r_i)
        !        mdot_edge(i + 1) (between r_i, r_{i+1})
        !
        ! mdot_i   = coef_i   * [ nu(i)*sqrt(ri)*Sigma(i)
        !                        -nu(i-1)*sqrt(ri-1)*Sigma(i-1)]
        !
        ! mdot_i+1 = coef_ip1 * [ nu(i+1)*sqrt(ri+1)*Sigma(i+1)
        !                        -nu(i)*sqrt(ri)*Sigma(i)      ]
        !---------------------------------------------
        ! mdot_i:
        dr_im1 = r(i) - r(i-1)
        r_face = 0.5_dp*(r(i-1) + r(i))
        coef   = -6.0_dp*pi*sqrt(r_face)/dr_im1

        A_L = coef * ( - nu(i-1)*sqrt(r(i-1)) )                 ! Sigma(i-1)
        A_C = coef * ( + nu(i  )*sqrt(r(i  )) )                 ! Sigma(i)

        ! mdot_{i+1}:
        dr_i   = r(i+1) - r(i)
        r_face = 0.5_dp*(r(i) + r(i+1))
        coef   = -6.0_dp*pi*sqrt(r_face)/dr_i

        B_C = coef * ( - nu(i  )*sqrt(r(i  )) )                 ! Sigma(i)
        B_R = coef * ( + nu(i+1)*sqrt(r(i+1)) )                 ! Sigma(i+1)

        A_cell = 2.0_dp*pi*r(i)*dr_cell(i)

        aL_i = A_L        / A_cell
        bL_i = (A_C-B_C)  / A_cell
        cL_i = - B_R      / A_cell

     end if

     !---------------------------------
     ! 3. Theta method (Crank–Nicolson)
     !    (I - theta*dt*L) Sigma^{n+1} = ...
     !---------------------------------
     a(i) = -theta*dt * aL_i
     b(i) =  1.0_dp  - theta*dt * bL_i
     c(i) = -theta*dt * cL_i

  end do

end subroutine build_tridiag_coeff

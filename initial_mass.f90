  module initial_mass_mod
  use kind_params, only : i4b, dp
  use mod_global, only : nr, r, sigmat
  implicit none

  contains
  subroutine initial_mass(M0)
     real(dp), intent(out) :: M0
     integer(i4b) :: i
     real(dp) :: r_edge(nr+1), dr_cell(nr)
     real(dp), parameter :: pi = 3.14159265358979323846_dp

     ! Construct edges exactly as in evolve.f90
     r_edge(1) = r(1) - 0.5_dp * (r(2) - r(1))
     do i = 2, nr
        r_edge(i) = 0.5_dp * (r(i-1) + r(i))
     end do
     r_edge(nr+1) = r(nr) + 0.5_dp * (r(nr) - r(nr-1))

     do i = 1, nr
        dr_cell(i) = r_edge(i+1) - r_edge(i)
     end do

     M0 = 0.0_dp
     do i = 1, nr
        M0 = M0 + sigmat(1,i) * 2.0_dp * pi * r(i) * dr_cell(i)
     end do
  end subroutine initial_mass

  end module initial_mass_mod

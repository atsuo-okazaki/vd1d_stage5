!==========================================================
! solve_tridiag
!
! Solves the tridiagonal system:
!   a(i)*x(i-1) + b(i)*x(i) + c(i)*x(i+1) = d(i)
!
! INPUT:
!   a(1..n)   -- lower diagonal (a(1) must be 0)
!   b(1..n)   -- main diagonal
!   c(1..n)   -- upper diagonal (c(n) must be 0)
!   d(1..n)   -- right-hand-side
!
! OUTPUT:
!   x(1..n)   -- solution vector
!==========================================================
subroutine solve_tridiag(n, a, b, c, d, x)
  use kind_params, only : dp, i4b
  implicit none

  integer(i4b), intent(in) :: n
  real(dp), intent(in)  :: a(n), b(n), c(n), d(n)
  real(dp), intent(out) :: x(n)

  real(dp) :: cp(n), dpv(n)
  integer(i4b) :: i
  real(dp) :: denom

  !-----------------------------------------
  ! Forward elimination
  !-----------------------------------------
  cp(1) = c(1) / b(1)
  dpv(1) = d(1) / b(1)

  do i = 2, n
     denom = b(i) - a(i)*cp(i-1)
     cp(i) = c(i) / denom
     dpv(i) = (d(i) - a(i)*dpv(i-1)) / denom
  end do

  !-----------------------------------------
  ! Back substitution
  !-----------------------------------------
  x(n) = dpv(n)
  do i = n-1, 1, -1
     x(i) = dpv(i) - cp(i)*x(i+1)
  end do

end subroutine solve_tridiag

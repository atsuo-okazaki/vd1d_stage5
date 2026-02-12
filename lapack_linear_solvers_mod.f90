module lapack_linear_solvers
  use kind_params, only: dp, i4b
  implicit none

  interface
     subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
       integer :: n, nrhs, lda, ldb, info
       integer :: ipiv(*)
       double precision :: a(lda,*), b(ldb,*)
     end subroutine dgesv
  end interface

contains

  subroutine solve_dgesv(n, A_in, b_in, x, info)
    implicit none
    integer(i4b), intent(in)  :: n
    real(dp),     intent(in)  :: A_in(n,n)
    real(dp),     intent(in)  :: b_in(n)
    real(dp),     intent(out) :: x(n)
    integer(i4b), intent(out) :: info

    real(dp) :: A(n,n), b(n,1)
    integer(i4b) :: ipiv(n)

    A = A_in
    b(:,1) = b_in

    call dgesv(n, 1, A, n, ipiv, b, n, info)

    if (info == 0) then
       x = b(:,1)
    else
       x = 0.0_dp
    end if
  end subroutine solve_dgesv

end module lapack_linear_solvers

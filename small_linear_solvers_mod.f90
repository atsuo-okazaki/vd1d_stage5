module small_linear_solvers
  use kind_params, only: dp, i4b
  implicit none
contains

  subroutine solve_small_linear_scaled(n, A_in, b_in, x, info)
    implicit none
    integer(i4b), intent(in)  :: n
    real(dp),     intent(in)  :: A_in(n,n)
    real(dp),     intent(in)  :: b_in(n)
    real(dp),     intent(out) :: x(n)
    integer(i4b), intent(out) :: info

    real(dp) :: A(n,n), b(n)
    real(dp) :: s(n)                 ! row scales
    integer(i4b) :: i, j, k, piv
    real(dp) :: eps, rmax, ratio, best, tmp, factor

    A = A_in
    b = b_in
    x = 0.0_dp
    info = 0
    eps = epsilon(1.0_dp)

    ! Build row scaling factors: s(i) = max_j |A(i,j)|
    do i = 1, n
       rmax = 0.0_dp
       do j = 1, n
          rmax = max(rmax, abs(A(i,j)))
       end do
       s(i) = rmax
    end do

    ! Forward elimination with scaled partial pivoting
    do k = 1, n-1

       piv  = k
       best = -1.0_dp
       do i = k, n
          if (s(i) > 0.0_dp) then
             ratio = abs(A(i,k)) / s(i)
          else
             ratio = 0.0_dp
          end if
          if (ratio > best) then
             best = ratio
             piv  = i
          end if
       end do

       ! If best is tiny, the column is (effectively) singular at this step
       if (best <= eps) then
          info = k
          return
       end if

       ! Swap rows k and piv (A, b, and s)
       if (piv /= k) then
          do j = k, n
             tmp = A(k,j); A(k,j) = A(piv,j); A(piv,j) = tmp
          end do
          tmp = b(k); b(k) = b(piv); b(piv) = tmp
          tmp = s(k); s(k) = s(piv); s(piv) = tmp
       end if

       ! Pivot sanity check (relative to its row scale)
       if (abs(A(k,k)) <= eps * max(1.0_dp, s(k))) then
          info = k
          return
       end if

       ! Eliminate below pivot
       do i = k+1, n
          factor = A(i,k) / A(k,k)
          A(i,k) = 0.0_dp
          do j = k+1, n
             A(i,j) = A(i,j) - factor * A(k,j)
          end do
          b(i) = b(i) - factor * b(k)
       end do

    end do

    ! Final pivot check
    if (abs(A(n,n)) <= eps * max(1.0_dp, s(n), 1.0_dp)) then
      info = n
      return
    end if

    ! Back substitution
    x(n) = b(n) / A(n,n)
    do i = n-1, 1, -1
       tmp = b(i)
       do j = i+1, n
          tmp = tmp - A(i,j) * x(j)
       end do
       if (abs(A(i,i)) <= eps * max(1.0_dp, s(i), 1.0_dp)) then
          info = i
          return
       end if
       x(i) = tmp / A(i,i)
    end do

  end subroutine solve_small_linear_scaled

end module small_linear_solvers

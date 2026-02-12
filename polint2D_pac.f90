!***********************
   MODULE polint2D_pac
!***********************
!  This module (=2D polynomial-interpolation package) collects
!  subroutines and functions used in 2D polynomial interpolation
!  procedure.

     use kind_params, only : dp, i4b, lgt
     implicit none

   CONTAINS
!     ********************************************
           SUBROUTINE        hunt(xx,n,x,jlo)
!     ********************************************
!       Given an array xx(1:n), and given a value x,
!       returns a value jlo such that x is between xx(jlo)
!       and xx(jlo+1). xx(1:n) must be monotonic, either
!       increasing or decreasing. jlo=0 or jlo=n is returned
!       to indicate that x is out of range. jlo on input is
!       taken as the guess for jlo on output. [A bug-fixed 
!       f90 version for the f77 version of 'Numerical
!       Recipes 2nd Ed.', p.112 by Press et al. (1992),
!       which returns a wrong value jlo=0 when x=xx(1).]

      implicit none

      integer, intent(in) ::  n
      real(DP), intent(in) :: x
      real(DP), dimension(n), intent(in) :: xx
      integer, intent(inout) :: jlo
      integer :: inc, jhi, jm
      logical :: ascnd

      !n = size(xx)
      ascnd = (xx(n) >= xx(1))
      if (jlo <= 0 .or. jlo > n) then
         jlo=0
         jhi=n+1
      else
         inc=1
         if (x >= xx(jlo) .eqv. ascnd) then
            do
               jhi=jlo+inc
               if (jhi > n) then
                  jhi=n+1
                  exit
               else
                  if (x < xx(jhi) .eqv. ascnd) exit
                  jlo=jhi
                  inc=inc+inc
               end if
            end do
         else
            jhi=jlo
            do
               jlo=jhi-inc
               if (jlo < 1) then
                  jlo=0
                  exit
               else
                  if (x >= xx(jlo) .eqv. ascnd) exit
                  jhi=jlo
                  inc=inc+inc
               end if
            end do
         end if
      end if
      do
         if (jhi-jlo <= 1) then
            if (x == xx(n)) jlo=n-1
            if (x == xx(1)) jlo=1
            exit
         else
            jm=(jhi+jlo)/2
            if (x >= xx(jm) .eqv. ascnd) then
               jlo=jm
            else
               jhi=jm
            end if
         end if
      end do

      END SUBROUTINE hunt

!     ************************************************************
           SUBROUTINE      polin2mod(x1a, x2a, ya, &
                                     jx1, jx2, x1, x2, y, dy)
!     ************************************************************
!       Given arrays x1a(1:m) and x2a(1:n) of independent
!       variables, and an m*n array of function values ya(1:m,1:n),
!       tabulated at the grid points x1a and x2a, and given values
!       x1 and x2 of the independent variables, this routine returns
!       a function value y interpolated with a polynomial of degree
!       mm-1, and an accuracy indication dy (based only on the 
!       interpolation in the x1 direction, however).

      implicit none

      integer (i4b), parameter :: mm = 3 ! mm-1: degree of the polynomial
      integer (i4b), intent(in) :: jx1, jx2
      integer (i4b) :: m, n, mm1, i, j, i1, i2, j1, j2
      real (dp), dimension(:), intent(in) :: x1a
      real (dp), dimension(:), intent(in) :: x2a
      real (dp), dimension(:,:), intent(in) :: ya
      real (dp), intent(in) :: x1, x2
      real (dp), intent(out) :: y, dy
      real (dp), dimension(mm) :: x1asub, f1tmp
      real (dp), dimension(mm) :: x2asub, f2tmp

      IF (size(x1a) == size(ya,1)) THEN
         m = size(x1a)
      ELSE
         WRITE (*,*) 'size(x1a)=', size(x1a), &
              ' /= size(ya,1)=', size(ya,1), ' in polin2mod'
         STOP
      END IF
      
      IF (size(x2a) == size(ya,2)) THEN
         n = size(x2a)
      ELSE
         WRITE (*,*) 'size(x2a)=', size(x2a), &
              ' /= size(ya,2)=', size(ya,2), ' in polin2mod'
         STOP
      END IF

      IF (mm > MIN(m, n)) THEN
         WRITE (6, "('###  Number of points used (', i3, &
              & ') is too large for a (', i3, ',', i3, ') array', &
              & ' ###')") mm, m, n
         STOP
      END IF

      mm1 = mm - 1
      i1 = MIN(MAX(jx1 - mm1/2, 1), m - mm1)
      i2 = i1 + mm1
      j1 = MIN(MAX(jx2 - mm1/2, 1), n - mm1)
      j2 = j1 + mm1

      DO i = i1, i2
         x1asub(i-i1+1) = x1a(i)
      END DO
      DO j = j1, j2
         x2asub(j-j1+1) = x2a(j)
      END DO

      DO j = j1, j2
         DO i = i1, i2
            f1tmp(i-i1+1) = ya(i,j)
         END DO
         CALL polint(x1asub, f1tmp, x1, f2tmp(j-j1+1), dy)
      END DO
      CALL polint(x2asub, f2tmp, x2, y, dy)

      END SUBROUTINE polin2mod

!     **********************************************
          SUBROUTINE    polint(xa, ya, x, y, dy)
!     **********************************************
!     Given arrays xa and ya of length N, and given a value x,
!     this routine returns a value y, and an error estimate dy.
!     If P(x) is the polynomial of degree N − 1 such that
!     P(xai) = yai, i = 1, . . . ,N, then the returned value y = P(x).

      implicit none

      real (dp), dimension(:), intent(in) :: xa, ya
      real (dp), intent(in) :: x
      real (dp), intent(out) :: y, dy
      real (dp), dimension(size(xa)) :: c, d, den, ho
      integer (i4b), dimension(1) :: imin
      integer (i4b) :: m, n, ns

      IF (size(xa) == size(ya)) THEN
         n = size(xa)
      ELSE
         WRITE (*,*) 'size(xa)=',size(xa),' /= size(ya)=',size(ya), &
              ' in polint'
         STOP
      END IF

      !-- Initialize the tableau of c’s and d’s.
      c = ya
      d = ya
      ho = xa - x

      !-- Find index ns of closest table entry.
      imin = MINLOC(ABS(x-xa))
      ns = imin(1)

      y = ya(ns)  ! This is the initial approximation to y.
      ns = ns-1

      !-- For each column of the tableau,
      !   we loop over the current c’s and d’s and update them.
      DO m = 1,n-1
         den(1:n-m) = ho(1:n-m) - ho(1+m:n)
         !-- The followin error occurs only if two input xa’s are
         !   (to within roundoff) identical.
         IF (any(den(1:n-m) == 0.0)) THEN
            WRITE (*,*) 'polint: calculation failure'
            STOP
         END IF
         den(1:n-m) = (c(2:n-m+1) - d(1:n-m)) / den(1:n-m)

         !-- Update c's and d's.
         d(1:n-m) = ho(1+m:n) * den(1:n-m)
         c(1:n-m) = ho(1:n-m) * den(1:n-m)

         !-- After each column in the tableau is completed, we decide
         ! which correction, c or d, we want to add to our accumulating
         ! value of y, i.e., which path to take through
         ! the tableau—forking up or down. We do this in such a
         ! way as to take the most “straight line” route through the
         ! tableau to its apex, updating ns accordingly to keep track
         ! of where we are. This route keeps the partial approximations
         ! centered (insofar as possible) on the target x. The
         ! last dy added is thus the error indication.
         IF (2*ns < n-m) THEN
            dy = c(ns+1)
         ELSE
            dy = d(ns)
            ns = ns - 1
         END IF
         y = y + dy
      END DO

    END SUBROUTINE polint

 END MODULE polint2D_pac

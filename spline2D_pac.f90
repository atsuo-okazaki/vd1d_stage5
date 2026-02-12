!  *************************
       MODULE spline2D_pac
!  *************************
!  This module (=2D spline-interpolation package) collects
!  subroutines and functions used in 2D spline interpolation
!  procedure.

     use kind_params, only : dp, i4b, lgt
     implicit none

    CONTAINS
!     *********************************************
          SUBROUTINE splie2 (x1a, x2a, ya, y2a)
!     *********************************************
!     Given an MxN tabulated function ya, and N tabulated
!     independent variables x2a, this routine constructs
!     one-dimensional natural cubic splines of the rows of
!     ya and returns the second derivatives in the MxN array y2a.
!     (The array x1a is included in the argument list merely
!     for consistency with routine splin2.)

      implicit none
 
      real (dp), dimension(:), intent(in) :: x1a, x2a
      real (dp), dimension(:,:), intent(in) :: ya
      real (dp), dimension(:,:), intent(out) :: y2a
      integer (i4b)  :: j, m, ndum

      IF (size(x1a) == size(ya,1) .AND. &
          size(ya,1) == size(y2a,1)) THEN
         m = size(x1a)
      ELSE
         WRITE (*,*) 'Error in splie2: First array sizes don''t coincide'
         STOP
      END IF

      IF (size(x2a) == size(ya,2) .AND. &
          size(ya,2) == size(y2a,2)) THEN
         ndum = size(x2a)
      ELSE
         WRITE (*,*) 'Error in splie2: Second array sizes don''t coincide'
         STOP
      END IF

      DO j = 1, m
         CALL spline(x2a, ya(j,:), 1.0e30_dp, 1.0e30_dp, y2a(j,:))
      END DO

      END SUBROUTINE splie2

!     **************************************************
          FUNCTION splin2(x1a, x2a, ya, y2a, x1, x2)
!     **************************************************
!     Given x1a, x2a, ya as described in splie2 and y2a as
!     produced by that routine; and given a desired interpolating
!     point x1, x2; this routine returns an interpolated function
!     value by bicubic spline interpolation.

      implicit none

      real (dp), dimension(:), intent(in) :: x1a, x2a
      real (dp), dimension(:,:), intent(in) :: ya, y2a
      real (dp), intent(in) :: x1, x2
      real (dp) :: splin2
      integer (i4b) :: j, m, ndum
      real (dp), dimension(size(x1a)) :: yytmp, y2tmp2

      !m=assert_eq(size(x1a),size(ya,1),size(y2a,1),'splin2: m')
      !ndum=assert_eq(size(x2a),size(ya,2),size(y2a,2),'splin2: ndum')
      IF (size(x1a) == size(ya,1) .AND. &
          size(ya,1) == size(y2a,1)) THEN
         m = size(x1a)
      ELSE
         WRITE (*,*) 'Error in splin2: First array sizes don''t coincide'
         STOP
      END IF

      IF (size(x2a) == size(ya,2) .AND. &
          size(ya,2) == size(y2a,2)) THEN
         ndum = size(x2a)
      ELSE
         WRITE (*,*) 'Error in splin2: Second array sizes don''t coincide'
         STOP
      END IF

      DO j=1,m
         yytmp(j) = splint(x2a, ya(j,:), y2a(j,:), x2)
      END DO

      CALL spline(x1a, yytmp, 1.0e30_dp, 1.0e30_dp, y2tmp2)

      splin2 = splint(x1a, yytmp, y2tmp2, x1)

      END FUNCTION splin2

!     **************************************************
          SUBROUTINE    spline(x, y, yp1, ypn, y2)
!     **************************************************
!     Given arrays x(1:n) and y(1:n) containing a tabulated
!     function, i.e., yi = f(xi), with x1 < x2 < ::: < xN,
!     and given values yp1 and ypn for the first derivative
!     of the interpolating function at points 1 and n, respectively,
!     this routine returns an array y2(1:n) of length n which
!     contains the second derivatives of the interpolating function
!     at the tabulated points xi. If yp1 and/or ypn are equal to
!     1.0e30 or larger, the routine is signaled to set the
!     corresponding boundary condition for a natural spline,
!     with zero second derivative on that boundary.

      implicit none

      real(dp), dimension(:), intent(in) :: x, y
      real(dp), intent(in) :: yp1, ypn
      real(dp), dimension(:), intent(out) :: y2
      integer(i4b) :: n
      real(dp), dimension(size(x)) :: a, b, c, r

      !n = assert_eq(size(x),size(y),size(y2),'spline')
      IF (size(x) == size(y) .AND. size(x) == size(y2)) THEN
         n = size(x)
      ELSE
         WRITE (*,*) 'size(x)=', size(x), &
              ' /= size(y)=', size(y), &
              ' or size(y2)=', size(y2), ' in spline'
         STOP
      END IF

      c(1:n-1) = x(2:n)-x(1:n-1)
      r(1:n-1) = 6.0_dp*((y(2:n)-y(1:n-1))/c(1:n-1))
      r(2:n-1) = r(2:n-1)-r(1:n-2)
      a(2:n-1) = c(1:n-2)
      b(2:n-1) = 2.0_dp*(c(2:n-1)+a(2:n-1))
      b(1) = 1.0
      b(n) = 1.0
      if (yp1 > 0.99e30_dp) then
         r(1) = 0.0
         c(1) = 0.0
      else
         r(1) = (3.0_dp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
         c(1) = 0.5
      end if
      if (ypn > 0.99e30_dp) then
         r(n) = 0.0
         a(n) = 0.0
      else
         r(n) = (-3.0_dp/(x(n)-x(n-1)))*((y(n)-y(n-1)) &
                /(x(n)-x(n-1))-ypn)
         a(n) = 0.5
      end if
      call tridag(a(2:n), b(1:n), c(1:n-1), r(1:n), y2(1:n))

      END SUBROUTINE spline

!     ***************************************
          FUNCTION splint(xa, ya, y2a, x)
!     ***************************************
!     Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate
!     a function (with the xai's in order), and given the array
!     y2a(1:n), which is the output from spline above, and given
!     a value of x, this routine returns a cubic-spline interpolated
!     value y.            

      implicit none

      real(dp), dimension(:), intent(in) :: xa, ya, y2a
      real(dp), intent(in) :: x
      real(dp) :: splint
      integer(i4b) :: khi, klo, n
      real(dp) :: a, b, h

      !n = assert_eq(size(xa),size(ya),size(y2a),'splint')
      IF (size(xa) == size(ya) .AND. size(xa) == size(y2a)) THEN
         n = size(xa)
      ELSE
         WRITE (*,*) 'size(xa)=', size(xa), &
              ' /= size(ya)=', size(ya), &
              ' or size(y2a)=', size(y2a), ' in splint'
         STOP
      END IF

      klo = max(min(locate(xa,x),n-1),1)
      khi = klo+1
      h = xa(khi)-xa(klo)
      IF (h == 0.0) THEN
         print *, 'bad xa input in splint'
         STOP
      END IF
      a = (xa(khi)-x)/h
      b = (x-xa(klo))/h
      splint = a*ya(klo) + b*ya(khi) &
           + ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi)) &
           * (h**2)/6.0_dp

      END FUNCTION splint

!     *******************************************
          SUBROUTINE    tridag(a, b, c, r, u)
!     *******************************************
      implicit none

      real(dp), dimension(:), intent(in) :: a,b,c,r
      real(dp), dimension(:), intent(out) :: u
      real(dp), dimension(size(b)) :: gam
      integer(i4b) :: n,j
      real(dp) :: bet

      !n = assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/), &
      !     'tridag')
      IF (size(a)+1 == size(b) .AND. size(a) == size(c) &
           .AND. size(b) == size(r) .AND. size(b) == size(u)) THEN
         n = size(a) + 1
      ELSE
         WRITE (*,*) 'Array sized don''t coincide'
         STOP
      END IF
      
      bet = b(1)
      if (bet == 0.0) then
         print *, 'tridag_ser: Error at code stage 1'
         stop
      end if

      u(1) = r(1) / bet
      do j = 2, n
         gam(j) = c(j-1) / bet
         bet = b(j) - a(j-1) * gam(j)
         if (bet == 0.0) then
            print *, 'tridag_ser: Error at code stage 2'
            stop
         end if
         u(j) = (r(j) - a(j-1) * u(j-1)) / bet
      end do

      do j = n-1, 1, -1
         u(j) = u(j) - gam(j+1) * u(j+1)
      end do

      END SUBROUTINE tridag
    
!     *********************************
          FUNCTION    locate(xx, x)
!     *********************************
      implicit none

      real(dp), dimension(:), intent(in) :: xx
      real(dp), intent(in) :: x
      integer(i4b) :: locate
      integer(i4b) :: n,jl,jm,ju
      logical(lgt) :: ascnd

      n = size(xx)
      !ascnd = (xx(n) >= xx(1))
      IF (xx(n) >= xx(1)) THEN
         ascnd = .true.
      ELSE
         ascnd = .false.
      END IF
        
      jl = 0
      ju = n+1
      do
         if (ju-jl <= 1) exit
         jm = (ju+jl)/2
         if (ascnd .eqv. (x >= xx(jm))) then
            jl = jm
         else
            ju = jm
         end if
      end do
      if (x == xx(1)) then
         locate = 1
      else if (x == xx(n)) then
         locate = n-1
      else
         locate = jl
      end if

      END FUNCTION locate

  END MODULE spline2D_pac

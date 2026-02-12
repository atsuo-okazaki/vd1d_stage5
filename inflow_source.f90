module inflow_source_mod
  use kind_params,  only : dp, i4b
  use constants,    only : pi
  use mod_global,   only : use_inflow, t_in_start, t_in_end, &
                           rinj_min, rinj_max, mdot_inj_edd
  use units_disk_mod, only : mdot_nd_from_edd
  implicit none
contains

  !------------------------------------------
  ! Dimensionless inflow rate (PDE unit)
  !------------------------------------------
  pure function mdot_inflow_nd(t_nd) result(mdot_nd)
    real(dp), intent(in) :: t_nd
    real(dp)             :: mdot_nd

    if (.not. use_inflow) then
      mdot_nd = 0.0_dp
    else if (t_nd < t_in_start .or. t_nd > t_in_end) then
      mdot_nd = 0.0_dp
    else
      ! Convert on the fly: input is mdot_inj_edd (Edd unit)
      mdot_nd = mdot_nd_from_edd(mdot_inj_edd)
    end if
  end function mdot_inflow_nd


  !------------------------------------------
  ! Dimensionless source term:
  !   source(i) = dSigma'(i)/dt'
  !
  ! Normalization:
  !   mdot_inj_nd = ∫ 2π r' source(r') dr'
  !------------------------------------------
  subroutine compute_source(t_nd, r, dr, source)
    real(dp), intent(in)  :: t_nd
    real(dp), intent(in)  :: r(:)
    real(dp), intent(in)  :: dr(:)
    real(dp), intent(out) :: source(:)

    integer(i4b) :: i, n
    real(dp)     :: mdot_nd
    real(dp)     :: area_tot_nd
    real(dp)     :: r_left, r_right
    real(dp)     :: overlap, area_i_nd, weight_i
    real(dp)     :: coeff

    n      = size(r)
    source = 0.0_dp

    mdot_nd = mdot_inflow_nd(t_nd)
    if (mdot_nd == 0.0_dp) return

    ! 1) total "injection area" in dimensionless sense
    area_tot_nd = 0.0_dp
!$omp parallel do default(shared) private(i,    &
!$omp&  r_left, r_right, overlap, area_i_nd)    &
!$omp&  reduction(+: area_tot_nd)
    do i = 1, n
       r_left  = r(i) - 0.5_dp*dr(i)
       r_right = r(i) + 0.5_dp*dr(i)

       overlap = max(0.0_dp, min(r_right, rinj_max) - max(r_left, rinj_min))
       if (overlap > 0.0_dp) then
          area_i_nd   = 2.0_dp * pi * r(i) * overlap
          area_tot_nd = area_tot_nd + area_i_nd
       end if
    end do
!$omp end parallel do

    if (area_tot_nd <= 0.0_dp) return

    ! 2) Prefactor so that integral gives mdot_nd:
    coeff = mdot_nd

!$omp parallel do default(shared) private(i,  &
!$omp& r_left, r_right, overlap, area_i_nd, weight_i)
    do i = 1, n
       r_left  = r(i) - 0.5_dp*dr(i)
       r_right = r(i) + 0.5_dp*dr(i)

       overlap = max(0.0_dp, min(r_right, rinj_max) - max(r_left, rinj_min))
       if (overlap > 0.0_dp .and. r(i) > 0.0_dp .and. dr(i) > 0.0_dp) then
          area_i_nd = 2.0_dp * pi * r(i) * overlap
          weight_i  = area_i_nd / area_tot_nd
          source(i) = coeff * weight_i / (2.0_dp*pi*r(i)*dr(i))
       else
          source(i) = 0.0_dp
       end if
    end do
!$omp end parallel do
    !write (*, '("t_nd =", 1pe12.4, ": source(n-3), ..., source(n) =", 1p4e12.4)') &
    !      t_nd, (source(i), i=n-3, n)

  end subroutine compute_source

end module inflow_source_mod

!==============================================================
! wind_mod.f90
!
! Wind model for Be disk truncation.
! - Wind parameters are taken from mod_global.
! - wind_velocity(:) and wind_density(:) are precomputed arrays.
!==============================================================
module wind_mod
  use kind_params, only : dp, i4b
  use constants,   only : pi, msun, year
  use mod_global,  only : use_wind_truncation, mdot_w_msunyr, &
                          vinf_w, beta_w, f_rho_wind, R_star
  implicit none

  logical :: initialized = .false.

  real(dp), allocatable :: rho_wind(:)  !! wind density profile [g/cm^3]
  real(dp), allocatable :: vw_wind(:)   !! wind velocity profile [cm/s]

  real(dp), parameter :: vw_floor = 1.0e6_dp  !! minimum wind speed [cm/s]

contains

  !================================================================
  subroutine init_wind_model(nr, r_cgs)
    !! Precompute wind velocity and density over the radial grid.
    integer(i4b), intent(in) :: nr
    real(dp),     intent(in) :: r_cgs(nr)

    real(dp) :: mdot_cgs, x
    integer(i4b), save :: iu_w = -1
    integer(i4b) :: i

    if (allocated(rho_wind)) deallocate(rho_wind)
    if (allocated(vw_wind))  deallocate(vw_wind)

    allocate(rho_wind(nr))
    allocate(vw_wind(nr))

    ! Convert wind mass-loss rate from Msun/yr to g/s
    mdot_cgs = mdot_w_msunyr * msun / year

    do i = 1, nr
       if (r_cgs(i) <= 0.0_dp) then
          vw_wind(i)  = vw_floor
          rho_wind(i) = 0.0_dp
       else
          !---------------------------------------------------
          ! 1) Wind velocity (beta-law)
          !    v(r) = v_inf * (1 - R_star / r)^beta
          !---------------------------------------------------
          x = 1.0_dp - R_star / r_cgs(i)
          x = max(0.0_dp, min(1.0_dp, x))
          vw_wind(i) = max( vinf_w * x**beta_w , vw_floor )

          !---------------------------------------------------
          ! 2) Wind density from mass continuity
          !    rho_w = Mdot_w / (4π r^2 v_w)
          !---------------------------------------------------
          rho_wind(i) = mdot_cgs / (4.0_dp * pi * r_cgs(i)**2 * vw_wind(i))
       end if
       !print "('r (cm) =', 1pe13.5, ': rho_wind (g/cm^3) =', 1pe13.5)", &
       !     r_cgs(i), rho_wind(i)
    end do

    open(newunit=iu_w, file='wind_density.dat', form='formatted')
    write (iu_w, "(a, 4x, a, 6x, a, 3x, a, 1x, a)") '#', 'r(cm)', &
       'r/R_star', 'vw_wind(cm/s)', 'rho_wind(g/cm^3)'
    do i = 1, nr
       write (iu_w, "(1p, 4e13.5)") r_cgs(i), r_cgs(i)/R_star, vw_wind(i), &
             rho_wind(i)
    end do
    close (iu_w)

    initialized = .true.
  end subroutine init_wind_model
  !================================================================

end module wind_mod

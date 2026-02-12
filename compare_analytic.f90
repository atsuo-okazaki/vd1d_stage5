module compare_analytic_mod
  use kind_params,    only : dp, i4b
  use mod_global,     only : nr, r, sigmat, dt, nu0_nd, r0
  use analytic_lbp_mod, only : sigma_lbp_ring
  implicit none
contains

  subroutine write_profile_with_analytic(it, M0)
    integer(i4b), intent(in) :: it
    real(dp),     intent(in) :: M0

    integer(i4b) :: i
    real(dp)     :: t, sigma_num, sigma_ana
    character(len=64) :: fname
    integer(i4b) :: unit

    t = dt * real(it-1, dp)

    write(fname,'("profile_t",i7.7,".dat")') it
    open(newunit=unit, file=fname, status='replace')

    do i = 1, nr
       sigma_num = sigmat(it, i)
       sigma_ana = sigma_lbp_ring(r(i), t, r0, nu0_nd, M0)
       write(unit,'(1p,3e15.7)') r(i), sigma_num, sigma_ana
    end do

    close(unit)
  end subroutine write_profile_with_analytic

end module compare_analytic_mod

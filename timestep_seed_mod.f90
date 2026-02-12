module timestep_seed_mod
  use kind_params, only : dp
  implicit none
  private
  public :: reset_dt_seed, get_dt_seed, set_dt_seed

  real(dp), save :: dt_seed_nd = 0.0_dp

contains

  subroutine reset_dt_seed(dt0)
    real(dp), intent(in) :: dt0
    dt_seed_nd = max(dt0, 0.0_dp)
  end subroutine reset_dt_seed

  real(dp) function get_dt_seed()
    get_dt_seed = dt_seed_nd
  end function get_dt_seed

  subroutine set_dt_seed(dt_used)
    real(dp), intent(in) :: dt_used
    if (dt_used > 0.0_dp) dt_seed_nd = dt_used
  end subroutine set_dt_seed

end module timestep_seed_mod

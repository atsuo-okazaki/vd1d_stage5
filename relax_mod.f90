module relax_mod
  use kind_params, only: dp
  implicit none
contains
  pure function ur_update(x_old, x_trial, omega) result(x_new)
    ! Under-relaxation: x_new = x_old + omega*(x_trial - x_old)
    real(dp), intent(in) :: x_old, x_trial, omega
    real(dp) :: x_new

    x_new = x_old + omega * (x_trial - x_old)
    
  end function ur_update
end module relax_mod

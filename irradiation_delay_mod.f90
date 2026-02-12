!===============================================================
! irradiation_delay_mod.f90
!
! Arbitrary time delay for inner-edge mass accretion rate Mdot_in.
! Uses a ring buffer of (t, mdot) and linear interpolation.
!
! Time unit:
!   - This module expects dimensionless time t_nd (t/t0).
!   - tau_lag is also given in dimensionless units.
!===============================================================
module irradiation_delay_mod
  !! Ring buffer for delayed irradiation driven by inner mdot history.
  !!
  !! Stores (t_nd, mdot_nd) and returns mdot_nd(t - tau_lag).
  !!
  !! - tau_lag_nd may vary with timestep.
  !! - If requested history is older than buffer, the oldest value is returned.

  use kind_params, only : dp, i4b
  implicit none

  type :: mdot_delay_buffer
     integer(i4b) :: nbuf = 0
     integer(i4b) :: head = 0
     integer(i4b) :: count = 0
     real(dp), allocatable :: t_nd(:)
     real(dp), allocatable :: mdot_nd(:)
  end type mdot_delay_buffer

  public :: mdot_delay_init, mdot_delay_push, mdot_delay_get

contains

  subroutine mdot_delay_init(buf, nbuf_in)
    type(mdot_delay_buffer), intent(inout) :: buf
    integer(i4b), intent(in) :: nbuf_in

    buf%nbuf  = nbuf_in
    buf%head  = 0
    buf%count = 0

    if (allocated(buf%t_nd))    deallocate(buf%t_nd)
    if (allocated(buf%mdot_nd)) deallocate(buf%mdot_nd)

    allocate(buf%t_nd(nbuf_in))
    allocate(buf%mdot_nd(nbuf_in))

    buf%t_nd(:)    = 0.0_dp
    buf%mdot_nd(:) = 0.0_dp
  end subroutine mdot_delay_init

  subroutine mdot_delay_push(buf, t_now_nd, mdot_now_nd)
    type(mdot_delay_buffer), intent(inout) :: buf
    real(dp), intent(in) :: t_now_nd, mdot_now_nd

    buf%head = mod(buf%head, buf%nbuf) + 1
    buf%t_nd(buf%head)    = t_now_nd
    buf%mdot_nd(buf%head) = mdot_now_nd

    if (buf%count < buf%nbuf) buf%count = buf%count + 1
  end subroutine mdot_delay_push

  function mdot_delay_get(buf, t_now_nd, tau_lag_nd) result(mdot_out)
    type(mdot_delay_buffer), intent(in) :: buf
    real(dp), intent(in) :: t_now_nd, tau_lag_nd
    real(dp) :: mdot_out

    integer(i4b) :: i, idx
    real(dp) :: t_target, dt_best, dt

    mdot_out = 0.0_dp
    if (buf%count == 0) return

    t_target = t_now_nd - max(tau_lag_nd, 0.0_dp)

    ! Search backward from newest to oldest
    dt_best = huge(1.0_dp)
    do i = 0, buf%count-1
       idx = buf%head - i
       if (idx <= 0) idx = idx + buf%nbuf

       dt = abs(buf%t_nd(idx) - t_target)
       if (dt < dt_best) then
          dt_best = dt
          mdot_out = buf%mdot_nd(idx)
       end if
    end do
  end function mdot_delay_get

end module irradiation_delay_mod

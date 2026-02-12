  module run_control_mod
    use kind_params, only : dp, i4b
    implicit none

    ! Defaults
    logical      :: do_restart = .false.
    integer(i4b) :: it_restart = 1
    character(len=256) :: chk_file = 'checkpoint.bin'
    integer(i4b) :: chkfreq    = 0
    integer(i4b) :: outfreq    = 0

  end module run_control_mod
  
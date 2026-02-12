module checkpoint_util_mod
  use kind_params, only : i4b
  implicit none
contains
  function checkpoint_name(base, it) result(fname)
    character(len=*), intent(in) :: base
    integer(i4b),     intent(in) :: it
    character(len=256) :: fname
    character(len=32)  :: s

    write(s,'("it", i9.9)') it
    fname = trim(base)//'_'//trim(adjustl(s))
  end function checkpoint_name
end module checkpoint_util_mod

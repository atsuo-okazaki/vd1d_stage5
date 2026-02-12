!==============================================================
! kind_params.f90: Define kind parameters
!==============================================================
 MODULE kind_params
    implicit none
    save

    integer, parameter :: sp = selected_real_kind(p=6,r=37)
    integer, parameter :: dp = selected_real_kind(p=13,r=200)
    integer, parameter :: i1b = selected_int_kind(r=2)
    integer, parameter :: i2b = selected_int_kind(r=4)
    integer, parameter :: i4b = selected_int_kind(r=9)
    integer, parameter :: i8b = selected_int_kind(r=18)
    integer, parameter :: lgt = kind(.true.)

 END MODULE kind_params

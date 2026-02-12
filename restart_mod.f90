module restart_mod
  use kind_params, only : dp, i4b
  use mod_global,  only : nr, nt, r, sigmat, nu_conv, Tmid, H, rho, kappaR, tauR, &
                          Qirr, dYdXi, k_iter, i_edge, r_edge, t_nd
  implicit none
  private
  public :: checkpoint_append, restart_read_to_it

  integer(i4b), parameter :: MAGIC = 31415926_i4b
  integer(i4b), parameter :: VERSION = 1_i4b

contains

  subroutine checkpoint_append(fname, it_save, t_save)
    ! Append one checkpoint record to a single stream/unformatted file.
    character(len=*), intent(in) :: fname
    integer(i4b),     intent(in) :: it_save
    real(dp),         intent(in) :: t_save

    integer(i4b), save :: iu_ap = -1
    integer(i4b) :: ios
    logical :: ex

    inquire(file=fname, exist=ex)

    if (.not. ex) then
       open(newfile=iu_ap, file=fname, access='stream', form='unformatted', status='new', &
            action='write', iostat=ios)
       if (ios /= 0) stop 'checkpoint_append: cannot open file'
    end if

    ! Header
    write(iu_ap) MAGIC
    write(iu_ap) VERSION
    write(iu_ap) nr
    write(iu_ap) nt
    write(iu_ap) it_save
    write(iu_ap) t_save

    ! Grid (optional but recommended for consistency checks)
    write(iu_ap) r(:)

    ! State at it_save
    write(iu_ap) sigmat(it_save, :)
    write(iu_ap) nu_conv(it_save, :)
    write(iu_ap) Tmid(it_save, :)
    write(iu_ap) H(it_save, :)
    write(iu_ap) rho(it_save, :)
    write(iu_ap) kappaR(it_save, :)
    write(iu_ap) tauR(it_save, :)
    write(iu_ap) Qirr(it_save, :)
    write(iu_ap) dYdXi(it_save, :)

    ! Diagnostics / metadata arrays (scalar per it)
    write(iu_ap) k_iter(it_save)
    write(iu_ap) i_edge(it_save)
    write(iu_ap) r_edge(it_save)
  end subroutine checkpoint_append


  subroutine restart_read_to_it(fname, it_target, it_loaded, t_loaded)
    ! Read checkpoints sequentially and load the record with it_save == it_target.
    ! If not found, this subroutine stops with an error.
    character(len=*), intent(in)  :: fname
    integer(i4b),     intent(in)  :: it_target
    integer(i4b),     intent(out) :: it_loaded
    real(dp),         intent(out) :: t_loaded

    integer(i4b) :: iu_re = -1
    integer(i4b) :: ios, magic_in, ver_in, nr_in, nt_in
    integer(i4b) :: it_save
    real(dp)     :: t_save
    real(dp), allocatable :: r_in(:)
    real(dp), allocatable :: buf(:)

    open(newunit=iu_re, file=fname, access='stream', form='unformatted', status='old', &
         action='read', iostat=ios)
    if (ios /= 0) stop 'restart_read_to_it: cannot open file'

    allocate(r_in(nr))
    allocate(buf(nr))

    it_loaded = -1
    t_loaded  = 0.0_dp

    do
       ! Try to read MAGIC; if EOF -> exit
       read(iu_re, iostat=ios) magic_in
       if (ios /= 0) exit

       read(iu_re) ver_in
       read(iu_re) nr_in
       read(iu_re) nt_in
       read(iu_re) it_save
       read(iu_re) t_save

       if (magic_in /= MAGIC) stop 'restart_read_to_it: bad MAGIC'
       if (ver_in   /= VERSION) stop 'restart_read_to_it: unsupported VERSION'
       if (nr_in    /= nr) stop 'restart_read_to_it: nr mismatch'
       if (nt_in    /= nt) then
          ! nt mismatch is not always fatal, but it typically indicates a different run setup.
          !stop 'restart_read_to_it: nt mismatch'
          write (*,*) 'restart_read_to_it: nt mismatch'
       end if

       read(iu_re) r_in(:)

       ! Read state
       read(iu_re) buf(:); if (it_save == it_target) sigmat(it_save,:) = buf(:)
       read(iu_re) buf(:); if (it_save == it_target) nu_conv(it_save,:) = buf(:)
       read(iu_re) buf(:); if (it_save == it_target) Tmid(it_save,:) = buf(:)
       read(iu_re) buf(:); if (it_save == it_target) H(it_save,:) = buf(:)
       read(iu_re) buf(:); if (it_save == it_target) rho(it_save,:) = buf(:)
       read(iu_re) buf(:); if (it_save == it_target) kappaR(it_save,:) = buf(:)
       read(iu_re) buf(:); if (it_save == it_target) tauR(it_save,:) = buf(:)
       read(iu_re) buf(:); if (it_save == it_target) Qirr(it_save,:) = buf(:)
       read(iu_re) buf(:); if (it_save == it_target) dYdXi(it_save,:) = buf(:)

       if (it_save == it_target) then
          read(iu_re) k_iter(it_save)
          read(iu_re) i_edge(it_save)
          read(iu_re) r_edge(it_save)

          it_loaded = it_save
          t_loaded  = t_save
          exit
       else
          ! Skip the scalars for this record
          call skip_int(iu_re)
          call skip_int(iu_re)
          call skip_real(iu_re)
       end if
    end do

    close(iu_re)

    if (it_loaded < 0) stop 'restart_read_to_it: requested it not found'

    ! Restore global time and ensure main can continue
    t_nd = t_loaded

    deallocate(r_in, buf)
  contains
    subroutine skip_int(iu)
      integer(i4b), intent(in) :: iu
      integer(i4b) :: tmp
      read(iu) tmp
    end subroutine skip_int
    subroutine skip_real(iu)
      integer(i4b), intent(in) :: iu
      real(dp) :: tmp
      read(iu) tmp
    end subroutine skip_real
  end subroutine restart_read_to_it

end module restart_mod

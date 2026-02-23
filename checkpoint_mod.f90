!==============================================================
! checkpoint_mod.f90
!   Binary checkpoint I/O for 1D viscous disk code (2-slice).
!
!   - Stores two time slices: (it-1) and (it) for theta=0.5 safety
!   - Uses unformatted stream I/O (portable record layout)
!   - Stores minimal run metadata + r-grid + state arrays
!
!   Notes:
!     * Call write_checkpoint2(it) AFTER sigmat(it,:) etc. are finalized.
!     * Call read_checkpoint2(it_restart) BEFORE entering the time loop;
!       it will fill (it_restart-1) and (it_restart) in global arrays,
!       and set t_nd = t_nd_chk (stored in file).
!
!   Author: (template)
!==============================================================
module checkpoint_mod
  use kind_params, only : dp, i4b
  use mod_global,  only : nr, r, sigmat, nu_conv, t_nd, dt, alpha, &
                          Tmid, H, rho, kappaR, kappa_planck, tauR, Qvis, Qrad, Qirr, dYdXi, is_shadow, &
                          use_irradiation_delay
  use irradiation_mod, only : save_irradiation_buffer, load_irradiation_buffer
  implicit none

  integer(i4b), parameter :: CKPT_VER = 5
  integer(i4b), parameter :: CKPT_MAGIC = int(z'434B5054', i4b)  ! 'CKPT'

contains

  subroutine write_checkpoint2(fname, it_chk)
    !! Write a single checkpoint by replacing the file (atomic via temp+rename).
    character(len=*), intent(in) :: fname
    integer(i4b),     intent(in) :: it_chk

    character(len=512) :: tmp
    integer(i4b), save :: iu_cp = -1
    integer(i4b) :: ios
    logical :: ex

    tmp = trim(fname)//'.tmp'

    open(newunit=iu_cp, file=tmp, access='stream', form='unformatted', status='replace', &
         action='write', iostat=ios)
    if (ios /= 0) then
       write(*,*) 'ERROR: cannot open checkpoint tmp file: ', trim(tmp)
       stop
    end if

    ! --- Header
    write(iu_cp) CKPT_MAGIC
    write(iu_cp) CKPT_VER
    write(iu_cp) nr
    write(iu_cp) it_chk
    write(iu_cp) t_nd
    write(iu_cp) dt
    write(iu_cp) alpha
    write(iu_cp) r(:)

    ! --- 2-slice fields needed for theta=0.5 restart
    ! Sigma_{n-1}, Sigma_n
    write(iu_cp) sigmat(it_chk-1, :)
    write(iu_cp) sigmat(it_chk,   :)

    ! nu_{n-1}, nu_n  (store nu_conv as the canonical time-slice)
    write(iu_cp) nu_conv(it_chk-1, :)
    write(iu_cp) nu_conv(it_chk,   :)

    ! --- Optional but recommended: structure 2-slice (versioned)
    write(iu_cp) Tmid(it_chk-1, :)
    write(iu_cp) Tmid(it_chk,   :)
    write(iu_cp) H(it_chk-1, :)
    write(iu_cp) H(it_chk,   :)
    write(iu_cp) rho(it_chk-1, :)
    write(iu_cp) rho(it_chk,   :)
    write(iu_cp) kappaR(it_chk-1, :)
    write(iu_cp) kappaR(it_chk,   :)
    write(iu_cp) kappa_planck(it_chk-1, :)
    write(iu_cp) kappa_planck(it_chk,   :)
    write(iu_cp) tauR(it_chk-1, :)
    write(iu_cp) tauR(it_chk,   :)
    write(iu_cp) Qirr(it_chk-1, :)
    write(iu_cp) Qirr(it_chk,   :)
    write(iu_cp) dYdXi(it_chk-1, :)
    write(iu_cp) dYdXi(it_chk,   :)

    ! V3: Qvis, Qrad, is_shadow (needed for proper restart with irradiation)
    write(iu_cp) Qvis(it_chk-1, :)
    write(iu_cp) Qvis(it_chk,   :)
    write(iu_cp) Qrad(it_chk-1, :)
    write(iu_cp) Qrad(it_chk,   :)
    write(iu_cp) is_shadow(it_chk-1, :)
    write(iu_cp) is_shadow(it_chk,   :)

    ! V4: irradiation delay buffer (use_irradiation_delay only)
    call save_irradiation_buffer(iu_cp)

    close(iu_cp)

    ! --- Atomic-ish replace: remove old then rename temp to final
    inquire(file=fname, exist=ex)
    if (ex) then
       open(iu_cp, file=fname, status='old', iostat=ios)
       if (ios == 0) close(iu_cp, status='delete')
    end if
    call rename(trim(tmp), trim(fname))
  end subroutine write_checkpoint2


  subroutine read_checkpoint2(fname, it_restart, strict_dt)
    !! Read checkpoint and populate the last two time-slices in globals.
    character(len=*), intent(in) :: fname
    integer(i4b),     intent(in) :: it_restart
    logical,          intent(in) :: strict_dt

    integer(i4b), save :: iu_re = -1
    integer(i4b) :: ios
    integer(i4b) :: magic, ver, nr_file, it_file
    real(dp) :: t_file, dt_file, alpha_file
    real(dp), allocatable :: r_file(:)
    real(dp), allocatable :: sig_nm1(:), sig_n(:), nu_nm1(:), nu_n(:)
    real(dp), allocatable :: tmp(:)

    open(newunit=iu_re, file=fname, access='stream', form='unformatted', status='old', &
         action='read', iostat=ios)
    if (ios /= 0) then
       write(*,*) 'ERROR: cannot open checkpoint file: ', trim(fname)
       stop
    end if

    read(iu_re, iostat=ios) magic
    if (ios /= 0 .or. magic /= CKPT_MAGIC) then
       write(*,*) 'ERROR: invalid checkpoint magic.'
       stop
    end if

    read(iu_re) ver
    if (ver < 2 .or. ver > CKPT_VER) then
       write(*,*) 'ERROR: unsupported checkpoint version =', ver, ' (supported: 2-', CKPT_VER, ')'
       stop
    end if

    read(iu_re) nr_file
    if (nr_file /= nr) then
      write(*,*) 'ERROR: nr mismatch in checkpoint. nr_file=', nr_file, ' nr=', nr
      stop
    end if

    read(iu_re) it_file
    if (it_file /= it_restart) then
      write(*,*) 'ERROR: it mismatch in checkpoint. it_file=', it_file, ' it_restart=', it_restart
      stop
    end if

    read(iu_re) t_file
    read(iu_re) dt_file
    read(iu_re) alpha_file

    if (strict_dt) then
      if (abs(dt_file - dt) > 0.0_dp) then
        write(*,*) 'ERROR: dt mismatch (strict). dt_file=', dt_file, ' dt=', dt
        stop
      end if
    end if

    allocate(r_file(nr))
    read(iu_re) r_file(:)

    ! Enforce grid equality (recommended to be strict)
    if (maxval(abs(r_file(:) - r(:))) > 0.0_dp) then
      write(*,*) 'ERROR: r grid mismatch in checkpoint.'
      stop
    end if
    deallocate(r_file)

    allocate(sig_nm1(nr), sig_n(nr), nu_nm1(nr), nu_n(nr))
    read(iu_re) sig_nm1(:)
    read(iu_re) sig_n(:)
    read(iu_re) nu_nm1(:)
    read(iu_re) nu_n(:)

    ! Populate globals: we only guarantee slices it_restart-1 and it_restart.
    sigmat(:, :) = 0.0_dp
    nu_conv(:, :) = 0.0_dp

    sigmat(it_restart-1, :) = sig_nm1(:)
    sigmat(it_restart,   :) = sig_n(:)
    nu_conv(it_restart-1, :) = nu_nm1(:)
    nu_conv(it_restart,   :) = nu_n(:)

    ! Restore time at the checkpoint slice
    t_nd = t_file

    ! Optional structure (must exist in ver=2 format)
    allocate(tmp(nr))
    read(iu_re) tmp(:); Tmid(it_restart-1, :) = tmp(:)
    read(iu_re) tmp(:); Tmid(it_restart,   :) = tmp(:)
    read(iu_re) tmp(:); H(it_restart-1, :) = tmp(:)
    read(iu_re) tmp(:); H(it_restart,   :) = tmp(:)
    read(iu_re) tmp(:); rho(it_restart-1, :) = tmp(:)
    read(iu_re) tmp(:); rho(it_restart,   :) = tmp(:)
    read(iu_re) tmp(:); kappaR(it_restart-1, :) = tmp(:)
    read(iu_re) tmp(:); kappaR(it_restart,   :) = tmp(:)
    if (ver >= 5) then
       read(iu_re) tmp(:); kappa_planck(it_restart-1, :) = tmp(:)
       read(iu_re) tmp(:); kappa_planck(it_restart,   :) = tmp(:)
    else
       kappa_planck(it_restart-1, :) = kappaR(it_restart-1, :)
       kappa_planck(it_restart,   :) = kappaR(it_restart,   :)
    end if
    read(iu_re) tmp(:); tauR(it_restart-1, :) = tmp(:)
    read(iu_re) tmp(:); tauR(it_restart,   :) = tmp(:)
    read(iu_re) tmp(:); Qirr(it_restart-1, :) = tmp(:)
    read(iu_re) tmp(:); Qirr(it_restart,   :) = tmp(:)
    read(iu_re) tmp(:); dYdXi(it_restart-1, :) = tmp(:)
    read(iu_re) tmp(:); dYdXi(it_restart,   :) = tmp(:)

    ! V3: Qvis, Qrad, is_shadow (optional for ver 2)
    if (ver >= 3) then
       read(iu_re) tmp(:); Qvis(it_restart-1, :) = tmp(:)
       read(iu_re) tmp(:); Qvis(it_restart,   :) = tmp(:)
       read(iu_re) tmp(:); Qrad(it_restart-1, :) = tmp(:)
       read(iu_re) tmp(:); Qrad(it_restart,   :) = tmp(:)
       block
         logical, allocatable :: shad_tmp(:)
         allocate(shad_tmp(nr))
         read(iu_re) shad_tmp(:); is_shadow(it_restart-1, :) = shad_tmp(:)
         read(iu_re) shad_tmp(:); is_shadow(it_restart,   :) = shad_tmp(:)
         deallocate(shad_tmp)
       end block
    else
       ! Ver 2: Qvis, Qrad, is_shadow not in file; leave as-is or zero
       Qvis(it_restart-1:it_restart, :) = 0.0_dp
       Qrad(it_restart-1:it_restart, :) = 0.0_dp
       is_shadow(it_restart-1:it_restart, :) = .false.
    end if
    deallocate(tmp)

    ! V4: irradiation delay buffer
    if (ver >= 4) call load_irradiation_buffer(iu_re)

    close(iu_re)

    deallocate(sig_nm1, sig_n, nu_nm1, nu_n)
  end subroutine read_checkpoint2

end module checkpoint_mod

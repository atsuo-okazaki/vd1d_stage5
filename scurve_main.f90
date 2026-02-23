!======================================================================
! scurve_main.f90
!
! Standalone program to generate S-curves (Sigma vs T at thermal equilibrium)
! for thermal instability analysis. Use when opacity transitions near
! T ~ 10^4 K create thermally unstable regions.
!
! Qirr options:
!   qirr_mode = 'constant'  : use Qirr_in from namelist (0 or any value)
!   qirr_mode = 'from_file'  : read Qirr from disk_t*.dat at r_cgs
!
! Usage:
!   1. Create run with same setup as main simulation (ad1d.in, etc.)
!   2. Compile: make scurve
!   3. Run: ./scurve
!   4. Edit scurve.nml for r_cgs, Sigma range, Qirr mode, etc.
!
! Output: scurve_*.dat (Sigma, T, kappa, Qvis, Qrad, dQ+/dT, dQ-/dT, stability)
!======================================================================
program scurve_main
  use kind_params,   only : dp, i4b
  use constants,     only : gg, pi
  use mod_global,   only : M_star, R0, sigma_init, alphaSS
  use setup_mod,     only : setup
  use opacity_table_mod, only : init_opacity_tables
  use thermal_stability_analysis_mod, only : generate_scurve_at_radius
  use units_disk_mod, only : r_dim, omegaK_dim
  implicit none

  real(dp) :: r_cgs, r_nd
  real(dp) :: Sigma_min, Sigma_max
  real(dp) :: T_min, T_max
  real(dp) :: Sigma_sim, Tmid_sim
  real(dp) :: Qirr_in
  character(len=16) :: qirr_mode
  character(len=128) :: qirr_file
  integer(i4b) :: nSigma
  logical :: shadow
  character(len=64) :: outfile
  integer(i4b) :: iu
  namelist /scurve/ r_cgs, r_nd, Sigma_min, Sigma_max, nSigma, &
       T_min, T_max, Qirr_in, qirr_mode, qirr_file, shadow, outfile

  ! Defaults
  r_cgs     = 1.0e12_dp    ! [cm] e.g. ~7 AU for solar-type star
  r_nd      = -1.0_dp      ! if > 0, overrides r_cgs via r_cgs = r_dim(r_nd)
  Sigma_min = 1.0e-2_dp    ! [g/cm^2]
  Sigma_max = 1.0e4_dp     ! [g/cm^2]
  nSigma    = 200
  T_min     = 30.0_dp      ! [K] temperature range for root search
  T_max     = 1.0e5_dp     ! [K]
  Qirr_in   = 0.0_dp       ! [erg/cm^2/s] used when qirr_mode='constant'
  qirr_mode = 'constant'   ! 'constant' or 'from_file'
  qirr_file = 'disk_t00000100.dat'
  shadow    = .false.
  outfile   = 'scurve.dat'

  ! Setup builds grid, allocates globals, init_units
  call setup()

  ! Init opacity tables (required for heating_cooling_cell)
  call init_opacity_tables()

  ! Read namelist if present
  block
    logical :: ex
    inquire(file='scurve.nml', exist=ex)
    if (ex) then
       open(newunit=iu, file='scurve.nml', status='old', action='read')
       read(iu, nml=scurve)
       close(iu)
    end if
  end block

  if (r_nd > 0.0_dp) r_cgs = r_dim(r_nd)

  ! Set Qirr, Sigma_sim, Tmid_sim: constant or read from simulation output
  Sigma_sim = -1.0_dp
  Tmid_sim  = -1.0_dp
  if (trim(qirr_mode) == 'from_file') then
     call read_from_disk_file(trim(qirr_file), r_cgs, Qirr_in, Sigma_sim, Tmid_sim)
  end if

  write(*,'(a,1pe14.6,a)') 'Generating S-curve at r = ', r_cgs, ' [cm]'
  write(*,'(a,1pe14.6,a,1pe14.6,a)') 'Sigma range: ', Sigma_min, ' - ', Sigma_max, ' [g/cm^2]'
  write(*,'(a,1pe14.6,a,1pe14.6,a)') 'T range: ', T_min, ' - ', T_max, ' [K]'
  write(*,'(a,i0)') 'nSigma = ', nSigma
  write(*,'(a,a,a)') 'Qirr mode: ', trim(qirr_mode), '  => '
  if (trim(qirr_mode) == 'from_file') then
     write(*,'(a,a,a,1pe14.6,a)') '  from ', trim(qirr_file), ' => Qirr = ', Qirr_in, ' [erg/cm^2/s]'
  else
     write(*,'(a,1pe14.6,a)') '  Qirr = ', Qirr_in, ' [erg/cm^2/s]'
  end if
  write(*,'(a)') 'Output: '//trim(outfile)

  call generate_scurve_at_radius(r_cgs, omegaK_dim(r_cgs / R0), shadow, Qirr_in, &
       Sigma_min, Sigma_max, nSigma, trim(outfile), T_min=T_min, T_max=T_max, &
       Sigma_sim=Sigma_sim, Tmid_sim=Tmid_sim)

  write(*,'(a)') 'Done. Plot with: gnuplot> plot "scurve.dat" u 2:1 with lines'

contains

  !------------------------------------------------------------
  ! Read Qirr, Sigma, Tmid at r_target from disk structure output (disk_t*.dat).
  ! File format: # comment lines, then 15 columns per row:
  !   xi, r_cgs, sigma, H, H/r, dYdXi, rho, nu, Qvis, Qirr, Qrad, Tmid, Tsurf, kappaR, tauR
  ! Linear interpolation in r_cgs.
  !------------------------------------------------------------
  subroutine read_from_disk_file(fname, r_target, Qirr_out, Sigma_out, Tmid_out)
    character(len=*), intent(in) :: fname
    real(dp), intent(in) :: r_target
    real(dp), intent(out) :: Qirr_out, Sigma_out, Tmid_out

    integer(i4b) :: iu, ios, n
    character(len=256) :: line
    real(dp) :: xi, r_cgs, sigma, h_val, yr, dydx, rho, nu
    real(dp) :: qvis, qirr, qrad, tmid, tsurf, kap, tau
    real(dp), allocatable :: r_arr(:), qirr_arr(:), sigma_arr(:), tmid_arr(:)
    integer(i4b) :: i, i1, i2

    open(newunit=iu, file=trim(fname), status='old', action='read', iostat=ios)
    if (ios /= 0) then
       write(*,'(a,a,a)') 'ERROR: Cannot open file "', trim(fname), '"'
       stop
    end if

    n = 0
    do
       read(iu, '(a)', iostat=ios) line
       if (ios /= 0) exit
       if (len_trim(line) < 2 .or. line(1:1) == '#') cycle
       n = n + 1
    end do

    rewind(iu)
    allocate(r_arr(n), qirr_arr(n), sigma_arr(n), tmid_arr(n))

    n = 0
    do
       read(iu, '(a)', iostat=ios) line
       if (ios /= 0) exit
       if (len_trim(line) < 2 .or. line(1:1) == '#') cycle
       n = n + 1
       ! Columns: xi, r_cgs, sigma, H, H/r, dYdXi, rho, nu, Qvis, Qirr, Qrad, Tmid, Tsurf, kappaR, tauR
       read(line, *, iostat=ios) xi, r_cgs, sigma, h_val, yr, dydx, rho, nu, qvis, qirr, qrad, tmid, tsurf, kap, tau
       if (ios /= 0) then
          write(*,'(a,i0,a)') 'ERROR: Cannot parse line ', n, ' in file'
          stop
       end if
       r_arr(n) = r_cgs
       qirr_arr(n) = qirr
       sigma_arr(n) = sigma
       tmid_arr(n) = tmid
    end do
    close(iu)

    if (n == 0) then
       write(*,'(a)') 'ERROR: No data found in file'
       stop
    end if

    ! Linear interpolation in r_cgs
    if (r_target <= r_arr(1)) then
       Qirr_out  = qirr_arr(1)
       Sigma_out = sigma_arr(1)
       Tmid_out  = tmid_arr(1)
    else if (r_target >= r_arr(n)) then
       Qirr_out  = qirr_arr(n)
       Sigma_out = sigma_arr(n)
       Tmid_out  = tmid_arr(n)
    else
       Qirr_out = qirr_arr(1)
       Sigma_out = sigma_arr(1)
       Tmid_out = tmid_arr(1)
       do i = 1, n - 1
          if (r_arr(i) <= r_target .and. r_target <= r_arr(i+1)) then
             i1 = i
             i2 = i + 1
             Qirr_out = qirr_arr(i1) + (qirr_arr(i2) - qirr_arr(i1)) * &
                  (r_target - r_arr(i1)) / (r_arr(i2) - r_arr(i1))
             Sigma_out = sigma_arr(i1) + (sigma_arr(i2) - sigma_arr(i1)) * &
                  (r_target - r_arr(i1)) / (r_arr(i2) - r_arr(i1))
             Tmid_out = tmid_arr(i1) + (tmid_arr(i2) - tmid_arr(i1)) * &
                  (r_target - r_arr(i1)) / (r_arr(i2) - r_arr(i1))
             exit
          end if
       end do
    end if

    deallocate(r_arr, qirr_arr, sigma_arr, tmid_arr)

  end subroutine read_from_disk_file

end program scurve_main

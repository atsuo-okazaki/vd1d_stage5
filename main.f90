!======================================================================
! main.f90
!   1D viscous disk evolution code (theta-method diffusion solver)
!
!   Modules used:
!     kind_params   – kind definitions (i4b, dp)
!     mod_global    – global parameters and arrays
!     setup_mod     – grid setup and initial conditions
!     evolve_mod    – diffusion solver (time stepping)
!     output_mod    – write Sigma(r,t) data to text files (for Gnuplot)
!     inflow_params – mass inflow source term
!
!   Author : Atsuo Okazaki (2025, revised)
!======================================================================
program viscous_disk_main
  use kind_params,   only : i4b, dp
  use mod_global,    only : nr, nt, t_nd, dt, r, sigmat,          &
                            use_energy_balance, use_be_decretion, &
                            t_sim_end_eff
  use inflow_source_mod, only : compute_source
  use setup_mod,     only : setup, init_initial_conditions
  use opacity_table_mod, only : init_opacity_tables
  use evolve_try_mod,    only : evolve_try_to_target
  use state_io_mod,       only : load_state_from_history, store_state_to_history
  !use evolve_subcycle_mod,only : evolve_to_target
  use output_mod,    only : output_units_summary, output_mdot_profile, &
                            compute_disk_mass, output_disk_structure_dtl, &
                            output_edge_radius
  use output_run_summary_mod, only : output_run_summary
  use checkpoint_mod, only : read_checkpoint2, write_checkpoint2
  use checkpoint_util_mod, only : checkpoint_name
  use input_file_mod, only : infile
  use run_control_mod, only : do_restart, it_restart, chk_file, chkfreq, outfreq
  use timestep_seed_mod, only : reset_dt_seed
  use irradiation_mod, only : init_irradiation_buffer
  use disk_energy_mod, only : build_initial_disk_structure
  use hot_region_metrics_mod, only : init_hot_region_metrics
  implicit none

  namelist /run_control/ do_restart, it_restart, chk_file, chkfreq, outfreq

  integer(i4b) :: it, i, it_start
  real(dp), allocatable :: dr(:)
  real(dp), allocatable :: nu_new(:)
  real(dp) :: dt_out
  logical :: success

  ! 1) Setup always builds grid/dt/nt/allocates arrays, but does NOT create ICs.
  call setup()
  call init_irradiation_buffer(nt)
  call init_hot_region_metrics(Tcrit_in=1.0e4_dp, filename='hot_region_metrics.dat')
  call output_units_summary()
  
  ! 2) Read run_control AFTER setup (your requested policy)
  call try_read_run_control('run_control.nml', do_restart, it_restart, &
                            chk_file, chkfreq, outfreq)

  allocate(dr(nr))
  allocate(nu_new(nr))

  do i = 1, nr
     if (i == 1) then
        dr(i) = r(2) - r(1)
     else if (i == nr) then
        dr(i) = r(nr) - r(nr-1)
     else
        dr(i) = 0.5_dp * (r(i+1) - r(i-1))
     end if
  end do

  if (outfreq <= 0) outfreq = max(1_i4b, nt / 1000_i4b)
  if (chkfreq <= 0) chkfreq = outfreq

  ! 3) Initialize state (restart or fresh)
  if (do_restart) then
     call read_checkpoint2(chk_file, it_restart, strict_dt=.false.)
     it_start = it_restart
  else
     call init_initial_conditions()
     it_start = 1
  end if

  !------------------------------------------------------------
  ! 3) For a fresh run, build a consistent initial structure:
  !     Tmid(1,:), H(1,:), nu(:), Qirr(1,:) ...
  !     This avoids a violent transient at the first PDE step.
  !------------------------------------------------------------
  if (.not. do_restart) then
     ! Here nu_new is just a work array; update_nu_and_temp will also fill
     ! the global diagnostic fields at itp1=1 (Tmid,H,rho,kappaR,tauR,Qirr,...)
     call build_initial_disk_structure(nr, r, sigmat(1,:), 1)
  end if

  call output_run_summary('run_summary.txt', trim(infile), &
                          'run_control.nml', do_restart, it_restart)

  if (use_energy_balance) call init_opacity_tables()

  call reset_dt_seed(dt)

  ! Output cadence (uniform): always advance by dt_out in nd units
  dt_out = dt

  write(*,*) 'Starting diffusion evolution...'
  do it = it_start, nt-1

     call load_state_from_history(it)

     success = evolve_try_to_target(it, t_nd, dt_out, dr)

     if (.not. success) then
       write(*,'(a,i9,a,1pe12.4)') 'FATAL: time step failed at it=', it, ', t_nd=', t_nd
       stop
     end if

     t_nd = t_nd + dt_out

     call store_state_to_history(it+1)

     if (outfreq > 0 .and. mod(it, outfreq) == 0) then
        write(*,'(a,i9,a,1pe12.3,a,1pe12.3)') &
             'Step =', it, '  Time =', t_nd, ' / ', t_sim_end_eff
        call output_disk_structure_dtl(it+1, t_nd)
        !call output_disk_structure(it+1, t_nd)
        call output_mdot_profile(it+1, t_nd)
        call compute_disk_mass(it+1, t_nd)
        if (use_be_decretion) call output_edge_radius(it+1, t_nd)
     end if

     if (chkfreq > 0 .and. mod(it, chkfreq) == 0) then
        call write_checkpoint2(chk_file, it+1)
        call write_checkpoint2(checkpoint_name(chk_file, it+1), it+1)
     end if

  end do

  deallocate(dr)

contains

  subroutine try_read_run_control(fname, do_restart, it_restart, chk_file, chkfreq, outfreq)
    character(len=*), intent(in)    :: fname
    logical,          intent(inout) :: do_restart
    integer(i4b),     intent(inout) :: it_restart, chkfreq, outfreq
    character(len=*), intent(inout) :: chk_file
    integer(i4b), save :: iu_run = -1
    integer(i4b) :: ios
    logical :: ex

    inquire(file=fname, exist=ex)
    if (.not. ex) return

    open(newunit=iu_run, file=fname, status='old', action='read', iostat=ios)
    if (ios /= 0) return
    read(iu_run, nml=run_control, iostat=ios)
    close(iu_run)
  end subroutine try_read_run_control

end program viscous_disk_main

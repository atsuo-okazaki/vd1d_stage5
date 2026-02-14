  module run_control_mod
    use kind_params, only : dp, i4b
    implicit none

    ! Defaults
    logical      :: do_restart = .false.
    integer(i4b) :: it_restart = 1
    character(len=256) :: chk_file = 'checkpoint.bin'
    integer(i4b) :: chkfreq    = 0       ! checkpoint frequency; 0 = auto (reduced)
    integer(i4b) :: outfreq    = 0       ! disk structure (disk_t*.dat) for animation
    integer(i4b) :: thermal_stability_freq = 0  ! thermal stability output; 0 = auto (reduced)
    logical      :: do_thermal_stability_output = .false.  ! thermal instability analysis
    logical      :: do_thermal_stability_all_roots = .false.  ! find ALL roots at each r

    ! Debug / diagnostics (off by default for performance)
    logical      :: do_hot_region_metrics = .false.  ! output T>=10^4 K extent to hot_region_metrics.dat
    integer(i4b) :: iprint_irr = 0  ! m_irr convergence diagnostics; >=1 to enable

  end module run_control_mod
  
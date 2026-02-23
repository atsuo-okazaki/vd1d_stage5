module output_run_summary_mod
  use kind_params, only : dp, i4b
  implicit none
contains

  subroutine output_run_summary(fname, input_file, run_control_file, do_restart, it_restart)
    use constants, only : year, msun
    use mod_global, only : nr, nt, dt, alpha, t_nd, t_sim_start, t_sim_end, t_sim_end_eff, &
                           use_energy_balance, use_be_decretion, use_irradiation, &
                           use_irradiation_delay, use_finite_irradiation_source, use_energy_pde, &
                           inner_bc_type, outer_bc_type, p_nu_isothermal, &
                           M_star, R_star, Teff_star, q, R0, L_star, eta_acc, f_edd_cap, &
                           delta, alphaSS, hr0, nu0_nd, nu0_dim, fwhm, sigma_init, &
                           use_wind_truncation, mdot_w_msunyr, vinf_w, beta_w, f_rho_wind, &
                           use_inflow, mdot_inj_edd, mdot_inj_phys, mdot_inj_msunyr, mdot_inj_nd, &
                           rinj_min, rinj_max, t_in_start, t_in_end
    use irradiation_mod, only : L_irr
    use run_control_mod, only : do_thermal_stability_output, do_thermal_stability_all_roots, &
                                do_hot_region_metrics
    implicit none

    character(len=*), intent(in) :: fname
    character(len=*), intent(in) :: input_file
    character(len=*), intent(in) :: run_control_file
    logical,          intent(in) :: do_restart
    integer(i4b),     intent(in) :: it_restart

    integer(i4b) :: iu = -1
    integer(i4b) :: ios
    logical :: ex

    open(newunit=iu, file=trim(fname), status='replace', action='write', iostat=ios)
    if (ios /= 0) then
       write(*,*) 'ERROR: cannot open run summary file: ', trim(fname)
       return
    end if

    call write_header(iu)
    call write_scalar_summary(iu, do_restart, it_restart)

    call write_separator(iu, 'INPUT FILE (verbatim copy)')
    call copy_text_file(iu, input_file)

    call write_separator(iu, 'RUN CONTROL (verbatim copy, if exists)')
    call copy_text_file_if_exists(iu, run_control_file)

    call write_separator(iu, 'DERIVED / RUNTIME VALUES (selected)')
    write(iu,'(a)') '--- logical parameters (all) ---'
    write(iu,'(a,l1)') 'do_restart                  = ', do_restart
    write(iu,'(a,l1)') 'do_thermal_stability_output = ', do_thermal_stability_output
    write(iu,'(a,l1)') 'do_thermal_stability_all_roots = ', do_thermal_stability_all_roots
    write(iu,'(a,l1)') 'do_hot_region_metrics       = ', do_hot_region_metrics
    write(iu,'(a,l1)') 'use_energy_balance          = ', use_energy_balance
    write(iu,'(a,l1)') 'use_be_decretion            = ', use_be_decretion
    write(iu,'(a,l1)') 'use_irradiation             = ', use_irradiation
    write(iu,'(a,l1)') 'use_irradiation_delay       = ', use_irradiation_delay
    write(iu,'(a,l1)') 'use_finite_irradiation_source= ', use_finite_irradiation_source
    write(iu,'(a,l1)') 'use_energy_pde              = ', use_energy_pde
    write(iu,'(a,l1)') 'use_inflow                  = ', use_inflow
    write(iu,'(a,l1)') 'use_wind_truncation         = ', use_wind_truncation
    write(iu,'(a)') ''
    write(iu,'(a,1pe20.12)') 't_sim_start        = ', t_sim_start
    write(iu,'(a,1pe20.12)') 't_sim_end          = ', t_sim_end
    write(iu,'(a,1pe20.12)') 't_sim_end_eff      = ', t_sim_end_eff
    write(iu,'(a,1pe20.12)') 'dt (dimensionless) = ', dt
    write(iu,'(a,i12)')      'nt                 = ', nt
    write(iu,'(a,1pe20.12)') 'alpha (theta)      = ', alpha
    write(iu,'(a,1pe20.12)') 't_nd (current)     = ', t_nd

    write(iu,'(a)') '--- disk_mode ---'
    write(iu,'(a,l1)') 'use_energy_balance = ', use_energy_balance
    write(iu,'(a,l1)') 'use_be_decretion   = ', use_be_decretion
    write(iu,'(a,l1)') 'use_irradiation    = ', use_irradiation
    write(iu,'(a,i12)') 'inner_bc_type      = ', inner_bc_type
    write(iu,'(a,i12)') 'outer_bc_type      = ', outer_bc_type
    write(iu,'(a,1pe20.12)') 'p_nu_isothermal    = ', p_nu_isothermal

    write(iu,'(a)') '--- star / scale ---'
    write(iu,'(a,1pe20.12)') 'M_star [g]         = ', M_star
    write(iu,'(a,1pe20.12)') 'R_star [cm]        = ', R_star
    write(iu,'(a,1pe20.12)') 'Teff_star [K]      = ', Teff_star
    write(iu,'(a,1pe20.12)') 'q                  = ', q
    write(iu,'(a,1pe20.12)') 'R0 [cm]            = ', R0
    write(iu,'(a,1pe20.12)') 'L_star [erg/s]     = ', L_star
    write(iu,'(a,1pe20.12)') 'L_irr  [erg/s]     = ', L_irr
    write(iu,'(a,1pe20.12)') 'eta_acc            = ', eta_acc
    write(iu,'(a,1pe20.12)') 'f_edd_cap          = ', f_edd_cap

    write(iu,'(a)') '--- viscosity / init ---'
    write(iu,'(a,1pe20.12)') 'delta              = ', delta
    write(iu,'(a,1pe20.12)') 'alphaSS            = ', alphaSS
    write(iu,'(a,1pe20.12)') 'hr0                = ', hr0
    write(iu,'(a,1pe20.12)') 'nu0_nd             = ', nu0_nd
    write(iu,'(a,1pe20.12)') 'nu0_dim            = ', nu0_dim
    write(iu,'(a,1pe20.12)') 'fwhm               = ', fwhm
    write(iu,'(a,1pe20.12)') 'sigma_init         = ', sigma_init

    write(iu,'(a)') '--- inflow ---'
    write(iu,'(a,l1)')        'use_inflow         = ', use_inflow
    write(iu,'(a,1pe20.12)')  'mdot_inj_edd       = ', mdot_inj_edd
    write(iu,'(a,1pe20.12)')  'mdot_inj_phys [g/s]= ', mdot_inj_phys
    write(iu,'(a,1pe20.12)')  'mdot_inj_msunyr    = ', mdot_inj_msunyr
    write(iu,'(a,1pe20.12)')  'mdot_inj_nd        = ', mdot_inj_nd
    write(iu,'(a,1pe20.12)')  'rinj_min           = ', rinj_min
    write(iu,'(a,1pe20.12)')  'rinj_max           = ', rinj_max
    write(iu,'(a,1pe20.12)')  't_in_start         = ', t_in_start
    write(iu,'(a,1pe20.12)')  't_in_end           = ', t_in_end

    write(iu,'(a)') '--- wind ---'
    write(iu,'(a,l1)')        'use_wind_truncation= ', use_wind_truncation
    write(iu,'(a,1pe20.12)')  'mdot_w_msunyr      = ', mdot_w_msunyr
    write(iu,'(a,1pe20.12)')  'vinf_w [cm/s]      = ', vinf_w
    write(iu,'(a,1pe20.12)')  'beta_w             = ', beta_w
    write(iu,'(a,1pe20.12)')  'f_rho_wind         = ', f_rho_wind

    close(iu)
  end subroutine output_run_summary

  subroutine write_header(iu)
    integer(i4b), intent(in) :: iu
    write(iu,'(a)') '#============================================================'
    write(iu,'(a)') '# run_summary.txt (auto-generated)'
    write(iu,'(a)') '#============================================================'
  end subroutine write_header

  subroutine write_scalar_summary(iu, do_restart, it_restart)
    integer(i4b), intent(in) :: iu
    logical,      intent(in) :: do_restart
    integer(i4b), intent(in) :: it_restart
    write(iu,'(a)') ''
    write(iu,'(a)') '[RUN MODE]'
    write(iu,'(a,l1)') 'do_restart = ', do_restart
    write(iu,'(a,i12)') 'it_restart = ', it_restart
  end subroutine write_scalar_summary

  subroutine write_separator(iu, title)
    integer(i4b), intent(in) :: iu
    character(len=*), intent(in) :: title
    write(iu,'(a)') ''
    write(iu,'(a)') '#------------------------------------------------------------'
    write(iu,'(a)') '# ' // trim(title)
    write(iu,'(a)') '#------------------------------------------------------------'
  end subroutine write_separator

  subroutine copy_text_file(iu_out, fname)
    integer(i4b), intent(in) :: iu_out
    character(len=*), intent(in) :: fname
    integer(i4b) :: iu_in = -1
    integer(i4b) :: ios
    character(len=512) :: line

    open(newunit=iu_in, file=trim(fname), status='old', action='read', iostat=ios)
    if (ios /= 0) then
       write(iu_out,'(a)') '# [copy failed] cannot open: ' // trim(fname)
       return
    end if

    do
       read(iu_in,'(a)', iostat=ios) line
       if (ios /= 0) exit
       write(iu_out,'(a)') trim(line)
    end do
    close(iu_in)
  end subroutine copy_text_file

  subroutine copy_text_file_if_exists(iu_out, fname)
    integer(i4b), intent(in) :: iu_out
    character(len=*), intent(in) :: fname
    logical :: ex
    inquire(file=trim(fname), exist=ex)
    if (.not. ex) then
       write(iu_out,'(a)') '# [not found] ' // trim(fname)
       return
    end if
    call copy_text_file(iu_out, fname)
  end subroutine copy_text_file_if_exists

end module output_run_summary_mod

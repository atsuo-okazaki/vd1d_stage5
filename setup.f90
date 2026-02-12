!==============================================================
! setup_mod.f90
!   - Introduce t_sim_start/t_sim_end in /disk_time/
!   - Compute dt first using temporary r_tmp/nu_tmp arrays
!   - Compute nt from (t_sim_end - t_sim_start)/dt
!   - Allocate globals after nt is known
!==============================================================
module setup_mod
  use kind_params,  only : dp, i4b
  use constants, only : gg, msun, rsun, cc, pi, year, kb, sbc, mp, mu
  use mod_global,   only : nr, nt, r, nu, sigmat, Tmid, dt, alpha, t_nd, &
       t_sim_start, t_sim_end, t_sim_end_eff, &
       M_star, R_star, Teff_star, q, Temp0, L_star, R0, &
       fwhm, sigma_init, delta, alphaSS, hr0, nu0_nd, nu0_dim, &
       init_units, allocate_global, alloc_cur_state, print_global_info, &
       use_energy_balance, use_energy_pde, use_be_decretion, use_irradiation, &
       inner_bc_type, outer_bc_type, mdot_inj_nd, mdot_inj_phys, &
       mdot_inj_edd, mdot_inj_msunyr, eta_acc, f_edd_cap, &
       kappa_es, t0, p_nu_isothermal, rho_cut, r_edge, i_edge, &
       use_wind_truncation, mdot_w_msunyr, vinf_w, beta_w, f_rho_wind, &
       tau_irr_lag_nd
  use star_params_mod, only : &
       model_type, M_star_msun, R_star_rsun, Teff_star_in, &
       r0_mode, r0_over_Rstar, r0_over_Rg, R0_input,      &
       star_params, scale_params, set_default_star_params
  use inflow_source_mod, only : use_inflow, rinj_min, rinj_max, &
       t_in_start, t_in_end
  use input_file_mod, only : infile
  use units_disk_mod, only : t_dim, mdot_edd_unit, mdot_pde_unit, &
       mdot_from_edd, mdot_from_nd
  use irradiation_mod, only : L_irr
  use disk_energy_pde_mod, only : rebuild_structure_from_current_T
  implicit none

contains

  subroutine setup()

   integer(i4b) :: i, ios
    real(dp) :: rin_value, rout_value
    character(len=8) :: r_unit
    real(dp) :: tvisc_min, f_dt
    real(dp), parameter :: ln2  = 0.6931471805599453_dp
    real(dp), parameter :: tiny = 1.0e-30_dp

    ! Local namelist variables
    real(dp)     :: alpha_scheme
    !real(dp)     :: sigma_init
    real(dp)     :: dr_min, nu_max, C_stab

    ! Simulation time window (dimensionless)

    ! inflow (outer supply)

    ! disk_mode

    ! wind_params

    integer(i4b) :: narg

    ! Temporary arrays to compute dt and nt before allocate_global()
    real(dp), allocatable :: r_tmp(:), nu_tmp(:)

    !------------------------------------------------------------
    ! Namelists (nt removed; disk_time now defines run window)
    !------------------------------------------------------------
    namelist /disk_grid/ nr, rin_value, rout_value, r_unit
    namelist /disk_time/ alpha_scheme, t_sim_start, t_sim_end
    namelist /disk_init/ fwhm, sigma_init
    namelist /disk_visc/ delta, alphaSS, hr0
    namelist /disk_mode/ use_energy_balance, use_energy_pde, use_be_decretion, &
                     use_irradiation, inner_bc_type, outer_bc_type, &
                     p_nu_isothermal, tau_irr_lag_nd
    namelist /star_params/ model_type, M_star_msun, R_star_rsun, Teff_star_in, q
    namelist /scale_params/ r0_mode, r0_over_Rstar, r0_over_Rg, R0_input
    namelist /inflow/ use_inflow, mdot_inj_edd, rinj_min, rinj_max, &
                      t_in_start, t_in_end
    namelist /wind_params/ use_wind_truncation, mdot_w_msunyr, &
                           vinf_w, beta_w, f_rho_wind

    !------------------------------------------------------------
    ! 1. Defaults
    !------------------------------------------------------------
    nr           = 600
    rin_value    = 1.0_dp
    rout_value   = 100.0_dp
    r_unit       = 'dimless'

    alpha_scheme = 0.5_dp

    t_sim_start  = 0.0_dp
    t_sim_end    = 100.0_dp   ! default run duration in t' units

    fwhm         = 0.2_dp
    !sigma_init   = 0.0_dp

    delta        = 0.0_dp
    alphaSS      = 0.1_dp
    hr0          = 0.01_dp

    use_energy_balance = .true.
    use_be_decretion   = .false.
    use_irradiation    = .false.
    inner_bc_type      = 0
    outer_bc_type      = 0
    p_nu_isothermal    = 1.5_dp

    ! inflow window defaults (only used when use_inflow = .true.)
    use_inflow   = .false.
    mdot_inj_edd = 0.0_dp
    rinj_min     = 0.0_dp
    rinj_max     = 0.0_dp
    t_in_start   = 0.0_dp
    t_in_end     = 1.0e2_dp

    ! wind parameter defaults (only used when use_wind_truncation = .true.)
    use_wind_truncation = .false.
    mdot_w_msunyr = 0.0_dp
    vinf_w        = 2.5e8_dp
    beta_w        = 1.0_dp
    f_rho_wind    = 1.0_dp

    tau_irr_lag_nd = 0.0_dp   ! default: no lag

    call set_default_star_params()

    !------------------------------------------------------------
    ! 1b. Input file name
    !------------------------------------------------------------
    narg = command_argument_count()
    if (narg >= 1) call get_command_argument(1, infile)

    !------------------------------------------------------------
    ! 2. Read input file
    !------------------------------------------------------------
    open(unit=10, file=trim(infile), status='old', action='read', iostat=ios)
    if (ios == 0) then
       read(10, nml=disk_grid, iostat=ios)
       if (rin_value <= 0.0_dp) then
          write(*,*) '+++ rin_value must be positive. +++'
          stop
       end if
       if (r_unit /= 'dimless') then
          write(*,*) '+++ r_unit=''dimless'' is the only available unit. +++'
          stop
       end if

       read(10, nml=disk_time, iostat=ios)
       read(10, nml=disk_init, iostat=ios)
       read(10, nml=disk_visc, iostat=ios)

       read(10, nml=disk_mode, iostat=ios)
       if (ios /= 0) then
          write(*,*) 'Warning: namelist /disk_mode/ not found. Using defaults.'
          ios = 0
       end if

       read(10, nml=star_params, iostat=ios)
       if (ios /= 0) then
          write(*,*) 'Warning: namelist /star_params/ not found. Using defaults.'
          ios = 0
       end if

       read(10, nml=scale_params, iostat=ios)
       if (ios /= 0) then
          write(*,*) 'Warning: namelist /scale_params/ not found. Using defaults.'
          ios = 0
       end if

       read(10, nml=inflow, iostat=ios)
       if (ios /= 0) then
          write(*,*) 'Warning: namelist /inflow/ not found. Using defaults.'
          ios = 0
       end if

       read(10, nml=wind_params, iostat=ios)
       if (ios /= 0) then
          write(*,*) 'Warning: namelist /wind_params/ not found. Using defaults.'
          ios = 0
       end if

       close(10)
    else
       write(*,*) 'Warning: input file "', trim(infile), '" not found. Using defaults.'
    end if

    !------------------------------------------------------------
    ! 3. Copy scalar params to globals (no self-assignments)
    !------------------------------------------------------------
    alpha            = alpha_scheme

    !------------------------------------------------------------
    ! 4. Physical stellar parameters and R0, then init_units()
    !------------------------------------------------------------
    M_star    = M_star_msun * msun
    R_star    = R_star_rsun * rsun

    select case (trim(r0_mode))
    case ('star')
      R0 = r0_over_Rstar * R_star
    case ('rg')
      R0 = r0_over_Rg * (2.0_dp * gg * M_star / (cc*cc))
    case ('explicit')
      if (R0_input > 0.0_dp) then
        R0 = R0_input
      else
        write(*,*) 'Warning: r0_mode=''explicit'' but R0_input <= 0. Using R0 = R_star.'
        R0 = R_star
      end if
    case default
      write(*,*) 'Warning: unknown r0_mode. Using R0 = R_star.'
      R0 = R_star
    end select

    if (model_type == 'bh') then
       ! Teff_star is the virial temperature at R0
       Teff_star = mu * mp / kb * gg * M_star / R0
    else 
       Teff_star = Teff_star_in
    end if
    ! Just in case temperature has be normalized
    Temp0 = Teff_star
 
    L_star = 4.0_dp * pi * sbc * R0 * R0 * Teff_star**4

    call init_units()

    !------------------------------------------------------------
    ! 5. Mass-flow conversions (single source of truth)
    !------------------------------------------------------------
    mdot_inj_phys   = mdot_from_edd(mdot_inj_edd)     ! [g/s]
    mdot_inj_msunyr = mdot_inj_phys / msun * year
    mdot_inj_nd     = mdot_inj_phys / mdot_pde_unit() ! [-]

    !------------------------------------------------------------
    ! 6. Build temporary r and nu to compute dt before allocation
    !------------------------------------------------------------
    allocate(r_tmp(nr), nu_tmp(nr))

    do i = 1, nr
       r_tmp(i) = rin_value * (rout_value/rin_value)** &
                  ( real(i-1,dp) / real(nr-1,dp) )
    end do

    nu0_nd = 2.0_dp / 3.0_dp * alphaSS * hr0**2
    do i = 1, nr
       nu_tmp(i) = nu0_nd * r_tmp(i)**(delta + 1.5_dp)
    end do

    !------------------------------------------------------------
    ! 7. Compute dt (dimensionless)
    !------------------------------------------------------------
    if (alpha < 0.5_dp) then
       dr_min = 1.0e99_dp
       do i = 1, nr-1
          dr_min = min(dr_min, r_tmp(i+1) - r_tmp(i))
       end do
       nu_max = maxval(nu_tmp)
       C_stab = 0.1_dp
       dt = C_stab * dr_min * dr_min / max(nu_max, 1.0e-30_dp)
    else
       tvisc_min = 1.0e99_dp
       do i = 1, nr
          tvisc_min = min(tvisc_min, r_tmp(i)**2 / max(nu_tmp(i), 1.0e-30_dp))
       end do
       if (use_be_decretion) then
          f_dt = 1.0e-4_dp
       else
          f_dt = 1.0e-2_dp
       end if
       dt = f_dt * tvisc_min
    end if

    if (dt <= 0.0_dp) then
       write(*,*) 'ERROR: dt <= 0 after estimation.'
       stop
    end if

    !------------------------------------------------------------
    ! 8. Compute nt from (t_sim_end - t_sim_start) and dt
    !------------------------------------------------------------
    if (t_sim_end <= t_sim_start) then
       write(*,*) 'ERROR: t_sim_end must be > t_sim_start.'
       stop
    end if

    nt = int( (t_sim_end - t_sim_start)/dt ) + 1
    if (nt < 2) nt = 2

    t_sim_end_eff = t_sim_start + dt * real(nt-1, dp)

    !------------------------------------------------------------
    ! 9. Allocate globals and copy tmp arrays
    !------------------------------------------------------------
    call allocate_global(nr, nt)
    call alloc_cur_state(nr)

    r(:)  = r_tmp(:)
    nu(:) = nu_tmp(:)

    deallocate(r_tmp, nu_tmp)

    !------------------------------------------------------------
    ! 10. After allocation, keep arrays in a safe state.
    ! Initial conditions are constructed in init_initial_conditions()
    ! (fresh run) or overwritten by checkpoint (restart).
    !------------------------------------------------------------
    sigmat(:, :) = 0.0_dp
    t_nd         = t_sim_start

    !------------------------------------------------------------
    ! 11. Debug prints
    !------------------------------------------------------------
    write(*,*) '----------------------------------------------------------'
    write(*,*) 'Setup summary (dimensionless)'
    write(*,'(a)')            'Input file: '//trim(infile)
    write(*,'(a,1pe12.4)')    't_sim_start  = ', t_sim_start
    write(*,'(a,1pe12.4)')    't_sim_end    = ', t_sim_end
    write(*,'(a,1pe12.4)')    't_sim_end_eff= ', t_sim_end_eff
    write(*,'(a,i9)')         'nr = ', nr
    write(*,'(a,i9)')         'nt = ', nt
    write(*,'(a,1pe12.4)')    'dt'' = ', dt
    write(*,'(a,1pe12.4,a,1pe12.4)') 'r'' in [', r(1), ', ', r(nr), ']'
    write(*,'(a,1pe12.4)')    'mdot_inj_edd = ', mdot_inj_edd
    write(*,'(a,1pe12.4)')    'mdot_inj_msunyr = ', mdot_inj_msunyr
    write(*,'(a,1pe12.4)')    'mdot_inj_nd  = ', mdot_inj_nd
    write(*,*) '----------------------------------------------------------'

    call print_global_info()

  end subroutine setup

  subroutine init_initial_conditions()
    use kind_params, only : dp, i4b
    use constants,   only : pi
    use mod_global,  only : nr, nt, r, sigmat, t_nd, t_sim_start, &
                            use_inflow, use_be_decretion, fwhm, sigma_init
    implicit none

    integer(i4b) :: i
    real(dp), parameter :: ln2  = 0.6931471805599453_dp
    real(dp), parameter :: tiny = 1.0e-30_dp

    ! Clear all slices
    sigmat(:, :) = 0.0_dp

    ! Build Sigma(r, t=t_sim_start)
    if (use_inflow) then
       sigmat(1, :) = 0.0_dp
    else if (use_be_decretion) then
       sigmat(1, :) = 0.0_dp
    else
       do i = 1, nr
          sigmat(1, i) = sigma_init * exp( -4.0_dp*ln2 * ((r(i) - 1.0_dp)/fwhm)**2 )
          if (sigmat(1, i) < tiny) sigmat(1, i) = 0.0_dp
       end do
    end if

    ! Initialize Tmid at t=1
    !Tmid(1,:) = Temp0   ! or Teff_star, or a power-law seed
    Tmid(1,:) = 0.0_dp   ! or Teff_star, or a power-law seed
    !nu_tmp(:) = nu(:)   ! local array if needed; or just reuse nu(:)

    call rebuild_structure_from_current_T(1, sigmat(1,:), nu(:))

    t_nd = t_sim_start

    ! Keep remaining slices defined (optional but safe)
    if (nt > 1) sigmat(2:nt, :) = 0.0_dp

end subroutine init_initial_conditions

end module setup_mod

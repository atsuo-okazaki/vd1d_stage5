!==============================================================
! output_mode.f90
!   Write disk surface-density evolution to text files
!   (for Gnuplot or other plotting tools)
!==============================================================
module output_mod
  use kind_params,  only : dp, i4b, lgt
  use constants,    only : pi, kb, mp, mu, msun, year, sbc
  use mod_global,   only : nr, r, nu, sigmat, Tmid,     &
                           H, rho, kappaR, kappa_planck, tauR, dt, fwhm,    &
                           delta, alpha, kappa_es, kappa0, &
                           M_star, R_star, Teff_star, L_star, R0, &
                           t0, sigma_init, Temp0, dYdxi, &
                           alphaSS, hr0, nu0_nd, nu0_dim, &
                           nu_conv, k_iter, m_iter, r_edge, i_edge, &
                           mdot_inj_edd, mdot_inj_phys, mdot_inj_msunyr, &
                           mdot_inj_nd, Qvis, Qirr, Qrad, &
                           mdot_w_msunyr, vinf_w, beta_w, f_rho_wind, &
                           outer_bc_type, use_energy_balance, &
                           use_irradiation, is_shadow
  use setup_mod, only : model_type
  use units_disk_mod, only : t_dim, r_dim, sigma_dim, omegaK_dim
  use mdot_units_mod, only : mdot_edd_unit, mdot_pde_unit

contains

  !------------------------------------------------------------
  !  output_units_summary
  !
  !  Writes a summary of the physical and dimensionless units
  !  used in the current run.
  !
  !  File: units_summary.dat
  !
  !  Assumes:
  !    - M_star, R0, sigma_init, Temp0, alphaSS, hr0, delta, nu0_nd
  !      have been set (and init_units() has been called, so
  !      t0 corresponds to the Keplerian time at R0).
  !    - r(:) and dt are the dimensionless grid and time step.
  !------------------------------------------------------------
  subroutine output_units_summary()
    integer(i4b), save :: iu_sum = -1
    real(dp) :: Omega0, P0, P0_days
    real(dp) :: R0_AU
    real(dp) :: M_star_Msun
    real(dp) :: Mdot_Edd, Mdot_Edd_Msun_yr
    real(dp) :: rmin_nd, rmax_nd, dt_nd
    real(dp), parameter :: au = 1.495978707e13_dp  ! [cm]

    ! Basic derived quantities
    Omega0 = 1.0_dp / t0
    P0     = 2.0_dp * pi / Omega0
    P0_days = P0 / (86400.0_dp)

    R0_AU      = R0 / au
    M_star_Msun = M_star / msun

    Mdot_Edd       = mdot_edd_unit()           ! [g/s]
    Mdot_Edd_Msun_yr = Mdot_Edd / msun * year  ! [Msun/yr]

    rmin_nd = r(1)
    rmax_nd = r(nr)
    dt_nd   = dt

    open(newunit=iu_sum, file='units_summary.dat', status='replace', action='write')

    write(iu_sum,'(a)') '# ===================================================='
    write(iu_sum,'(a)') '# Units summary for 1D viscous disk simulation'
    write(iu_sum,'(a)') '# ===================================================='
    write(iu_sum,'(a)') '# BASIC PHYSICAL SCALES'
    write(iu_sum,'(a,1pe15.7)') '#   M_star [g]      = ', M_star
    write(iu_sum,'(a,1pe15.7)') '#   M_star [Msun]   = ', M_star_Msun
    write(iu_sum,'(a,1pe15.7)') '#   R0 [cm]         = ', R0
    write(iu_sum,'(a,1pe15.7)') '#   R0 [AU]         = ', R0_AU
    write(iu_sum,'(a,1pe15.7)') '#   t0 [s]          = ', t0
    write(iu_sum,'(a,1pe15.7)') '#   sigma_init [g/cm^2] = ', sigma_init
    write(iu_sum,'(a,1pe15.7)') '#   Temp0 [K]       = ', Temp0

    write(iu_sum,'(a)') '# '
    write(iu_sum,'(a)') '# DYNAMICAL QUANTITIES AT R0'
    write(iu_sum,'(a,1pe15.7)') '#   Omega0 [1/s]    = ', Omega0
    write(iu_sum,'(a,1pe15.7)') '#   P0 [s]          = ', P0
    write(iu_sum,'(a,1pe15.7)') '#   P0 [days]       = ', P0_days

    write(iu_sum,'(a)') '# '
    write(iu_sum,'(a)') '# ALPHA-DISK PARAMETERS'
    write(iu_sum,'(a,1pe15.7)') '#   alphaSS            = ', alphaSS
    write(iu_sum,'(a,1pe15.7)') '#   hr0 = H/R at R0    = ', hr0
    write(iu_sum,'(a,1pe15.7)') '#   delta (nu index)   = ', delta
    write(iu_sum,'(a,1pe15.7)') '#   nu0_nd (dimless)   = ', nu0_nd

    write(iu_sum,'(a)') '# '
    write(iu_sum,'(a)') '# OPACITY PARAMETERS'
    write(iu_sum,'(a,1pe15.7)') '#   kappa_es [cm^2/g] = ', kappa_es
    write(iu_sum,'(a,1pe15.7)') '#   kappa0   [cgs]    = ', kappa0

    write(iu_sum,'(a)') '# '
    write(iu_sum,'(a)') '# EDDINGTON-BASED ACCRETION SCALE'
    write(iu_sum,'(a,1pe15.7)') '#   Mdot_unit [g/s]      = ', Mdot_Edd
    write(iu_sum,'(a,1pe15.7)') '#   Mdot_unit [Msun/yr]  = ', Mdot_Edd_Msun_yr
    write(iu_sum,'(a)') '# '
    write(iu_sum,'(a)') '# MASS INJECTION RATE'
    write(iu_sum,'(a,1pe15.7)') '#   mdot_inj_phys [g/s]       = ', mdot_inj_phys
    write(iu_sum,'(a,1pe15.7)') '#   mdot_inj_msunyr [Msun/yr] = ', mdot_inj_msunyr
    write(iu_sum,'(a,1pe15.7)') '#   mdot_inj_edd [Edd scale]  = ', mdot_inj_edd

    write(iu_sum,'(a)') '# '
    write(iu_sum,'(a)') '# DIMENSIONLESS GRID AND TIME'
    write(iu_sum,'(a,1pe15.7,a,1pe15.7,a)') '#   r'' in [', rmin_nd, ', ', rmax_nd, ']'
    write(iu_sum,'(a,1pe15.7)') '#   dt'' (dimensionless) = ', dt_nd

    write(iu_sum,'(a)') '# ===================================================='

    close(iu_sum)

  end subroutine output_units_summary
  
  !------------------------------------------------------------
  !  output_sigmat
  !
  !  Parameters:
  !    nt_out  : number of time steps that have been computed
  !    outfreq : output interval in time steps
  !
  !  Output file:
  !    disk_evolution.dat
  !
  !  Format:uic
  !    Header with model parameters
  !    Then columns:
  !      r   Sigma(t=0)  Sigma(t=outfreq)  Sigma(t=2*outfreq) ...
  !------------------------------------------------------------
  subroutine output_sigmat(nt_out, outfreq)
    integer(i4b), intent(in) :: nt_out, outfreq
    integer(i4b), save :: idisk1 = -1
    integer(i4b) :: i, it, nstep
    real(dp)     :: t
    character(len=*), parameter :: outfile = 'disk_evolution.dat'

    if (nt_out <= 0) return
    if (outfreq <= 0) return

    nstep = (nt_out - 1) / outfreq + 1

    open(newunit=idisk1, file = outfile, status = 'replace', &
         action = 'write')

    write(idisk1, "(a)") '# 1D viscous disk evolution'
    write (idisk1, '(a, 1pe11.3, a, 1pe11.3, a, 1pe11.3)') &
          '# M_star =', M_star, ', R_star =', R_star, &
          ', Teff_star =', Teff_star
    write(idisk1, "(a,1pe11.3,a,1pe11.3,a,1pe11.3)") &
         '# r_min =', r(1), ', r_max =', r(nr), ', r0 =', r0
    write(idisk1, "(a,1pe11.3)") &
         '# FWHM (in r0 units) =', fwhm
    write(idisk1, "(a,1pe11.3)") &
         '# nu0_nd =', nu0_nd
    write(idisk1, "(a,f7.3)") &
         '# delta =', delta
    write(idisk1, "(a,1pe11.3)") &
         '# alpha (time integration parameter) =', alpha
    write(idisk1, "(a,1pe11.3)") &
         '# dt =', dt
    write(idisk1, "(a,i6,a,i7,a,i6)") &
         '# nr =', nr, ', nt_out =', nt_out, ', outfreq =', outfreq

    write(idisk1, "(a)", advance='no') '# sample times:'
    do it = 1, min(nt_out, 5*outfreq), outfreq
       t = dt * real(it-1, dp)
       write(idisk1, "(1x,1pe11.3)", advance='no') t
    end do
    write(idisk1, "('  ...')")

    write(idisk1, "(a)") '# r   Sigma(t=0)  Sigma(t=outfreq)  Sigma(t=2*outfreq) ... in CGS units'

    do i = 1, nr
       write(idisk1, "(1p,1e13.5)", advance='no') r_dim(r(i))
       do it = 1, nt_out, outfreq
          write(idisk1, "(1x,1p,1e13.5)", advance='no') &
               sigma_dim(sigmat(it, i))
       end do
       write(idisk1, *)
    end do

    close(idisk1)
  end subroutine output_sigmat


  !------------------------------------------------------------
  !  output_time_history
  !
  !  (Optional) Write simple time history of total mass, etc.
  !  For now this is a placeholder. You can fill in Mdisk(t)
  !  using a trapezoidal integrator over sigmat(it,:).
  !------------------------------------------------------------
  subroutine output_time_history()
    ! This routine can be implemented once you have
    ! a mass_trapz function and an array of saved times.
  end subroutine output_time_history

  
  subroutine output_mdot_profile(it, t)
    use mod_global, only : mdot_inner_nd, mdot_inner_phys
    use irradiation_mod, only : L_irr
    integer(i4b), intent(in) :: it
    real(dp),    intent(in)  :: t

    integer(i4b), save :: iu_his = -1
    logical(lgt), save :: first = .true.

    ! Accretion / zero-torque BC: recover flux from G-gradient
    ! Assume G = nu*Sigma*sqrt(r) = 0 at the inner edge.
    
    if (first) then
       open(newunit=iu_his, file='mdot_inner.dat', status='replace', &
            action='write')
       write (iu_his, "('# t0 =', 1pe12.5, ', Ledd/c^2 =', 1pe12.5)") &
            t0, mdot_edd_unit()
       write (iu_his, '(a, 1pe11.3, a, 1pe11.3, a, 1pe11.3, a, 1pe11.3)') &
          '# M_star =', M_star, ', R_star =', R_star, &
          ', Teff_star =', Teff_star, ', R0 =', R0
       !if (use_irradiation) then
          write (iu_his, '("#", 5x, "t/t0", 12x, "t[s]", 6x, "mdot_inner_nd", &
                 & 1x, "mdot_inner[g/s]", 2x, "L_irr[erg/s]")')
       !else
       !   write (iu_his, '("#", 5x, "t/t0", 12x, "t[s]", 6x, "mdot_inner_nd", &
       !          & 1x, "mdot_inner[g/s]")')
       !end if
       first = .false.
    end if
    !if (use_irradiation) then
       write(iu_his,'(1p,5e15.7)') t, t_dim(t), mdot_inner_nd, mdot_inner_phys, L_irr
    !else
    !   write(iu_his,'(1p,4e15.7)') t, t_dim(t), mdot_inner_nd, mdot_inner_phys
    !end if

  end subroutine output_mdot_profile

  subroutine compute_disk_mass(it, t)
    integer(i4b), intent(in) :: it
    real(dp),    intent(in)  :: t

    integer(i4b) :: i
    real(dp)     :: mass, r1_phys, r2_phys, dr_phys, &
         sigmat1_phys, sigmat2_phys
    integer(i4b), save :: iu_md = -1
    logical(lgt), save :: first = .true.

    mass = 0.0_dp
    do i = 1, nr-1
       r1_phys = r_dim(r(i))
       r2_phys = r_dim(r(i+1))
       dr_phys = r2_phys - r1_phys
       sigmat1_phys = sigma_dim(sigmat(it,i))
       sigmat2_phys = sigma_dim(sigmat(it,i+1))
       mass = mass + 0.5_dp * (2.0_dp * pi) * &
           (r1_phys * sigmat1_phys  + &
            r2_phys * sigmat2_phys) * dr_phys
    end do

    if (first) then
       open(newunit=iu_md, file='Mdisk_vs_t.dat', status='replace', &
            action='write')
       write (iu_md, "('# t0 =', 1pe12.5, ', Ledd/c^2 =', 1pe12.5)") &
            t0, mdot_edd_unit()
       write (iu_md, "('#', 5x, 't/t0', 12x, 't[s]', 8x, 'Mdisk[g]', &
            & 5x, 'Mdisk/(Ledd/c^2)[s]')")
       first = .false.
    end if
    write(iu_md,'(1p,4e15.7)') t, t_dim(t), mass, mass/mdot_edd_unit()

  end subroutine compute_disk_mass

  subroutine output_edge_radius(it, t)
    integer(i4b), intent(in) :: it
    real(dp),    intent(in)  :: t
    integer(i4b), save :: iu_ed = -1
    logical(lgt), save :: first = .true.

    if (first) then
       open(newunit=iu_ed, file='r_edge_vs_t.dat', status='replace', &
            action='write')
       write (iu_ed, '(a)') '# Radius to be truncated by wind-driven ablation'
       write (iu_ed, '(a, 1pe12.4, a, 1pe12.4, a, 1pe12.5, a)') &
              '# M_star =', M_star, '[g], R_star =', R_star, &
              '[cm], Teff_star =', Teff_star, '[K]'
       write (iu_ed, '(a, 1pe12.4, a, 1pe12.4, a, 1pe12.5, a)') &
              '# Mdot_wind =', mdot_w_msunyr, '[Msun/yr], beta =', &
              beta_w, ', v_inf =', vinf_w, '[cm/s]'
       write (iu_ed, '(a, 1pe12.4, a, 1pe12.4)') '# t0 =', t0, ', R0 =', R0
       write (iu_ed, '(a, 5x, a, 12x, a, 8x, a, 3x, a)') '#', &
              't/t0', 't[s]', 'r_edge[cm]', 'r_edge/R_star'
       first = .false.
    end if
    write(iu_ed,'(1p,4e15.7)') t, t_dim(t), r_edge(it), r_edge(it)/R_star
    
  end subroutine output_edge_radius

  subroutine output_disk_structure_dtl(it, t, rQ_rms, rQ_p95, rQ_max, &
                                       rdt_max, rdt_rms, rdt_p95, &
                                       qbal_max, qbal_p95)
    use mod_global, only : mdot_inner_phys
    use irradiation_mod, only : L_irr
    integer(i4b), intent(in) :: it
    real(dp),    intent(in)  :: t
    real(dp),     intent(in) :: rQ_rms, rQ_p95, rQ_max
    real(dp),     intent(in) :: rdt_max, rdt_rms, rdt_p95
    real(dp),     intent(in) :: qbal_max, qbal_p95

    integer(i4b), save :: iu_st = -1
    integer(i4b) :: i
    character(len=64) :: fname
    real(dp) :: Qrad_cell(nr), Qvis_cell(nr), Qirr_cell(nr)
    real(dp) :: t_cgs, r_cgs, sigma_cgs, omegaK_cgs
    real(dp) :: nu_dim_loc, T_loc
    real(dp) :: kappaR_loc, tauR_loc, kappa_out, tau_out
    real(dp) :: Qvis_out, Qrad_out
    real(dp) :: T_emit
    real(dp) :: cs2
    real(dp), parameter :: tiny_sigma = 1.0e-8_dp

    ! ---- SED cell work arrays (reused to avoid recomputation)
    real(dp) :: Sigma_cgs_cell(nr), Tspec_cell(nr), kapP_cell(nr)

    t_cgs = t_dim(t)
    
    write(fname,'("disk_t",i8.8,".dat")') it
    open(newunit=iu_st, file=fname, status='replace')

    write(iu_st,'(a, i8, a, 1pe13.5, a)') '# it =', it, ', t =', t_cgs, '[s]'
    write(iu_st,'(A, 1pe12.4, A, 1pe12.4, A, 1pe12.4)') &
        &  '# irr_fixedpoint_residual (full step aggregate): rms =', rQ_rms, &
        &  '  p95 =', rQ_p95, '  max =', rQ_max
    write(iu_st,'(A, 1pe12.4, A, 1pe12.4, A, 1pe12.4)') &
        &  '# Time-integration residuals: rms =', rdt_rms, &
        &  '  p95 =', rdt_p95, '  max =', rdt_max
    write(iu_st,'(A, 1pe12.4, A, 1pe12.4)') &
        &  '# Energy-balance residuals n time-integration eq.: p95 =', qbal_p95, &
        &  '  max =', qbal_max
    if (model_type == 'star') then
       write (iu_st, '(a, 1pe13.5, a, 1pe13.5, a, 1pe13.5, a, 1pe13.5, a)') &
          '# M_star =', M_star, '[g], R_star =', R_star, &
          '[cm], Teff_star =', Teff_star, '[K], R0 =', R0, '[cm]'
    else
       ! mdoel_type = 'bh'
       write (iu_st, '(a, 1pe13.5, a, 1pe13.5, a, 1pe13.5, a, 1pe13.5, a)') &
          '# M_star =', M_star, '[g], R0 =', R0, '[cm], Mdot(R0) =', &
          mdot_inner_phys, '[g/s], L_irr =', L_irr, '[erg/s]'
    end if
    write(iu_st,'(a)') '# kappa: kappaR if tau>=1, kappaP if tau<1. tau: tauR if tau>=1, tauP if tau<1. T_emit: (Qrad/2sbc)^0.25 if tau>=1, Tmid if tau<1 (SED-ready)'
    write(iu_st,'(a, 2x, a, 6x, a, 4x, a, 4x, a, 7x, a, 7x, a, 4x, a, 2x, a, 6x, &
         & a, 9x, a, 9x, a, 8x, a, 5x, a, 3x, a, 3x, a)') '#', 'xi(=r/R0)', 'r(cm)', &
         & 'sigma(g/cm^2)', 'H(cm)', 'Y(=H/r)', 'dY/dxi', &
         & 'rho(g/cm^3)', 'nu(cm^2/s)', &
         & 'Qvis', 'Qirr', 'Qrad', 'Tmid(K)', &
         & 'T_emit(K)', 'kappa(cm^2/g)', 'tau'

    ! Initialization of SED cell work arrays
    Sigma_cgs_cell(:) = 0.0_dp
    Tspec_cell(:)     = 0.0_dp
    kapP_cell(:)      = 0.0_dp
    Qrad_cell(:) = 0.0_dp
    Qvis_cell(:) = 0.0_dp
    Qirr_cell(:) = 0.0_dp

    do i = 1, nr
       r_cgs    = r_dim(r(i))
       Sigma_cgs = sigma_dim(sigmat(it, i))
       OmegaK_cgs = omegaK_dim(r(i))
       T_loc = Tmid(it, i)
       nu_dim_loc = nu_conv(it,i) * nu0_dim / nu0_nd
       kappaR_loc = kappaR(it, i)
       tauR_loc = tauR(it, i)

       if (use_energy_balance) then
          !-----------------------------------------------
          ! If Sigma is (numerically) zero, treat as empty
          !-----------------------------------------------
          if (sigma_cgs <= tiny_sigma .or. Tmid(it,i) <= 0.0_dp) then
             H(it, i)      = 0.0_dp
             rho(it, i)    = 0.0_dp
             kappaR(it, i) = 0.0_dp
             tauR(it, i)   = 0.0_dp
          
             kappa_out = 0.0_dp
             tau_out   = 0.0_dp
             write(iu_st,'(1p,15e13.5)') r(i), r_cgs, sigma_cgs, H(it,i), &
                  H(it,i)/r_cgs, 0.0_dp, rho(it,i), nu_dim_loc, &
                  0.0_dp, 0.0_dp, 0.0_dp, &  ! Qvis, Qirr, Qrad
                  0.0_dp, 0.0_dp, 0.0_dp, &  ! Tmid, T_emit, kappa
                  0.0_dp                     ! tau

             Sigma_cgs_cell(i) = 0.0_dp
             kapP_cell(i)      = 0.0_dp
             Tspec_cell(i)     = 0.0_dp
             Qrad_cell(i)      = 0.0_dp
             Qvis_cell(i)      = 0.0_dp
             Qirr_cell(i)      = 0.0_dp
             
             cycle
          end if

          ! For SED: tau>=1 use kappaR/tauR (optically thick); tau<1 use kappaP/tauP (optically thin)
          if (tauR_loc >= 1.0_dp) then
             kappa_out = kappaR(it, i)
             tau_out   = tauR(it, i)
          else
             kappa_out = kappa_planck(it, i)
             tau_out   = 0.5_dp * kappa_planck(it, i) * Sigma_cgs
          end if

          ! Use saved Qvis/Qrad (consistent with thermal solver)
          Qvis_out = Qvis(it, i)
          Qrad_out = Qrad(it, i)
          Qvis_cell(i) = Qvis_out
          Qrad_cell(i) = Qrad_out
          Qirr_cell(i) = Qirr(it, i)
          ! For SED: tau>=1 use effective surface T from flux; tau<1 use Tmid (optically thin)
          if (tauR_loc >= 1.0_dp) then
             T_emit = (Qrad_out / (2.0_dp * sbc))**0.25_dp
          else
             T_emit = Tmid(it, i)
          end if

          ! ---- Cache SED inputs for reuse (cell-centered)
          Sigma_cgs_cell(i) = Sigma_cgs
          kapP_cell(i)      = kappa_planck(it, i)   ! absorption-only Planck mean assumed
          Tspec_cell(i)     = T_emit                ! SED-ready temperature definition

       else
          ! Isothermal-disk quantities (no thermal balance; Qvis/Qrad not defined)
          Qvis_out = 0.0_dp
          Qrad_out = 0.0_dp
          if (sigma_cgs <= tiny_sigma) then
             nu_conv(it, i) = nu(i)
             H(it, i)     = 0.0_dp
             rho(it, i)   = 0.0_dp
             Tmid(it, i)  = 0.0_dp
             T_emit       = 0.0_dp

             Sigma_cgs_cell(i) = 0.0_dp
             kapP_cell(i)      = 0.0_dp
             Tspec_cell(i)     = 0.0_dp
             Qrad_cell(i)      = 0.0_dp
             Qvis_cell(i)      = 0.0_dp
             Qirr_cell(i)      = 0.0_dp
          else
            nu_conv(it, i) = nu(i)
            cs2            = 1.5_dp * nu_dim_loc * omegaK_cgs / alphaSS
            H(it, i)       = sqrt(cs2) / omegaK_cgs
            rho(it, i)     = sigma_cgs / (2.0_dp * H(it, i))
            Tmid(it, i)    = mu * mp * cs2 / kb
            T_emit        = Tmid(it, i)
            !write (*, '(a, 1pe121.3, a, 1pe11.3, a, 1pe11.3, a, 1pe11.3)') &
            !   'cs2=', cs2, ', H=', H(it, i), ', rho=', rho(it, i), &
            !   ', Tmid=', Tmid(it, i)
          end if
          kappaR(it, i) = 0.0_dp
          tauR(it, i) = 0.0_dp
          kappa_out = 0.0_dp
          tau_out   = 0.0_dp

          ! ---- Cache SED inputs for reuse (cell-centered)
          Sigma_cgs_cell(i) = Sigma_cgs
          kapP_cell(i)      = kappa_planck(it, i)   ! absorption-only Planck mean assumed
          Tspec_cell(i)     = T_emit                ! SED-ready temperature definition
          Qrad_cell(i) = Qrad(it, i)                ! Qrad from both faces (zero in this case)
          Qvis_cell(i) = Qvis(it, i)
          Qirr_cell(i) = Qirr(it, i)
       end if
 
       write(iu_st,'(1p,15e13.5)') r(i), r_cgs, sigma_cgs, H(it, i), &
             H(it, i)/r_cgs, dYdXi(it, i), rho(it, i), &
             nu_conv(it, i) * nu0_dim / nu0_nd, &
             Qvis_out, Qirr(it, i), Qrad_out, &
             Tmid(it, i), T_emit, kappa_out, tau_out
    end do

    close(iu_st)

    if (use_energy_balance) then
       call output_sed_dtl_from_cells(it, t, Sigma_cgs_cell, Tspec_cell, kapP_cell, &
                                      Qrad_cell, Qvis_cell, Qirr_cell)
    end if

  end subroutine output_disk_structure_dtl

  subroutine output_sed_dtl_from_cells(it, t, Sigma_cgs_cell, Tspec_cell, kapP_cell, &
                                       Qrad_cell, Qvis_cell, Qirr_cell)
   !-----------------------------------------------------------------------
   ! Output CBD SED snapshot synchronized with disk_t*.dat.
   !
   ! Input arrays are cell-centered values cached in output_disk_structure_dtl:
   !   Sigma_cgs_cell(i) : surface density [g/cm^2]
   !   Tspec_cell(i)     : SED-ready temperature [K]
   !                       (thick: Teff from Qrad; thin: Tmid)
   !   kapP_cell(i)      : Planck mean absorption opacity [cm^2/g]
   !
   ! Gray slab SED:
   !   dLnu = (2*pi*r*dr) * (2*pi*Bnu(T)) * (1 - exp(-tauP))
   !   tauP = 0.5 * kapP * Sigma   (one-side)
   !
   ! Output columns:
   !   nu [Hz], L_nu [erg/s/Hz], nu*L_nu [erg/s]
   !
   ! Parameters are read via run_control.nml through run_control_mod:
   !   sed_nnu, sed_nu_min_fac, sed_nu_max_fac
   !-----------------------------------------------------------------------
   use kind_params, only : i4b, dp
   use constants,   only : pi, hpl, kb, cc
   use mod_global,  only : nr, r, M_star, R0
   use run_control_mod, only : sed_nnu, sed_nu_min_fac, sed_nu_max_fac
   use mod_global,  only : mdot_inner_phys
   use star_params_mod, only : model_type
   use irradiation_mod, only : L_irr
   use units_disk_mod, only : t_dim, r_dim
   implicit none
 
   integer(i4b), intent(in) :: it
   real(dp),     intent(in) :: t
   real(dp),     intent(in) :: Sigma_cgs_cell(nr), Tspec_cell(nr), kapP_cell(nr)
   real(dp),     intent(in) :: Qrad_cell(nr), Qvis_cell(nr), Qirr_cell(nr)

   integer(i4b), save :: iu_sed = -1
   integer(i4b) :: i, j, nnu
   character(len=64) :: fname_sed
   real(dp) :: t_cgs
   real(dp), allocatable :: nu_grid(:), Lnu(:), nuLnu(:)
   real(dp), allocatable :: Lnu_vis(:), Lnu_irr(:), Lnu_tr(:)
   real(dp), allocatable :: nuLnu_vis(:), nuLnu_irr(:), nuLnu_tr(:)
   real(dp) :: Tmin, Tmax, nu_min, nu_max, lmin, lmax, dl
   real(dp) :: r_i, r_ip1, dr_cgs, rmid_cgs, dA
   real(dp) :: tauP_sed, fac, Bnu
   real(dp) :: Lbol_sed, Lbol_qrad, Lbol_vis, Lbol_irr, Lbol_tr, relerr, ltot_err
   real(dp) :: dnu, dLnu
   real(dp) :: qradp, f_vis, f_irr, f_tr
   real(dp), parameter :: tiny_q = 1.0e-99_dp

   t_cgs = t_dim(t)
 
   ! ---- Validate controls
   nnu = max(10_i4b, sed_nnu)
 
   allocate(nu_grid(nnu), Lnu(nnu), nuLnu(nnu))
   allocate(Lnu_vis(nnu), Lnu_irr(nnu), Lnu_tr(nnu))
   allocate(nuLnu_vis(nnu), nuLnu_irr(nnu), nuLnu_tr(nnu))
   Lnu_vis(:)   = 0.0_dp
   Lnu_irr(:)   = 0.0_dp
   Lnu_tr(:)    = 0.0_dp

   !--------------------------------------------
   ! 1) Determine frequency range from Tspec(min/max)
   !--------------------------------------------
   Tmin = huge(1.0_dp)
   Tmax = 0.0_dp
   do i = 1, nr
     if (Sigma_cgs_cell(i) <= 0.0_dp) cycle
     if (Tspec_cell(i) <= 0.0_dp) cycle
     Tmin = min(Tmin, Tspec_cell(i))
     Tmax = max(Tmax, Tspec_cell(i))
   end do
 
   if (Tmax <= 0.0_dp .or. Tmin >= huge(1.0_dp)) then
     ! Empty disk: write an empty SED safely
     Tmin = 10.0_dp
     Tmax = 10.0_dp
   end if
 
   ! Guard against pathological user inputs
   if (sed_nu_min_fac <= 0.0_dp) then
     nu_min = 0.1_dp * kb*Tmin / hpl
   else
     nu_min = sed_nu_min_fac * kb*Tmin / hpl
   end if
   if (sed_nu_max_fac <= 0.0_dp) then
     nu_max = 10.0_dp * kb*Tmax / hpl
   else
     nu_max = sed_nu_max_fac * kb*Tmax / hpl
   end if
 
   if (nu_max <= nu_min) then
     nu_max = 10.0_dp * nu_min
   end if
 
   !--------------------------------------------
   ! 2) Build log-spaced frequency grid
   !--------------------------------------------
   lmin = log(nu_min)
   lmax = log(nu_max)
   dl   = (lmax - lmin) / real(nnu-1, dp)
   do j = 1, nnu
     nu_grid(j) = exp(lmin + dl*real(j-1, dp))
   end do
 
   !--------------------------------------------
   ! 3) Annulus integration
   !--------------------------------------------
   Lnu(:) = 0.0_dp
   Lbol_sed  = 0.0_dp
   Lbol_qrad = 0.0_dp
   do i = 1, nr-1
     if (Sigma_cgs_cell(i) <= 0.0_dp) cycle
     if (Tspec_cell(i) <= 0.0_dp) cycle
 
     r_i   = r_dim(r(i))
     r_ip1 = r_dim(r(i+1))
     dr_cgs = r_ip1 - r_i
     if (dr_cgs <= 0.0_dp) cycle
 
     rmid_cgs = 0.5_dp*(r_ip1 + r_i)
     dA       = 2.0_dp*pi*rmid_cgs*dr_cgs
 
     ! Gray absorption optical depth (one-side)
     tauP_sed = 0.5_dp * max(kapP_cell(i), 0.0_dp) * max(Sigma_cgs_cell(i), 0.0_dp)
 
     if (tauP_sed > 50.0_dp) then
       fac = 1.0_dp
     else
       fac = 1.0_dp - exp(-max(tauP_sed, 0.0_dp))
     end if

     qradp = max(Qrad_cell(i), tiny_q)

     ! Fractions based on local heating / radiative loss
     f_vis = max(Qvis_cell(i), 0.0_dp) / qradp
     f_irr = max(Qirr_cell(i), 0.0_dp) / qradp

     ! Clamp to [0,1] to avoid pathological overshoots
     f_vis = min(max(f_vis, 0.0_dp), 1.0_dp)
     f_irr = min(max(f_irr, 0.0_dp), 1.0_dp)

     ! Residual/trans component (radial transport, time-dependence, numerical residuals, etc.)
     f_tr  = 1.0_dp - (f_vis + f_irr)
     if (f_tr < 0.0_dp) then
        if (f_vis + f_irr > 0.0_dp) then
           f_vis = f_vis / (f_vis + f_irr)
           f_irr = 1.0_dp - f_vis
        else
           f_vis = 0.0_dp
           f_irr = 0.0_dp
        end if
        f_tr = 0.0_dp
     end if

     do j = 1, nnu
       Bnu  = Bnu_planck(nu_grid(j), Tspec_cell(i))
    
       ! Both faces: Fnu = 2*pi*Bnu*(1-exp(-tau))
       dLnu = dA * (2.0_dp*pi) * Bnu * fac
    
       Lnu(j)     = Lnu(j)     + dLnu
       Lnu_vis(j) = Lnu_vis(j) + dLnu * f_vis
       Lnu_irr(j) = Lnu_irr(j) + dLnu * f_irr
       Lnu_tr(j)  = Lnu_tr(j)  + dLnu * f_tr
     end do

     ! ---- Bolometric luminosity from Qrad: L = ∫ (2*pi*r*dr) * Qrad
     Lbol_qrad = Lbol_qrad + dA * max(Qrad_cell(i), 0.0_dp)

     if (i == 1) then
       write(*,'(a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4)') 'DEBUG f_vis=', f_vis, ' f_irr=', f_irr, ' f_tr=', f_tr, ' qrad=', qradp
     end if
   end do
 
   do j = 1, nnu
     nuLnu(j)     = nu_grid(j) * Lnu(j)
     nuLnu_vis(j) = nu_grid(j) * Lnu_vis(j)
     nuLnu_irr(j) = nu_grid(j) * Lnu_irr(j)
     nuLnu_tr(j)  = nu_grid(j) * Lnu_tr(j)
   end do

   ! ---- Bolometric luminosity from SED: integrate over log nu for robustness
   ! Lbol = ∫ Lnu dnu = ∫ (nu*Lnu) d(ln nu)
   Lbol_sed = 0.0_dp
   Lbol_vis = 0.0_dp
   Lbol_irr = 0.0_dp
   Lbol_tr  = 0.0_dp
   do j = 1, nnu-1
      dnu = log(nu_grid(j+1)) - log(nu_grid(j))
      Lbol_sed = Lbol_sed + 0.5_dp * (nu_grid(j)*Lnu(j) + nu_grid(j+1)*Lnu(j+1)) * dnu
      Lbol_vis = Lbol_vis + 0.5_dp * (nu_grid(j)*Lnu_vis(j) + nu_grid(j+1)*Lnu_vis(j+1)) * dnu
      Lbol_irr = Lbol_irr + 0.5_dp * (nu_grid(j)*Lnu_irr(j) + nu_grid(j+1)*Lnu_irr(j+1)) * dnu
      Lbol_tr  = Lbol_tr + 0.5_dp * (nu_grid(j)*Lnu_tr(j) + nu_grid(j+1)*Lnu_tr(j+1)) * dnu
   end do

   ! ---- Relative error between SED and Qrad bolometric luminosities
   if (Lbol_qrad > 0.0_dp) then
      relerr = (Lbol_sed - Lbol_qrad) / Lbol_qrad
   else if (Lbol_sed > 0.0_dp) then
      relerr = 1.0_dp
   else
      relerr = 0.0_dp
   end if

   ! ---- Relative error between SED and Qvis bolometric luminosities
   if (Lbol_sed > 0.0_dp) then
     ltot_err = (Lbol_vis + Lbol_irr + Lbol_tr - Lbol_sed) / Lbol_sed
   else
     ltot_err = 0.0_dp
   end if

   !--------------------------------------------
   ! 4) Write file
   !--------------------------------------------
   write(fname_sed,'("sed_t",i8.8,".dat")') it
   open(newunit=iu_sed, file=fname_sed, status='replace', action='write')
 
   write(iu_sed,'(a, i8, a, 1pe13.5, a)') '# it =', it, ', t =', t_cgs, '[s]'
   if (model_type == 'bh') then
     write(iu_sed,'(a,1pe13.5,a,1pe13.5,a,1pe13.5,a,1pe13.5,a)') &
       '# M_star =', M_star, '[g], R0 =', R0, '[cm], Mdot(R0) =', &
       mdot_inner_phys, '[g/s], L_irr =', L_irr, '[erg/s]'
   else
     write(iu_sed,'(a,1pe13.5,a,1pe13.5,a)') '# M_star =', M_star, '[g], R0 =', R0, '[cm]'
   end if
   write(iu_sed,'(a, i0, a, 1pe13.5, a, 1pe13.5)') '# sed_nnu = ', nnu, &
        ', sed_nu_min_fac =', sed_nu_min_fac, ', sed_nu_max_fac =', sed_nu_max_fac
   write(iu_sed,'(a, 1pe13.5, a, 1pe13.5, a)') '# nu_min = ', nu_min, &
        ' [Hz], nu_max =', nu_max, ' [Hz]'
   write(iu_sed,'(a,1pe13.5,a,1pe13.5,a,1pe13.5)') &
        '# Lbol = ', Lbol_sed, ' [erg/s], Integration of Qrad = ', Lbol_qrad, &
        ' [erg/s] -> error(relative) = ', relerr
   write(iu_sed,'(a,1pe11.3,a,1pe11.3,a,1pe11.3,a, 1pe11.3)') &
        '# Lbol_vis =', Lbol_vis, ' [erg/s], Lbol_irr =', Lbol_irr, &
        ' [erg/s], Lbol_tr =', Lbol_tr, &
        ' [erg/s] vs. Lbol -> error(relative) =', ltot_err
   write(iu_sed,'(a, 1pe13.5, a, 1pe13.5, a, 1pe13.5)') &
        '# Fractions (bolometric): Lbol_vis/Lbol =', Lbol_vis/Lbol_sed, &
        ', Lbol_irr/Lbol =', Lbol_irr/Lbol_sed, &
        ', Lbol_tr/Lbol =', Lbol_tr/Lbol_sed
   write(iu_sed,'(a)') '#'
   write(iu_sed,'(a, 5x, a, 9x, a, 5x, a, 5x, a, 5x, a, 5x, a, 5x, a, 5x, a, 6x, a)') '#', &
         'nu', 'Lnu_tot', 'nuLnu_tot', 'Lnu_vis', 'nuLnu_vis', 'Lnu_irr', 'nuLnu_irr', &
         'Lnu_tr', 'nuLnu_tr'
   write(iu_sed,'(a, 4x, a, 7x, a, 4x, a, 5x, a, 4x, a, 5x, a, 4x, a, 5x, a, 4x, a)') '#', &
         '[Hz]', '[erg/s/Hz]', '[erg/s]', '[erg/s/Hz]', '[erg/s]', '[erg/s/Hz]', &
         '[erg/s]', '[erg/s/Hz]', '[erg/s]'
   do j = 1, nnu
     write(iu_sed,'(1p9e13.5)') nu_grid(j), Lnu(j), nuLnu(j), &
           Lnu_vis(j), nuLnu_vis(j), Lnu_irr(j), nuLnu_irr(j), Lnu_tr(j), nuLnu_tr(j)
   end do
   close(iu_sed)
 
   deallocate(nu_grid, Lnu, nuLnu)
 
 contains
 
   pure function Bnu_planck(nu, T) result(B)
     use kind_params, only : dp
     use constants,   only : hpl, kb, cc
     real(dp), intent(in) :: nu, T
     real(dp) :: B, x, ex
 
     if (T <= 0.0_dp .or. nu <= 0.0_dp) then
       B = 0.0_dp
       return
     end if
 
     x = hpl*nu/(kb*T)
 
     if (x > 700.0_dp) then
       B = 0.0_dp
     else
       ex = exp(x)
       B  = (2.0_dp*hpl*nu**3 / cc**2) / (ex - 1.0_dp)
     end if
   end function Bnu_planck
 
  end subroutine output_sed_dtl_from_cells
   
end module output_mod

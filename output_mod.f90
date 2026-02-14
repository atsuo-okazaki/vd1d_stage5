!==============================================================
! output_mode.f90
!   Write disk surface-density evolution to text files
!   (for Gnuplot or other plotting tools)
!==============================================================
module output_mod
  use kind_params,  only : dp, i4b, lgt
  use constants,    only : pi, kb, mp, mu, msun, year, sbc
  use mod_global,   only : nr, r, nu, sigmat, Tmid,     &
                           H, rho, kappaR, tauR, dt, fwhm,    &
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
  use disk_thermal_mod, only : heating_cooling_cell, Qvis_and_Qrad
  use opacity_table_mod, only : get_opacity_Planck_rhoT

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

  subroutine output_disk_structure_dtl(it, t)
    use mod_global, only : mdot_inner_phys
    use irradiation_mod, only : L_irr
    integer(i4b), intent(in) :: it
    real(dp),    intent(in)  :: t

    integer(i4b), save :: iu_st = -1
    integer(i4b) :: i, ierr_planck
    character(len=64) :: fname
    real(dp) :: t_cgs, r_cgs, sigma_cgs, omegaK_cgs
    real(dp) :: nu_dim_loc, T_loc
    real(dp) :: kappaR_loc, tauR_loc, kappaP_loc, Qvis_loc, Qrad_loc
    real(dp) :: Tsurf
    real(dp) :: cs2
    real(dp), parameter :: tiny_sigma = 1.0e-8_dp

    t_cgs = t_dim(t)
    
    write(fname,'("disk_t",i8.8,".dat")') it
    open(newunit=iu_st, file=fname, status='replace')

    if (k_iter(it) > 0) then
       if (use_irradiation) then
          if (m_iter(it) > 0) then
             write(iu_st,'(a, i8, a, 1pe13.5, a, i3, a, i3, a)') '# it =', it, &
                  ', t =', t_cgs, '[s] (converged after', k_iter(it), &
                  ' iterations; consistent Tmid and Qirr after', &
                  m_iter(it), ' iterations)'
          else
             write(iu_st,'(a, i8, a, 1pe13.5, a, i3, a)') '# it =', it, &
                  ', t =', t_cgs, '[s] (converged after', k_iter(it), &
                  ' iterations; no consistent Tmid and Qirr)'
          end if
       else
          write(iu_st,'(a, i8, a, 1pe13.5, a, i3, a)') '# it =', it, &
               ', t =', t_cgs, '[s] (converged after', k_iter(it), ' iterations)'
       end if
    else if (k_iter(it) == 0) then
       if (use_irradiation) then
          if (m_iter(it) > 0) then
             write(iu_st,'(a, i8, a, 1pe13.5, a, i3, a, i3, a)') '# it =', it, &
                  ', t =', t_cgs, '[s] (Consistent Tmid and Qirr after', &
                  m_iter(it), ' iterations)'
          else
             write(iu_st,'(a, i8, a, 1pe13.5, a, i3, a)') '# it =', it, &
                  ', t =', t_cgs, '[s] (No consistent Tmid and Qirr)'
          end if
       else
          write(iu_st,'(a, i8, a, 1pe13.5, a, i3, a)') '# it =', it, &
               ', t =', t_cgs, '[s]'
       end if
    else
       if (use_irradiation) then
          if (m_iter(it) > 0) then
             write(iu_st,'(a, i8, a, 1pe13.5, a, i3, a)') '# it =', it, &
                  ', t =', t_cgs, ' [s] (not converged; consistent Tmid and Qirr after', &
                  m_iter(it), ' iterations)'
          else
             write(iu_st,'(a, i8, a, 1pe13.5, a)') '# it =', it, &
                  ', t =', t_cgs, '[s] (not converged; no consistent Tmid and Qirr)'
          end if
       else
          write(iu_st,'(a, i8, a, 1pe13.5, a)') '# it =', it, &
               ', t =', t_cgs, '[s] (not converged)'
      end if
    end if
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
    write(iu_st,'(a, 2x, a, 6x, a, 4x, a, 4x, a, 7x, a, 7x, a, 4x, a, 2x, a, 6x, &
         & a, 9x, a, 9x, a, 8x, a, 5x, a, 3x, a, 3x, a)') '#', 'xi(=r/R0)', 'r(cm)', &
         & 'sigma(g/cm^2)', 'H(cm)', 'Y(=H/r)', 'dY/dxi', &
         & 'rho(g/cm^3)', 'nu(cm^2/s)', &
         & 'Qvis', 'Qirr', 'Qrad', 'Tmid(K)', &
         & 'Tsurf(K)', 'kappaR(cm^2/g)', 'tauR'

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
          
             write(iu_st,'(1p,15e13.5)') r(i), r_cgs, sigma_cgs, H(it,i), &
                  H(it,i)/r_cgs, 0.0_dp, rho(it,i), nu_dim_loc, &
                  0.0_dp, 0.0_dp, 0.0_dp, &  ! Qvis, Qirr, Qrad
                  0.0_dp, 0.0_dp, 0.0_dp, &  ! Tmid, Tsurf, kappaR
                  0.0_dp                     ! tauR
             cycle
          end if

          call get_opacity_Planck_rhoT(rho(it, i), T_loc, kappaP_loc, ierr_planck)
          if (ierr_planck /= 0) kappaP_loc = kappaR_loc
          call Qvis_and_Qrad(r_cgs, Sigma_cgs, OmegaK_cgs, nu_dim_loc, T_loc, &
                             kappaR_loc, tauR_loc, Qvis_loc, Qrad_loc, kappaP_loc=kappaP_loc)
   
          !! Qrad_loc; radiative flux from both surfaces
          !! Qrad_loc = 2 * sbc *Tsurf**4
          Tsurf = (Qrad_loc / (2.0_dp * sbc))**0.25_dp
       else
          ! Isotherma-disk quantities
          if (sigma_cgs <= tiny_sigma) then
             nu_conv(it, i) = nu(i)
             H(it, i)     = 0.0_dp
             rho(it, i)   = 0.0_dp
             Tmid(it, i)  = 0.0_dp
             Tsurf        = 0.0_dp
          else
            nu_conv(it, i) = nu(i)
            cs2            = 1.5_dp * nu_dim_loc * omegaK_cgs / alphaSS
            H(it, i)       = sqrt(cs2) / omegaK_cgs
            rho(it, i)     = sigma_cgs / (2.0_dp * H(it, i))
            Tmid(it, i)    = mu * mp * cs2 / kb
            Tsurf          = Tmid(it, i)
            !write (*, '(a, 1pe121.3, a, 1pe11.3, a, 1pe11.3, a, 1pe11.3)') &
            !   'cs2=', cs2, ', H=', H(it, i), ', rho=', rho(it, i), &
            !   ', Tmid=', Tmid(it, i)
          end if
          Qvis_loc = 0.0_dp
          Qrad_loc = 0.0_dp
          kappaR(it, i) = 0.0_dp
          tauR(it, i) = 0.0_dp
       end if
 
       write(iu_st,'(1p,15e13.5)') r(i), r_cgs, sigma_cgs, H(it, i), &
             H(it, i)/r_cgs, dYdXi(it, i), rho(it, i), &
             nu_conv(it, i) * nu0_dim / nu0_nd, &
             Qvis_loc, Qirr(it, i), Qrad_loc, &
             Tmid(it, i), Tsurf, kappaR(it, i), tauR(it, i)
    end do

    close(iu_st)

  end subroutine output_disk_structure_dtl

  subroutine output_disk_structure(it, t)
    use mod_global, only : mdot_inner_phys
    use irradiation_mod, only : L_irr
    integer(i4b), intent(in) :: it
    real(dp),    intent(in)  :: t

    integer(i4b), save :: iu_dtl = -1
    integer(i4b) :: i, ierr_planck
    character(len=64) :: fname
    real(dp) :: t_cgs, r_cgs, sigma_cgs, omegaK_cgs
    real(dp) :: nu_dim_loc, T_loc
    real(dp) :: kappaR_loc, tauR_loc, kappaP_loc, Qvis_loc, Qrad_loc
    real(dp) :: Tsurf
    real(dp) :: cs2
    real(dp), parameter :: tiny_sigma = 1.0e-8_dp

    t_cgs = t_dim(t)
    
    write(fname,'("disk_t",i8.8,".dat")') it
    open(newunit=iu_dtl, file=fname, status='replace')

    if (k_iter(it) > 0) then
       if (use_irradiation) then
          if (m_iter(it) > 0) then
             write(iu_dtl,'(a, i8, a, 1pe13.5, a, i3, a, i3, a)') '# it =', it, &
                  ', t =', t_cgs, '[s] (converged after', k_iter(it), &
                  ' iterations; consistent Tmid and Qs after', &
                  m_iter(it), ' iterations)'
          else
             write(iu_dtl,'(a, i8, a, 1pe13.5, a, i3, a)') '# it =', it, &
                  ', t =', t_cgs, '[s] (converged after', k_iter(it), &
                  ' iterations; no consistent Tmid and Qs)'
          end if
       else
          write(iu_dtl,'(a, i8, a, 1pe13.5, a, i3, a)') '# it =', it, &
               ', t =', t_cgs, '[s] (converged after', k_iter(it), ' iterations)'
       end if
    else if (k_iter(it) == 0) then
       write(iu_dtl,'(a, i8, a, 1pe13.5, a, i3, a)') '# it =', it, &
            ', t =', t_cgs, '[s]'
    else
       if (use_irradiation) then
          if (m_iter(it) > 0) then
             write(iu_dtl,'(a, i8, a, 1pe13.5, a, i3, a)') '# it =', it, &
                  ', t =', t_cgs, ' [s] (not converged; consistent Tmid and Qs after', &
                  m_iter(it), ' iterations)'
          else
             write(iu_dtl,'(a, i8, a, 1pe13.5, a)') '# it =', it, &
                  ', t =', t_cgs, '[s] (not converged; no consistent Tmid and Qs)'
          end if
       else
          write(iu_dtl,'(a, i8, a, 1pe13.5, a)') '# it =', it, &
               ', t =', t_cgs, '[s] (not converged)'
      end if
    end if
    if (model_type == 'star') then
       write (iu_dtl, '(a, 1pe13.5, a, 1pe13.5, a, 1pe13.5, a, 1pe13.5, a)') &
          '# M_star =', M_star, '[g], R_star =', R_star, &
          '[cm], Teff_star =', Teff_star, '[K], R0 =', R0, '[cm]'
    else
       ! mdoel_type = 'bh'
       write (iu_dtl, '(a, 1pe13.5, a, 1pe13.5, a, 1pe13.5, a, 1pe13.5, a)') &
          '# M_star =', M_star, '[g], R0 =', R0, '[cm], Mdot(R0) =', &
          mdot_inner_phys, '[g/s], L_irr =', L_irr, '[erg/s]'
    end if
    write(iu_dtl,'(a, 4x, a, 4x, a, 4x, a, 5x, a, 2x, a, 6x, &
         & a, 9x, a, 9x, a, 8x, a, 5x, a, 3x, a, 3x, a)') '#', 'r(cm)', &
         & 'sigma(g/cm^2)', 'H(cm)', 'rho(g/cm^3)', 'nu(cm^2/s)', &
         & 'Qvis', 'Qirr', 'Qrad', 'Tmid(K)', &
         & 'Tsurf(K)', 'kappaR(cm^2/g)', 'tauR'

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
          
             write(iu_dtl,'(1p,13e13.5)') r_cgs, sigma_cgs, H(it,i), &
                  rho(it,i), nu_dim_loc, &
                  0.0_dp, 0.0_dp, 0.0_dp, &  ! Qvis, Qirr, Qrad
                  0.0_dp, 0.0_dp, 0.0_dp, &  ! Tmid, Tsurf, kappaR
                  0.0_dp                     ! tauR
             cycle
          end if

          call get_opacity_Planck_rhoT(rho(it, i), T_loc, kappaP_loc, ierr_planck)
          if (ierr_planck /= 0) kappaP_loc = kappaR_loc
          call Qvis_and_Qrad(r_cgs, Sigma_cgs, OmegaK_cgs, nu_dim_loc, T_loc, &
                             kappaR_loc, tauR_loc, Qvis_loc, Qrad_loc, kappaP_loc=kappaP_loc)
   
          !! Qrad_loc; radiative flux from both surfaces
          !! Qrad_loc = 2 * sbc *Tsurf**4
          Tsurf = (Qrad_loc / (2.0_dp * sbc))**0.25_dp
       else
          ! Isotherma-disk quantities
          if (sigma_cgs <= tiny_sigma) then
             nu_conv(it, i) = nu(i)
             H(it, i)     = 0.0_dp
             rho(it, i)   = 0.0_dp
             Tmid(it, i)  = 0.0_dp
             Tsurf        = 0.0_dp
          else
            nu_conv(it, i) = nu(i)
            cs2            = 1.5_dp * nu_dim_loc * omegaK_cgs / alphaSS
            H(it, i)       = sqrt(cs2) / omegaK_cgs
            rho(it, i)     = sigma_cgs / (2.0_dp * H(it, i))
            Tmid(it, i)    = mu * mp * cs2 / kb
            Tsurf          = Tmid(it, i)
            !write (*, '(a, 1pe121.3, a, 1pe11.3, a, 1pe11.3, a, 1pe11.3)') &
            !   'cs2=', cs2, ', H=', H(it, i), ', rho=', rho(it, i), &
            !   ', Tmid=', Tmid(it, i)
          end if
          Qvis_loc = 0.0_dp
          Qrad_loc = 0.0_dp
          kappaR(it, i) = 0.0_dp
          tauR(it, i) = 0.0_dp
       end if
 
       write(iu_dtl,'(1p,12e13.5)') r_cgs, sigma_cgs, H(it, i), &
             rho(it, i), nu_conv(it, i) * nu0_dim / nu0_nd, &
             Qvis_loc, Qirr(it, i), Qrad_loc, &
             Tmid(it, i), Tsurf, kappaR(it, i), tauR(it, i)
    end do

    close(iu_dtl)

  end subroutine output_disk_structure
  
end module output_mod

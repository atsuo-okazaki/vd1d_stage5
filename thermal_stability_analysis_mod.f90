!==============================================================
! thermal_stability_analysis_mod.f90
!
! Thermal instability analysis for accretion disks.
!
! Background: Around Tmid ~ 10^4 K, opacity tables (OP-Ferguson,
! Ferguson-Semenov transitions) change abruptly. This can create
! thermally unstable regions where dQ-/dT < dQ+/dT at equilibrium.
!
! Thermal stability criterion:
!   - Stable:   dQminus/dT > dQplus/dT  (perturbation damps)
!   - Unstable: dQminus/dT < dQplus/dT  (perturbation grows)
!
! This module provides:
!   - compute_thermal_derivatives: dQ+/dT, dQ-/dT via finite difference
!   - thermal_stability_flag: stable (1) or unstable (-1)
!   - output_thermal_stability_profile: append to disk structure output
!   - output_thermal_stability_all_roots: find ALL equilibrium solutions at each r
!   - generate_scurve_at_radius: S-curve (Sigma vs T) at fixed r
!==============================================================
module thermal_stability_analysis_mod
  use kind_params,  only : dp, i4b
  use constants,    only : kb, mp, mu
  use mod_global,   only : t0, R0
  use disk_thermal_mod, only : heating_cooling_cell
  use units_disk_mod, only : r_dim, sigma_dim, omegaK_dim, t_dim
  use radiation_params_mod, only : T_floor, T_ceiling
  use opacity_table_mod, only : thigh_tab, tlow_tab, logRmin_tab, rhomin_tab
  use ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
  implicit none

  ! Relative perturbation for finite-difference derivatives
  real(dp), parameter :: dlogT_pert = 1.0e-4_dp
  ! Max number of thermal equilibrium roots per (r, Sigma)
  integer(i4b), parameter :: max_roots = 20

contains

  !------------------------------------------------------------
  ! compute_thermal_derivatives
  !
  ! Compute d(Qplus)/dT and d(Qminus)/dT at given (r, Sigma, T).
  ! Qplus = Qvis + Qirr,  Qminus = Qrad.
  ! Uses central finite difference in log(T).
  !
  ! Output:
  !   dQplus_dT, dQminus_dT : [erg/cm^2/s/K]
  !   thermally_stable      : .true. if dQminus_dT > dQplus_dT
  !------------------------------------------------------------
  subroutine compute_thermal_derivatives(r_cgs, Sigma_cgs, OmegaK_cgs, &
       shadow, Tmid_val, Qirr_in, dQplus_dT, dQminus_dT, thermally_stable)
    real(dp), intent(in)  :: r_cgs, Sigma_cgs, OmegaK_cgs
    logical,  intent(in)  :: shadow
    real(dp), intent(in)  :: Tmid_val, Qirr_in
    real(dp), intent(out) :: dQplus_dT, dQminus_dT
    logical,  intent(out) :: thermally_stable

    real(dp) :: T_lo, T_hi
    real(dp) :: H_tmp, rho_tmp, nu_tmp, kappa_tmp, kappa_planck_tmp, tau_tmp
    real(dp) :: Qvis_lo, Qirr_lo, Qrad_lo
    real(dp) :: Qvis_hi, Qirr_hi, Qrad_hi
    real(dp) :: Qplus_lo, Qminus_lo, Qplus_hi, Qminus_hi
    real(dp) :: dT
    real(dp), parameter :: T_tiny = 1.0e-99_dp

    T_lo = Tmid_val * (10.0_dp**(-dlogT_pert))
    T_hi = Tmid_val * (10.0_dp**(+dlogT_pert))

    call heating_cooling_cell(r_cgs, Sigma_cgs, OmegaK_cgs, shadow, T_lo, &
         H_tmp, rho_tmp, nu_tmp, kappa_tmp, kappa_planck_tmp, tau_tmp, &
         Qvis_lo, Qirr_lo, Qrad_lo, Qirr_in=Qirr_in)
    Qplus_lo = Qvis_lo + Qirr_lo
    Qminus_lo = Qrad_lo

    call heating_cooling_cell(r_cgs, Sigma_cgs, OmegaK_cgs, shadow, T_hi, &
         H_tmp, rho_tmp, nu_tmp, kappa_tmp, kappa_planck_tmp, tau_tmp, &
         Qvis_hi, Qirr_hi, Qrad_hi, Qirr_in=Qirr_in)
    Qplus_hi = Qvis_hi + Qirr_hi
    Qminus_hi = Qrad_hi

    dT = T_hi - T_lo
    if (abs(dT) < T_tiny) then
       dQplus_dT  = 0.0_dp
       dQminus_dT = 0.0_dp
       thermally_stable = .true.
       return
    end if

    dQplus_dT  = (Qplus_hi  - Qplus_lo)  / dT
    dQminus_dT = (Qminus_hi - Qminus_lo) / dT
    thermally_stable = (dQminus_dT > dQplus_dT)

  end subroutine compute_thermal_derivatives

  !------------------------------------------------------------
  ! thermal_stability_flag
  !
  ! Returns +1 (stable) or -1 (unstable) for plotting.
  !------------------------------------------------------------
  integer(i4b) function thermal_stability_flag(dQplus_dT, dQminus_dT) result(flag)
    real(dp), intent(in) :: dQplus_dT, dQminus_dT
    if (dQminus_dT > dQplus_dT) then
       flag = 1   ! stable
    else
       flag = -1  ! unstable
    end if
  end function thermal_stability_flag

  !------------------------------------------------------------
  ! output_thermal_stability_profile
  !
  ! Append thermal stability diagnostics to disk structure output.
  ! Call this from output_mod when use_energy_balance is true.
  !
  ! Columns appended: dQplus_dT  dQminus_dT  stability_flag
  !   stability_flag: 1 = stable, -1 = unstable
  !------------------------------------------------------------
  subroutine output_thermal_stability_profile(nr_out, r_nd, sigma_nd, &
       Tmid_arr, Qirr_arr, is_shadow_arr, it, t_nd, fname_base)
    integer(i4b), intent(in) :: nr_out
    real(dp),     intent(in) :: r_nd(nr_out), sigma_nd(nr_out)
    real(dp),     intent(in) :: Tmid_arr(nr_out), Qirr_arr(nr_out)
    logical,      intent(in) :: is_shadow_arr(nr_out)
    integer(i4b), intent(in) :: it
    real(dp),     intent(in) :: t_nd
    character(len=*), intent(in), optional :: fname_base

    integer(i4b) :: iu, i
    real(dp) :: r_cgs, Sigma_cgs, OmegaK_cgs
    real(dp) :: dQplus_dT, dQminus_dT
    logical  :: thermally_stable
    character(len=64) :: fname
    real(dp), parameter :: tiny_sigma = 1.0e-8_dp
    real(dp), parameter :: tiny_T = 1.0_dp

    if (present(fname_base)) then
       write(fname, '(a,"_thermal_stability_",i8.8,".dat")') trim(fname_base), it
    else
       write(fname, '("thermal_stability_",i8.8,".dat")') it
    end if

    open(newunit=iu, file=trim(fname), status='replace', action='write')
    write(iu, '(a, i8, a, 1pe14.6, a, 1pe14.6, a)') '# it = ', it, ', t_nd = ', t_nd, ' (t = ', t_dim(t_nd), ' [s])'
    write(iu, '(a, 1pe14.6, a, 1pe14.6, a)') '# t0 = ', t0, ', R0 = ', R0, ' [cm]'
    write(iu, '(a)') '# Thermal instability analysis: dQ+/dT, dQ-/dT [erg/cm^2/s/K]'
    write(iu, '(a)') '# stability: 1=stable, -1=unstable, 0=off-equilibrium(|R| too large)'
    write(iu, '(a)') '# resR = |(Qvis+Qirr)-Qrad| / max(Qrad,tiny)'
    write(iu, '(a)') '# near_turning: 1 if |(dF/dlogSigma)/Qrad| is small at this (Sigma,T)'
    write(iu, '(a)') '#    r_cgs        Sigma        Tmid         kappa        ' &
         &     //'Qvis         Qirr         Qrad       dQplus_dT   dQminus_dT ' &
         &     //'stability   resR  near_turning dF_dlogSigma_norm'

    block
      integer(i4b) :: n_unstable
      real(dp) :: r_unstable_min, r_unstable_max
      n_unstable = 0
      r_unstable_min = 1.0e99_dp
      r_unstable_max = -1.0_dp

    do i = 1, nr_out
       r_cgs     = r_dim(r_nd(i))
       Sigma_cgs = sigma_dim(sigma_nd(i))
       OmegaK_cgs = omegaK_dim(r_nd(i))

       if (Sigma_cgs <= tiny_sigma .or. Tmid_arr(i) <= tiny_T) then
         write(iu, '(1p,9e13.5, 3x, i3)') r_cgs, Sigma_cgs, 0.0_dp, 0.0_dp, &
               0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1
          cycle
       end if

       ! Evaluate Q's and opacity at the current (Sigma,Tmid) first
       block
         real(dp) :: H_tmp, rho_tmp, nu_tmp, tau_tmp
         real(dp) :: Qplus, Qminus, f, resR
         real(dp) :: kappa_tmp, kappa_planck_tmp, Qvis_tmp, Qirr_tmp, Qrad_tmp
         real(dp) :: dF_dlogSigma_norm
         integer(i4b) :: stab_flag, near_turning

         real(dp), parameter :: tiny_Q = 1.0e-99_dp
         real(dp), parameter :: resR_crit = 1.0e-3_dp   ! <-- tune if needed
         real(dp), parameter :: turn_crit = 1.0e-2_dp   ! <-- |dF/dlogSigma|/Qrad threshold

         call heating_cooling_cell(r_cgs, Sigma_cgs, OmegaK_cgs, is_shadow_arr(i), &
              Tmid_arr(i), H_tmp, rho_tmp, nu_tmp, kappa_tmp, kappa_planck_tmp, tau_tmp, &
              Qvis_tmp, Qirr_tmp, Qrad_tmp, Qirr_in=Qirr_arr(i))

         Qplus  = Qvis_tmp + Qirr_tmp
         Qminus = Qrad_tmp
         f      = Qplus - Qminus
         resR   = abs(f) / max(abs(Qminus), tiny_Q)

         ! Default outputs
         dQplus_dT  = 0.0_dp
         dQminus_dT = 0.0_dp
         stab_flag  = 0
         near_turning = 0
         dF_dlogSigma_norm = 0.0_dp

         if (resR <= resR_crit) then
            ! On/near equilibrium: stability is meaningful
            call compute_thermal_derivatives(r_cgs, Sigma_cgs, OmegaK_cgs, &
                 is_shadow_arr(i), Tmid_arr(i), Qirr_arr(i), &
                 dQplus_dT, dQminus_dT, thermally_stable)

            if (thermally_stable) then
               stab_flag = 1
            else
               stab_flag = -1
            end if

            ! Turning-point proximity (cheap local diagnostic):
            ! turning in T(Sigma) occurs near dF/dlogSigma ~ 0 at fixed T
            call turning_proximity_flag(r_cgs, Sigma_cgs, OmegaK_cgs, is_shadow_arr(i), &
                 Tmid_arr(i), Qirr_arr(i), dF_dlogSigma_norm, near_turning, turn_crit)
         else
            ! Off equilibrium: keep stab_flag=0, derivatives=0, near_turning=0
         end if

         if (stab_flag == -1) then
            n_unstable = n_unstable + 1
            r_unstable_min = min(r_unstable_min, r_cgs)
            r_unstable_max = max(r_unstable_max, r_cgs)
         end if

         write(iu, '(1p,9e13.5, 3x, i3, 1x, e13.5, 3x, i3, 6x, e13.5)') &
              r_cgs, Sigma_cgs, Tmid_arr(i), kappa_tmp, &
              Qvis_tmp, Qirr_tmp, Qrad_tmp, dQplus_dT, dQminus_dT, &
              stab_flag, resR, near_turning, dF_dlogSigma_norm
       end block
    end do

    if (n_unstable > 0) then
       write(iu, '(a, i0, a, 1pe14.6, a, 1pe14.6, a)') &
            '# Summary: n_unstable = ', n_unstable, &
            ', r_unstable [cm] = [', r_unstable_min, ', ', r_unstable_max, ']'
    else
       write(iu, '(a)') '# Summary: n_unstable = 0'
    end if
    end block
    close(iu)

  end subroutine output_thermal_stability_profile

  subroutine turning_proximity_flag(r_cgs, Sigma_cgs, OmegaK_cgs, shadow, Tmid, Qirr_in, &
                                   dF_dlogSigma_norm, near_turning, turn_crit)
    ! Estimate proximity to a turning point of the equilibrium manifold F(T,Sigma)=0
    ! using a local finite-difference of F with respect to log(Sigma) at fixed T.
    !
    ! At a turning point in T(Sigma), dT/dSigma diverges -> (∂F/∂Sigma)_T ~ 0.
    !
    ! We compute:
    !   F = (Qvis+Qirr) - Qrad
    !   dF/dlogSigma ~ (F_hi - F_lo) / (logSigma_hi - logSigma_lo)
    ! and normalize by Qrad for robustness.

    real(dp), intent(in)  :: r_cgs, Sigma_cgs, OmegaK_cgs
    logical,  intent(in)  :: shadow
    real(dp), intent(in)  :: Tmid, Qirr_in
    real(dp), intent(out) :: dF_dlogSigma_norm
    integer(i4b), intent(out) :: near_turning
    real(dp), intent(in)  :: turn_crit

    real(dp) :: Sig_lo, Sig_hi, logSig_lo, logSig_hi
    real(dp) :: H_t, rho_t, nu_t, kap_t, kapP_t, tau_t
    real(dp) :: Qv, Qi, Qr, F_lo, F_hi
    real(dp) :: Qrad_ref
    real(dp), parameter :: dlogSigma = 1.0e-3_dp  ! ~0.23% in Sigma
    real(dp), parameter :: tiny_Q = 1.0e-99_dp

    near_turning = 0
    dF_dlogSigma_norm = 0.0_dp

    Sig_lo = Sigma_cgs * (10.0_dp**(-dlogSigma))
    Sig_hi = Sigma_cgs * (10.0_dp**(+dlogSigma))

    logSig_lo = log10(max(Sig_lo, 1.0e-99_dp))
    logSig_hi = log10(max(Sig_hi, 1.0e-99_dp))

    call heating_cooling_cell(r_cgs, Sig_lo, OmegaK_cgs, shadow, Tmid, &
                H_t, rho_t, nu_t, kap_t, kapP_t, tau_t, Qv, Qi, Qr, Qirr_in=Qirr_in)
                F_lo = (Qv + Qi) - Qr

    call heating_cooling_cell(r_cgs, Sig_hi, OmegaK_cgs, shadow, Tmid, &
                H_t, rho_t, nu_t, kap_t, kapP_t, tau_t, Qv, Qi, Qr, Qirr_in=Qirr_in)

    F_hi = (Qv + Qi) - Qr
    Qrad_ref = max(abs(Qr), tiny_Q)

    dF_dlogSigma_norm = abs((F_hi - F_lo) / max((logSig_hi - logSig_lo), 1.0e-99_dp)) / Qrad_ref

    if (dF_dlogSigma_norm <= turn_crit) then
       near_turning = 1
    end if
  end subroutine turning_proximity_flag

  !------------------------------------------------------------
  ! find_all_thermal_roots
  !
  ! For fixed (r, Sigma, Qirr, shadow), find ALL T such that
  ! Qvis + Qirr = Qrad (thermal equilibrium).
  ! Scans T from T_floor to T_ceiling, detects sign changes in F(T),
  ! bisects each bracket. Returns T_roots(1:n_roots) sorted ascending.
  !------------------------------------------------------------
  subroutine find_all_thermal_roots(r_cgs, Sigma_cgs, OmegaK_cgs, shadow, &
       Qirr_in, T_roots, n_roots)
    real(dp), intent(in)  :: r_cgs, Sigma_cgs, OmegaK_cgs
    logical,  intent(in)  :: shadow
    real(dp), intent(in)  :: Qirr_in
    real(dp), intent(out) :: T_roots(max_roots)
    integer(i4b), intent(out) :: n_roots

    integer(i4b), parameter :: Nscan = 500
    integer(i4b), parameter :: iter_max = 50
    real(dp), parameter :: tol_rel = 1.0e-6_dp

    integer(i4b) :: i, j, iL
    integer(i4b), allocatable :: iCand(:), ieee_ok(:)
    real(dp), allocatable :: logTgrid(:), Fgrid(:)
    real(dp) :: logT_lo, logT_hi, dlogT, logT_mid
    real(dp) :: T_trial, cs2, H_trial, rho_trial, logR_trial
    real(dp) :: H_tmp, rho_tmp, nu_tmp, kappa_tmp, kappa_planck_tmp, tau_tmp
    real(dp) :: Qvis_tmp, Qirr_tmp, Qrad_tmp
    real(dp) :: f_lo, f_hi, f_mid

    n_roots = 0
    T_roots = 0.0_dp

    if (Sigma_cgs <= 1.0e-20_dp) return

    logT_lo = log10(max(T_floor, tlow_tab))
    logT_hi = log10(min(T_ceiling, thigh_tab))
    dlogT = (logT_hi - logT_lo) / real(Nscan, dp)

    allocate(logTgrid(Nscan+1), Fgrid(Nscan+1), ieee_ok(Nscan+1))
    ieee_ok(:) = 0

    do i = 1, Nscan + 1
       logTgrid(i) = logT_lo + dlogT * real(i-1, dp)
       T_trial = 10.0_dp**logTgrid(i)
       if (T_trial > thigh_tab .or. T_trial < tlow_tab) then
          ieee_ok(i) = 1
          cycle
       end if
       cs2 = kb * T_trial / (mu * mp)
       H_trial = sqrt(cs2) / OmegaK_cgs
       rho_trial = Sigma_cgs / (2.0_dp * H_trial)
       logR_trial = log10(rho_trial) - 3.0_dp*log10(T_trial) + 18.0_dp
       if (logR_trial < logRmin_tab .or. rho_trial < rhomin_tab) then
          ieee_ok(i) = 1
          cycle
       end if

       call heating_cooling_cell(r_cgs, Sigma_cgs, OmegaK_cgs, shadow, T_trial, &
            H_tmp, rho_tmp, nu_tmp, kappa_tmp, kappa_planck_tmp, tau_tmp, &
            Qvis_tmp, Qirr_tmp, Qrad_tmp, Qirr_in=Qirr_in)
       Fgrid(i) = (Qvis_tmp + Qirr_tmp) - Qrad_tmp
       if (ieee_is_nan(Fgrid(i)) .or. .not. ieee_is_finite(Fgrid(i))) ieee_ok(i) = 1
    end do

    allocate(iCand(Nscan))
    do i = 1, Nscan
       if (ieee_ok(i) /= 0 .or. ieee_ok(i+1) /= 0) cycle
       if (Fgrid(i) * Fgrid(i+1) < 0.0_dp) then
          if (n_roots < max_roots) then
             n_roots = n_roots + 1
             iCand(n_roots) = i
          end if
       end if
    end do

    do j = 1, n_roots
       iL = iCand(j)
       logT_lo = logTgrid(iL)
       logT_hi = logTgrid(iL+1)

       T_trial = 10.0_dp**logT_lo
       call heating_cooling_cell(r_cgs, Sigma_cgs, OmegaK_cgs, shadow, T_trial, &
            H_tmp, rho_tmp, nu_tmp, kappa_tmp, kappa_planck_tmp, tau_tmp, &
            Qvis_tmp, Qirr_tmp, Qrad_tmp, Qirr_in=Qirr_in)
       f_lo = (Qvis_tmp + Qirr_tmp) - Qrad_tmp

       T_trial = 10.0_dp**logT_hi
       call heating_cooling_cell(r_cgs, Sigma_cgs, OmegaK_cgs, shadow, T_trial, &
            H_tmp, rho_tmp, nu_tmp, kappa_tmp, kappa_planck_tmp, tau_tmp, &
            Qvis_tmp, Qirr_tmp, Qrad_tmp, Qirr_in=Qirr_in)
       f_hi = (Qvis_tmp + Qirr_tmp) - Qrad_tmp

       do i = 1, iter_max
          logT_mid = 0.5_dp * (logT_lo + logT_hi)
          T_trial = 10.0_dp**logT_mid
          call heating_cooling_cell(r_cgs, Sigma_cgs, OmegaK_cgs, shadow, T_trial, &
               H_tmp, rho_tmp, nu_tmp, kappa_tmp, kappa_planck_tmp, tau_tmp, &
               Qvis_tmp, Qirr_tmp, Qrad_tmp, Qirr_in=Qirr_in)
          f_mid = (Qvis_tmp + Qirr_tmp) - Qrad_tmp
          if (f_lo * f_mid <= 0.0_dp) then
             logT_hi = logT_mid
             f_hi = f_mid
          else
             logT_lo = logT_mid
             f_lo = f_mid
          end if
          if (abs(logT_hi - logT_lo) / max(abs(logT_hi), 1.0e-99_dp) < tol_rel) then
             T_roots(j) = 10.0_dp**(0.5_dp * (logT_lo + logT_hi))
             exit
          end if
       end do
       T_roots(j) = 10.0_dp**(0.5_dp * (logT_lo + logT_hi))
    end do

    deallocate(logTgrid, Fgrid, ieee_ok, iCand)

  end subroutine find_all_thermal_roots

  !------------------------------------------------------------
  ! output_thermal_stability_all_roots
  !
  ! For each radial cell, find ALL equilibrium solutions (Q+=Q-)
  ! and output each with its stability. Format: one line per root.
  ! Tmid_arr: simulation's adopted T at each radius (to mark adopted root).
  ! Columns: r_cgs  Sigma  i_root  n_roots  T  kappa  ...  stability  adopted
  !   adopted: 1 = this root was selected by simulation, 0 = not
  !------------------------------------------------------------
  subroutine output_thermal_stability_all_roots(nr_out, r_nd, sigma_nd, Tmid_arr, &
       Qirr_arr, is_shadow_arr, it, t_nd, fname_base)
    integer(i4b), intent(in) :: nr_out
    real(dp),     intent(in) :: r_nd(nr_out), sigma_nd(nr_out), Tmid_arr(nr_out)
    real(dp),     intent(in) :: Qirr_arr(nr_out)
    logical,      intent(in) :: is_shadow_arr(nr_out)
    integer(i4b), intent(in) :: it
    real(dp),     intent(in) :: t_nd
    character(len=*), intent(in), optional :: fname_base

    integer(i4b) :: iu, i, j, n_roots, j_adopted
    real(dp) :: r_cgs, Sigma_cgs, OmegaK_cgs
    real(dp) :: T_roots(max_roots)
    real(dp) :: dQplus_dT, dQminus_dT
    logical  :: thermally_stable
    real(dp) :: H_tmp, rho_tmp, nu_tmp, kappa_tmp, kappa_planck_tmp, tau_tmp
    real(dp) :: Qvis_tmp, Qirr_tmp, Qrad_tmp
    character(len=64) :: fname
    real(dp), parameter :: tiny_sigma = 1.0e-8_dp
    real(dp), parameter :: tiny_T = 1.0_dp

    if (present(fname_base)) then
       write(fname, '(a,"_thermal_stability_all_roots_",i8.8,".dat")') trim(fname_base), it
    else
       write(fname, '("thermal_stability_all_roots_",i8.8,".dat")') it
    end if

    open(newunit=iu, file=trim(fname), status='replace', action='write')
    write(iu, '(a, i8, a, 1pe14.6, a, 1pe14.6, a)') '# it = ', it, ', t_nd = ', t_nd, ' (t = ', t_dim(t_nd), ' [s])'
    write(iu, '(a, 1pe14.6, a, 1pe14.6, a)') '# t0 = ', t0, ', R0 = ', R0, ' [cm]'
    write(iu, '(a)') '# All equilibrium solutions at each (r, Sigma)'
    write(iu, '(a)') '# stability: 1 = stable, -1 = unstable'
    write(iu, '(a)') '# adopted: 1 = solution selected by simulation, 0 = not'
    write(iu, '(a)') '# r_cgs  Sigma  i_root  n_roots  T  kappa  Qvis  Qirr  Qrad  dQplus_dT  dQminus_dT  stability  adopted'
    write(iu, '(a)') '#'

    do i = 1, nr_out
       r_cgs = r_dim(r_nd(i))
       Sigma_cgs = sigma_dim(sigma_nd(i))
       OmegaK_cgs = omegaK_dim(r_nd(i))

       if (Sigma_cgs <= tiny_sigma) cycle

       call find_all_thermal_roots(r_cgs, Sigma_cgs, OmegaK_cgs, &
            is_shadow_arr(i), Qirr_arr(i), T_roots, n_roots)

       ! Which root is adopted? (closest to Tmid_arr)
       j_adopted = 0
       if (Tmid_arr(i) > tiny_T .and. n_roots > 0) then
          j_adopted = 1
          do j = 2, n_roots
             if (abs(T_roots(j) - Tmid_arr(i)) < abs(T_roots(j_adopted) - Tmid_arr(i))) &
                  j_adopted = j
          end do
       end if

       do j = 1, n_roots
          call heating_cooling_cell(r_cgs, Sigma_cgs, OmegaK_cgs, is_shadow_arr(i), &
               T_roots(j), H_tmp, rho_tmp, nu_tmp, kappa_tmp, kappa_planck_tmp, tau_tmp, &
               Qvis_tmp, Qirr_tmp, Qrad_tmp, Qirr_in=Qirr_arr(i))
          call compute_thermal_derivatives(r_cgs, Sigma_cgs, OmegaK_cgs, &
               is_shadow_arr(i), T_roots(j), Qirr_arr(i), &
               dQplus_dT, dQminus_dT, thermally_stable)
          if (thermally_stable) then
             write(iu, '(1p,7e13.5, 2i4, 1p,2e13.5, 2i4)') r_cgs, Sigma_cgs, &
                  T_roots(j), kappa_tmp, Qvis_tmp, Qirr_tmp, Qrad_tmp, &
                  j, n_roots, dQplus_dT, dQminus_dT, 1, merge(1, 0, j == j_adopted)
          else
             write(iu, '(1p,7e13.5, 2i4, 1p,2e13.5, 2i4)') r_cgs, Sigma_cgs, &
                  T_roots(j), kappa_tmp, Qvis_tmp, Qirr_tmp, Qrad_tmp, &
                  j, n_roots, dQplus_dT, dQminus_dT, -1, merge(1, 0, j == j_adopted)
          end if
       end do
    end do
    close(iu)

  end subroutine output_thermal_stability_all_roots

  !------------------------------------------------------------
  ! generate_scurve_at_radius
  !
  ! Generate S-curve: thermal equilibrium Sigma vs T at fixed r.
  ! For each Sigma in [Sigma_min, Sigma_max], find T such that
  ! Qvis + Qirr = Qrad. Optionally compute stability at each point.
  !
  ! Output file: Sigma, T, kappa, Qvis, Qrad, dQ+/dT, dQ-/dT, stability
  !------------------------------------------------------------
  subroutine generate_scurve_at_radius(r_cgs, OmegaK_cgs, shadow, Qirr_in, &
       Sigma_min, Sigma_max, nSigma, filename, T_min, T_max, Sigma_sim, Tmid_sim)
    real(dp), intent(in) :: r_cgs, OmegaK_cgs
    logical,  intent(in) :: shadow
    real(dp), intent(in) :: Qirr_in
    real(dp), intent(in) :: Sigma_min, Sigma_max
    integer(i4b), intent(in) :: nSigma
    character(len=*), intent(in) :: filename
    real(dp), intent(in), optional :: T_min, T_max, Sigma_sim, Tmid_sim

    integer(i4b) :: i, iu, nroots, ierr, j
    real(dp) :: Sigma_cgs, T_root, logT_lo, logT_hi, dlogT
    real(dp) :: logSigma_min, logSigma_max, logSigma
    real(dp) :: logT_min, logT_max
    real(dp) :: H_loc, rho_loc, nu_dim, kappa_loc, kappa_planck, tau_loc
    real(dp) :: Qvis_tmp, Qirr_tmp, Qrad_tmp
    real(dp) :: Tmid_roots(20)
    real(dp) :: dQplus_dT, dQminus_dT
    logical  :: thermally_stable
    integer(i4b), parameter :: Nscan = 400

    ! Temperature range for root search [K]; defaults ~30 K - 10^5 K
    if (present(T_min) .and. present(T_max)) then
       logT_min = log10(max(T_min, 1.0e-30_dp))
       logT_max = log10(max(T_max, 1.0e-30_dp))
    else
       logT_min = 1.5_dp
       logT_max = 5.0_dp
    end if

    open(newunit=iu, file=trim(filename), status='replace', action='write')
    write(iu, '(a, 1pe14.6, a, 1pe14.6, a)') '# S-curve at r = ', r_cgs, ' [cm] ( r/R0 = ', r_cgs/R0, ' )'
    if (present(Sigma_sim) .and. present(Tmid_sim)) then
       if (Sigma_sim > 0.0_dp .and. Tmid_sim > 0.0_dp) then
          write(iu, '(a, 1pe14.6, a, 1pe14.6, a)') &
               '# Sigma_cgs = ', Sigma_sim, ' [g/cm^2], Tmid = ', Tmid_sim, ' [K]'
       end if
    end if
    write(iu, '(a, 1pe14.6, a)') '# Qirr = ', Qirr_in, ' [erg/cm^2/s]'
    write(iu, '(a)') '# Sigma  T  branch  kappa  Qvis  Qrad  dQplus_dT  dQminus_dT  stability'

    ! Log scale for Sigma: uniform in log10(Sigma) for fine resolution over wide range
    logSigma_min = log10(max(Sigma_min, 1.0e-30_dp))
    logSigma_max = log10(max(Sigma_max, 1.0e-30_dp))

    do i = 0, nSigma
       logSigma = logSigma_min + (logSigma_max - logSigma_min) * real(i, dp) / max(1, nSigma)
       Sigma_cgs = 10.0_dp**logSigma

       if (Sigma_cgs <= 1.0e-20_dp) cycle

       ! Scan T to find ALL roots of F(T) = (Qvis+Qirr) - Qrad = 0
       nroots = 0
       dlogT  = (logT_max - logT_min) / real(Nscan, dp)

       do j = 1, Nscan
          logT_lo = logT_min + dlogT * real(j - 1, dp)
          logT_hi = logT_min + dlogT * real(j, dp)

          call find_root_bracket(r_cgs, Sigma_cgs, OmegaK_cgs, shadow, Qirr_in, &
               10.0_dp**logT_lo, 10.0_dp**logT_hi, T_root, ierr)

          if (ierr == 0) then
             call push_root_unique(T_root, Tmid_roots, nroots, 0.02_dp)  ! 2% in logT
          end if
       end do

       if (nroots == 0) cycle

       ! For each root, compute structure and stability and write out
       do j = 1, nroots
          T_root = Tmid_roots(j)

          call heating_cooling_cell(r_cgs, Sigma_cgs, OmegaK_cgs, shadow, T_root, &
               H_loc, rho_loc, nu_dim, kappa_loc, kappa_planck, tau_loc, &
               Qvis_tmp, Qirr_tmp, Qrad_tmp, Qirr_in=Qirr_in)

          call compute_thermal_derivatives(r_cgs, Sigma_cgs, OmegaK_cgs, &
               shadow, T_root, Qirr_in, dQplus_dT, dQminus_dT, thermally_stable)

          if (thermally_stable) then
             write(iu, '(1p,2e14.6,1x,i3,1x,5e14.6,i4)') Sigma_cgs, T_root, j, kappa_loc, &
                  Qvis_tmp, Qrad_tmp, dQplus_dT, dQminus_dT, 1
          else
             write(iu, '(1p,2e14.6,1x,i3,1x,5e14.6,i4)') Sigma_cgs, T_root, j, kappa_loc, &
                  Qvis_tmp, Qrad_tmp, dQplus_dT, dQminus_dT, -1
          end if
       end do

    end do
    close(iu)

  contains
    subroutine find_root_bracket(r_c, Sig_c, OmK_c, shad, Qirr, T_lo, T_hi, &
         T_out, ierr)
      real(dp), intent(in)  :: r_c, Sig_c, OmK_c, T_lo, T_hi
      logical,  intent(in)  :: shad
      real(dp), intent(in)  :: Qirr
      real(dp), intent(out) :: T_out
      integer(i4b), intent(out) :: ierr

      real(dp) :: f_lo, f_hi, f_mid, logT_lo, logT_hi, logT_mid
      real(dp) :: H_t, rho_t, nu_t, kap_t, kapP_t, tau_t, Qv_t, Qi_t, Qr_t
      integer(i4b) :: iter
      integer(i4b), parameter :: iter_max = 50
      real(dp), parameter :: tol_rel = 1.0e-6_dp

      call heating_cooling_cell(r_c, Sig_c, OmK_c, shad, T_lo, &
           H_t, rho_t, nu_t, kap_t, kapP_t, tau_t, Qv_t, Qi_t, Qr_t, Qirr_in=Qirr)
      f_lo = (Qv_t + Qi_t) - Qr_t

      call heating_cooling_cell(r_c, Sig_c, OmK_c, shad, T_hi, &
           H_t, rho_t, nu_t, kap_t, kapP_t, tau_t, Qv_t, Qi_t, Qr_t, Qirr_in=Qirr)
      f_hi = (Qv_t + Qi_t) - Qr_t

      if (f_lo * f_hi >= 0.0_dp) then
         ierr = 1
         return
      end if

      logT_lo = log10(T_lo)
      logT_hi = log10(T_hi)
      do iter = 1, iter_max
         logT_mid = 0.5_dp * (logT_lo + logT_hi)
         T_out = 10.0_dp**logT_mid
         call heating_cooling_cell(r_c, Sig_c, OmK_c, shad, T_out, &
              H_t, rho_t, nu_t, kap_t, kapP_t, tau_t, Qv_t, Qi_t, Qr_t, Qirr_in=Qirr)
         f_mid = (Qv_t + Qi_t) - Qr_t
         if (f_lo * f_mid <= 0.0_dp) then
            logT_hi = logT_mid
            f_hi = f_mid
         else
            logT_lo = logT_mid
            f_lo = f_mid
         end if
         if (abs(logT_hi - logT_lo) / max(abs(logT_hi), 1.0e-99_dp) < tol_rel) then
            ierr = 0
            T_out = 10.0_dp**(0.5_dp * (logT_lo + logT_hi))
            return
         end if
      end do
      ierr = 0
      T_out = 10.0_dp**(0.5_dp * (logT_lo + logT_hi))
    end subroutine find_root_bracket
  end subroutine generate_scurve_at_radius

  subroutine push_root_unique(Tnew, Troots, nroots, dlogT_merge)
   ! Store a root if it is not a near-duplicate of existing roots.
   ! Duplicate merging is done in log10(T) space for robustness.
   !
   ! Inputs:
   !   Tnew        : candidate root temperature [K]
   !   dlogT_merge : merge threshold in log10(T) (e.g., 0.02 ~ 5%)
   !
   ! In/Out:
   !   Troots(:)   : stored roots
   !   nroots      : number of stored roots
   real(dp), intent(in)    :: Tnew, dlogT_merge
   real(dp), intent(inout) :: Troots(:)
   integer(i4b), intent(inout) :: nroots

   integer(i4b) :: k
   real(dp) :: logTnew, logTk

   if (Tnew <= 0.0_dp) return

   logTnew = log10(Tnew)

   do k = 1, nroots
      logTk = log10(max(Troots(k), 1.0e-99_dp))
      if (abs(logTnew - logTk) <= dlogT_merge) then
         ! Keep the one closer to the center of the bracket is not tracked here;
         ! just ignore near-duplicates.
         return
      end if
   end do

   if (nroots + 1 > size(Troots)) then
      ! Too many roots; increase Troots size if needed.
      return
   end if

   nroots = nroots + 1
   Troots(nroots) = Tnew

   call sort_roots_ascending(Troots, nroots)
 end subroutine push_root_unique


 subroutine sort_roots_ascending(Troots, nroots)
   ! Simple insertion sort for small nroots.
   real(dp), intent(inout) :: Troots(:)
   integer(i4b), intent(in) :: nroots
   integer(i4b) :: i, j
   real(dp) :: key

   do i = 2, nroots
      key = Troots(i)
      j = i - 1
      do while (j >= 1 .and. Troots(j) > key)
         Troots(j+1) = Troots(j)
         j = j - 1
      end do
      Troots(j+1) = key
   end do
 end subroutine sort_roots_ascending

end module thermal_stability_analysis_mod

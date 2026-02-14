module evolve_try_mod
  !! Adaptive timestep wrapper for one diffusion+thermal step.
  !!
  !! Key design:
  !! - TRY phase must NOT write to global state at it+1 (Tmid/Qvis/Qrad/etc.).
  !! - We compute trial arrays locally and validate them.
  !! - Only after acceptance, we commit everything to mod_global(it+1,:) and nu(:).
  !!
  use kind_params, only : dp, i4b
  use mod_global,  only : nr, t0, Tmid, sigmat, nu, &
                          Tmid, H, rho, kappaR, tauR, Qvis, Qirr, dYdXi, Qrad, &
                          use_be_decretion, use_irradiation, is_shadow
  use units_disk_mod, only : r_dim, sigma_dim, omegaK_dim
  public :: evolve_try

contains

logical function evolve_try_to_target(it_out, t_nd0, dt_out, dr) result(ok)
  use kind_params, only : dp, i4b
  use mod_global,  only : nr, r
  use inflow_source_mod, only : compute_source
  use timestep_seed_mod, only : get_dt_seed, set_dt_seed
  implicit none
  integer(i4b), intent(in) :: it_out
  real(dp),     intent(in) :: t_nd0, dt_out
  real(dp),     intent(in) :: dr(nr)

  real(dp) :: t_local, remaining
  real(dp) :: dt_seed, dt_try, dt_used
  real(dp) :: source_n(nr), source_np1(nr)
  integer(i4b) :: attempt
  integer(i4b), parameter :: max_retry_dt = 12
  real(dp),    parameter :: dt_shrink = 0.5_dp
  real(dp),    parameter :: dt_min = 1.0e-12_dp  ! nd safeguard (tune)
  real(dp),    parameter :: dt_grow = 1.2_dp      ! optional mild growth

  ok = .false.

  t_local   = t_nd0
  remaining = dt_out

  !--------------------------------------------
  ! Seed dt_try from last success (bounded)
  !--------------------------------------------
  dt_seed = get_dt_seed()              ! module-saved value
  if (dt_seed <= 0.0_dp) dt_seed = dt_out
  dt_try  = min(dt_out, dt_seed)

  do while (remaining > 0.0_dp)

     ! never step beyond the output interval
     dt_try = min(dt_try, remaining)

     ! source at the beginning of this substep
     call compute_source(t_local, r, dr, source_n)

     !------------------------------------
     ! retry loop for THIS substep only
     !------------------------------------
     dt_used = dt_try
     do attempt = 1, max_retry_dt

        call compute_source(t_local + dt_used, r, dr, source_np1)

        ok = substep_try_and_commit(it_out, t_local, dt_used, source_n, source_np1)
        if (ok) exit

        dt_used = dt_used * dt_shrink
        if (dt_used < dt_min) then
           ok = .false.
           return
        end if
     end do

     if (.not. ok) return

     !------------------------------------
     ! substep succeeded and is already committed
     !------------------------------------
     t_local   = t_local + dt_used
     remaining = remaining - dt_used

     ! remember success for the NEXT call (store a *grown* seed)
     call set_dt_seed( min(dt_out, dt_used * dt_grow) )

     ! next try inside this same output interval
     dt_try = min( max(dt_used * dt_grow, dt_min), remaining )

  end do

  ok = .true.
end function evolve_try_to_target


  subroutine evolve_physics_one_substep_try( dt_loc, &
                               sigma_in, nu_in, T_in, &
                               source_n, source_np1, &
                               Qirr_prof_in, dYdXi_in, shadow_in, &
                               sigma_out, nu_out, &
                               T_out, H_out, rho_out, kappa_out, tau_out, &
                               Qvis_out, Qrad_out, Qirr_out, dYdXi_out, shadow_out, &
                               kiter_out, miter_out )
    use mod_global,  only : nr, r, nu, sigmat, alpha, &
                            use_energy_balance, use_energy_pde, &
                            use_be_decretion, use_wind_truncation, &
                            r_edge, i_edge, nu0_dim, nu0_nd, nu_conv
    use disk_energy_mod, only : solve_structure_from_sigma
    use irradiation_mod, only : L_irr, rin_cgs, A1, L1, beta1, beta2, Q12, &
         set_zero_irradiation, &
         set_irradiation_luminosity_from_arrays, &
         compute_shadow_loggrid_hyst, &
         build_Qirr_profile_eq16_with_raw, &
         fix_Qirr_negative_spikes
    implicit none

    real(dp),  intent(in) :: dt_loc
    real(dp),  intent(in) :: sigma_in(nr), nu_in(nr), T_in(nr)
    real(dp),  intent(in) :: source_n(nr), source_np1(nr)

    real(dp), intent(out) :: sigma_out(nr), nu_out(nr)
    real(dp), intent(out) :: T_out(nr), H_out(nr), rho_out(nr), kappa_out(nr), tau_out(nr)
    real(dp), intent(out) :: Qvis_out(nr), Qrad_out(nr), Qirr_out(nr), dYdXi_out(nr)
    logical,  intent(out) :: shadow_out(nr)
    integer(i4b), intent(out) :: kiter_out
    integer(i4b), intent(out) :: miter_out

    real(dp) :: sigma_old(nr)
    real(dp) :: sigma_star(nr), sigma_commit(nr)
    real(dp) :: src_theta(nr)
    real(dp) :: sigma_prev(nr), sigma_new(nr)
    real(dp) :: H_prev(nr), rho_prev(nr), kappa_prev(nr), tau_prev(nr),  &
                Qvis_prev(nr), Qrad_prev(nr), Qirr_prev(nr), dYdXi_prev(nr)
    real(dp) :: Qirr_prof_in(nr), dYdXi_in(nr)
    real(dp) :: T_prev(nr), T_trial(nr), T_relaxed(nr)
    real(dp) :: nu_new(nr), nu_trial(nr), nu_relaxed(nr)

    real(dp) :: H_trial(nr), rho_trial(nr)
    real(dp) :: kappa_trial(nr), tau_trial(nr)

    real(dp) :: Qvis_trial(nr), Qrad_trial(nr), Qirr_trial(nr)
    real(dp) :: dYdXi_trial(nr)
    logical  :: shadow_in(nr), shadow_prev(nr), shadow_trial(nr)

    real(dp) :: dt_phys
    real(dp) :: r_cgs(nr), Sigma_cgs(nr), OmegaK_cgs(nr)
    real(dp) :: Qirr_prof(nr)
    logical :: shadow_try(nr)
    logical :: converged

    real(dp) :: err_sig, err_T, err_nu, err_dQ
    real(dp) :: sig2, T2, nu2, dQ2, err2, err_rms
    real(dp), parameter :: tau_freeze = 1.0e30_dp
    real(dp), parameter :: tiny = 1.0e-12_dp
    real(dp), parameter :: omega_T  = 0.3_dp
    real (dp), parameter :: omega_nu = 0.3_dp

    real(dp), parameter :: eps_sigma = 1.0e-3_dp
    real(dp), parameter :: eps_T     = 1.0e-3_dp
    real(dp), parameter :: eps_nu    = 1.0e-3_dp
    real(dp), parameter :: eps_Q     = 1.0e-3_dp

    integer(i4b) :: i, k
    integer(i4b) :: iedge_local
    real(dp)     :: redge_cgs
    integer(i4b) :: i_maxresT, i_maxresnu
    real(dp) :: resT_max, resnu_max, resT, resnu

    integer(i4b), parameter :: iter_max = 100
    integer(i4b), save :: iu_res = -1, iu_conv = -1
    logical, save :: firstA = .true., firstB = .true.
    character(len=70) :: file91A, file91B

    real(dp) :: xi(nr), Ytmp(nr), dYdXi_tmp(nr)

    integer(i4b), parameter :: mmax_irr = 6
    real(dp),     parameter :: eps_rQ   = 1.0e-2_dp
    real(dp),     parameter :: w_init_irr = 0.3_dp
    real(dp),     parameter :: w_min_irr  = 0.02_dp
    real(dp),     parameter :: w_shrink_irr = 0.5_dp
    real(dp),     parameter :: grow_tol_irr = 1.05_dp
    integer(i4b), parameter :: blow_max_irr = 3

    ! shadow hysteresis parameters (same as substep_try_and_commit)
    real(dp),     parameter :: eps_on  = 3.0e-3_dp
    real(dp),     parameter :: eps_off = 1.5e-3_dp

    ! Optional: "soft shadow" (more physical: finite source + scattering)
    logical,      parameter :: use_soft_shadow = .false.
    real(dp),     parameter :: soft_f_mid = 1.0_dp - 0.5_dp*(eps_on + eps_off)
    real(dp),     parameter :: soft_df    = 0.5_dp*(eps_on - eps_off)   ! transition width

    real(dp) :: Qirr_calc(nr), Qirr_raw(nr), Y_prof(nr), dYdXi_calc(nr)
    logical  :: shadow_new(nr)
    real(dp) :: hoverr_raw(nr), hoverr_sm(nr), Hshadow_cgs(nr)
    real(dp) :: wloc_irr, rQ, numr, denr, relprev_rQ
    integer(i4b) :: m_irr, n_blow_irr

    file91A = 'residuals.txt'
    file91B = 'convergence_history.txt'

    sigma_old(:) = sigma_in(:)

    src_theta(:) = (1.0_dp - alpha) * source_n(:) + alpha * source_np1(:)

    T_prev(:) = T_in(:)

    if (.not. use_energy_balance) then
       ! Case A: isothermal (nu is fixed here)
       call diffusion_theta_step(nr, r, nu, sigma_old, dt_loc, alpha, sigma_new)

       sigma_star(:)   = sigma_new(:)
       sigma_commit(:) = max(0.0_dp, sigma_star(:) + dt_loc * src_theta(:))

       sigma_out(:) = sigma_commit(:)
       nu_out(:)    = nu(:)          ! unchanged in isothermal
       kiter_out    = 0

       ! For consistency, produce "dummy" structure outputs (zeros)
       T_out(:) = 0.0_dp; H_out(:)=0.0_dp; rho_out(:)=0.0_dp
       kappa_out(:)=0.0_dp; tau_out(:)=0.0_dp
       Qvis_out(:)=0.0_dp; Qrad_out(:)=0.0_dp
       Qirr_out(:)=0.0_dp; dYdXi_out(:)=0.0_dp
       shadow_out(:)=.true.

    else
       if (use_energy_pde) then
          !=========================================================
          ! Case: use_energy_pde  (TRY mode)
          !   - Advance Sigma first (theta diffusion + source)
          !   - Then advance T by local energy ODE (implicit Euler)
          !   - Recompute structure/nu/Q terms consistently from new T
          !   - NO fixed-point iteration here
          !=========================================================
          call diffusion_theta_step(nr, r, nu, sigma_old, dt_loc, alpha, sigma_new)

          sigma_star(:)   = sigma_new(:)
          sigma_commit(:) = max(0.0_dp, sigma_star(:) + dt_loc * src_theta(:))

          sigma_out(:) = sigma_commit(:)

          !---------------------------------------------
          ! Build CGS arrays for this TRY step
          !---------------------------------------------
          dt_phys = dt_loc * t0

          do i = 1, nr
             r_cgs(i)     = r_dim(r(i))
             Sigma_cgs(i) = sigma_dim(sigma_commit(i))
             if (r_cgs(i) > 0.0_dp) then
                OmegaK_cgs(i) = omegaK_dim(r(i))
             else
                OmegaK_cgs(i) = 0.0_dp
             end if
          end do

          !---------------------------------------------
          ! Irradiation / shadow + self-consistent fixed-point iteration
          !   Goal: make (T,H) consistent with Qirr(H) and shadow(H) within THIS substep
          !   Sigma_commit is held fixed (diffusion already done).
          !---------------------------------------------
          ! Initialize irradiation inputs for this TRY step
          if (.not. use_irradiation) then
             call set_zero_irradiation()
             Qirr_prof(:)  = 0.0_dp
             dYdXi_out(:)  = 0.0_dp
             shadow_try(:) = .false.
          else
             Qirr_prof(:)  = Qirr_prof_in(:)
             dYdXi_out(:)  = dYdXi_in(:)
             shadow_try(:) = shadow_in(:)
          end if

          ! Fixed-point loop: (T,H) -> Qirr(H) + shadow(H) -> update Qirr -> repeat
          wloc_irr     = w_init_irr
          relprev_rQ   = huge(1.0_dp)
          n_blow_irr   = 0

          do m_irr = 1, mmax_irr

             ! Solve energy PDE with the CURRENT (shadow_try, Qirr_prof)
             call energy_pde_step_try( nr, r_cgs, Sigma_cgs, OmegaK_cgs, shadow_try, Qirr_prof, &
                  dt_phys, T_prev, &
                  T_out, H_out, rho_out, kappa_out, tau_out, &
                  Qvis_out, Qrad_out, Qirr_out, &
                  nu_out, converged )

             if (.not. converged) then
                kiter_out    = -99
                miter_out    = -99
                nu_out(:)    = 0.0_dp
                T_out(:)     = 0.0_dp
                H_out(:)     = 0.0_dp
                rho_out(:)   = 0.0_dp
                kappa_out(:) = 0.0_dp
                tau_out(:)   = 0.0_dp
                Qvis_out(:)  = 0.0_dp
                Qrad_out(:)  = 0.0_dp
                Qirr_out(:)  = 0.0_dp
                shadow_out(:)= .true.
                return
             end if

             if (.not. use_irradiation) exit  ! nothing else to iterate

             !------------------------------------------------------
             ! Recompute shadow(H_out) using the SAME hysteresis logic
             !------------------------------------------------------
             do i = 1, nr
                hoverr_raw(i) = H_out(i) / max(r_cgs(i), 1.0e-99_dp)
             end do
             call smooth_profile_logr(n=nr, r_cgs=r_cgs, y_in=hoverr_raw, halfwin=5, y_out=hoverr_sm)
             do i = 1, nr
                Hshadow_cgs(i) = hoverr_sm(i) * r_cgs(i)
             end do
             call compute_shadow_loggrid_hyst(nr, r_cgs, Hshadow_cgs, halfwin=5, &
                  eps_on=eps_on, eps_off=eps_off, shadow=shadow_new)

             !------------------------------------------------------
             ! Recompute Qirr(H_out) using LOH24 Eq.(16) wrapper
             !------------------------------------------------------
             call build_Qirr_profile_eq16_with_raw(nr, r_cgs, H_out, Qirr_raw, Qirr_calc, &
                  dYdXi_calc, Y_prof)

             call fix_Qirr_negative_spikes(n=nr, r_cgs=r_cgs, shadow=shadow_new, &
                  Qirr_raw=Qirr_raw, Qirr_prof=Qirr_calc, &
                  max_run=3, Qirr_floor_abs=0.0_dp, use_log_r=.true.)

             if (use_soft_shadow) then
                ! Soft transition around the shadow boundary (more physical, smoother Jacobian)
                do i = 1, nr
                   ! f = (H/r)_sm / envelope is hidden inside compute_shadow_loggrid_hyst,
                   ! so here we approximate smoothing by using the local hoverr_sm ratio to max so far.
                   ! For a stricter implementation, you can re-evaluate the envelope and f(i) explicitly.
                   ! Simple robust option: keep hard shadow as default (use_soft_shadow=.false.).
                   if (shadow_new(i)) then
                      Qirr_calc(i) = 0.0_dp
                   end if
                end do
             else
                ! Hard mask (your current model)
                where (shadow_new)
                   Qirr_calc = 0.0_dp
                end where
             end if

             !------------------------------------------------------
             ! True fixed-point residual: ||Qirr_calc - Qirr_prof|| / ||Qirr_prof||
             !------------------------------------------------------
             numr = 0.0_dp; denr = 0.0_dp
!$omp parallel do default(shared) private(i) reduction(+:numr,denr)
             do i = 1, nr
                if (Sigma_cgs(i) <= 1.0e-8_dp) cycle
                numr = numr + (Qirr_calc(i) - Qirr_prof(i))**2
                denr = denr + max(abs(Qirr_prof(i)), 1.0e-99_dp)**2
             end do
!$omp end parallel do
             rQ = sqrt(numr / max(denr, 1.0e-99_dp))

             write (*, '("m_irr =", i3, ": rQ =", 1pe12.4)') m_irr, rQ

             if (rQ < eps_rQ) then
                ! Accept: update outputs so caller gets consistent arrays
                shadow_try(:) = shadow_new(:)
                dYdXi_out(:)  = dYdXi_calc(:)
                miter_out     = m_irr
                exit
             end if

             ! If residual is growing, shrink relaxation
             if (rQ > grow_tol_irr*relprev_rQ) then
                n_blow_irr = n_blow_irr + 1
                if (wloc_irr > w_min_irr) wloc_irr = max(w_min_irr, wloc_irr*w_shrink_irr)
             else
                n_blow_irr = 0
             end if
             if (n_blow_irr >= blow_max_irr) then
                ! Give up early; keep the latest converged PDE solution but do not claim self-consistency
                exit
             end if
             relprev_rQ = rQ

             ! Under-relax Qirr profile (fixed-point mixing)
             Qirr_prof(:)  = (1.0_dp - wloc_irr)*Qirr_prof(:) + wloc_irr*Qirr_calc(:)
             shadow_try(:) = shadow_new(:)
             dYdXi_out(:)  = dYdXi_calc(:)

             ! Better initial guess for the next PDE solve
             T_prev(:) = T_out(:)

          end do  ! m_irr

          shadow_out(:) = shadow_try(:)
          kiter_out     = 0
          miter_out     = -99

       else
          !=======================================================
          ! Legacy fixed-point iteration with UNDER-RELAXATION
          !   Unknowns: Sigma, Tmid, nu
          !=======================================================

          call diffusion_theta_step(nr, r, nu, sigma_old, dt_loc, alpha, sigma_new)

          sigma_prev(:) = sigma_new(:)
          nu_new(:)     = nu(:)

          ! ---------------------------------------------------------
          ! Seed "previous" state for relaxation WITHOUT time-history arrays.
          ! Use the substep inputs as the previous state.
          ! ---------------------------------------------------------
          T_prev(:) = T_in(:)
          nu_new(:) = nu_in(:)

          ! Build a consistent previous structure at (sigma_new, T_prev).
          ! This avoids starting H/rho/kappa/tau/Q arrays from zeros.
          call solve_structure_from_sigma( nr, r, sigma_new, T_prev, &
                                           nu_new, T_prev, H_prev, rho_prev, &
                                           kappa_prev, tau_prev, &
                                           Qvis_prev, Qrad_prev, Qirr_prev, &
                                           dYdXi_prev, shadow_prev, miter_out )

          kiter_out = 0

          do k = 1, iter_max
 
             sigma_prev(:) = sigma_new(:)

             !----------------------------------------------------
             ! (1) Solve thermal / structure problem at fixed Sigma
             !----------------------------------------------------
             !call solve_structure_from_sigma_freeze( nr, r, sigma_prev, T_prev,  &
             !     nu_new, H_prev, rho_prev, kappa_prev, tau_prev,                &
             !     Qvis_prev, Qrad_prev, Qirr_prev, dYdXi_prev, shadow_prev,      &
             !     tau_freeze,                                                    &
             !     nu_trial, T_trial, H_trial, rho_trial, kappa_trial, tau_trial, &
             !     Qvis_trial, Qrad_trial, Qirr_trial, dYdXi_trial, shadow_trial )
             call solve_structure_from_sigma( nr, r, sigma_prev, T_prev, &
                                    nu_trial, T_trial, H_trial, rho_trial, &
                                    kappa_trial, tau_trial,                &
                                    Qvis_trial, Qrad_trial, Qirr_trial,    &
                                    dYdXi_trial, shadow_trial, miter_out )

             !----------------------------------------------------
             ! (2) Under-relaxation in log(T) and linear nu
             !----------------------------------------------------
             do i = 1, nr
                if (T_prev(i) > 0.0_dp .and. T_trial(i) > 0.0_dp) then
                   T_relaxed(i) = exp((1.0_dp-omega_T)*log(T_prev(i)) &
                            + omega_T*log(T_trial(i)))
                else
                   T_relaxed(i) = T_trial(i)
                end if
                nu_relaxed(i) = (1.0_dp-omega_nu)*nu_new(i) + omega_nu*nu_trial(i)
             end do

             !----------------------------------------------------
             ! (3) Diffuse Sigma again with relaxed nu
             !----------------------------------------------------
             call diffusion_theta_step(nr, r, nu_relaxed, sigma_old, dt_loc, alpha, sigma_new)
             sigma_new(:) = max(0.0_dp, sigma_new(:) + dt_loc * src_theta(:))
             !----------------------------------------------------
             ! (4) Convergence diagnostics
             !----------------------------------------------------
             err_sig = 0.0_dp; sig2 = 0.0_dp
             err_T   = 0.0_dp; T2   = 0.0_dp
             err_nu  = 0.0_dp; nu2  = 0.0_dp
             err_dQ = 0.0_dp; dQ2  = 0.0_dp
             resT_max = 0.0_dp; resnu_max = 0.0_dp
             i_maxresT = 1; i_maxresnu = 1
             do i = 1, nr
                if (tau_trial(i) > 1.0_dp) then
                   err_sig = err_sig + (sigma_new(i) - sigma_prev(i))**2
                   sig2    = sig2    + max(abs(sigma_prev(i)), tiny)**2

                   resT = T_relaxed(i) - T_prev(i)
                   if (abs(resT) >= resT_max) then
                      i_maxresT = i
                      resT_max = resT
                   end if
                   err_T   = err_T   + resT*resT
                   T2      = T2      + max(abs(T_prev(i)), tiny)**2

                   resnu = nu_relaxed(i) - nu_new(i)
                   if (abs(resnu) >= resnu_max) then
                      i_maxresnu = i
                      resnu_max = resnu
                   end if
                   err_nu  = err_nu  + resnu*resnu
                   nu2     = nu2     + max(abs(nu_new(i)), tiny)**2

                   err_dQ  = err_dQ + (Qvis_trial(i) + Qirr_trial(i) - Qrad_trial(i))**2
                   dQ2     = dQ2    + max(Qvis_trial(i), Qirr_trial(i), Qrad_trial(i), tiny)**2
                end if
             end do

             err_sig = sqrt(err_sig / max(sig2, tiny))
             err_T   = sqrt(err_T   / max(T2,   tiny))
             err_nu  = sqrt(err_nu  / max(nu2,  tiny))
             err_dQ  = sqrt(err_dQ  / max(dQ2,  tiny))

             !-------------------------
             if (firstA) then
                open(newunit=iu_res, file=trim(file91A), status='replace', &
                     action='write')
                write (iu_res, '("  k  T/nu i_max  Max(residual)     tauR        Sigma        Tmid          nu")')
                firstA = .false.
             end if
             write (iu_res, '(i3, 1x, a3, 2x, i5, 2x, 1p5e13.5)') k, ' T', i_maxresT, &
                   resT_max, tau_trial(i_maxresT), sigma_new(i_maxresT), &
                   T_relaxed(i_maxresT), nu_relaxed(i_maxresT)
             write (iu_res, '(i3, 1x, a3, 2x, i5, 2x, 1p5e13.5)') k, 'nu', i_maxresnu, &
                   resnu_max, tau_trial(i_maxresnu), sigma_new(i_maxresnu), &
                   T_relaxed(i_maxresnu), nu_relaxed(i_maxresnu)
             !-------------------------

             !----------------------------------------------------
             ! (5) Convergence test (ALL variables)
             !----------------------------------------------------
             !if (err_sig < eps_sigma .and. err_T < eps_T .and. &
             !    err_nu < eps_nu .and. err_dQ < eps_Q) then
             if (err_sig < eps_sigma .and. &
                 err_nu < eps_nu .and. err_dQ < eps_Q) then
                kiter_out = k
                if (firstB) then
                   open(newunit=iu_conv, file=trim(file91B), status='replace', &
                        action='write')
                   firstB = .false.
                end if
                write (iu_conv, '("Converged at k =", i3)') k
                write (iu_conv, '("err_sig =", 1pe12.4, " err_T =", 1pe12.4, &
                       ", err_nu =", 1pe12.4, ", err_dQ =", 1pe12.4)') &
                       err_sig, err_T, err_nu, err_dQ
               exit
             end if

             !----------------------------------------------------
             ! (6) Prepare next iteration
             !----------------------------------------------------
             T_prev(:) = T_relaxed(:)
             nu_new(:) = nu_relaxed(:)

          end do

          if (kiter_out == 0) then
             kiter_out = -99
             if (firstB) then
                open(newunit=iu_conv, file=trim(file91B), status='replace', &
                     action='write')
                firstB = .false.
             end if
             write (iu_conv, '("Not converged afer", i4, " iterations")') &
                   iter_max
             write (iu_conv, '("err_sig =", 1pe12.4, " err_T =", 1pe12.4, &
                    ", err_nu =", 1pe12.4, ", err_dQ =", 1pe12.4)') &
                    err_sig, err_T, err_nu, err_dQ
             write (iu_conv, '("       r       sigma_new   sigma_old   T_relaxed    T_prev", &
                 & "    nu_relaxed    nu_new", &
                 & "    Qvis_trial  Qirr_trial  Qrad_trial", &
                 & "   Q^+ - Q^-")')
             do i = 1, nr
                write (iu_conv, '(1p11e12.4)') &
                     r(i), sigma_new(i), sigma_prev(i), T_relaxed(i), T_prev(i), &
                     nu_relaxed(i), nu_new(i), Qvis_trial(i), Qirr_trial(i), &
                     Qrad_trial(i), Qvis_trial(i) + Qirr_trial(i) - Qrad_trial(i)
             end do
          end if

          !------------------------------------------------------------
          ! After the iteration loop (legacy branch), ALWAYS fill outputs
          !------------------------------------------------------------

          ! 1) Add the source ONCE (time-centered) to obtain committed Sigma
          sigma_star(:)   = sigma_new(:)
          sigma_commit(:) = max(0.0_dp, sigma_star(:) + dt_loc * src_theta(:))

          sigma_out(:) = sigma_commit(:)

          if (kiter_out > 0) then
             ! 2) Final structure solve at committed Sigma to populate ALL output arrays
             !    Use the latest relaxed temperature as the seed (T_prev holds it)
             !call solve_structure_from_sigma_freeze( nr, r, sigma_commit, T_prev, &
             !     nu_new, H_prev, rho_prev, kappa_prev, tau_prev,                 &
             !     Qvis_prev, Qrad_prev, Qirr_prev, dYdXi_prev, shadow_prev,       &
             !     tau_freeze,                                                     &
             !     nu_out, T_out, H_out, rho_out, kappa_out, tau_out,              &
             !     Qvis_out, Qrad_out, Qirr_out, dYdXi_out, shadow_out )
             call solve_structure_from_sigma( nr, r, sigma_commit, T_prev, &
                                    nu_out, T_out, H_out, rho_out, &
                                    kappa_out, tau_out,            &
                                    Qvis_out, Qrad_out, Qirr_out,  &
                                    dYdXi_out, shadow_out, miter_out )
          else
             ! Not converged -> return "safe" zeros so accept_step_try can reject
             nu_out(:)     = 0.0_dp
             T_out(:)      = 0.0_dp
             H_out(:)      = 0.0_dp
             rho_out(:)    = 0.0_dp
             kappa_out(:)  = 0.0_dp
             tau_out(:)    = 0.0_dp
             Qvis_out(:)   = 0.0_dp
             Qrad_out(:)   = 0.0_dp
             Qirr_out(:)   = 0.0_dp
             dYdXi_out(:)  = 0.0_dp
             shadow_out(:) = .true.
             nu_out(:)     = 0.0_dp
          end if

       end if
    end if

    ! Optional wind truncation diagnosis SHOULD be done on accepted step only.
    ! If you absolutely need it for diagnostics during try, you can compute it here
    ! but do NOT write to global i_edge/r_edge until commit_step.

  end subroutine evolve_physics_one_substep_try


logical function substep_try_and_commit(it_out, t_local, dt_nd, source_n, source_np1) result(ok)
  use mod_global, only : nr, r, &
                         sigma_cur, nu_cur, Tmid_cur, H_cur, rho_cur, kappa_cur, tau_cur, &
                         Qvis_cur, Qrad_cur, Qirr_cur, dYdXi_cur, shadow_cur, &
                         use_energy_pde, use_irradiation, mdot_inner_nd, mdot_inner_phys
  use units_disk_mod, only : r_dim, mdot_edd_unit
  use irradiation_mod, only : set_zero_irradiation, &
                              set_irradiation_luminosity_from_arrays, &
                              compute_shadow_loggrid_hyst, &
                              build_Qirr_profile_eq16_with_raw, & ! a wrapper routine
                              fix_Qirr_negative_spikes
  use hot_region_metrics_mod, only : update_hot_region_metrics_after_commit
  implicit none
  integer(i4b), intent(in) :: it_out
  real(dp),     intent(in) :: t_local, dt_nd
  real(dp),     intent(in) :: source_n(nr), source_np1(nr)

  ! trial arrays
  real(dp) :: sigma_try(nr), nu_try(nr)
  real(dp) :: T_try(nr), H_try(nr), rho_try(nr), kappa_try(nr), tau_try(nr)
  real(dp) :: Qvis_try(nr), Qrad_try(nr), Qirr_try(nr), dYdXi_try(nr)
  logical  :: shadow_try(nr)
  integer(i4b) :: kiter_try, miter_try

  ! irradiation profile for THIS substep (computed at substep start geometry)
  real(dp) :: Qirr_prof(nr), dYdXi_prof(nr)
  real(dp) :: Qirr_raw(nr), Y_prof(nr)
  real(dp) :: hoverr(nr)
  integer(i4b), save :: iu_dbg = -1
  logical  :: shadow_prof(nr)
  logical, save :: dbg_first = .true.
  character(len=*), parameter :: dbgfile = 'diag_shadow_qirr.dat'

  ! Work array: H in CGS for LOH24
  real(dp) :: H_cgs(nr), r_cgs(nr)
  real(dp) :: mdot_tmp
  integer(i4b) :: i

  ! Arrays for smoothing H/r
  real(dp) :: Hshadow_cgs(nr)
  real(dp) :: hoverr_raw(nr), hoverr_sm(nr)

  !--------------------------------------------
  ! Irradiation update at the SUBSTEP START time
  !--------------------------------------------
  if (.not. use_irradiation) then
     call set_zero_irradiation()
     Qirr_prof(:)  = 0.0_dp
     dYdXi_prof(:) = 0.0_dp
     shadow_prof(:)= .false.
  else
     !--------------------------------------------
     ! Shadowing at SUBSTEP START geometry (no delay)
     !--------------------------------------------
     do i = 1, nr
        r_cgs(i) = r_dim(r(i)) 
     end do

     ! Build a smoothed H/r ONLY for shadowing to suppress single-cell spikes.
     do i = 1, nr
        hoverr_raw(i) = H_cur(i) / max(r_cgs(i), 1.0e-99_dp)   ! dimensionless
     end do

     call smooth_profile_logr(n=nr, r_cgs=r_cgs, y_in=hoverr_raw, halfwin=5, &
                             y_out=hoverr_sm)

     do i = 1, nr
        Hshadow_cgs(i) = hoverr_sm(i) * r_cgs(i)
     end do

     !!call compute_shadow_from_HoverR(nr, r_cgs, Hshadow_cgs, shadow_prof)
     !call compute_shadow_loggrid_smooth(nr, r_cgs, Hshadow_cgs, halfwin=5, &
     !                                  eps_shadow=3.0e-3_dp, shadow=shadow_prof)
     call compute_shadow_loggrid_hyst(nr, r_cgs, Hshadow_cgs, halfwin=5, &
                                     eps_on=3.0e-3_dp, eps_off=1.5e-3_dp, &
                                     shadow=shadow_prof)

     ! Update L_irr and LOH24 parameters using CURRENT state (cur arrays)
     call set_irradiation_luminosity_from_arrays(t_nd_now=t_local, sigma_nd=sigma_cur, nu_nd=nu_cur, do_push=.false.)

     ! Build Qirr(r) using CURRENT geometry (H_cur) at substep start
     do i = 1, nr
        H_cgs(i) = H_cur(i)   ! ensure CGS conversion if needed
     end do

     call build_Qirr_profile_eq16_with_raw(nr, r_cgs, H_cgs, Qirr_raw, Qirr_prof, &
                                          dYdXi_prof, Y_prof)

     ! Post-process: fill short negative spikes instead of hard zeroing
     call fix_Qirr_negative_spikes(n=nr, r_cgs=r_cgs, shadow=shadow_prof, &
                                  Qirr_raw=Qirr_raw, Qirr_prof=Qirr_prof, &
                                  max_run=3, Qirr_floor_abs=0.0_dp, use_log_r=.true.)

     ! Enforce shadow mask (keep this in ONE place: here)
     where (shadow_prof)
        Qirr_prof  = 0.0_dp
        ! dYdXi_prof = 0.0_dp   ! usually unnecessary; keep only if your solver assumes it
     end where

     !-----------------------------
     ! Diagnostics output (optional)
     !-----------------------------
     ! Example: dump every accepted output step or every N steps
     if (mod(it_out, 1000) == 0) then
        if (dbg_first) then
           open(newunit=iu_dbg, file=dbgfile, status='replace', action='write')
           write(iu_dbg,'(a)') '# it  t_nd  i  r_cgs  H_over_r  shadow  Y  dYdXi  Qirr_raw  Qirr_prof'
           dbg_first = .false.
        end if

        do i = 1, nr
           hoverr(i) = H_cur(i) / max(r_cgs(i), 1.0e-99_dp)
           write(iu_dbg,'(i8,1x,1pe14.6,1x,i5,1x,1pe14.6,1x,1pe14.6,1x,i1,1x,1pe14.6,1x,1pe14.6,1x,1pe14.6,1x,1pe14.6)') &
                it_out, t_local, i, r_cgs(i), hoverr(i), merge(1,0,shadow_prof(i)), Y_prof(i), dYdXi_prof(i), Qirr_raw(i), Qirr_prof(i)
        end do

        write(iu_dbg,'(a)') ''   ! blank line as a block separator
     end if
  end if

  !--------------------------------------------
  ! TRY step (uses Qirr_prof + shadow_prof)
  !--------------------------------------------
  call evolve_physics_one_substep_try( &
       dt_nd, &
       sigma_cur, nu_cur, Tmid_cur, &                   ! input = cur
       source_n, source_np1, &
       Qirr_prof, dYdXi_prof, shadow_prof, &            ! ★ NEW
       sigma_try, nu_try, T_try, H_try, rho_try, kappa_try, tau_try, &
       Qvis_try, Qrad_try, Qirr_try, dYdXi_try, shadow_try, &
       kiter_try, miter_try )

  ok = accept_step_try(nr, kiter_try, sigma_try, nu_try, T_try, Qvis_try, Qrad_try, Qirr_try)
  if (.not. ok) return

  call commit_step_cur( sigma_try, nu_try, &
       T_try, H_try, rho_try, kappa_try, tau_try, &
       Qvis_try, Qrad_try, Qirr_try, dYdXi_try, shadow_try, &
       kiter_try, miter_try )

  ! ---- after commit of the substep (ACCEPTED) ----
  ! After commit: update hot-region metrics once (accepted state)
  call update_hot_region_metrics_after_commit(t_nd=t_local+dt_nd, r_nd=r, Tmid=T_try)

  if (use_irradiation) then
     call set_irradiation_luminosity_from_arrays( t_nd_now=t_local + dt_nd, &
             sigma_nd=sigma_try, nu_nd=nu_try, do_push=.true. )
  end if
  ok = .true.
end function substep_try_and_commit


subroutine smooth_profile_logr(n, r_cgs, y_in, halfwin, y_out)
  ! Smooth y(r) on a log-r grid using a robust two-stage filter:
  !  (1) 3-point median prefilter (kills isolated spikes)
  !  (2) moving average over +/- halfwin (gentle smoothing)
  use kind_params, only: dp, i4b
  implicit none
  integer(i4b), intent(in) :: n, halfwin
  real(dp),     intent(in) :: r_cgs(n), y_in(n)
  real(dp),     intent(out):: y_out(n)

  integer(i4b) :: i, j, j1, j2, cnt
  real(dp) :: y_med(n), sumv

  ! 3-point median prefilter
  y_med(:) = y_in(:)
  if (n >= 3) then
     do i = 2, n-1
        y_med(i) = median3(y_in(i-1), y_in(i), y_in(i+1))
     end do
     y_med(1) = y_in(1)
     y_med(n) = y_in(n)
  end if

  ! Moving average (index window; for log grid this is effectively log-r smoothing)
  do i = 1, n
     j1 = max(1, i-halfwin)
     j2 = min(n, i+halfwin)
     sumv = 0.0_dp
     cnt  = 0
     do j = j1, j2
        sumv = sumv + y_med(j)
        cnt  = cnt + 1
     end do
     y_out(i) = sumv / max(1, cnt)
  end do

contains
  pure real(dp) function median3(a, b, c)
    real(dp), intent(in) :: a, b, c
    if ((a <= b .and. b <= c) .or. (c <= b .and. b <= a)) then
       median3 = b
    else if ((b <= a .and. a <= c) .or. (c <= a .and. a <= b)) then
       median3 = a
    else
       median3 = c
    end if
  end function median3
end subroutine smooth_profile_logr


  logical function accept_step_try(nr, kiter_try, sigma_try, nu_try, T_try, &
                                   Qvis_try, Qrad_try, Qirr_try) result(ok)
    use, intrinsic :: ieee_arithmetic
    implicit none
    integer(i4b), intent(in) :: nr
    integer(i4b), intent(in) :: kiter_try
    real(dp),     intent(in) :: sigma_try(nr), nu_try(nr), T_try(nr), &
                                Qvis_try(nr), Qrad_try(nr), Qirr_try(nr)

    integer(i4b) :: i
    real(dp) :: denom, qscale, qfloor, qactive, res
    real(dp), parameter :: tiny = 1.0e-60_dp
    real(dp), parameter :: eps_Q = 1.0e-2_dp
    real(dp), parameter :: Tmin  = 1.0_dp
    real(dp), parameter :: Tmax  = 1.0e7_dp

    ok = .true.

    if (kiter_try == -99) ok = .false.

    do i = 1, nr
       if (ieee_is_nan(T_try(i)) .or. .not. ieee_is_finite(T_try(i))) ok = .false.
       if (ieee_is_nan(nu_try(i)) .or. .not. ieee_is_finite(nu_try(i))) ok = .false.
    end do

    ! Reject a "zero structure" solution if Sigma is non-trivial
    if (maxval(sigma_try) > tiny) then
       if (maxval(T_try) <= 0.0_dp) ok = .false.
       if (maxval(nu_try) <= 0.0_dp) ok = .false.
    end if

    !-----------------------------------------------------------------
    ! 3) thermal balance residual (robust)
    !   Evaluate only where heating/cooling is "active".
    !-----------------------------------------------------------------
    ! A reference floor to avoid divide-by-tiny in near-vacuum cells.
    ! You can start conservative and tune:
    qfloor  = 1.0e-30_dp   ! [same units as Q*]; will be overridden by max below

    ! "Active" threshold: below this, skip the relative-residual test.
    ! Recommended: tie it to something that scales with your run.
    ! Here we use a fraction of global maxima in this try-step arrays:
    qactive = 1.0e-10_dp * max( maxval(abs(Qrad_try(:))), &
              maxval(abs(Qvis_try(:) + Qirr_try(:))) )

    do i = 1, nr

       qscale = max( abs(Qvis_try(i) + Qirr_try(i)), abs(Qrad_try(i)) )

       ! If the cell is essentially inactive (nearly zero heating/cooling),
       ! skip the relative residual check (otherwise it always blows up).
       if (qscale < qactive) cycle

       ! Use a robust denominator: representative thermal scale (not Qvis alone).
       denom = max(qscale, qfloor)

       res = abs( (Qvis_try(i) + Qirr_try(i)) - Qrad_try(i) ) / denom

       if (res > eps_Q) then
          !!ok = .false.
          !write(*,'("i=",i4,"  res=",1pe12.4,"  Qvis=",1pe12.4, &
          !        "  Qirr=",1pe12.4,"  Qrad=",1pe12.4)') i, res, &
          !        Qvis_try(i), Qirr_try(i), Qrad_try(i)
       end if
    end do

  end function accept_step_try


  subroutine commit_step(it, sigma_try, nu_try, &
                         T_try, H_try, rho_try, kappa_try, tau_try, &
                         Qvis_try, Qrad_try, Qirr_try, dYdXi_try, shadow_try, &
                         kiter_try, miter_try)
    use mod_global, only : nr, sigmat, nu, nu_conv, k_iter, m_iter
    use disk_energy_mod, only : commit_structure_to_global
    implicit none

    integer(i4b), intent(in) :: it
    real(dp),     intent(in) :: sigma_try(nr), nu_try(nr)
    real(dp),     intent(in) :: T_try(nr), H_try(nr), rho_try(nr), kappa_try(nr), tau_try(nr)
    real(dp),     intent(in) :: Qvis_try(nr), Qrad_try(nr), Qirr_try(nr), dYdXi_try(nr)
    logical,      intent(in) :: shadow_try(nr)
    integer(i4b), intent(in) :: kiter_try, miter_try

    integer(i4b) :: itp1
    itp1 = it + 1

    sigmat(itp1, :) = sigma_try(:)
    Tmid(itp1,:)    = T_try(:)
    nu(:)           = nu_try(:)
    nu_conv(itp1,:) = nu_try(:)
    k_iter(itp1)    = kiter_try
    m_iter(itp1)    = miter_try

    call commit_structure_to_global( itp1, &
         T_try, H_try, rho_try, kappa_try, tau_try, &
         Qvis_try, Qrad_try, Qirr_try, dYdXi_try, shadow_try, nu_try, miter_try )

  end subroutine commit_step

  subroutine commit_step_cur(sigma_try, nu_try, &
                           T_try, H_try, rho_try, kappa_try, tau_try, &
                           Qvis_try, Qrad_try, Qirr_try, dYdXi_try, shadow_try, &
                           kiter_try, miter_try)
    use mod_global, only : nr, sigma_cur, Tmid_cur, nu, nu_cur, k_iter_cur, &
                         m_iter_cur, H_cur, rho_cur, kappa_cur, tau_cur, &
                         Qvis_cur, Qrad_cur, Qirr_cur, dYdXi_cur, shadow_cur
    implicit none

    real(dp), intent(in) :: sigma_try(nr), nu_try(nr)
    real(dp), intent(in) :: T_try(nr), H_try(nr), rho_try(nr), kappa_try(nr), tau_try(nr)
    real(dp), intent(in) :: Qvis_try(nr), Qrad_try(nr), Qirr_try(nr), dYdXi_try(nr)
    logical,  intent(in) :: shadow_try(nr)
    integer(i4b), intent(in) :: kiter_try, miter_try

    sigma_cur(:) = sigma_try(:)
    Tmid_cur(:)  = T_try(:)

    nu(:)        = nu_try(:)
    nu_cur(:)    = nu_try(:)
    k_iter_cur   = kiter_try
    m_iter_cur   = miter_try

    H_cur(:)        = H_try(:)
    rho_cur(:)      = rho_try(:)
    kappa_cur(:)    = kappa_try(:)
    tau_cur(:)      = tau_try(:)
    Qvis_cur(:)     = Qvis_try(:)
    Qrad_cur(:)     = Qrad_try(:)
    Qirr_cur(:)     = Qirr_try(:)
    dYdXi_cur(:)    = dYdXi_try(:)
    shadow_cur(:)   = shadow_try(:)

   end subroutine commit_step_cur


  subroutine energy_pde_step_try( n, r_cgs, Sigma_cgs, OmegaK_cgs, shadow, Qirr_prof, &
                                dt_phys, T_old, &
                                T_new, H_new, rho_new, kappa_new, tau_new, &
                                Qvis_new, Qrad_new, Qirr_new, &
                                nu_nd_new, converged )
    use kind_params, only : dp, i4b
    use constants,   only : kb, mp, mu
    use mod_global,  only : nu0_dim, nu0_nd
    use radiation_params_mod, only : T_floor, T_ceiling
    use disk_thermal_mod, only : heating_cooling_cell
    use disk_energy_pde_mod, only : apply_radial_Tdiff_implicit
    use, intrinsic :: ieee_arithmetic
    implicit none

    integer(i4b), intent(in)  :: n
    real(dp),     intent(in)  :: r_cgs(n), Sigma_cgs(n), OmegaK_cgs(n)
    logical,      intent(in)  :: shadow(n)
    real(dp),     intent(in)  :: Qirr_prof(n)
    real(dp),     intent(in)  :: dt_phys
    real(dp),     intent(in)  :: T_old(n)

    real(dp),     intent(out) :: T_new(n), H_new(n), rho_new(n), kappa_new(n), tau_new(n)
    real(dp),     intent(out) :: Qvis_new(n), Qrad_new(n), Qirr_new(n)
    real(dp),     intent(out) :: nu_nd_new(n)
    logical,      intent(out) :: converged

    integer(i4b) :: i, iter
    real(dp) :: cv, sigcv, Told, T, Hloc, rholoc, nudim, kap, tau
    real(dp) :: Qv, Qi, Qm
    real(dp) :: R, Rp, dRdT, dT, eps
    real(dp) :: Rdt_over_T, dTrel, Qbal
    real(dp) :: T_tmp(n), H_tmp(n), rho_tmp(n), kappa_tmp(n)
    logical  :: active(n)
    logical  :: cell_ok

    integer(i4b), parameter :: iter_max = 30
    real(dp),    parameter :: tiny   = 1.0e-99_dp
    real(dp),    parameter :: tiny_sigma = 1.0e-8_dp
    real(dp),    parameter :: tiny_tau = 1.0e-6_dp
    real(dp), parameter :: tol_dT = 1.0e-3_dp
    real(dp), parameter :: tol_Rdt = 1.0e-3_dp   ! tolerance for |R|*dt/T

    ! For ideal monoatomic gas gamma=5/3:
    ! cV = (3/2) kB / (mu mp)
    cv = 1.5_dp * kb / (mu * mp)

    converged = .true.

    T_new(:)     = 0.0_dp
    H_new(:)     = 0.0_dp
    rho_new(:)   = 0.0_dp
    kappa_new(:) = 0.0_dp
    tau_new(:)   = 0.0_dp
    Qvis_new(:)  = 0.0_dp
    Qrad_new(:)  = 0.0_dp
    Qirr_new(:)  = 0.0_dp
    nu_nd_new(:) = 0.0_dp

    do i = 1, n

       if (Sigma_cgs(i) <= tiny_sigma .or. OmegaK_cgs(i) <= 0.0_dp .or. r_cgs(i) <= 0.0_dp) then
          T_tmp(i)     = 0.0_dp
          H_tmp(i)     = 0.0_dp
          rho_tmp(i)   = 0.0_dp
          kappa_tmp(i) = 0.0_dp
          active(i)    = .false.
          cycle
       end if

       sigcv = max(Sigma_cgs(i) * cv, tiny)

       Told = T_old(i)
       if (Told <= 0.0_dp) Told = T_floor
       Told = max(T_floor, min(Told, T_ceiling))

       !-----------------------------------------
       ! Newton for implicit Euler residual:
       ! R(T) = (T - Told)/dt - (Qv(T) + Qirr - Qm(T)) / (Sigma*cV)
       !-----------------------------------------
       T = Told

       cell_ok = .false.
       do iter = 1, iter_max

          call heating_cooling_cell( r_cgs(i), Sigma_cgs(i), OmegaK_cgs(i), shadow(i), T, &
                                 Hloc, rholoc, nudim, kap, tau, Qv, Qi, Qm, &
                                 Qirr_in=Qirr_prof(i) )

          if (tau < tiny_tau) then
             !-----------------------------------------
             ! Optically thin fallback:
             ! Do NOT treat this as Newton convergence.
             ! Freeze temperature to a safe value and
             ! accept this cell without Newton iteration.
             !-----------------------------------------
             T = max(T_floor, Told)
             cell_ok = .true.
             exit
          end if

          R  = (T - Told)/max(dt_phys,tiny) - ( (Qv + Qi - Qm) / sigcv )
          eps  = max(1.0e-2_dp * max(T, T_floor), 1.0_dp)

          call heating_cooling_cell( r_cgs(i), Sigma_cgs(i), OmegaK_cgs(i), shadow(i), T+eps, &
                                 Hloc, rholoc, nudim, kap, tau, Qv, Qi, Qm, &
                                 Qirr_in=Qirr_prof(i) )

          Rp = ( (T+eps) - Told)/max(dt_phys,tiny) - ( (Qv + Qi - Qm) / sigcv )

          dRdT = (Rp - R) / max(eps, tiny)
          if (abs(dRdT) < tiny) exit   ! derivative failure -> treat as not converged

          dT = -R / dRdT

          ! Simple line-search damping: keep the step bounded and inside temperature limits
          dT = max(min(dT, 0.5_dp*max(T, T_floor)), -0.5_dp*max(T, T_floor))
          T  = T + dT
          T  = max(T_floor, min(T, T_ceiling))

          if (.not. ieee_is_finite(T)) exit

          Rdt_over_T = abs(R) * dt_phys / max(abs(T), 1.0_dp)
          dTrel = abs(dT)  /max(abs(T), 1.0_dp)
          Qbal = Qv + Qi - Qm
          !if (i > n - 30) then
          !   write (*, '("i =", i4, ", iter =", i3, ":")') i, iter
          !   write (*, '(3x, "Rdt_over_T =", 1pe12.4, ", dTrel =", 1pe12.4, &
          !          ", Qbal =", 1pe12.4)') Rdt_over_T, dTrel, Qbal
          !end if
          if ( abs(dT)/max(abs(T),1.0_dp) < tol_dT .and. Rdt_over_T < tol_Rdt ) then
             cell_ok = .true.
             !if (i > n - 30) then
             !   write (*, '(3x, "-> Converged after ", i3, " iterations")') iter
             !end if
             exit
          end if

       end do

       if (.not. cell_ok) then
          converged = .false.
          return
       end if

       T_tmp(i) = T

       ! Call this subroutine to obtain coefficients of diffusion constant
       call heating_cooling_cell( r_cgs(i), Sigma_cgs(i), OmegaK_cgs(i), &
                                shadow(i), T_tmp(i), &
                                Hloc, rholoc, nudim, kap, tau, Qv, Qi, Qm, &
                                Qirr_in=Qirr_prof(i) )

       H_tmp(i)     = Hloc
       rho_tmp(i)   = rholoc
       kappa_tmp(i) = kap

       ! Cells to inhibit diffusion
       if (Sigma_cgs(i) <= tiny_sigma) then
          H_tmp(i)     = 0.0_dp
          kappa_tmp(i) = 0.0_dp
       end if
       if (tau < tiny_tau) then
          H_tmp(i)     = 0.0_dp
          kappa_tmp(i) = 0.0_dp
       end if
    end do

    call apply_radial_Tdiff_implicit( &
                    n         = n,       &
                    dt        = dt_phys, &
                    r_cgs     = r_cgs,   &
                    Sigma_cgs = Sigma_cgs, &
                    H_cgs     = H_tmp,   &
                    kappa     = kappa_tmp, &
                    cv        = cv,      &
                    T_inout   = T_tmp)

    do i = 1, n
       if (Sigma_cgs(i) <= tiny_sigma .or. OmegaK_cgs(i) <= 0.0_dp .or. r_cgs(i) <= 0.0_dp) then
          T_new(i)     = 0.0_dp
          H_new(i)     = 0.0_dp
          rho_new(i)   = 0.0_dp
          kappa_new(i) = 0.0_dp
          tau_new(i)   = 0.0_dp
          Qvis_new(i)  = 0.0_dp
          Qrad_new(i)  = 0.0_dp
          Qirr_new(i)  = 0.0_dp
          nu_nd_new(i) = 0.0_dp
          cycle
       end if

       T_new(i) = max(T_floor, min(T_tmp(i), T_ceiling))

       call heating_cooling_cell( r_cgs(i), Sigma_cgs(i), OmegaK_cgs(i), shadow(i), T_new(i), &
                             Hloc, rholoc, nudim, kap, tau, Qv, Qi, Qm, &
                             Qirr_in=Qirr_prof(i) )

       H_new(i)     = Hloc
       rho_new(i)   = rholoc
       kappa_new(i) = kap
       tau_new(i)   = tau
       Qvis_new(i)  = Qv
       Qrad_new(i)  = Qm
       Qirr_new(i)  = Qi

       if (nu0_dim > 0.0_dp) then
          nu_nd_new(i) = (nudim / nu0_dim) * nu0_nd
       else
          nu_nd_new(i) = 0.0_dp
       end if

    end do

  end subroutine energy_pde_step_try


  !-----------------------------------------------------------------
  subroutine compute_local_structure(r_cgs, Sigma_cgs, OmegaK_cgs, TT, &
                                     H_out, rho_out, nu_out, kappa_out, tau_out)
    !! Fallback structure computation at given TT (no root find).
    !!
    !! This is your compute_local_structure_be, but in generic naming.
    use constants, only : kb, mu, mp
    use mod_global, only : kappa0, kappa_es, alphaSS
    use opacity_table_mod, only : opacity_tables, findkappa_OP_AES_S03, &
                                findkappa_OP_S03, findkappa_OP_F05_S03
    use radiation_params_mod, only : T_floor, T_ceiling
    implicit none

    real(dp), intent(in)  :: r_cgs, Sigma_cgs, OmegaK_cgs, TT
    real(dp), intent(out) :: H_out, rho_out, nu_out, kappa_out, tau_out
    real(dp) :: cs2, rho_mid, kappa_ff, kappa_tab, T_loc
    integer(i4b) :: ierror_kap

    if (OmegaK_cgs <= 0.0_dp .or. Sigma_cgs <= 0.0_dp) then
        H_out     = 0.0_dp
       rho_out   = 0.0_dp
       kappa_out = 0.0_dp
       nu_out    = 0.0_dp
       tau_out   = 0.0_dp
       return
    end if

    T_loc   = max( T_floor, min(TT, T_ceiling) )
    cs2     = kb * T_loc / (mu * mp)
    H_out   = sqrt(cs2) / OmegaK_cgs
    rho_mid = Sigma_cgs / (2.0_dp * H_out)

    select case (opacity_tables)
    case ('OP+AES+S03')
       call findkappa_OP_AES_S03(rho_mid, T_loc, kappa_tab, ierror_kap)
    case ('OP+F05+S03')
       call findkappa_OP_F05_S03(rho_mid, T_loc, kappa_tab, ierror_kap)
    !case ('OP+AES')
    !   call findkappa_OP_AES(rho_mid, T_loc, kappa_tab, ierror_kap)
    case ('OP+S03')
       call findkappa_OP_S03(rho_mid, T_loc, kappa_tab, ierror_kap)
    end select

    if (ierror_kap /= 0) then
       kappa_ff  = kappa0 * rho_mid * T_loc**(-3.5_dp)
       kappa_out = kappa_es + kappa_ff
    else
       kappa_out = kappa_tab
    end if

    tau_out = 0.5_dp * kappa_out * Sigma_cgs
    nu_out  = alphaSS * cs2 / OmegaK_cgs
  end subroutine compute_local_structure


  !================================================================
  subroutine diagnose_wind_edge(itp1, sigma_in, iedge_local, redge_cgs)
  !================================================================
    use kind_params, only : dp, i4b
    use mod_global,  only : nr, r, H, rho, M_star, R0
    use constants,   only : gg
    use units_disk_mod, only : omegaK_dim, sigma_dim
    use wind_mod, only : initialized, init_wind_model, rho_wind, vw_wind
    implicit none

    integer(i4b), intent(in)  :: itp1
    real(dp),     intent(in)  :: sigma_in(nr)
    integer(i4b), intent(out) :: iedge_local
    real(dp),     intent(out) :: redge_cgs

    integer(i4b) :: i
    real(dp) :: r_cgs_i, H_i, rho_i, rho_w, vw_rel, vphiK
    real(dp) :: gz, f_stab, k, lambda_fac
    real(dp), parameter :: tiny_sig = 1.0e-20_dp

    if (.not. initialized) call init_wind_model(nr, r*R0)

    lambda_fac  = 1.0_dp
    iedge_local = nr + 1
    redge_cgs   = 0.0_dp

    do i = 1, nr
       if (sigma_dim(sigma_in(i)) <= tiny_sig) cycle

       H_i   = H(itp1, i)
       rho_i = rho(itp1, i)
       rho_w = rho_wind(i)
       if (H_i <= 0.0_dp .or. rho_i <= 0.0_dp .or. rho_w <= 0.0_dp) cycle

       r_cgs_i = r(i) * R0
       vphiK   = r_cgs_i * omegaK_dim(r(i))
       vw_rel  = sqrt(vw_wind(i)**2 + vphiK**2)
       if (vw_rel <= 0.0_dp) cycle

       gz = gg * M_star * H_i / ((r_cgs_i*r_cgs_i + H_i*H_i)**1.5_dp)
       k  = 1.0_dp / H_i

       f_stab = gz * (rho_i**2 - rho_w**2) / (rho_i*rho_w*vw_rel**2) - k

       if (f_stab < 0.0_dp) then
          iedge_local = i
          redge_cgs   = r_cgs_i
          exit
       end if
    end do

  end subroutine diagnose_wind_edge

  subroutine apply_outer_isothermal_cap(nr, itp1, r_nd, sigma_new, nu_new)
    use kind_params,  only : dp, i4b
    use constants,    only : kb, mp, mu
    use mod_global,   only : Teff_star, Tmid, R_star, H, rho, &
                             kappaR, tauR, nu0_dim, nu0_nd, nu_conv, &
                             use_irradiation, use_be_decretion
    use units_disk_mod, only : r_dim, sigma_dim, omegaK_dim, sigma_nd
    implicit none

    integer(i4b), intent(in)    :: nr, itp1
    real(dp),     intent(in)    :: r_nd(nr), sigma_new(nr)
    real(dp),     intent(inout) :: nu_new(nr)

    integer(i4b) :: i, i_cap_start
    real(dp) :: T_iso
    real(dp) :: r_cgs, Sigma_cgs, OmegaK_cgs
    real(dp) :: H_loc, rho_loc, kappa_loc, tau_loc, nu_dim_loc
    real(dp), parameter :: tiny_sigma = 1.0e-8_dp  ! keep consistent with caller
    real(dp), parameter :: cap_tau_thresh = 4.0_dp / 3.0_dp
    real(dp), parameter :: cap_Teff_factor = 0.6_dp
    logical, parameter  :: use_outer_iso_cap = .true.

    if (.not. use_be_decretion) return
    if (.not. use_outer_iso_cap) return

    T_iso = cap_Teff_factor * Teff_star
    i_cap_start = nr + 1

    ! Find the first radius where the disk becomes optically thin and "too hot"
    do i = 1, nr
       Sigma_cgs = sigma_dim(sigma_new(i))
       if (Sigma_cgs <= tiny_sigma) cycle

       if (tauR(itp1,i) <= cap_tau_thresh .and. Tmid(itp1,i) >= T_iso) then
          i_cap_start = i
          exit
       end if
    end do

    if (i_cap_start > nr) return

    ! Apply the cap outward and recompute structure consistently at T_iso
    do i = i_cap_start, nr
       Sigma_cgs = sigma_dim(sigma_new(i))
       if (Sigma_cgs <= tiny_sigma) cycle

       r_cgs     = r_dim(r_nd(i))
       OmegaK_cgs = omegaK_dim(r_nd(i))
       if (OmegaK_cgs <= 0.0_dp) cycle

       Tmid(itp1,i) = T_iso

       call compute_local_structure(r_cgs, Sigma_cgs, OmegaK_cgs, T_iso, &
                                  H_loc, rho_loc, nu_dim_loc, kappa_loc, tau_loc)

       H(itp1,i)      = H_loc
       rho(itp1,i)    = rho_loc
       kappaR(itp1,i) = kappa_loc
       tauR(itp1,i)   = tau_loc

       if (nu0_dim > 0.0_dp) then
          nu_new(i) = (nu_dim_loc / nu0_dim) * nu0_nd
       else
          !nu_new(i) = 0.0_dp
          nu_new(i) = nu_conv(itp1-1,i)
       end if
       nu_conv(itp1,i) = nu_new(i)
    end do
  end subroutine apply_outer_isothermal_cap

end module evolve_try_mod

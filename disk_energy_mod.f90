module disk_energy_mod
  !! Energy/viscosity update for a 1D viscous disk.
  !!
  !! This module is split into:
  !!  (1) solve_structure_from_sigma: side-effect-free solver for TRY step
  !!  (2) commit_structure_to_global: commit results to mod_global after acceptance
  !!
  use kind_params,  only : dp, i4b
  use constants,    only : kb, mp, pi, mu
  use mod_global,   only : M_star, q, R0, sigma_init, nu0_dim, nu0_nd, nu_conv, &
                           Tmid, H, rho, kappaR, tauR, Qvis, Qirr, dYdXi, Qrad, &
                           use_be_decretion, use_irradiation, is_shadow
  use units_disk_mod,       only : r_dim, omegaK_dim, sigma_dim
  use disk_thermal_mod,     only : heating_cooling_cell, build_shadow_flags
  use relax_mod,            only : ur_update
  use irradiation_mod, only : L_irr, Lirr_floor, rin_cgs, A1, L1, beta1, beta2, Q12
  use, intrinsic :: ieee_arithmetic
  implicit none

  private
  public :: build_initial_disk_structure
  public :: solve_structure_from_sigma
  public :: commit_structure_to_global

contains

  subroutine build_initial_disk_structure(nr, r_nd, sigma_in, itp1)
    use kind_params, only : dp, i4b
    use mod_global,  only : nu, nu0_dim, nu0_nd, Tmid
    implicit none

    integer(i4b), intent(in) :: nr, itp1
    real(dp),     intent(in) :: r_nd(nr), sigma_in(nr)

    ! ---- local work arrays (structure at itp1) ----
    real(dp) :: T_prev(nr)
    real(dp) :: nu_new(nr)          ! dimensionless nu (this will be updated)
    real(dp) :: Tmid_new(nr), H_loc(nr), rho_loc(nr), kappa_loc(nr), tau_loc(nr)
    real(dp) :: Qvis_loc(nr), Qrad_loc(nr), Qirr_prof(nr), dYdXi_loc(nr)
    logical  :: shadow_loc(nr)

    integer(i4b) :: mm_iter
    integer(i4b) :: i

    !------------------------------------------------------------
    ! 0) Prepare T_prev and initial guess nu_new
    !    - For itp1=1 in a fresh run, you may not have Tmid(0,:).
    !    - Use zeros or a floor as your convention, or reuse Tmid(1,:) if preset.
    !------------------------------------------------------------
    if (itp1 > 1) then
       T_prev(:) = Tmid(itp1-1, :)
    else
       T_prev(:) = 0.0_dp
    end if

    nu_new(:) = nu(:)   ! current dimensionless viscosity as initial guess

    !------------------------------------------------------------
    ! 1) Solve structure from Sigma (NO global commit here)
    !    Signature must match exactly what you pasted.
    !------------------------------------------------------------
    call solve_structure_from_sigma( nr, r_nd, sigma_in, T_prev, &
                                  nu_new, Tmid_new, H_loc, rho_loc, kappa_loc, tau_loc, &
                                  Qvis_loc, Qrad_loc, Qirr_prof, dYdXi_loc, shadow_loc, &
                                  mm_iter )

    !------------------------------------------------------------
    ! 2) Commit structure arrays to global at itp1
    !    (This writes Tmid(itp1,:), H(itp1,:), ... etc.)
    !------------------------------------------------------------
    call commit_structure_to_global( itp1, &
                                   Tmid_new, H_loc, rho_loc, kappa_loc, tau_loc, &
                                   Qvis_loc, Qrad_loc, Qirr_prof, dYdXi_loc, shadow_loc, &
                                   nu_new, mm_iter )

    !------------------------------------------------------------
    ! 3) Commit viscosity (dimensionless) to global nu(:)
    !    Here nu_new is already dimensionless, per your code convention.
    !------------------------------------------------------------
    nu(:) = nu_new(:)

  end subroutine build_initial_disk_structure

  subroutine solve_structure_from_sigma( nr, r_nd, sigma_new, T_prev, &
                                        nu_new, Tmid_new, H_loc, rho_loc, kappa_loc, tau_loc, &
                                        Qvis_loc, Qrad_loc, Qirr_prof, dYdXi_loc, shadow_loc, &
                                        mm_iter )
    !! Side-effect-free structure solver.
    !!
    !! Inputs:
    !!   sigma_new(:) : trial surface density (dimensionless)
    !!   T_prev(:)    : previous committed temperature used as root selector / UR seed
    !!
    !! Outputs:
    !!   nu_new(:)    : trial viscosity (dimensionless, consistent with nu0_nd)
    !!   (Tmid_new, H_loc, rho_loc, kappa_loc, tau_loc, Qvis_loc, Qrad_loc, Qirr_prof, dYdXi_loc, shadow_loc)
    !!
    integer(i4b), intent(in)  :: nr
    real(dp),     intent(in)  :: r_nd(nr), sigma_new(nr), T_prev(nr)
    real(dp),     intent(out) :: nu_new(nr)
    real(dp),     intent(out) :: Tmid_new(nr), H_loc(nr), rho_loc(nr), kappa_loc(nr), tau_loc(nr)
    real(dp),     intent(out) :: Qvis_loc(nr), Qrad_loc(nr), Qirr_prof(nr), dYdXi_loc(nr)
    logical,      intent(out) :: shadow_loc(nr)
    integer(i4b), intent(out) :: mm_iter

    real(dp) :: r_cgs(nr), Sigma_cgs(nr), OmegaK_cgs(nr)
    real(dp) :: nu_dim(nr)

    integer(i4b) :: i
    logical :: converged
    real(dp), parameter :: tiny_sigma = 1.0e-8_dp

    ! 1) Units + OmegaK
!$omp parallel do default(shared) private(i)
    do i = 1, nr
       r_cgs(i)     = r_dim(r_nd(i))
       Sigma_cgs(i) = sigma_dim(sigma_new(i))
       if (r_cgs(i) > 0.0_dp) then
          OmegaK_cgs(i) = omegaK_dim(r_nd(i))
       else
          OmegaK_cgs(i) = 0.0_dp
       end if
    end do
!$omp end parallel do

    ! Default outputs
    Tmid_new(:)  = 0.0_dp
    H_loc(:)     = 0.0_dp
    rho_loc(:)   = 0.0_dp
    kappa_loc(:) = 0.0_dp
    tau_loc(:)   = 0.0_dp
    Qvis_loc(:)  = 0.0_dp
    Qrad_loc(:)  = 0.0_dp
    Qirr_prof(:) = 0.0_dp
    dYdXi_loc(:) = 0.0_dp
    shadow_loc(:)= .true.
    nu_dim(:)    = 0.0_dp
    nu_new(:)    = 0.0_dp

    if (.not. use_irradiation .or. L_irr <= Lirr_floor) then
       shadow_loc(:) = .true.
       Qirr_prof(:)  = 0.0_dp
       dYdXi_loc(:)  = 0.0_dp

       call thermal_solve_stageA_try( nr, r_nd, r_cgs, Sigma_cgs, OmegaK_cgs, &
                                      shadow_loc, Qirr_prof, T_prev, &
                                      Tmid_new, H_loc, rho_loc, kappa_loc, tau_loc, nu_dim, &
                                      Qvis_loc, Qrad_loc )

    else
       ! Shadow flags computed from current geometry state.
       ! If you want TRY-consistent shadow, you must compute it from H_loc iteratively.
       call build_shadow_flags(nr, shadow_loc)

       call iterate_structure_irradiation_try( nr, r_nd, r_cgs, Sigma_cgs, OmegaK_cgs, &
                                              shadow_loc, rin_cgs, A1, L1, Q12, beta1, beta2, &
                                              T_prev, &
                                              Tmid_new, H_loc, rho_loc, kappa_loc, tau_loc, nu_dim, &
                                              Qvis_loc, Qrad_loc, Qirr_prof, dYdXi_loc, &
                                              mm_iter )
    end if

    ! Map nu_dim -> nu_new (dimensionless)
!$omp parallel do default(shared) private(i)
    do i = 1, nr
       if (Sigma_cgs(i) <= tiny_sigma) then
          ! Keep previous nu if available via nu_conv in caller (recommended),
          ! but here we just set 0; caller can overwrite if desired.
          nu_new(i) = 0.0_dp
       else
          if (nu0_dim > 0.0_dp) then
             nu_new(i) = (nu_dim(i) / nu0_dim) * nu0_nd
          else
             nu_new(i) = 0.0_dp
          end if
       end if
    end do
!$omp end parallel do
    !write(*,'("STRUCT: max Sigma_nd=",1pe12.4," max T=",1pe12.4," max nu=",1pe12.4," nTpos=",i8)') &
    !      maxval(sigma_new), maxval(Tmid_new), maxval(nu_new), count(Tmid_new > 0.0_dp)
  end subroutine solve_structure_from_sigma


  subroutine commit_structure_to_global( itp1, &
       Tmid_new, H_loc, rho_loc, kappa_loc, tau_loc, &
       Qvis_loc, Qrad_loc, Qirr_loc, dYdXi_loc, shadow_loc, &
       nu_new, mm_iter )
    use kind_params, only : dp, i4b
    use mod_global,  only : nr, Tmid, H, rho, kappaR, tauR, &
                            Qvis, Qrad, Qirr, dYdXi, is_shadow, &
                            nu_conv
    implicit none

    integer(i4b), intent(in) :: itp1, mm_iter
    real(dp), intent(in) :: Tmid_new(nr), H_loc(nr), rho_loc(nr)
    real(dp), intent(in) :: kappa_loc(nr), tau_loc(nr)
    real(dp), intent(in) :: Qvis_loc(nr), Qrad_loc(nr), Qirr_loc(nr)
    real(dp), intent(in) :: dYdXi_loc(nr)
    logical, intent(in)  :: shadow_loc(nr)
    real(dp), intent(in) :: nu_new(nr)

    integer(i4b) :: i

!$omp parallel do default(shared) private(i)
    do i = 1, nr
       Tmid(itp1, i)   = Tmid_new(i)
       H(itp1, i)      = H_loc(i)
       rho(itp1, i)    = rho_loc(i)
       kappaR(itp1, i) = kappa_loc(i)
       tauR(itp1, i)   = tau_loc(i)

       Qvis(itp1, i)   = Qvis_loc(i)
       Qrad(itp1, i)   = Qrad_loc(i)
       Qirr(itp1, i)   = Qirr_loc(i)

       dYdXi(itp1, i)  = dYdXi_loc(i)
       is_shadow(itp1,i)    = shadow_loc(i)

       nu_conv(itp1, i) = nu_new(i)

    end do
!$omp end parallel do

  end subroutine commit_structure_to_global

  !-----------------------------------------------------------------
  subroutine thermal_solve_stageA_try( nr, r_nd, r_cgs, Sigma_cgs, OmegaK_cgs, shadow, &
                                      Qirr_prof, T_prev, &
                                      Tmid_new, H_loc, rho_loc, kappa_loc, tau_loc, nu_dim, &
                                      Qvis_loc, Qrad_loc )
    !! Thermal balance solve for each cell (no global reference).
    use radiation_params_mod, only: T_floor
    implicit none
    integer(i4b), intent(in) :: nr
    real(dp), intent(in)     :: r_nd(nr), r_cgs(nr), Sigma_cgs(nr), OmegaK_cgs(nr)
    logical, intent(in)      :: shadow(nr)
    real(dp), intent(in)     :: Qirr_prof(nr)
    real(dp), intent(in)     :: T_prev(nr)

    real(dp), intent(out) :: Tmid_new(nr), H_loc(nr), rho_loc(nr), kappa_loc(nr), tau_loc(nr), nu_dim(nr)
    real(dp), intent(out) :: Qvis_loc(nr), Qrad_loc(nr)

    real(dp) :: Tmid_roots(20)
    real(dp) :: T_old, T_root, logT_new
    real(dp) :: Qplus_visc_dummy, Qplus_irr_dummy, Qminus_dummy
    integer(i4b) :: nroots, ierr, i
    real(dp), parameter :: tiny_sigma = 1.0e-8_dp
    real(dp), parameter :: omega_T = 0.3_dp

    nu_dim(:)   = 0.0_dp
    Qvis_loc(:) = 0.0_dp
    Qrad_loc(:) = 0.0_dp

!$omp parallel do default(shared) private(i, Tmid_roots, T_old, T_root, nroots, ierr, &
!$omp& Qplus_visc_dummy, Qplus_irr_dummy, Qminus_dummy, logT_new)
    do i = 1, nr
       Tmid_new(i)  = 0.0_dp
       H_loc(i)     = 0.0_dp
       rho_loc(i)   = 0.0_dp
       kappa_loc(i) = 0.0_dp
       tau_loc(i)   = 0.0_dp

       if (Sigma_cgs(i) <= tiny_sigma) cycle

       T_old = max(0.0_dp, T_prev(i))

       call solve_Tmid_cell( r_cgs(i), Sigma_cgs(i), OmegaK_cgs(i), shadow(i), &
                             T_old, Tmid_roots, H_loc(i), rho_loc(i), nu_dim(i), &
                             kappa_loc(i), tau_loc(i), nroots, ierr, &
                             Qvis_loc(i), Qirr_prof(i), Qrad_loc(i) )

       if (ierr == 0 .and. nroots > 0 .and. Tmid_roots(1) > 0.0_dp) then
          T_root = Tmid_roots(1)

          !if (T_old > 0.0_dp) then
          !   logT_new    = ur_update(log(T_old), log(T_root), omega_T)
          !   Tmid_new(i) = exp(logT_new)
          !else
             Tmid_new(i) = T_root
          !end if

          call heating_cooling_cell( r_cgs(i), Sigma_cgs(i), OmegaK_cgs(i), shadow(i), &
                                     Tmid_new(i), H_loc(i), rho_loc(i), nu_dim(i), &
                                     kappa_loc(i), tau_loc(i), &
                                     Qplus_visc_dummy, Qplus_irr_dummy, Qminus_dummy, &
                                     Qirr_in=Qirr_prof(i) )
       else
          ! If no root, keep zeros. Caller will likely reject by accept_step_try.
         write (*, '("+++ i =", i4, ": No root was found. Stop. +++")') i
         stop
       end if
    end do
!$omp end parallel do

  end subroutine thermal_solve_stageA_try


  !-----------------------------------------------------------------
  subroutine iterate_structure_irradiation_try( nr, r_nd, r_cgs, Sigma_cgs, OmegaK_cgs, &
                                               shadow, rin_cgs_in, A1_in, L1_in, Q12_in, beta1_in, beta2_in, &
                                               T_prev, &
                                               Tmid_new, H_loc, rho_loc, kappa_loc, tau_loc, nu_dim, &
                                               Qvis_loc, Qrad_loc, Qirr_prof, dYdXi_loc, &
                                               mm_iter )
    !! TRY-safe wrapper around your existing fixed-point irradiation iteration.
    !!
    !! IMPORTANT:
    !! - This routine must NOT write to mod_global arrays.
    !! - If your existing iterate_structure_irradiation writes globals, you must
    !!   create a "_try" variant that returns arrays instead.
    !!
    !! Here we assume you already have a side-effect-free implementation
    !! (or you will refactor it similarly).
    implicit none
    integer(i4b), intent(in) :: nr
    real(dp), intent(in) :: r_nd(nr), r_cgs(nr), Sigma_cgs(nr), OmegaK_cgs(nr)
    logical,  intent(in) :: shadow(nr)
    real(dp), intent(in) :: rin_cgs_in, A1_in, L1_in, Q12_in, beta1_in, beta2_in
    real(dp), intent(in) :: T_prev(nr)

    real(dp), intent(out) :: Tmid_new(nr), H_loc(nr), rho_loc(nr), kappa_loc(nr), tau_loc(nr), nu_dim(nr)
    real(dp), intent(out) :: Qvis_loc(nr), Qrad_loc(nr), Qirr_prof(nr), dYdXi_loc(nr)
    integer(i4b), intent(out) :: mm_iter

    integer(i4b), parameter :: mmax = 20
    real(dp),     parameter :: epsQ = 1.0e-3_dp
    logical :: converged

    call iterate_structure_irradiation( nr, r_nd, r_cgs, Sigma_cgs, OmegaK_cgs, &
                                        0, shadow, rin_cgs_in, A1_in, L1_in, Q12_in, &
                                        beta1_in, beta2_in, mmax, epsQ, mm_iter, &
                                        Tmid_new, H_loc, rho_loc, kappa_loc, tau_loc, nu_dim, &
                                        Qirr_prof, dYdXi_loc )

    ! Qvis/Qrad are not returned by your current iterate_structure_irradiation interface.
    ! If you want accept_step to check thermal balance, you must compute them consistently
    ! for the final (Tmid_new, ...) solution:
    call build_Qvis_Qrad_from_solution( nr, r_cgs, Sigma_cgs, OmegaK_cgs, shadow, &
                                        Tmid_new, Qirr_prof, &
                                        Qvis_loc, Qrad_loc )

  end subroutine iterate_structure_irradiation_try

  !-----------------------------------------------------------------
  subroutine iterate_structure_irradiation( nr, r_nd, r_cgs, Sigma_cgs, OmegaK_cgs,  &
                               itp1, is_shadow,                        &
                               rin_cgs, A1, L1, Q12, beta1, beta2,     &
                               mmax, epsQ, mm_iter,                    &
                               Tmid_new, H_loc, rho_loc, kappa_loc, tau_loc, nu_dim, &
                               Qirr_new, dYdXi_loc)
    use kind_params, only : dp, i4b
    implicit none
    integer(i4b), intent(in) :: nr, itp1, mmax
    real(dp),     intent(in) :: r_nd(nr), r_cgs(nr), Sigma_cgs(nr), OmegaK_cgs(nr)
    logical,      intent(in) :: is_shadow(nr)
    real(dp),     intent(in) :: rin_cgs, A1, L1, Q12, beta1, beta2, epsQ
    real(dp),     intent(out):: Tmid_new(nr), H_loc(nr), rho_loc(nr), kappa_loc(nr), tau_loc(nr)
    real(dp),     intent(out):: nu_dim(nr), Qirr_new(nr), dYdXi_loc(nr)
    !real(dp),     intent(out):: Qvis_loc(nr), Qrad_loc(nr)
    integer(i4b), intent(out) :: mm_iter

    real(dp) :: Qirr_old(nr), Qirr_in(nr), Qirr_calc(nr), Y(nr)
    integer(i4b) :: m, i
    real(dp) :: Tmid_old(nr)
    real(dp), parameter :: tinyQ = 1.0e-99_dp
    real(dp), parameter :: tinyT = 1.0e-99_dp
    real(dp) :: num, den, relQ, relT, rel
    real(dp), parameter :: tiny_sigma = 1.0e-8_dp
    real(dp) :: rel_prev, rel_best
    real(dp) :: wloc
    real(dp) :: rQ
    real(dp) :: dotp, n1, n2, costh
    real(dp) :: Qirr_prev(nr), Qirr_prev2(nr)

    integer(i4b) :: n_stall, n_blow
    integer(i4b), parameter :: stall_max = 6     ! how many non-improving iters allowed
    integer(i4b), parameter :: blow_max  = 3     ! how many increasing iters allowed
    real(dp), parameter :: w_init = 0.3_dp
    real(dp), parameter :: w_min  = 0.02_dp
    real(dp), parameter :: w_shrink = 0.5_dp    ! w <- w*w_shrink when needed
    real(dp), parameter :: grow_tol = 1.05_dp   ! consider "growing" if rel > 1.05*rel_prev
    real(dp), parameter :: imp_tol  = 0.98_dp   ! consider "improving" if rel < 0.98*rel_prev

    write(*,'(A,I0,A,I0)') 'ENTER iterate_structure_irradiation: itp1=', itp1, ' mmax=', mmax
    ! Initial guess for stageA input
    if (itp1 > 1) then
       Qirr_old(:) = Qirr(itp1-1, :)
    else
       Qirr_old(:) = 0.0_dp
    end if

    ! First structure solve with the initial irradiation guess
    call thermal_solve_stageA(nr, r_nd, r_cgs, Sigma_cgs, OmegaK_cgs, is_shadow, &
                          Qirr_old, itp1, &
                          Tmid_new, H_loc, rho_loc, kappa_loc, tau_loc, nu_dim)

    wloc     = w_init
    rel_prev = huge(1.0_dp)
    rel_best = huge(1.0_dp)
    n_stall  = 0
    n_blow   = 0
    Qirr_prev(:)  = Qirr_old(:)
    Qirr_prev2(:) = Qirr_old(:)

    do m = 1, mmax

       Tmid_old(:) = Tmid_new(:)

       call irradiation_stageB_eq16(nr, r_nd, r_cgs, H_loc, &
                                rin_cgs, A1, L1, Q12, beta1, beta2, &
                                Y, dYdXi_loc, Qirr_calc)

       ! Under-relaxation with adaptive weight
       Qirr_in(:) = (1.0_dp - wloc)*Qirr_old(:) + wloc*Qirr_calc(:)

       call thermal_solve_stageA(nr, r_nd, r_cgs, Sigma_cgs, OmegaK_cgs, is_shadow, &
                             Qirr_in, itp1, &
                             Tmid_new, H_loc, rho_loc, kappa_loc, tau_loc, nu_dim)

       !------------------------------------------------------------
       ! Convergence check: combined relative change in Qirr_in and Tmid
       ! relQ = ||Qirr_in - Qirr_old|| / ||Qirr_in||
       ! relT = ||Tmid_new - Tmid_old|| / ||Tmid_new||
       ! rel  = max(relQ, relT)
       !------------------------------------------------------------
       ! We compute RMS-type relative norms over cells with meaningful Sigma
       num = 0.0_dp; den = 0.0_dp
!$omp parallel do default(shared) private(i) reduction(+:num,den)
       do i = 1, nr
          if (Sigma_cgs(i) <= tiny_sigma) cycle
          num = num + (Qirr_in(i) - Qirr_old(i))**2
          den = den + max(abs(Qirr_in(i)), tinyQ)**2
       end do
!$omp end parallel do
       relQ = sqrt(num / max(den, tinyQ))

       ! --- true fixed-point residual for irradiation (not step size) ---
       num = 0.0_dp; den = 0.0_dp
!$omp parallel do default(shared) private(i) reduction(+:num,den)
       do i = 1, nr
          if (Sigma_cgs(i) <= tiny_sigma) cycle
          num = num + (Qirr_calc(i) - Qirr_in(i))**2
          den = den + max(abs(Qirr_in(i)), tinyQ)**2
       end do
!$omp end parallel do
       rQ = sqrt(num / max(den, tinyQ))

       ! --- oscillation detector using successive update directions ---
       dotp = 0.0_dp; n1 = 0.0_dp; n2 = 0.0_dp
!$omp parallel do default(shared) private(i) reduction(+:dotp,n1,n2)
       do i = 1, nr
          if (Sigma_cgs(i) <= tiny_sigma) cycle
          dotp = dotp + (Qirr_in(i) - Qirr_prev(i)) * (Qirr_prev(i) - Qirr_prev2(i))
          n1   = n1   + (Qirr_in(i) - Qirr_prev(i))**2
          n2   = n2   + (Qirr_prev(i) - Qirr_prev2(i))**2
       end do
!$omp end parallel do
       costh = dotp / sqrt(max(n1, tinyQ)*max(n2, tinyQ))

       num = 0.0_dp; den = 0.0_dp
!$omp parallel do default(shared) private(i) reduction(+:num,den)
       do i = 1, nr
          if (Sigma_cgs(i) <= tiny_sigma) cycle
          num = num + (Tmid_new(i) - Tmid_old(i))**2
          den = den + max(abs(Tmid_new(i)), tinyT)**2
       end do
!$omp end parallel do
       relT = sqrt(num / max(den, tinyT))

       rel = max(relQ, relT)

       write (*, '("m =", i3, ": relQ =", 1pe12.4, ", relT =", 1pe12.4, &
             ", rel =", 1pe12.4, ", rQ =", 1pe12.4, ", costh =", 1pe12.4, &
             ", wloc =", 1pe12.4)') m, relQ, relT, rel, rQ, costh, wloc

       Qirr_prev2(:) = Qirr_prev(:)
       Qirr_prev(:)  = Qirr_in(:)
       Qirr_old(:) = Qirr_in(:)

       if (rel < rel_best) rel_best = rel

       if (rel < epsQ) then
          mm_iter = m
          write(*,'(A,I0,A,1PE12.4)') 'CONVERGED itp1=', itp1, ' rel=', rel
          exit
       end if

       ! --- convergence diagnostics / adaptation ---
       if (rel > grow_tol*rel_prev) then
          n_blow = n_blow + 1
          ! If it's growing, shrink relaxation weight
          if (wloc > w_min) wloc = max(w_min, wloc*w_shrink)
       else
          n_blow = 0
       end if

       if (rel < imp_tol*rel_prev) then
          n_stall = 0
       else
          n_stall = n_stall + 1
       end if

       ! Early exit if it clearly stalls or blows up
       if (n_blow >= blow_max .or. n_stall >= stall_max) then
          mm_iter = -1   ! non-converged, early-aborted
          write(*,'(A,I0,A,I0,A,1PE12.4)') 'EARLY-ABORT itp1=', itp1, ' m=', m, ' rel=', rel
          exit
       end if

       rel_prev = rel

       if (m == mmax) mm_iter = -99
       write(*,'(A,I0,A,I0,A,1PE12.4)') 'ABORT itp1=', itp1, ' m=', m, ' rel=', rel
    end do

    Qirr_new(:) = Qirr_old(:)
    
  end subroutine iterate_structure_irradiation

  subroutine build_Qvis_Qrad_from_solution( nr, r_cgs, Sigma_cgs, OmegaK_cgs, shadow, &
                                            Tmid_new, Qirr_prof, Qvis_loc, Qrad_loc )
    !! Compute Qvis/Qrad consistently for acceptance checks.
    implicit none
    integer(i4b), intent(in) :: nr
    real(dp), intent(in) :: r_cgs(nr), Sigma_cgs(nr), OmegaK_cgs(nr)
    logical, intent(in)  :: shadow(nr)
    real(dp), intent(in) :: Tmid_new(nr), Qirr_prof(nr)
    real(dp), intent(out):: Qvis_loc(nr), Qrad_loc(nr)

    integer(i4b) :: i
    real(dp) :: H_tmp, rho_tmp, nu_tmp, kappa_tmp, tau_tmp
    real(dp) :: Qvis_tmp, Qirr_tmp, Qrad_tmp

!$omp parallel do default(shared) private(i, H_tmp, rho_tmp, nu_tmp, kappa_tmp, tau_tmp, Qvis_tmp, Qirr_tmp, Qrad_tmp)
    do i = 1, nr
       Qvis_loc(i) = 0.0_dp
       Qrad_loc(i) = 0.0_dp
       if (Sigma_cgs(i) <= 0.0_dp) cycle
       if (OmegaK_cgs(i) <= 0.0_dp) cycle
       if (Tmid_new(i) <= 0.0_dp) cycle

       call heating_cooling_cell( r_cgs(i), Sigma_cgs(i), OmegaK_cgs(i), shadow(i), Tmid_new(i), &
                                  H_tmp, rho_tmp, nu_tmp, kappa_tmp, tau_tmp, &
                                  Qvis_tmp, Qirr_tmp, Qrad_tmp, Qirr_in=Qirr_prof(i) )
       Qvis_loc(i) = Qvis_tmp
       Qrad_loc(i) = Qrad_tmp
    end do
!$omp end parallel do

  end subroutine build_Qvis_Qrad_from_solution

  !-----------------------------------------------------------------
  subroutine irradiation_stageB_eq16(nr, xi, r_cgs, H_cgs, &
                                     rin_cgs, A1, L1, Q12, beta1, beta2, &
                                     Y, dYdXi, Qirr_prof)
    use kind_params, only: dp, i4b
    use irradiation_mod, only: compute_Y_dYdXi, compute_Y_dYdXi_sg, compute_Qirr_eq16
    implicit none
    integer(i4b), intent(in) :: nr
    real(dp), intent(in) :: xi(nr), r_cgs(nr), H_cgs(nr)
    real(dp), intent(in) :: rin_cgs, A1, L1, Q12, beta1, beta2
    real(dp), intent(out):: Y(nr), dYdXi(nr), Qirr_prof(nr)
    integer(i4b) :: halfwin, poly_order

    ! No smoothing of Y and dY/dXi
    !call compute_Y_dYdXi(nr, xi, r_cgs, H_cgs, Y, dYdXi)
    ! Smoothing of Y and dYdXi by Savitzky–Golay method
    call compute_Y_dYdXi_sg(n=nr, xi=xi, r_cgs=r_cgs, H_cgs=H_cgs, &
                        halfwin=7, poly_order=2, Y=Y, dYdXi=dYdXi)
    call compute_Qirr_eq16(n=nr, xi=xi, Y=Y, dYdXi=dYdXi, rin_cgs=rin_cgs, &
                        A1=A1, L1=L1, Q12=Q12, beta1=beta1, beta2=beta2, &
                        Qirr_out=Qirr_prof)
    !call compute_Y_dYdXi_sg(nr, xi, r_cgs, H_cgs, halfwin, poly_order, Y, dYdXi)
    !call compute_Qirr_eq16(nr, xi, Y, dYdXi, rin_cgs, A1, L1, Q12, beta1, beta2, Qirr_prof)

  end subroutine irradiation_stageB_eq16

  !-----------------------------------------------------------------
  subroutine thermal_solve_stageA( nr, r_nd, r_cgs, Sigma_cgs, OmegaK_cgs, is_shadow, &
                               Qirr_prof, itp1, &
                               Tmid_new, H_loc, rho_loc, kappa_loc, tau_loc, nu_dim)
    use kind_params, only: dp, i4b
    use radiation_params_mod, only: T_floor
    implicit none
    integer(i4b), intent(in) :: nr, itp1
    real(dp), intent(in) :: r_nd(nr), r_cgs(nr), Sigma_cgs(nr), OmegaK_cgs(nr)
    logical, intent(in) :: is_shadow(nr)
    real(dp), intent(in) :: Qirr_prof(nr)

    real(dp), intent(out) :: Tmid_new(nr), H_loc(nr), rho_loc(nr), &
                             kappa_loc(nr), tau_loc(nr), nu_dim(nr)
    real(dp) :: Qvis_loc(nr), Qrad_loc(nr)

    real(dp) :: Tmid_roots(20), T_old, T_root, logT_new
    real(dp) :: Qplus_visc_dummy, Qplus_irr_dummy, Qminus_dummy
    integer(i4b) :: nroots, ierr, i
    real(dp), parameter :: tiny_sigma = 1.0e-8_dp
    real(dp), parameter :: omega_T = 0.3_dp   ! start conservative: 0.05-0.3

    !write (*, '("thermal_solve_stageA: Qirr_prof(1) =", 1pe12.4)') Qirr_prof(1)
    !write (*, '("thermal_solve_stageA: Qirr_prof(nr) =", 1pe12.4)') Qirr_prof(nr)

    nu_dim(:) = 0.0_dp

!$omp parallel do default(shared) private(i, Tmid_roots, &
!$omp& T_old, T_root, nroots, ierr, &
!$omp& Qplus_visc_dummy, Qplus_irr_dummy, Qminus_dummy)
    do i = 1, nr
       Tmid_new(i)  = 0.0_dp
       H_loc(i)     = 0.0_dp
       rho_loc(i)   = 0.0_dp
       kappa_loc(i) = 0.0_dp
       tau_loc(i)   = 0.0_dp

       if (Sigma_cgs(i) <= tiny_sigma) cycle

       T_old = 0.0_dp
       if (itp1 > 1) T_old = Tmid(itp1-1, i)

       !write (*, '("Calling solve_Tmid_cell: Sigma_cgs =", 1pe12.4, " at i =", i4)') &
       !       Sigma_cgs(i), i
       call solve_Tmid_cell( r_cgs(i), Sigma_cgs(i), OmegaK_cgs(i), &
                            is_shadow(i), T_old, Tmid_roots,       &
                            H_loc(i), rho_loc(i), nu_dim(i),       &
                            kappa_loc(i), tau_loc(i),              &
                            nroots, ierr, Qvis_loc(i), Qirr_prof(i), Qrad_loc(i) )

       if (ierr == 0 .and. nroots > 0 .and. Tmid_roots(1) > 0.0_dp) then
          T_root       = Tmid_roots(1)
          Tmid_new(i)  = max(T_floor, T_root)
          !if (itp1 > 1 .and. T_old > 0.0_dp) then
          !   ! Under-relax temperature across time steps
          !   logT_new = ur_update(log(T_old), log(T_root), omega_T)
          !   Tmid_new(i) = exp(logT_new)
          !else
          !   Tmid_new(i) = T_root
          !end if

          ! Recompute consistent structure with external Qirr (LOH24 Eq.16)
          call heating_cooling_cell( r_cgs(i), Sigma_cgs(i), OmegaK_cgs(i), is_shadow(i), &
                                   Tmid_new(i), H_loc(i), rho_loc(i), nu_dim(i),        &
                                   kappa_loc(i), tau_loc(i),                            &
                                   Qplus_visc_dummy, Qplus_irr_dummy, Qminus_dummy,     &
                                   Qirr_in=Qirr_prof(i) )
       else
          ! Fallback: use T_floor and compute consistent structure
          Tmid_new(i) = T_floor
          call heating_cooling_cell( r_cgs(i), Sigma_cgs(i), OmegaK_cgs(i), is_shadow(i), &
                                   Tmid_new(i), H_loc(i), rho_loc(i), nu_dim(i),        &
                                   kappa_loc(i), tau_loc(i),                            &
                                   Qplus_visc_dummy, Qplus_irr_dummy, Qminus_dummy,     &
                                   Qirr_in=Qirr_prof(i) )
          ! Warn only once per cell per run (itp1<=1 avoids spam in iteration)
          ! if (itp1 <= 1) write (*, '("Warning: No thermal root at i =", i4)') i
       end if
    end do
!$omp end parallel do

  end subroutine thermal_solve_stageA

  !-----------------------------------------------------------------
  subroutine solve_Tmid_cell(r_cgs, Sigma_cgs, OmegaK_cgs, shadow, &
                           Tmid_old, Tmid_roots, H_loc, rho_loc, &
                           nu_dim, kappa_loc, tau_loc, nroots, ierr, &
                           Qvis_tmp, Qirr_in, Qrad_tmp )
    !! Wrapper around heating_cooling_cell(...) to find thermal-balance roots.
    use ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
    use radiation_params_mod, only: T_floor, T_ceiling
    use opacity_table_mod, only : thigh_tab, tlow_tab, logRmin_tab, rhomin_tab
    implicit none

    real(dp), intent(in)  :: r_cgs, Sigma_cgs, OmegaK_cgs
    logical,  intent(in)  :: shadow
    real(dp), intent(in)  :: Tmid_old
    real(dp), intent(in)  :: Qirr_in
    !integer(i4b), intent(in) :: itp1

    real(dp), intent(out) :: Tmid_roots(20)
    real(dp), intent(out) :: H_loc, rho_loc, kappa_loc, tau_loc, nu_dim
    real(dp), intent(out) :: Qvis_tmp, Qrad_tmp
    integer(i4b), intent(out) :: nroots, ierr

    integer(i4b) :: Nscan
    integer(i4b), parameter :: iter_max = 30
    integer(i4b), parameter :: Nscan0 = 200
    integer(i4b), parameter :: scan_max = 3
    real(dp), parameter :: fTrange = 1.7_dp
    real(dp),    parameter :: tol_rel  = 1.0e-6_dp

    integer(i4b) :: i, j, jbest, iter, iL, scan_loop
    integer(i4b), allocatable :: iCand(:), ieee_ok(:)
    real(dp),     allocatable :: logTgrid(:), Fgrid(:)

    real(dp) :: logT_lo, logT_hi, dlogT
    real(dp) :: f_lo, f_hi, f_mid, logT_mid
    real(dp) :: dmin
    real(dp) :: T_trial, H_trial, rho_trial, cs2, logR_trial, &
                H_tmp, rho_tmp, kappa_tmp, nu_tmp
    real(dp) :: tauR_tmp, Qirr_tmp
    logical, parameter :: debug_on = .false.
    logical, save :: first = .true.
    integer(i4b), save :: iu_ro = -1
    character(len=70), parameter :: file_Troots = 'T_roots.txt'


    ierr       = 0
    nroots     = 0
    Tmid_roots = 0.0_dp

    if (Sigma_cgs <= 0.0_dp .or. OmegaK_cgs <= 0.0_dp) then
       ierr = -1
       return
    end if

    !write (*, '("solve_Tmid_cell: Qirr_in =", 1pe12.4)') Qirr_in
    
    do scan_loop = 1, scan_max
       !------------------------------------------
       ! Choose scan range
       !------------------------------------------
       !if (Tmid_old > 0.0_dp .and. scan_loop == 1) then
       if (Tmid_old > 0.0_dp) then
          logT_lo = log10(max(T_floor,   Tmid_old / fTrange**scan_loop))
          logT_hi = log10(min(T_ceiling, Tmid_old * fTrange**scan_loop))
          if (logT_hi <= logT_lo + 1.0e-6_dp) then
             logT_lo = log10(T_floor)
             logT_hi = log10(T_ceiling)
          end if
          Nscan = Nscan0 * scan_loop
       else
          logT_lo = log10(T_floor)
          logT_hi = log10(T_ceiling)
          !logT_lo = log10(max(T_floor,   Tmid_old / fTrange**scan_loop))
          !logT_hi = log10(min(T_ceiling, Tmid_old))
          Nscan = Nscan0 * scan_max
       end if

       dlogT = (logT_hi - logT_lo) / real(Nscan, dp)

       allocate(logTgrid(Nscan+1), Fgrid(Nscan+1), ieee_ok(Nscan+1))
       ieee_ok(:) = 0

       ! Scan (high -> low) is fine, but keep consistent
       do i = 1, Nscan + 1
          ! logTgrid has to be in the ascending order
          logTgrid(i) = logT_lo + dlogT*real(i-1,dp)
          T_trial     = 10.0_dp**(logTgrid(i))
          if (T_trial > thigh_tab .or. T_trial < tlow_tab) then
             ieee_ok(i) = 1
             cycle
          else
             cs2     = kb * T_trial / (mu * mp)
             H_trial   = sqrt(cs2) / OmegaK_cgs
             rho_trial = Sigma_cgs / (2.0_dp * H_trial)
             logR_trial = log10(rho_trial) - 3.0_dp*log10(T_trial) + 18.0_dp
             if (logR_trial < logRmin_tab .or. rho_trial < rhomin_tab) then
                ieee_ok(i) = 1
                cycle
             end if
          end if

          call heating_cooling_cell( r_cgs, Sigma_cgs, OmegaK_cgs, shadow, T_trial, &
                                H_tmp, rho_tmp, nu_tmp, kappa_tmp, tauR_tmp,   &
                                Qvis_tmp, Qirr_tmp, Qrad_tmp,                  &
                                Qirr_in=Qirr_in )

          Fgrid(i) = (Qvis_tmp + Qirr_tmp) - Qrad_tmp

          if (ieee_is_nan(Fgrid(i)) .or. .not. ieee_is_finite(Fgrid(i))) then
             ieee_ok(i) = 1
          end if
       end do

       allocate(iCand(Nscan))
       nroots = 0
       do i = 1, Nscan
          if (ieee_ok(i) /= 0 .or. ieee_ok(i+1) /= 0) cycle
          if (Fgrid(i) * Fgrid(i+1) < 0.0_dp) then
             if (nroots < size(Tmid_roots)) then
                nroots = nroots + 1
                iCand(nroots) = i
             end if
          end if
       end do

       if (nroots > 0) then
          ierr = 0
          exit
       else
          ierr = 1
          deallocate(logTgrid, Fgrid, ieee_ok, iCand)
       end if
    end do

    if (nroots == 0) return

    ! Refine each bracket by bisection in logT
    do j = 1, nroots
       iL = iCand(j)
       logT_lo = logTgrid(iL)
       logT_hi = logTgrid(iL+1)

       T_trial = 10.0_dp**(logT_lo)
       call heating_cooling_cell( r_cgs, Sigma_cgs, OmegaK_cgs, shadow, T_trial, &
                                H_tmp, rho_tmp, nu_tmp, kappa_tmp, tauR_tmp,   &
                                Qvis_tmp, Qirr_tmp, Qrad_tmp,                  &
                                Qirr_in=Qirr_in )
       f_lo = (Qvis_tmp + Qirr_tmp) - Qrad_tmp

       T_trial = 10.0_dp**(logT_hi)
       call heating_cooling_cell( r_cgs, Sigma_cgs, OmegaK_cgs, shadow, T_trial, &
                                H_tmp, rho_tmp, nu_tmp, kappa_tmp, tauR_tmp,   &
                                Qvis_tmp, Qirr_tmp, Qrad_tmp,                  &
                                Qirr_in=Qirr_in )
       f_hi = (Qvis_tmp + Qirr_tmp) - Qrad_tmp

       do iter = 1, iter_max
          logT_mid = 0.5_dp * (logT_lo + logT_hi)
          T_trial  = 10.0_dp**(logT_mid)

          call heating_cooling_cell( r_cgs, Sigma_cgs, OmegaK_cgs, shadow, T_trial, &
                                   H_tmp, rho_tmp, nu_tmp, kappa_tmp, tauR_tmp,   &
                                   Qvis_tmp, Qirr_tmp, Qrad_tmp,                  &
                                   Qirr_in=Qirr_in )
          f_mid = (Qvis_tmp + Qirr_tmp) - Qrad_tmp

          if (f_lo * f_mid <= 0.0_dp) then
             logT_hi = logT_mid
             f_hi    = f_mid
          else
             logT_lo = logT_mid
             f_lo    = f_mid
          end if

          if (abs(logT_hi - logT_lo) / max(abs(logT_hi), 1.0e-99_dp) < tol_rel) then
             !if (j <= size(Tmid_roots)) Tmid_roots(j) = 10.0_dp**(0.5_dp*(logT_lo + logT_hi))
             Tmid_roots(j) = 10.0_dp**(0.5_dp*(logT_lo + logT_hi))
             exit
          end if
       end do

       ! --- after Tmid_roots(j) is set for each j
       if (debug_on) then
          if (r_cgs > 1.6e12_dp .and. r_cgs < 2.0e12_dp) then
             T_trial = Tmid_roots(j)
             call heating_cooling_cell(r_cgs, Sigma_cgs, OmegaK_cgs, shadow, T_trial, &
                  H_tmp, rho_tmp, nu_tmp, kappa_tmp, tauR_tmp, Qvis_tmp, Qirr_tmp, Qrad_tmp, &
                  Qirr_in=Qirr_in)
             write(*,'(a,1pe12.4,a,i2,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4)') &
               'root: T=',T_trial,' j=',j,' kappa=',kappa_tmp,' tau=',tauR_tmp, &
               ' Qvis=',Qvis_tmp,' Qirr=',Qirr_tmp,' Qrad=',Qrad_tmp, &
               ' F=',(Qvis_tmp+Qirr_tmp-Qrad_tmp)
          end if
       end if
    end do

    !-------------------------
    if (r_cgs > 2.15*R0 .and. r_cgs < 2.20*R0) then
       if (first) then
          open(newunit=iu_ro, file=trim(file_Troots), status='replace', &
               action='write')
          write (iu_ro, '("#(roots)   Tmid''s")')
          first = .false.
       end if
       write (iu_ro, '(i7, 2x, 10(1pe12.4,:))') nroots, &
             (Tmid_roots(j), j=1, nroots)
    end if
    !-------------------------

    deallocate(logTgrid, Fgrid, ieee_ok, iCand)

    ! Choose best root
    if (Tmid_old > 0.0_dp) then
       dmin  = 1.0e99_dp
       jbest = 1
       !do j = 1, min(nroots, size(Tmid_roots))
       do j = 1, nroots
          if (Tmid_roots(j) <= 0.0_dp) cycle
          if (abs(Tmid_roots(j) - Tmid_old) < dmin) then
             dmin  = abs(Tmid_roots(j) - Tmid_old)
             jbest = j
          end if
       end do
       Tmid_roots(1) = Tmid_roots(jbest)
       nroots        = 1
    else
       Tmid_roots(1) = minval(Tmid_roots(1:nroots))
       !Tmid_roots(1) = maxval(Tmid_roots(1:nroots))
       nroots        = 1
    end if
    if (r_cgs > 2.15*R0 .and. r_cgs < 2.20*R0) then
       write (iu_ro, '(8x, "-> Tmid(selected) =", 1pe12.4)') Tmid_roots(1)
    end if

    ! Return consistent structure at selected root (including external irradiation)
    T_trial = Tmid_roots(1)
    call heating_cooling_cell( r_cgs, Sigma_cgs, OmegaK_cgs, shadow, T_trial, &
                             H_loc, rho_loc, nu_dim, kappa_loc, tau_loc,     &
                             Qvis_tmp, Qirr_tmp, Qrad_tmp,                  &
                             Qirr_in=Qirr_in )

    ierr = 0
  end subroutine solve_Tmid_cell

end module disk_energy_mod

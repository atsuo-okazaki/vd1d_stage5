module evolve_substep_mod
  use kind_params, only : dp, i4b
  use mod_global,  only : nr, r, &
                         sigma_cur, nu_cur, Tmid_cur, &
                         H_cur, rho_cur, kappa_cur, tau_cur, &
                         Qvis_cur, Qrad_cur, Qirr_cur, dYdXi_cur, shadow_cur, &
                         use_irradiation
  use evolve_try_mod, only : accept_step_try, evolve_physics_one_substep_try
  use irradiation_mod, only : L_irr, Lirr_floor, set_zero_irradiation, &
                              set_irradiation_luminosity_from_arrays, &
                              build_Qirr_profile_eq16
  implicit none
contains

  logical function substep_try_and_commit(it_out, t_local, dt_nd, source_n, source_np1) result(ok)
    integer(i4b), intent(in) :: it_out
    real(dp),     intent(in) :: t_local, dt_nd
    real(dp),     intent(in) :: source_n(nr), source_np1(nr)

    ! Trial arrays
    real(dp) :: sigma_try(nr), nu_try(nr)
    real(dp) :: T_try(nr), H_try(nr), rho_try(nr), kappa_try(nr), tau_try(nr)
    real(dp) :: Qvis_try(nr), Qrad_try(nr), Qirr_try(nr), dYdXi_try(nr)
    logical  :: shadow_try(nr)
    integer(i4b) :: kiter_try, miter_try

    ! Irradiation profile for THIS substep
    real(dp) :: Qirr_prof(nr), dYdXi_prof(nr)
    logical  :: shadow_prof(nr)

    !--------------------------------------------
    ! Build irradiation profile at SUBSTEP START
    !--------------------------------------------
    if (.not. use_irradiation .or. L_irr <= Lirr_floor) then
    !if (.not. use_irradiation) then
       call set_zero_irradiation()
       Qirr_prof(:)  = 0.0_dp
       dYdXi_prof(:) = 0.0_dp
       shadow_prof(:)= .false.
    else
       ! shadowing not implemented yet
       shadow_prof(:)= .false.

       ! Update L_irr and LOH24 params using current state arrays
       call set_irradiation_luminosity_from_arrays(t_nd_now=t_local, sigma_nd=sigma_cur, nu_nd=nu_cur)

       ! Build Qirr(r) from CURRENT geometry (H_cur is already CGS in your case)
       call build_Qirr_profile_eq16(nr, r, H_cur, Qirr_prof, dYdXi_prof)
    end if

    !--------------------------------------------
    ! 1) Build trial step from current state
    !    IMPORTANT: pass profiles as INPUT
    !--------------------------------------------
    call evolve_physics_one_substep_try( &
         dt_loc      = dt_nd, &
         sigma_in    = sigma_cur, &
         nu_in       = nu_cur, &
         T_in        = Tmid_cur, &
         source_n    = source_n, &
         source_np1  = source_np1, &
         Qirr_prof_in= Qirr_prof, &
         dYdXi_in    = dYdXi_prof, &
         shadow_in   = shadow_prof, &
         sigma_out   = sigma_try, &
         nu_out      = nu_try, &
         T_out       = T_try, &
         H_out       = H_try, &
         rho_out     = rho_try, &
         kappa_out   = kappa_try, &
         tau_out     = tau_try, &
         Qvis_out    = Qvis_try, &
         Qrad_out    = Qrad_try, &
         Qirr_out    = Qirr_try, &
         dYdXi_out   = dYdXi_try, &
         shadow_out  = shadow_try, &
         kiter_out   = kiter_try, &
         miter_out   = miter_try )

    ! 2) Accept / reject
    ok = accept_step_try(nr, kiter_try, sigma_try, nu_try, T_try, Qvis_try, Qrad_try, Qirr_try)
    if (.not. ok) return

    ! 3) Commit trial -> current
    sigma_cur(:)  = sigma_try(:)
    nu_cur(:)     = nu_try(:)
    Tmid_cur(:)   = T_try(:)
    H_cur(:)      = H_try(:)
    rho_cur(:)    = rho_try(:)
    kappa_cur(:)  = kappa_try(:)
    tau_cur(:)    = tau_try(:)
    Qvis_cur(:)   = Qvis_try(:)
    Qrad_cur(:)   = Qrad_try(:)
    Qirr_cur(:)   = Qirr_try(:)
    dYdXi_cur(:)  = dYdXi_try(:)
    shadow_cur(:) = shadow_try(:)

    ok = .true.
  end function substep_try_and_commit

end module evolve_substep_mod

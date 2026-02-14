module disk_thermal_mod
  !! Local thermal balance for a single radial cell in a Be disk.
  !!
  !! This module provides:
  !!   - heating_cooling_cell_be: given (r, Sigma, OmegaK, Tmid, shadow),
  !!     compute H, rho, kappaR, nu, tauR, Qplus_visc, Qplus_irr, Qminus.

  use kind_params,  only : dp, i4b
  use constants,    only : kb, mp, mu, sbc, pi
  use mod_global,   only : alphaSS, L_star, R_star, kappa0, kappa_es, &
                           use_irradiation
  use irradiation_mod, only : L_irr, cos_inc_min
  use radiation_params_mod, only : tau_min, T_floor, T_ceiling
  use opacity_table_mod, only : opacity_tables, findkappa_OP_AES_S03, &
                                findkappa_OP_S03, findkappa_OP_F05_S03, &
                                get_opacity_Planck_rhoT, &
                         Tcrit_OF_tab, &   ! Boundary temperature of OP-Ferguson
                         Tcrit_FS_tab, &   ! Boundary temperature of Ferguson-Semenov
                         Tcrit_OA_tab, &   ! Boundary temperature of OP-AESOPUS
                         Tcrit_AS_tab, &   ! Boundary temperature of AESOPUS-Semenov
                         thigh_tab, &      ! Max temperature (OP)
                         tlow_tab, &       ! Min temperature (Semenov)
                         logRmax_OF_tab, & ! Max R parameter (OP & Ferguson)
                         logRmax_A_tab, &  ! Max R parameter (AESOPUS)
                         logRmin_tab, &    ! Min R parameter (OP, AESOPUS, & Ferguson)
                         rhomax_tab, &     ! Max density (Semenov)
                         rhomin_tab        ! Min density (Semenov)
  implicit none

contains

  subroutine heating_cooling_cell( r_cgs, Sigma_cgs, OmegaK_cgs, shadow, Tmid, &
                                 H_loc, rho_loc, nu_dim_loc, kappa_loc,        &
                                 tau_loc, Qplus_visc, Qplus_irr, Qminus,       &
                                 Qirr_in )
    use,intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic

    implicit none
    real(dp), intent(in)  :: r_cgs, Sigma_cgs, OmegaK_cgs, Tmid
    real(dp), intent(in), optional :: Qirr_in
    logical,  intent(in)  :: shadow
    real(dp), intent(out) :: H_loc, rho_loc, nu_dim_loc, kappa_loc
    real(dp), intent(out) :: tau_loc, Qplus_visc, Qplus_irr, Qminus

    ! Local variables
    real(dp) :: T_loc, cs2, kappa_ff, kappa_tab, kappaP_loc
    real(dp) :: logR_loc, rho_use, T_use, logR
    real(dp) :: F_irr, cos_geom, cos_star, cos_inc
    real(dp) :: Q_thick, Q_thin, w
    real(dp), parameter :: kappa_floor = 1.0e-3_dp
    real(dp), parameter :: Tcrit_ff_es = 1.0e4_dp
    integer(i4b) :: ierror_kap, ierror

    logical :: irr_on   ! effective irradiation switch (robust)

    !-----------------------------
    ! 0. Initialization (for making debugging easyier)
    !-----------------------------
    H_loc      = transfer(-1_int64,0.0_real64)
    rho_loc    = transfer(-1_int64,0.0_real64)
    nu_dim_loc = transfer(-1_int64,0.0_real64)
    kappa_loc  = transfer(-1_int64,0.0_real64)
    tau_loc    = transfer(-1_int64,0.0_real64)
    Qplus_visc = transfer(-1_int64,0.0_real64)
    Qplus_irr  = transfer(-1_int64,0.0_real64)
    Qminus     = transfer(-1_int64,0.0_real64)

    !-----------------------------
    ! 1. Safety on Tmid
    !-----------------------------
    T_loc = max( T_floor, min(Tmid, T_ceiling) )

    if (Sigma_cgs <= 0.0_dp .or. OmegaK_cgs <= 0.0_dp) then
       H_loc      = 0.0_dp
       rho_loc    = 0.0_dp
       nu_dim_loc = 0.0_dp
       kappa_loc  = 0.0_dp
       tau_loc   = 0.0_dp
       Qplus_visc = 0.0_dp
       Qplus_irr  = 0.0_dp
       Qminus     = 0.0_dp
       return
    end if

    !-----------------------------
    ! 2. Vertical structure
    !-----------------------------
    cs2     = kb * T_loc / (mu * mp)
    H_loc   = sqrt(cs2) / OmegaK_cgs
    rho_loc = Sigma_cgs / (2.0_dp * H_loc)

    nu_dim_loc = 2.0_dp / 3.0_dp * alphaSS * cs2 / OmegaK_cgs

    !-----------------------------
    ! 2b. Opacity (table + fallback)
    !-----------------------------
    !!! For debug
    !!kappa_tab = 1.0e-1_dp
    !!ierror_kap = 0
    select case (opacity_tables)
    case ('OP+AES+S03')
       call findkappa_OP_AES_S03(rho_loc, T_loc, kappa_tab, ierror_kap)
    case ('OP+F05+S03')
       call findkappa_OP_F05_S03(rho_loc, T_loc, kappa_tab, ierror_kap)
    !case ('OP+AES')
    !   call findkappa_OP_AES(rho_loc, T_loc, kappa_tab, ierror_kap)
    case ('OP+S03')
       call findkappa_OP_S03(rho_loc, T_loc, kappa_tab, ierror_kap)
    end select

    if (ierror_kap == 0) then
       kappa_loc = kappa_tab
    else
       !write (*, '("r_cgs =", 1pe12.4, ": rho_loc =", 1pe12.4, ", T_loc =", 1pe12.4, &
       !      " -> ierror =", i2, ", kappa =", 1pe12.4)') r_cgs, rho_loc, T_loc, &
       !      ierror_kap, kappa_tab
       !if (T_loc > Tcrit_ff_es) then
       !   ! Free-free + electron scattering opacity (Kramers + constant)
       !   kappa_loc = kappa0 * rho_loc * T_loc**(-3.5_dp) + kappa_es
       !else
       !   ! Clamping kappa to a value when opacity table fails
          kappa_loc = kappa_tab
       !end if
    end if
    kappa_loc = max(kappa_loc, kappa_floor)
    tau_loc = 0.5_dp * kappa_loc * Sigma_cgs

    ! Planck mean for optically thin cooling (D'Alessio: Q_rad thin ∝ κ_P)
    call get_opacity_Planck_rhoT(rho_loc, T_loc, kappaP_loc, ierror)
    kappaP_loc = max(kappaP_loc, kappa_floor)

    !-----------------------------
    ! 3. Heating and cooling rates
    !-----------------------------
    ! Effective irradiation switch:
    ! - If use_irradiation is false, irradiation must be OFF regardless of "shadow".
    ! - If use_irradiation is true, shadow controls local blocking.
    irr_on = use_irradiation .and. (.not. shadow)

    if (present(Qirr_in)) then
       ! Use externally provided irradiation heating (e.g., LOH24 Eq.16)
       Qplus_irr = max(Qirr_in, 0.0_dp)

    else if (.not. irr_on) then
       Qplus_irr = 0.0_dp

    else
    !   ! Fallback to the legacy simple prescription
    !   if (r_cgs <= 0.0_dp .or. R_star <= 0.0_dp .or. L_star <= 0.0_dp) then
    !      Qplus_irr = 0.0_dp
    !   else
    !      cos_geom = H_loc / r_cgs
    !      cos_star = (2.0_dp / pi) * atan( R_star / r_cgs )
    !      cos_inc  = max( max(cos_geom, cos_star), cos_inc_min )
    !
    !      F_irr     = L_irr / (4.0_dp * pi * r_cgs**2) * cos_inc
    !      Qplus_irr = 2.0_dp * (1.0_dp - albedo) * F_irr
    !   end if
       write (*, '("heating_cooling_cell: Qirr_in not given at r =", 1pe12.4)') r_cgs
       stop

    end if

    ! Viscous heating rate and radiative cooling rate
    call Qvis_and_Qrad_KFM2008(r_cgs, Sigma_cgs, OmegaK_cgs, nu_dim_loc, T_loc, &
                       kappa_loc, kappaP_loc, tau_loc, Qplus_visc, Qminus)

    !-----------------------------
    ! 5. Sanity check
    !-----------------------------
    ierror = 0
    if (.not. ieee_is_finite(H_loc)) ierror = ierror + 1
    if (.not. ieee_is_finite(rho_loc)) ierror = ierror + 1
    if (.not. ieee_is_finite(nu_dim_loc)) ierror = ierror + 1
    if (.not. ieee_is_finite(kappa_loc)) ierror = ierror + 1
    if (.not. ieee_is_finite(tau_loc)) ierror = ierror + 1
    if (.not. ieee_is_finite(Qplus_visc)) ierror = ierror + 1
    if (.not. ieee_is_finite(Qplus_irr)) ierror = ierror + 1
    if (.not. ieee_is_finite(Qminus)) ierror = ierror + 1
    if (ierror /= 0) then
       write (*, '("+++", i2, " variables remain to be substituted their values")') ierror
       write (*, '(4x,"H_loc      =", 1pe12.4)') H_loc
       write (*, '(4x,"rho_loc    =", 1pe12.4)') rho_loc
       write (*, '(4x,"nu_dim_loc =", 1pe12.4)') nu_dim_loc
       write (*, '(4x,"kappa_loc  =", 1pe12.4)') kappa_loc
       write (*, '(4x,"tau_loc    =", 1pe12.4)') tau_loc
       write (*, '(4x,"Qplus_visc =", 1pe12.4)') Qplus_visc
       write (*, '(4x,"Qplus_irr  =", 1pe12.4)') Qplus_irr
       write (*, '(4x,"Qminus     =", 1pe12.4)') Qminus
       stop
    end if

  end subroutine heating_cooling_cell

  subroutine build_shadow_flags(nr, is_shadow)
    use kind_params, only: i4b
    use mod_global,  only: use_irradiation
    implicit none
    integer(i4b), intent(in) :: nr
    logical,      intent(out):: is_shadow(nr)

    if (.not. use_irradiation) then
      is_shadow(:) = .true.    ! irradiation OFF everywhere
    else
      is_shadow(:) = .false.   ! placeholder: no shadow anywhere
    end if
  end subroutine build_shadow_flags

subroutine Qvis_and_Qrad_KFM2008(r_cgs, Sigma_cgs, OmegaK_cgs, nu_dim_loc, T_loc, &
                                kappaR_loc, kappaP_loc, tauR_loc, Qvis_loc, Qrad_loc)
  !-----------------------------------------------------------------------
  ! Compute viscous heating Qvis and radiative cooling Qrad (per unit area)
  ! using coefficients consistent with Kato, Fukue & Mineshige (2008).
  !
  ! Conventions in THIS routine:
  ! - Qvis_loc : total viscous dissipation rate per unit disk surface area,
  !              consistent with KFM2008 Eq.(3.34) as you cited (9/4 factor).
  ! - Qrad_loc : total radiative cooling from BOTH disk faces,
  !              with optically-thick diffusion limit coefficient consistent
  !              with KFM2008 Eq.(3.38) as you cited (64/3 factor).
  ! - tauR_loc is assumed to be tau = kappaR_loc * Sigma_cgs / 2 (one-side),
  !   and is used ONLY for bridging/switching (not for setting coefficients).
  ! - Optically THICK: uses kappaR (Rosseland mean) for diffusion.
  ! - Optically THIN:  uses kappaP (Planck mean) for emission (D'Alessio et al.).
  !-----------------------------------------------------------------------
  use kind_params, only : dp
  use constants,  only : sbc
  implicit none

  real(dp), intent(in)  :: r_cgs, Sigma_cgs, OmegaK_cgs, nu_dim_loc, T_loc
  real(dp), intent(in)  :: kappaR_loc, kappaP_loc, tauR_loc
  real(dp), intent(out) :: Qvis_loc, Qrad_loc

  logical,  parameter   :: use_bridging = .true.
  real(dp), parameter   :: tau_thick    = 4.0_dp / 3.0_dp
  real(dp), parameter   :: tiny_sig     = 1.0e-99_dp
  real(dp), parameter   :: tiny_kap     = 1.0e-99_dp
  real(dp), parameter   :: tiny_q       = 1.0e-99_dp

  real(dp) :: sig_eff, kap_eff, kapP_eff, tau_eff
  real(dp) :: Q_thick, Q_thin, w

  !--- Safety: enforce non-negative / non-singular inputs locally
  sig_eff = max(Sigma_cgs, tiny_sig)
  kap_eff = max(kappaR_loc, tiny_kap)
  kapP_eff = max(kappaP_loc, tiny_kap)

  ! If tauR_loc is passed inconsistently, you can override with this:
  ! tau_eff = 0.5_dp * kap_eff * sig_eff
  tau_eff = max(tauR_loc, 0.0_dp)

  !--- Viscous heating (KFM2008 Eq. 3.34 as cited in your code)
  Qvis_loc = 9.0_dp / 4.0_dp * nu_dim_loc * sig_eff * OmegaK_cgs**2

  !--- Radiative cooling (KFM2008 Eq. 3.38 as cited in your code)
  ! Optically thick (diffusion-like; uses Rosseland mean)
  Q_thick = 64.0_dp * sbc * T_loc**4 / (3.0_dp * kap_eff * sig_eff)

  ! Optically thin limit: Planck mean for emission (D'Alessio: Q_rad ∝ κ_P Σ σT⁴)
  Q_thin  = 2.0_dp * kapP_eff * sig_eff * sbc * T_loc**4

  !--- Combine thick/thin (bridging or sharp switch)
  if (use_bridging) then
     ! Smooth transition weight based on one-side optical depth tau = kappa*Sigma/2
     w = tau_eff / (1.0_dp + tau_eff)

     ! Harmonic blend with weight (robust when one branch is very small)
     Qrad_loc = 1.0_dp / ( (1.0_dp - w) / max(Q_thin,  tiny_q) &
                         + (      w) / max(Q_thick, tiny_q) )
  else
     if (tau_eff >= tau_thick) then
        Qrad_loc = Q_thick
     else
        Qrad_loc = Q_thin
     end if
  end if

end subroutine Qvis_and_Qrad_KFM2008

  subroutine Qvis_and_Qrad(r_cgs, Sigma_cgs, OmegaK_cgs, nu_dim_loc, T_loc, &
                           kappaR_loc, tauR_loc, Qvis_loc, Qrad_loc, kappaP_loc)
    real(dp), intent(in)  :: r_cgs, Sigma_cgs, OmegaK_cgs, nu_dim_loc, T_loc, &
                             kappaR_loc, tauR_loc
    real(dp), intent(out) :: Qvis_loc, Qrad_loc
    real(dp), intent(in), optional :: kappaP_loc
    logical, parameter    :: use_bridging = .true.
    real(dp), parameter   :: tau_thick = 4.0_dp / 3.0_dp
    real(dp)              :: Q_thick, Q_thin, w, kapP_use

    ! Viscous heating rate
    !Qvis_loc = 9.0_dp / 8.0_dp * nu_dim_loc * Sigma_cgs * OmegaK_cgs**2
    ! Eq(3.34) of Kato+ (2008)
    Qvis_loc = 9.0_dp / 4.0_dp * nu_dim_loc * Sigma_cgs * OmegaK_cgs**2

    ! Radiative cooling rate
    ! Optically thick: Rosseland. Optically thin: Planck (when provided)
    kapP_use = kappaR_loc
    if (present(kappaP_loc)) kapP_use = kappaP_loc
    Q_thick = 64.0_dp * sbc * T_loc**4 / (3.0_dp * kappaR_loc * Sigma_cgs)
    Q_thin  = 2.0_dp * kapP_use * Sigma_cgs * sbc * T_loc**4

    if (use_bridging) then
       ! If you want a strict floor to avoid numerical singularities, enable this:
       w        = tauR_loc / (1.0_dp + tauR_loc)
       Qrad_loc = 1.0_dp / ( (1.0_dp - w) / Q_thin + w / Q_thick )
    else
       if (tauR_loc >= tau_thick) then
          Qrad_loc = Q_thick
       else
          Qrad_loc = Q_thin
       end if
    end if

  end subroutine Qvis_and_Qrad

end module disk_thermal_mod

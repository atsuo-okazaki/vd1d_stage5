!===============================================================
! temp_solver_mod.f90
!
! Temperature BE update per cell with damped Newton and bisection fallback.
!
! This version is integrated to a "mod_global" style codebase:
!   - nr from mod_global
!   - caller passes Sigma_new(:), T_old(:), Qirr_lag(:), dt
!   - heating/cooling functions are resolved inside this module
!     by "use heatingcooling_mod" (you will map names accordingly)
!===============================================================
module temp_solver_mod
  use kind_params, only : dp, i4b
  use mod_global,  only : nr
  implicit none
  private

  public :: temp_solver_params
  public :: update_temperature_BE_global

  type :: temp_solver_params
     integer(i4b) :: max_newton_iter = 30_i4b
     integer(i4b) :: max_bisect_iter = 80_i4b
     real(dp)     :: reltol = 1.0e-10_dp
     real(dp)     :: abstol = 1.0e-10_dp
     real(dp)     :: Tmin   = 1.0e1_dp
     real(dp)     :: Tmax   = 1.0e8_dp
     real(dp)     :: fd_eps = 1.0e-6_dp
  end type temp_solver_params

contains

  subroutine update_temperature_BE_global(Sigma_new, T_old, Qirr_lag, dt, pars, T_new, istat)
    real(dp), intent(in)  :: Sigma_new(nr)
    real(dp), intent(in)  :: T_old(nr)
    real(dp), intent(in)  :: Qirr_lag(nr)
    real(dp), intent(in)  :: dt
    type(temp_solver_params), intent(in) :: pars
    real(dp), intent(out) :: T_new(nr)
    integer(i4b), intent(out) :: istat(nr)

    integer(i4b) :: i
    real(dp) :: Ti

    do i = 1, nr
       Ti = clamp(T_old(i), pars%Tmin, pars%Tmax)
       call solve_cell(i, Sigma_new(i), T_old(i), Qirr_lag(i), dt, pars, Ti, istat(i))
       T_new(i) = Ti
    end do
  end subroutine update_temperature_BE_global


  subroutine solve_cell(i, Sigma, Told, Qirr, dt, pars, Tsol, istat)
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: Sigma, Told, Qirr, dt
    type(temp_solver_params), intent(in) :: pars
    real(dp), intent(inout) :: Tsol
    integer(i4b), intent(out) :: istat

    integer(i4b) :: it
    real(dp) :: F, dF, step, Ttrial, lambda
    logical :: conv

    istat = 0_i4b
    conv = .false.

    do it = 1, pars%max_newton_iter
       call cell_F_and_dF(i, Sigma, Told, Qirr, dt, pars, Tsol, F, dF)

       if (is_converged(F, Tsol, pars)) then
          conv = .true.
          exit
       end if

       if (.not. isfinite(dF) .or. abs(dF) <= 0.0_dp) exit

       step = -F / dF

       ! Damped update to ensure residual decreases and T stays physical.
       lambda = 1.0_dp
       do
          Ttrial = clamp(Tsol + lambda*step, pars%Tmin, pars%Tmax)
          if (Ttrial == Tsol) exit
          if (improves(i, Sigma, Told, Qirr, dt, Tsol, Ttrial)) then
             Tsol = Ttrial
             exit
          end if
          lambda = 0.5_dp * lambda
          if (lambda < 1.0e-6_dp) exit
       end do

       if (lambda < 1.0e-6_dp) exit
    end do

    if (conv) then
       istat = 0_i4b
       return
    end if

    call bisection_fallback(i, Sigma, Told, Qirr, dt, pars, Tsol, istat)
  end subroutine solve_cell


  subroutine cell_F_and_dF(i, Sigma, Told, Qirr, dt, pars, T, F, dF)
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: Sigma, Told, Qirr, dt
    type(temp_solver_params), intent(in) :: pars
    real(dp), intent(in) :: T
    real(dp), intent(out) :: F, dF

    real(dp) :: h, Tp, Tm, Fp, Fm

    F = cell_F_only(i, Sigma, Told, Qirr, dt, T)

    h  = pars%fd_eps * max(abs(T), 1.0_dp)
    Tp = clamp(T + h, pars%Tmin, pars%Tmax)
    Tm = clamp(T - h, pars%Tmin, pars%Tmax)

    if (Tp == Tm) then
       dF = 0.0_dp
       return
    end if

    Fp = cell_F_only(i, Sigma, Told, Qirr, dt, Tp)
    Fm = cell_F_only(i, Sigma, Told, Qirr, dt, Tm)
    dF = (Fp - Fm) / (Tp - Tm)
  end subroutine cell_F_and_dF


  function cell_F_only(i, Sigma, Told, Qirr, dt, T) result(F)
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: Sigma, Told, Qirr, dt, T
    real(dp) :: F

    real(dp) :: cv, qvis, qrad

    !-----------------------------------------------------------
    ! Map these three calls to your production code.
    ! Replace the bodies of cv_local / qvis_local / qrad_local
    ! at the bottom of this module.
    !-----------------------------------------------------------
    cv   = cv_local(i, T)
    qvis = qvis_local(i, Sigma, T)
    qrad = qrad_local(i, Sigma, T)

    F = Sigma * cv * (T - Told) / dt - qvis - Qirr + qrad
  end function cell_F_only


  subroutine bisection_fallback(i, Sigma, Told, Qirr, dt, pars, Tsol, istat)
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: Sigma, Told, Qirr, dt
    type(temp_solver_params), intent(in) :: pars
    real(dp), intent(inout) :: Tsol
    integer(i4b), intent(out) :: istat

    real(dp) :: TL, TU, FL, FU, TM, FM
    integer(i4b) :: it
    logical :: ok

    call find_bracket(i, Sigma, Told, Qirr, dt, pars, TL, TU, FL, FU, ok)
    if (.not. ok) then
      Tsol  = clamp(Tsol, pars%Tmin, pars%Tmax)
      istat = 2_i4b
      return
    end if

    do it = 1, pars%max_bisect_iter
      TM = 0.5_dp * (TL + TU)
      FM = cell_F_only(i, Sigma, Told, Qirr, dt, TM)

      if (abs(FM) <= max(pars%abstol, pars%reltol*max(abs(TM), 1.0_dp))) then
         Tsol = TM
         istat = 1_i4b
         return
      end if

      if (FL*FM <= 0.0_dp) then
         TU = TM; FU = FM
      else
         TL = TM; FL = FM
      end if

      if (abs(TU - TL) <= pars%reltol * max(abs(TM), 1.0_dp)) then
         Tsol = TM
         istat = 1_i4b
         return
      end if
    end do

    Tsol = 0.5_dp * (TL + TU)
    istat = 3_i4b
  end subroutine bisection_fallback


  subroutine find_bracket(i, Sigma, Told, Qirr, dt, pars, TL, TU, FL, FU, ok)
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: Sigma, Told, Qirr, dt
    type(temp_solver_params), intent(in) :: pars
    real(dp), intent(out) :: TL, TU, FL, FU
    logical, intent(out) :: ok

    integer(i4b) :: k
    real(dp) :: T0, facL, facU

    ok = .false.
    T0 = clamp(Told, pars%Tmin, pars%Tmax)

    TL = clamp(T0*0.8_dp, pars%Tmin, pars%Tmax)
    TU = clamp(T0*1.2_dp, pars%Tmin, pars%Tmax)
    FL = cell_F_only(i, Sigma, Told, Qirr, dt, TL)
    FU = cell_F_only(i, Sigma, Told, Qirr, dt, TU)
    if (FL*FU <= 0.0_dp) then
       ok = .true.
       return
    end if

    facL = 0.5_dp
    facU = 2.0_dp
    do k = 1, 40
      TL = clamp(T0*facL, pars%Tmin, pars%Tmax)
      TU = clamp(T0*facU, pars%Tmin, pars%Tmax)
      FL = cell_F_only(i, Sigma, Told, Qirr, dt, TL)
      FU = cell_F_only(i, Sigma, Told, Qirr, dt, TU)
      if (FL*FU <= 0.0_dp) then
         ok = .true.
         return
      end if
      if (TL <= pars%Tmin .and. TU >= pars%Tmax) exit
      facL = facL*0.7_dp
      facU = facU*1.3_dp
    end do
  end subroutine find_bracket


  pure logical function is_converged(F, T, pars)
    real(dp), intent(in) :: F, T
    type(temp_solver_params), intent(in) :: pars
    is_converged = (abs(F) <= max(pars%abstol, pars%reltol*max(abs(T), 1.0_dp)))
  end function is_converged

  logical function improves(i, Sigma, Told, Qirr, dt, ToldT, TnewT)
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: Sigma, Told, Qirr, dt
    real(dp), intent(in) :: ToldT, TnewT
    real(dp) :: Fold, Fnew

    Fold = cell_F_only(i, Sigma, Told, Qirr, dt, ToldT)
    Fnew = cell_F_only(i, Sigma, Told, Qirr, dt, TnewT)
    improves = (abs(Fnew) < abs(Fold))
  end function improves

  pure real(dp) function clamp(x, xmin, xmax)
    real(dp), intent(in) :: x, xmin, xmax
    clamp = max(xmin, min(xmax, x))
  end function clamp

  pure logical function isfinite(x)
    real(dp), intent(in) :: x
    isfinite = (x == x) .and. (abs(x) < huge(x))
  end function isfinite


  !=============================================================
  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !  Replace these three functions with your production hooks.
  !
  !  You can either:
  !    (A) "use your_heatingcooling_mod" and call your existing
  !        Qvis/Qrad/cV functions here, or
  !    (B) move these into a small adapter module and keep
  !        temp_solver_mod physics-agnostic.
  !=============================================================

  function cv_local(i, T) result(cv)
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: T
    real(dp) :: cv
    ! Placeholder: constant cV.
    cv = 1.0_dp
  end function cv_local

  function qvis_local(i, Sigma, T) result(qvis)
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: Sigma, T
    real(dp) :: qvis
    ! Placeholder: replace with your viscous heating.
    qvis = 0.0_dp
  end function qvis_local

  function qrad_local(i, Sigma, T) result(qrad)
    integer(i4b), intent(in) :: i
    real(dp), intent(in) :: Sigma, T
    real(dp) :: qrad
    ! Placeholder: replace with your radiative cooling.
    qrad = 0.0_dp
  end function qrad_local

end module temp_solver_mod

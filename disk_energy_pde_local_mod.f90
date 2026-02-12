module disk_energy_pde_mod
  !! Time-dependent temperature PDE:
  !!
  !!   Sigma*cV*(T^{n+1}-T^n)/dt = Qvis(T^{n+1}) + Qirr_lag - Qrad(T^{n+1})
  !!
  !! - cV is constant (radiation pressure ignored).
  !! - Qirr_lag is provided externally (already delayed).
  !! - Uses local Newton iteration per radial cell.

  use kind_params, only : dp, i4b
  use constants,   only : kb, mp, mu
  use disk_thermal_mod, only : heating_cooling_cell
  implicit none

  real(dp), parameter :: cV_const = kb / ((5.0_dp/3.0_dp - 1.0_dp)*mu*mp)
  real(dp), parameter :: tol_relT = 1.0e-3_dp
  integer(i4b), parameter :: newton_max = 20

  public :: update_temperature_pde

contains

  subroutine update_temperature_pde(nr, r_cgs, Sigma_cgs, OmegaK_cgs, shadow, &
                                    T_old, Qirr_lag, dt_phys, T_new)
    integer(i4b), intent(in) :: nr
    real(dp), intent(in) :: r_cgs(nr), Sigma_cgs(nr), OmegaK_cgs(nr)
    logical,  intent(in) :: shadow(nr)
    real(dp), intent(in) :: T_old(nr), Qirr_lag(nr), dt_phys
    real(dp), intent(out):: T_new(nr)

    integer(i4b) :: i, k
    real(dp) :: T, Told, F, dFdT
    real(dp) :: H,rho,nu,kappa,tau,Qvis,Qirr,Qrad
    real(dp) :: H2,rho2,nu2,kappa2,tau2,Qvis2,Qirr2,Qrad2
    real(dp) :: dT

    do i = 1, nr
       if (Sigma_cgs(i)<=0.0_dp .or. OmegaK_cgs(i)<=0.0_dp) then
          T_new(i)=0.0_dp
          cycle
       end if

       T = max(T_old(i), 1.0_dp)
       Told = T

       do k = 1, newton_max
          call heating_cooling_cell(r_cgs(i),Sigma_cgs(i),OmegaK_cgs(i), &
                                    shadow(i),T,H,rho,nu,kappa,tau, &
                                    Qvis,Qirr,Qrad, Qirr_in=Qirr_lag(i))

          F = Sigma_cgs(i)*cV_const*(T-Told)/dt_phys - (Qvis+Qirr-Qrad)

          call heating_cooling_cell(r_cgs(i),Sigma_cgs(i),OmegaK_cgs(i), &
                                    shadow(i),T*1.01_dp,H2,rho2,nu2,kappa2,tau2, &
                                    Qvis2,Qirr2,Qrad2, Qirr_in=Qirr_lag(i))

          dFdT = ( Sigma_cgs(i)*cV_const*(1.01_dp*T-Told)/dt_phys &
                   - (Qvis2+Qirr2-Qrad2) - F )/(0.01_dp*T)

          if (abs(dFdT)<=0.0_dp) exit

          dT = -F/dFdT
          T  = max(1.0_dp, T + dT)

          if (abs(T-Told)/max(T,1.0_dp) < tol_relT) exit
          Told = T
       end do

       T_new(i) = T
    end do
  end subroutine update_temperature_pde

end module disk_energy_pde_mod

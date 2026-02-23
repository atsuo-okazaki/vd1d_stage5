module disk_energy_init_mod
  !! Initial temperature construction consistent with local thermal balance.
  !!
  !! This module computes T(r, t0) by solving:
  !!   Qvis(T) + Qirr(t0) = Qrad(T)
  !!
  !! This avoids violent transients at the first PDE timestep.

  use kind_params, only : dp, i4b
  use disk_thermal_mod, only : heating_cooling_cell
  use radiation_params_mod, only : T_floor, T_ceiling
  implicit none

  public :: init_temperature_equilibrium

contains

  subroutine init_temperature_equilibrium(nr, r, Sigma, OmegaK, shadow, Qirr, T)
    integer(i4b), intent(in) :: nr
    real(dp), intent(in) :: r(nr), Sigma(nr), OmegaK(nr), Qirr(nr)
    logical,  intent(in) :: shadow(nr)
    real(dp), intent(out):: T(nr)

    integer(i4b) :: i, k
    real(dp) :: Tloc, H,rho,nu,kappa,kappaP,tau,Qvis,Qirr_loc,Qrad
    real(dp) :: F, dFdT, Ttrial

    do i = 1, nr
       if (Sigma(i) <= 0.0_dp) then
          T(i) = T_floor
          cycle
       end if

       Tloc = max(T_floor, min(T_ceiling, 0.5_dp*(T_floor+T_ceiling)))

       do k = 1, 50
          call heating_cooling_cell(r(i),Sigma(i),OmegaK(i),shadow(i),Tloc, &
                                    H,rho,nu,kappa,kappaP,tau,Qvis,Qirr_loc,Qrad, &
                                    Qirr_in=Qirr(i))
          F = Qvis + Qirr_loc - Qrad
          if (abs(F)/max(Qrad,1.0_dp) < 1.0e-4_dp) exit

          Ttrial = Tloc * 1.01_dp
          call heating_cooling_cell(r(i),Sigma(i),OmegaK(i),shadow(i),Ttrial, &
                                    H,rho,nu,kappa,kappaP,tau,Qvis,Qirr_loc,Qrad, &
                                    Qirr_in=Qirr(i))
          dFdT = (Qvis + Qirr_loc - Qrad - F)/(Ttrial-Tloc)

          if (abs(dFdT) <= 0.0_dp) exit
          Tloc = max(T_floor, min(T_ceiling, Tloc - F/dFdT))
       end do

       T(i) = Tloc
    end do

  end subroutine init_temperature_equilibrium

end module disk_energy_init_mod

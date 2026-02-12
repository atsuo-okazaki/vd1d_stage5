module units_disk_mod
  use kind_params, only : dp
  use constants,  only : gg, cc, pi
  use mod_global, only : M_star, R0, t0, sigma_init, Temp0, kappa_es
  implicit none

contains

  !===========================================================
  ! Units of mass accretion/decretion rate
  !===========================================================
  pure function mdot_edd_unit() result(mdot_u)
    ! Eddington mass accretion rate [g/s] = L_Edd / c^2
    real(dp) :: mdot_u
    mdot_u = 4.0_dp * pi * gg * M_star / (kappa_es * cc)
  end function mdot_edd_unit

  pure function mdot_pde_unit() result(mdot_u)
    ! PDE unit for mass flow [g/s] = sigma_init * R0^2 / t0
    real(dp) :: mdot_u
    mdot_u = sigma_init * R0*R0 / t0
  end function mdot_pde_unit

  pure function mdot_from_edd(mdot_edd) result(mdot_phys)
    real(dp), intent(in) :: mdot_edd
    real(dp)             :: mdot_phys
    mdot_phys = mdot_edd * mdot_edd_unit()
  end function mdot_from_edd

  pure function mdot_from_nd(mdot_nd) result(mdot_phys)
    real(dp), intent(in) :: mdot_nd
    real(dp)             :: mdot_phys
    mdot_phys = mdot_nd * mdot_pde_unit()
  end function mdot_from_nd

  pure function mdot_nd_from_edd(mdot_edd) result(mdot_nd)
    real(dp), intent(in) :: mdot_edd
    real(dp)             :: mdot_nd
    mdot_nd = (mdot_edd * mdot_edd_unit()) / mdot_pde_unit()
  end function mdot_nd_from_edd

  !===========================================================
  ! Radius and time
  !
  !   r  = R0 * r'
  !   t  = t0 * t'
  !===========================================================
  pure function r_dim(r_nd) result(r)
    ! Dimensionless radius r' -> physical radius r [cm]
    real(dp), intent(in) :: r_nd
    real(dp)             :: r
    r = r_nd * R0
  end function r_dim

  pure function t_dim(t_nd) result(t)
    ! Dimensionless time t' -> physical time t [s]
    real(dp), intent(in) :: t_nd
    real(dp)             :: t
    t = t_nd * t0
  end function t_dim

  pure function r_nd(r) result(r_nd_out)
    ! Physical radius r [cm] -> dimensionless radius r'
    real(dp), intent(in) :: r
    real(dp)             :: r_nd_out
    r_nd_out = r / R0
  end function r_nd

  pure function t_nd(t) result(t_nd_out)
    ! Physical time t [s] -> dimensionless time t'
    real(dp), intent(in) :: t
    real(dp)             :: t_nd_out
    t_nd_out = t / t0
  end function t_nd

  !===========================================================
  ! Keplerian angular frequency
  !
  !   Omega_K(r) = sqrt(G M_star / r^3)
  !
  ! Input is dimensionless radius r', internally converted.
  !===========================================================
  pure function omegaK_dim(r_nd) result(omegaK)
    ! Dimensionless radius r' -> Keplerian angular frequency [s^-1]
    real(dp), intent(in) :: r_nd
    real(dp)             :: r, omegaK
    r      = r_dim(r_nd)
    omegaK = sqrt( gg * M_star / r**3 )
  end function omegaK_dim

  !===========================================================
  ! Surface density
  !
  !   Sigma  = sigma_init * Sigma'
  !===========================================================
  pure function sigma_dim(sigma_nd) result(sigma)
    ! Dimensionless surface density Sigma' -> physical Sigma [g/cm^2]
    real(dp), intent(in) :: sigma_nd
    real(dp)             :: sigma
    sigma = sigma_nd * sigma_init
  end function sigma_dim

  pure function sigma_nd(sigma) result(sigma_nd_out)
    ! Physical Sigma [g/cm^2] -> dimensionless Sigma'
    real(dp), intent(in) :: sigma
    real(dp)             :: sigma_nd_out
    sigma_nd_out = sigma / sigma_init
  end function sigma_nd

  !===========================================================
  ! Temperature (optional)
  !
  ! For now, you mostly keep T in physical [K]. These helpers
  ! are provided in case you later want a dimensionless T' = T/Temp0.
  !===========================================================
  pure function Temp_dim(Temp_nd) result(Temp)
    ! Dimensionless temperature T' -> physical temperature T [K]
    real(dp), intent(in) :: Temp_nd
    real(dp)             :: Temp
    Temp = Temp_nd * Temp0
  end function Temp_dim

  pure function Temp_nd(Temp) result(Temp_nd_out)
    ! Physical temperature T [K] -> dimensionless temperature T'
    real(dp), intent(in) :: Temp
    real(dp)             :: Temp_nd_out
    Temp_nd_out = Temp / Temp0
  end function Temp_nd

end module units_disk_mod

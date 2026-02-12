module mdot_units_mod
  use kind_params, only : dp
  use constants,   only : pi, gg, cc
  use mod_global,  only : M_star, kappa_es, sigma_init, R0, t0
  implicit none
contains
  pure function mdot_edd_unit() result(mdot_u)
    real(dp) :: mdot_u
    mdot_u = 4.0_dp * pi * gg * M_star / (kappa_es * cc)
  end function mdot_edd_unit

  pure function mdot_pde_unit() result(mdot_u)
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
end module mdot_units_mod

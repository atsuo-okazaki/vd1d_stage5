module disk_flux_mod
  use kind_params, only : dp, i4b
  use constants,   only : pi
  use mod_global,  only : r
  use units_disk_mod, only : mdot_edd_unit
  implicit none

contains
  subroutine measure_mdot_inner_from_arrays(sigma_nd, nu_nd, mdot_out)
    ! Recover mdot at the inner boundary from G-gradient using the CURRENT arrays.
    real(dp), intent(in)  :: sigma_nd(:), nu_nd(:)
    real(dp), intent(out) :: mdot_out

    real(dp) :: r_edge1, g_left, g_right, dGdr, r_face
    real(dp) :: mdot_inner_nd, mdot_inner_phys

    mdot_out = 0.0_dp
    if (size(r) < 2) return

    r_edge1 = r(1) - 0.5_dp * ( r(2) - r(1) )
    r_face  = r_edge1

    g_left  = 0.0_dp
    g_right = nu_nd(1) * sigma_nd(1) * sqrt(max(r(1), 0.0_dp))

    dGdr = (g_right - g_left) / max(r(1) - r_edge1, 1.0e-99_dp)

    mdot_inner_nd   = 6.0_dp * pi * sqrt(max(r_face, 0.0_dp)) * dGdr
    mdot_inner_phys = mdot_inner_nd * mdot_edd_unit()

    mdot_out = mdot_inner_phys
  end subroutine measure_mdot_inner_from_arrays

end module disk_flux_mod
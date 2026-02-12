module disk_state_update_mod
  use kind_params, only: dp, i4b
  use mod_global,  only: nr, Tmid, H, rho, kappaR, tauR, nu_conv
  implicit none
contains

  subroutine compute_state_from_sigma(itp1, r_nd, sigma_nd,          &
                                      Tmid_out, H_out, rho_out,      &
                                      kappa_out, tau_out, nu_out,    &
                                      Qvis_out, Qirr_out, Qrad_out)
    ! Pure compute: NO writes to mod_global arrays (except reading past values).
    integer(i4b), intent(in) :: itp1
    real(dp), intent(in)  :: r_nd(:), sigma_nd(:)
    real(dp), intent(out) :: Tmid_out(:), H_out(:), rho_out(:)
    real(dp), intent(out) :: kappa_out(:), tau_out(:), nu_out(:)
    real(dp), intent(out), optional :: Qvis_out(:), Qirr_out(:), Qrad_out(:)
  end subroutine compute_state_from_sigma

  subroutine commit_state(itp1, Tmid_in, H_in, rho_in, kappa_in, tau_in, nu_in, &
                        Qvis_in, Qirr_in, Qrad_in)
    use mod_global, only: Tmid, H, rho, kappaR, tauR, nu_conv, Qvis, Qirr, Qrad
    integer(i4b), intent(in) :: itp1
    real(dp), intent(in) :: Tmid_in(:), H_in(:), rho_in(:), kappa_in(:), tau_in(:), nu_in(:)
    real(dp), intent(in), optional :: Qvis_in(:), Qirr_in(:), Qrad_in(:)

    Tmid(itp1,:)   = Tmid_in(:)
    H(itp1,:)      = H_in(:)
    rho(itp1,:)    = rho_in(:)
    kappaR(itp1,:) = kappa_in(:)
    tauR(itp1,:)   = tau_in(:)
    nu_conv(itp1,:) = nu_in(:)

    if (present(Qvis_in)) Qvis(itp1,:) = Qvis_in(:)
    if (present(Qirr_in)) Qirr(itp1,:) = Qirr_in(:)
    if (present(Qrad_in)) Qrad(itp1,:) = Qrad_in(:)
  end subroutine commit_state

end module disk_state_update_mod

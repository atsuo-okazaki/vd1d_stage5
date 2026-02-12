module state_io_mod
  use kind_params, only : dp, i4b
  use mod_global,  only : nr, &
                          sigmat, Tmid, H, rho, kappaR, tauR, Qvis, Qrad, Qirr, dYdXi, nu_conv, &
                          sigma_cur, nu_cur, Tmid_cur, H_cur, rho_cur, kappa_cur, tau_cur, &
                          Qvis_cur, Qrad_cur, Qirr_cur, dYdXi_cur, shadow_cur
  implicit none
contains

  subroutine load_state_from_history(it)
    use mod_global, only : nr, sigmat, Tmid, nu_conv, k_iter, m_iter, &
                         sigma_cur, Tmid_cur, nu, nu_cur, k_iter_cur, m_iter_cur, &
                         H, rho, kappaR, tauR, Qvis, Qrad, Qirr, dYdXi, is_shadow, &
                         H_cur, rho_cur, kappa_cur, tau_cur, Qvis_cur, Qrad_cur, &
                         Qirr_cur, dYdXi_cur, shadow_cur
    implicit none
    integer(i4b), intent(in) :: it

    sigma_cur(:) = sigmat(it,:)
    Tmid_cur(:)  = Tmid(it,:)

    nu(:)        = nu_conv(it,:)
    nu_cur(:)    = nu_conv(it,:)
    k_iter_cur   = k_iter(it)
    m_iter_cur   = m_iter(it)

    ! Restore current state from history
    H_cur(:)        = H(it,:)
    rho_cur(:)      = rho(it,:)
    kappa_cur(:)    = kappaR(it,:)
    tau_cur(:)      = tauR(it,:)
    Qvis_cur(:)     = Qvis(it,:)
    Qrad_cur(:)     = Qrad(it,:)
    Qirr_cur(:)     = Qirr(it,:)
    dYdXi_cur(:)    = dYdXi(it,:)
    !shadow_cur(:)  = is_shadow(it,:)
    shadow_cur(:)   = .false.

  end subroutine load_state_from_history


  subroutine store_state_to_history(itp1)
    use mod_global, only : nr, sigmat, Tmid, nu_conv, k_iter, m_iter, &
                         sigma_cur, Tmid_cur, nu, nu_cur, k_iter_cur, m_iter_cur, &
                         H, rho, kappaR, tauR, Qvis, Qrad, Qirr, dYdXi, is_shadow, &
                         H_cur, rho_cur, kappa_cur, tau_cur, Qvis_cur, Qrad_cur, &
                         Qirr_cur, dYdXi_cur, shadow_cur
    implicit none
    integer(i4b), intent(in) :: itp1

    sigmat(itp1,:) = sigma_cur(:)
    Tmid(itp1,:)   = Tmid_cur(:)

    nu_conv(itp1,:)= nu_cur(:)
    k_iter(itp1)   = k_iter_cur
    m_iter(itp1)   = m_iter_cur

    ! diagnostics を history に保存するなら
    H(itp1,:)       = H_cur(:)
    rho(itp1,:)     = rho_cur(:)
    kappaR(itp1,:)  = kappa_cur(:)
    tauR(itp1,:)    = tau_cur(:)
    Qvis(itp1,:)    = Qvis_cur(:)
    Qrad(itp1,:)    = Qrad_cur(:)
    Qirr(itp1,:)    = Qirr_cur(:)
    dYdXi(itp1,:)   = dYdXi_cur(:)
    is_shadow(itp1,:)= shadow_cur(:)

  end subroutine store_state_to_history

end module state_io_mod

subroutine evolve_physics_one_step(it, source_n, source_np1, dt_loc)
  use kind_params, only : dp, i4b
  use mod_global,  only : nr, r, nu, sigmat, alpha, &
                          use_energy_balance, use_energy_pde, &
                          use_be_decretion, use_wind_truncation, &
                          nu_conv, k_iter, r_edge, i_edge
  implicit none

  integer(i4b), intent(in) :: it
  real(dp),     intent(in) :: dt_loc          ! <<< dt を引数に
  real(dp),     intent(in) :: source_n(nr), source_np1(nr)

  real(dp) :: sigma_old(nr)
  real(dp) :: sigma_new(nr), sigma_prev(nr)
  real(dp) :: nu_new(nr)
  real(dp) :: sigma_star(nr), sigma_commit(nr)
  real(dp) :: src_theta(nr)

  real(dp) :: err2, sig2, err_rms
  real(dp), parameter :: tiny = 1.0e-12_dp
  integer(i4b) :: i, k
  integer(i4b) :: iedge_local
  real(dp)     :: redge_cgs

  integer(i4b), parameter :: iter_max = 10
  real(dp),     parameter :: eps_iter  = 1.0e-2_dp

  !---------------------------------------------
  ! Sigma^n
  !---------------------------------------------
  sigma_old(:) = sigmat(it, :)

  !---------------------------------------------
  ! Time-centered source (theta method)
  !---------------------------------------------
  src_theta(:) = (1.0_dp - alpha) * source_n(:) + &
                  alpha            * source_np1(:)

  if (.not. use_energy_balance) then
     !==========================================================
     ! Case A: isothermal
     !==========================================================
     call diffusion_theta_step(nr, r, nu, sigma_old, dt_loc, alpha, sigma_new)

     sigma_star(:)   = sigma_new(:)
     sigma_commit(:) = max(0.0_dp, sigma_star(:) + dt_loc * src_theta(:))

     sigmat(it+1, :)  = sigma_commit(:)
     nu_conv(it+1, :) = nu(:)
     k_iter(it+1)     = 0

  else
     !==========================================================
     ! Case B: energy balance
     !==========================================================

     if (use_energy_pde) then
        !-------------------------------------------------------
        ! (PDE) diffuse Sigma with lagged nu
        !-------------------------------------------------------
        call diffusion_theta_step(nr, r, nu, sigma_old, dt_loc, alpha, sigma_new)

        sigma_star(:)   = sigma_new(:)
        sigma_commit(:) = max(0.0_dp, sigma_star(:) + dt_loc * src_theta(:))

        if (use_be_decretion .and. use_wind_truncation) then
           call diagnose_wind_edge(it+1, sigma_commit, iedge_local, redge_cgs)
           i_edge(it+1) = iedge_local
           r_edge(it+1) = redge_cgs
        end if

        sigmat(it+1, :) = sigma_commit(:)

        call energy_pde_step(it, it+1, sigma_commit, compute_local_structure)

        nu_conv(it+1,:) = nu(:)
        k_iter(it+1)    = 0

        if (use_be_decretion) then
           call apply_outer_isothermal_cap(nr, it+1, r, sigma_commit, nu)
        end if

     else
        !=======================================================
        ! Legacy fixed-point iteration
        !=======================================================

        call diffusion_theta_step(nr, r, nu, sigma_old, dt_loc, alpha, sigma_new)
        nu_new(:)    = nu(:)
        k_iter(it+1) = 1

        do k = 1, iter_max
           sigma_prev(:) = sigma_new(:)

           call update_nu_and_temp(nr, r, sigma_prev, nu_new, it+1)

           call diffusion_theta_step(nr, r, nu_new, sigma_old, dt_loc, alpha, sigma_new)

           err2 = 0.0_dp
           sig2 = 0.0_dp
           do i = 1, nr
              err2 = err2 + (sigma_new(i) - sigma_prev(i))**2
              sig2 = sig2 + max(abs(sigma_prev(i)), tiny)**2
           end do

           err_rms = sqrt(err2 / max(sig2, 1.0e-99_dp))

           if (err_rms < eps_iter) then
              k_iter(it+1) = k
              exit
           end if

           if (k == iter_max) k_iter(it+1) = -99
        end do

        sigma_star(:)   = sigma_new(:)
        sigma_commit(:) = max(0.0_dp, sigma_star(:) + dt_loc * src_theta(:))

        if (use_be_decretion .and. use_wind_truncation) then
           call diagnose_wind_edge(it+1, sigma_commit, iedge_local, redge_cgs)
           i_edge(it+1) = iedge_local
           r_edge(it+1) = redge_cgs
        end if

        sigmat(it+1, :) = sigma_commit(:)

        call advance_temperature_and_structure(it+1, dt_loc, sigma_commit, nu_new)

        nu(:)            = nu_new(:)
        nu_conv(it+1, :) = nu_new(:)

        if (use_be_decretion) then
           call apply_outer_isothermal_cap(nr, it+1, r, sigma_commit, nu_new)
        end if

     end if
  end if

end subroutine evolve_physics_one_step

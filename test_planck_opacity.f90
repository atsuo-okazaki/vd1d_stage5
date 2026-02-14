! Quick test for get_opacity_Planck_rhoT
program test_planck_opacity
  use kind_params, only : dp
  use opacity_table_mod, only : init_opacity_tables, get_opacity_Planck_rhoT
  implicit none
  real(dp) :: rho, T, kappaP
  integer :: ierr

  call init_opacity_tables()

  ! T < 10^4 K: Semenov table
  rho = 1.0e-12_dp
  T = 1.0e3_dp
  call get_opacity_Planck_rhoT(rho, T, kappaP, ierr)
  print '(a,1p,e12.4,a,i0)', 'T=1e3 K, rho=1e-12: kappaP=', kappaP, ', ierr=', ierr

  ! T = 10^4 K: boundary (Semenov)
  T = 1.0e4_dp
  call get_opacity_Planck_rhoT(rho, T, kappaP, ierr)
  print '(a,1p,e12.4,a,i0)', 'T=1e4 K, rho=1e-12: kappaP=', kappaP, ', ierr=', ierr

  ! T > 10^4 K: Kramers + kappa_es
  T = 2.0e4_dp
  call get_opacity_Planck_rhoT(rho, T, kappaP, ierr)
  print '(a,1p,e12.4,a,i0)', 'T=2e4 K, rho=1e-12: kappaP=', kappaP, ', ierr=', ierr

  T = 1.0e5_dp
  call get_opacity_Planck_rhoT(rho, T, kappaP, ierr)
  print '(a,1p,e12.4,a,i0)', 'T=1e5 K, rho=1e-12: kappaP=', kappaP, ', ierr=', ierr

  print *, 'done.'
end program test_planck_opacity

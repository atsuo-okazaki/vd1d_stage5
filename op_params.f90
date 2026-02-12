!  *****************************
   module      op_params
!  *****************************
     use kind_params, only : i4b, dp
     implicit none
     save

     !----------------------------------------------------------
     ! Opacity tables:
     !    - OP table + Ferguson+ (2005) + Semenov+ (2003)
     !    - OP table + AESOPUS 2.1 + Semenov+ (2003)
     !    - OP table + Semenov+ (2003)
     !    - OP table + AESOPUS 2.1
     !    (Caution: AESOPUSTable.data is for gas only (no grains)
     !----------------------------------------------------------
     character(len=10), parameter :: opacity_tables = 'OP+AES+S03'
     !character(len=10), parameter :: opacity_tables = 'OP+F05+S03'
     !character(len=10), parameter :: opacity_tables = 'OP+S03'
     !character(len=10), parameter :: opacity_tables = 'OP+AES'
 
     ! Interpolation method applied to opacity tables
     !character(len=6), parameter :: intpol = 'poly'
     character(len=6),  parameter :: intpol  = 'spline'

     character(len=70) :: file_highT, file_midT, file_lowT
     !character(len=70), parameter :: file_highT = 'data/OPTable73.data', &
     !     file_midT  = 'data/FergusonA09Table.data', &
     !     file_lowT  = 'data/SemenovRosselandTable.data'
     integer(i4b) :: NThigh, NTmid, NTlow, NRhigh, NRmid, NRlow, Nrho
     real(dp), parameter :: kes = 0.34_dp, &
          kff0 = 6.2e22_dp ! kappa_ff = kff0 * rho * Tc**(-3.5)
     real(dp), allocatable, dimension(:) :: logThigh, logTmid, logTlow, &
          logRhigh, logRmid, logrho, logRlow
     real(dp), allocatable, dimension(:,:) :: op_highT, op_midT, op_lowT, &
          op2_highT, op2_midT, op2_lowT
     real(dp) :: eta, zeta

  end module op_params

!==============================================================
! opacity_table_mod.f90
!
! Unified interface to OPAL + Semenov opacity tables.
!
! - Uses the existing table reader  : optables
! - Uses the existing interpolation : findkappa(rho, T, kappa, ierror)
!
! Public API for other codes:
!
!   call init_opacity_tables()
!      -> read OPAL + Semenov tables (only once)
!
!   call get_opacity_rhoT(rho, T, kappaR, ierror)
!      -> rho [g/cm^3], T [K]  --> kappaR [cm^2/g]
!         ierror = 0: OK
!                > 0: outside table range, see findkappa
!
!   call get_opacity_Planck_rhoT(rho, T, kappaP, ierror)
!      -> rho [g/cm^3], T [K]  --> kappaP [cm^2/g]  (Planck mean)
!         T <= 10^4 K: Semenov et al. (2003) Planck table
!         T >  10^4 K: Kramers free-free + electron scattering
!
!==============================================================
module opacity_table_mod
  use kind_params,  only : dp, i4b, lgt
  use op_params,    only : opacity_tables, &
                           file_highT, file_midT, file_lowT, &
                           op_highT, op_midT, op_lowT, logThigh, logRhigh, &
                           op2_highT, op2_midT, op2_lowT, &
                           logTmid, logRmid, logTlow, logrho, &
                           NThigh, NTmid, NTlow, NRhigh, NRmid, Nrho, &
                           kes, intpol
  use run_control,  only : iprint
  use spline2D_pac
  use polint2D_pac
  implicit none

  real(dp), parameter :: Tcrit_OF_tab  = 10.0_dp**4.0_dp, &   ! OP-Ferguson
                         Tcrit_FS_tab  = 10.0_dp**3.0_dp, &   ! Ferguson-Semenov
                         Tcrit_OA_tab  = 10.0_dp**4.3_dp, &   ! OP-AESOPUS
                         Tcrit_AS_tab  = 10.0_dp**3.3_dp, &   ! AESOPUS-Semenov
                         thigh_tab   = 1.0e8_dp, &            ! OP
                         tlow_tab    = 1.0e1_dp, &            ! Semenov
                         logRmax_OF_tab = 1.0_dp, &           ! OP & Ferguson
                         logRmax_A_tab = 6.0_dp, &            ! AESOPUS
                         logRmin_tab = -8.0_dp, &             ! OP & Ferguson
                         rhomax_tab  = 1.0e-7_dp, &           ! Semenov
                         rhomin_tab  = 1.0e-17_dp             ! Semenov
  ! Planck mean opacity (Semenov table T <= 10^4 K; Kramers+kappa_es for T > 10^4 K)
  real(dp), allocatable :: op_planck_lowT(:,:)
  real(dp), allocatable :: op2_planck_lowT(:,:)
  logical :: planck_table_read = .false.
  real(dp), parameter :: T_Planck_Semenov_max = 1.0e4_dp
  ! Kramers: kappa_ff = kff_coeff * g_ff * (1-Z)(1+X) * rho * T^(-3.5)  [cm^2/g]
  real(dp), parameter :: kff_coeff = 3.68e22_dp
  real(dp), parameter :: g_ff = 1.0_dp
  real(dp), parameter :: X_abund_P = 0.7_dp, Z_abund_P = 0.02_dp  ! cosmic
  ! Flag to avoid reading the tables more than once
  logical :: opacity_initialized = .false.

contains

  !------------------------------------------------------------
  ! init_opacity_tables
  !
  ! Read OPAL + Semenov opacity tables using the existing
  ! reader subroutine "optables".
  !
  ! This routine is safe to call multiple times; the tables are
  ! actually read only once.
  !------------------------------------------------------------
  subroutine init_opacity_tables()
    implicit none

    character(len=70), parameter :: dir = '/Users/atsuo/Documents/programs/AD/AD1D/vd1d_stage5/data/'
    !external :: optables_OP_AES_S03, optables_OP_F05_S03, optables_OP_S03

    if (opacity_initialized) return

    select case (opacity_tables)
    case ('OP+AES+S03')
       file_highT = trim(dir)//'OPTable73.data'
       file_midT  = trim(dir)//'AESOPUSTable.data'
       file_lowT  = trim(dir)//'SemenovRosselandTable.data'
       call optables_OP_AES_S03()
    case ('OP+F05+S03')
       file_highT = trim(dir)//'OPTable73.data'
       file_midT  = trim(dir)//'FergusonA09Table.data'
       file_lowT  = trim(dir)//'SemenovRosselandTable.data'
       call optables_OP_F05_S03()
    !case ('OP+AES')
    !   file_highT = trim(dir)//'OPTable73.data'
    !   file_lowT  = trim(dir)//'AESOPUSTable.data'
    !   call optables_OP_AES()
    case ('OP+S03')
       file_highT = trim(dir)//'OPTable73.data'
       file_lowT  = trim(dir)//'SemenovRosselandTable.data'
       call optables_OP_S03()
    case default
       write (*, '("+++ opacity_tables = ", a, " is not availble. +++")') &
             opacity_tables
       stop
    end select

    call read_semenov_planck_table()
    opacity_initialized = .true.
  end subroutine init_opacity_tables


  !------------------------------------------------------------
  ! read_semenov_planck_table
  ! Read Planck mean opacity table (Semenov et al. 2003).
  ! Uses same (logT, logrho) grid as Semenov Rosseland; op_lowT,
  ! logTlow, logrho, NTlow, Nrho must already be set by optables.
  !------------------------------------------------------------
  subroutine read_semenov_planck_table()
    use kind_params, only : i4b, dp
    use op_params,   only : logTlow, logrho, NTlow, Nrho
    implicit none
    character(len=70), parameter :: dir = '/Users/atsuo/Documents/programs/AD/AD1D/vd1d_stage5/data/'
    character(len=100), parameter :: file_planck = trim(dir)//'SemenovPlanckTable.data'
    character(len=3000) line
    character(len=5)  :: a5
    character(len=7)  :: a7
    character(len=10) :: a10
    integer(i4b) :: iu, ios, i, j, alloc_err
    integer(i4b) :: nt_read, nrho_read

    if (planck_table_read) return
    allocate(op_planck_lowT(NTlow, Nrho), stat=alloc_err)
    if (alloc_err /= 0) then
       write (*,'("### Error allocating op_planck_lowT ###")')
       stop
    end if

    open(newunit=iu, file=file_planck, iostat=ios, status='old', &
         form='formatted', action='read')
    if (ios /= 0) then
       write (*,'("### Error opening ", a, " ###")') trim(file_planck)
       stop
    end if

    do
       read(iu,'(a)') line
       if (index(line,'# NT:') > 0) exit
    end do
    read(line,"(a5, i4)") a5, nt_read
    read(iu,'(a)') line
    read(line,"(a7, i4)") a7, nrho_read
    if (nt_read /= NTlow .or. nrho_read /= Nrho) then
       write (*,'("### Semenov Planck table dims mismatch: ", 2i5, " vs ", 2i5, " ###")') &
            nt_read, nrho_read, NTlow, Nrho
       stop
    end if

    do
       read(iu,"(a)") line
       if (index(line,'# log T(K)') > 0) exit
    end do
    read(line,"(a10, 301(f9.4,:))") a10, (logrho(j), j = 1, Nrho)
    do i = 1, NTlow
       read(iu,"(a)") line
       read(line, "(f9.4, 1x, 301(f9.4,:))") logTlow(i), &
            (op_planck_lowT(i, j), j = 1, Nrho)
    end do
    close(iu)
    planck_table_read = .true.
  end subroutine read_semenov_planck_table


  !------------------------------------------------------------
  ! get_opacity_Planck_rhoT
  ! Returns Planck mean opacity [cm^2/g] for given rho, T.
  ! T <= 10^4 K: Semenov et al. (2003) table interpolation
  ! T >  10^4 K: Kramers free-free + electron scattering
  !------------------------------------------------------------
  subroutine get_opacity_Planck_rhoT(rho, T, kappaP, ierror)
    use kind_params, only : i4b, dp, lgt
    use op_params,   only : logTlow, logrho, NTlow, Nrho, kes, intpol
    use spline2D_pac
    use polint2D_pac
    implicit none
    real(dp), intent(in)  :: rho, T
    real(dp), intent(out) :: kappaP
    integer(i4b), intent(out) :: ierror

    real(dp), parameter :: eps_edge = 1.0e-10_dp
    real(dp), parameter :: tiny_rho = 1.0e-99_dp
    real(dp), parameter :: logT_P_max = log10(T_Planck_Semenov_max)

    logical(lgt), save :: first_planck = .true.
    integer(i4b) :: alloc_err, jx1, jx2
    real(dp) :: logTx, logrhox, rhox_eff, logkapx, dkapx
    real(dp) :: kappa_ff, kappa_es

    ierror = 0

    if (.not. opacity_initialized) call init_opacity_tables()

    if (T <= T_Planck_Semenov_max) then
       ! Semenov Planck table (log10(kappa) stored)
       rhox_eff = max(rho, tiny_rho)
       if (rhox_eff > rhomax_tab .or. rhox_eff < rhomin_tab) then
          rhox_eff = clamp_interior(rhox_eff, rhomin_tab, rhomax_tab, eps_edge)
          ierror = 1
       end if
       logTx = log10(max(T, tlow_tab))
       logTx = clamp_interior(logTx, log10(tlow_tab), logT_P_max, eps_edge)
       logrhox = log10(rhox_eff)

       if (intpol == 'spline') then
          if (first_planck) then
             allocate(op2_planck_lowT(NTlow, Nrho), stat=alloc_err)
             call splie2(logTlow, logrho, op_planck_lowT, op2_planck_lowT)
             first_planck = .false.
          end if
          logkapx = splin2(logTlow, logrho, op_planck_lowT, op2_planck_lowT, logTx, logrhox)
       else
          jx1 = NTlow
          jx2 = Nrho
          call hunt(logTlow, NTlow, logTx, jx1)
          call hunt(logrho, Nrho, logrhox, jx2)
          call polin2mod(logTlow, logrho, op_planck_lowT, jx1, jx2, logTx, logrhox, logkapx, dkapx)
       end if
       kappaP = 10.0_dp**logkapx
    else
       ! T > 1e4 K: true absorption only (free-free; add bound-free if you have it)
       kappa_ff = kff_coeff * g_ff * (1.0_dp - Z_abund_P) * (1.0_dp + X_abund_P) &
                  * rho * T**(-3.5_dp) 
       kappaP = kappa_ff      ! <-- DO NOT add kappa_es here for thin cooling
    end if
  end subroutine get_opacity_Planck_rhoT


  pure real(dp) function clamp_interior(x, xmin, xmax, eps_edge) result(y)
    use kind_params, only : dp
    real(dp), intent(in) :: x, xmin, xmax, eps_edge
    real(dp) :: a, b
    a = xmin + eps_edge*(xmax - xmin)
    b = xmax - eps_edge*(xmax - xmin)
    y = min(max(x, a), b)
  end function clamp_interior

! ***************************************
  SUBROUTINE optables_OP_AES()
! ***************************************
      use kind_params, only : i4b, dp
      use op_params, only : file_highT, file_lowT, &
           op_highT, op_lowT, logThigh, logRhigh, &
           logTlow, logRlow, NThigh, NTlow, NRhigh, NRlow
      implicit none

      character(len=3000) line
      character(len=5)  :: a5
      character(len=10) :: a10
      character(len=18) :: a18
      character(len=23) :: a23
      integer(i4b), save :: iu_op = -1
      integer(i4b) :: ios, alloc_error, i, j

      !---------------------------------------------------
      ! Read the high-temperature opacity table: OP table
      !---------------------------------------------------
      OPEN (newunit=iu_op, FILE=trim(file_highT), iostat=ios, &
           status='old', form='formatted', action='read')
         IF (ios /= 0) THEN
            WRITE (*,"('### Error in opening ', A70)") file_highT
            STOP
         END IF

      DO
         READ(iu_op,'(A)') line
         if (index(line,'log(T) range') > 0) exit
      END DO
      READ(line,'(a18,i3)') a18, NThigh
      READ(iu_op,'(A)') line
      READ(line,'(a18,i3)') a18, NRhigh

      !-- Allocate global arrays for high-temperature opacities
      ALLOCATE (logThigh(NThigh), STAT=alloc_error)
      ALLOCATE (logRhigh(NRhigh), STAT=alloc_error)
      ALLOCATE (op_highT(NThigh, NRhigh), STAT=alloc_error)

      DO
         READ(iu_op,'(A)') line
         if (index(line,'#logT') > 0) exit
      END DO
      
      !-- R = density[g/cm**3]/T6**3, T6=1.e-6*T[degrees]
      !   logR = log10(R)
      READ (line,"(a5, 301(f5.1,:,2x))") a5, &
           (logRhigh(j), j = 1, NRhigh)
      DO i = 1, NThigh
         READ(iu_op,'(A)') line
         READ (line, "(f4.2, 301(f7.3,:))") logThigh(i), &
              (op_highT(i, j), j = 1, NRhigh)
      END DO
      CLOSE (iu_op)
      !-- Test output
      WRITE (*,"(a70)") file_highT
      WRITE (*, "('log R:', 50(f6.2,:))") (logRhigh(j), j = 1, NRhigh)
      DO i = 1, NThigh
         WRITE (*, "(51(f7.3,:))") logThigh(i), &
              (op_highT(i, j), j = 1, NRhigh)
      END DO

      !----------------------------------------------------
      ! Read the low-temperature opacity table: AESOPUS 2.1
      !----------------------------------------------------
      iu_op = -1
      OPEN (newunit=iu_op, FILE=trim(file_lowT), iostat=ios, &
           status='old', form='formatted', action='read')
         IF (ios /= 0) THEN
            WRITE (*,"('### Error in opening ', A70)") file_lowT
            STOP
         END IF

      DO
         READ(iu_op,'(A)') line
         if (index(line,'log10(T) range') > 0) exit
      END DO
      READ (line,"(a23, i3)") a23, NTlow

      DO
         READ(iu_op,'(A)') line
         if (index(line,'log10(R) range') > 0) exit
      END DO
      READ (line,"(a23, i3)") a23, NRlow

      !-- Allocate global arrays for low-temperature opacities
      ALLOCATE (logTlow(NTlow), STAT=alloc_error)
      ALLOCATE (logRlow(NRlow), STAT=alloc_error)
      ALLOCATE (op_lowT(NTlow, NRlow), STAT=alloc_error)

      DO
         READ(iu_op,'(A)') line
         if (index(line,'# log10(T)') > 0) exit
      END DO
      
      !-- R = density[g/cm**3]/T6**3, T6=1.e-6*T[degrees]
      !   logR = log10(R)
      READ (line,"(a10, 301(f9.4,:))") a10, &
           (logRlow(j), j = 1, NRlow)
      DO i = 1, NTlow
         READ(iu_op,'(A)') line
         READ (line, "(f5.3, 1x, 301(f9.4,:))") logTlow(i), &
              (op_lowT(i, j), j = 1, NRlow)
      END DO
      CLOSE (iu_op)
      !-- Test output
      WRITE (*,"(a70)") file_lowT
      WRITE (*, "('log R:', 50(f7.3,:))") (logRlow(j), j = 1, NRlow)
      DO i = 1, NTlow
         WRITE (*, "(51(f7.3,:))") logTlow(i), &
              (op_lowT(i, j), j = 1, NRlow)
      END DO

      !print *, 'size(logThigh) =', size(logThigh)
      !print *, 'size(logRhigh) =', size(logRhigh)
      !print *, 'size(op_highT,1) =', size(op_highT,1)
      !print *, 'size(op_highT,2) =', size(op_highT,2)
      !print *, 'size(logTlow) =', size(logTlow)
      !print *, 'size(logRlow) =', size(logRlow)
      !print *, 'size(op_lowT,1) =', size(op_lowT,1)
      !print *, 'size(op_lowT,2) =', size(op_lowT,2)
      
  END SUBROUTINE optables_OP_AES


! ***************************************
  SUBROUTINE optables_OP_S03()
! ***************************************
      use kind_params, only : i4b, dp
      use op_params, only : file_highT, file_lowT, &
           op_highT, op_lowT, logThigh, logRhigh, &
           logTlow, logrho, NThigh, NTlow, NRhigh, Nrho
      implicit none

      character(len=1) :: a1
      character(len=5) :: a5
      character(len=7) :: a7
      character(len=10) :: a10
      character(len=18) :: a18
      integer(i4b), save :: iu_op = -1
      integer(i4b) :: ios, alloc_error, i, j, nskip1, nskip2

      !-- Read the opacity table for T >= 10^4 K
      OPEN (newunit=iu_op, FILE=trim(file_highT), iostat=ios, &
           status='old', form='formatted', action='read')
         IF (ios /= 0) THEN
            WRITE (*,"('### Error in opening ', A70)") file_highT
            STOP
         END IF

      nskip1 = 0
      DO
         READ (iu_op,"(a18)") a18
         IF (a18 == '#    log(T) range:') THEN
            EXIT
         ELSE
            nskip1 = nskip1 + 1
         END IF
      END DO
      REWIND (iu_op)
      
      DO i = 1, nskip1
         READ (iu_op,"(a1)") a1
      END DO
      READ (iu_op,"(a18, i3)") a18, NThigh
      READ (iu_op,"(a18, i3)") a18, NRhigh

      !-------------------------------------------------------
      ! Allocate global arrays for high-temperature opacities
      ALLOCATE (logThigh(NThigh), STAT=alloc_error)
      ALLOCATE (logRhigh(NRhigh), STAT=alloc_error)
      ALLOCATE (op_highT(NThigh, NRhigh), STAT=alloc_error)
      !------------------------------------------------------

      nskip2 = nskip1 + 2
      DO
         READ (iu_op,"(a5)") a5
         IF (a5 == '#logT') THEN
            EXIT
         ELSE
            nskip2 = nskip2 + 1
         END IF
      END DO
      REWIND (iu_op)
      
      DO i = 1, nskip2
         READ (iu_op,"(a1)") a1
      END DO
      
      !-- R = density[g/cm**3]/T6**3, T6=1.e-6*T[degrees]
      !   logR = log R
      READ (iu_op,"(a5, 301(f5.1,:,2x))") a5, &
           (logRhigh(j), j = 1, NRhigh)
      DO i = 1, NThigh
         READ (iu_op, "(f4.2, 301(f7.3,:))") logThigh(i), &
              (op_highT(i, j), j = 1, NRhigh)
      END DO
      CLOSE (iu_op)
      !-- Test output
      !WRITE (*,"(a70)") file_highT
      !WRITE (*, "('log R:', 50(f6.2,:))") (logRhigh(j), j = 1, NRhigh)
      !DO i = 1, NThigh
      !   WRITE (*, "(51(f7.3,:))") logThigh(i), &
      !        (op_highT(i, j), j = 1, NRhigh)
      !END DO
      
      !-- Read the opacity le for T < 10^4 K
      iu_op = -1
      OPEN (newunit=iu_op, FILE=trim(file_lowT), iostat=ios, &
           status='old', form='formatted', action='read')
         IF (ios /= 0) THEN
            WRITE (*,"('### Error in opening ', A70)") file_lowT
            STOP
         END IF

      nskip1 = 0
      DO
         READ (iu_op,"(a5)") a5
         IF (a5 == '# NT:') THEN
            EXIT
         ELSE
            nskip1 = nskip1 + 1
         END IF
      END DO
      REWIND (iu_op)
      
      DO i = 1, nskip1
         READ (iu_op,"(a1)") a1
      END DO
      READ (iu_op,"(a5, i4)") a5, NTlow
      READ (iu_op,"(a7, i4)") a7, Nrho

      !-------------------------------------------------------
      ! Allocate global arrays for low-temperature opacities
      ALLOCATE (logTlow(NTlow), STAT=alloc_error)
      ALLOCATE (logrho(Nrho), STAT=alloc_error)
      !ALLOCATE (logRlow(NRlow), STAT=alloc_error)
      ALLOCATE (op_lowT(NTlow, Nrho), STAT=alloc_error)
      !------------------------------------------------------

      nskip2 = nskip1 + 2
      DO
         READ (iu_op,"(a10)") a10
         IF (a10 == '# log T(K)') THEN
            EXIT
         ELSE
            nskip2 = nskip2 + 1
         END IF
      END DO
      REWIND (iu_op)
      
      DO i = 1, nskip2
         READ (iu_op,"(a1)") a1
      END DO
      READ (iu_op,"(a10, 301(f9.4,:))") a10, (logrho(j), j = 1, Nrho)
      DO i = 1, NTlow
         READ (iu_op, "(f9.4, 1x, 301(f9.4,:))") logTlow(i), &
              (op_lowT(i, j), j = 1, Nrho)
      END DO
      CLOSE (iu_op)
      !-- Test output
      !WRITE (*,*)
      !WRITE (*,"(a70)") file_lowT
      !WRITE (*, "('log rho:', 301(f9.4,:))") (logrho(j), j = 1, Nrho)
      !DO i = 1, NTlow
      !   WRITE (*, "(f9.4, 1x, 301(f9.4,:))") logTlow(i), &
      !        (op_lowT(i, j), j = 1, Nrho)
      !END DO

      !print *, 'size(logThigh) =', size(logThigh)
      !print *, 'size(logRhigh) =', size(logRhigh)
      !print *, 'size(op_highT,1) =', size(op_highT,1)
      !print *, 'size(op_highT,2) =', size(op_highT,2)
      !print *, 'size(logTlow) =', size(logTlow)
      !print *, 'size(logrho) =', size(logrho)
      !print *, 'size(op_lowT,1) =', size(op_lowT,1)
      !print *, 'size(op_lowT,2) =', size(op_lowT,2)
      
  END SUBROUTINE optables_OP_S03


! ***************************************
  SUBROUTINE optables_OP_AES_S03()
! ***************************************
      use kind_params, only : i4b, dp
      use op_params, only : file_highT, file_midT, file_lowT, &
           op_highT, op_midT, op_lowT, logThigh, logRhigh, &
           logTmid, logRmid, logTlow, logrho, &
           NThigh, NTmid, NTlow, NRhigh, NRmid, Nrho
      implicit none

      character(len=3000) line
      character(len=5) :: a5
      character(len=7) :: a7
      character(len=10) :: a10
      character(len=18) :: a18
      character(len=23) :: a23
      integer(i4b), save :: iu_op = -1
      integer(i4b) :: ios, alloc_error, i, j

      !---------------------------------------------------
      ! Read the high-temperature opacity table: OP table
      !---------------------------------------------------
      OPEN (newunit=iu_op, FILE=trim(file_highT), iostat=ios, &
           status='old', form='formatted', action='read')
         IF (ios /= 0) THEN
            WRITE (*,"('### Error in opening ', A70)") file_highT
            STOP
         END IF

      DO
         READ(iu_op,'(A)') line
         if (index(line,'log(T) range') > 0) exit
      END DO
      READ(line,'(a18,i3)') a18, NThigh
      READ(iu_op,'(A)') line
      READ(line,'(a18,i3)') a18, NRhigh

      !-- Allocate global arrays for high-temperature opacities
      ALLOCATE (logThigh(NThigh), STAT=alloc_error)
      ALLOCATE (logRhigh(NRhigh), STAT=alloc_error)
      ALLOCATE (op_highT(NThigh, NRhigh), STAT=alloc_error)

      DO
         READ(iu_op,'(A)') line
         if (index(line,'#logT') > 0) exit
      END DO
      
      !-- R = density[g/cm**3]/T6**3, T6=1.e-6*T[degrees]
      !   logR = log10(R)
      READ (line,"(a5, 301(f5.1,:,2x))") a5, &
           (logRhigh(j), j = 1, NRhigh)
      DO i = 1, NThigh
         READ(iu_op,'(A)') line
         READ (line, "(f4.2, 301(f7.3,:))") logThigh(i), &
              (op_highT(i, j), j = 1, NRhigh)
      END DO
      CLOSE (iu_op)
      !-- Test output
      !WRITE (*,"(a70)") file_highT
      !WRITE (*, "('log R:', 50(f6.2,:))") (logRhigh(j), j = 1, NRhigh)
      !DO i = 1, NThigh
      !   WRITE (*, "(51(f7.3,:))") logThigh(i), &
      !        (op_highT(i, j), j = 1, NRhigh)
      !END DO

      !--------------------------------------------------------
      ! Read the mid-temperature opacity table: AESOPUS 2.1
      !--------------------------------------------------------
      iu_op = -1
      OPEN (newunit=iu_op, FILE=trim(file_midT), iostat=ios, &
           status='old', form='formatted', action='read')
         IF (ios /= 0) THEN
            WRITE (*,"('### Error in opening ', A70)") file_midT
            STOP
         END IF

      DO
         READ(iu_op,'(A)') line
         if (index(line,'log10(T) range') > 0) exit
      END DO
      READ (line,"(a23, i3)") a23, NTmid

      DO
         READ(iu_op,'(A)') line
         if (index(line,'log10(R) range') > 0) exit
      END DO
      READ (line,"(a23, i3)") a23, NRmid

      !-- Allocate global arrays for mid-temperature opacities
      ALLOCATE (logTmid(NTmid), STAT=alloc_error)
      ALLOCATE (logRmid(NRmid), STAT=alloc_error)
      ALLOCATE (op_midT(NTmid, NRmid), STAT=alloc_error)

      DO
         READ(iu_op,'(A)') line
         if (index(line,'# log10(T)') > 0) exit
      END DO
      
      !-- R = density[g/cm**3]/T6**3, T6=1.e-6*T[degrees]
      !   logR = log10(R)
      READ (line,"(a10, 301(f9.4,:))") a10, &
           (logRmid(j), j = 1, NRmid)
      DO i = 1, NTmid
         READ(iu_op,'(A)') line
         READ (line, "(f5.3, 1x, 301(f9.4,:))") logTmid(i), &
              (op_midT(i, j), j = 1, NRmid)
      END DO
      CLOSE (iu_op)
      !-- Test output
      !WRITE (*,"(a70)") file_midT
      !WRITE (*, "('log R:', 50(f7.3,:))") (logRmid(j), j = 1, NRmid)
      !DO i = 1, NTmid
      !   WRITE (*, "(51(f7.3,:))") logTmid(i), &
      !        (op_midT(i, j), j = 1, NRmid)
      !END DO
      
      !----------------------------------------------------
      ! Read the low-temperature opacity table: Semenov 03
      !----------------------------------------------------
      iu_op = -1
      OPEN (newunit=iu_op, FILE=trim(file_lowT), iostat=ios, &
           status='old', form='formatted', action='read')
         IF (ios /= 0) THEN
            WRITE (*,"('### Error in opening ', A70)") file_lowT
            STOP
         END IF

      DO
         READ(iu_op,'(A)') line
         if (index(line,'# NT:') > 0) exit
      END DO
      READ (line,"(a5, i4)") a5, NTlow
      READ(iu_op,'(A)') line
      READ (line,"(a7, i4)") a7, Nrho

      !-- Allocate global arrays for low-temperature opacities
      ALLOCATE (logTlow(NTlow), STAT=alloc_error)
      ALLOCATE (logrho(Nrho), STAT=alloc_error)
      !ALLOCATE (logRlow(NRlow), STAT=alloc_error)
      ALLOCATE (op_lowT(NTlow, Nrho), STAT=alloc_error)

      DO
         READ (iu_op,"(a)") line
         if (index(line,'# log T(K)') > 0) exit
      END DO
      READ (line,"(a10, 301(f9.4,:))") a10, (logrho(j), j = 1, Nrho)
      DO i = 1, NTlow
         READ (iu_op,"(a)") line
         READ (line, "(f9.4, 1x, 301(f9.4,:))") logTlow(i), &
              (op_lowT(i, j), j = 1, Nrho)
      END DO
      CLOSE (iu_op)
      !-- Test output
      !WRITE (*,*)
      !WRITE (*,"(a70)") file_lowT
      !WRITE (*, "('log rho:', 301(f9.4,:))") (logrho(j), j = 1, Nrho)
      !DO i = 1, NTlow
      !   WRITE (*, "(f9.4, 1x, 301(f9.4,:))") logTlow(i), &
      !        (op_lowT(i, j), j = 1, Nrho)
      !END DO

      !print *, 'size(logThigh) =', size(logThigh)
      !print *, 'size(logRhigh) =', size(logRhigh)
      !print *, 'size(op_highT,1) =', size(op_highT,1)
      !print *, 'size(op_highT,2) =', size(op_highT,2)
      !print *, 'size(logTmid) =', size(logTmid)
      !print *, 'size(logRmid) =', size(logRmid)
      !print *, 'size(op_midT,1) =', size(op_midT,1)
      !print *, 'size(op_midT,2) =', size(op_midT,2)
      !print *, 'size(logTlow) =', size(logTlow)
      !print *, 'size(logrho) =', size(logrho)
      !print *, 'size(op_lowT,1) =', size(op_lowT,1)
      !print *, 'size(op_lowT,2) =', size(op_lowT,2)
      
  END SUBROUTINE optables_OP_AES_S03


! ***************************************
  SUBROUTINE optables_OP_F05_S03()
! ***************************************
      use kind_params, only : i4b, dp
      !use op_params, only : file_highT, file_midT, file_lowT, &
      !     op_highT, op_midT, op_lowT, logThigh, logRhigh, &
      !     logTmid, logRmid, logTlow, logrho, &
      !     NThigh, NTmid, NTlow, NRhigh, NRmid, Nrho
      implicit none

      character(len=1) :: a1
      character(len=5) :: a5
      character(len=7) :: a7
      character(len=10) :: a10
      character(len=18) :: a18
      integer(i4b), save :: iu_op = -1
      integer(i4b) :: ios, alloc_error, i, j, nskip1, nskip2

      !----------------------------------------------
      ! Read the opacity table for high temperatures
      !----------------------------------------------
      OPEN (newunit=iu_op, FILE=file_highT, iostat=ios, &
           status='old', form='formatted', action='read')
         IF (ios /= 0) THEN
            WRITE (*,"('### Error in opening file_highT: ', A70)") file_highT
            STOP
         END IF

      nskip1 = 0
      DO
         READ (iu_op,"(a18)") a18
         IF (a18 == '#    log(T) range:') THEN
            EXIT
         ELSE
            nskip1 = nskip1 + 1
         END IF
      END DO
      REWIND (iu_op)
      
      DO i = 1, nskip1
         READ (iu_op,"(a1)") a1
      END DO
      READ (iu_op,"(a18, i3)") a18, NThigh
      READ (iu_op,"(a18, i3)") a18, NRhigh

      !-- Allocate global arrays for high-temperature opacities
      ALLOCATE (logThigh(NThigh), STAT=alloc_error)
      ALLOCATE (logRhigh(NRhigh), STAT=alloc_error)
      ALLOCATE (op_highT(NThigh, NRhigh), STAT=alloc_error)

      nskip2 = nskip1 + 2
      DO
         READ (iu_op,"(a5)") a5
         IF (a5 == '#logT') THEN
            EXIT
         ELSE
            nskip2 = nskip2 + 1
         END IF
      END DO
      REWIND (iu_op)
      
      DO i = 1, nskip2
         READ (iu_op,"(a1)") a1
      END DO
      
      !-- R = density[g/cm**3]/T6**3, T6=1.e-6*T[degrees]
      !   logR = log R
      READ (iu_op,"(a5, 301(f5.1,:,2x))") a5, &
           (logRhigh(j), j = 1, NRhigh)
      DO i = 1, NThigh
         READ (iu_op, "(f4.2, 301(f7.3,:))") logThigh(i), &
              (op_highT(i, j), j = 1, NRhigh)
      END DO
      CLOSE (iu_op)
      !-- Test output
      !WRITE (*,"(a70)") file_highT
      !WRITE (*, "('log R:', 50(f6.2,:))") (logRhigh(j), j = 1, NRhigh)
      !DO i = 1, NThigh
      !   WRITE (*, "(51(f7.3,:))") logThigh(i), &
      !        (op_highT(i, j), j = 1, NRhigh)
      !END DO

      !--------------------------------------------------------
      ! Read the opacity table for temperatures around 10**4 K
      !--------------------------------------------------------
      iu_op = -1
      OPEN (newunit=iu_op, FILE=file_midT, iostat=ios, &
           status='old', form='formatted', action='read')
         IF (ios /= 0) THEN
            WRITE (*,"('### Error in opening file_midT: ', A70)") file_midT
            STOP
         END IF

      nskip1 = 0
      DO
         READ (iu_op,"(a5)") a5
         IF (a5 == '# NT:') THEN
            EXIT
         ELSE
            nskip1 = nskip1 + 1
         END IF
      END DO
      REWIND (iu_op)
      
      DO i = 1, nskip1
         READ (iu_op,"(a1)") a1
      END DO
      READ (iu_op,"(a5, i3)") a5, NTmid
      READ (iu_op,"(a5, i3)") a5, NRmid

      !-- Allocate global arrays for mid-temperature opacities
      ALLOCATE (logTmid(NTmid), STAT=alloc_error)
      ALLOCATE (logRmid(NRmid), STAT=alloc_error)
      ALLOCATE (op_midT(NTmid, NRmid), STAT=alloc_error)

      nskip2 = nskip1 + 2
      DO
         READ (iu_op,"(a7)") a7
         IF (a7 == '#log T ') THEN
            EXIT
         ELSE
            nskip2 = nskip2 + 1
         END IF
      END DO
      REWIND (iu_op)
      
      DO i = 1, nskip2
         READ (iu_op,"(a1)") a1
      END DO
      
      !-- R = density[g/cm**3]/T6**3, T6=1.e-6*T[degrees]
      !   logR = log R
      READ (iu_op,"(a7, 301(f7.3,:))") a5, &
           (logRmid(j), j = 1, NRmid)
      DO i = 1, NTmid
         READ (iu_op, "(f5.3, 1x, 301(f7.3,:))") logTmid(i), &
              (op_midT(i, j), j = 1, NRmid)
      END DO
      CLOSE (iu_op)
      !-- Test output
      !WRITE (*,"(a70)") file_midT
      !WRITE (*, "('log R:', 50(f7.3,:))") (logRmid(j), j = 1, NRmid)
      !DO i = 1, NTmid
      !   WRITE (*, "(51(f7.3,:))") logTmid(i), &
      !        (op_midT(i, j), j = 1, NRmid)
      !END DO
      
      !----------------------------------------------
      ! Read the opacity table for low temperatures
      !----------------------------------------------
      iu_op = -1
      OPEN (newunit=iu_op, FILE=file_lowT, iostat=ios, &
           status='old', form='formatted', action='read')
         IF (ios /= 0) THEN
            WRITE (*,"('### Error in opening file_lowT: ', A70)") file_lowT
            STOP
         END IF

      nskip1 = 0
      DO
         READ (iu_op,"(a5)") a5
         IF (a5 == '# NT:') THEN
            EXIT
         ELSE
            nskip1 = nskip1 + 1
         END IF
      END DO
      REWIND (iu_op)
      
      DO i = 1, nskip1
         READ (iu_op,"(a1)") a1
      END DO
      READ (iu_op,"(a5, i4)") a5, NTlow
      READ (iu_op,"(a7, i4)") a7, Nrho

      !-- Allocate global arrays for low-temperature opacities
      ALLOCATE (logTlow(NTlow), STAT=alloc_error)
      ALLOCATE (logrho(Nrho), STAT=alloc_error)
      !ALLOCATE (logRlow(NRlow), STAT=alloc_error)
      ALLOCATE (op_lowT(NTlow, Nrho), STAT=alloc_error)

      nskip2 = nskip1 + 2
      DO
         READ (iu_op,"(a10)") a10
         IF (a10 == '# log T(K)') THEN
            EXIT
         ELSE
            nskip2 = nskip2 + 1
         END IF
      END DO
      REWIND (iu_op)
      
      DO i = 1, nskip2
         READ (iu_op,"(a1)") a1
      END DO
      READ (iu_op,"(a10, 301(f9.4,:))") a10, (logrho(j), j = 1, Nrho)
      DO i = 1, NTlow
         READ (iu_op, "(f9.4, 1x, 301(f9.4,:))") logTlow(i), &
              (op_lowT(i, j), j = 1, Nrho)
      END DO
      CLOSE (iu_op)
      !-- Test output
      !WRITE (*,*)
      !WRITE (*,"(a70)") file_lowT
      !WRITE (*, "('log rho:', 301(f9.4,:))") (logrho(j), j = 1, Nrho)
      !DO i = 1, NTlow
      !   WRITE (*, "(f9.4, 1x, 301(f9.4,:))") logTlow(i), &
      !        (op_lowT(i, j), j = 1, Nrho)
      !END DO

      !print *, 'size(logThigh) =', size(logThigh)
      !print *, 'size(logRhigh) =', size(logRhigh)
      !print *, 'size(op_highT,1) =', size(op_highT,1)
      !print *, 'size(op_highT,2) =', size(op_highT,2)
      !print *, 'size(logTmid) =', size(logTmid)
      !print *, 'size(logRmid) =', size(logRmid)
      !print *, 'size(op_midT,1) =', size(op_midT,1)
      !print *, 'size(op_midT,2) =', size(op_midT,2)
      !print *, 'size(logTlow) =', size(logTlow)
      !print *, 'size(logrho) =', size(logrho)
      !print *, 'size(op_lowT,1) =', size(op_lowT,1)
      !print *, 'size(op_lowT,2) =', size(op_lowT,2)
      
  END SUBROUTINE optables_OP_F05_S03


  !------------------------------------------------------------
  ! findkappa
  !
  ! Original CBD routine, slightly wrapped into this module.
  !
  !  Input:
  !    rhox   : density [g/cm^3]
  !    tempx  : temperature [K]
  !
  !  Output:
  !    kappax : Rosseland-mean opacity [cm^2/g]
  !    ierror : error flag
  !             0 : success
  !             1 : T below Semenov table range
  !             2 : rho outside Semenov table range
  !
  !  Internal logic:
  !    - High-T   (logT >= logTcrit_OA_tab + dlogT_tab): OPAL only
  !    - Low-T    (logT <= logTcrit_OA_tab - dlogT_tab): Semenov only
  !    - Middle-T (transition): smooth blend between the two
  !
  !  Units:
  !    - cgs (rho in g/cm^3, T in K, kappa in cm^2/g)
  !------------------------------------------------------------
  subroutine findkappa_OP_S03(rhox, tempx, kappax, ierror)
    use kind_params, only : i4b, dp, lgt
    ! Tables and parameters are imported from op_params at module level

    real(dp), intent(in)  :: tempx, rhox
    real(dp), intent(out) :: kappax
    integer(i4b), intent(out) :: ierror

    !character(len=3), parameter :: testing = 'on'
    character(len=3),  parameter :: testing = 'off'

    integer(i4b), save      :: n_warn = 0
    integer(i4b), parameter :: max_warn = 100

    integer(i4b) :: alloc_error
    integer(i4b) :: jx1, jx2

    ! Opacity table parameters
    real(dp), parameter :: dlogT_tab = 0.05_dp
    real(dp) :: logRx, logrhox, logTx, logkapx, dkapx
    real(dp) :: logTcrit_OA_tab, kap_high, kap_low, w

    logical(lgt), save :: first_highT = .true., &
         first_lowT  = .true.

    ierror = 0

    !---------------------------------------------------------
    ! Optional sanity check: force Thomson opacity for debugging.
    !---------------------------------------------------------
    if (testing == 'on') then
       kappax = kes
       return
    end if

    logrhox  = log10(rhox)
    logTx    = log10(tempx)
    logRx    = logrhox - 3.0_dp*logTx + 18.0_dp
    logTcrit_OA_tab = log10(Tcrit_OA_tab)

    !---------------------------------------------------------
    ! Temperature regions:
    !
    !   logT >= logTcrit_OA_tab + dlogT_tab : pure OP       (high-T)
    !   logT <= logTcrit_OA_tab - dlogT_tab : pure Semenov  (low-T)
    !   otherwise               : smooth blend across the transition.
    !---------------------------------------------------------

    if (logTx >= logTcrit_OA_tab + dlogT_tab) then
       if (tempx > thigh_tab) then
          ! Temperature is above the valid range of OP tables.
          if (iprint >= 1 .and. n_warn < max_warn) then
             write (*,"('Temperature (', 1pe13.5, &
             &         ') is too high to apply OP table.')") tempx
             n_warn = n_warn + 1
             if (n_warn == max_warn) then
                write (*,*) 'Further OP/Semenov opacity warnings suppressed...'
             end if
          end if
          ierror = 1
          return
       else if (logRx > logRmax_OF_tab .or. logRx < logRmin_tab) then
          ! R (= rhox/(tempx/1.0e6)^3) is outside the valid range of OP tables.
          if (iprint >= 1 .and. n_warn < max_warn) then
             write (*,'("R (=rho/T6^3)",1pe13.5,") is out of range to apply ", &
             &         "OP tables (rho =", 1pe13.5, ", T = ",1pe13.5,")")') &
             &         10.0_dp**logRx, rhox, tempx
             n_warn = n_warn + 1
             if (n_warn == max_warn) then
                write (*,*) 'Further OP/Semenov opacity warnings suppressed...'
             end if
          end if
          ierror = 2
          return
       end if

       !---------------- High-temperature branch: OPAL only ----------------
       if (intpol == 'spline') then
          if (first_highT) then
             ! Initialize spline for high-T OP table (one-time).
             allocate (op2_highT(NThigh, NRhigh), stat=alloc_error)
             call splie2(logThigh, logRhigh, op_highT, op2_highT)
             first_highT = .false.
          end if
          logkapx = splin2(logThigh, logRhigh, op_highT, op2_highT, &
                           logTx, logRx)
       else
          jx1 = NThigh
          jx2 = NRhigh
          call hunt(logThigh, NThigh, logTx, jx1)
          call hunt(logRhigh, NRhigh, logRx, jx2)
          call polin2mod(logThigh, logRhigh, op_highT, &
                         jx1, jx2, logTx, logRx, logkapx, dkapx)
       end if
       kappax = 10.0_dp**logkapx
       !print '(a, 1pe12.4, a, 1pe12.4, a, 1pe12.4)', 'rhox=', rhox, &
       !      ', tempx=', tempx, '-> kappax =', kappax

    else if (logTx <= logTcrit_OA_tab - dlogT_tab) then

       !---------------- Low-temperature branch: Semenov only ---------------
       !    tlow   <= T   <= Tcrit
       !    rhomin <= rho <= rhomax

       if (tempx < tlow_tab) then
          ! Temperature is below the valid range of Semenov's tables.
          if (iprint >= 1 .and. n_warn < max_warn) then
             write (*,"('Temperature (', 1pe13.5, &
             &         ') is too low to apply Semenov''s model.')") tempx
             n_warn = n_warn + 1
             if (n_warn == max_warn) then
                write (*,*) 'Further OP/Semenov opacity warnings suppressed...'
             end if
          end if
          ierror = 1
          return

       else if (rhox > rhomax_tab .or. rhox < rhomin_tab) then
          ! Density is outside the valid range of Semenov's tables.
          if (iprint >= 1 .and. n_warn < max_warn) then
             write (*,'("Density (",1pe13.5,") is out of range to apply ", &
             &         "Semenov''s model (T = ",1pe13.5,")")') rhox, tempx
             n_warn = n_warn + 1
             if (n_warn == max_warn) then
                write (*,*) 'Further OP/Semenov opacity warnings suppressed...'
             end if
          end if
          ierror = 2
          return
       end if

       if (intpol == 'spline') then
          if (first_lowT) then
             ! Initialize spline for low-T Semenov table (one-time).
             allocate (op2_lowT(NTlow, Nrho), stat=alloc_error)
             call splie2(logTlow, logrho, op_lowT, op2_lowT)
             first_lowT = .false.
          end if
          logkapx = splin2(logTlow, logrho, op_lowT, op2_lowT, &
                           logTx, logrhox)
       else
          jx1 = NTlow
          jx2 = Nrho
          call hunt(logTlow, NTlow, logTx, jx1)
          call hunt(logrho, Nrho, logrhox, jx2)
          call polin2mod(logTlow, logrho, op_lowT, &
                         jx1, jx2, logTx, logrhox, logkapx, dkapx)
       end if
       kappax = 10.0_dp**logkapx

    else

       !---------------- Transition region: smooth blend --------------------
       !   logTcrit_OA_tab - dlogT_tab < logT < logTcrit_OA_tab + dlogT_tab
       !   logT = logTcrit_OA_tab - dlogT_tab -> pure Semenov (w = 1)
       !   logT = logTcrit_OA_tab + dlogT_tab -> pure OPAL    (w = 0)

       ! High-T opacity from OPAL.
       if (intpol == 'spline') then
          if (first_highT) then
             allocate (op2_highT(NThigh, NRhigh), stat=alloc_error)
             call splie2(logThigh, logRhigh, op_highT, op2_highT)
             first_highT = .false.
          end if
          logkapx = splin2(logThigh, logRhigh, op_highT, op2_highT, &
                           logTx, logRx)
       else
          jx1 = NThigh
          jx2 = NRhigh
          call hunt(logThigh, NThigh, logTx, jx1)
          call hunt(logRhigh, NRhigh, logRx, jx2)
          call polin2mod(logThigh, logRhigh, op_highT, &
                         jx1, jx2, logTx, logRx, logkapx, dkapx)
       end if
       kap_high = 10.0_dp**logkapx

       ! Low-T opacity from Semenov (with the same validity checks).
       if (tempx < tlow_tab) then
          if (iprint >= 1 .and. n_warn < max_warn) then
             write (*,"('Temperature (', 1pe13.5, &
               &       ') is too low to apply Semenov''s model.')") tempx
             n_warn = n_warn + 1
             if (n_warn == max_warn) then
                write (*,*) 'Further OP/Semenov opacity warnings suppressed...'
             end if
          end if
          ierror = 1
          return

       else if (rhox > rhomax_tab .or. rhox < rhomin_tab) then
          if (iprint >= 1 .and. n_warn < max_warn) then
             write (*,'("Density (",1pe13.5,") is out of range to apply ", &
               &       "Semenov''s model (T = ",1pe13.5,")")') rhox, tempx
             n_warn = n_warn + 1
             if (n_warn == max_warn) then
                write (*,*) 'Further OP/Semenov opacity warnings suppressed...'
             end if
          end if
          ierror = 2
          return
       end if

       if (intpol == 'spline') then
          if (first_lowT) then
             allocate (op2_lowT(NTlow, Nrho), stat=alloc_error)
             call splie2(logTlow, logrho, op_lowT, op2_lowT)
             first_lowT = .false.
          end if
          logkapx = splin2(logTlow, logrho, op_lowT, op2_lowT, &
                           logTx, logrhox)
       else
          jx1 = NTlow
          jx2 = Nrho
          call hunt(logTlow, NTlow, logTx, jx1)
          call hunt(logrho, Nrho, logrhox, jx2)
          call polin2mod(logTlow, logrho, op_lowT, &
                         jx1, jx2, logTx, logrhox, logkapx, dkapx)
       end if
       kap_low = 10.0_dp**logkapx

       ! Linear weight in logT across the transition region
       w = (logTcrit_OA_tab + dlogT_tab - logTx)/(2.0_dp*dlogT_tab)
       w = max(0.0_dp, min(1.0_dp, w))

       kappax = w*kap_low + (1.0_dp - w)*kap_high

    end if

  end subroutine findkappa_OP_S03


! ***********************************************************
SUBROUTINE findkappa_OP_AES_S03(rhox, tempx, kappa, ierror)
! ***********************************************************
  use kind_params, only : i4b, dp, lgt
  use op_params, only : op_highT, op_midT, op_lowT, &
       op2_highT, op2_midT, op2_lowT, &
       logThigh, logRhigh, logTmid, logRmid, &
       logTlow, logrho, NThigh, NTmid, NTlow, &
       NRhigh, NRmid, Nrho, kes
  use spline2D_pac
  use polint2D_pac
  implicit none

  !character(len=3),  parameter :: testing = 'on'
  character(len=3),  parameter :: testing = 'off'

  integer(i4b), intent(out) :: ierror
  integer(i4b)              :: alloc_error
  integer(i4b)              :: jx1, jx2
  real(dp), intent(in)      :: tempx, rhox
  real(dp), intent(out)     :: kappa

  ! Opacity table parameters
  real(dp), parameter :: dlogT_OA_tab = 0.05_dp, dlogT_AS_tab = 0.05_dp
  real(dp), parameter :: eps_edge = 1.0e-10_dp
  real(dp), parameter :: tiny_rho = 1.0e-99_dp

  integer(i4b), save      :: n_warn = 0
  integer(i4b), parameter :: max_warn = 100

  real(dp) :: logTx_eff, logrho_phys, logR_phys
  real(dp) :: logTcrit_OA_tab, logTcrit_AS_tab
  real(dp) :: kap_high, kap_mid, kap_low, w
  real(dp) :: logkapx, dkapx

  real(dp) :: rhox_eff, logrho_eff
  real(dp) :: tempx_eff
  integer(i4b) :: ierrT, ierrX

  real(dp) :: logR_OP, logR_A

  logical(lgt), save :: first_highT = .true., &
                        first_midT  = .true., &
                        first_lowT  = .true.

  !-- Optional sanity check: force Thomson opacity for debugging.
  if (testing == 'on') then
     kappa  = kes
     ierror = 0
     return
  end if

  ierror = 0
  ierrT  = 0
  ierrX  = 0

  !---------------------------------------------------------
  ! 0) Precompute transition temperatures in log10 space
  !---------------------------------------------------------
  logTcrit_OA_tab = log10(Tcrit_OA_tab)
  logTcrit_AS_tab = log10(Tcrit_AS_tab)

  !---------------------------------------------------------
  ! 1) Global temperature clamp (only temperature)
  !---------------------------------------------------------
  tempx_eff = tempx
  if (tempx_eff > thigh_tab) then
     tempx_eff = thigh_tab
     ierrT = 1
  else if (tempx_eff < tlow_tab) then
     tempx_eff = tlow_tab
     ierrT = 1
  end if
  tempx_eff = clamp_interior(tempx_eff, tlow_tab, thigh_tab, eps_edge)

  logTx_eff = log10(tempx_eff)

  !---------------------------------------------------------
  ! 2) Compute "physical" logR for OP/AESOPUS from raw rho
  !    (do NOT use Semenov-clamped rho here)
  !---------------------------------------------------------
  logrho_phys = log10(max(rhox, tiny_rho))
  logR_phys   = logrho_phys - 3.0_dp*logTx_eff + 18.0_dp

  !-----------------------------------------------------------------
  ! Temperature regions (ALL comparisons use logTx_eff)
  !
  !   logT >= logTcrit_OA_tab + dlogT_OA_tab : pure OPAL        (high-T)
  !   logT <= logTcrit_AS_tab - dlogT_AS_tab : pure Semenov     (low-T)
  !   logTcrit_AS_tab + dlogT_AS_tab <= logT <= logTcrit_OA_tab - dlogT_OA_tab : pure AESOPUS
  !   otherwise : smooth blends across the two transition regions.
  !-----------------------------------------------------------------

  if (logTx_eff >= logTcrit_OA_tab + dlogT_OA_tab) then
     !---------------- High-temperature branch: OPAL only ----------------
     logR_OP = logR_phys
     if (logR_OP > logRmax_OF_tab .or. logR_OP < logRmin_tab) then
        logR_OP = clamp_interior(logR_OP, logRmin_tab, logRmax_OF_tab, eps_edge)
        ierrX = 1
     end if

     if (intpol == 'spline') then
        if (first_highT) then
           allocate(op2_highT(NThigh, NRhigh), stat=alloc_error)
           call splie2(logThigh, logRhigh, op_highT, op2_highT)
           first_highT = .false.
        end if
        logkapx = splin2(logThigh, logRhigh, op_highT, op2_highT, logTx_eff, logR_OP)
     else
        jx1 = NThigh
        jx2 = NRhigh
        call hunt(logThigh, NThigh, logTx_eff, jx1)
        call hunt(logRhigh, NRhigh, logR_OP, jx2)
        call polin2mod(logThigh, logRhigh, op_highT, jx1, jx2, logTx_eff, logR_OP, logkapx, dkapx)
     end if
     kappa = 10.0_dp**logkapx

     ! For debug
     !write (*, '("tempx =", 1pe12.4, " -> logTx_eff =", 1pe12.4, " -> ", a)') &
     !      tempx, logTx_eff, 'OP'

  else if (logTx_eff <= logTcrit_AS_tab - dlogT_AS_tab) then
     !---------------- Low-temperature branch: Semenov only ---------------
     ! Semenov uses (logT, logrho): clamp rho ONLY here
     rhox_eff = rhox
     if (rhox_eff > rhomax_tab .or. rhox_eff < rhomin_tab) then
        rhox_eff = clamp_interior(rhox_eff, rhomin_tab, rhomax_tab, eps_edge)
        ierrX = 1
     end if
     logrho_eff = log10(max(rhox_eff, tiny_rho))

     if (intpol == 'spline') then
        if (first_lowT) then
           allocate(op2_lowT(NTlow, Nrho), stat=alloc_error)
           call splie2(logTlow, logrho, op_lowT, op2_lowT)
           first_lowT = .false.
        end if
        logkapx = splin2(logTlow, logrho, op_lowT, op2_lowT, logTx_eff, logrho_eff)
     else
        jx1 = NTlow
        jx2 = Nrho
        call hunt(logTlow, NTlow, logTx_eff, jx1)
        call hunt(logrho, Nrho, logrho_eff, jx2)
        call polin2mod(logTlow, logrho, op_lowT, jx1, jx2, logTx_eff, logrho_eff, logkapx, dkapx)
     end if
     kappa = 10.0_dp**logkapx

     ! For debug
     !write (*, '("tempx =", 1pe12.4, " -> logTx_eff =", 1pe12.4, " -> ", a)') &
     !      tempx, logTx_eff, 'Semenov'

  else if (logTx_eff >= logTcrit_AS_tab + dlogT_AS_tab .and. &
           logTx_eff <= logTcrit_OA_tab - dlogT_OA_tab) then
     !---------------- Intermediate branch: AESOPUS only ------------------
     logR_A = logR_phys
     if (logR_A > logRmax_A_tab .or. logR_A < logRmin_tab) then
        logR_A = clamp_interior(logR_A, logRmin_tab, logRmax_A_tab, eps_edge)
        ierrX = 1
     end if

     if (intpol == 'spline') then
        if (first_midT) then
           allocate(op2_midT(NTmid, NRmid), stat=alloc_error)
           call splie2(logTmid, logRmid, op_midT, op2_midT)
           first_midT = .false.
        end if
        logkapx = splin2(logTmid, logRmid, op_midT, op2_midT, logTx_eff, logR_A)
     else
        jx1 = NTmid
        jx2 = NRmid
        call hunt(logTmid, NTmid, logTx_eff, jx1)
        call hunt(logRmid, NRmid, logR_A, jx2)
        call polin2mod(logTmid, logRmid, op_midT, jx1, jx2, logTx_eff, logR_A, logkapx, dkapx)
     end if
     kappa = 10.0_dp**logkapx

     ! For debug
     !write (*, '("tempx =", 1pe12.4, " -> logTx_eff =", 1pe12.4, " -> ", a)') &
     !      tempx, logTx_eff, 'AESOPUS'

  else
     !---------------- Transition regions: smooth blends ------------------
     if (logTx_eff <= logTcrit_AS_tab + dlogT_AS_tab) then
        ! Blend Semenov (low) and AESOPUS (mid) around Tcrit_AS_tab.
        ! logT = logTcrit_AS_tab - dlogT_AS_tab -> pure Semenov (w = 1)
        ! logT = logTcrit_AS_tab + dlogT_AS_tab -> pure AESOPUS (w = 0)

        ! --- Semenov side (low) ---
        rhox_eff = rhox
        if (rhox_eff > rhomax_tab .or. rhox_eff < rhomin_tab) then
           rhox_eff = clamp_interior(rhox_eff, rhomin_tab, rhomax_tab, eps_edge)
           ierrX = 1
        end if
        logrho_eff = log10(max(rhox_eff, tiny_rho))

        if (intpol == 'spline') then
           if (first_lowT) then
              allocate(op2_lowT(NTlow, Nrho), stat=alloc_error)
              call splie2(logTlow, logrho, op_lowT, op2_lowT)
              first_lowT = .false.
           end if
           logkapx = splin2(logTlow, logrho, op_lowT, op2_lowT, logTx_eff, logrho_eff)
        else
           jx1 = NTlow
           jx2 = Nrho
           call hunt(logTlow, NTlow, logTx_eff, jx1)
           call hunt(logrho, Nrho, logrho_eff, jx2)
           call polin2mod(logTlow, logrho, op_lowT, jx1, jx2, logTx_eff, logrho_eff, logkapx, dkapx)
        end if
        kap_low = 10.0_dp**logkapx

        ! --- AESOPUS side (mid) ---
        logR_A = logR_phys
        if (logR_A > logRmax_A_tab .or. logR_A < logRmin_tab) then
           logR_A = clamp_interior(logR_A, logRmin_tab, logRmax_A_tab, eps_edge)
           ierrX = 1
        end if

        if (intpol == 'spline') then
           if (first_midT) then
              allocate(op2_midT(NTmid, NRmid), stat=alloc_error)
              call splie2(logTmid, logRmid, op_midT, op2_midT)
              first_midT = .false.
           end if
           logkapx = splin2(logTmid, logRmid, op_midT, op2_midT, logTx_eff, logR_A)
        else
           jx1 = NTmid
           jx2 = NRmid
           call hunt(logTmid, NTmid, logTx_eff, jx1)
           call hunt(logRmid, NRmid, logR_A, jx2)
           call polin2mod(logTmid, logRmid, op_midT, jx1, jx2, logTx_eff, logR_A, logkapx, dkapx)
        end if
        kap_mid = 10.0_dp**logkapx

        w = (logTcrit_AS_tab + dlogT_AS_tab - logTx_eff) / (2.0_dp*dlogT_AS_tab)
        w = max(0.0_dp, min(1.0_dp, w))
        kappa = w*kap_low + (1.0_dp - w)*kap_mid

        ! For debug
        !write (*, '("tempx =", 1pe12.4, " -> logTx_eff =", 1pe12.4, " -> ", a)') &
        !      tempx, logTx_eff, 'AESOPUS-Semenov'

     else
        ! Blend AESOPUS (mid) and OPAL (high) around Tcrit_OA_tab.
        ! logT = logTcrit_OA_tab - dlogT_OA_tab -> pure AESOPUS (w = 1)
        ! logT = logTcrit_OA_tab + dlogT_OA_tab -> pure OPAL    (w = 0)

        ! --- AESOPUS side (mid) ---
        logR_A = logR_phys
        if (logR_A > logRmax_A_tab .or. logR_A < logRmin_tab) then
           logR_A = clamp_interior(logR_A, logRmin_tab, logRmax_A_tab, eps_edge)
           ierrX = 1
        end if

        if (intpol == 'spline') then
           if (first_midT) then
              allocate(op2_midT(NTmid, NRmid), stat=alloc_error)
              call splie2(logTmid, logRmid, op_midT, op2_midT)
              first_midT = .false.
           end if
           logkapx = splin2(logTmid, logRmid, op_midT, op2_midT, logTx_eff, logR_A)
        else
           jx1 = NTmid
           jx2 = NRmid
           call hunt(logTmid, NTmid, logTx_eff, jx1)
           call hunt(logRmid, NRmid, logR_A, jx2)
           call polin2mod(logTmid, logRmid, op_midT, jx1, jx2, logTx_eff, logR_A, logkapx, dkapx)
        end if
        kap_mid = 10.0_dp**logkapx

        ! --- OPAL side (high) ---
        logR_OP = logR_phys
        if (logR_OP > logRmax_OF_tab .or. logR_OP < logRmin_tab) then
           logR_OP = clamp_interior(logR_OP, logRmin_tab, logRmax_OF_tab, eps_edge)
           ierrX = 1
        end if

        if (intpol == 'spline') then
           if (first_highT) then
              allocate(op2_highT(NThigh, NRhigh), stat=alloc_error)
              call splie2(logThigh, logRhigh, op_highT, op2_highT)
              first_highT = .false.
           end if
           logkapx = splin2(logThigh, logRhigh, op_highT, op2_highT, logTx_eff, logR_OP)
        else
           jx1 = NThigh
           jx2 = NRhigh
           call hunt(logThigh, NThigh, logTx_eff, jx1)
           call hunt(logRhigh, NRhigh, logR_OP, jx2)
           call polin2mod(logThigh, logRhigh, op_highT, jx1, jx2, logTx_eff, logR_OP, logkapx, dkapx)
        end if
        kap_high = 10.0_dp**logkapx

        w = (logTcrit_OA_tab + dlogT_OA_tab - logTx_eff) / (2.0_dp*dlogT_OA_tab)
        w = max(0.0_dp, min(1.0_dp, w))
        kappa = w*kap_mid + (1.0_dp - w)*kap_high

        ! For debug
        !write (*, '("tempx =", 1pe12.4, " -> logTx_eff =", 1pe12.4, " -> ", a)') &
        !         tempx, logTx_eff, 'OP-AESOPUS'
     end if
  end if

  !---------------------------------------------------------
  ! 3) Report error code (your convention)
  !---------------------------------------------------------
  if (ierrT == 1 .and. ierrX == 1) then
     ierror = 3
  else if (ierrT == 1) then
     ierror = 1
  else if (ierrX == 1) then
     ierror = 2
  else
     ierror = 0
  end if

END SUBROUTINE findkappa_OP_AES_S03

  
! ***********************************************************
  SUBROUTINE findkappa_OP_F05_S03(rhox, tempx, kappa, ierror)
! ***********************************************************
      use kind_params, only : i4b, dp, lgt
      use op_params, only : op_highT, op_midT, op_lowT, &
           op2_highT, op2_midT, op2_lowT, &
           logThigh, logRhigh, logTmid, logRmid, &
           logTlow, logrho, NThigh, NTmid, NTlow, &
           NRhigh, NRmid, Nrho, kes
      use spline2D_pac
      use polint2D_pac
      implicit none

      !character(len=3), parameter :: testing = 'on'
      character(len=3),  parameter :: testing = 'off'

      integer(i4b), intent(out) :: ierror
      integer(i4b)              :: alloc_error
      integer(i4b)              :: jx1, jx2
      real(dp), intent(in)      :: tempx, rhox
      real(dp), intent(out)     :: kappa

      ! Opacity table parameters
      real(dp), parameter :: dlogT_OF_tab = 0.05_dp, dlogT_FS_tab = 0.05_dp

      integer(i4b), save      :: n_warn = 0
      integer(i4b), parameter :: max_warn = 100

      real(dp) :: logRx, logrhox, logTx, logkapx, dkapx
      real(dp) :: logTcrit_OF_tab, logTcrit_FS_tab
      real(dp) :: kap_high, kap_mid, kap_low, w

      logical(lgt), save :: first_highT = .true., &
           first_midT  = .true., &
           first_lowT  = .true.

      ierror = 0

      !-- Optional sanity check: force Thomson opacity for debugging.
      IF (testing == 'on') THEN
         kappa = kes
         RETURN
      END IF

      logrhox   = LOG10(rhox)
      logTx     = LOG10(tempx)
      logRx     = logrhox - 3.0_dp*logTx + 18.0_dp
      logTcrit_OF_tab = LOG10(Tcrit_OF_tab)
      logTcrit_FS_tab = LOG10(Tcrit_FS_tab)

      !-----------------------------------------------------------------
      ! Temperature regions:
      !   logT >= logTcrit_OF_tab + dlogT_OF_tab : pure OPAL        (high-T)
      !   logT <= logTcrit_FS_tab - dlogT_FS_tab : pure Semenov     (low-T)
      !   logTcrit_FS_atb + dlogT_FS_tab <= logT <= logTcrit_OF_tab - dlogT_OF_tab : pure Ferguson
      !   otherwise : smooth blends across the two transition regions.
      !-----------------------------------------------------------------

      IF (logTx >= logTcrit_OF_tab + dlogT_OF_tab) THEN
         !---------------- High-temperature branch: OPAL only ----------------
         if (tempx > thigh_tab) then
            ! Temperature is above the valid range of OP tables.
            if (iprint >= 1 .and. n_warn < max_warn) then
               write (*,"('Temperature (', 1pe13.5, &
                 &       ') is too high to apply OP table.')") tempx
               n_warn = n_warn + 1
               if (n_warn == max_warn) then
                  write (*,*) 'Further OP/Ferguson/Semenov opacity warnings suppressed...'
               end if
            end if
            ierror = 1
            return
         else if (logRx > logRmax_OF_tab .or. logRx < logRmin_tab) then
            ! R (= rhox/(tempx/1.0e6)^3) is outside the valid range of OP tables.
            if (iprint >= 1 .and. n_warn < max_warn) then
               write (*,'("R (=rho/T6^3)",1pe13.5,") is out of range to apply ", &
                 &       "OP tables (rho =", 1pe13.5, ", T = ",1pe13.5,")")') &
                        10.0_dp**logRx, rhox, tempx
               n_warn = n_warn + 1
               if (n_warn == max_warn) then
                  write (*,*) 'Further OP/Ferguson/Semenov opacity warnings suppressed...'
               end if
            end if
            ierror = 2
            return
         end if

         IF (intpol == 'spline') THEN
            IF (first_highT) THEN
               ! Initialize spline for high-T OPAL table (one-time).
               ALLOCATE (op2_highT(NThigh, NRhigh), STAT=alloc_error)
               CALL splie2(logThigh, logRhigh, op_highT, op2_highT)
               first_highT = .false.
            END IF
            logkapx = splin2(logThigh, logRhigh, op_highT, op2_highT, &
                 logTx, logRx)
         ELSE
            jx1 = NThigh
            jx2 = NRhigh
            CALL hunt(logThigh, NThigh, logTx, jx1)
            CALL hunt(logRhigh, NRhigh, logRx, jx2)
            CALL polin2mod(logThigh, logRhigh, op_highT, &
                 jx1, jx2, logTx, logRx, logkapx, dkapx)
         END IF
         kappa = 10.0_dp**logkapx

      ELSE IF (logTx <= logTcrit_FS_tab - dlogT_FS_tab) THEN
         !---------------- Low-temperature branch: Semenov only ---------------
         IF (tempx < tlow_tab) THEN
            if (iprint >= 1 .and. n_warn < max_warn) then
               WRITE (*,"('Temperature (', 1pe13.5, ') is too low to ', &
                    & 'apply Semenov''s model.')") tempx
               n_warn = n_warn + 1
               if (n_warn == max_warn) then
                  write (*,*) 'Further OP/Ferguson/Semenov opacity warnings suppressed...'
               end if
            end if
            ierror = 1
            !-- temporary setting for plotting log kappa
            !kappa = 1.0e-99_dp
            RETURN
         ELSE IF (rhox > rhomax_tab .OR. rhox < rhomin_tab) THEN
            if (iprint >= 1 .and. n_warn < max_warn) then
               WRITE (*,"('Density (', 1pe12.5, ') is out of range to ', &
                    & 'apply Semenov''s model (T =', 1pe12.5, ')')") &
                    rhox, tempx
               n_warn = n_warn + 1
               if (n_warn == max_warn) then
                  write (*,*) 'Further OP/Ferguson/Semenov opacity warnings suppressed...'
               end if
            end if
            ierror = 2
            !-- temporary setting for plotting log kappa
            !kappa = 1.0e-99_dp
            RETURN
         END IF

         IF (intpol == 'spline') THEN
            IF (first_lowT) THEN
               ! Initialize spline for low-T Semenov table (one-time).
               ALLOCATE (op2_lowT(NTlow, Nrho), STAT=alloc_error)
               CALL splie2(logTlow, logrho, op_lowT, op2_lowT)
               first_lowT = .false.
            END IF
            logkapx = splin2(logTlow, logrho, op_lowT, op2_lowT, &
                 logTx, logrhox)
         ELSE
            jx1 = NTlow
            jx2 = Nrho
            CALL hunt(logTlow, NTlow, logTx, jx1)
            CALL hunt(logrho, Nrho, logrhox, jx2)
            CALL polin2mod(logTlow, logrho, op_lowT, &
                 jx1, jx2, logTx, logrhox, logkapx, dkapx)
         END IF
         kappa = 10.0_dp**logkapx

      ELSE IF (logTx >= logTcrit_FS_tab + dlogT_FS_tab .AND. &
               logTx <= logTcrit_OF_tab - dlogT_OF_tab) THEN
         !---------------- Intermediate branch: Ferguson only -----------------
         if (logRx > logRmax_OF_tab .or. logRx < logRmin_tab) then
            ! R (= rhox/(tempx/1.0e6)^3) is outside the valid range of OP tables.
            if (iprint >= 1 .and. n_warn < max_warn) then
               write (*,'("R (=rho/T6^3)",1pe13.5,") is out of range to apply ", &
                 &       "Ferguson table (rho =", 1pe13.5, ", T = ",1pe13.5,")")') &
                        10.0_dp**logRx, rhox, tempx
               n_warn = n_warn + 1
               if (n_warn == max_warn) then
                  write (*,*) 'Further OP/Ferguson/Semenov opacity warnings suppressed...'
               end if
            end if
            ierror = 2
            return
         end if

         IF (intpol == 'spline') THEN
            IF (first_midT) THEN
               ! Initialize spline for intermediate-T Ferguson table (one-time).
               ALLOCATE (op2_midT(NTmid, NRmid), STAT=alloc_error)
               CALL splie2(logTmid, logRmid, op_midT, op2_midT)
               first_midT = .false.
            END IF
            logkapx = splin2(logTmid, logRmid, op_midT, op2_midT, &
                 logTx, logRx)
         ELSE
            jx1 = NTmid
            jx2 = NRmid
            CALL hunt(logTmid, NTmid, logTx, jx1)
            CALL hunt(logRmid, NRmid, logRx, jx2)
            CALL polin2mod(logTmid, logRmid, op_midT, &
                 jx1, jx2, logTx, logRx, logkapx, dkapx)
         END IF
         kappa = 10.0_dp**logkapx

      ELSE

         !---------------- Transition regions: smooth blends ------------------
         IF (logTx <= logTcrit_FS_tab + dlogT_FS_tab) THEN

            !---- Blend Semenov (low) and Ferguson (mid) around Tcrit_FS_tab.
            !     logT = logTcrit_FS_tab - dlogT_FS_tab -> pure Semenov (w = 1)
            !     logT = logTcrit_FS_tab + dlogT_FS_tab -> pure Ferguson (w = 0)

            ! Low-T opacity from Semenov.
            IF (tempx < tlow_tab) THEN
               if (iprint >= 1 .and. n_warn < max_warn) then
                  WRITE (*,"('Temperature (', 1pe13.5, ') is too low to ', &
                       & 'apply Semenov''s model.')") tempx
                  n_warn = n_warn + 1
                  if (n_warn == max_warn) then
                     write (*,*) 'Further OP/Ferguson/Semenov opacity warnings suppressed...'
                  end if
               end if
               ierror = 1
               RETURN
            ELSE IF (rhox > rhomax_tab .OR. rhox < rhomin_tab) THEN
               if (iprint >= 1 .and. n_warn < max_warn) then
                  WRITE (*,"('Density (', 1pe12.5, ') is out of range to ', &
                      & 'apply Semenov''s model (T =', 1pe12.5, ')')") &
                      rhox, tempx
                  n_warn = n_warn + 1
                  if (n_warn == max_warn) then
                     write (*,*) 'Further OP/Ferguson/Semenov opacity warnings suppressed...'
                  end if
               end if
               ierror = 2
               RETURN
            END IF

            IF (intpol == 'spline') THEN
               IF (first_lowT) THEN
                  ! Initialize spline for low-T Semenov table (one-time).
                  ALLOCATE (op2_lowT(NTlow, Nrho), STAT=alloc_error)
                  CALL splie2(logTlow, logrho, op_lowT, op2_lowT)
                  first_lowT = .false.
               END IF
               logkapx = splin2(logTlow, logrho, op_lowT, op2_lowT, &
                    logTx, logrhox)
            ELSE
               jx1 = NTlow
               jx2 = Nrho
               CALL hunt(logTlow, NTlow, logTx, jx1)
               CALL hunt(logrho, Nrho, logrhox, jx2)
               CALL polin2mod(logTlow, logrho, op_lowT, &
                    jx1, jx2, logTx, logrhox, logkapx, dkapx)
            END IF
            kap_low = 10.0_dp**logkapx

            ! Intermediate-T opacity from Ferguson.
            IF (intpol == 'spline') THEN
               IF (first_midT) THEN
                  ! Initialize spline for intermediate-T Ferguson table (one-time).
                  ALLOCATE (op2_midT(NTmid, NRmid), STAT=alloc_error)
                  CALL splie2(logTmid, logRmid, op_midT, op2_midT)
                  first_midT = .false.
               END IF
               logkapx = splin2(logTmid, logRmid, op_midT, op2_midT, &
                    logTx, logRx)
            ELSE
               jx1 = NTmid
               jx2 = NRmid
               CALL hunt(logTmid, NTmid, logTx, jx1)
               CALL hunt(logRmid, NRmid, logRx, jx2)
               CALL polin2mod(logTmid, logRmid, op_midT, &
                    jx1, jx2, logTx, logRx, logkapx, dkapx)
            END IF
            kap_mid = 10.0_dp**logkapx

            w = (logTcrit_FS_tab + dlogT_FS_tab - logTx)/(2.0_dp*dlogT_FS_tab)
            w = MAX(0.0_dp, MIN(1.0_dp, w))
            kappa = w*kap_low + (1.0_dp - w)*kap_mid

         ELSE

            !---- Blend Ferguson (mid) and OPAL (high) around Tcrit_OF.
            !     logT = logTcrit_OF - dlogT_OF_tab -> pure Ferguson (w = 1)
            !     logT = logTcrit_OF + dlogT_OF_tab -> pure OPAL    (w = 0)

            ! Intermediate-T opacity from Ferguson.
            IF (intpol == 'spline') THEN
               IF (first_midT) THEN
                  ! Initialize spline for intermediate-T Ferguson table (one-time).
                  ALLOCATE (op2_midT(NTmid, NRmid), STAT=alloc_error)
                  CALL splie2(logTmid, logRmid, op_midT, op2_midT)
                  first_midT = .false.
               END IF
               logkapx = splin2(logTmid, logRmid, op_midT, op2_midT, &
                    logTx, logRx)
            ELSE
               jx1 = NTmid
               jx2 = NRmid
               CALL hunt(logTmid, NTmid, logTx, jx1)
               CALL hunt(logRmid, NRmid, logRx, jx2)
               CALL polin2mod(logTmid, logRmid, op_midT, &
                    jx1, jx2, logTx, logRx, logkapx, dkapx)
            END IF
            kap_mid = 10.0_dp**logkapx

            ! High-T opacity from OPAL.
            IF (intpol == 'spline') THEN
               IF (first_highT) THEN
                  ! Initialize spline for high-T OPAL table (one-time).
                  ALLOCATE (op2_highT(NThigh, NRhigh), STAT=alloc_error)
                  CALL splie2(logThigh, logRhigh, op_highT, op2_highT)
                  first_highT = .false.
               END IF
               logkapx = splin2(logThigh, logRhigh, op_highT, op2_highT, &
                    logTx, logRx)
            ELSE
               jx1 = NThigh
               jx2 = NRhigh
               CALL hunt(logThigh, NThigh, logTx, jx1)
               CALL hunt(logRhigh, NRhigh, logRx, jx2)
               CALL polin2mod(logThigh, logRhigh, op_highT, &
                    jx1, jx2, logTx, logRx, logkapx, dkapx)
            END IF
            kap_high = 10.0_dp**logkapx

            w = (logTcrit_OF_tab + dlogT_OF_tab - logTx)/(2.0_dp*dlogT_OF_tab)
            w = MAX(0.0_dp, MIN(1.0_dp, w))
            kappa = w*kap_mid + (1.0_dp - w)*kap_high

         END IF

      END IF
      
      END SUBROUTINE findkappa_OP_F05_S03

end module opacity_table_mod

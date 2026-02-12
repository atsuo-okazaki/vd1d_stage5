!***********************************************************************
!>  MODULE: run_control
!-----------------------------------------------------------------------
!  Purpose:
!     Provide run-time control parameters such as iprint and runid
!
!  Public entities:
!     (Add exported types, procedures, and parameters)
!
!  Notes:
!     (Add any remarks or design decisions)
!-----------------------------------------------------------------------
!  Author: A. Okazaki
!  Created: 2025-11-09
!***********************************************************************

!**********************************************
      MODULE run_control
!**********************************************
  use kind_params, only : i4b
  implicit none
  save

  !-------------------------------------------------
  ! Verbosity level for run-time messages
  !
  !  0 : minimal (fatal errors + tiny summary)
  !  1 : normal (summary + basic info)
  !  2 : verbose (diagnostics, detailed warnings)
  !-------------------------------------------------
  integer(i4b) :: iprint = 1

  !-------------------------------------------------
  ! Run identifier that is attached to output files.
  !
  ! Examples:
  !   runid = 'M100_q1.0'
  !   runid = 'test'
  !
  ! If runid is empty, it is simply skipped in filenames.
  !-------------------------------------------------
  character(len=64) :: runid = 'default'

 END MODULE run_control

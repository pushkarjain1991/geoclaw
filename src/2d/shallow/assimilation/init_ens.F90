!$Id: init_ens.F90 1443 2013-10-04 10:52:09Z lnerger $
!BOP
!
! !ROUTINE: init_ens --- Initialize ensemble
!
! !INTERFACE:
SUBROUTINE init_ens(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! The routine is called when the filter is
! initialized in PDAF\_filter\_init.  It has
! to initialize an ensemble of dim\_ens states.
! Typically, the ensemble will be directly read from files.
!
! The routine is called by all filter processes and
! initializes the ensemble for the PE-local domain.
!
! Implementation for the 2D online example
! without parallelization.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code based on offline_1D
! Later revisions - see svn log
!
! !USES:
  USE mod_model, ONLY: reshaped_recv_ic, field
  use mod_assimilation, only: type_ensinit
  

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: filtertype              ! Type of filter to initialize
  INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens                 ! Size of ensemble
  REAL, INTENT(inout) :: state_p(dim_p)          ! PE-local model state
  ! It is not necessary to initialize the array 'state_p' for SEIK.
  ! It is available here only for convenience and can be used freely.
  REAL, INTENT(inout) :: Uinv(dim_ens-1,dim_ens-1) ! Array not referenced for SEIK
  REAL, INTENT(out)   :: ens_p(dim_p, dim_ens)   ! PE-local state ensemble
  INTEGER, INTENT(inout) :: flag                 ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: PDAF_filter_init    (as U_ens_init)
!EOP

! *** local variables ***
  INTEGER :: i, j, member  ! Counters
  INTEGER, SAVE :: allocflag = 0      ! Flag for memory counting
  CHARACTER(len=3) :: ensstr          ! String for ensemble member


! **********************
! *** INITIALIZATION ***
! **********************

  ! *** Generate full ensemble on filter-PE 0 ***
  WRITE (*, '(/9x, a)') 'Initialize state ensemble'
  WRITE (*, '(9x, a)') '--- read ensemble from files'
  WRITE (*, '(9x, a, i5)') '--- Ensemble size:  ', dim_ens
!
  
  state_p = field
  IF (TRIM(type_ensinit) == 'eof') THEN
     ! Initialize by 2nd-order exact sampling from EOFs
     print *, "type_ensinit is EOF"
     CALL init_ens_eof(dim_p, dim_ens, state_p, ens_p, flag)
!  ELSE IF (TRIM(type_ensinit) == 'rnd') THEN
!     ! Initialize by random sampling from state trajectory
!     CALL init_ens_rnd(dim_p, dim_ens, state_p, ens_p, flag)
   ELSE 
      stop "No type_ensinit specified"
  END IF


!DO member = 1, dim_ens
!    ens_p(:, member) = reshaped_recv_ic(:,member)
!    print *, "shape of ens_p", shape(ens_p)
!enddo
    
! ****************
! *** clean up ***
! ****************
!deallocate(reshaped_recv_ic)

END SUBROUTINE init_ens

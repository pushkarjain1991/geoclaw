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
  USE mod_model, &
       ONLY: nx, ny

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
  REAL, ALLOCATABLE :: field(:)     ! global model field
  CHARACTER(len=2) :: ensstr          ! String for ensemble member
  INTEGER :: index_2d_pdaf(2)
  !INTEGER :: index_2d_row, index_2d_col
  INTEGER :: index_1d_pdaf


! **********************
! *** INITIALIZATION ***
! **********************

  ! *** Generate full ensemble on filter-PE 0 ***
  WRITE (*, '(/9x, a)') 'Initialize state ensemble'
  WRITE (*, '(9x, a)') '--- read ensemble from files'
  WRITE (*, '(9x, a, i5)') '--- Ensemble size:  ', dim_ens
  
  ! allocate memory for temporary fields
  ALLOCATE(field(ny*nx))
 !print *, nx, ny

! ********************************
! *** Read ensemble from files ***
! ********************************

  DO member = 1, dim_ens
     WRITE (ensstr, '(i2)') member
     OPEN(24, file = '../ens_'//TRIM(ADJUSTL(ensstr))//'.txt', status='old')
     !do i=1,4
         !do j =1,nx*ny
     !    do j = 1,nx1*ny1
     !        read(24,*) field((i-1)*nx1*ny1 + j)
     !    enddo
     !enddo
     read(24,*)field
     !print *,field
     CLOSE(24)
     ens_p(:,member) = field(:)
     !DO i = 1, ny
     !   READ (20, *) field(i, :)
     !END DO
     !DO i = 1, ny
     !   ens_p(1 + (i-1)*nx : i*nx, member) = field(i, 1:nx)
     !END DO

  END DO


! ****************
! *** clean up ***
! ****************

  DEALLOCATE(field)

  !Just testing printing
!  call oned_to_twod(3200,2,2,index_2d_pdaf)
!  print *, "2d domain is", index_2d_pdaf
!  call oned_to_twod(5445,2,2,index_2d_pdaf)
!  print *, "2d domain is", index_2d_pdaf
!  call oned_to_twod(5720,2,2,index_2d_pdaf)
!  print *, "2d domain is", index_2d_pdaf
!  call oned_to_twod(10000,2,2,index_2d_pdaf)
!  print *, "2d domain is", index_2d_pdaf

!  call twod_to_oned(14,100,index_1d_pdaf)
!  print *, "1d domain is", index_1d_pdaf
!  call twod_to_oned(59,45,index_1d_pdaf)
!  print *, "1d domain is", index_1d_pdaf
!  call twod_to_oned(65,20,index_1d_pdaf)
!  print *, "1d domain is", index_1d_pdaf
!  call twod_to_oned(100,100,index_1d_pdaf)
!  print *, "1d domain is", index_1d_pdaf

END SUBROUTINE init_ens

!$Id: init_ens_eof.F90 1606 2016-05-30 15:23:06Z lnerger $
!BOP
!
! !ROUTINE: init_ens_eof --- Initialize ensemble from EOF decomposition
!
! !INTERFACE:
SUBROUTINE init_ens_eof(dim, dim_ens, state, ens, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF (SEIK):
!
! The routine is called by init_seik. It 
! initializes an ensemble of dim\_ens states
! by exact 2nd order sampling.
! State vectors of the form
!   $x_i = x + sqrt(FAC) eofV (\Omega C^{-1})^T$
! fulfill the condition
!   $P = 1/(FAC)  \sum_{i=1}^{dim\_ens} (x_i - x)(x_i - x)^T$
! The matrix is initialized in the form of
! singular values and singular vectors.
!
! This version is for the Lorenz96 model
! without parallelization.
!
! !REVISION HISTORY:
! 2009-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
use mod_model, only: field
  IMPLICIT NONE


! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim                 ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens             ! Size of ensemble
  REAL, INTENT(inout) :: state(dim)          ! PE-local model state
  ! It is not necessary to initialize the array 'state' for SEIK. 
  ! It is available here only for convenience and can be used freely.
  REAL, INTENT(out)   :: ens(dim, dim_ens)   ! PE-local state ensemble
  INTEGER, INTENT(inout) :: flag             ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: init_seik
! Calls: seik_omega
! Calls: timeit
! Calls: memcount
! Calls: dgemm (BLAS)
!EOP

! *** local variables ***
  INTEGER :: i, s, row, col       ! counters
  INTEGER, SAVE :: allocflag = 0  ! Flag for memory counting
  REAL, ALLOCATABLE :: eofV(:,:)  ! matrix of eigenvectors V 
  REAL, ALLOCATABLE :: svals(:)   ! singular values
  REAL, ALLOCATABLE :: omega(:,:) ! Matrix Omega
  REAL :: fac                     ! Square-root of dim_eof+1 or dim_eof
  INTEGER :: dim_file             ! State dimension in file
  INTEGER :: rank                 ! Rank of approximated covariance matrix
  INTEGER :: rank_file            ! Rank of covariance matrix stored in file

! *** Just for debugging ***
  REAL :: zero_mean(dim)
  REAL :: obs_matrix(dim, dim_ens)   ! PE-local state ensemble
  CHARACTER(len=3) :: ensstr

! **********************
! *** INITIALIZATION ***
! **********************

  ! *** Rank of matrix is ensemble size minus one
  rank = dim_ens - 1
  
  ! *** Generate full ensemble on filter-PE 0 ***
  WRITE (*, '(9x, a)') '--- generate ensemble from covariance matrix'
  WRITE (*, '(9x, a)') &
       '--- use rank reduction and 2nd order exact sampling (SEIK type)'
  WRITE (*, '(9x, a, i5)') '--- Ensemble size:  ', dim_ens
  WRITE (*, '(9x, a, i5)') '--- number of EOFs: ', rank

  ! allocate memory for temporary fields
  ALLOCATE(eofV(dim, rank))
  ALLOCATE(svals(rank))
  ALLOCATE(omega(rank + 1, rank))


! *************************************************
! *** Initialize initial state and covar matrix ***
! *************************************************
  
  WRITE(*,'(9x,a,a)') '--- Reading covariance information'


  ! Read singular values
  open(22, file='../singular_values',status='old')
  read(22,*) svals(:)
  close(22)


!  ! Read mean state
!  open(22, file='mean_state',status='old')
!  read(22, *) meanstate(:)
!  close(22)
  state(:) = field(:)

  ! *** Read singular vectors

  open(22, file='../singular_vector', status='old')
  readvectors: DO i = 1, rank
    ! Read singular vector
    read(22,*) eofV(:,i)
  END DO readvectors
  close(22)


! *****************************************
! *** Generate ensemble of model states ***
! *****************************************

  WRITE (*,'(9x, a)') '--- generate state ensemble'
     
  ! Use PDAF routine to generate ensemble from covariance matrix
  CALL PDAF_SampleEns(dim, dim_ens, eofV, svals, state, ens, flag)
     
! *****************************************
! *** JUST FOR DEBUGGING ***
! *** Checking magnitude of perturbations ***
! *****************************************
  zero_mean(:) = 0.d0
  CALL PDAF_SampleEns(dim, dim_ens, eofV, svals, zero_mean, obs_matrix, flag)
  open(22, file='check_obs.txt', status='replace')
  do i=1,dim_ens
    write(22,*) obs_matrix(:,i)
  enddo
  close(22)

  do i=1,dim_ens
      write(ensstr,'(i3.1)') i
      open(22, file="init_ens" // trim(adjustl(ensstr)),status='new')
      write(22,*) ens(:,i)
      close(22)
  enddo
  

! ****************
! *** clean up ***
! ****************

  DEALLOCATE(svals, eofV, omega)

END SUBROUTINE init_ens_eof

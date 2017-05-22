!$Id: init_dim_obs_pdaf.F90 1443 2013-10-04 10:52:09Z lnerger $
!BOP
!
! !ROUTINE: init_dim_obs_pdaf --- Compute number of observations
!
! !INTERFACE:
SUBROUTINE init_dim_obs_pdaf(step, dim_obs_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/ETKF/ESTKF
!
! The routine is called at the beginning of each
! analysis step.  It has to initialize the size of 
! the observation vector according to the current 
! time step for the PE-local domain.
!
! Implementation for the 2D online example
! without parallelization.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE

! ARGUMENTS:
  INTEGER, INTENT(in)  :: step      ! Current time step
  INTEGER, INTENT(inout) :: dim_obs_p ! Dimension of full observation vector
! LOCAL variables
  CHARACTER(len=20) :: fname
  CHARACTER(len=4) :: stepstr         ! String for time step
  REAL(KIND=8):: x,y
  REAL(KIND=8),ALLOCATABLE :: xx(:),yy(:),qq(:)
  INTEGER :: i,IO
  INTEGER :: cnt0


  INTERFACE
      SUBROUTINE get_obs(xx,yy,qq,dim_obs_p)
          REAL(KIND=8), intent(in) :: xx(:),yy(:),qq(:)
          INTEGER, intent(inout) :: dim_obs_p
      END SUBROUTINE get_obs
  END INTERFACE



! !CALLING SEQUENCE:
! Called by: PDAF_seek_analysis    (as U_init_dim_obs)
! Called by: PDAF_seik_analysis, PDAF_seik_analysis_newT
! Called by: PDAF_enkf_analysis_rlm, PDAF_enkf_analysis_rsm
! Called by: PDAF_etkf_analysis, PDAF_etkf_analysis_T
! Called by: PDAF_estkf_analysis, PDAF_estkf_analysis_fixed
!EOP

   print *, "Running init_dim_obs_pdaf"

   WRITE (stepstr, '(i4)') step
   fname='../obs_step'//TRIM(ADJUSTL(stepstr))//'.txt'
   !fname=trim(adjustl(fname))

   ! Counting the number of observation values
   cnt0=0
   OPEN (12, file= fname, status='old')
   DO
      READ(12,*,IOSTAT=IO) x,y
      IF (IO /= 0) exit 
      cnt0=cnt0+1
   ENDDO 
   if (IO>0) then
       print *,"error: error reading the file ",fname
       stop
   endif
   print *, "Obs file ", fname, " has ", cnt0, " observations"

   ! Allocate space for x,y coordinates and observations q
   allocate(xx(cnt0))
   allocate(yy(cnt0))
   allocate(qq(cnt0))
   rewind(12)
   do i=1,cnt0
       read(12,*) xx(i),yy(i),qq(i)
   enddo
   close(12)

   call get_obs(xx,yy,qq,dim_obs_p) 
   print *, "dim_obs_p = ", dim_obs_p

   print *, "Finished running init_dim_obs_pdaf"

END SUBROUTINE init_dim_obs_pdaf

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
  INTEGER :: dim_obs_p ! Dimension of full observation vector
! LOCAL variables
  CHARACTER(len=20) :: fname
  CHARACTER(len=4) :: stepstr         ! String for time step
  REAL(KIND=8):: x,y
  REAL(KIND=8),ALLOCATABLE :: xx(:),yy(:),qq(:)
  INTEGER :: i,IO
  INTEGER :: cnt0


  INTERFACE
      SUBROUTINE get_obs(xx,yy,qq,dim_obs_p)
          REAL(KIND=8),ALLOCATABLE :: xx(:),yy(:),qq(:)
          INTEGER :: dim_obs_p
      END SUBROUTINE get_obs
  END INTERFACE



! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_init_dim_obs)
! Called by: PDAF_lestkf_update  (as U_init_dim_obs)
! Called by: PDAF_letkf_update   (as U_init_dim_obs)
!EOP

   WRITE (stepstr, '(i4)') step
   fname='../obs_step'//TRIM(ADJUSTL(stepstr))//'.txt'
   fname=trim(adjustl(fname))
   cnt0=0
   OPEN (12, file= fname,status='old')
 
   DO
      READ(12,*,IOSTAT=IO) x,y
      IF (IO /= 0) exit 
      cnt0=cnt0+1
   ENDDO 
   if (IO>0) then
       print *,"error: error reading the file ",fname
       stop
   endif

   allocate(xx(cnt0))
   allocate(yy(cnt0))
   allocate(qq(cnt0))
   rewind(12)
   do i=1,cnt0
       read(12,*) xx(i),yy(i),qq(i)
   enddo
   close(12)
!   print *,cnt0
!   print *,xx(1:10),yy(1:10),qq(1:10)
   print *, "init_dim_obs_pdaf1"
   call get_obs(xx,yy,qq,dim_obs_p) 
   print *, "init_dim_obs_pdaf2"





END SUBROUTINE init_dim_obs_pdaf

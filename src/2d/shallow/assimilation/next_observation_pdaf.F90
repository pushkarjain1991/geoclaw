!$Id: next_observation_pdaf.F90 1437 2013-10-04 07:20:38Z lnerger $
!BOP
!
! !ROUTINE: next_observation_pdaf --- Initialize information on next observation
!
! !INTERFACE:
SUBROUTINE next_observation_pdaf(stepnow, nsteps_pdaf, doexit, time)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! The subroutine is called before each forecast phase
! by PDAF\_get\_state. It has to initialize the number 
! of time steps until the next available observation 
! (nsteps_pdaf) and the current model time (time). In 
! addition the exit flag (exit) has to be initialized.
! It indicates if the data assimilation process is 
! completed such that the ensemble loop in the model 
! routine can be exited.
!
! The routine is called by all processes. 
!         
! Version for the dummy model. Identical for
! mode- and domain-decomposition .
!
! !REVISION HISTORY:
! 2004-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, ONLY: delt_obs
  USE mod_parallel, &
       ONLY: mype_world
  USE mod_model, &
       ONLY: total_steps

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: stepnow  ! Number of the current time step
  INTEGER, INTENT(out) :: nsteps_pdaf   ! Number of time steps until next obs
  INTEGER, INTENT(out) :: doexit   ! Whether to exit forecasting (1 for exit)
  REAL, INTENT(out)    :: time     ! Current model (physical) time

! !CALLING SEQUENCE:
! Called by: PDAF_get_state   (as U_next_obs)
!EOP



! *******************************************************
! *** Set number of time steps until next observation ***
! *******************************************************


  IF (stepnow + nsteps_pdaf <= total_steps) THEN
     ! Next obbservation ahead
     nsteps_pdaf = delt_obs
     print *, "delt_obs =", delt_obs
     print *, "stepnow = ", stepnow
     !have_obs = .TRUE.
     doexit = 0          ! Do not exit assimilation
     if (mype_world==0) WRITE (*, '(i7, 3x, a, i7)') &
          stepnow, 'Next observation at time step', stepnow + nsteps_pdaf
  
  ELSE
     nsteps_pdaf = 0
     doexit = 1
  !have_obs = .FALSE.
     
     ! End of assimilation
     if (mype_world==0) WRITE (*, '(i7, 3x, a)') &
     stepnow, 'No more observations - end assimilation'

  END IF

END SUBROUTINE next_observation_pdaf

!$Id: mod_parallel_pdaf.F90 1415 2013-09-25 14:33:26Z lnerger $
!BOP
!
! !MODULE:
MODULE mod_parallel_ens_gen

! !DESCRIPTION:
! This modules provides variables for the MPI parallelization
! to be shared between model-related routines. The are variables
! that are used in the model, even without PDAF and additional
! variables that are only used, if data assimialtion with PDAF
! is performed.
! In addition methods to initialize and finalize MPI are provided.
! The initialization routine is only for the model itself, the 
! more complex initialization of communicators for xecution with
! PDAF is peformed in init\_parallel\_pdaf.
!
! !REVISION HISTORY:
! 2004-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE
  SAVE 

  INCLUDE 'mpif.h'

! !PUBLIC DATA MEMBERS:
  ! Basic variables for model state integrations

  ! Additional variables for use with PDAF
  INTEGER :: dim_ens, mype_world   ! # PEs and PE rank in MPI_COMM_WORLD
  INTEGER :: MPIerr      ! Error flag for MPI
  INTEGER :: MPIstatus(MPI_STATUS_SIZE)       ! Status array for MPI
  LOGICAL :: second_valout=.false.
!EOP
  
CONTAINS
!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_parallel - Initialize MPI
!
! !INTERFACE:
  SUBROUTINE init_parallel()

! !DESCRIPTION:
! Routine to initialize MPI, the number of PEs
! (npes\_world) and the rank of a PE (mype\_world).
! The model is executed within the scope of the
! communicator Comm_model. It is also initialized
! here together with its size (npes\_model) and 
! the rank of a PE (mype\_model) within Comm_model.
!EOP

    IMPLICIT NONE

    INTEGER :: i
  
    CALL MPI_INIT(i);
    CALL MPI_Comm_Size(MPI_COMM_WORLD,dim_ens,i)
    CALL MPI_Comm_Rank(MPI_COMM_WORLD,mype_world,i)

  END SUBROUTINE init_parallel
!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: finalize_parallel - Finalize MPI
!
! !INTERFACE:
  SUBROUTINE finalize_parallel()

! !DESCRIPTION:
! Routine to finalize MPI
!EOP

    IMPLICIT NONE
    
    CALL  MPI_Barrier(MPI_COMM_WORLD,MPIerr)
    CALL  MPI_Finalize(MPIerr)

  END SUBROUTINE finalize_parallel
!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: abort_parallel - Abort MPI
!
! !INTERFACE:
  SUBROUTINE abort_parallel()

! !DESCRIPTION:
! Routine to abort MPI program
!EOP

    IMPLICIT NONE
    
    CALL  MPI_Abort(MPI_COMM_WORLD, 1, MPIerr)

  END SUBROUTINE abort_parallel

END MODULE mod_parallel_ens_gen

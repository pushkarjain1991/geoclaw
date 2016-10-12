!$Id: mod_model.F90 1409 2013-09-25 11:47:03Z lnerger $
!BOP
!
! !MODULE:
MODULE mod_model

! !DESCRIPTION:
! This module provides variables needed for the
! 2-dimensional tutorial model without parallelization.
!
! !REVISION HISTORY:
! 2013-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE
  SAVE
!EOP


! *** Variables specific for 2D tutorial model ***

  !INTEGER :: nx, ny               ! Size of 2D grid
  !REAL(KIND=8):: xlow,xhigh
  !REAL(KIND=8):: ylow,yhigh
  INTEGER :: total_steps          ! Total number of time steps
  !REAL(KIND=8), ALLOCATABLE :: field(:,:) ! Model field
  REAL(KIND=8), ALLOCATABLE :: field(:) ! Model field
  REAL(KIND=8), ALLOCATABLE :: reshaped_recv_ic(:,:) ! Model field
!  REAL :: dtinit                  ! Time step size
!  REAL :: time                    ! Model time

!  INTEGER :: num_grids
!  INTEGER, ALLOCATABLE :: nx_local(:)
!  INTEGER, ALLOCATABLE :: ny_local(:)
!  INTEGER :: numcells_pdaf
!  REAL(KIND=8),ALLOCATABLE :: rnode_pdaf(:,:)
!  INTEGER,ALLOCATABLE :: node_pdaf(:,:)
  !INTEGER,ALLOCATABLE :: mptr_array(:)
  !INTEGER,ALLOCATABLE :: ordered_mptr_array(:)
  
END MODULE mod_model

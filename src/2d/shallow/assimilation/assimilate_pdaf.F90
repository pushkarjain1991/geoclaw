!$Id: assimilate_pdaf.F90 1415 2013-09-25 14:33:26Z lnerger $
!BOP
!
! !ROUTINE: assimilate_pdaf - Routine to control perform analysis step
!
! !INTERFACE:
SUBROUTINE assimilate_pdaf(nvar,naux,mxlevel,time_geoclaw)

! !DESCRIPTION:
! This routine is called during the model integrations at each time 
! step. It check whether the forecast phase is completed. If so, 
! PDAF_put_state_X is called to perform the analysis step.
!
! !REVISION HISTORY:
! 2013-08 - Lars Nerger - Initial code for NEMO
! Later revisions - see svn log
!
! !USES:
!  use mod_model, only: field
!  use mod_assimilation, only: dim_ens, filename,
!     & filtertype, first_assimilation
!  use PDAF_mod_filter, only: cnt_steps,nsteps,type_filter
  use mod_parallel, only: mype_world,abort_parallel,mpi_comm_world,mpierr
  USE mod_assimilation, &      ! Variables for assimilation
       ONLY: filtertype, stepnow_pdaf,assimilate_step,ncycle_pdaf
   use mod_model,only: field

  IMPLICIT NONE

! !CALLING SEQUENCE:
! Called by: step
! CAlls: PDAF_assimilate_X
!EOP

! Local variables
  INTEGER,INTENT(IN)::nvar
  INTEGER,INTENT(IN)::naux
  INTEGER,INTENT(IN)::mxlevel
  INTEGER :: status_pdaf       ! PDAF status flag
  INTEGER :: steps
  INTEGER :: doexit
  INTEGER :: ii,jj
  REAL(8) :: time_geoclaw
  REAL(8) :: time_pdaf

  CHARACTER(len=3) :: cyclestr
  CHARACTER(len=30) :: fname
  LOGICAL :: there
  LOGICAL :: isoutput
!  INTEGER :: assimilate_level(2)


  ! External subroutines
  EXTERNAL :: prodRinvA_pdaf,&
             prodRinvA_l_pdaf,&
             next_observation_pdaf,&
      ! Provide time step, model time and dimension of next observation
             distribute_state_pdaf,&
      ! Routine to distribute a state vector to model fields
             prepoststep_ens_pdaf,&
      ! Routine to collect a state vector from model fields
             collect_state_pdaf,&
      ! Initialize dimension of observation vector
             init_dim_obs_pdaf,&
             init_dim_obs_f_pdaf,&
             init_dim_obs_l_pdaf,&
      ! Implementation of observation operator
             obs_op_pdaf,&
             obs_op_f_pdaf,&
      ! Routine to provide vector of measurements
             init_obs_pdaf,&
             init_obs_f_pdaf,&
             init_obs_l_pdaf,&
      ! Add obs. error covariance R to HPH in EnKF
             add_obs_error_pdaf,&
      ! Initialize obs error covar R in EnKF
             init_obscovar_pdaf,&
      !
             init_n_domains_pdaf,&
             init_dim_l_pdaf,&
             g2l_state_pdaf,&
             l2g_state_pdaf,&
             g2l_obs_pdaf,&
             init_obsvar_pdaf,&
             init_obsvar_l_pdaf
! *********************************
! *** Call assimilation routine ***
! *********************************
!  assimilate_level(1)=mxlevel;assimilate_level(2)=1
   !stepnow_pdaf=stepnow_pdaf+1 !Moving to tick so that can update flags at
   !assimilation
   ncycle_pdaf=ncycle_pdaf+1
!   call alloc2field(nvar,naux)
   if (stepnow_pdaf==assimilate_step) then
       

       call alloc2field(nvar,naux)
       !CALL  MPI_Barrier(mpi_comm_world,MPIerr)
       !print *, "mypeyhmm", mype_world
       !CALL  MPI_Barrier(mpi_comm_world,MPIerr)
       !stop
       !print *,'(1,1) before get state is',field(25311)
        if(mype_world==0) write(*,'(a,5x,a)') 'PDAF','Perform assimilation with PDAF '
        if (filtertype==2) then
             CALL PDAF_put_state_enkf(collect_state_pdaf,&
                 init_dim_obs_pdaf, obs_op_pdaf, init_obs_pdaf,&
                 prepoststep_ens_pdaf, add_obs_error_pdaf,&
                 init_obscovar_pdaf,status_pdaf)
        else if (filtertype==1) then
             CALL PDAF_put_state_seik(collect_state_pdaf,&
                 init_dim_obs_pdaf, obs_op_pdaf, init_obs_pdaf,&
                 prepoststep_ens_pdaf, prodRinvA_pdaf,&
                 init_obsvar_pdaf, status_pdaf)
        else if (filtertype==4) then
             CALL PDAF_put_state_etkf(collect_state_pdaf,&
                 init_dim_obs_pdaf, obs_op_pdaf, init_obs_pdaf,&
                 prepoststep_ens_pdaf, prodRinvA_pdaf,&
                 init_obsvar_pdaf,status_pdaf)
        else if (filtertype==0) then
            CALL PDAF_assimilate_seek(collect_state_pdaf,&
                init_dim_obs_pdaf, obs_op_pdaf, init_obs_pdaf,&
                prepoststep_ens_pdaf, prodRinvA_pdaf, status_pdaf)
        else if (filtertype==6) then
            CALL PDAF_put_state_estkf(collect_state_pdaf,&
                init_dim_obs_pdaf, obs_op_pdaf, init_obs_pdaf,&
                prepoststep_ens_pdaf, prodRinvA_pdaf,&
                init_obsvar_pdaf, status_pdaf)
        !else if (filtertype==7) then
        !    CALL PDAF_put_state_lestkf(collect_state_pdaf,&
        !        init_dim_obs_f_pdaf,obs_op_f_pdaf,init_obs_f_pdaf,&
        !        init_obs_l_pdaf,prepoststep_ens_pdaf,&
        !        prodRinvA_l_pdaf,init_n_domains_pdaf,&
        !        init_dim_l_pdaf,init_dim_obs_l_pdaf,g2l_state_pdaf,&
        !        l2g_state_pdaf, g2l_obs_pdaf, init_obsvar_pdaf,&
        !        init_obsvar_l_pdaf,status_pdaf)
        !else if (filtertype==3) then
        !    CALL PDAF_put_state_lseik(collect_state_pdaf,&
        !        init_dim_obs_f_pdaf,obs_op_f_pdaf,init_obs_f_pdaf,&
        !        init_obs_l_pdaf,prepoststep_ens_pdaf,&
        !        prodRinvA_l_pdaf,init_n_domains_pdaf,&
        !        init_dim_l_pdaf,init_dim_obs_l_pdaf,g2l_state_pdaf,&
        !        l2g_state_pdaf, g2l_obs_pdaf, init_obsvar_pdaf,&
        !        init_obsvar_l_pdaf,status_pdaf)
        !else if (filtertype==5) then
        !    CALL PDAF_put_state_letkf(collect_state_pdaf,&
        !        init_dim_obs_f_pdaf,obs_op_f_pdaf,init_obs_f_pdaf,&
        !        init_obs_l_pdaf,prepoststep_ens_pdaf,&
        !        prodRinvA_l_pdaf,init_n_domains_pdaf,&
        !        init_dim_l_pdaf,init_dim_obs_l_pdaf,g2l_state_pdaf,&
        !        l2g_state_pdaf, g2l_obs_pdaf, init_obsvar_pdaf,&
        !        init_obsvar_l_pdaf,status_pdaf)
        else
             if (mype_world==0) print *,"invalid filter type"
             CALL  abort_parallel()
        endif
!      call MPI_barrier(mpi_comm_world,mpierr)
        IF (status_pdaf==0) then
#ifdef dasd
    output : if (mype_world==0) then! output for mype=0
                 write(cyclestr,'(i3)') ncycle_pdaf
                 fname='state_step'//trim(adjustl(cyclestr))//'_ana.txt'
                 fname=trim(adjustl(fname))
!                 print *, 'the size of field is ', size(field)
!                 INQUIRE(file='state_step'//trim(adjustl(cyclestr))//'_ana.txt',exist=there)
                 INQUIRE(file=fname,exist=there)
                 if (there) then
                     open(15,file=fname,status='old')
                     read(15,*) field
                     close(15)
                 else
                     print *,'Cannot find assimilate file: ',fname
!                     call abort_parallel()
                     stop
                 endif
                 call field2alloc(nvar,naux,.true.)
                 if (mxlevel /=1) then
                     do ii=mxlevel-1,1,-1!!!!!!!!!!!!!!!!!!!!!!!!!!!
                       call update(ii,nvar,naux)
                     enddo
                 endif
                 call valout(1,mxlevel,time_geoclaw,nvar,naux)
             endif output
#endif
!             print *,'before get state'
             call MPI_barrier(mpi_comm_world,mpierr)
!             print *, "I am process", mype_world
!             call MPI_barrier(mpi_comm_world,mpierr)
             call PDAF_get_state(steps,time_pdaf,doexit,next_observation_pdaf,&
                 distribute_state_pdaf,prepoststep_ens_pdaf,status_pdaf)
             print *,'(1,1) after get state is',field(25311)
!             print *, "I am process", mype_world
!             call MPI_barrier(mpi_comm_world,mpierr)
             stepnow_pdaf=0
!             call field2alloc(nvar,naux)
             call field2alloc(nvar,naux,.false.)!not for output purpose
             

             if (mxlevel /=1) then
                 do ii=mxlevel -1,1!!!!!!!!!!!!
                   call update(ii,nvar,naux)
                 enddo
             endif
        Endif

     else ! not in assimilation time
         status_pdaf=0
     endif
  ! Check for errors during execution of PDAF

!  call field2alloc(nvar,naux)
  IF (status_pdaf /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in PDAF_put_state - stopping! (PE ', mype_world,')'
     CALL  abort_parallel()
  END IF

END SUBROUTINE assimilate_pdaf

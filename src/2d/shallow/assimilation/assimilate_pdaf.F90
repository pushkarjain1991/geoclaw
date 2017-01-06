!$Id: assimilate_pdaf.F90 1415 2013-09-25 14:33:26Z lnerger $
!BOP
!
! !ROUTINE: assimilate_pdaf - Routine to control perform analysis step
!
! !INTERFACE:
SUBROUTINE assimilate_pdaf()

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
  use mod_parallel, only: mype_world,abort_parallel
  USE mod_assimilation, &      ! Variables for assimilation
       ONLY: filtertype

  IMPLICIT NONE

! !CALLING SEQUENCE:
! Called by: step
! CAlls: PDAF_assimilate_X
!EOP

! Local variables
  INTEGER :: status_pdaf       ! PDAF status flag


  ! External subroutines
  EXTERNAL :: collect_state_pdaf, & ! Routine to collect a state vector from model fields
       init_dim_obs_pdaf, &         ! Initialize Dimension Of Observation Vector
       obs_op_pdaf, &               ! Implementation of the Observation operator
       init_obs_pdaf, &             ! Routine to provide vector of measurements
       prepoststep_ens_pdaf, &      ! User supplied pre/poststep routine
       prodRinvA_pdaf, &            ! Provide product R^-1 A for some matrix A
       init_obsvar_pdaf, &          ! Initialize mean observation error variance
       next_observation_pdaf, &     ! Provide time step, model time, &
                                    ! and dimension of next observation
       distribute_state_pdaf, &     ! Routine to distribute a state vector to model fields
       add_obs_error_pdaf, &        ! Add obs. error covariance R to HPH in EnKF
       init_obscovar_pdaf           ! Initialize obs error covar R in EnKF
  EXTERNAL :: init_n_domains_pdaf, &   ! Provide number of local analysis domains
       init_dim_l_pdaf, &              ! Initialize state dimension for local ana. domain
       init_dim_obs_l_pdaf,&           ! Initialize dim. of obs. vector for local ana. domain
       g2l_state_pdaf, &               ! Get state on local ana. domain from global state
       l2g_state_pdaf, &               ! Init global state from state on local analysis domain
       g2l_obs_pdaf, &                 ! Restrict a global obs. vector to local analysis domain
       init_obs_l_pdaf, &              ! Provide vector of measurements for local ana. domain
       prodRinvA_l_pdaf, &             ! Provide product R^-1 A for some local matrix A
       init_obsvar_l_pdaf, &           ! Initialize local mean observation error variance
       init_obs_f_pdaf, &              ! Provide full vector of measurements for PE-local domain
       obs_op_f_pdaf, &                ! Obs. operator for full obs. vector for PE-local domain
       init_dim_obs_f_pdaf             ! Get dimension of full obs. vector for PE-local domain


! *********************************
! *** Call assimilation routine ***
! *********************************
      if(mype_world==0) write(*,'(a,5x,a)') 'PDAF','Perform assimilation with PDAF '
      if (filtertype==0) then
        CALL PDAF_assimilate_seek(collect_state_pdaf,&
                 distribute_state_pdaf, init_dim_obs_pdaf,&
                 obs_op_pdaf, init_obs_pdaf,&
                 prepoststep_ens_pdaf, prodRinvA_pdaf,&
                 next_observation_pdaf, status_pdaf)
      else if (filtertype==1) then
        CALL PDAF_assimilate_seik(collect_state_pdaf,&
                 distribute_state_pdaf, init_dim_obs_pdaf,&
                 obs_op_pdaf, init_obs_pdaf,&
                 prepoststep_ens_pdaf, prodRinvA_pdaf,&
                 init_obscovar_pdaf,next_observation_pdaf,&
                 status_pdaf)
      else if (filtertype==2) then
        CALL PDAF_assimilate_enkf(collect_state_pdaf,&
                 distribute_state_pdaf, init_dim_obs_pdaf,&
                 obs_op_pdaf, init_obs_pdaf,&
                 prepoststep_ens_pdaf, add_obs_error_pdaf,&
                 init_obscovar_pdaf,next_observation_pdaf,&
                 status_pdaf)
        else if (filtertype==5) then
        CALL PDAF_assimilate_letkf(collect_state_pdaf,&
                 distribute_state_pdaf, init_dim_obs_f_pdaf,&
                 obs_op_f_pdaf, init_obs_f_pdaf,&
                 init_obs_l_pdaf, prepoststep_ens_pdaf, &
                 prodRinvA_l_pdaf,init_n_domains_pdaf,&
                 init_dim_l_pdaf, init_dim_obs_l_pdaf, &
                 g2l_state_pdaf, l2g_state_pdaf, &
                 g2l_obs_pdaf, init_obsvar_pdaf, &
                 init_obsvar_l_pdaf, next_observation_pdaf,&
                 status_pdaf)
      else if (filtertype==4) then
        CALL PDAF_assimilate_etkf(collect_state_pdaf,&
                 distribute_state_pdaf, init_dim_obs_pdaf,&
                 obs_op_pdaf, init_obs_pdaf,&
                 prepoststep_ens_pdaf, prodRinvA_pdaf,&
                 init_obscovar_pdaf,next_observation_pdaf,&
                 status_pdaf)
      else if (filtertype==6) then
        CALL PDAF_assimilate_estkf(collect_state_pdaf,&
                 distribute_state_pdaf, init_dim_obs_pdaf,&
                 obs_op_pdaf, init_obs_pdaf,&
                 prepoststep_ens_pdaf, prodRinvA_pdaf,&
                 init_obscovar_pdaf,next_observation_pdaf,&
                 status_pdaf)
      !else if (filtertype==8) then
      !  CALL PDAF_assimilate_netf(collect_state_pdaf,&
      !           distribute_state_pdaf, init_dim_obs_pdaf,&
      !           obs_op_pdaf, init_obs_pdaf,&
      !           prepoststep_ens_pdaf, prodRinvA_pdaf,&
      !           init_obscovar_pdaf,next_observation_pdaf,&
      !           status_pdaf)
        else if (filtertype==7) then
        CALL PDAF_assimilate_lestkf(collect_state_pdaf,&
                 distribute_state_pdaf, init_dim_obs_f_pdaf,&
                 obs_op_f_pdaf, init_obs_f_pdaf,&
                 init_obs_l_pdaf, prepoststep_ens_pdaf, &
                 prodRinvA_l_pdaf,init_n_domains_pdaf,&
                 init_dim_l_pdaf, init_dim_obs_l_pdaf, &
                 g2l_state_pdaf, l2g_state_pdaf,g2l_obs_pdaf,&
                 init_obsvar_pdaf,init_obsvar_l_pdaf, &
                 next_observation_pdaf, status_pdaf)
      else
        if (mype_world==0) print *,"invalid filter type"
        CALL  abort_parallel()
      endif

      if (status_pdaf /= 0) then
        WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
        'ERROR ', status_pdaf, &
            ' in PDAF_put_state - stopping! (PE ', mype_world,')'
        CALL  abort_parallel()
      END IF

!    else ! not in assimilation time
!        status_pdaf=0
!    endif
END SUBROUTINE assimilate_pdaf

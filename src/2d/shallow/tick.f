c
c  -------------------------------------------------------------
c
      subroutine tick(nvar,cut,nstart,vtime,time,naux,start_time,
     &                rest,dt_max)
c
      use geoclaw_module
      use sortarr
      use mapdomain
      use refinement_module, only: varRefTime
      use amr_module
      use topo_module, only: dt_max_dtopo, num_dtopo, topo_finalized,
     &                       aux_finalized, topo0work
      use gauges_module, only: setbestsrc, num_gauges
      use gauges_module, only: print_gauges_and_reset_nextLoc

      use storm_module, only: landfall, display_landfall_time

#ifdef USE_PDAF
      use mod_model, only: field
      use mod_assimilation, only: dim_ens, filename, filtertype
#endif

      implicit double precision (a-h,o-z)

      logical vtime,dumpout/.false./,dumpchk/.false./,rest,dump_final
      dimension dtnew(maxlv), ntogo(maxlv), tlevel(maxlv)
      integer clock_start, clock_finish, clock_rate
 
#ifdef USE_PDAF

      EXTERNAL :: prodRinvA_pdaf,
     &             next_observation_pdaf, 
     & ! Provide time step, model time and dimension of next observation
     &            distribute_state_pdaf, 
     & ! Routine to distribute a state vector to model fields
     &             prepoststep_ens_pdaf,
     & ! Routine to collect a state vector from model fields
     &             collect_state_pdaf,
     & ! Initialize dimension of observation vector
     &             init_dim_obs_pdaf, 
     & ! Implementation of observation operator
     &             obs_op_pdaf, 
     & ! Routine to provide vector of measurements
     &             init_obs_pdaf,
     & ! Add obs. error covariance R to HPH in EnKF
     &             add_obs_error_pdaf, 
     & ! Initialize obs error covar R in EnKF
     &             init_obscovar_pdaf
      
      integer :: index_2d_pdaf(2)

      integer,PARAMETER :: nx1 = 100
      integer, PARAMETER :: ny1 = 100
      integer i_pkj
      integer j_pkj
      integer dim_counter ! Variable to count the ensemble number. Currently it
                          ! hardcoded to 3. In future has to be automated
      !integer cycle_counter
      integer :: nsteps_pdaf
      integer :: doexit
      integer :: status_pdaf
      real(kind=8) :: h_pkj, hu_pkj, hv_pkj, eta_pkj 
      CHARACTER(len=3) :: ncycle_str
      CHARACTER(len=21) :: assim_filestr
      LOGICAL :: there
      !INTEGER, PARAMETER :: numofens = 5
      integer :: input_counter
      integer :: i1,i2,i_mod,j_mod
      integer, allocatable :: mptr_array(:)
      integer, allocatable :: ordered_mptr_array(:)
      integer :: Ntot
      real(kind=8),allocatable :: corner_array(:,:)
      real(kind=8), allocatable :: sample_array(:,:)
      integer fieldval, corner_counter

                real(kind=8), allocatable :: total_array(:)
                integer, allocatable :: global_to_totalarray_map(:)
                integer, allocatable :: ptr_location(:) !Pointer

      counter_func(i,j) = j + (i-1)*100
!      neighbors(i,j) = total_array(global_to_totalarray_map(counter_func
!     .(i,j)))
      neighbors(i,j) = global_to_totalarray_map(j + 
     . (i-1)*100)
      neighbors2(i,j) = total_array(neighbors(i,j))
      iadd(ivar,i,j)  = loc + ivar - 1 + nvar*((j-1)* mitot+i-1)
     
      iaddaux(iaux,i,j) = locaux + iaux-1 + naux*(i-1) +
     .                                      naux*mitot*(j-1)

      !REAL :: timenow
      !integer nsteps1 = 73
      type ensstate
          integer :: ens_number
          !REAL(KIND=8), DIMENSION(2, nx1*ny1) :: mom = 0.0d0
          REAL(KIND=8), ALLOCATABLE :: mom(:,:)
      end type ensstate
      type(ensstate) :: ens_num(dim_ens)
#endif
#ifndef USE_PDAF
      integer,PARAMETER :: nx1 = 50
      integer, PARAMETER :: ny1 = 50
      integer :: i1,i2,i_mod,j_mod
      integer, allocatable :: mptr_array(:)
      integer, allocatable :: ordered_mptr_array(:)
      integer :: Ntot
      real(kind=8),allocatable :: corner_array(:,:)
      integer fieldval, corner_counter
      real(kind=8) :: ximc,xim,x,xip,xipc,yjmc,yjm,y,yjp,yjpc,dq
      REAL(KIND=8) :: field(nx1*ny1) = 0.0d0
      iadd(ivar,i,j)  = loc + ivar - 1 + nvar*((j-1)* mitot+i-1)
      iaddaux(iaux,i,j) = locaux + iaux-1 + naux*(i-1) +
     .                                      naux*mitot*(j-1)
      fielddata(i,j) = loc+nvar*((j-1)*mitot+i-1)
#endif
      character(len=128) :: time_format
      real(kind=8) cpu_start,cpu_finish

c
c :::::::::::::::::::::::::::: TICK :::::::::::::::::::::::::::::
c  main driver routine.  controls:
c        integration  of all grids.
c        error estimation / regridding
c        output counting
c        updating of fine to coarse grids

c  parameters:
c     nstop   = # of coarse grid time steps to be taken
c     iout    = output interval every 'iout' coarse time steps
c               (if 0, not used - set to inf.)
c     vtime   = true for variable timestep, calculated each coarse step
c
c  integration strategy is to advance a fine grid until it catches
c  up to the coarse grid. this strategy is applied recursively.
c  coarse grid goes first.
c
c  nsteps: used to count how number steps left for a level to be
c          integrated before it catches up with the next coarser level.
c  ncycle: counts number of coarse grid steps = # cycles.
c
c  icheck: counts the number of steps (incrementing by 1
c          each step) to keep track of when that level should
c          have its error estimated and finer levels should be regridded.
c ::::::::::::::::::::::::::::::::::::;::::::::::::::::::::::::::
c
      ncycle         = nstart
      call setbestsrc()     ! need at very start of run, including restart
      if (iout .eq. 0) then
c        # output_style 1 or 2
         iout  = iinfinity
         nextout = 0
         if (nout .gt. 0) then
            nextout = 1
            if (nstart .gt. 0) then
c              # restart: make sure output times start after restart time
               do ii = 1, nout
                 if (tout(ii) .gt. time) then
                   nextout = ii
                   go to 2
                 endif
               end do
  2         continue
            endif
         endif
      endif


      nextchk = 1
      if ((nstart .gt. 0) .and. (abs(checkpt_style).eq.2)) then
c        if this is a restart, make sure chkpt times start after restart time
         do ii = 1, nchkpt
           if (tchk(ii) .gt. time) then
              nextchk = ii
              go to 3
              endif
           enddo
  3      continue
         endif

      tlevel(1)      = time

      do 5 i       = 2, mxnest
       tlevel(i) = tlevel(1)
 5     continue

#ifndef USE_PDAF
         print *,"numgrids = ",numgrids(1)
         !-------------------
         !Get the order of mptr 
         !-------------------
         corner_counter = numgrids(1)
         print *,"num of corners",corner_counter

!         allocate(mptr_array(corner_counter))
!         allocate(ordered_mptr_array(corner_counter))
!         allocate(corner_array(corner_counter,2))
!         
!         mptr = lstart(1)
!         do i_pkj =1,corner_counter
!             mptr_array(i_pkj) = mptr    
!             corner_array(i_pkj,1) = rnode(cornxlo,mptr)
!             corner_array(i_pkj,2) = rnode(cornylo,mptr)
!             mptr = node(levelptr, mptr)
!         enddo
!         print *,"mptr_array = ", mptr_array
!      print *,"corner array = ",(corner_array(i1,:),i1=1,corner_counter)
!         
!      call set_global(corner_array,mptr_array,ordered_mptr_array)
       call get_ordered_array(mptr_array,ordered_mptr_array)


          OPEN(24,file="../ens_1.txt", STATUS="old")
          do i2 = 1,size(ordered_mptr_array)
                  mptr1 = ordered_mptr_array(i2)
                  nx = node(ndihi,mptr1) - node(ndilo, mptr1) + 1
                  ny = node(ndjhi,mptr1) - node(ndjlo, mptr1) + 1
                  Ntot = nx*ny
                  DO i_pkj = 1,Ntot
                      READ(24,*) field((i2-1)*Ntot+i_pkj)
                  !    print *,i_pkj
                  ENDDO
          ENDDO
          CLOSE(24)
          print *,"Read ens_1.tx"
!         level = 1
!         lend = lfine
!         ngrids=0
! 65      if (level .gt. lend) go to 91

         mptr = lstart(1)
! 71      if (mptr .eq. 0) go to 91
          do i2 = 1,size(ordered_mptr_array)
                  mptr1 = ordered_mptr_array(i2)
                  nx = node(ndihi,mptr1) - node(ndilo, mptr1) + 1
                  ny = node(ndjhi,mptr1) - node(ndjlo, mptr1) + 1
                  mitot   = nx + 2*nghost
                  mjtot   = ny + 2*nghost
                  corn1 = rnode(cornxlo,mptr1)
                  corn2 = rnode(cornylo,mptr1)
                  loc     = node(store1, mptr1)
                  locaux  = node(storeaux,mptr1)
                  Ntot = nx*ny
                  
                  print *,"loc=",loc
                  print *,"mptr1=",mptr1
                  print *, "nx = ",nx
                  print *, "ny = ",ny
                  print *, "corners = ",corn1,corn2


                !From left to right. Bsically from 1 to nx
                  !CLOSE(24)

                  do j_pkj = nghost+1, mjtot-nghost
                      do i_pkj = nghost+1, mitot-nghost
                          do ivar=1,nvar
                  if (abs(alloc(iadd(ivar,i_pkj,j_pkj))) < 1d-90) then
                                  alloc(iadd(ivar,i_pkj,j_pkj)) = 0.d0
                              endif
                          enddo
                        ! Extract depth and momenta
                        i_mod = i_pkj-nghost
                        j_mod = j_pkj - nghost
                        !print *,i_mod,j_mod
!                        alloc(iadd(1,i_pkj,j_pkj)) = 
!     .                            field(i_mod+nx*(j_mod-1)) 
!     .                    - alloc(iaddaux(1,i_pkj, j_pkj))
                        alloc(iadd(1,i_pkj,j_pkj)) =
     .              field((i2-1)*Ntot+i_mod+nx*(j_mod-1))
     .                    - alloc(iaddaux(1,i_pkj, j_pkj))
                        alloc(iadd(2,i_pkj, j_pkj)) = 0.d0
                        alloc(iadd(3,i_pkj, j_pkj)) = 0.d0
                  !   endif
                        
                  !      if (alloc(iadd(1,i_pkj,j_pkj)) < 0.d0) then
                  !              alloc(iadd(1,i_pkj,j_pkj)) = 0.d0
                  !      endif
                  
                  !        do ivar=1,nvar
                  !if (abs(alloc(iadd(ivar,i_pkj,j_pkj))) < 1d-90) then
                  !                alloc(iadd(ivar,i_pkj,j_pkj)) = 0.d0
                  !            endif
                  !        enddo
                        !h_pkj = alloc(iadd(1,i_pkj,j_pkj)) 
                        !hu_pkj = alloc(iadd(2,i_pkj,j_pkj))
                        !hv_pkj = alloc(iadd(3,i_pkj,j_pkj))
                        !eta_pkj = h_pkj + alloc(iaddaux(1,i_pkj,j_pkj))
                        !print *,field(i_pkj, j_pkj)
                        
                        !if (abs(eta_pkj) < 1d-90) then
                        !   eta_pkj = 0.d0
                        !end if

                        !print *,h_pkj, hu_pkj, hv_pkj, eta_pkj
               !         print *,alloc(iadd(1,i_pkj,j_pkj)), 
     .         !                 alloc(iadd(2,i_pkj,j_pkj)),
     .         !                 alloc(iadd(3,i_pkj,j_pkj)),  
     .         !                 alloc(iaddaux(1,i_pkj,j_pkj))
                     enddo
                  enddo
             !mptr = node(levelptr, mptr)
             enddo
!             if (mptr .ne. 0) go to 71
!             go to 71
! 82         level = level +1
!             go to 65
! 91          continue 
#endif

#ifdef USE_PDAF
      dim_counter = 0

      DO i_pkj=1,dim_ens
          ALLOCATE(ens_num(i_pkj)%mom(2,numcells(1)))
          ens_num(i_pkj)%mom(:,:) = 0.0d0
      enddo

      print *, "Number of ensemble members is ",dim_ens
      
      !Start of the model loop. This loop runs from 1:num_ens
      pdaf_modelloop: DO

         ens_num(dim_counter + 1)%ens_number = dim_counter + 1
         print *,"Executing forecast of ens_number", 
     &    ens_num(dim_counter + 1)%ens_number

         ! Gets the current time (timenow) and the number of steps (nsteps_pdaf) to forecast.
         ! next_observation - Check next_observation_pdaf.F
         ! distribute_state_pdaf - Copy wse from state_p to field
         call PDAF_get_state(nsteps_pdaf, timenow, doexit,
     &         next_observation_pdaf,distribute_state_pdaf, 
     &         prepoststep_ens_pdaf, status_pdaf)

         !--------------------------
         !Printing the received PDAF state
         !print *,"Just got pdaf_get_state"
         !do i_pkj = 1,50
         !    do  j_pkj = 1,50
         !        print *,field(i_pkj, j_pkj)
         !    enddo
         !    print *,''
         !enddo
         !--------------------------

         print *,"I am being executed"
         !print *,doexit, status_pdaf

          !Check if forecast has to be performed
          !if (time.lt.obstime .and. time + possk(1) .ge. obstime) then
          checkforecast: if (doexit /=1 .and. status_pdaf == 0) then
              
              ! ***Forecast ensemble state ***
              forecast_ensemble: if (nsteps_pdaf>0) then 
                  
                  !Initialize current model time 
                  print *,"Time1 = ",time
                  time = timenow
                  !ncycle = nsteps_pdaf
                  ncycle = time/possk(1)
                  tlevel(1) = time !same as 20 lines above. PDAF loop can start there
                  do 6 i       = 2, mxnest
                      tlevel(i) = tlevel(1)
 6                continue
                   
                  print *,"Time2 = ",time
                  print *,"ncycle = ",ncycle
                  !cycle_counter = 0
        



         !-------------------
         !Get the order of mptr 
         !-------------------
         corner_counter = numgrids(1)


         print *,"num of corners",corner_counter

!         allocate(mptr_array(corner_counter))
!         allocate(ordered_mptr_array(corner_counter))
!         allocate(corner_array(corner_counter,2))
!         
!         mptr = lstart(1)
!         do i_pkj =1,corner_counter
!             mptr_array(i_pkj) = mptr    
!             corner_array(i_pkj,1) = rnode(cornxlo,mptr)
!             corner_array(i_pkj,2) = rnode(cornylo,mptr)
!             mptr = node(levelptr, mptr)
!         enddo
!         print *,"mptr_array = ", mptr_array
!      print *,"corner array = ",(corner_array(i1,:),i1=1,corner_counter)
!         
!      call set_global(corner_array,mptr_array,ordered_mptr_array)
      call get_ordered_array(mptr_array,ordered_mptr_array)
      
!      allocate(sample_array(100,100))
!      !open(unit=23,file="../inputs_online/masked_topo.txt")
!      open(unit=23,file="../yoyo.txt")
!      do i_pkj=1,100
!          read(23,*) sample_array(i_pkj,:)
!      !print *,i_pkj
!     end do
!      close(23)

!      call traverse_global(sample_array, ordered_mptr_array
!     &        , total_array,global_to_totalarray_map, ptr_location)
!
!      print *,"mik1 = ", total_array(global_to_totalarray_map(63 + 
!     . (51-1)*100))
!      print *,"mik2 = ", total_array(neighbors(51,63))
!      print *, "mik3 = ", sample_array(51,63)
!      print *,"2d to 1d is = ", neighbors(14,100), neighbors(59,45)
!     ., neighbors(65,20), neighbors(100,100)
!      
!      print *,"yolo"
!      call oned_to_twod(3200,2,2,index_2d_pdaf)
!      print *, "2d domain is", index_2d_pdaf
!
!      deallocate(sample_array)
      

         mptr = lstart(1)
! 71      if (mptr .eq. 0) go to 91
          do i2 = 1,size(ordered_mptr_array)
                  mptr1 = ordered_mptr_array(i2)
                  nx = node(ndihi,mptr1) - node(ndilo, mptr1) + 1
                  ny = node(ndjhi,mptr1) - node(ndjlo, mptr1) + 1
                  mitot   = nx + 2*nghost
                  mjtot   = ny + 2*nghost
                  corn1 = rnode(cornxlo,mptr1)
                  corn2 = rnode(cornylo,mptr1)
                  loc     = node(store1, mptr1)
                  locaux  = node(storeaux,mptr1)
                  
                  print *,"loc=",loc
                  print *,"mptr1 1=",mptr1
                  print *, "nx = ",nx
                  print *, "ny = ",ny
                  print *, "corners = ",corn1,corn2

                  Ntot = nx*ny
                  do j_pkj = nghost+1, mjtot-nghost
                      do i_pkj = nghost+1, mitot-nghost
                          do ivar=1,nvar
                  if (abs(alloc(iadd(ivar,i_pkj,j_pkj))) < 1d-90) then
                                  alloc(iadd(ivar,i_pkj,j_pkj)) = 0.d0
                              endif
                          enddo
                        ! Extract depth and momenta
                        i_mod = i_pkj-nghost
                        j_mod = j_pkj - nghost
!                        alloc(iadd(1,i_pkj,j_pkj)) = 
!     .                            field(i_mod+nx*(j_mod-1)) 
!     .                    - alloc(iaddaux(1,i_pkj, j_pkj))
                        alloc(iadd(1,i_pkj,j_pkj)) =
     .              field((i2-1)*Ntot+i_mod+nx*(j_mod-1))
     .                    - alloc(iaddaux(1,i_pkj, j_pkj))
                        alloc(iadd(2,i_pkj, j_pkj)) = 
     .   ens_num(dim_counter + 1)%mom(1, (i2-1)*Ntot+i_mod+nx*(j_mod-1))
                    alloc(iadd(3,i_pkj, j_pkj)) = 
     .   ens_num(dim_counter + 1)%mom(2, (i2-1)*Ntot+i_mod+nx*(j_mod-1))
                     enddo
                  enddo
             enddo
         deallocate(mptr_array)
         deallocate(ordered_mptr_array)
!         deallocate(corner_array)


                  !Advance pdaf_nsteps of forward model 
                  !CALL integrate(time, nsteps_pdaf)
#endif
c
c  ------ start of coarse grid integration loop. ------------------
!
 20   if (ncycle .ge. nstop .or. time .ge. tfinal) goto 999

      if (nout .gt. 0) then
          if (nextout  .le. nout) then
             outtime       = tout(nextout)
          else
             outtime       = rinfinity
          endif
      else
          outtime = tfinal
      endif

      if (nextchk  .le. nchkpt) then
         chktime       = tchk(nextchk)
      else
         chktime       = rinfinity
      endif

      dumpout = .false.  !# may be reset below


      if (time.lt.outtime .and. time+1.001*possk(1) .ge. outtime) then
c        ## adjust time step  to hit outtime exactly, and make output
c        #  apr 2010 mjb: modified so allow slightly larger timestep to
c        #  hit output time exactly, instead of taking minuscule timestep
c        #  should still be stable since increase dt in only 3rd digit.
         oldposs = possk(1)
         possk(1) = outtime - time
c        write(*,*)" old possk is ", possk(1)
         diffdt = oldposs - possk(1)  ! if positive new step is smaller


         if (.false.) then  
            write(*,122) diffdt,outtime  ! notify of change
 122        format(" Adjusting timestep by ",e10.3,
     .             " to hit output time of ",e12.6)
c           write(*,*)" new possk is ", possk(1)
            if (diffdt .lt. 0.) then ! new step is slightly larger
              pctIncrease = -100.*diffdt/oldposs   ! minus sign to make whole expr. positive
              write(*,123) pctIncrease
 123          format(" New step is ",e8.2," % larger.",
     .               "  Should still be stable")
              endif
            endif


         do i = 2, mxnest
            possk(i) = possk(i-1) / kratio(i-1)
            enddo
         if (nout .gt. 0) then
            nextout = nextout + 1
            dumpout = .true.
            endif
      endif


      if (time.lt.chktime .and. time + possk(1) .ge. chktime) then
c        ## adjust time step  to hit chktime exactly, and do checkpointing
         possk(1) = chktime - time
         do 13 i = 2, mxnest
 13         possk(i) = possk(i-1) / kratio(i-1)
         nextchk = nextchk + 1
        dumpchk = .true.
      else
        dumpchk = .false.
      endif

c
      level        = 1
      ntogo(level) = 1
      dtnew(1:maxlv) = rinfinity
C       do i = 1, maxlv
C          dtnew(i)  = rinfinity
C       enddo

c     We should take at least one step on all levels after any
c     moving topography (dtopo) has been finalized to insure that
c     all aux arrays are consistent with the final topography.
c     The variable aux_finalized is incremented so that we can check
c     if this is true by checking if aux_finalized == 2 elsewhere in code.

      if (aux_finalized .eq. 1 .and. num_dtopo > 0) then
c         # this is only true once, and only if there was moving topo
          deallocate(topo0work)
          endif 
      if (topo_finalized .and. (aux_finalized .lt. 2)) then
          aux_finalized = aux_finalized + 1
          endif

c
c     ------------- regridding  time?  ---------
c
c check if either
c   (i)  this level should have its error estimated before being advanced
c   (ii) this level needs to provide boundary values for either of
c        next 2 finer levels to have their error estimated.
c        this only affects two grid levels higher, occurs because
c        previous time step needs boundary vals for giant step.
c  no error estimation on finest possible grid level
c
 60       continue
          if (icheck(level) .ge. kcheck) then
               lbase = level
          else if (level+1 .ge. mxnest) then
               go to 90
          else if (icheck(level+1) .ge. kcheck) then
               lbase = level+1
          else if (level+2 .ge. mxnest) then
               go to 90
          else if (icheck(level+2) .ge. kcheck) then
               lbase = level+2
          else
               go to 90
          endif
          if (lbase .eq. mxnest .or. lbase .gt. lfine) go to 70
c
c regrid level 'lbase+1' up to finest level.
c level 'lbase' stays fixed.
c
          if (rprint) write(outunit,101) lbase
101       format(8h  level ,i5,32h  stays fixed during regridding )

          call system_clock(clock_start,clock_rate)
          call cpu_time(cpu_start)
          call regrid(nvar,lbase,cut,naux,start_time)
          call system_clock(clock_finish,clock_rate)
          call cpu_time(cpu_finish)
          timeRegridding = timeRegridding + clock_finish - clock_start
          timeRegriddingCPU=timeRegriddingCPU+cpu_finish-cpu_start

!#ifdef USE_PDAF
!          if (dim_counter == 0) then
!              print *,"hello123"
!              call setbestsrc()     ! need at every grid change
!          endif
!#endif

#ifndef USE_PDAF
          call setbestsrc()     ! need at every grid change
#endif
c         call conck(1,nvar,naux,time,rest)
c         call outtre(lstart(lbase+1),.true.,nvar,naux)
c note negative time to signal regridding output in plots
c         call valout(lbase,lfine,-tlevel(lbase),nvar,naux)
c
c  maybe finest level in existence has changed. reset counters.
c
!#ifdef USE_PDAF
!          if(dim_counter == 0) then
!          if (rprint .and. lbase .lt. lfine) then
!             call outtre(lstart(lbase+1),.false.,nvar,naux)
!          endif
!          endif
!#endif

#ifndef USE_PDAF
          if (rprint .and. lbase .lt. lfine) then
             call outtre(lstart(lbase+1),.false.,nvar,naux)
          endif
#endif

70       continue
          do 80  i  = lbase, lfine
 80          icheck(i) = 0
          do 81  i  = lbase+1, lfine
 81          tlevel(i) = tlevel(lbase)
c
c          MJB: modified to check level where new grids start, which is lbase+1
          if (verbosity_regrid.ge.lbase+1) then
                 do levnew = lbase+1,lfine
                     write(6,1006) intratx(levnew-1),intraty(levnew-1),
     &                             kratio(levnew-1),levnew
 1006                format('   Refinement ratios...  in x:', i3, 
     &                 '  in y:',i3,'  in t:',i3,' for level ',i4)
                 end do

              endif

c  ------- done regridding --------------------
c
c integrate all grids at level 'level'.
c
 90       continue


          call advanc(level,nvar,dtlevnew,vtime,naux)

c Output time info
          timenew = tlevel(level)+possk(level)
          time_format = "(' AMRCLAW: level ',i2,'  CFL = ',e8.3," //
     &                  "'  dt = ',e10.4,  '  final t = ',e12.6)"
          if (display_landfall_time) then
            timenew = (timenew - landfall) / (3.6d3 * 24d0)
            time_format = "(' AMRCLAW: level ',i2,'  CFL = ',e8.3," //
     &                  "'  dt = ',e10.4,  '  final t = ', f5.2)"
          end if
          if (tprint) then
              write(outunit, time_format) level, cfl_level, 
     &                                    possk(level), timenew
          endif
          if (method(4).ge.level) then
              print time_format, level, cfl_level, possk(level), timenew
          endif

c        # to debug individual grid updates...
c        call valout(level,level,time,nvar,naux)
c
c done with a level of integration. update counts, decide who next.
c
          ntogo(level)  = ntogo(level) - 1
          dtnew(level)  = dmin1(dtnew(level),dtlevnew)
          tlevel(level) = tlevel(level) + possk(level)
          icheck(level) = icheck(level) + 1
c
          if (level .lt. lfine) then
             level = level + 1
c            #  check if should adjust finer grid time step to start wtih
             if (((possk(level-1) - dtnew(level-1))/dtnew(level-1)) .gt.
     .            .05) then
                dttemp = dtnew(level-1)/kratio(level-1)
                ntogo(level) = (tlevel(level-1)-tlevel(level))/dttemp+.9
              else
                ntogo(level) = kratio(level-1)
              endif
             possk(level) = possk(level-1)/ntogo(level)
             go to 60
          endif
c
 105      if (level .eq. 1) go to 110
              if (ntogo(level) .gt. 0) then
c                same level goes again. check for ok time step
 106             if ((possk(level)-dtnew(level))/dtnew(level)
     .                .gt. .05)  then

                    write(6,601) level, time
 601                format(" ***adjusting timestep for level ", i3,
     &                     " at t = ",d16.6)
                    print *,"    old ntogo dt",ntogo(level),possk(level)

c                   adjust time steps for this and finer levels
                    ntogo(level) = ntogo(level) + 1
                    possk(level) = (tlevel(level-1)-tlevel(level))/
     .                             ntogo(level)
                    if (varRefTime) then
                       kratio(level-1) = ceiling(possk(level-1) /
     .                                           possk(level))
                    endif
                    print *,"    new ntogo dt ",ntogo(level),
     &                      possk(level)
                    go to 106
                 endif
                 if (ntogo(level) .gt. 100) then
                     write(6,*) "**** Too many dt reductions ****"
                     write(6,*) "**** Stopping calculation   ****"
                     write(6,*) "**** ntogo = ",ntogo(level)
                     write(6,1006) intratx(level-1),intraty(level-1),
     &                             kratio(level-1),level
                     write(6,*) "Writing checkpoint file at t = ",time
                     call check(ncycle,time,nvar,naux)
                     if (num_gauges .gt. 0) then
                        do ii = 1, num_gauges
                           call print_gauges_and_reset_nextLoc(ii, nvar)
                        end do
                     endif
                     stop
                 endif

                 go to 60
              else
                 level = level - 1
                 call system_clock(clock_start,clock_rate)
                 call update(level,nvar,naux)
                 call system_clock(clock_finish,clock_rate)
                 timeUpdating=timeUpdating+clock_finish-clock_start
              endif
          go to 105
c
c      time for output?  done with the whole thing?
c
 110      continue
          

21        time    = time   + possk(1)
          ncycle  = ncycle + 1
          call conck(1,nvar,naux,time,rest)

#ifndef USE_PDAF

      if ( .not.vtime) goto 201

        ! Adjust time steps if variable time step and/or variable
        ! refinement ratios in time
        if (.not. varRefTime) then
          ! find new dt for next cycle (passed back from integration routine).
           do 115 i = 2, lfine
             ii = lfine+1-i
             dtnew(ii) = min(dtnew(ii),dtnew(ii+1)*kratio(ii))
 115       continue
           possk(1) = dtnew(1)
           do 120 i = 2, mxnest
 120         possk(i) = possk(i-1) / kratio(i-1)
        else  ! since refinement ratio in time can change need to set new timesteps in different order
c             ! use same alg. as when setting refinement when first make new fine grids
          dtnew(1) = min(dtnew(1),dt_max)
          if ((num_dtopo>0).and.(topo_finalized.eqv..false.)) then
              dtnew(1) = min(dtnew(1),dt_max_dtopo)
          endif

          possk(1) = dtnew(1)
          do 125 i = 2, lfine
             if (dtnew(i)  .gt. possk(i-1)) then
               kratio(i-1) = 1  ! cant have larger timestep than parent level
               possk(i)    = possk(i-1)
            else
               kratio(i-1) = ceiling(possk(i-1)/dtnew(i))  ! round up for stable integer ratio
               possk(i)    = possk(i-1)/kratio(i-1)        ! set exact timestep on this level
           endif
 125    continue


      endif

 201  if ((abs(checkpt_style).eq.3 .and. 
     &      mod(ncycle,checkpt_interval).eq.0) .or. dumpchk) then
                call check(ncycle,time,nvar,naux)
                dumpchk = .true.
               if (num_gauges .gt. 0) then
                  do ii = 1, num_gauges
                     call print_gauges_and_reset_nextLoc(ii, nvar)
                  end do
               endif
       endif

       if ((mod(ncycle,iout).eq.0) .or. dumpout) then
         call valout(1,lfine,time,nvar,naux)
         if (printout) call outtre(mstart,.true.,nvar,naux)
         if (num_gauges .gt. 0) then
            do ii = 1, num_gauges
               call print_gauges_and_reset_nextLoc(ii, nvar)
            end do
         endif
       endif
#endif

#ifdef USE_PDAF
          !cycle_counter = cycle_counter+1
          !print *,cycle_counter
          !if (cycle_counter==nsteps_pdaf) then
          if (mod(ncycle,nsteps_pdaf) == 0) then
              goto 22
          endif
#endif

      go to 20

#ifdef USE_PDAF
!22            continue
22              end if forecast_ensemble
                  
                  ! Calculating eta after forecast step which causes
                  ! alloc(iadd(1,i,j)) to modify.
                  ! eta = h + topo

                  print *,"Done forecasting for ens_num", dim_counter+1




         !-------------------
         !Get the order of mptr 
         !-------------------
!         mptr = lstart(1)
         !if (mptr .eq. 0) go to 89
!         corner_counter = 0
!94         mptr = node(levelptr, mptr)
!         corner_counter = corner_counter + 1
!         if (mptr .ne. 0) go to 94
!         print *,"num of corners",corner_counter

!         allocate(mptr_array(corner_counter))
!         allocate(ordered_mptr_array(corner_counter))
!         allocate(corner_array(corner_counter,2))
!         
!         mptr = lstart(1)
!         do i_pkj =1,corner_counter
!             mptr_array(i_pkj) = mptr    
!             corner_array(i_pkj,1) = rnode(cornxlo,mptr)
!             corner_array(i_pkj,2) = rnode(cornylo,mptr)
!             mptr = node(levelptr, mptr)
!         enddo
!         print *,"mptr_array = ", mptr_array
!      print *,"corner array = ",(corner_array(i1,:),i1=1,corner_counter)
         
!      call set_global(corner_array,mptr_array,ordered_mptr_array)
      call get_ordered_array(mptr_array,ordered_mptr_array)

         mptr = lstart(1)
! 71      if (mptr .eq. 0) go to 91
          do i2 = 1,size(ordered_mptr_array)
                  mptr1 = ordered_mptr_array(i2)
                  nx = node(ndihi,mptr1) - node(ndilo, mptr1) + 1
                  ny = node(ndjhi,mptr1) - node(ndjlo, mptr1) + 1
                  mitot   = nx + 2*nghost
                  mjtot   = ny + 2*nghost
                  corn1 = rnode(cornxlo,mptr1)
                  corn2 = rnode(cornylo,mptr1)
                  loc     = node(store1, mptr1)
                  locaux  = node(storeaux,mptr1)
                  Ntot = nx*ny 
                  print *,"loc=",loc
                  print *,"mptr1 2=",mptr1
                  print *, "nx = ",nx
                  print *, "ny = ",ny
                  print *, "corners = ",corn1,corn2


                  do j_pkj = nghost+1, mjtot-nghost
                      do i_pkj = nghost+1, mitot-nghost
                          do ivar=1,nvar
                  if (abs(alloc(iadd(ivar,i_pkj,j_pkj))) < 1d-90) then
                                  alloc(iadd(ivar,i_pkj,j_pkj)) = 0.d0
                              endif
                          enddo
                        ! Extract depth and momenta
                        i_mod = i_pkj-nghost
                        j_mod = j_pkj - nghost
                        field((i2-1)*Ntot+i_mod+nx*(j_mod-1)) = 
     .      alloc(iadd(1,i_pkj,j_pkj)) + alloc(iaddaux(1,i_pkj, j_pkj))
         ens_num(dim_counter + 1)%mom(1,(i2-1)*Ntot+i_mod+nx*(j_mod-1))=
     .                      alloc(iadd(2,i_pkj, j_pkj))

         ens_num(dim_counter + 1)%mom(2,(i2-1)*Ntot+i_mod+nx*(j_mod-1))=
     .                  alloc(iadd(3,i_pkj, j_pkj)) 
                     enddo
                  enddo
             enddo
         deallocate(mptr_array)
         deallocate(ordered_mptr_array)
         !deallocate(corner_array)




              ! *** PDAF: Send State forecast to filter;
              ! *** PDAF: Perform assimilation if ensemble forecast is completed
              ! field is the eta here. Check collect_state
              if (filtertype==2) then
                  CALL PDAF_put_state_enkf(collect_state_pdaf, 
     &             init_dim_obs_pdaf, obs_op_pdaf, init_obs_pdaf, 
     &             prepoststep_ens_pdaf, add_obs_error_pdaf, 
     &             init_obscovar_pdaf, status_pdaf)
              else if (filtertype==1) then
                  CALL PDAF_put_state_seik(collect_state_pdaf, 
     &             init_dim_obs_pdaf, obs_op_pdaf, init_obs_pdaf, 
     &             prepoststep_ens_pdaf, prodRinvA_pdaf, 
     &             init_obsvar_pdaf, status_pdaf)
              else if (filtertype==0) then
                  CALL PDAF_put_state_seek(collect_state_pdaf, 
     &             init_dim_obs_pdaf, obs_op_pdaf, init_obs_pdaf, 
     &             prepoststep_ens_pdaf, prodRinvA_pdaf, status_pdaf)
              else if (filtertype==6) then
                  CALL PDAF_put_state_estkf(collect_state_pdaf, 
     &             init_dim_obs_pdaf, obs_op_pdaf, init_obs_pdaf, 
     &             prepoststep_ens_pdaf, prodRinvA_pdaf, 
     &             init_obsvar_pdaf, status_pdaf)
              else
                      print *,"nothing"
              endif
              
              dim_counter = dim_counter + 1

              !init_obsvar_pdaf, status_pdaf)

              ! After assimilation step, alloc(iadd()) must contain assimilated
              ! value for appropriate valout
              !WRITE(ncycle_str,'(i2.2)') ncycle
              WRITE(ncycle_str,'(i3)') ncycle
              print *,'dim_counter = ',dim_counter
         
             !-----------------------------------
              ! This is just for valout purpose     
             !-----------------------------------
              dimcounter: if (dim_counter == dim_ens) then
                          !dim_counter = 0
                   !Read Analysis file and overwrite iadd(1,i,j)
                WRITE(assim_filestr,'(A21)') 'state_step'//
     &        TRIM(ADJUSTL(ncycle_str))//'_ana.txt'
            assim_filestr = TRIM(ADJUSTL(assim_filestr))
            INQUIRE(FILE=assim_filestr, EXIST=there)
            IF (THERE) THEN
                PRINT *,assim_filestr,' exists'
                OPEN(20, file = assim_filestr, status = 'old')
                READ(20,*) field
                
                !DO i_pkj = 1,ny1
                    !READ(20,*) field(:,i_pkj)
                !ENDDO
                        
                CLOSE(20)
            ELSE 
                PRINT *,"Cannot find assimilated file", assim_filestr
                EXIT PDAF_MODELLOOP
            ENDIF
            !print *,'field read from PDAF output = ', field 
                
         !-------------------
         !Get the order of mptr 
         !-------------------
!         mptr = lstart(1)
         !if (mptr .eq. 0) go to 89
!         corner_counter = 0
!95         mptr = node(levelptr, mptr)
!         corner_counter = corner_counter + 1
!         if (mptr .ne. 0) go to 95
!         print *,"num of corners",corner_counter

!         allocate(mptr_array(corner_counter))
!         allocate(ordered_mptr_array(corner_counter))
!         allocate(corner_array(corner_counter,2))
         
!         mptr = lstart(1)
!         do i_pkj =1,corner_counter
!             mptr_array(i_pkj) = mptr    
!             corner_array(i_pkj,1) = rnode(cornxlo,mptr)
!             corner_array(i_pkj,2) = rnode(cornylo,mptr)
!             mptr = node(levelptr, mptr)
!         enddo
!         print *,"mptr_array = ", mptr_array
!      print *,"corner array = ",(corner_array(i1,:),i1=1,corner_counter)
         
!      call set_global(corner_array,mptr_array,ordered_mptr_array)
       call get_ordered_array(mptr_array,ordered_mptr_array)

         mptr = lstart(1)
          do i2 = 1,size(ordered_mptr_array)
                  mptr1 = ordered_mptr_array(i2)
                  nx = node(ndihi,mptr1) - node(ndilo, mptr1) + 1
                  ny = node(ndjhi,mptr1) - node(ndjlo, mptr1) + 1
                  mitot   = nx + 2*nghost
                  mjtot   = ny + 2*nghost
                  corn1 = rnode(cornxlo,mptr1)
                  corn2 = rnode(cornylo,mptr1)
                  loc     = node(store1, mptr1)
                  locaux  = node(storeaux,mptr1)
                  Ntot = nx*ny 
                  do j_pkj = nghost+1, mjtot-nghost
                      do i_pkj = nghost+1, mitot-nghost
                        i_mod = i_pkj-nghost
                        j_mod = j_pkj - nghost
                          !READ(20,*) field(i_pkj, j_pkj)
                          !print *, field(i_pkj-nghost, j_pkj-nghost)
                        alloc(iadd(1,i_pkj,j_pkj)) = 
     .                            field((i2-1)*Ntot+i_mod+nx*(j_mod-1)) 
     .                    - alloc(iaddaux(1,i_pkj, j_pkj))
!                        alloc(iadd(1,i_pkj,j_pkj)) = 
!     .                            field(i_mod+nx*(j_mod-1)) 
!     .                    - alloc(iaddaux(1,i_pkj, j_pkj))
!!                          alloc(iadd(1,i_pkj,j_pkj)) = 
!!     .                            field(i_pkj-nghost,j_pkj-nghost) 
!!     .                            - alloc(iaddaux(1,i_pkj, j_pkj))
                          !print *,alloc(iadd(1,i_pkj,j_pkj))
                      enddo
                  !print *,''
                  enddo
            enddo
         deallocate(mptr_array)
         deallocate(ordered_mptr_array)
!         deallocate(corner_array)
                          
                   
                   if ( .not.vtime) goto 202

                       ! Adjust time steps if variable time step and/or variable
                       ! refinement ratios in time
                       if (.not. varRefTime) then
                           ! find new dt for next cycle (passed back from integration routine).
                           do 115 i = 2, lfine
                               ii = lfine+1-i
                               dtnew(ii) = min(dtnew(ii),dtnew(ii+1)*
     &                         kratio(ii))
 115                       continue
                           possk(1) = dtnew(1)
                           do 120 i = 2, mxnest
 120                           possk(i) = possk(i-1) / kratio(i-1)
                       else  ! since refinement ratio in time can change need to set new timesteps in different order
c                          ! use same alg. as when setting refinement when first make new fine grids
                           dtnew(1) = min(dtnew(1),dt_max)
                           if ((num_dtopo>0).and.(topo_finalized.eqv.
     &                     .false.)) then
                               dtnew(1) = min(dtnew(1),dt_max_dtopo)
                           endif

                           possk(1) = dtnew(1)
                           do 125 i = 2, lfine
                               if (dtnew(i)  .gt. possk(i-1)) then
                                   kratio(i-1) = 1  ! cant have larger timestep than parent level
                                   possk(i)    = possk(i-1)
                               else
                                   kratio(i-1) = ceiling(possk(i-1)/
     &                             dtnew(i))
                                   ! round up for stable integer ratio
                                   possk(i)    = possk(i-1)/kratio(i-1) 
                                   ! set exact timestep on this level
                               endif
 125                       continue


                       endif

 202              if ((checkpt_style.eq.3 .and. 
     &                 mod(ncycle,checkpt_interval).eq.0) 
     &                 .or. dumpchk) then
                       call check(ncycle,time,nvar,naux)
                       dumpchk = .true.
                   endif

                   if ((mod(ncycle,iout).eq.0) .or. dumpout) then
                       call valout(1,lfine,time,nvar,naux)
                      if (printout) call outtre(mstart,.true.,nvar,naux)
                   endif

              endif dimcounter
              
              if (dim_counter == dim_ens) then
                  dim_counter = 0
              endif
              !call MPI_barrier(COMM_model, MPIERR)
              !call finalize_pdaf()

          else checkforecast
              !No more assimilation work; exit loop
              EXIT pdaf_modelloop
      
          end if checkforecast
          

      ENDDO pdaf_modelloop
      
      !call finalize_pdaf()

#endif
c
999   continue
              

c
c  # computation is complete to final time or requested number of steps
c
       if (ncycle .ge. nstop .and. tfinal .lt. rinfinity) then
c         # warn the user that calculation finished prematurely
          write(outunit,102) nstop
          write(6,102) nstop
  102     format('*** Computation halted after nv(1) = ',i8,
     &           '  steps on coarse grid')
          endif
c
c  # final output (unless we just did it above)
c
      dump_final = ((iout.lt.iinfinity) .and. (mod(ncycle,iout).ne.0))
      if (.not. dumpout) then
          if (nout > 0) then
              dump_final = (tout(nout).eq.tfinal)
              endif
          endif
      
      if (dump_final) then
           call valout(1,lfine,time,nvar,naux)
           if (printout) call outtre(mstart,.true.,nvar,naux)
           if (num_gauges .gt. 0) then
              do ii = 1, num_gauges
                 call print_gauges_and_reset_nextLoc(ii, nvar)
              end do
           endif
      endif

c  # checkpoint everything for possible future restart
c  # (unless we just did it based on dumpchk)
c

c
      if (checkpt_style .ne. 0) then  ! want a chckpt
         ! check if just did it so dont do it twice
         if (.not. dumpchk) call check(ncycle,time,nvar,naux)
      endif
      if (num_gauges .gt. 0) then
         do ii = 1, num_gauges
            call print_gauges_and_reset_nextLoc(ii, nvar)
         end do
      endif

      write(6,*) "Done integrating to time ",time
      return
      end

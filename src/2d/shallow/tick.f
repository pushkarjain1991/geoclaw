c
c  -------------------------------------------------------------
c
      subroutine tick(nvar,cut,nstart,vtime,time,naux,start_time,
     &                rest,dt_max)
c
      use geoclaw_module
#ifdef USE_PDAF
      use mod_parallel, only: mype_world, mpierr, mpi_comm_world
      use mod_assimilation, only: stepnow_pdaf, assimilate_step,
     & regrid_assim, second_valout, assimilation_time, 
     & global_coordinate
      use mod_model, only: field
      use common_level, only: regions_assim
#endif
#ifdef CHILE_GEN_ENS
      use mod_parallel_ens_gen
#endif
      use refinement_module, only: varRefTime
      use amr_module
      use topo_module, only: dt_max_dtopo, num_dtopo, topo_finalized,
     &                       aux_finalized, topo0work
      use gauges_module, only: setbestsrc, num_gauges
      use gauges_module, only: print_gauges_and_reset_nextLoc

      use storm_module, only: landfall, display_landfall_time


      implicit double precision (a-h,o-z)

      logical vtime,dumpout/.false./,dumpchk/.false./,rest,dump_final
      dimension dtnew(maxlv), ntogo(maxlv), tlevel(maxlv)
      integer clock_start, clock_finish, clock_rate
      integer tick_clock_start, tick_clock_finish, tick_clock_rate
      character(len=128) :: time_format
      real(kind=8) cpu_start,cpu_finish
#ifdef USE_PDAF
      character(len=3) :: stepstr1
      character(len=2) :: ensstr2
      logical :: dir_exists
      real(kind=8), allocatable :: temp_field(:)
      integer :: curr_tot_step
      integer :: num_ex = 1
      integer :: field_size
      integer :: i
#endif

#ifdef CHILE_GEN_ENS
      character(len=2) :: ensstr2
      logical :: dir_exists
      integer :: i
#endif

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
      call system_clock(tick_clock_start,tick_clock_rate)
      call cpu_time(tick_cpu_start)


      ncycle         = nstart
#ifdef USE_PDAF
      if(mype_world == 0) then
          call setbestsrc()     ! need at very start of run, including restart
      endif
#else
      call setbestsrc()     ! need at very start of run, including restart
#endif

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

c
c  ------ start of coarse grid integration loop. ------------------
c
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
     .             " to hit output time of ",e13.6)
c           write(*,*)" new possk is ", possk(1)
            if (diffdt .lt. 0.) then ! new step is slightly larger
              pctIncrease = -100.*diffdt/oldposs   ! minus sign to make whole expr. positive
              write(*,123) pctIncrease
 123          format(" New step is ",e9.2," % larger.",
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
#ifdef USE_PDAF
              print *, mype_world, "Regridding non-assimilation.... 
     &       lbase = ", lbase
          regrid_assim = .false.
#endif
          call regrid(nvar,lbase,cut,naux,start_time)
#ifdef USE_PDAF
          !call mpi_barrier(mpi_comm_world, mpierr)
          call print_num_cells(nvar, naux)
          !call mpi_barrier(mpi_comm_world, mpierr)
#endif
          call system_clock(clock_finish,clock_rate)
          call cpu_time(cpu_finish)
          timeRegridding = timeRegridding + clock_finish - clock_start
          timeRegriddingCPU=timeRegriddingCPU+cpu_finish-cpu_start

#ifdef USE_PDAF
      if(mype_world == 0) then
          call setbestsrc()     ! need at very start of run, including restart
      endif
#else
      call setbestsrc()     ! need at every grid change
#endif

c         call conck(1,nvar,naux,time,rest)
c         call outtre(lstart(lbase+1),.true.,nvar,naux)
c note negative time to signal regridding output in plots
c         call valout(lbase,lfine,-tlevel(lbase),nvar,naux)
c
c  maybe finest level in existence has changed. reset counters.
c
          if (rprint .and. lbase .lt. lfine) then
             call outtre(lstart(lbase+1),.false.,nvar,naux)
          endif
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
#ifdef USE_PDAF
          time_format = "(' Ens_num ',i1, ' AMRCLAW: level ',i2,' 
     & CFL = ',e8.3," //"'  dt = ',e10.4,  '  final t = ',e12.6)"
#elif CHILE_GEN_ENS
          time_format = "(' Ens_num ',i1, ' AMRCLAW: level ',i2,' 
     & CFL = ',e8.3," //"'  dt = ',e10.4,  '  final t = ',e12.6)"
#else
          time_format = "(' AMRCLAW: level ',i2,'  CFL = ',e8.3," //
     &                  "'  dt = ',e10.4,  '  final t = ',e12.6)"
#endif
          if (display_landfall_time) then
            timenew = (timenew - landfall) / (3.6d3 * 24d0)
            time_format = "(' AMRCLAW: level ',i2,'  CFL = ',e8.3," //
     &                  "'  dt = ',e10.4,  '  final t = ', f5.2)"
          end if
          if (tprint) then
#ifdef USE_PDAF
              write(outunit, time_format) mype_world, level, cfl_level, 
     &                                    possk(level), timenew
#elif CHILE_GEN_ENS
              write(outunit, time_format) mype_world, level, cfl_level, 
     &                                    possk(level), timenew
#else
              write(outunit, time_format) level, cfl_level, 
     &                                    possk(level), timenew
#endif
          endif
          if (method(4).ge.level) then
#ifdef USE_PDAF
              print time_format, mype_world, level, cfl_level, 
     & possk(level), timenew
#elif CHILE_GEN_ENS
              print time_format, mype_world, level, cfl_level, 
     & possk(level), timenew
#else
              print time_format, level, cfl_level, possk(level), timenew
#endif
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
                           call print_gauges_and_reset_nextLoc(ii)
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
c  --------------one complete coarse grid integration cycle done. -----
c
c      time for output?  done with the whole thing?
c
 110      continue
          time    = time   + possk(1)
          ncycle  = ncycle + 1
          call conck(1,nvar,naux,time,rest)

          print *, "timetime ", time
#ifdef USE_PDAF
          assimilation_time = outtime
          stepnow_pdaf = stepnow_pdaf + 1
          
          !Ready for assimilation
          if (abs(time - assimilation_time) < 1.0D-06) then
              !Adhoc thingy. Should be replaced
              stepnow_pdaf = assimilate_step
          
              print *, "Assimilating data at time = ", assimilation_time
              regrid_assim=.true.

              ! Perform regridding with regrid_assim True
              ! This is to get the ensembles obtain same flagging 
              ! and hence same patches
              !call mpi_barrier(mpi_comm_world, mpierr)
              print *, "Regridding first time ya", mype_world 
              call print_num_cells(nvar, naux)
              print *, "yo123"
              do i = 1, mxnest-1
                call extract_regions()
                call regrid(nvar,1,0.7,naux,start_time)
                if(mype_world == 0) then
                    call setbestsrc()     ! need at every grid change
                endif
              enddo
              call print_num_cells(nvar, naux)
              call check_ensemble_grid()
              print *, "yo1234"
              
          
          if (rprint .and. lbase .lt. lfine) then
             call outtre(lstart(lbase+1),.false.,nvar,naux)
          endif

              if (mxnest <= 3) then
                if (lfine /=1) then
                  do ii=lfine-1,1,-1
                      call update(1, nvar, naux)
                  enddo
                  print *, "Update done"
                endif
              endif

              ! Put the geoclaw alloc values to PDAF field
              ! Next step is to use the field for assimilation
              !call mpi_barrier(mpi_comm_world, mpierr)
              call alloc2field(nvar,naux)
              print *, "field size = ", size(field)
              !call mpi_barrier(mpi_comm_world, mpierr)
              
              ! Update dim_state_p and state
              !call mpi_barrier(mpi_comm_world, mpierr)
              call update_dim_state_p()
              !call mpi_barrier(mpi_comm_world, mpierr)

              field_size = size(field)
              !Set global coordinate - For localization only
              if(allocated(global_coordinate)) then
                deallocate(global_coordinate)
              endif
              allocate(global_coordinate(2,field_size))
       call set_global_coordinate_array(field_size, global_coordinate)

          
          ! Perform assimilation
          ! 1. alloc2field
          ! 2. Put state
          ! 3. Update
          ! 4. field2alloc
              print *, "time =" , time
              call assimilate_pdaf(time)
              call mpi_barrier(mpi_comm_world, mpierr)

            
              ! Put the assimilated values from field to alloc
              call field2alloc(nvar,naux)!not for output purpose
              print *, "Running field2alloc after assimilation"
              !call mpi_barrier(mpi_comm_world, mpierr)
              deallocate(field) 

              !if (mxnest <= 3) then
               !if (lfine /=1) then
               !   do ii=lfine -1, 1, -1
                      call put2zero(nvar,naux)
                      call update(1,nvar,naux)
               !   enddo
               ! endif
              !endif
              !call mpi_barrier(mpi_comm_world, mpierr)
              
          endif
          
#endif
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

#ifndef USE_PDAF
 201  if ((abs(checkpt_style).eq.3 .and. 
     &      mod(ncycle,checkpt_interval).eq.0) .or. dumpchk) then
                call check(ncycle,time,nvar,naux)
                dumpchk = .true.
               if (num_gauges .gt. 0) then
                  do ii = 1, num_gauges
                     call print_gauges_and_reset_nextLoc(ii)
                  end do
               endif
       endif
#endif

#ifdef USE_PDAF

       ! Output state_ana
201    if (mype_world == 0) then
           if ((mod(ncycle,iout).eq.0) .or. dumpout) then
               
               ! Copy original alloc to field
               call alloc2field(nvar, naux)
               call update_dim_state_p(nvar, naux)
               field_size = size(field)

               ! Copy original field to temporary field
               if(allocated(temp_field)) deallocate(temp_field)
               allocate(temp_field(field_size))
               temp_field(:) = field(:)

               ! Read state_step_ana
               curr_tot_step = num_ex * assimilate_step 
               num_ex = num_ex + 1
               write(stepstr1, '(i3.1)') curr_tot_step
               OPEN(20, file='state_step'//TRIM(ADJUSTL(stepstr1))//
     &         '_ana.txt', status = 'old')
               
               ! Assign state_step_ana to field
               print *, "field has size = ", field_size
               if(allocated(field)) deallocate(field)
               allocate(field(field_size))
               read(20,"(e26.16)") field(:)
               close(20)

               ! Perform field to alloc
               call field2alloc(nvar, naux)
               
               ! Perform valout
               print *, "valout state_ana "
               call valout(1,lfine,time,nvar,naux)
               ! Perform gauge output
               if (printout) call outtre(mstart,.true.,nvar,naux)
               if (num_gauges .gt. 0) then
                   do ii = 1, num_gauges
                       call print_gauges_and_reset_nextLoc(ii)
                   end do
               endif

               ! Assign temporary field to field
               if(allocated(field)) deallocate(field)
               allocate(field(field_size))
               field(:) = temp_field(:)

               ! Assign field to alloc
               call field2alloc(nvar, naux)

               ! Deallocate temp_field
               deallocate(temp_field)
               deallocate(field)

           endif
       endif
       !call mpi_barrier(mpi_comm_world, mpierr)
          
!       ! Forecast output for every ensemble  
!       if ((mod(ncycle,iout).eq.0) .or. dumpout) then
!           !call mpi_barrier(mpi_comm_world, mpierr)
!           print *, "Forecast valout for ens ", mype_world
!           !call mpi_barrier(mpi_comm_world, mpierr)
!           
!           write(ensstr2, '(i2.1)') mype_world
!           inquire(file="_output_"//trim(adjustl(ensstr2))//"_for", 
!     &     exist=dir_exists)
!           if(dir_exists .eqv. .false.) then
!               call system('mkdir _output_'// 
!     &     trim(adjustl(ensstr2))//"_for")
!           endif
!           call chdir("_output_"//trim(adjustl(ensstr2))//"_for")
!           second_valout = .true.
!           call valout(1,lfine,time,nvar,naux)
!           second_valout = .false.
!           call chdir("../")
!      endif


       ! Analysis output for every ensemble  
       if ((mod(ncycle,iout).eq.0) .or. dumpout) then
           !call mpi_barrier(mpi_comm_world, mpierr)
           print *, "Analysis valout for ens ", mype_world
           !call mpi_barrier(mpi_comm_world, mpierr)
           
           write(ensstr2, '(i2.1)') mype_world
           inquire(file="_output_"//trim(adjustl(ensstr2))//"_ana", 
     &     exist=dir_exists)
           if(dir_exists .eqv. .false.) then
               call system('mkdir _output_'// 
     &         trim(adjustl(ensstr2))//"_ana")
           endif
           call chdir("_output_"//trim(adjustl(ensstr2))//"_ana")
           second_valout = .true.
           call valout(1,lfine,time,nvar,naux)
           second_valout = .false.
           call chdir("../")
       endif
          
#elif CHILE_GEN_ENS
       ! Forecast output for every ensemble  
       if ((mod(ncycle,iout).eq.0) .or. dumpout) then
           !call mpi_barrier(mpi_comm_world, mpierr)
           print *, "Analysis valout for ens ", mype_world
           !call mpi_barrier(mpi_comm_world, mpierr)
           
           write(ensstr2, '(i2.1)') mype_world
           inquire(file="_output_"//trim(adjustl(ensstr2))//"_for", 
     &     exist=dir_exists)
           if(dir_exists .eqv. .false.) then
               call system('mkdir _output_'// 
     &     trim(adjustl(ensstr2))//"_for")
           endif
           call chdir("_output_"//trim(adjustl(ensstr2))//"_for")
           second_valout = .true.
           call valout(1,lfine,time,nvar,naux)
           second_valout = .false.
           call chdir("../")
      endif
#endif

#ifndef USE_PDAF ! CHECK THIS PLEASE IF IT IS NEEDED
       if ((mod(ncycle,iout).eq.0) .or. dumpout) then
         call valout(1,lfine,time,nvar,naux)
         if (printout) call outtre(mstart,.true.,nvar,naux)
         if (num_gauges .gt. 0) then
            do ii = 1, num_gauges
               call print_gauges_and_reset_nextLoc(ii)
            end do
         endif
       endif
#endif

#ifdef USE_PDAF
       !VERY VERY IMPORTANT
       !if (stepnow_pdaf == assimilate_step) then
       !if (time == assimilation_time) then
       if (abs(time - assimilation_time) < 1.0D-06) then
           stepnow_pdaf = 0
       endif
       print *, "reached yo2", mype_world, stepnow_pdaf
       !call mpi_barrier(mpi_comm_world, mpierr)
#endif


      go to 20
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
      
#ifndef USE_PDAF
      if (dump_final) then
           call valout(1,lfine,time,nvar,naux)
           if (printout) call outtre(mstart,.true.,nvar,naux)
           if (num_gauges .gt. 0) then
              do ii = 1, num_gauges
                 call print_gauges_and_reset_nextLoc(ii)
              end do
           endif
      endif
#endif

c  # checkpoint everything for possible future restart
c  # (unless we just did it based on dumpchk)
c
      call system_clock(tick_clock_finish,tick_clock_rate)
      call cpu_time(tick_cpu_finish)
      timeTick = timeTick + tick_clock_finish - tick_clock_start 
      timeTickCPU = timeTickCPU + tick_cpu_finish - tick_cpu_start 


c  # checkpoint everything for possible future restart
c  # (unless we just did it based on dumpchk)
c
#ifndef USE_PDAF
      if (checkpt_style .ne. 0) then  ! want a chckpt
         ! check if just did it so dont do it twice
         if (.not. dumpchk) call check(ncycle,time,nvar,naux)
      endif

      if (num_gauges .gt. 0) then
         do ii = 1, num_gauges
            call print_gauges_and_reset_nextLoc(ii)
         end do
      endif
#endif

      write(6,*) "Done integrating to time ",time
      return
      end

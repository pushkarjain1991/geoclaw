!
!  Use adaptive mesh refinement to solve the hyperbolic 2-d equation:
!
!              u  +  f(u)    + g(u)   = 0
!               t         x        y
!
! or the more general non-conservation law form:
!              u  +  A u     + B u    = 0
!               t         x        y
!
!  using the wave propagation method as in CLAWPACK in combination
!  with the locally uniform embedded grids of AMR.
!
!  Estimate error with Richardson extrap. (in errest.f)
!  + gradient checking (in errsp.f).  Initial conditions set
!  in (qinit.f), b.c.'s in (physbd.f).
!
!  Specify rectangular domain from
!           (xlower,ylower) to (xupper,yupper).
!
!
! =========================================================================
!
!  This software is made available for research and instructional use only.
!  You may copy and use this software without charge for these non-commercial
!  purposes, provided that the copyright notice and associated text is
!  reproduced on all copies.  For all other uses (including distribution of
!  modified versions), please contact the author at the address given below.
!
!  *** This software is made available "as is" without any assurance that it
!  *** will work for your purposes.  The software may in fact have defects, so
!  *** use the software at your own risk.
!
!  --------------------------------------
!    AMRCLAW Version 5.0,  2012
!    Homepage: http://www.clawpack.org
!  --------------------------------------
!
!   Authors:
!
!             Marsha J. Berger
!             Courant Institute of Mathematical Sciences
!             New York University
!             251 Mercer St.
!             New York, NY 10012
!             berger@cims.nyu.edu
!
!             Randall J. LeVeque
!             Applied Mathematics
!             Box 352420
!             University of Washington,
!             Seattle, WA 98195-2420
!             rjl@uw.edu
!
! =========================================================================
program amr2

    use amr_module, only: dbugunit, parmunit, outunit, inunit, matunit
    use amr_module, only: mxnest, rinfinity, iinfinity
    use amr_module, only: xupper, yupper, xlower, ylower
    use amr_module, only: hxposs, hyposs, intratx, intraty, kratio
    use amr_module, only: cfl, cflv1, cflmax, evol

    use amr_module, only: checkpt_style, checkpt_interval, tchk, nchkpt
    use amr_module, only: rstfile

    use amr_module, only: max1d, maxvar, maxlv

    use amr_module, only: method, mthlim, use_fwaves, numgrids
    use amr_module, only: nghost, mwaves, mcapa, auxtype
    use amr_module, only: tol, tolsp, flag_richardson, flag_gradient

    use amr_module, only: nghost, mthbc
    use amr_module, only: xperdom, yperdom, spheredom

    use amr_module, only: nstop, nout, iout, tfinal, tout, output_style
    use amr_module, only: output_format, printout, verbosity_regrid
    use amr_module, only: output_q_components, output_aux_components
    use amr_module, only: output_aux_onlyonce, matlabu

    use amr_module, only: lfine, lentot, iregridcount, avenumgrids
    use amr_module, only: tvoll, tvollCPU, rvoll, rvol, mstart, possk, ibuff
    use amr_module, only: timeRegridding,timeUpdating, timeValout
    use amr_module, only: timeBound,timeStepgrid, timeFlagger,timeBufnst,timeFilvalTot
    use amr_module, only: timeBoundCPU,timeStepGridCPU,timeSetauxCPU,timeRegriddingCPU
    use amr_module, only: timeSetaux, timeSetauxCPU, timeValoutCPU
    use amr_module, only: kcheck, iorder, lendim, lenmax

    use amr_module, only: dprint, eprint, edebug, gprint, nprint, pprint
    use amr_module, only: rprint, sprint, tprint, uprint

    use amr_module, only: t0, tstart_thisrun

    ! Data modules
    use geoclaw_module, only: set_geo
#ifdef USE_PDAF_CHILE2
    use topo_module, only: read_topo_settings, read_dtopo_settings, &
    dtopowork
#else
    use topo_module, only: read_topo_settings, read_dtopo_settings
#endif

    use qinit_module, only: set_qinit
    use fixedgrids_module, only: set_fixed_grids
    use refinement_module, only: set_refinement
    use storm_module, only: set_storm
    use friction_module, only: setup_variable_friction
    use gauges_module, only: set_gauges, num_gauges
    use regions_module, only: set_regions, regions, num_regions
    use fgmax_module, only: set_fgmax, FG_num_fgrids
    use multilayer_module, only: set_multilayer

#ifdef USE_PDAF
    use mpi, only: mpi_real8
    use mod_model, only : field, total_steps, reshaped_recv_ic, &
    random_pert
    use mod_parallel, only: MPIerr, mpi_comm_world, mype_world, &
    n_modeltasks
    use mod_assimilation, only: dim_state_p, dim_ens, regrid_assim
    use common_level, only: set_original_regions
    use gauges_module, only: setbestsrc
#endif

    implicit none

    ! Local variables
    integer :: i, iaux, mw, level
    integer :: ndim, nvar, naux, mcapa1, mindim, dimensional_split
    integer :: nstart, nsteps, nv1, nx, ny, lentotsave, num_gauge_SAVE
    integer :: omp_get_max_threads, maxthreads
    real(kind=8) :: time, ratmet, cut, dtinit, dt_max
    logical :: vtime, rest, output_t0    

    ! Timing variables
    integer :: clock_start, clock_finish, clock_rate, ttotal
    real(kind=8) :: cpu_start, cpu_finish, ttotalcpu

    ! Common block variables
    real(kind=8) :: dxmin, dymin

    common /comfine/ dxmin,dymin

    character(len=364) :: format_string
    character(len=*), parameter :: clawfile = 'claw.data'
    character(len=*), parameter :: amrfile = 'amr.data'
    character(len=*), parameter :: outfile = 'fort.amr'
    character(len=*), parameter :: dbugfile = 'fort.debug'
    character(len=*), parameter :: matfile = 'fort.nplot'
    character(len=*), parameter :: parmfile = 'fort.parameters'
#ifdef USE_PDAF
    character(len=3):: ensstr1
    character(len=3):: ensstr3
    character(len=255) :: cwd
    real(kind=8), allocatable :: recv_ic(:)
    integer, parameter :: root = 0
    logical :: dir_exists1
#endif
    
#ifdef USE_PDAF_CHILE
    integer :: n, clock
    real(kind=8) :: pert_mu
    real(kind=8),allocatable :: rand_pert(:)
    integer(kind=4) :: seed
    real(kind=8) :: pert_sigma
    !real(kind=8) :: r8_normal_ab
    !character(len=1) :: ensstr
    real(kind=8) :: recv_random_pert = 0.0
    integer, parameter :: root1 = 0
    logical :: mpierr1
#endif

#ifdef USE_PDAF_CHILE_PERT_RESTART
    real(kind=8), allocatable :: temp_field(:)
    integer :: ii
#endif


#ifdef USE_PDAF
    ! Add parallelization
    CALL init_parallel_pdaf(0,2)
#endif
    ! Open parameter and debug files
    open(dbugunit,file=dbugfile,status='unknown',form='formatted')
    open(parmunit,file=parmfile,status='unknown',form='formatted')

    maxthreads = 1    !! default, if no openmp


    ! Open AMRClaw primary parameter file
    call opendatafile(inunit,clawfile)

    ! Number of space dimensions, not really a parameter but we read it in and
    ! check to make sure everyone is on the same page. 
    read(inunit,"(i1)") ndim  
    if (ndim /= 2) then
        print *,'Error ***   ndim = 2 is required,  ndim = ',ndim
        print *,'*** Are you sure input has been converted'
        print *,'*** to Clawpack 5.x form?'
        stop
    endif
          
    ! Domain variables
    read(inunit,*) xlower, ylower
    read(inunit,*) xupper, yupper
    read(inunit,*) nx, ny
    read(inunit,*) nvar    ! meqn
    read(inunit,*) mwaves
    read(inunit,*) naux
    read(inunit,*) t0

    ! ==========================================================================
    ! Output Options
    ! Output style
    read(inunit,*) output_style
    if (output_style == 1) then
        read(inunit,*) nout
        read(inunit,*) tfinal
        read(inunit,*) output_t0

        iout = 0
    else if (output_style == 2) then
        read(inunit,*) nout
        allocate(tout(nout))
        read(inunit,*) (tout(i), i=1,nout)
        output_t0 = (abs(tout(1) - t0) < 1e-15)
        ! Move output times down one index
        if (output_t0) then
            nout = nout - 1
            do i=1,nout
                tout(i) = tout(i+1)
            enddo
        endif
        iout = 0
        tfinal = tout(nout)
    else if (output_style == 3) then
        read(inunit,*) iout
        read(inunit,*) nstop
        read(inunit,*) output_t0
        nout = 0
        tfinal = rinfinity
    else
        stop "Error ***   Invalid output style."
    endif

    ! Error checking
    if ((output_style == 1) .and. (nout > 0)) then
        allocate(tout(nout))
        do i=1,nout
            tout(i) = t0 + i * (tfinal - t0) / real(nout,kind=8)
        enddo
    endif

    ! What and how to output
    read(inunit,*) output_format
    allocate(output_q_components(nvar))
    read(inunit,*) (output_q_components(i),i=1,nvar)
    if (naux > 0) then
        allocate(output_aux_components(naux))
        read(inunit,*) (output_aux_components(i),i=1,naux)
        read(inunit,*) output_aux_onlyonce
    endif
    ! ==========================================================================

    ! ==========================================================================
    !  Algorithm parameters

    read(inunit,*) possk(1)   ! dt_initial
    read(inunit,*) dt_max     ! largest allowable dt
    read(inunit,*) cflv1      ! cfl_max
    read(inunit,*) cfl        ! clf_desired
    read(inunit,*) nv1        ! steps_max
      
    if (output_style /= 3) then
        nstop = nv1
    endif

    read(inunit,*) vtime      ! dt_variable
    if (vtime) then
        method(1) = 2
    else
        method(1) = 1
    endif

    read(inunit,*) method(2)  ! order
    iorder = method(2)
    read(inunit,*) method(3)  ! order_trans

    read(inunit,*) dimensional_split
    if (dimensional_split > 0) then
        print *, '*** ERROR ***  dimensional_split = ', dimensional_split
        print *, ' dimensional splitting not supported in amrclaw'
        stop
    endif

    read(inunit,*) method(4)   ! verbosity
    read(inunit,*) method(5)   ! src_split
    read(inunit,*) mcapa1
    
    read(inunit,*) use_fwaves
    allocate(mthlim(mwaves))
    read(inunit,*) (mthlim(mw), mw=1,mwaves)

    ! Boundary conditions
    read(inunit,*) nghost
    read(inunit,*) mthbc(1),mthbc(3)
    read(inunit,*) mthbc(2),mthbc(4)

    ! 1 = left, 2 = right 3 = bottom 4 = top boundary
    xperdom = (mthbc(1) == 2 .and. mthbc(2) == 2)
    yperdom =  (mthbc(3) == 2 .and. mthbc(4) == 2)
    spheredom =  (mthbc(3) == 5 .and. mthbc(4) == 5)

    if ((mthbc(1).eq.2 .and. mthbc(2).ne.2) .or. &
        (mthbc(2).eq.2 .and. mthbc(1).ne.2)) then
        
        print *, '*** ERROR ***  periodic boundary conditions: '
        print *, '  mthbc(1) and mthbc(2) must BOTH be set to 2'
        stop
    endif

    if ((mthbc(3).eq.2 .and. mthbc(4).ne.2) .or. &
        (mthbc(4).eq.2 .and. mthbc(3).ne.2)) then

        print *, '*** ERROR ***  periodic boundary conditions: '
        print *, '  mthbc(3) and mthbc(4) must BOTH be set to 2'
        stop
    endif

    if ((mthbc(3).eq.5 .and. mthbc(4).ne.5) .or. &
        (mthbc(4).eq.5 .and. mthbc(3).ne.5)) then
    
        print *, '*** ERROR ***  sphere bcs at top and bottom: '
        print *, '  mthbc(3) and mthbc(4) must BOTH be set to 5'
        stop
    endif

    if (spheredom .and. .not. xperdom) then

        print *,'*** ERROR ***  sphere bcs at top and bottom: '
        print *,'must go with periodic bcs at left and right  '
        stop
    endif

    ! ==========================================================================
    !  Restart and Checkpointing

    read(inunit,*) rest
    read(inunit,*) rstfile

    read(inunit,*) checkpt_style
    if (checkpt_style == 0) then
        ! Never checkpoint:
        checkpt_interval = iinfinity

    else if (abs(checkpt_style) == 2) then
        read(inunit,*) nchkpt
        allocate(tchk(nchkpt))
        read(inunit,*) (tchk(i), i=1,nchkpt)

    else if (abs(checkpt_style) == 3) then
        ! Checkpoint every checkpt_interval steps on coarse grid
        read(inunit,*) checkpt_interval
    endif

    close(inunit)

    ! ==========================================================================
    !  Refinement Control
    call opendatafile(inunit, amrfile)

    read(inunit,*) mxnest
    if (mxnest <= 0) then
        stop 'Error ***   mxnest (amrlevels_max) <= 0 not allowed'
    endif
          
    if (mxnest > maxlv) then
        stop 'Error ***   mxnest > max. allowable levels (maxlv) in common'
    endif
      
    ! Anisotropic refinement always allowed in 5.x:
    read(inunit,*) (intratx(i),i=1,max(1,mxnest-1))
    read(inunit,*) (intraty(i),i=1,max(1,mxnest-1))
    read(inunit,*) (kratio(i), i=1,max(1,mxnest-1))
    read(inunit,*)

    do i=1,mxnest-1
        if ((intratx(i) > max1d) .or. (intraty(i) > max1d)) then 
            print *, ""
            format_string = "(' *** Error: Refinement ratios must be no " // &
                            "larger than max1d = ',i5,/,'     (set max1d" // &
                            " in amr_module.f90)')"
            print format_string, max1d
            stop
        endif
    enddo

    if (naux > 0) then
        allocate(auxtype(naux))
        read(inunit,*) (auxtype(iaux), iaux=1,naux)
    endif
    read(inunit,*)
              
    read(inunit,*) flag_richardson
    read(inunit,*) tol            ! for richardson
    read(inunit,*) flag_gradient
    read(inunit,*) tolsp          ! for gradient
    read(inunit,*) kcheck
    read(inunit,*) ibuff
    read(inunit,*) cut
    read(inunit,*) verbosity_regrid

    ! read verbose/debugging flags
    read(inunit,*) dprint
    read(inunit,*) eprint
    read(inunit,*) edebug
    read(inunit,*) gprint
    read(inunit,*) nprint
    read(inunit,*) pprint
    read(inunit,*) rprint
    read(inunit,*) sprint
    read(inunit,*) tprint
    read(inunit,*) uprint

    close(inunit)
    ! Finished with reading in parameters
    ! ==========================================================================

    ! Look for capacity function via auxtypes:
    mcapa = 0
    do iaux = 1, naux
        if (auxtype(iaux) == "capacity") then
            if (mcapa /= 0) then
                print *, " only 1 capacity allowed"
                stop
            else
                mcapa = iaux
            endif
        endif

        ! Change to Version 4.1 terminology:
        if (auxtype(iaux) == "leftface") auxtype(iaux) = "xleft"
        if (auxtype(iaux) == "bottomface") auxtype(iaux) = "yleft"
        if (.not. (auxtype(iaux) .eq."xleft" .or. &
                   auxtype(iaux) .eq. "yleft".or. &
                   auxtype(iaux) .eq. "capacity".or. &
                   auxtype(iaux) .eq. "center"))  then
            print *," unknown type for auxiliary variables"
            print *," use  xleft/yleft/center/capacity"
            stop
        endif
    enddo

    ! Error checking of input data
    if (mcapa /= mcapa1) then
        stop 'Error ***  mcapa does not agree with auxtype'
    endif
    if (nvar > maxvar) then
        stop 'Error ***   nvar > maxvar in common'
    endif
    if (2*nghost > min(nx,ny) .and. ny /= 1) then
        mindim = 2 * nghost
        print *, 'Error ***   need finer domain >', mindim, ' cells'
        stop
    endif
    if (mcapa > naux) then     
        stop 'Error ***   mcapa > naux in input file'
    endif

    if (.not. vtime .and. nout /= 0) then
        print *,        ' cannot specify output times with fixed dt'
        stop
    endif


    ! Write out parameters
    write(parmunit,*) ' '
    write(parmunit,*) 'Running amrclaw with parameter values:'
    write(parmunit,*) ' '


    print *, ' '
    print *, 'Running amrclaw ...  '
    print *, ' '

    hxposs(1) = (xupper - xlower) / nx
    hyposs(1) = (yupper - ylower) / ny

    ! initialize frame number for output.  
    ! Note: might be reset in restrt if this is a restart
    if (output_t0) then
        matlabu   = 0
    else
        matlabu   = 1
    endif

    if (rest) then

        open(outunit, file=outfile, status='unknown', position='append', &
                      form='formatted')

        call restrt(nsteps,time,nvar,naux)
        nstart  = nsteps
        tstart_thisrun = time
        print *, ' '
        print *, 'Restarting from previous run'
        print *, '   at time = ',time
        print *, ' '


        ! Call user routine to set up problem parameters:
        call setprob()

        ! Non-user defined setup routine
        call set_geo()                    ! sets basic parameters g and coord system
        call set_refinement()             ! sets refinement control parameters
        call read_dtopo_settings()        ! specifies file with dtopo from earthquake
        call read_topo_settings()         ! specifies topography (bathymetry) files
        call set_qinit()                  ! specifies file with dh if this used instead
        call set_fixed_grids()            ! Fixed grid settings
        call setup_variable_friction()    ! Variable friction parameter
        call set_multilayer()             ! Set multilayer SWE parameters
        call set_storm()                  ! Set storm parameters
        call set_regions()                ! Set refinement regions
        call set_gauges(rest, nvar)       ! Set gauge output
        call set_fgmax()

#ifdef USE_PDAF
        regrid_assim = .false.
        call set_original_regions()
#endif

    else
!This is for setting same level at initial time t=0
! regrid_Assim will be switched off after setting initial domain
        open(outunit, file=outfile, status='unknown', form='formatted')

        tstart_thisrun = t0

        ! Call user routine to set up problem parameters:
        call setprob()

        ! Non-user defined setup routine
        call set_geo()                    ! sets basic parameters g and coord system
        call set_refinement()             ! sets refinement control parameters
        call read_dtopo_settings()        ! specifies file with dtopo from earthquake

#ifdef USE_PDAF_CHILE2
        !Initial perturbation in topo itself
        pert_sigma = .005D+00
        pert_mu = 0.0D+00
        if(mype_world == 0) then
            seed = clock
            !seed=123456789

            allocate(rand_pert(n_modeltasks))
            call r8vec_normal_ab(n_modeltasks, pert_mu, pert_sigma, seed, rand_pert)
            !print *, "printing random number", r, mype_world
            call r8vec_print(n_modeltasks, rand_pert, ' Vector of Normal AB values:')
        endif

        call mpi_scatter(rand_pert, 1, mpi_real8, recv_random_pert,1, mpi_real8, &
        root1, mpi_comm_world, mpierr)

        dtopowork(:) = dtopowork(:)*(1 + recv_random_pert) !dtopo is negative
        !dtopowork(:) = dtopowork(:) + recv_random_pert !dtopo is negative
        !print *, dtopowork(:)
#endif
        call read_topo_settings()         ! specifies topography (bathymetry) files
        call set_qinit()                  ! specifies file with dh if this used instead
        call set_fixed_grids()            ! Fixed grid settings
        call setup_variable_friction()    ! Variable friction parameter
        call set_multilayer()             ! Set multilayer SWE parameters
        call set_storm()                  ! Set storm parameters
        call set_regions()                ! Set refinement regions
        call set_gauges(rest, nvar)       ! Set gauge output
!#ifdef USE_PDAF_CHILE
!        write(ensstr3, '(i3.1)') mype_world
!        call getcwd(cwd)
!        print *, cwd
!        inquire(file="_output_"//trim(adjustl(ensstr3))//"_ana",&     
!        exist=dir_exists1)
!        if(dir_exists1 .eqv. .false.) then
!          call system('mkdir _output_'//&
!          trim(adjustl(ensstr3))//'_ana')
!          call chdir('_output_'//&
!          trim(adjustl(ensstr3))//'_ana')
!          call system('cp ../gauges.data .')
!          call set_gauges(rest, nvar)       ! Set gauge output
!          open(12,file='gauge32412.txt',status='new')
!          close(12)
!          call chdir('../')
!        endif
!#endif
        call set_fgmax()

        cflmax = 0.d0   ! otherwise use previously heckpointed val

        lentot = 0
        lenmax = 0
        lendim = 0
        rvol   = 0.0d0
        do i   = 1, mxnest
            rvoll(i) = 0.0d0
        enddo
        evol = 0.0d0
        call stst1()


        ! changed 4/24/09: store dxmin,dymin for setaux before
        ! grids are made, in order to average up from finest grid.
        dxmin = hxposs(mxnest)
        dymin = hyposs(mxnest)

        call domain(nvar,vtime,nx,ny,naux,t0)

        ! Hold off on gauges until grids are set. 
        ! The fake call to advance at the very first timestep 
        ! looks at the gauge array but it is not yet built
        num_gauge_SAVE = num_gauges
        num_gauges = 0
        call setgrd(nvar,cut,naux,dtinit,t0)
        num_gauges = num_gauge_SAVE


! commented out to match 4-x version
!!$        if (possk(1) .gt. dtinit*cflv1/cfl .and. vtime) then
!!$            ! initial time step was too large. reset to dt from setgrd
!!$            print *, "*** Initial time step reset for desired cfl"
!!$            possk(1) = dtinit
!!$            do i = 2, mxnest-1
!!$                possk(i) = possk(i-1)*kratio(i-1)
!!$            end do
!!$        endif

        time = t0
        nstart = 0
#ifdef USE_PDAF
        regrid_assim = .false.
        if(num_regions /= 0) then 
           call set_original_regions()
        endif
#endif
    endif


    write(parmunit,*) ' '
    write(parmunit,*) '--------------------------------------------'
    write(parmunit,*) ' '
    write(parmunit,*) '   rest = ', rest, '   (restart?)'
    write(parmunit,*) '   start time = ',time
    write(parmunit,*) ' '

!$   maxthreads = omp_get_max_threads() 
#ifdef USE_PDAF
     if (mype_world==0) then
         write(outunit,*)" max threads set to ",maxthreads
         print *," max threads set to ",maxthreads
     endif
#else
     write(outunit,*)" max threads set to ",maxthreads
     print *," max threads set to ",maxthreads
#endif
    
    !
    !  print out program parameters for this run
    !
    format_string = "(/' amrclaw parameters:',//," // &
                      "' error tol            ',e12.5,/," // &
                      "' spatial error tol    ',e12.5,/," // &
                      "' order of integrator     ',i9,/," // &
                      "' error checking interval ',i9,/," // &
                      "' buffer zone size        ',i9,/," // &
                      "' nghost                  ',i9,/," // &
                      "' volume ratio cutoff  ',e12.5,/," // &
                      "' max. refinement level   ',i9,/," // &
                      "' user sub. calling times ',i9,/," // &
                      "' cfl # (if var. delt) ',e12.5,/)"
    write(outunit,format_string) tol,tolsp,iorder,kcheck,ibuff,nghost,cut, &
                                 mxnest,checkpt_interval,cfl
    format_string = "(' xupper(upper corner) ',e12.5,/," // &
                     "' yupper(upper corner) ',e12.5,/," // &
                     "' xlower(lower corner) ',e12.5,/," // &
                     "' ylower(lower corner) ',e12.5,/," // &
                     "' nx = no. cells in x dir.',i9,/," // &
                     "' ny = no. cells in y dir.',i9,/," // &
                     "' refinement ratios       ',6i5,/,/)"
    write(outunit,format_string) xupper,yupper,xlower,ylower,nx,ny
    write(outunit,"(' refinement ratios:       ',6i5,/)"  ) &
                                                        (intratx(i),i=1,mxnest)
    write(outunit,"(' refinement ratios:       ',6i5,/)"  ) &
                                                        (intraty(i),i=1,mxnest)
    write(outunit,"(' no. auxiliary vars.     ',i9)") naux
    write(outunit,"('       var ',i5,' of type ', a10)") &
                                                (iaux,auxtype(iaux),iaux=1,naux)
    if (mcapa > 0) write(outunit,"(' capacity fn. is aux. var',i9)") mcapa

#ifdef USE_PDAF
    if (mype_world==0) then
        print *, ' '
        print *, 'Done reading data, starting computation ...  '
        print *, ' '
    endif
#else
        print *, ' '
        print *, 'Done reading data, starting computation ...  '
        print *, ' '
#endif


    call outtre (mstart,printout,nvar,naux)
    write(outunit,*) "  original total mass ..."
    call conck(1,nvar,naux,time,rest)

    ! Timing
    call system_clock(clock_start,clock_rate)
    call cpu_time(cpu_start)

    if (output_t0) then
#ifdef USE_PDAF_CHILE
     if(mype_world == 0) then
        call valout(1,lfine,time,nvar,naux)
    endif
#else
        call valout(1,lfine,time,nvar,naux)
#endif
    endif
    close(parmunit)


#ifdef USE_PDAF_CHILE_PERT_RESTART
        pert_sigma = .01D+00
        pert_mu = 0.0D+00
        if(mype_world == 0) then
            seed = clock
            !seed=123456789

            allocate(rand_pert(n_modeltasks))
            call r8vec_normal_ab(n_modeltasks, pert_mu, pert_sigma, seed, rand_pert)
            !print *, "printing random number", r, mype_world
            call r8vec_print(n_modeltasks, rand_pert, ' Vector of Normal AB values:')
        endif

        call mpi_scatter(rand_pert, 1, mpi_real8, recv_random_pert,1, mpi_real8, &
        root1, mpi_comm_world, mpierr)
        call alloc2field(nvar, naux)
        !temp_field = field + 0.01*mype_world
        temp_field = field + recv_random_pert
        where((field > 0.1 .and. field < 0.3) .or. (field > -0.3 .and. field<-0.1))
            field = temp_field
        endwhere
        !field(:) = temp_field(:)
        call field2alloc(nvar, naux)
        
        !regrid_assim = .true.
        !call extract_regions(t0, regrid_assim)
        !call regrid(nvar,1,cut,naux,t0)
        !call setbestsrc()     ! need at every grid change
        
        !print *, "lfine intiially = ", lfine
        !if (lfine /=1) then
        !  do ii=lfine-1,1
        !    call update(ii, nvar, naux)
        !  enddo
        !  print *, "Update done initially"
        !endif
#endif

#ifdef USE_PDAF
    total_steps = nstop
    
    call alloc2field(nvar, naux)
    !call update_dim_state_p(nvar, naux)
    dim_state_p = size(field)
    
    !call mpi_barrier(mpi_comm_world, mpierr)
    print *, "Initial dim_state_p = ", mype_world, dim_state_p
    !call mpi_barrier(mpi_comm_world, mpierr)
    
    write(ensstr1, '(i3.1)') mype_world
    print *, "writing init_ens# file"
    open(unit=57, file='init_ens'//trim(adjustl(ensstr1)), status='replace')
    do i=1,size(field) 
        write(57,*) field(i)
    enddo
    close(57)
    print *, "Done writing init_ens# file"
    !call mpi_barrier(mpi_comm_world, mpierr)
    
    !Gather the initial conditions to root processor
    print *, "Reading IC, ", mype_world
    print *, "Allocating space for recv_ic"
    if (mype_world == 0) allocate(recv_ic(dim_state_p*n_modeltasks))
    !!PARALLEL_DEBUG
    print *, "Done allocating space for recv_ic"
    !call mpi_barrier(mpi_comm_world, mpierr)
    print *, "Performing initial condition gathering..."
    call mpi_gather(field, dim_state_p, mpi_real8, recv_ic, dim_state_p,&
    mpi_real8, root, mpi_comm_world, mpierr) 
    print *, "Done Performing initial condition gathering..."
    print *, "dim_state_p, n_modeltasks = ", dim_state_p, n_modeltasks
    
    if (mype_world == 0) then
        if(allocated(reshaped_recv_ic)) deallocate(reshaped_recv_ic)
        allocate(reshaped_recv_ic(dim_state_p, n_modeltasks))
        reshaped_recv_ic(:,:) = reshape(recv_ic, (/dim_state_p, n_modeltasks/))
        !print *, "reshaped_recv_ic = ", reshaped_recv_ic(:,:)
        deallocate(recv_ic)
    endif

    if (mype_world==0) then
        WRITE (*, '(1x, a)') 'INITIALIZE GEOCLAW with PDAF'
        WRITE (*, '(10x,a,i4,1x,a1,1x,i4)') 'Grid size:', nx, 'x', ny
        WRITE (*, '(10x,a,i4)') 'Time steps', total_steps
    endif
    CALL init_pdaf()
#endif

    ! --------------------------------------------------------
    !  Tick is the main routine which drives the computation:
    ! --------------------------------------------------------
    
    call tick(nvar,cut,nstart,vtime,time,naux,t0,rest,dt_max)
    ! --------------------------------------------------------

    ! Done with computation, cleanup:

    ! Print out the fgmax files
    if (FG_num_fgrids > 0) call fgmax_finalize()
    
    

    call system_clock(clock_finish,clock_rate)
    call cpu_time(cpu_finish)
    
    !output timing data
    write(*,*)
    write(outunit,*)
    format_string="('============================== Timing Data ==============================')"
    write(outunit,format_string)
    write(*,format_string)
    
    write(*,*)
    write(outunit,*)
    
    !Integration time
    format_string="('Integration Time (stepgrid + BC + overhead)')"
    write(outunit,format_string)
    write(*,format_string)
    
    !Advanc time
    format_string="('Level           Wall Time (seconds)    CPU Time (seconds)   Total Cell Updates')"
    write(outunit,format_string)
    write(*,format_string)
    ttotalcpu=0.d0
    ttotal=0
    do level=1,mxnest
        format_string="(i3,'           ',1f15.3,'        ',1f15.3,'    ', e17.3)"
        write(outunit,format_string) level, &
             real(tvoll(level),kind=8) / real(clock_rate,kind=8), tvollCPU(level), rvoll(level)
        write(*,format_string) level, &
             real(tvoll(level),kind=8) / real(clock_rate,kind=8), tvollCPU(level), rvoll(level)
        ttotalcpu=ttotalcpu+tvollCPU(level)
        ttotal=ttotal+tvoll(level)
    end do
    
    format_string="('total         ',1f15.3,'        ',1f15.3,'    ', e17.3)"
    write(outunit,format_string) &
             real(ttotal,kind=8) / real(clock_rate,kind=8), ttotalCPU, rvol
    write(*,format_string) &
             real(ttotal,kind=8) / real(clock_rate,kind=8), ttotalCPU, rvol
    
    write(*,*)
    write(outunit,*)
    
    
    format_string="('All levels:')"
    write(*,format_string)
    write(outunit,format_string)
    
    
    
    !stepgrid
    format_string="('stepgrid      ',1f15.3,'        ',1f15.3,'    ',e17.3)"
    write(outunit,format_string) &
         real(timeStepgrid,kind=8) / real(clock_rate,kind=8), timeStepgridCPU
    write(*,format_string) &
         real(timeStepgrid,kind=8) / real(clock_rate,kind=8), timeStepgridCPU
    
    !bound
    format_string="('BC/ghost cells',1f15.3,'        ',1f15.3)"
    write(outunit,format_string) &
         real(timeBound,kind=8) / real(clock_rate,kind=8), timeBoundCPU
    write(*,format_string) &
         real(timeBound,kind=8) / real(clock_rate,kind=8), timeBoundCPU
    
    !regridding time
    format_string="('Regridding    ',1f15.3,'        ',1f15.3,'  ')"
    write(outunit,format_string) &
            real(timeRegridding,kind=8) / real(clock_rate,kind=8), timeRegriddingCPU
    write(*,format_string) &
            real(timeRegridding,kind=8) / real(clock_rate,kind=8), timeRegriddingCPU
    
    !output time
    format_string="('Output (valout)',1f14.3,'        ',1f15.3,'  ')"
    write(outunit,format_string) &
            real(timeValout,kind=8) / real(clock_rate,kind=8), timeValoutCPU
    write(*,format_string) &
            real(timeValout,kind=8) / real(clock_rate,kind=8), timeValoutCPU
    
    write(*,*)
    write(outunit,*)
    
    !Total Time
    format_string="('Total time:   ',1f15.3,'        ',1f15.3,'  ')"
    write(outunit,format_string) &
            real(clock_finish - clock_start,kind=8) / real(clock_rate,kind=8), &
            cpu_finish-cpu_start
    write(*,format_string) &
            real(clock_finish - clock_start,kind=8) / real(clock_rate,kind=8), &
            cpu_finish-cpu_start
    
    format_string="('Using',i3,' thread(s)')"
    write(outunit,format_string) maxthreads
    write(*,format_string) maxthreads
    
    
    write(*,*)
    write(outunit,*)
    
    
    write(*,"('Note: The CPU times are summed over all threads.')")
    write(outunit,"('Note: The CPU times are summed over all threads.')")
    write(*,"('      Total time includes more than the subroutines listed above')")
    write(outunit,"('      Total time includes more than the subroutines listed above')")
    
    
    !end of timing data
    write(*,*)
    write(outunit,*)
    format_string="('=========================================================================')"
    write(outunit,format_string)
    write(*,format_string)
    write(*,*)
    write(outunit,*)
    
    
    
    
    !format_string = "('Total time to solution = ',1f16.8,' s, using', i3,' threads')"
    !write(outunit,format_string) &
    !        real(clock_finish - clock_start,kind=8) / real(clock_rate,kind=8), maxthreads
    !write(*,format_string) &
    !        real(clock_finish - clock_start,kind=8) / real(clock_rate,kind=8), maxthreads

    !do level = 1, mxnest            
    !  format_string = "('Total advanc time on level',i3,':',1f16.8,' s')"
    !  write(outunit,format_string) level, &
    !         real(tvoll(level),kind=8) / real(clock_rate,kind=8)
    !  write(*,format_string) level, &
    !         real(tvoll(level),kind=8) / real(clock_rate,kind=8)
    !end do
    !format_string = "('Total updating time:          ',1f16.8,' s')"
    !write(outunit,format_string)  real(timeUpdating,kind=8) / real(clock_rate,kind=8)
    !write(*,format_string) real(timeUpdating,kind=8) / real(clock_rate,kind=8)
    !format_string = "('Total valout time:            ',1f16.8,' s')"
    !write(outunit,format_string)  real(timeValout,kind=8) / real(clock_rate,kind=8)
    !write(*,format_string) real(timeValout,kind=8) / real(clock_rate,kind=8)

    !write(*,*)
    !format_string = "('Total regridding time (clock):',1f16.8,' s')"
    !write(outunit,format_string)  real(timeRegridding,kind=8) / real(clock_rate,kind=8)
    !write(*,format_string) real(timeRegridding,kind=8) / real(clock_rate,kind=8)
 
    !format_string = "('   Total Grdfit     time (clock):    ',1f16.8,' s')"
    !write(outunit,format_string)  real(timeGrdfitAll,kind=8) / real(clock_rate,kind=8)
    !write(*,format_string) real(timeGrdfitAll,kind=8) / real(clock_rate,kind=8)
    !format_string = "('      Total Flglvl  time                ',1f16.8,' s')"
    !write(outunit,format_string)  real(timeFlglvl,kind=8) / real(clock_rate,kind=8)
    !write(*,format_string) real(timeFlglvl,kind=8) / real(clock_rate,kind=8)
    !format_string = "('         Total    Flagger     time (clock)     ',1f16.8,' s')"
    !write(outunit,format_string)  real(timeFlagger,kind=8) / real(clock_rate,kind=8)
    !write(*,format_string) real(timeFlagger,kind=8) / real(clock_rate,kind=8)
    !format_string = "('         Total    Bufnst     time (clock)     ',1f16.8,' s')"
    !write(outunit,format_string)  real(timeBufnst,kind=8) / real(clock_rate,kind=8)
    !write(*,format_string) real(timeBufnst,kind=8) / real(clock_rate,kind=8)

    !format_string = "('   Total gfixup     time  (clock)   ',1f16.8,' s')"
    !write(outunit,format_string)  real(timeGrdfit2,kind=8) / real(clock_rate,kind=8)
    !write(*,format_string) real(timeGrdfit2,kind=8) / real(clock_rate,kind=8)
    !format_string = "('      Total filval (wall) time      ',1f16.8,' s')"
    !write(outunit,format_string)  real(timeFilvalTot,kind=8) / real(clock_rate,kind=8)
    !write(*,format_string) real(timeFilvalTot,kind=8) / real(clock_rate,kind=8)

    !format_string = "('          Total filval (all cores) time     ',1f16.8,' s')"
    !write(outunit,format_string)  real(timeFilval,kind=8) / real(clock_rate,kind=8)
    !write(*,format_string) real(timeFilval,kind=8) / real(clock_rate,kind=8)

    !write(*,*)
    !format_string = "('Total setaux (all cores) time:',1f16.8,' s')"
    !write(outunit,format_string)  real(timeSetaux,kind=8) / real(clock_rate,kind=8)
    !write(*,format_string) real(timeSetaux,kind=8) / real(clock_rate,kind=8)


    !write(*,*)
    !write(*,*)"Integration time, still not counting saveqc"
    !format_string = "('Total Bound (all levels) wall time:   ',1f16.8,' s')"
    !write(outunit,format_string)  real(timeBound,kind=8) / real(clock_rate,kind=8)
    !write(*,format_string) real(timeBound,kind=8) / real(clock_rate,kind=8)
    !format_string = "('Total stepgrid (all levels) wall time:',1f16.8,' s')"
    !write(outunit,format_string)  real(timeStepgrid,kind=8) / real(clock_rate,kind=8)
    !write(*,format_string) real(timeStepgrid,kind=8) / real(clock_rate,kind=8)

    ! Done with computation, cleanup:
    lentotsave = lentot
    call cleanup(nvar,naux)
    if (lentot /= 0) then
        write(outunit,*) lentot," words not accounted for in memory cleanup"
        print *,         lentot," words not accounted for in memory cleanup"
    endif
    
    !
    ! report on statistics
    !
    open(matunit,file=matfile,status='unknown',form='formatted')
    write(matunit,*) matlabu-1
    write(matunit,*) mxnest
    close(matunit)

    write(outunit,*)
    write(outunit,*)
    do i = 1, mxnest
      if (iregridcount(i) > 0) then
        write(outunit,801) i,avenumgrids(i)/iregridcount(i),iregridcount(i)
 801    format("for level ",i3, " average num. grids = ",f10.2," over ",i10,  &
               " regridding steps")
        write(outunit,802) i,numgrids(i)
 802    format("for level ",i3,"  current num. grids = ",i7)
      endif
    end do

    write(outunit,*)
    write(outunit,*)
    write(outunit,"('current  space usage = ',i12)") lentotsave
    write(outunit,"('maximum  space usage = ',i12)") lenmax
    write(outunit,"('need space dimension = ',i12,/)") lendim

    write(outunit,"('number of cells advanced for time integration = ',f20.6)")&
                    rvol
    do level = 1,mxnest
        write(outunit,"(3x,'# cells advanced on level ',i4,' = ',f20.2)") &
                    level, rvoll(level)
    enddo

    write(outunit,"('number of cells advanced for error estimation = ',f20.6,/)") &
                     evol
    if (evol + rvol > 0.d0) then
        ratmet = rvol / (evol + rvol) * 100.0d0
    else
        ratmet = 0.0d0
    endif
    write(outunit,"(' percentage of cells advanced in time  = ', f10.2)") ratmet
    write(outunit,"(' maximum Courant number seen = ', f10.2)") cflmax

    write(outunit,"(//,' ------  end of AMRCLAW integration --------  ')")

    ! Close output and debug files.
    close(outunit)
    close(dbugunit)

#ifdef USE_PDAF
    call mpi_finalize(mpierr)
#endif

end program amr2

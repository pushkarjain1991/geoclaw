! ::::::::::::::::::::: flag2refine ::::::::::::::::::::::::::::::::::
!
! User routine to control flagging of points for refinement.
!
! Specific for GeoClaw for tsunami applications and related problems
!
!
! The logical function allowflag(x,y,t) is called to
! check whether further refinement at this level is allowed in this cell
! at this time.
!
!    q   = grid values including ghost cells (bndry vals at specified
!          time have already been set, so can use ghost cell values too)
!
!  aux   = aux array on this grid patch
!
! amrflags  = array to be flagged with either the value
!             DONTFLAG (no refinement needed)  or
!             DOFLAG   (refinement desired)
!
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine flag2refine2(mx,my,mbc,mbuff,meqn,maux,xlower,ylower,dx,dy,t,level, &
                       tolsp,q,aux,amrflags,DONTFLAG,DOFLAG)
    
    use amr_module, only: mxnest, t0
    use geoclaw_module, only:dry_tolerance, sea_level
    use geoclaw_module, only: spherical_distance, coordinate_system

    use topo_module, only: tlowtopo,thitopo,xlowtopo,xhitopo,ylowtopo,yhitopo
    use topo_module, only: minleveltopo,mtopofiles

    use topo_module, only: tfdtopo,xlowdtopo,xhidtopo,ylowdtopo,yhidtopo
    use topo_module, only: minleveldtopo,num_dtopo

    use qinit_module, only: x_low_qinit,x_hi_qinit,y_low_qinit,y_hi_qinit
    use qinit_module, only: min_level_qinit,qinit_type

    use storm_module, only: storm_type, wind_refine, R_refine, storm_location
    use storm_module, only: wind_forcing, wind_index, wind_refine

    use regions_module, only: num_regions, regions
    use refinement_module

#ifdef USE_PDAF
      use mod_parallel, only: mype_world, MPIerr, mpi_comm_world, n_modeltasks
      use mpi, only: mpi_real8
      use mod_assimilation, only: stepnow_pdaf, assimilate_step, &
      dim_ens, regrid_assim
      !USE PDAF_mod_filtermpi, only: MPI_REALTYPE
#endif

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: mx,my,mbc,meqn,maux,level,mbuff
    real(kind=8), intent(in) :: xlower,ylower,dx,dy,t,tolsp

    real(kind=8), intent(in) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Flagging
    real(kind=8), intent(in out) :: amrflags(1-mbuff:mx+mbuff,1-mbuff:my+mbuff)
    real(kind=8), intent(in) :: DONTFLAG
    real(kind=8), intent(in) :: DOFLAG

    logical :: allowflag
    external allowflag

    ! Generic locals
    integer :: i,j,m
    real(kind=8) :: x_c,y_c,x_low,y_low,x_hi,y_hi
    real(kind=8) :: speed, eta, ds

    ! Storm specific variables
    real(kind=8) :: R_eye(2), wind_speed

#ifdef USE_PDAF
    integer :: gsize
    !real(kind=8), dimension(5000) :: recvbuf = 0
    !real(kind=8), dimension(2500) :: sendbuf = 0
    real(kind=8), allocatable :: recvbuf(:)
    real(kind=8), allocatable, dimension(:,:) :: reshaped_buffer
    real(kind=8), allocatable :: sendbuf(:)
    real(kind=8), allocatable :: flagunion(:)
    logical, allocatable :: mask(:)
    integer :: alloc_stat
    
    integer :: sendcount, recvcount
    integer, parameter :: root = 0
    integer :: recvsize
    character(len=1) :: ensstr
    character(len=7) :: file1
    character(len=7) :: file2
    character(len=17) :: filename
    integer :: csize
#endif

    ! Initialize flags
    amrflags = DONTFLAG


    ! Loop over interior points on this grid
    ! (i,j) grid cell is [x_low,x_hi] x [y_low,y_hi], cell center at (x_c,y_c)
    y_loop: do j=1,my
        y_low = ylower + (j - 1) * dy
        y_c = ylower + (j - 0.5d0) * dy
        y_hi = ylower + j * dy

        x_loop: do i = 1,mx
            x_low = xlower + (i - 1) * dx
            x_c = xlower + (i - 0.5d0) * dx
            x_hi = xlower + i * dx

            ! The following conditions are only checked in the horizontal and
            ! override the allowflag routine

            ! ************* Storm Based Refinement ****************
            ! Check to see if we are some specified distance from the eye of
            ! the storm and refine if we are
            if (storm_type > 0) then
                R_eye = storm_location(t)
                do m=1,size(R_refine,1)
                    if (coordinate_system == 2) then
                        ds = spherical_distance(x_c, y_c, R_eye(1), R_eye(2))
                    else
                        ds = sqrt((x_c - R_eye(1))**2 + (y_c - R_eye(2))**2)
                    end if
                    
                    if ( ds < R_refine(m) .and. level <= m ) then
                        amrflags(i,j) = DOFLAG
                        cycle x_loop
                    endif
                enddo
                
                ! Refine based on wind speed
                if (wind_forcing) then
                    wind_speed = sqrt(aux(wind_index,i,j)**2 + aux(wind_index+1,i,j)**2)
                    do m=1,size(wind_refine,1)
                        if ((wind_speed > wind_refine(m)) .and. (level <= m)) then
                            amrflags(i,j) = DOFLAG
                            cycle x_loop
                        endif
                    enddo
                endif
            endif
            ! *****************************************************

            ! Check to see if refinement is forced in any topography file region:
            do m=1,mtopofiles
                if (level < minleveltopo(m) .and. t >= tlowtopo(m) .and. t <= thitopo(m)) then
                    if (  x_hi > xlowtopo(m) .and. x_low < xhitopo(m) .and. &
                          y_hi > ylowtopo(m) .and. y_low < yhitopo(m) ) then

                        amrflags(i,j) = DOFLAG
                        cycle x_loop
                    endif
                endif
            enddo

            ! Check to see if refinement is forced in any other region:
            do m=1,num_regions
                if (level < regions(m)%min_level .and. &
                    t >= regions(m)%t_low .and. t <= regions(m)%t_hi) then
                    if (x_hi > regions(m)%x_low .and. x_low < regions(m)%x_hi .and. &
                        y_hi > regions(m)%y_low .and. y_low < regions(m)%y_hi ) then

                        amrflags(i,j) = DOFLAG
                        cycle x_loop
                    endif
                endif
            enddo

            ! Check if we're in the dtopo region and need to refine:
            ! force refinement to level minleveldtopo
            do m = 1,num_dtopo
                if (level < minleveldtopo(m).and. &
                    t <= tfdtopo(m) .and. & !t.ge.t0dtopo(m).and.
                    x_hi > xlowdtopo(m) .and. x_low < xhidtopo(m).and. &
                    y_hi > ylowdtopo(m) .and. y_low < yhidtopo(m)) then

                    amrflags(i,j) = DOFLAG
                    cycle x_loop
                endif
            enddo

            ! Check if we're in the region where initial perturbation is
            ! specified and need to force refinement:
            ! This assumes that t0 = 0.d0, should really be t0 but we do
            ! not have access to that parameter in this routine
            if (qinit_type > 0 .and. t == t0) then
                if (level < min_level_qinit .and. &
                    x_hi > x_low_qinit .and. x_low < x_hi_qinit .and. &
                    y_hi > y_low_qinit .and. y_low < y_hi_qinit) then

                    amrflags(i,j) = DOFLAG
                    cycle x_loop
                endif
            endif

            ! -----------------------------------------------------------------
            ! Refinement not forced, so check if it is allowed and if so,
            ! check if there is a reason to flag this point:
            if (allowflag(x_c,y_c,t,level)) then

                if (q(1,i,j) > dry_tolerance) then
                    eta = q(1,i,j) + aux(1,i,j)

                    ! Check wave criteria
                    if (abs(eta - sea_level) > wave_tolerance) then
                        ! Check to see if we are near shore
                        if (q(1,i,j) < deep_depth) then
                            amrflags(i,j) = DOFLAG
                            cycle x_loop
                        ! Check if we are allowed to flag in deep water
                        ! anyway
                        else if (level < max_level_deep) then
                            amrflags(i,j) = DOFLAG
                            cycle x_loop
                        endif
                    endif

                    ! Check speed criteria, note that it might be useful to
                    ! also have a per layer criteria since this is not
                    ! gradient based
                    speed = sqrt(q(2,i,j)**2 + q(3,i,j)**2) / q(1,i,j)
                    do m=1,min(size(speed_tolerance),mxnest)
                        if (speed > speed_tolerance(m) .and. level <= m) then
                            amrflags(i,j) = DOFLAG
                            cycle x_loop
                        endif
                    enddo
                endif
            endif

        enddo x_loop
    enddo y_loop

#ifdef USE_PDAF
    !write(ensstr,*) mype_world
    if(regrid_assim .eqv. .true.) then
    if (stepnow_pdaf == assimilate_step .or. t == t0) then
        print *, "reached here123 ", mype_world
        if (level == 1) then
            
            sendcount = mx*my
            recvcount = mx*my
            gsize = n_modeltasks
            recvsize = gsize*sendcount
            call MPI_barrier(mpi_comm_world, mpierr)
            print *, "flag2refine debugging"
            print *, "mype = ", mype_world, "level = ", level, "nx,ny = ", mx, my
            call MPI_barrier(mpi_comm_world, mpierr)
            if(allocated(sendbuf)) deallocate(sendbuf)
            allocate(sendbuf(sendcount))

            if (mype_world == 0) then
                allocate(recvbuf(recvsize))
            endif
            do i = 1,mx
                do j = 1,my
                    sendbuf(i+(j-1)*mx) = amrflags(i,j)
                enddo
            enddo
            

            if (xlower==0.0) then
                if (ylower==0.0) then
                    if (mype_world == 0) then
                        open(unit=46, file='file0', status='replace')
                        do j = 1,mx*my
                            write(46,*) sendbuf(j)
                        enddo
                        close(46)
                    endif
                    
                    if (mype_world == 1) then
                        open(unit=47, file='file1', status='replace')
                        do i = 1,mx*my
                            write(47,*) sendbuf(i)
                        enddo
                        close(47)
                    endif
                endif
            endif
            
            print *, "Gathering... ", mype_world
            call mpi_gather(amrflags(1:mx, 1:my), sendcount, mpi_real8, recvbuf, recvcount, mpi_real8, root, mpi_comm_world, mpierr)
            print *, "done gathering... ", mype_world
            !call MPI_barrier(mpi_comm_world, mpierr)

            !Comparing flags of different ensembles
            if (mype_world == 0) then
                if (xlower==0.0) then
                    if (ylower==0.0) then
                        !print *, "size", size(recvbuf), xlower, ylower
                        !print *, "amrflags = ", amrflags(:,:)
                        open(unit=45, file='file_combined', status='replace')
                        do i=1,recvsize
                            write(45,*) recvbuf(i)
                        enddo
                        !write(45,*) recvbuf
                        close(45)
                    endif
                endif
            endif

            !RESHAPING recvbuf
            if (mype_world == 0) then
                !Allocating memory for reshaped buffer
                allocate(reshaped_buffer(mx*my, gsize))

                !Reshaping recvbuf
                reshaped_buffer = reshape(recvbuf,(/mx*my, gsize/))

                !Writing to file for debugging
                if (xlower==0.0) then
                    if (ylower==0.0) then
                        print *, "recvbuf_shape = ", shape(recvbuf)
                        print *, "reshaped_buffer_shape = ", shape(reshaped_buffer)
                        open(unit=48, file='reshaped_buffer', status='replace')
                        do i = 1,mx*my
                            write(48,*) reshaped_buffer(i,:)
                        enddo
                        close(48)
                    endif
                endif
            endif

            if(allocated(flagunion)) deallocate(flagunion)
            allocate(flagunion(mx*my))

            !Perform boolean operation
            if (mype_world == 0) then
                flagunion = 0.0

                if(allocated(mask)) deallocate(mask)
                allocate(mask(mx*my))
                mask = any(reshaped_buffer .eq. 2.0, 2)
                where (mask .eqv. .true.) 
                    flagunion = 2.0
                endwhere
                
                if (xlower==0.0) then
                    if (ylower==0.0) then
                        open(unit=48, file='flag_union', status='replace')
                        do i = 1,mx*my
                            write(48,*) flagunion(i)
                        enddo
                        close(48)
                    endif
                endif
            endif

            !Broadcast flagunion
            print *, "Broadcasting..."
            call mpi_bcast(flagunion, mx*my, mpi_real8, root, mpi_comm_world, mpierr)
                
            if (mype_world == 1) then
                if (xlower==0.0) then
                    if (ylower==0.0) then
                        open(unit=48, file='check_broadcast', status='replace')
                        do i = 1,mx*my
                            write(48,*) flagunion(i)
                        enddo
                        close(48)
                    endif
                endif
            endif

            !Write flagunion to amrflags
            amrflags(1:mx, 1:my) = reshape(flagunion,(/mx,my/))
            !print *, "new amrflag = ", amrflags(1:mx, 1:my)

             
            !if (mype_world == 0) then
            !    allocate(reshaped_amrflag(2500,2), stat=alloc_stat)
            !    if (alloc_stat /= 0) then
            !            print *, "Memory not allocated"
            !    endif
            !!    reshaped_amrflag = reshape(amrflags,(/2500,2/))
            !!    print *, "reshaped_amrflag", reshaped_amrflag
            !    deallocate(reshaped_amrflag)
            !endif
            !
            !if (mype_world == 0) then
            !    deallocate(recvbuf)
            !endif
        endif
    endif
    endif
#endif
end subroutine flag2refine2

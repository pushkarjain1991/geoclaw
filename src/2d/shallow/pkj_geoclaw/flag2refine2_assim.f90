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
subroutine flag2refine2_assim(mx,my,mbc,mbuff,meqn,maux,xlower,ylower,dx,dy,t,level, &
                       tolsp,q,aux,amrflags,DONTFLAG,DOFLAG)
    
    use common_level, only: num_regions_assim, regions_assim
    use mod_assimilation, only: same_final_grid

    use refinement_module
      
    use mod_parallel, only: mype_world, MPIerr, mpi_comm_world, n_modeltasks
    use mpi, only: mpi_real8


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


    ! Generic locals
    integer :: i,j,m
    real(kind=8) :: x_c,y_c,x_low,y_low,x_hi,y_hi
    
    ! Delete from here
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
    ! Delete till here

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

            ! Check to see if refinement is forced in any other region:
            do m=1,num_regions_assim
                !if (level < regions_assim(m)%min_level .and. &
                !    t >= regions_assim(m)%t_low .and. t <= regions(m)%t_hi) then
                if (level < regions_assim(m)%min_level) then
                    if (x_hi > regions_assim(m)%x_low .and. x_low < regions_assim(m)%x_hi .and. &
                        y_hi > regions_assim(m)%y_low .and. y_low < regions_assim(m)%y_hi ) then

                        amrflags(i,j) = DOFLAG
                        cycle x_loop
                    endif
                endif
            enddo

        enddo x_loop
    enddo y_loop

!            print *, "mype = ", mype_world, "level = ", level, "nx,ny = ", mx, my
    !write(ensstr,*) mype_world
    if ((level == 2) .and. (same_final_grid .eqv. .true.)) then
            
      sendcount = mx*my
      recvcount = mx*my
      gsize = n_modeltasks
      recvsize = gsize*sendcount
            !print *, "flag2refine debugging"
            !print *, "mype = ", mype_world, "level = ", level, "nx,ny = ", mx, my
            
      if(allocated(sendbuf)) deallocate(sendbuf)
      allocate(sendbuf(sendcount))

      if (mype_world == 0) allocate(recvbuf(recvsize))
            
      do i = 1,mx
        do j = 1,my
          sendbuf(i+(j-1)*mx) = amrflags(i,j)
        enddo
      enddo
            
      print *, "Gathering... ", mype_world
      call mpi_gather(amrflags(1:mx, 1:my), sendcount, mpi_real8, &
           recvbuf, recvcount, mpi_real8, root, mpi_comm_world, &
           mpierr)
      print *, "done gathering... ", mype_world

      if (mype_world == 0) then
      !Allocating memory for reshaped buffer
        allocate(reshaped_buffer(mx*my, gsize))
 
        !Reshaping recvbuf
        reshaped_buffer = reshape(recvbuf,(/mx*my, gsize/))

      endif

      if(allocated(flagunion)) deallocate(flagunion)
      allocate(flagunion(mx*my))

      !Perform boolean operation
      if (mype_world == 0) then
        flagunion = 0.0
        if(allocated(mask)) deallocate(mask)
        allocate(mask(mx*my))
        !mask = any(reshaped_buffer .eq. 2.0, 2)
        mask = any(reshaped_buffer .eq. DOFLAG, 2)
        where (mask .eqv. .true.) 
        !flagunion = 2.0
          flagunion = DOFLAG
        endwhere
      endif

      !Broadcast flagunion
      print *, "Broadcasting..."
      call mpi_bcast(flagunion, mx*my, mpi_real8, root, mpi_comm_world, mpierr)
                
      !Write flagunion to amrflags
      amrflags(1:mx, 1:my) = reshape(flagunion,(/mx,my/))
      !print *, "new amrflag = ", amrflags(1:mx, 1:my)
    endif
    
end subroutine flag2refine2_assim

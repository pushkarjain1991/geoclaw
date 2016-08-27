module qinit_module

    use amr_module, only: rinfinity

    implicit none
    save

    logical :: module_setup = .false.
    
    ! Type of q initialization
    integer, public :: qinit_type
    
    ! Work array
    real(kind=8), private, allocatable :: qinit(:)

    ! Geometry
    real(kind=8) :: x_low_qinit
    real(kind=8) :: y_low_qinit
    real(kind=8) :: t_low_qinit
    real(kind=8) :: x_hi_qinit
    real(kind=8) :: y_hi_qinit
    real(kind=8) :: t_hi_qinit
    real(kind=8) :: dx_qinit
    real(kind=8) :: dy_qinit
    
    integer, private :: mx_qinit
    integer, private :: my_qinit
    integer :: min_level_qinit
    integer :: max_level_qinit

contains

    subroutine set_qinit(fname)
    
        use geoclaw_module, only: GEO_PARM_UNIT
    
        implicit none
        
        ! Subroutine arguments
        character(len=*), optional, intent(in) :: fname
        
        ! File handling
        integer, parameter :: unit = 7
        character(len=150) :: qinit_fname
        
        if (.not.module_setup) then
            write(GEO_PARM_UNIT,*) ' '
            write(GEO_PARM_UNIT,*) '--------------------------------------------'
            write(GEO_PARM_UNIT,*) 'SETQINIT:'
            write(GEO_PARM_UNIT,*) '-------------'
            
            ! Open the data file
            if (present(fname)) then
                call opendatafile(unit,fname)
            else
                call opendatafile(unit,"qinit.data")
            endif
            
            read(unit,"(i1)") qinit_type
            if (qinit_type == 0) then
                ! No perturbation specified
                write(GEO_PARM_UNIT,*)  '  qinit_type = 0, no perturbation'
                print *,'  qinit_type = 0, no perturbation'
                return
            endif
            read(unit,*) qinit_fname
            read(unit,"(2i2)") min_level_qinit, max_level_qinit

            write(GEO_PARM_UNIT,*) '   min_level, max_level, qinit_fname:'
            write(GEO_PARM_UNIT,*)  min_level_qinit, max_level_qinit, qinit_fname
            
            call read_qinit(qinit_fname)

            module_setup = .true.
        end if
    
    end subroutine set_qinit


    subroutine add_perturbation(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    
        use geoclaw_module, only: sea_level, coordinate_system
        use amr_module, only: mcapa
    
        implicit none
    
        ! Subroutine arguments
        integer, intent(in) :: meqn,mbc,mx,my,maux
        real(kind=8), intent(in) :: xlower,ylower,dx,dy
        real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
        real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
        
        ! Local
        integer :: i,j
        real(kind=8) :: ximc,xim,x,xip,xipc,yjmc,yjm,y,yjp,yjpc,dq
        
        ! Topography integral function
        real(kind=8) :: topointegral
        
        if (qinit_type > 0) then
            do i=1-mbc,mx+mbc
                x = xlower + (i-0.5d0)*dx
                xim = x - 0.5d0*dx
                xip = x + 0.5d0*dx
                do j=1-mbc,my+mbc
                    y = ylower + (j-0.5d0)*dy
                    yjm = y - 0.5d0*dy
                    yjp = y + 0.5d0*dy

                    ! Check to see if we are in the qinit region at this grid point
                    if ((xip > x_low_qinit).and.(xim < x_hi_qinit).and.  &
                        (yjp > y_low_qinit).and.(yjm < y_hi_qinit)) then

                        xipc=min(xip,x_hi_qinit)
                        ximc=max(xim,x_low_qinit)

                        yjpc=min(yjp,y_hi_qinit)
                        yjmc=max(yjm,y_low_qinit)

                        dq = topointegral(ximc,xipc,yjmc,yjpc,x_low_qinit, &
                                          y_low_qinit,dx_qinit,dy_qinit,mx_qinit, &
                                          my_qinit,qinit,1)
                        if (coordinate_system == 2) then
                            dq = dq / ((xipc-ximc)*(yjpc-yjmc)*aux(mcapa,i,j))
                        else
                            dq = dq / ((xipc-ximc)*(yjpc-yjmc))
                        endif 

                        if (qinit_type < 4) then 
                            if (aux(1,i,j) <= sea_level) then
                                q(qinit_type,i,j) = q(qinit_type,i,j) + dq
                            endif
                        else if (qinit_type == 4) then
                            q(1,i,j) = max(dq-aux(1,i,j),0.d0)
                        endif
                    endif
                enddo
            enddo
        endif
        
    end subroutine add_perturbation

        
    ! currently only supports one file type:
    ! x,y,z values, one per line in standard order from NW corner to SE
    ! z is perturbation from standard depth h,hu,hv set in qinit_geo,
    ! if iqinit = 1,2, or 3 respectively.
    ! if iqinit = 4, the z column corresponds to the definition of the 
    ! surface elevation eta. The depth is then set as q(i,j,1)=max(eta-b,0)
    subroutine read_qinit(fname)
    
        use geoclaw_module, only: GEO_PARM_UNIT
        
        implicit none
        
        ! Subroutine arguments
        character(len=150) :: fname
        
        ! Data file opening
        integer, parameter :: unit = 19
        integer :: i,num_points,status
        double precision :: x,y
        
        print *,'  '
        print *,'Reading qinit data from file  ', fname
        print *,'  '

        write(GEO_PARM_UNIT,*) '  '
        write(GEO_PARM_UNIT,*) 'Reading qinit data from'
        write(GEO_PARM_UNIT,*) fname
        write(GEO_PARM_UNIT,*) '  '
        
        open(unit=unit, file=fname, iostat=status, status="unknown", &
             form='formatted',action="read")
        if ( status /= 0 ) then
            print *,"Error opening file", fname
            stop
        endif
        
        ! Initialize counters
        num_points = 0
        mx_qinit = 0
        
        ! Read in first values, determines x_low and y_hi
        read(unit,*) x_low_qinit,y_hi_qinit
        num_points = num_points + 1
        mx_qinit = mx_qinit + 1
        
        ! Sweep through first row figuring out mx
        y = y_hi_qinit
        do while (y_hi_qinit == y)
            read(unit,*) x,y
            num_points = num_points + 1
            mx_qinit = mx_qinit + 1
        enddo
        ! We over count by one in the above loop
        mx_qinit = mx_qinit - 1
        
        ! Continue to count the rest of the lines
        do
            read(unit,*,iostat=status) x,y
            if (status /= 0) exit
            num_points = num_points + 1
        enddo
        if (status > 0) then
            print *,"ERROR:  Error reading qinit file ",fname
            stop
        endif
        
        ! Extract rest of geometry
        x_hi_qinit = x
        y_low_qinit = y
        my_qinit = num_points / mx_qinit
        dx_qinit = (x_hi_qinit - x_low_qinit) / (mx_qinit-1)
        dy_qinit = (y_hi_qinit - y_low_qinit) / (my_qinit-1)
        
        rewind(unit)
        allocate(qinit(num_points))
        
        ! Read and store the data this time
        do i=1,num_points
            read(unit,*) x,y,qinit(i)
        enddo
        close(unit)
        
    end subroutine read_qinit
    
    subroutine add_momentum(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    
        !use geoclaw_module, only: sea_level, coordinate_system
        !use amr_module, only: mcapa
    
        implicit none
    
        ! Subroutine arguments
        logical there
        logical ens_tracker_exist
        integer, intent(in) :: meqn,mbc,mx,my,maux
        real(kind=8), intent(in) :: xlower,ylower,dx,dy
        real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
        real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
        real(kind=8) :: height, xveldata, yveldata, etadata
        ! Local
        integer :: i,j, p, Reason
        real(kind=8) :: xim,x,xip,yjm,y,yjp
        integer next_ens_number, ens_number
        character(*), parameter :: fileplace = "/h2/pkjain/Desktop/Pushkar/clawpack/geoclaw&
        /examples/tsunami/PDAF-D_V1.11.1/testsuite/&
        simplified_offline_geoclaw/python_script/"
        character(145) :: fileplace2
        character(100) :: totallength
        character(2) :: str_ens_number
        character(LEN=*), PARAMETER :: FMT1 = "(E26.16, E26.16, E26.16, E26.16)"
        !--------------------------------------------------------------!
        !ENSEMBLE TRACKER
        !--------------------------------------------------------------!
        !Open and read ens_tracker
        inquire(file=fileplace//"ens_tracker", exist=ens_tracker_exist)
        if(.NOT. ens_tracker_exist) then
          return
        endif
        open(unit = 45, FILE=fileplace//"ens_tracker")
        read(45,*)ens_number
        print *, "ens_number read by geoclaw is ", ens_number
        close(45)
        next_ens_number = ens_number + 1 ! Calculate next ensemble number
        ! Overwrite ens_tracker with next ensemble number
        open(unit = 46, FILE=fileplace//"ens_tracker", status='replace')
        write(46,*)next_ens_number
        close(46)

        write(str_ens_number,'(I2)')ens_number
        fileplace2 = adjustl(fileplace//"fort.q0012"//"_ens_"//TRIM(ADJUSTL(str_ens_number)))
        !open(unit=2, FILE=fileplace2!)
        !close(2)
      
        inquire(file=fileplace2, exist=there)
        PRINT *,fileplace2,"FILE IS",there
        if ( there ) then
            open(unit=2, FILE=fileplace2)
            do p=1,9
                read(2,*)
             enddo
        
            if (qinit_type > 0) then
                !print *,x_low_qinit, x_hi_qinit, y_low_qinit, y_hi_qinit, mbc
                !print *,xlower, ylower
                !do i=1-mbc,mx+mbc
                !    x = xlower + (i-0.5d0)*dx
                !    xim = x - 0.5d0*dx
                !    xip = x + 0.5d0*dx
                !    do j=1-mbc,my+mbc
                !        y = ylower + (j-0.5d0)*dy
                !        yjm = y - 0.5d0*dy
                !        yjp = y + 0.5d0*dy
                do j=1-mbc,my+mbc
                    y = ylower + (j-0.5d0)*dy
                    yjm = y - 0.5d0*dy
                    yjp = y + 0.5d0*dy
                    do i=1-mbc,mx+mbc
                        x = xlower + (i-0.5d0)*dx
                        xim = x - 0.5d0*dx
                        xip = x + 0.5d0*dx
                    ! Check to see if we are in the qinit region at this grid point
                        if ((xip > x_low_qinit).and.(xim < x_hi_qinit).and.  &
                            (yjp > y_low_qinit).and.(yjm < y_hi_qinit)) then
                            read(2,*) q(1,i,j),q(2,i,j), q(3,i,j), height
                            !read(2,*) height,q(3,i,j), q(2,i,j), etadata
                            !print *,etadata, height
                            !print *,xip,xim,yjp,yjm,height
                        endif
                    enddo
                enddo
            endif
         ! Write rest of data to dummy file
         !open(unit=3,FILE=fileplace//"dummy")
         !READ(2,"(A)",IOSTAT=Reason)totallength
         !DO
         !    READ(2,"(A)",IOSTAT=Reason)totallength
         !    
         !    IF (Reason > 0) THEN
         !        print *,"something went wrong"
         !    ELSE IF (Reason < 0) THEN
         !        print *, "END OF FILE REACHED"
         !        close(3)
         !        exit
         !    ELSE
         !        !print *,totallength
         !        write(3,*)totallength
         !    ENDIF
         !ENDDO
         !   !close(2,status='delete')
         !   close(2)
        !CALL RENAME(fileplace//"fort.q0012",fileplace//"fort.orig")
        !close(3)
        !CALL RENAME(fileplace//"dummy",fileplace//"fort.q0012")
        endif
        
    end subroutine add_momentum
    

    subroutine add_momentum2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    
        !use geoclaw_module, only: sea_level, coordinate_system
        !use amr_module, only: mcapa
    
        implicit none
    
        ! Subroutine arguments
        logical there
        logical ens_tracker_exist
        integer, intent(in) :: meqn,mbc,mx,my,maux
        real(kind=8), intent(in) :: xlower,ylower,dx,dy
        real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
        real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
        real(kind=8) :: height, xveldata, yveldata, etadata
        ! Local
        integer :: i,j, p, Reason
        real(kind=8) :: xim,x,xip,yjm,y,yjp
        integer next_ens_number, ens_number
        character(*), parameter :: fileplace = "/h2/pkjain/Desktop/Pushkar/clawpack/geoclaw&
        /examples/tsunami/PDAF-D_V1.11.1/testsuite/&
        simplified_offline_geoclaw/python_script/"
        character(144) :: fileplace2
        !character(100) :: totallength
        character(1) :: str_ens_number
        character(LEN=*), PARAMETER :: FMT1 = "(E26.16, E26.16, E26.16, E26.16)"
        !--------------------------------------------------------------!
        !ENSEMBLE TRACKER
        !--------------------------------------------------------------!
        !Open and read ens_tracker
        inquire(file=fileplace//"ens_tracker", exist=ens_tracker_exist)
        if(.NOT. ens_tracker_exist) then
          return
        endif
        open(unit = 45, FILE=fileplace//"ens_tracker")
        read(45,*)ens_number
        print *, "ens_number read by geoclaw is ", ens_number
        close(45)
        next_ens_number = ens_number + 1 ! Calculate next ensemble number
        ! Overwrite ens_tracker with next ensemble number
        open(unit = 46, FILE=fileplace//"ens_tracker", status='replace')
        write(46,*)next_ens_number
        close(46)

        write(str_ens_number,'(I2)')ens_number
        !fileplace2 = adjustl(fileplace//"fort.q0012"//"_ens_"//str_ens_number)
        fileplace2 = adjustl(fileplace//"fort.q0012last")
        !open(unit=2, FILE=fileplace2)
        !close(2)
      
        inquire(file=fileplace2, exist=there)
        PRINT *,fileplace2,"FILE IS",there
        if ( .NOT. there ) then
            return
        else
            open(unit=2, FILE=fileplace2)
            !open(unit=21, FILE=fileplace//"yoyo"//str_ens_number)
            do p=1,9
                read(2,*)
             enddo
        
            if (qinit_type > 0) then
                !print *,x_low_qinit, x_hi_qinit, y_low_qinit, y_hi_qinit, mbc
                !print *,xlower, ylower
                !do i=1-mbc,mx+mbc
                !    x = xlower + (i-0.5d0)*dx
                !    xim = x - 0.5d0*dx
                !    xip = x + 0.5d0*dx
                !    do j=1-mbc,my+mbc
                !        y = ylower + (j-0.5d0)*dy
                !        yjm = y - 0.5d0*dy
                !        yjp = y + 0.5d0*dy
                do j=1-mbc,my+mbc
                    y = ylower + (j-0.5d0)*dy
                    yjm = y - 0.5d0*dy
                    yjp = y + 0.5d0*dy
                    do i=1-mbc,mx+mbc
                        x = xlower + (i-0.5d0)*dx
                        xim = x - 0.5d0*dx
                        xip = x + 0.5d0*dx
                    ! Check to see if we are in the qinit region at this grid point
                        if ((xip > x_low_qinit).and.(xim < x_hi_qinit).and.  &
                            (yjp > y_low_qinit).and.(yjm < y_hi_qinit)) then
                                
                                read(2,*) q(1,i,j),q(2,i,j), q(3,i,j), height                           
                                !write(21,FMT1) q(1,i,j),q(3,i,j), q(2,i,j), height                           
                                print *, q(1,i,j),q(2,i,j), q(3,i,j), height
                        endif
                    enddo
                enddo
            endif
            close(2)
            close(21)
         endif
         
    end subroutine add_momentum2

end module qinit_module

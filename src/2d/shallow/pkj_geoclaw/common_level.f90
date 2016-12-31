!*************************************************************************
!   > File Name: common_level.f90
!    > Author: Pushkar Kumar Jain
!    > Mail: pushkarjain1991@utexas.edu
! ************************************************************************/
module common_level
    use amr_module
    use regions_module, only: region_type
    use mod_parallel, only: mype_world, mpi_comm_world, &
            mpierr, n_modeltasks
    use mpi, only: mpi_int, mpi_sum, mpi_real8
    use checker_same_column
    implicit none
    
    real(kind=8),allocatable :: x_low_array(:)
    real(kind=8),allocatable :: y_low_array(:)
    real(kind=8),allocatable :: x_hi_array(:)
    real(kind=8),allocatable :: y_hi_array(:)
    integer, allocatable :: max_level_array(:)
    integer, allocatable :: min_level_array(:)
    integer :: region_counter
    
    real(kind=8),allocatable :: gather_x_low_array(:)
    real(kind=8),allocatable :: gather_y_low_array(:)
    real(kind=8),allocatable :: gather_x_hi_array(:)
    real(kind=8),allocatable :: gather_y_hi_array(:)
    integer, allocatable :: gather_max_level_array(:)
    integer, allocatable :: gather_min_level_array(:)
    
    type(region_type), allocatable :: regions_assim(:)
    integer :: num_regions_assim
  
    !interface check_same_column
    !  module procedure check_same_column_R
    !end interface

    contains

    subroutine num_region_count()
        integer :: level
        integer :: mptr
         
         level = 1
         region_counter = 0
 66      if (level .gt. lfine) go to 91
            mptr = lstart(level)
 71         if (mptr .eq. 0) go to 81
                region_counter = region_counter + 1
                mptr = node(levelptr, mptr)
             go to 71
 81      level = level + 1
         go to 66
 91      continue
         !print *, "num_regions = ", region_counter
    end subroutine num_region_count 

    subroutine set_regions_assim_arrays()

        integer :: ii,mptr,Ntot,Ntot_l,j_pkj,i_pkj
        integer :: i_mod,j_mod,loc,locaux,mitot,mjtot,nx,ny 
        integer :: ivar,i,j,iaux
        integer :: level
        integer :: region_putter
        real(kind=8) :: dx,dy,xlow,ylow
 
         call num_region_count()
         print *, "num_regions in set_regions_array= ", region_counter,mype_world
    
         if(allocated(x_low_array)) deallocate(x_low_array)
         if(allocated(y_low_array)) deallocate(y_low_array)
         if(allocated(x_hi_array)) deallocate(x_hi_array)
         if(allocated(y_hi_array)) deallocate(y_hi_array)
         if(allocated(max_level_array)) deallocate(max_level_array)
         if(allocated(min_level_array)) deallocate(min_level_array)
         allocate(x_low_array(region_counter))
         allocate(x_hi_array(region_counter))
         allocate(y_low_array(region_counter))
         allocate(y_hi_array(region_counter))
         allocate(max_level_array(region_counter))
         allocate(min_level_array(region_counter))

         region_putter = 1
         level = 1
 67      if (level .gt. lfine) go to 92
            mptr = lstart(level)
 72         if (mptr .eq. 0) go to 82

             nx = node(ndihi,mptr) - node(ndilo, mptr) + 1
             ny = node(ndjhi,mptr) - node(ndjlo, mptr) + 1
             xlow = rnode(cornxlo,mptr)
             ylow = rnode(cornylo,mptr)
             dx = hxposs(level)
             dy = hyposs(level)
             min_level_array(region_putter) = level
             max_level_array(region_putter) = level
             x_low_array(region_putter) = xlow
             y_low_array(region_putter) = ylow
             x_hi_array(region_putter) = xlow + nx*dx
             y_hi_array(region_putter) = ylow + ny*dy
               
             mptr = node(levelptr, mptr)
             region_putter = region_putter + 1
             go to 72
 82      level = level + 1
         go to 67
 92      continue

    end subroutine set_regions_assim_arrays

    subroutine create_gatherv_array()

    integer :: sum_ensemble_num_regions
    integer :: num_regions_array(n_modeltasks)
    integer :: root = 0
    integer :: sendcount
    integer :: i1
    integer, allocatable :: gather_displs(:)
    
    !Gather sizes of num_regions to root
    call mpi_gather(region_counter, 1, mpi_int, num_regions_array, &
            1, mpi_int, root, mpi_comm_world, mpierr)
    if(mype_world == 0) then
        print *, "Gathered ", num_regions_array
    endif
    
    !Get total number of regions at root
    !This will decide the size of the output array of the domain
    call mpi_allreduce(region_counter, sum_ensemble_num_regions,1, &
    mpi_int, mpi_sum, mpi_comm_world, mpierr)
    !sum_ensemble_num_regions = sum(num_regions_array)
    if(mype_world == 0) then
        print *, "Total regions = ", sum_ensemble_num_regions
    endif
    
    !Set gatherv parameters
    if (mype_world == root) then
        if(allocated(gather_displs)) deallocate(gather_displs)
        allocate(gather_displs(n_modeltasks))
        gather_displs(1) = 0
        do i1 = 2, n_modeltasks, 1
            gather_displs(i1) = gather_displs(i1-1) + num_regions_array(i1-1)
        end do
    endif

        if(allocated(gather_x_low_array)) deallocate(gather_x_low_array)
        if(allocated(gather_y_low_array)) deallocate(gather_y_low_array)
        if(allocated(gather_x_hi_array)) deallocate(gather_x_hi_array)
        if(allocated(gather_y_hi_array)) deallocate(gather_y_hi_array)
        if(allocated(gather_max_level_array)) deallocate(gather_max_level_array)
        if(allocated(gather_min_level_array)) deallocate(gather_min_level_array)
        
        allocate(gather_x_low_array(sum_ensemble_num_regions))
        allocate(gather_y_low_array(sum_ensemble_num_regions))
        allocate(gather_x_hi_array(sum_ensemble_num_regions))
        allocate(gather_y_hi_array(sum_ensemble_num_regions))
        allocate(gather_max_level_array(sum_ensemble_num_regions))
        allocate(gather_min_level_array(sum_ensemble_num_regions))

    print *, "Performing gahterv operation"
    call mpi_gatherv(x_low_array, region_counter, mpi_real8, &
            gather_x_low_array, num_regions_array, gather_displs, mpi_real8, &
            0, mpi_comm_world, mpierr)
    call mpi_gatherv(y_low_array, region_counter, mpi_real8, &
            gather_y_low_array, num_regions_array, gather_displs, mpi_real8, &
            0, mpi_comm_world, mpierr)
    call mpi_gatherv(x_hi_array, region_counter, mpi_real8, &
            gather_x_hi_array, num_regions_array, gather_displs, mpi_real8, &
            0, mpi_comm_world, mpierr)
    call mpi_gatherv(y_hi_array, region_counter, mpi_real8, &
            gather_y_hi_array, num_regions_array, gather_displs, mpi_real8, &
            0, mpi_comm_world, mpierr)
    call mpi_gatherv(max_level_array, region_counter, mpi_int, &
            gather_max_level_array, num_regions_array, gather_displs, mpi_int, &
            0, mpi_comm_world, mpierr)
    call mpi_gatherv(min_level_array, region_counter, mpi_int, &
            gather_min_level_array, num_regions_array, gather_displs, mpi_int, &
            0, mpi_comm_world, mpierr)
    print *, "Done gahterv operation"

    !At this stage, it is required to do rectangulation (Future work)
    ! Updated all the gather_arrays

    print *, "Performing broadcast operation"
    !Broadcast back the updated array to processors
    sendcount = size(gather_x_low_array)
    call mpi_bcast(gather_x_low_array, sendcount, mpi_real8, root, &
     mpi_comm_world, mpierr)
    call mpi_bcast(gather_y_low_array, sendcount, mpi_real8, root, &
     mpi_comm_world, mpierr)
    call mpi_bcast(gather_x_hi_array, sendcount, mpi_real8, root, &
     mpi_comm_world, mpierr)
    call mpi_bcast(gather_y_hi_array, sendcount, mpi_real8, root, &
     mpi_comm_world, mpierr)
    call mpi_bcast(gather_max_level_array, sendcount, mpi_int, root, &
     mpi_comm_world, mpierr)
    call mpi_bcast(gather_min_level_array, sendcount, mpi_int, root, &
     mpi_comm_world, mpierr)
     !print *, "bcasted field = ", gather_x_low_array(:)
     print *, "Done broadcast operation"
    end subroutine create_gatherv_array
    




    subroutine check_gatherv_array()

      integer :: sum_ensemble_num_regions
      integer :: num_regions_array(n_modeltasks)
      integer :: root = 0
      integer :: sendcount
      integer :: i, i1
      integer :: local_left, local_right
      integer, allocatable :: gather_displs(:)
    
      real(kind=8),allocatable :: reshaped_gather_x_low_array(:,:)
      real(kind=8),allocatable :: reshaped_gather_y_low_array(:,:)
      real(kind=8),allocatable :: reshaped_gather_x_hi_array(:,:)
      real(kind=8),allocatable :: reshaped_gather_y_hi_array(:,:)
      integer, allocatable :: reshaped_gather_max_level_array(:,:)
      integer, allocatable :: reshaped_gather_min_level_array(:,:)

      !Gather sizes of num_regions to root
      call mpi_gather(region_counter, 1, mpi_int, num_regions_array, &
            1, mpi_int, root, mpi_comm_world, mpierr)
      if(mype_world == 0) then
        print *, "Gathered ", num_regions_array
        if(all(num_regions_array /= num_regions_array(1))) then
          stop "Unequal regions in check_gatherv"
        endif
      endif
    
      !Get total number of regions at root
      !This will decide the size of the output array of the domain
      call mpi_allreduce(region_counter, sum_ensemble_num_regions,1, &
      mpi_int, mpi_sum, mpi_comm_world, mpierr)
      !sum_ensemble_num_regions = sum(num_regions_array)
      if(mype_world == 0) then
        print *, "Total regions = ", sum_ensemble_num_regions
      endif
    
      !Set gatherv parameters
      if (mype_world == 0) then
        if(allocated(gather_displs)) deallocate(gather_displs)
        allocate(gather_displs(n_modeltasks))
        gather_displs(1) = 0
        do i1 = 2, n_modeltasks, 1
          gather_displs(i1) = gather_displs(i1-1) + num_regions_array(i1-1)
        end do
      endif

      if(allocated(gather_x_low_array)) deallocate(gather_x_low_array)
      if(allocated(gather_y_low_array)) deallocate(gather_y_low_array)
      if(allocated(gather_x_hi_array)) deallocate(gather_x_hi_array)
      if(allocated(gather_y_hi_array)) deallocate(gather_y_hi_array)
      if(allocated(gather_max_level_array)) deallocate(gather_max_level_array)
      if(allocated(gather_min_level_array)) deallocate(gather_min_level_array)
      
      allocate(gather_x_low_array(sum_ensemble_num_regions))
      allocate(gather_y_low_array(sum_ensemble_num_regions))
      allocate(gather_x_hi_array(sum_ensemble_num_regions))
      allocate(gather_y_hi_array(sum_ensemble_num_regions))
      allocate(gather_max_level_array(sum_ensemble_num_regions))
      allocate(gather_min_level_array(sum_ensemble_num_regions))

      print *, "Performing gatherv operation"
      call mpi_gatherv(x_low_array, region_counter, mpi_real8, &
            gather_x_low_array, num_regions_array, gather_displs, mpi_real8, &
            0, mpi_comm_world, mpierr)
      call mpi_gatherv(y_low_array, region_counter, mpi_real8, &
            gather_y_low_array, num_regions_array, gather_displs, mpi_real8, &
            0, mpi_comm_world, mpierr)
      call mpi_gatherv(x_hi_array, region_counter, mpi_real8, &
            gather_x_hi_array, num_regions_array, gather_displs, mpi_real8, &
            0, mpi_comm_world, mpierr)
      call mpi_gatherv(y_hi_array, region_counter, mpi_real8, &
            gather_y_hi_array, num_regions_array, gather_displs, mpi_real8, &
            0, mpi_comm_world, mpierr)
      call mpi_gatherv(max_level_array, region_counter, mpi_int, &
            gather_max_level_array, num_regions_array, gather_displs, mpi_int, &
            0, mpi_comm_world, mpierr)
      call mpi_gatherv(min_level_array, region_counter, mpi_int, &
            gather_min_level_array, num_regions_array, gather_displs, mpi_int, &
            0, mpi_comm_world, mpierr)
      print *, "Done gatherv operation"

      if(mype_world == 0) then
        if(allocated(reshaped_gather_x_low_array)) deallocate(reshaped_gather_x_low_array)
        if(allocated(reshaped_gather_y_low_array)) deallocate(reshaped_gather_y_low_array)
        if(allocated(reshaped_gather_x_hi_array)) deallocate(reshaped_gather_x_hi_array)
        if(allocated(reshaped_gather_y_hi_array)) deallocate(reshaped_gather_y_hi_array)
        if(allocated(reshaped_gather_max_level_array)) deallocate(reshaped_gather_max_level_array)
        if(allocated(reshaped_gather_min_level_array)) deallocate(reshaped_gather_min_level_array)
      
        allocate(reshaped_gather_x_low_array(num_regions_array(1),n_modeltasks))
        reshaped_gather_x_low_array = &
        reshape(gather_x_low_array, (/num_regions_array(1), n_modeltasks/))

        allocate(reshaped_gather_y_low_array(num_regions_array(1),n_modeltasks))
        reshaped_gather_y_low_array = &
        reshape(gather_y_low_array, (/num_regions_array(1), n_modeltasks/))
        
        allocate(reshaped_gather_x_hi_array(num_regions_array(1),n_modeltasks))
        reshaped_gather_x_hi_array = &
        reshape(gather_x_hi_array, (/num_regions_array(1), n_modeltasks/))
        
        allocate(reshaped_gather_y_hi_array(num_regions_array(1),n_modeltasks))
        reshaped_gather_y_hi_array = &
        reshape(gather_y_hi_array, (/num_regions_array(1), n_modeltasks/))
        
        allocate(reshaped_gather_max_level_array(num_regions_array(1),n_modeltasks))
        reshaped_gather_max_level_array = &
        reshape(gather_max_level_array, (/num_regions_array(1), n_modeltasks/))
        
        allocate(reshaped_gather_min_level_array(num_regions_array(1),n_modeltasks))
        reshaped_gather_min_level_array = &
        reshape(gather_min_level_array, (/num_regions_array(1), n_modeltasks/))
        
        !if(check_same_column(reshaped_gather_x_low_array) .eqv. .false.) stop
        call check_same_column(reshaped_gather_x_low_array)
        call check_same_column(reshaped_gather_y_low_array)
        call check_same_column(reshaped_gather_x_hi_array)
        call check_same_column(reshaped_gather_y_hi_array)
        call check_same_column(reshaped_gather_max_level_array)
        call check_same_column(reshaped_gather_min_level_array)
        
        deallocate(reshaped_gather_x_low_array)
        deallocate(reshaped_gather_y_low_array)
        deallocate(reshaped_gather_x_hi_array)
        deallocate(reshaped_gather_y_hi_array)
        deallocate(reshaped_gather_max_level_array)
        deallocate(reshaped_gather_min_level_array)
      endif

    end subroutine check_gatherv_array


    subroutine update_regions()
        
        !real(kind=8), intent(in) :: time
        !logical, intent(in) :: regrid_assim

        integer :: total_num_regions
        integer :: i
        type(region_type), allocatable :: temp_regions(:)

          !Get current regions data
          total_num_regions = size(gather_x_low_array)
          print *, "Total regions = ", total_num_regions
          if(allocated(temp_regions)) deallocate(temp_regions)
          allocate(temp_regions(total_num_regions))

          !Update Geoclaw regions
          do i = 1, total_num_regions
            temp_regions(i)%min_level = gather_min_level_array(i)
            temp_regions(i)%max_level = gather_max_level_array(i)
            temp_regions(i)%t_low = 0.02
            temp_regions(i)%t_hi = 1.0E+10
            temp_regions(i)%x_low = gather_x_low_array(i)
            temp_regions(i)%x_hi = gather_x_hi_array(i)
            temp_regions(i)%y_low = gather_y_low_array(i)
            temp_regions(i)%y_hi = gather_y_hi_array(i)
          enddo

          if(allocated(regions_assim)) deallocate(regions_assim)
          allocate(regions_assim(total_num_regions))
          regions_assim(:) = temp_regions(:)
          num_regions_assim = total_num_regions

          deallocate(temp_regions)
        
    end subroutine update_regions
end module


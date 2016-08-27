module mapdomain
        use amr_module, only: numgrids, rnode, node, lstart, ndihi, ndilo, &
        ndjhi, ndjlo, mstart, cornxlo, cornylo, levelptr
        use sortarr

        implicit none

        !Local integers
        integer :: a
        integer :: total_nx, total_ny
        
        !real(kind=8) :: total_array(10000) = 0.0
        !integer :: global_to_totalarray_map(10000) = 0
        !integer :: ptr_location(4) !Pointer

        !numgrids(1) - no. of domains
        !numcells(1) - total number of cells
        !real(kind=8),allocatable :: corner_array(:,:)
        !real(kind=8),allocatable :: mptr_array(:)


contains
        ! This subroutine gives the ordered mptr array for the domains going
        ! from left to right and upwards
        subroutine get_ordered_array(mptr_array,ordered_mptr_array)

                integer,allocatable, intent(out) :: mptr_array(:)
                integer, allocatable, intent(out) :: ordered_mptr_array(:)

                real(kind=8), allocatable :: corner_array(:,:)
                integer :: i_pkj
                integer :: mptr

                allocate(mptr_array(numgrids(1)))
                allocate(ordered_mptr_array(numgrids(1)))
                allocate(corner_array(numgrids(1),2))

                mptr = lstart(1)

                do i_pkj=1,numgrids(1)
                mptr_array(i_pkj)=mptr
                corner_array(i_pkj,1) = rnode(cornxlo,mptr)
                corner_array(i_pkj,2) = rnode(cornylo,mptr)
                mptr = node(levelptr,mptr)
                enddo

                call set_global(corner_array,mptr_array,ordered_mptr_array)

        end subroutine get_ordered_array

        !subroutine allocate_space(ptr_location, ordered_mptr_array)
        !Allocates space for nx*ny cells and sets pointers for each domain 
        subroutine allocate_space(ordered_mptr_array, ptr_location)
                integer, intent(in) :: ordered_mptr_array(:)
                integer,intent(inout) :: ptr_location(:)

                !Local variables
                integer :: mptr1,i, nx, ny

                !allocate(ptr_location(numgrids(1)))
                ptr_location(1)=1

                do i = 2,size(ordered_mptr_array)
                    mptr1 = ordered_mptr_array(i-1)
                    nx = node(ndihi,mptr1) - node(ndilo,mptr1) + 1
                    ny = node(ndjhi,mptr1) - node(ndjlo,mptr1) + 1
                    ptr_location(i) = ptr_location(i-1) + nx*ny
                end do
                print *,"ptr_location = ",ptr_location

        end subroutine allocate_space

        !subroutine check_domain(i_global, j_global,a)
        subroutine check_domain(i_global, j_global, domain_num)
                integer, intent(in) :: i_global
                integer, intent(in) :: j_global
                integer, intent(out) :: domain_num

                integer :: domain_num_x
                integer :: domain_num_y
                
                !To make it generic, need a good subroutine to identify
                !domain_num_y and domain_num_x from any number of cells in the
                !local domains. Basically 50 must be replaced somehow
                domain_num_y = ((j_global-1)/50) + 1
                domain_num_x = ((i_global-1)/50) + 1
                domain_num = domain_num_y + (domain_num_x-1)*2
                print *,"Global_coord = ", j_global, i_global
                print *, "domain coordinate - ",domain_num_x, domain_num_y
                print *,"Domain number = ", domain_num
        end subroutine check_domain

        subroutine traverse_global(field_array, ordered_mptr_array, &
                        total_array,global_to_totalarray_map, ptr_location)
                !State values in the physical domain
                real(kind=8), intent(in) :: field_array(:,:)
                integer, intent(in) :: ordered_mptr_array(:)
                real(kind=8), allocatable, intent(out) :: total_array(:)
                integer, allocatable, intent(out) :: global_to_totalarray_map(:)
                integer, allocatable, intent(out) :: ptr_location(:) !Pointer

                ! Local variables
                integer :: i_global, j_global
                integer :: counter = 0
                integer :: domain_num
                real(kind=8) :: field_value

                allocate(total_array(10000))
                allocate(global_to_totalarray_map(10000))
                allocate(ptr_location(numgrids(1)))
                call allocate_space(ordered_mptr_array, ptr_location)
                !For making it generic, 100 has to be replaced with total number
                !of cells in x and y direction in the domain
                do i_global = 1,100
                    do j_global = 1,100

                        counter = counter + 1
                        field_value = field_array(i_global, j_global)
                        call check_domain(i_global, j_global, domain_num)

                        call fill_value(domain_num, field_value, counter, &
                                ptr_location, total_array, &
                                global_to_totalarray_map)
                    end do
                end do
                print *,"muhaha"
                total_nx = node(ndihi,mstart) - node(ndilo,mstart) + 1
                total_ny = node(ndjhi,mstart) - node(ndjlo,mstart) + 1
                print *,mstart, total_nx, total_ny
                print *,"global2totalarray = ", global_to_totalarray_map
                open(unit=24,file="total_array")
                do i_global = 1,10000
                    write(24,*) total_array(i_global)
                end do
                close(24)
        end subroutine traverse_global

        subroutine fill_value(domain_num, field_val, counter, &
                        ptr_location, total_array, &
                        global_to_totalarray_map)
                integer, intent(in) :: domain_num
                real(kind=8),intent(in) :: field_val
                integer, intent(in) :: counter
                integer,intent(inout) :: ptr_location(:)
                real(kind=8), intent(inout) :: total_array(:)
                integer, intent(inout) :: global_to_totalarray_map(:)
                
                !real(kind=8), intent(inout) :: total_array(:)
                global_to_totalarray_map(counter) = ptr_location(domain_num)
                total_array(ptr_location(domain_num)) = field_val
                ptr_location(domain_num) = ptr_location(domain_num) + 1
                print *,domain_num,ptr_location(domain_num)

        end subroutine fill_value

!        subroutine obtain_local_domain_width(ordered_mptr_array)
!          integer, intent(in) :: ordered_mptr_array(:)
!
!          integer :: i2, mptr1
!
!          do i2 = 1,size(ordered_mptr_array)
!                  mptr1 = ordered_mptr_array(i2)
!                  nx = node(ndihi,mptr1) - node(ndilo, mptr1) + 1
!                  ny = node(ndjhi,mptr1) - node(ndjlo, mptr1) + 1
!                  Ntot = nx*ny
!                  DO i_pkj = 1,Ntot
!                      READ(24,*) field((i2-1)*Ntot+i_pkj)
!                  !    print *,i_pkj
!                  ENDDO
!          ENDDO
!        end subroutine obtain_local_domain_width


end module mapdomain

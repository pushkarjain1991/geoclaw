!*************************************************************************
!   > File Name: set_global_coordinate_array.f90
!    > Author: Pushkar Kumar Jain
!    > Mail: pushkarjain1991@utexas.edu
! ************************************************************************/

subroutine set_global_coordinate_array(array_size, global_coordinate_array)
        use amr_module

        implicit none
  
        integer, intent(in) :: array_size
        real(kind=8), intent(inout) :: global_coordinate_array(2, array_size)
        integer :: nx, ny, level, mptr
        real(kind=8) :: xlow, ylow, dx,dy
        integer :: cell_cnt
        integer :: i,j


         !Count the total number of cells
         !The value will be used to allocate field
         cell_cnt = 1
         level = 1
 66      if (level .gt. lfine) go to 91
            mptr = lstart(level)
 71         if (mptr .eq. 0) go to 81
                nx = node(ndihi,mptr) - node(ndilo, mptr) + 1
                ny = node(ndjhi,mptr) - node(ndjlo, mptr) + 1
                xlow = rnode(cornxlo,mptr)
                ylow = rnode(cornylo,mptr)
                dx = hxposs(level)
                dy = hyposs(level)

                do j=1,ny
                  do i=1,nx
                  global_coordinate_array(1, cell_cnt) = xlow + (i-0.5)*dx
                  global_coordinate_array(2, cell_cnt) = ylow + (j-0.5)*dy
                  !global_coordinate(1, cell_cnt) = xlow + (i-0.5)*dx
                  !global_coordinate(2, cell_cnt) = ylow + (j-0.5)*dy
                  cell_cnt = cell_cnt + 1
                  enddo
                enddo

                mptr = node(levelptr, mptr)
             go to 71
 81         level = level + 1
         go to 66
 91      continue

end subroutine set_global_coordinate_array

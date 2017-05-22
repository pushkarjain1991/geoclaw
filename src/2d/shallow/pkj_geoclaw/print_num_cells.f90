!*************************************************************************
!   > File Name: alloc2field.f90
!    > Author: Shawn Lin
!    > Mail: lin_yuxiang@utexas.edu
!    > Created Time: Mon 01 Aug 2016 04:04:44 PM CDT
! ************************************************************************/

subroutine print_num_cells()
    use amr_module
    use mod_parallel, only: mype_world

        implicit none

        integer :: mptr, nx, ny, level
        integer :: i,j
        integer :: ncells_temp, level_num_cells

         !Count the total number of cells
         !The value will be used to allocate field
         ncells_temp=0
         level_num_cells = 0
         level = 1
 66      if (level .gt. lfine) go to 91
            mptr = lstart(level)
 71         if (mptr .eq. 0) go to 81
                nx = node(ndihi,mptr) - node(ndilo, mptr) + 1
                ny = node(ndjhi,mptr) - node(ndjlo, mptr) + 1

                ncells_temp = ncells_temp + nx*ny
                level_num_cells = level_num_cells + nx*ny
                mptr = node(levelptr, mptr)
             go to 71
 81      print *, mype_world, "level ", level, level_num_cells, " cells"
         level_num_cells = 0
         level = level + 1
         go to 66
 91      continue

         print *, "Total_cells = ", mype_world, ncells_temp, lfine

    end subroutine


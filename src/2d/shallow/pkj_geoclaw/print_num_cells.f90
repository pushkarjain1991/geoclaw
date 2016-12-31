!*************************************************************************
!   > File Name: alloc2field.f90
!    > Author: Shawn Lin
!    > Mail: lin_yuxiang@utexas.edu
!    > Created Time: Mon 01 Aug 2016 04:04:44 PM CDT
! ************************************************************************/

subroutine print_num_cells(nvar,naux)
    use amr_module
    use geoclaw_module, only: dry_tolerance
    use mod_parallel, only: mype_world

        implicit none

        integer,intent(in)::nvar
        integer,intent(in)::naux
        integer :: ii,mptr,Ntot,Ntot_l,j_pkj,i_pkj
        integer :: i_mod,j_mod,loc,locaux,mitot,mjtot,nx,ny 
        integer :: ivar,i,j,iaux
        integer :: iadd,iaddaux
        integer :: cell_cnt
        integer :: level
        integer :: ncells_temp, level_num_cells
        real(kind=8),allocatable :: temp_field(:)
        integer,allocatable :: temp_wet_cell_index(:)
        integer :: wet_cell_cnt
        logical :: field_is_total_height

        iadd(ivar,i,j)=loc+ivar-1+nvar*((j-1)*mitot+i-1)
        iaddaux(iaux,i,j)=locaux+iaux-1+naux*(i-1)+naux*mitot*(j-1)

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
                mitot   = nx + 2*nghost
                mjtot   = ny + 2*nghost
                loc     = node(store1, mptr)
                locaux  = node(storeaux,mptr)

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


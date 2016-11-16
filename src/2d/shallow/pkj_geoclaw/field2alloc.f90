!*************************************************************************
!   > File Name: field2alloc.f90
!    > Author: Shawn Lin
!    > Mail: lin_yuxiang@utexas.edu
!    > Created Time: Mon 01 Aug 2016 04:07:06 PM CDT
! ************************************************************************/

subroutine field2alloc(nvar,naux,analyze_water)
        use amr_module

#ifdef USE_PDAF
        use mod_model,only:field, wet_cell_index
        use geoclaw_module, only: dry_tolerance
!        use mod_assimilation,only: first_assimilation
!       use sortarr
!        use mapdomain,only: get_ordered_array
!        use mod_assimilation, only: ordered_mptr_array

        implicit none
        
        integer,intent(in) :: nvar
        integer,intent(in) :: naux
        logical,intent(in) :: analyze_water

!        integer,allocatable :: mptr_array(:)
!        integer,allocatable :: ordered_mptr_array(:)
        integer :: mptr, j_pkj, i_pkj
        integer :: i_mod, j_mod, loc, locaux
        integer :: mitot, mjtot,nx,ny 
        integer :: ivar,i,j,iaux
        integer :: iadd,iaddaux
        integer :: cell_cnt
        integer :: level
        integer :: wet_index_ptr
        logical :: first_water_analysis=.true.

        iadd(ivar,i,j)=loc+ivar-1+nvar*((j-1)*mitot+i-1)
        iaddaux(iaux,i,j)=locaux+iaux-1+naux*(i-1)+naux*mitot*(j-1)

        if(analyze_water) then
          cell_cnt = 1
          level = 1
          wet_index_ptr = 1
 65       if (level .gt. lfine) go to 90
            mptr = lstart(level)
 70         if (mptr .eq. 0) go to 80
              nx = node(ndihi,mptr) - node(ndilo, mptr) + 1
              ny = node(ndjhi,mptr) - node(ndjlo, mptr) + 1
              mitot   = nx + 2*nghost
              mjtot   = ny + 2*nghost
              loc     = node(store1, mptr)
              locaux  = node(storeaux,mptr)

              do j_pkj = nghost+1, mjtot-nghost
                do i_pkj = nghost+1, mitot-nghost
                  if (cell_cnt == wet_cell_index(wet_index_ptr)) then   
                    alloc(iadd(1,i_pkj,j_pkj)) = &
                    field(wet_index_ptr) - alloc(iaddaux(1,i_pkj, j_pkj))
                    wet_index_ptr = wet_index_ptr + 1
                  endif
                  cell_cnt = cell_cnt + 1

!                     !wse update might cause total hieght /= 0 
!                     if (abs(alloc(iadd(1,i_pkj,j_pkj))) < 1.0) then
!                         print *, "alloc non negative activated"
!                         alloc(iadd(1,i_pkj,j_pkj)) = 0.d0
!                     endif
                             
                enddo
              enddo
              mptr = node(levelptr, mptr)
            go to 70
 80         level = level + 1
          go to 65

 90       continue
     
        else
          cell_cnt = 1
          level = 1
 66       if (level .gt. lfine) go to 91
            mptr = lstart(level)
 71         if (mptr .eq. 0) go to 81
              nx = node(ndihi,mptr) - node(ndilo, mptr) + 1
              ny = node(ndjhi,mptr) - node(ndjlo, mptr) + 1
              mitot   = nx + 2*nghost
              mjtot   = ny + 2*nghost
              loc     = node(store1, mptr)
              locaux  = node(storeaux,mptr)

              do j_pkj = nghost+1, mjtot-nghost
                do i_pkj = nghost+1, mitot-nghost
                  alloc(iadd(1,i_pkj,j_pkj)) = &
                  field(cell_cnt) - alloc(iaddaux(1,i_pkj, j_pkj))
                  cell_cnt = cell_cnt + 1

                  if(first_water_analysis) then
                     !wse update might cause total hieght /= 0 
                     if (abs(alloc(iadd(1,i_pkj,j_pkj))) < dry_tolerance) then
                         print *, "alloc non negative activated"
                         alloc(iadd(1,i_pkj,j_pkj)) = 0.d0
                     endif
                  endif
                             
                enddo
              enddo
              mptr = node(levelptr, mptr)
            go to 71
 81         level = level + 1
          go to 66

 91     continue
!        first_water_analysis=.false.

     endif
#endif
     end subroutine field2alloc

!*************************************************************************
!   > File Name: field2alloc.f90
!    > Author: Shawn Lin
!    > Mail: lin_yuxiang@utexas.edu
!    > Created Time: Mon 01 Aug 2016 04:07:06 PM CDT
! ************************************************************************/

subroutine field2alloc(nvar,naux,isoutput)
        use amr_module

#ifdef USE_PDAF
        use mod_model,only:field,mptr_array, ordered_mptr_array
!        use mod_assimilation,only: first_assimilation
!       use sortarr
        use mapdomain,only: get_ordered_array
!        use mod_assimilation, only: ordered_mptr_array

        implicit none
        
        integer,intent(in) :: nvar
        integer,intent(in) :: naux
        logical,intent(in) :: isoutput

!        integer,allocatable :: mptr_array(:)
!        integer,allocatable :: ordered_mptr_array(:)
        integer :: ii,mptr,Ntot,Ntot_l,j_pkj,i_pkj
        integer :: i_mod,j_mod,loc,locaux,mitot,mjtot,nx,ny 
        integer :: ivar,i,j,iaux
        integer :: iadd,iaddaux
        integer :: cell_cnt = 1
        integer :: level
        iadd(ivar,i,j)=loc+ivar-1+nvar*((j-1)*mitot+i-1)

        iaddaux(iaux,i,j)=locaux+iaux-1+naux*(i-1)+naux*mitot*(j-1)
         level = 1
 65      if (level .gt. lfine) go to 90
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
                        
                alloc(iadd(1,i_pkj,j_pkj)) = &
                field(cell_cnt)-alloc(iaddaux(1,i_pkj, j_pkj))
                cell_cnt = cell_cnt + 1
            enddo
         enddo

            mptr = node(levelptr, mptr)
            go to 70
 80      level = level + 1
         go to 65

 90     continue
#endif
     end subroutine field2alloc

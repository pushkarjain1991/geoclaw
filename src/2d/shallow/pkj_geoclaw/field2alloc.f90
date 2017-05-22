!*************************************************************************
!   > File Name: field2alloc.f90
!    > Author: Shawn Lin
!    > Mail: lin_yuxiang@utexas.edu
!    > Created Time: Mon 01 Aug 2016 04:07:06 PM CDT
! ************************************************************************/

subroutine field2alloc(nvar,naux)
        use amr_module
        use mod_model,only:field

        !implicit none
        implicit double precision (a-h,o-z)
        
        integer,intent(in) :: nvar
        integer,intent(in) :: naux

        integer :: j_pkj, i_pkj
        !integer :: mptr
        !integer :: loc, locaux
        !integer :: mitot, mjtot,nx,ny 
        !integer :: ivar,i,j,iaux
        !integer :: iadd,iaddaux
        integer :: cell_cnt
        !integer :: level

        iadd(ivar,i,j)=loc+ivar-1+nvar*((j-1)*mitot+i-1)
        iaddaux(iaux,i,j)=locaux+iaux-1+naux*(i-1)+naux*mitot*(j-1)

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
                  !alloc(iadd(1,i_pkj,j_pkj)) = &
                  !field(cell_cnt)
                  alloc(iadd(1,i_pkj,j_pkj)) = &
                  field(cell_cnt) - alloc(iaddaux(1,i_pkj, j_pkj))
                  cell_cnt = cell_cnt + 1
                enddo
              enddo
              mptr = node(levelptr, mptr)
            go to 71
 81         level = level + 1
          go to 66

 91     continue

     end subroutine field2alloc

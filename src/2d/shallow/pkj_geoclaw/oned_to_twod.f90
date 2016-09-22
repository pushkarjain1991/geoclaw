subroutine ind1d_to_coord2d(index_1d, coord_2d)
#ifdef USE_PDAF
      use amr_module
      use mod_model,only: ordered_mptr_array
!      use mapdomain,only:get_ordered_array
      implicit none

      integer, intent(in) :: index_1d
!      integer, intent(in) :: numdomains_in_x, numdomains_in_y
      real(kind=8), intent(inout) :: coord_2d(2)

      !Local variables
!      integer, allocatable :: mptr_array(:), ordered_mptr_array(:)
      integer :: ii,mptr,nx,ny,row_l,column_l,level
      real(kind=8) :: xlow,ylow,dx,dy
      integer :: cellnum , remaining_cells
!      integer :: domain_num_row, domain_num_column
!      integer :: local_row, local_column
! :::::::::::1d to 2d :::::::::::::::::::::
! get the 2d coordinate of 1d index
! ::::::::::::::::::::::::::::::::::::::::::
!       print *, "executing the oned to twod"
!       call get_ordered_array(mptr_array,ordered_mptr_array)
       cellnum=0
!       ind=minloc(numgrids,1)

!       dx=hxposs(ind)
!       dy=hyposs(ind)
       do ii=1,size(ordered_mptr_array)
           mptr=ordered_mptr_array(ii)
           nx=node(ndihi,mptr)-node(ndilo,mptr)+1
           ny=node(ndjhi,mptr)-node(ndjlo,mptr)+1
           level=node(nestlevel,mptr)
           dx=hxposs(level)
           dy=hyposs(level)

           cellnum=cellnum+nx*ny
           if (index_1d<=cellnum) then
               xlow=rnode(cornxlo,mptr)
               ylow=rnode(cornylo,mptr)
               remaining_cells=index_1d-(cellnum-nx*ny)
               column_l=(remaining_cells -1)/nx+1
               row_l=mod(remaining_cells -1,nx)+1
               coord_2d(1)=xlow+(row_l-0.5)*dx
               coord_2d(2)=ylow+(column_l-0.5)*dy
               return
           endif
       enddo

#endif








!      numcells = 0
!      do i = 1,4
!          !numcells = numcells +  nx(i)*ny(i)
!          numcells = numcells +  50*50
!
!          if (index_1d <= numcells) then
!                  !remaining_cells = index_1d - (numcells - nx(i)*ny(i))
!                  remaining_cells = index_1d - (numcells - 50*50)
!!                  print *,"remaining cells = ", remaining_cells
!
!                  !local_row = remaining_cells/nx(i)
!                  local_row = (remaining_cells-1)/50 + 1
!                  !local_column = remaining_cells - local_row*nx(i)
!                  !local_column = remaining_cells - local_row*50
!                  local_column = mod(remaining_cells-1, 50) + 1
!!                  print *,"local row = ", local_row
!!                  print *,"local column = ", local_column
!!                  print *, "i = ", i
!
!                  domain_num_row = (i-1)/numdomains_in_x + 1
!                  !domain_num_column = i - (domain_num_row-1)*numdomains_in_x
!                  domain_num_column = mod(i-1, numdomains_in_x)+1
!!                  print *,"hello"
!!                  print *,"domain_num_row = ", domain_num_row
!!                  print *,"domain_num_col = ", domain_num_column
!                  !domain_num_column = mod(i, numdomains_in_x)
!
!                  !index_2d(1) = (domain_num_row-1)*ny(i) + local_row
!                  index_2d(1) = (domain_num_row-1)*50 + local_row
!                  !print *,"xcoord = ", index_2d(1)
!                  !index_2d(2) = (domain_num_column-1)*nx(i) +&
!                  index_2d(2) = (domain_num_column-1)*50 + local_column
!                  !print *,"ycoord = ", index_2d(2)
!                  return
!          endif
!
!      end do

end subroutine


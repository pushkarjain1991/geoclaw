subroutine oned_to_twod(index_1d, numdomains_in_x,&
numdomains_in_y, index_2d)
      implicit none
      
      integer, intent(in) :: index_1d
      integer, intent(in) :: numdomains_in_x, numdomains_in_y
      integer, intent(inout) :: index_2d(2)

      !Local variables
      integer :: i
      integer :: numcells, remaining_cells 
      integer :: domain_num_row, domain_num_column
      integer :: local_row, local_column


      numcells = 0
      do i = 1,4
          !numcells = numcells +  nx(i)*ny(i)
          numcells = numcells +  50*50

          if (index_1d <= numcells) then
                  !remaining_cells = index_1d - (numcells - nx(i)*ny(i))
                  remaining_cells = index_1d - (numcells - 50*50)
                  print *,"remaining cells = ", remaining_cells

                  !local_row = remaining_cells/nx(i)
                  local_row = (remaining_cells-1)/50 + 1
                  !local_column = remaining_cells - local_row*nx(i)
                  !local_column = remaining_cells - local_row*50
                  local_column = mod(remaining_cells-1, 50) + 1
                  print *,"local row = ", local_row
                  print *,"local column = ", local_column
                  print *, "i = ", i

                  domain_num_row = (i-1)/numdomains_in_x + 1
                  domain_num_column = i - (domain_num_row-1)*numdomains_in_x
                  !domain_num_column = mod(i, numdomains_in_x)+1
                  print *,"hello"
                  print *,"domain_num_row = ", domain_num_row
                  print *,"domain_num_col = ", domain_num_column
                  !domain_num_column = mod(i, numdomains_in_x)

                  !index_2d(1) = (domain_num_row-1)*ny(i) + local_row
                  index_2d(1) = (domain_num_row-1)*50 + local_row
                  !print *,"xcoord = ", index_2d(1)
                  !index_2d(2) = (domain_num_column-1)*nx(i) +&
                  index_2d(2) = (domain_num_column-1)*50 + local_column
                  !print *,"ycoord = ", index_2d(2)
                  return
          endif
         
      end do

end subroutine 
      

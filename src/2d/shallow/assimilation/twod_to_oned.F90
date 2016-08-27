subroutine twod_to_oned(index_2d_row,index_2d_col, index_1d)
      implicit none

      integer, intent(in) :: index_2d_row, index_2d_col
      integer, intent(out) :: index_1d

      !Local variables
      integer :: domain_num, i
      integer :: cells_temp
      integer :: local_row, local_col
      integer :: domain_num_row, domain_num_col

      call check_domain2(index_2d_row, index_2d_col, domain_num)

      cells_temp = 0
      if (domain_num /= 1) then
              do i = 1,domain_num-1
              !cells_temp = cells_temp + nx(i)*ny(i)
              cells_temp = cells_temp + 50*50
              enddo
      end if

      domain_num_col = ((index_2d_col-1)/50) + 1
      domain_num_row = ((index_2d_row-1)/50) + 1
     
      !local_row = index_2d_row - (domain_num_row - 1)*50
      !local_col = index_2d_col - (domain_num_col - 1)*50
      local_row = mod(index_2d_row - 1, 50) + 1
      local_col = mod(index_2d_col - 1, 50) + 1
!index_1d = cells_temp + (index_2d(1) - 1)*nx(domain_num) + index_2d(2)
      index_1d = cells_temp + (local_row - 1)*50 + local_col

end subroutine twod_to_oned

!subroutine check_domain(i_global, j_global,a)
subroutine check_domain2(i_global, j_global, domain_num)
      implicit none
      integer, intent(in) :: i_global
      integer, intent(in) :: j_global
      integer, intent(inout) :: domain_num

      integer :: domain_num_x
      integer :: domain_num_y

      !To make it generic, need a good subroutine to identify
      !domain_num_y and domain_num_x from any number of cells in the
      !local domains. Basically 50 must be replaced somehow
      domain_num_y = ((j_global-1)/50) + 1
      domain_num_x = ((i_global-1)/50) + 1
      domain_num = domain_num_y + (domain_num_x-1)*2
      print *,"Twod to oned"
      print *,"Global_coord = ", j_global, i_global
      print *, "domain coordinate - ",domain_num_x, domain_num_y
      print *,"Domain number = ", domain_num
end subroutine check_domain2

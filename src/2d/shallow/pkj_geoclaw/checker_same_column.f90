module checker_same_column
! Checks if the 2d array has same columns

  interface check_same_column
    module procedure check_same_column_R, check_same_column_I 
  end interface


contains      

  subroutine check_same_column_R(two_d_array)
    
    implicit none

    real(kind=8), intent(in) :: two_d_array(:,:)

    integer :: i
    integer :: num_rows, num_cols

    num_rows = size(two_d_array, 1)
    num_cols = size(two_d_array, 2)

    do i =1, num_cols
      if(.not. (all(abs(two_d_array(:,i) &
                  - two_d_array(:,1)) <= 1.0d-06))) then
        print *, two_d_array
        stop "Matrix column not same in checker_same_column"
        return
      endif
    enddo

    return
  end subroutine check_same_column_R

  subroutine check_same_column_I(two_d_array)
    implicit none

    integer, intent(in) :: two_d_array(:,:)

    integer :: i
    integer :: num_rows, num_cols

    num_rows = size(two_d_array, 1)
    num_cols = size(two_d_array, 2)

    do i =1, num_cols
      if(.not. (all(abs(two_d_array(:,i) &
                  - two_d_array(:,1)) == 0))) then
        print *, two_d_array
        stop "Matrix column not same in checker_same_column"
        return
      endif
    enddo

    return
  end subroutine check_same_column_I

end module checker_same_column

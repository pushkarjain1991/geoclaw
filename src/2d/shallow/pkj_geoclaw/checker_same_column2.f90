module checker_same_column

  interface check_same_column
   
    logical function check_same_column_R(two_d_array)
      real(kind=8), intent(in) :: two_d_array(:,:)
    end function check_same_column_R

    logical function check_same_column_I(two_d_array)
      integer, intent(in) :: two_d_array(:,:)
    end function check_same_column_I
  end interface

end module checker_same_column

      

  logical function check_same_column_R(two_d_array)
    
    implicit none

    real(kind=8), intent(in) :: two_d_array(:,:)

    integer :: i
    integer :: num_rows, num_cols

    num_rows = size(two_d_array, 1)
    num_cols = size(two_d_array, 2)

    check_same_column_R = .true.
    do i =1, num_cols
      if(.not. (all(abs(two_d_array(:,i) &
                  - two_d_array(:,1)) <= 1.0d-06))) then
        print *, two_d_array
        check_same_column_R = .false.
        return
      endif
    enddo

    return
  end function check_same_column_R

  logical function check_same_column_I(two_d_array)
    implicit none

    integer, intent(in) :: two_d_array(:,:)

    integer :: i
    integer :: num_rows, num_cols

    num_rows = size(two_d_array, 1)
    num_cols = size(two_d_array, 2)

    check_same_column_I = .true.
    do i =1, num_cols
      if(.not. (all(abs(two_d_array(:,i) &
                  - two_d_array(:,1)) == 0))) then
        print *, two_d_array
        check_same_column_I = .false.
        return
      endif
    enddo

    return
  end function check_same_column_I

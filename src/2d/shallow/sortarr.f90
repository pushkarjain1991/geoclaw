MODULE sortarr
    implicit none
     
      contains

      subroutine set_global(A,inputid,outputid) 
      real(kind=8),intent(in) :: A(:,:)
      integer,intent(in) :: inputid(:)
      integer,intent(out) :: outputid(:)

      integer :: num_points
      real(kind=8) :: input_point(2)
      real(kind=8) :: output_point(2)
      real(kind=8),allocatable :: output_points(:,:)
      integer :: i,j
      integer :: m,n
      logical :: firsttime = .true.
      logical :: yesno
      integer :: temp

      num_points = size(A,2)
      !allocate(outputid(num_points))
      allocate(output_points(num_points,2))
     
      output_points = A
  
      do i = 1,size(inputid)
          input_point = A(i,:)
          if (firsttime .eqv. .true.) then
              outputid(1) = inputid(1)
              output_points(i,:) = A(i,:) 
              firsttime = .false.
          else
              outputid(i) = inputid(i)
              !Check the position of input point in outputid array
              do j = 1,i-1
                  output_point = output_points(j,:)
                  call shouldbebefore(input_point,output_point,yesno) 
                  if (yesno .eqv. .true.) then
                      call swap(outputid(j),outputid(i))
                     call swaparr(output_points(j,:),output_points(i,:))
                  else
                      output_points(i,:) = input_point
                  endif 
              enddo
          endif
      enddo
             
      !print *,"input array",(A(i,:),i=1,num_points)  
      do i=1,4
      print *,"input array",(A(i,j),j=1,2)  
      enddo
      print *,"input id",inputid
      do i=1,4
      print *,"output array",(output_points(i,j),j=1,2)
      enddo
      print *,"outputid", outputid

      !deallocate(outputid)
      deallocate(output_points)

      end subroutine set_global

      subroutine swap(a,b)
          implicit none
          integer, intent(inout) :: a,b
          integer :: c
          c = a
          a = b
          b = c
      end subroutine swap
      
      subroutine swaparr(a,b)
          implicit none
          real(kind=8), intent(inout) :: a(2),b(2)
          real(kind=8) :: c(2)
          c = a
          a = b
          b = c
      end subroutine swaparr

      subroutine shouldbebefore(a,b,yup)
          real(kind=8),intent(in) :: a(2),b(2)
          logical,intent(out) :: yup
          if (a(2)<b(2)) then
              yup = .true.
          else if (a(2)==b(2)) then
              if (a(1) < b(1)) then
                  yup = .true.
              else
                  yup = .false.
              endif
          else
              yup = .false.
          endif


      end subroutine shouldbebefore
END MODULE sortarr

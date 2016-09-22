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
        iadd(ivar,i,j)=loc+ivar-1+nvar*((j-1)*mitot+i-1)

        iaddaux(iaux,i,j)=locaux+iaux-1+naux*(i-1)+naux*mitot*(j-1)
!        call get_ordered_array()   
!#ifdef USE_PDAF
!        if (first_assimilation) call field_conversion()
!#endif
!        print *,field(2500)
!        print *,ordered_mptr_array
         Ntot_l=0
 !        print *,ordered_mptr_array
        do ii = 1,size(ordered_mptr_array)
           mptr = ordered_mptr_array(ii)
           nx = node(ndihi,mptr) - node(ndilo, mptr) + 1
           ny = node(ndjhi,mptr) - node(ndjlo, mptr) + 1
           mitot   = nx + 2*nghost
           mjtot   = ny + 2*nghost
           loc     = node(store1, mptr)
           locaux  = node(storeaux,mptr)

!          print *,"loc=",loc
!          print *,"mptr =",mptr
!          print *, "nx = ",nx
!          print *, "ny = ",ny
!          print *, "corners = ",corn1,corn2
!          print *, "mptr= ",ordered_mptr_array

          Ntot_l = Ntot_l+nx*ny
          do j_pkj = nghost+1, mjtot-nghost
            do i_pkj = nghost+1, mitot-nghost
              do ivar=1,nvar
                if (abs(alloc(iadd(ivar,i_pkj,j_pkj))) < 1d-90) then
                    alloc(iadd(ivar,i_pkj,j_pkj)) = 0.d0
                endif
              enddo
                        
                i_mod = i_pkj - nghost
                j_mod = j_pkj - nghost
!                ind_alloc=iadd(1,i_pkj,j_pkj)                
!               print *,ind_alloc
                alloc(iadd(1,i_pkj,j_pkj)) = &
                field(Ntot_l-nx*ny+i_mod+nx*(j_mod-1))-alloc(iaddaux(1,i_pkj, j_pkj))
            enddo
         enddo
       enddo
!       deallocate(field)
!      print *,'finish field to alloc'
      if (.not. isoutput) then ! only for output, not need to deallcate
          deallocate(mptr_array)
          deallocate(ordered_mptr_array)
      endif 
!      print *,alloc(iadd(1,3,3))
!        if (first_assimilation) deallocate(field)
#endif
     end subroutine field2alloc

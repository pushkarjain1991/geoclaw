!*************************************************************************
!   > File Name: alloc_field_exchange.f90
!    > Author: Shawn Lin
!    > Mail: lin_yuxiang@utexas.edu
!    > Created Time: Fri 29 Jul 2016 04:29:37 PM CDT
! ************************************************************************

module alloc_field_exchange
        use amr_module

#ifdef USE_PDAF
        use mod_model,only:field
        use mod_assimilation,only: first_assimilation
#endif
        use sortarr
        use mapdomain,only: get_ordered_array,field_conversion

        implicit none
        integer,allocatable :: mptr_array(:)
        integer,allocatable :: ordered_mptr_array(:)
        integer :: ii,mptr,Ntot,Ntot_l,j_pkj,i_pkj
        integer :: i_mod,j_mod,loc,locaux,mitot,mjtot,nx,ny 
        

contains

!    function iadd(ivar,i,j)
!        integer :: ivar,i,j,iadd
!        iadd =loc+ivar-1+nvar*((j-1)*mitot+i-1)
!    end function iadd
!    
!    function iaddaux(naux,i,j)
!        integer :: naux,i,j,iaddaux
!        iaddaux=locaux+iaux-1+naux*(i-1)+naux*mitot*(j-1)
!    end function iaddaux

    
    subroutine field2alloc(nvar,naux)

        implicit none
        
        integer,intent(in)::nvar
        integer,intent(in)::naux
        integer :: ivar,i,j,iaux
        integer :: iadd,iaddaux
        iadd(ivar,i,j)=loc+ivar-1+nvar*((j-1)*mitot+i-1)

        iaddaux(iaux,i,j)=locaux+iaux-1+naux*(i-1)+naux*mitot*(j-1)
        call get_ordered_array(mptr_array,ordered_mptr_array)   
#ifdef USE_PDAF
        if (first_assimilation) call field_conversion()
#endif
!        print *,field(5050)
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

          Ntot_l = nx*ny
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
                field((ii-1)*Ntot_l+i_mod+nx*(j_mod-1))-alloc(iaddaux(1,i_pkj, j_pkj))
            enddo
         enddo
       enddo
       deallocate(mptr_array)
       deallocate(ordered_mptr_array)
#ifdef USE_PDAF
        if (first_assimilation) deallocate(field)
#endif
     end subroutine field2alloc
     
     
     subroutine alloc2field(nvar,naux)


        implicit none
        integer,intent(in)::nvar
        integer,intent(in)::naux
        integer :: ivar,i,j,iaux
        integer :: iadd,iaddaux
        iadd(ivar,i,j)=loc+ivar-1+nvar*((j-1)*mitot+i-1)

        iaddaux(iaux,i,j)=locaux+iaux-1+naux*(i-1)+naux*mitot*(j-1)

           Ntot=0
           if (.not. allocated(field)) then
               call get_ordered_array(mptr_array,ordered_mptr_array)   
               do ii = 1,size(ordered_mptr_array)
                   mptr = ordered_mptr_array(ii)
                   nx = node(ndihi,mptr) - node(ndilo, mptr) + 1
                   ny = node(ndjhi,mptr) - node(ndjlo, mptr) + 1
                   Ntot=Ntot+nx*ny
               enddo
               allocate(field(Ntot))  !need to figure out the dimension of new field
               print *,"allocated field"
          endif 
           do ii = 1,size(ordered_mptr_array)
           mptr = ordered_mptr_array(ii)
           nx = node(ndihi,mptr) - node(ndilo, mptr) + 1
           ny = node(ndjhi,mptr) - node(ndjlo, mptr) + 1
           mitot   = nx + 2*nghost
           mjtot   = ny + 2*nghost
           loc     = node(store1, mptr)
           locaux  = node(storeaux,mptr)

!                  print *,"loc=",loc
!                  print *,"mptr =",mptr
!                  print *, "nx = ",nx
!                  print *, "ny = ",ny
!                  print *, "corners = ",corn1,corn2

          Ntot_l = nx*ny
          do j_pkj = nghost+1, mjtot-nghost
            do i_pkj = nghost+1, mitot-nghost
              do ivar=1,nvar
                if (abs(alloc(iadd(ivar,i_pkj,j_pkj))) < 1d-90) then
                    alloc(iadd(ivar,i_pkj,j_pkj)) = 0.d0
                endif
              enddo
                        
                i_mod = i_pkj - nghost
                j_mod = j_pkj - nghost
                field((ii-1)*Ntot+i_mod+nx*(j_mod-1)) =&
                alloc(iadd(1,i_pkj,j_pkj)) +alloc(iaddaux(1,i_pkj, j_pkj))
            enddo
          enddo
        enddo
       deallocate(mptr_array)
       deallocate(ordered_mptr_array)

       end subroutine alloc2field

end module alloc_field_exchange

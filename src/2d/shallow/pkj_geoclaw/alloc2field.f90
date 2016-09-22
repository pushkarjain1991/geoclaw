!*************************************************************************
!   > File Name: alloc2field.f90
!    > Author: Shawn Lin
!    > Mail: lin_yuxiang@utexas.edu
!    > Created Time: Mon 01 Aug 2016 04:04:44 PM CDT
! ************************************************************************/

subroutine alloc2field(nvar,naux)
    use amr_module

#ifdef USE_PDAF
        use mod_model,only:field, ordered_mptr_array
!        use mod_assimilation,only: first_assimilation
        use mod_assimilation,only: dim_state_p
        use mod_parallel,only: mype_world,mpi_comm_world,mpierr
        use sortarr
        use mapdomain,only: get_ordered_array


        implicit none
        integer,intent(in)::nvar
        integer,intent(in)::naux
!        integer,intent(in)::level
!        integer,allocatable :: mptr_array(:)
!        integer,allocatable :: ordered_mptr_array(:)
        integer :: ii,mptr,Ntot,Ntot_l,j_pkj,i_pkj
        integer :: i_mod,j_mod,loc,locaux,mitot,mjtot,nx,ny 
        integer :: ivar,i,j,iaux
        integer :: iadd,iaddaux
        iadd(ivar,i,j)=loc+ivar-1+nvar*((j-1)*mitot+i-1)

        iaddaux(iaux,i,j)=locaux+iaux-1+naux*(i-1)+naux*mitot*(j-1)

!        mptr=4
!        loc=node(store1,mptr)
!        locaux=node(storeaux,mptr)
!        mitot=54
!        print *,alloc(iadd(1,3,3))



!           Ntot=0
           call get_ordered_array()   
!           print *,ordered_mptr_array
           if (.not. allocated(field)) allocate(field(dim_state_p))
!           allocate(field(dim_state_p))
!           if (mype_world==0) print *,"allocated field"
!           call mpi_barrier(mpi_comm_world,mpierr)
!               print *,mptr_array
!               do ii = 1,size(ordered_mptr_array)
!                   mptr = ordered_mptr_array(ii)
!                   nx = node(ndihi,mptr) - node(ndilo, mptr) + 1
!                   ny = node(ndjhi,mptr) - node(ndjlo, mptr) + 1
!                   Ntot=Ntot+nx*ny
!               enddo
!               allocate(field(Ntot))  !need to figure out the dimension of new field
!              print *,Ntot
!          endif 
          print *,ordered_mptr_array
          Ntot_l=0
!          print *,alloc(8416)
        do ii = 1,size(ordered_mptr_array)
           mptr = ordered_mptr_array(ii)
!           if (mype_world==0) print *,ii,mptr
           nx = node(ndihi,mptr) - node(ndilo, mptr) + 1
           ny = node(ndjhi,mptr) - node(ndjlo, mptr) + 1
           mitot   = nx + 2*nghost
           mjtot   = ny + 2*nghost
           loc     = node(store1, mptr)
           locaux  = node(storeaux,mptr)
!           if (ii==11) print *, alloc(iadd(1,3,3)) 

!                  print *,"loc=",loc
!                  print *,"mptr =",mptr
!                  print *, "nx = ",nx
!                  print *, "ny = ",ny
!                  print *, "corners = ",corn1,corn2

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
                field(Ntot_l-nx*ny+i_mod+nx*(j_mod-1)) =&
                alloc(iadd(1,i_pkj,j_pkj)) +alloc(iaddaux(1,i_pkj, j_pkj))
            enddo
          enddo
        enddo
!        print *,'finish alloc to field'
!       deallocate(mptr_array)
!       deallocate(ordered_mptr_array)
!       print *,alloc(iadd(1,3,3))
!       print *,"finished"
!       print *,'field=',field(2500)

#endif
end subroutine alloc2field

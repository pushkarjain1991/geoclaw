!/*************************************************************************
!   > File Name: get_dim_state.f90
!    > Author: Shawn Lin
!    > Mail: lin_yuxiang@utexas.edu
!    > Created Time: Fri 05 Aug 2016 12:25:50 PM CDT
! ************************************************************************/

      subroutine get_dim_state()
          !use amr_module,only: node,ndihi,ndilo,ndjhi,ndjlo
          use amr_module,only: numcells
#ifdef USE_PDAF
!          use mod_model,only: mptr_array, ordered_mptr_array
          use mod_assimilation, only: dim_state_p
!          use mapdomain,only: get_ordered_array
!!          use PDAF_mod_filter,only:dim_p
!          use mod_assimilation,only: global_assimilation!include level 1 mesh or not

          implicit none
!          integer,intent(in):: level
!          integer :: Ntot
!          integer :: i,nx,ny,mptr
!          integer,allocatable :: mptr_array(:)
!          integer,allocatable :: ordered_mptr_array(:)
          
!          Ntot=0
!          call get_ordered_array()
!!          print *,ordered_mptr_array
!          do i=1,size(mptr_array)
!              mptr=mptr_array(i)
!              nx=node(ndihi,mptr)-node(ndilo,mptr)+1
!              ny=node(ndjhi,mptr)-node(ndjlo,mptr)+1
!              Ntot=Ntot+nx*ny
!          enddo
          print *, "numcells" 
          print *, numcells
          dim_state_p=numcells(1)
!          dim_p=Ntot
!          print *,Ntot
!          deallocate(mptr_array)
!          deallocate(ordered_mptr_array)
#endif
      end subroutine get_dim_state


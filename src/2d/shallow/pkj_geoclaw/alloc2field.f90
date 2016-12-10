!*************************************************************************
!   > File Name: alloc2field.f90
!    > Author: Shawn Lin
!    > Mail: lin_yuxiang@utexas.edu
!    > Created Time: Mon 01 Aug 2016 04:04:44 PM CDT
! ************************************************************************/

subroutine alloc2field(nvar,naux, analyze_water, get_total_height)
    use amr_module
    use geoclaw_module, only: dry_tolerance

#ifdef USE_PDAF
        use mod_model, only: field, wet_cell_index
        use mod_assimilation,only: dim_state_p
        use mod_parallel,only: mype_world, mpi_comm_world, mpierr
!        use sortarr
!        use mapdomain,only: get_ordered_array

        implicit none

        integer,intent(in)::nvar
        integer,intent(in)::naux
        logical,intent(in)::analyze_water
        logical, optional::get_total_height
        integer :: ii,mptr,Ntot,Ntot_l,j_pkj,i_pkj
        integer :: i_mod,j_mod,loc,locaux,mitot,mjtot,nx,ny 
        integer :: ivar,i,j,iaux
        integer :: iadd,iaddaux
        integer :: cell_cnt
        integer :: level
        integer :: ncells_temp
        real(kind=8),allocatable :: temp_field(:)
        integer,allocatable :: temp_wet_cell_index(:)
        integer :: wet_cell_cnt
        logical :: field_is_total_height

        iadd(ivar,i,j)=loc+ivar-1+nvar*((j-1)*mitot+i-1)
        iaddaux(iaux,i,j)=locaux+iaux-1+naux*(i-1)+naux*mitot*(j-1)

        !if(present(get_total_height)) then
        !  field_is_total_height = get_total_height
        !else
        !  field_is_total_height = .false.
        !endif

         !Count the total number of cells
         !The value will be used to allocate field
         ncells_temp=0
         level = 1
 66      if (level .gt. lfine) go to 91
            mptr = lstart(level)
 71         if (mptr .eq. 0) go to 81
                nx = node(ndihi,mptr) - node(ndilo, mptr) + 1
                ny = node(ndjhi,mptr) - node(ndjlo, mptr) + 1
                mitot   = nx + 2*nghost
                mjtot   = ny + 2*nghost
                loc     = node(store1, mptr)
                locaux  = node(storeaux,mptr)

                ncells_temp = ncells_temp + nx*ny
                mptr = node(levelptr, mptr)
             go to 71
 81      level = level + 1
         go to 66
 91      continue

           print *, "Running alloc2field ", mype_world
           if (allocated(field)) deallocate (field)
           allocate(field(ncells_temp))

           ! Values from alloc to field
           cell_cnt=1
           level = 1
 67        if (level .gt. lfine) go to 92
             mptr = lstart(level)
 72          if (mptr .eq. 0) go to 82
               nx = node(ndihi,mptr) - node(ndilo, mptr) + 1
               ny = node(ndjhi,mptr) - node(ndjlo, mptr) + 1
               mitot   = nx + 2*nghost
               mjtot   = ny + 2*nghost
               loc     = node(store1, mptr)
               locaux  = node(storeaux,mptr)

               do j_pkj = nghost+1, mjtot-nghost
                 do i_pkj = nghost+1, mitot-nghost

                   !field(cell_cnt) =&
                   !alloc(iadd(1,i_pkj,j_pkj))
                   field(cell_cnt) =&
                   alloc(iadd(1,i_pkj,j_pkj)) +alloc(iaddaux(1,i_pkj, j_pkj))
                   cell_cnt = cell_cnt + 1
                 enddo
               enddo
               mptr = node(levelptr, mptr)
               go to 72
 82          level = level + 1
           go to 67

 92        continue

#endif


end subroutine alloc2field

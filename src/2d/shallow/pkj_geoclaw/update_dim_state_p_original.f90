!*************************************************************************
!   > File Name: alloc2field.f90
!    > Author: Shawn Lin
!    > Mail: lin_yuxiang@utexas.edu
!    > Created Time: Mon 01 Aug 2016 04:04:44 PM CDT
! ************************************************************************/

subroutine update_dim_state_p()

        use mod_model, only: field
        use PDAF_mod_filter,only: dim_p, state, eofv, dim_ens
        use mod_parallel,only: mype_world, mpi_comm_world, mpierr
        use pdaf_mod_filtermpi, only: filterpe, dim_ens_l, &
                mpi_comm_null,comm_couple

        implicit none


         dim_p = size(field)
         if(filterpe) then
             
             if(allocated(state)) deallocate(state)
             allocate(state(dim_p))

             if(allocated(eofv)) deallocate(eofv)
             allocate(eofv(dim_p, dim_ens))
         else
             if (comm_couple /= mpi_comm_null) then
                 if(allocated(eofv)) deallocate(eofv)
                 allocate(eofv(dim_p, dim_ens_l))
             endif
         endif
         print *, "updated_dim_state_p = ", dim_p, mype_world

end subroutine update_dim_state_p

!*************************************************************************
!   > File Name: alloc2field.f90
!    > Author: Shawn Lin
!    > Mail: lin_yuxiang@utexas.edu
!    > Created Time: Mon 01 Aug 2016 04:04:44 PM CDT
! ************************************************************************/

subroutine update_dim_state_p()

        use mod_model, only: field
        USE PDAF_mod_filter, &
        ONLY: screen, dim_ens, rank, dim_p, dim_bias_p, &
              state, state_inc, eofU, eofV, &
              sens, bias, dim_lag, flag, filterstr

        use mod_assimilation, only: subtype, dim_state_p
        !use mod_parallel,only: mype_world, mpi_comm_world, mpierr
        !use pdaf_mod_filtermpi, only: filterpe, dim_ens_l, &
        !        mpi_comm_null,comm_couple

        implicit none


         dim_p = size(field)
         dim_state_p = size(field)

         if(allocated(state)) deallocate(state)
         if(allocated(state_inc)) deallocate(state_inc)
         if(allocated(eofU)) deallocate(eofU)
         if(allocated(eofV)) deallocate(eofV)
         if(allocated(sens)) deallocate(sens)
         if(allocated(bias)) deallocate(bias)

         flag=0
         call PDAF_alloc_filters(filterstr, subtype, flag)

!         if(filterpe) then
!             
!             if(allocated(state)) deallocate(state)
!             allocate(state(dim_p))
!
!             if(allocated(eofv)) deallocate(eofv)
!             allocate(eofv(dim_p, dim_ens))
!         else
!             if (comm_couple /= mpi_comm_null) then
!                 if(allocated(eofv)) deallocate(eofv)
!                 allocate(eofv(dim_p, dim_ens_l))
!             endif
!         endif
!         print *, "updated_dim_state_p = ", dim_p, mype_world

end subroutine update_dim_state_p

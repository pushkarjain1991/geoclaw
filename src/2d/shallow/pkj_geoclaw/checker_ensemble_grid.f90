subroutine check_ensemble_grid
      ! This routine is run by every ensemble member

      use common_level
      implicit none

      ! Check if same numbe of regions exist
      call num_region_count()

      ! Obtain the grids 
      call set_regions_assim_arrays()

      ! Create gatherv arrays
      call check_gatherv_array()


end subroutine check_ensemble_grid


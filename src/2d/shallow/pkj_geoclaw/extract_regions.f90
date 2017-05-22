subroutine extract_regions()
    use common_level
    implicit none

    !real(kind=8), intent(in) :: time
    !logical, intent(in) :: regrid_assim

    !if(regrid_assim) then
      
      !Get number of regions in processor pe
      call num_region_count()
      print *, "Rank ", mype_world, "has ", region_counter, "regions"

      print *, "Performing set_regions_arrays at Rank", mype_world
      call set_regions_assim_arrays()
      print *, "Done set_regions_arrays at Rank", mype_world

      call create_gatherv_array()
    !endif

    !update regions structure
    !call update_regions(time, regrid_assim)
    call update_regions()

end subroutine

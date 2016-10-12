!/*************************************************************************
!   > File Name: get_obs.f90
!    > Author: Shawn Lin
!    > Mail: lin_yuxiang@utexas.edu
!    > Created Time: Sun 14 Aug 2016 10:48:47 PM CDT
! ************************************************************************/

      subroutine get_obs(x,y,q,cnt)
#ifdef USE_PDAF
          use amr_module
          !use mod_assimilation,only: obs,obs_index, coords_obs
          use mod_assimilation, only: obs, obs_index
          use mod_parallel, only: mpi_comm_world, mpierr, mype_world

          implicit none
          real(kind=8),intent(in) :: x(:)
          real(kind=8),intent(in) :: y(:)
          real(kind=8),intent(in) :: q(:)
          integer :: cnt
!         Local variables          
          integer :: i,j,mptr,nx,ny,level,row,coln,cnt0,cnt1
          real(8) :: dx,dy,xlow,ylow,left,right,up,down
          real(8) :: temp_coord_obs_2d(2)
          real(8), allocatable :: q_local(:)
          logical :: obs_flag(size(x))
          logical :: obs_in_domain
          cnt = 0 

          print *, "Running get_obs"

         ! Count number of observations at finest mesh
         level = 1
 65      if (level .gt. lfine) go to 90
            mptr = lstart(level)
 70         if (mptr .eq. 0) go to 80
              nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
              ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
              xlow = rnode(cornxlo,mptr)
              ylow = rnode(cornylo,mptr)
              dx = hxposs(level)
              dy = hyposs(level)

              do j=1,nx*ny
                  row=(j-1)/nx+1
                  coln=mod(j-1,nx)+1
                  left=xlow+(coln-1)*dx;right=xlow+coln*dx
                  down=ylow+(row -1)*dy;  up =ylow+row*dy

                  obs_in_domain = any((x>=left).and.(x<right) .and. &
                  (y>=down) .and. (y<up))
                  if (obs_in_domain) then 
                      cnt=cnt+1
                  endif 
              enddo
  
            mptr = node(levelptr, mptr)
            go to 70
 80       level = level + 1
          go to 65
 90     continue
        print *, "Size of observations vector - ", cnt

        if (allocated(obs_index)) deallocate(obs_index)
        allocate(obs_index(cnt))
        if (allocated(obs)) deallocate(obs)
        allocate(obs(cnt))
          !if (.not. allocated(obs_index)) allocate(obs_index(cnt))
          !if (.not. allocated(obs)) allocate(obs(cnt))
          !if (.not. allocated(coords_obs)) allocate(coords_obs(2,cnt))

        cnt0=0
        cnt1=0
         
        level = 1
 66     if (level .gt. lfine) go to 91
          mptr = lstart(level)
 71       if (mptr .eq. 0) go to 81
            nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
            ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
            xlow = rnode(cornxlo,mptr)
            ylow = rnode(cornylo,mptr)
            dx=hxposs(level);dy=hyposs(level)

            do j=1,nx*ny
              cnt0=cnt0+1
              row=(j-1)/nx+1
              coln=mod(j-1,nx)+1
              left=xlow+(coln-1)*dx
              right=xlow+coln*dx
              down=ylow+(row -1)*dy  
              up =ylow+row*dy
              obs_flag=(x>=left).and.(x<right).and.(y>=down)&
              .and.(y<up)
              obs_in_domain = any(obs_flag)
              if (obs_in_domain) then 
                cnt1=cnt1+1
                obs_index(cnt1)=cnt0

                q_local=pack(q,obs_flag)
                obs(cnt1)=sum(q_local)/size(q_local)
                !obs(cnt1)=0.01


                !call ind1d_to_coord2d(cnt0,temp_coord_obs_2d)
                !coords_obs(1,cnt1)=temp_coord_obs_2d(1)
                !coords_obs(2,cnt1)=temp_coord_obs_2d(2)
                      
!               print *,left,right,down,up,size(q_local),&
!               obs(cnt1),coords_obs(1,cnt1),&
!               coords_obs(2,cnt1),cnt0,cnt1
              endif
            enddo

            mptr = node(levelptr, mptr)
            go to 71
 81       level = level + 1
        go to 66

 91     continue
        !print *, "obs after get_obs = ", obs
        !print *, mype_world, " obs_index = ", obs_index
        print *, "Finished running get_obs"
#endif
      end subroutine get_obs












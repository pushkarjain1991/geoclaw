!/*************************************************************************
!   > File Name: get_obs.f90
!    > Author: Shawn Lin
!    > Mail: lin_yuxiang@utexas.edu
!    > Created Time: Sun 14 Aug 2016 10:48:47 PM CDT
! ************************************************************************/

      subroutine get_obs(x,y,q,cnt)
          use amr_module
          use mod_assimilation, only: obs, obs_index, coords_obs
          use mod_parallel, only: mpi_comm_world, mpierr, mype_world
          implicit none

          real(kind=8),intent(in) :: x(:)
          real(kind=8),intent(in) :: y(:)
          real(kind=8),intent(in) :: q(:)
          integer, intent(inout) :: cnt
!         Local variables          
          integer :: mptr,nx,ny,level
          integer :: cnt0, cnt1
          real(8) :: dx,dy,xlow,ylow,left,right,up,down
          logical :: obs_flag(size(x))
          logical :: obs_in_domain
          integer :: i,j
           

          print *, "Running get_obs"

         cnt = 0
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

              do j=1,ny
                do i=1,nx
                  !row=(j-1)/nx+1
                  !coln=mod(j-1,nx)+1

                  left = xlow + (i-1)*dx
                  right = xlow + i*dx
                  down = ylow + (j -1)*dy
                  up = ylow + j*dy

                  obs_in_domain = any((x>=left).and.(x<right) .and. &
                  (y>=down) .and. (y<up))
                  if (obs_in_domain) then 
                      cnt=cnt+1
                  endif 
                enddo
              enddo
  
            mptr = node(levelptr, mptr)
            go to 70
 80       level = level + 1
          go to 65
 90     continue
        print *, "Size of observations vector - ", cnt

        if(allocated(obs_index)) deallocate(obs_index)
        allocate(obs_index(cnt))
        if(allocated(obs)) deallocate(obs)
        allocate(obs(cnt))
        if(allocated(coords_obs)) deallocate(coords_obs)
        allocate(coords_obs(2,cnt))

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
            dx=hxposs(level)
            dy=hyposs(level)

            do j=1,ny
              do i = 1,nx
                cnt0=cnt0+1
                !row=(j-1)/nx+1
                !coln=mod(j-1,nx)+1
                left = xlow + (i-1)*dx
                right = xlow + i*dx
                down = ylow + (j-1)*dy  
                up = ylow + j*dy
                obs_flag=(x>=left).and.(x<right).and.(y>=down)&
                .and.(y<up)
                obs_in_domain = any(obs_flag)
                if (obs_in_domain) then 
                  cnt1=cnt1+1
                  obs_index(cnt1)=cnt0

                  obs(cnt1)=sum(pack(q, obs_flag))/&
                  size(pack(q, obs_flag))
                  coords_obs(1,cnt1)=sum(pack(x, obs_flag))/&
                  size(pack(x, obs_flag))
                  coords_obs(2,cnt1)=sum(pack(y, obs_flag))/&
                  size(pack(y, obs_flag))
                endif
              enddo
            enddo

            mptr = node(levelptr, mptr)
            go to 71
 81       continue
        level = level + 1
        go to 66
 91     continue
        !print *, "obs after get_obs = ", obs
        !print *, mype_world, " obs_index = ", obs_index
        print *, "Finished running get_obs"
      end subroutine get_obs












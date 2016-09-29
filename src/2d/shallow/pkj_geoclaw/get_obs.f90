!/*************************************************************************
!   > File Name: get_obs.f90
!    > Author: Shawn Lin
!    > Mail: lin_yuxiang@utexas.edu
!    > Created Time: Sun 14 Aug 2016 10:48:47 PM CDT
! ************************************************************************/

      subroutine get_obs(x,y,q,cnt)
#ifdef USE_PDAF
          use amr_module
          use mod_model,only:ordered_mptr_array
          use mod_assimilation,only:obs,obs_index,coords_obs

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
          logical :: have_obs
          cnt=0 
          ! Count number of observations at finest mesh
          do i=1,size(ordered_mptr_array)
              mptr=ordered_mptr_array(i)
              level=node(nestlevel,mptr)
              nx=node(ndihi,mptr)-node(ndilo,mptr)+1
              ny=node(ndjhi,mptr)-node(ndjlo,mptr)+1
              xlow=rnode(cornxlo,mptr)
              ylow=rnode(cornylo,mptr)
              dx=hxposs(level);dy=hyposs(level)
              !if (level /= mxnest) cycle
              do j=1,nx*ny
                  row=(j-1)/nx+1
                  coln=mod(j-1,nx)+1
                  left=xlow+(coln-1)*dx;right=xlow+coln*dx
                  down=ylow+(row -1)*dy;  up =ylow+row*dy
!                  print *,level,row,coln,dx,dy
!                  print *,left,right,up,down
!                  print *,x,y,q

                  have_obs=any((x>=left).and.(x<right).and.(y>=down)&
                      .and.(y<up))
!                  have_obs=any(obs_flag)
!                  print *,have_obs                  
                  if (have_obs) then 
                      cnt=cnt+1
!                      print *,left,right,down,up 
                  endif 
              enddo
          enddo

          if (.not. allocated(obs_index)) allocate(obs_index(cnt))
          if (.not. allocated(obs)) allocate(obs(cnt))
          if (.not. allocated(coords_obs)) allocate(coords_obs(2,cnt))

          cnt0=0
          cnt1=0
          do i=1,size(ordered_mptr_array)
              mptr=ordered_mptr_array(i)
              level=node(nestlevel,mptr)
              nx=node(ndihi,mptr)-node(ndilo,mptr)+1
              ny=node(ndjhi,mptr)-node(ndjlo,mptr)+1
              xlow=rnode(cornxlo,mptr)
              ylow=rnode(cornylo,mptr)
              dx=hxposs(level);dy=hyposs(level)
              do j=1,nx*ny
                cnt0=cnt0+1
                !if (level==mxnest) then
                    row=(j-1)/nx+1
                    coln=mod(j-1,nx)+1
                    left=xlow+(coln-1)*dx;right=xlow+coln*dx
                    down=ylow+(row -1)*dy;  up =ylow+row*dy

                    obs_flag=(x>=left).and.(x<right).and.(y>=down)&
                      .and.(y<up)
                    have_obs=any(obs_flag)
                    if (have_obs) then 
                        cnt1=cnt1+1
                        obs_index(cnt1)=cnt0
                        q_local=pack(q,obs_flag)! for 1d array will pick
!                        out elements from mask obs_flag
                        obs(cnt1)=sum(q_local)/size(q_local)


                        call ind1d_to_coord2d(cnt0,temp_coord_obs_2d)
                        coords_obs(1,cnt1)=temp_coord_obs_2d(1)
                        coords_obs(2,cnt1)=temp_coord_obs_2d(2)
                        
!                        print *,left,right,down,up,size(q_local),&
!                            obs(cnt1),coords_obs(1,cnt1),&
!                            coords_obs(2,cnt1),cnt0,cnt1
                    endif
                !endif
              enddo
          enddo
!          print *,obs_index(100:130)
!          print *,coords_obs(1,100:130)
!          print *,coords_obs(2,100:130)
!          print *,obs(100:130)
#endif     
      end subroutine get_obs












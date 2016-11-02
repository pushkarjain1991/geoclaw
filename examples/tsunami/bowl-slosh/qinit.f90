! qinit routine for parabolic bowl problem, only single layer
subroutine qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

    use geoclaw_module, only: grav

#ifdef USE_PDAF
    use mod_parallel, only: mype_world
#endif

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: meqn,mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Parameters for problem
    real(kind=8), parameter :: a = 1.d0
    real(kind=8), parameter :: sigma = 0.5d0
    real(kind=8) :: h0 = 0.1d0

    ! Other storage
    integer :: i,j
    real(kind=8) :: omega,x,y,eta
    
#ifdef USE_PDAF
    logical :: first_pert = .true.
    integer :: n, clock
    real(kind=8),save :: r
    real(kind=8) :: r8_normal_ab
    real(kind=8), parameter :: pert_mu = 0.0D+00
    real(kind=8), parameter :: pert_sigma = 0.05D+00
    integer,allocatable :: seed(:)
    !character(len=1) :: ensstr

    if(first_pert .eqv. .true.) then
        !write(ensstr, '(i1)') mype_world+1
        call random_number(r)
        allocate(seed(n))
        !seed = clock + 37*(/(i-1, i=1,n)/)
        seed = 12345 + 10*mype_world

        r = r8_normal_ab(pert_mu, pert_sigma, seed)
        print *, "printing random number", r, mype_world

        first_pert = .false.
    endif
    
    h0 = h0+r
#endif

    omega = sqrt(2.d0 * grav * h0) / a
    
    do i=1-mbc,mx+mbc
        x = xlower + (i - 0.5d0)*dx
        do j=1-mbc,my+mbc
            y = ylower + (j - 0.5d0) * dx
            eta = sigma * h0 / a**2 * (2.d0 * x - sigma)
            
            q(1,i,j) = max(0.d0,eta - aux(1,i,j))
            q(2,i,j) = 0.d0
            q(3,i,j) = sigma * omega * q(1,i,j)
        enddo
    enddo
    
end subroutine qinit

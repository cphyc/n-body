program n_body
  use initial_conditions
  use constants
  use physics
  use integration
  
  implicit none

  real(kind = xp), dimension(:), allocatable :: m
  real(kind = xp), dimension(:, :), allocatable :: r, v, a
  integer :: npoints
  integer :: iter = 0
  integer :: maxiter = 10000

  
  !---------------------------------------------
  ! Read number of points
  !---------------------------------------------
  open(unit=10, file='initial_conditions.dat')
  call read_npoints(10, npoints)

  !---------------------------------------------
  ! Allocations
  !---------------------------------------------
  allocate(m(npoints))
  allocate(r(npoints, 3))
  allocate(v(npoints, 3))
  allocate(a(npoints, 3))


  !---------------------------------------------
  ! Read initial positions
  !---------------------------------------------
  call read_mxyz(10, m, r, npoints)
  close(unit=10)

  !---------------------------------------------
  ! Compute initial speeds
  !---------------------------------------------
  call initial_speeds(r, npoints, v)

  !---------------------------------------------
  ! Loop over time
  !---------------------------------------------
  do while (iter < maxiter)
     iter = iter + 1
     !-------------------------
     ! Compute acceleration
     !-------------------------

     call compute_force(m, r, npoints, a)

     call integration(v, a, dt/2)
     call integration(r, v, dt)
  end do

  !---------------------------------------------
  ! Deallocations
  !---------------------------------------------
  deallocate(r)
  deallocate(v)
  deallocate(a)

end program n_body

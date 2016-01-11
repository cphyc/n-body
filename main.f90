program n_body
  use io_tools
  use constants
  use physics
  
  implicit none

  real(kind = xp), dimension(:), allocatable :: m
  real(kind = xp), dimension(:, :), allocatable :: r, v, a
  integer :: iter = 0
  integer :: dump_freq = 1
  integer :: maxiter = 10000

  
  !---------------------------------------------
  ! Read number of points
  !---------------------------------------------
  open(unit=10, file='initial_conditions.dat')
  call read_npoints(10)

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
  call read_mpos(10, m, r)
  close(unit=10)

  !---------------------------------------------
  ! Compute initial speeds
  !---------------------------------------------
  call initial_speeds(r, v)

  !---------------------------------------------
  ! Loop over time
  !---------------------------------------------
  open(unit=5, file='output.dat')
  do while (iter < maxiter)
     iter = iter + 1
     !-------------------------
     ! Compute acceleration
     !-------------------------

     call compute_force(m, r, a)

     call integration(v, a, dt/2)
     call integration(r, v, dt)

     if (mod(dump_freq, iter) == 0) then
        call write_dump(5, r, v, npoints)
     end if
  end do
  close(unit=5)
  !---------------------------------------------
  ! Deallocations
  !---------------------------------------------
  deallocate(r)
  deallocate(v)
  deallocate(a)

end program n_body

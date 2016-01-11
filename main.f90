program n_body
  use io_tools
  use constants
  use physics
  
  implicit none

  real(kind = xp), dimension(:), allocatable :: m
  real(kind = xp), dimension(:, :), allocatable :: r, v, a
  real(kind = xp) :: t = 0._xp, dt, Ec, Ep, E
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
  ! Compute initial speeds, accelerations
  !---------------------------------------------
  call initial_speeds(r, v)
  call compute_force(m, r, a)

  !---------------------------------------------
  ! Compute time step
  !---------------------------------------------
  dt = 1.e-2_xp * minval(abs(norm2(v,2)/norm2(a,2)))
  !dt = 1.e-2_xp * minval(abs(r/v))
  print *, dt

  !---------------------------------------------
  ! Loop over time
  !---------------------------------------------
  open(unit=5, file='output.dat')
  do while (iter < maxiter)
     iter = iter + 1
     !-------------------------
     ! Compute acceleration
     !-------------------------
     call integrate(v, a, dt/2)
     call integrate(r, v, dt)
     call compute_force(m, r, a)
     call integrate(v, a, dt/2)
     call compute_energy(m, r, v, Ec, Ep, E)

     t = t + dt
     if (mod(iter, dump_freq) == 0) then
        print*, 'Dump!'
        call write_dump(5, iter, Ec, Ep, E, t, r, v)
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

program n_body
  use initial_conditions
  use constants
  use physics
  
  implicit none

  real(kind = xp), dimension(:), allocatable :: m, x, y, z, vx, vy, vz
  integer :: npoints
  
  !---------------------------------------------
  ! Read number of points
  !---------------------------------------------
  open(unit=10, file='initial_conditions.dat')
  call read_npoints(10, npoints)

  !---------------------------------------------
  ! Allocations
  !---------------------------------------------
  allocate(m(npoints))
  allocate(x(npoints))
  allocate(y(npoints))
  allocate(z(npoints))
  allocate(vx(npoints))
  allocate(vy(npoints))
  allocate(vz(npoints))

  !---------------------------------------------
  ! Read initial positions
  !---------------------------------------------
  call read_mxyz(10, m, x, y, z, npoints)
  close(unit=10)

  !---------------------------------------------
  ! Compute initial speeds
  !---------------------------------------------
  call initial_speeds(x, y, z, npoints, vx, vy, vz)

  !---------------------------------------------
  ! Deallocations
  !---------------------------------------------
  deallocate(x)
  deallocate(y)
  deallocate(z)
  deallocate(vx)
  deallocate(vy)
  deallocate(vz)

end program n_body

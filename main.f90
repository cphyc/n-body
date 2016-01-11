program n_body
  use initial_conditions
  use constants
  
  implicit none

  real(kind = xp), dimension(:), allocatable :: r, mu, phi
  integer :: npoints
  
  !---------------------------------------------
  ! Read initial conditions file
  !---------------------------------------------
  call read_initial_conditions('initial_conditions.dat', r, mu, phi, npoints)

  !---------------------------------------------
  ! 
  !---------------------------------------------

end program n_body

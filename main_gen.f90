program gen

  use constants

  integer(kind=4)   :: ios     ! I/O test variable
  integer           :: i       ! Points iteration variable
  character(len=50) :: line    ! String reading variable

  real(xp) :: x, y, z
  real(xp) :: radius
  real(xp) :: mass

  ! Open file
  open(unit=11, file="./input_parameters.dat", action="read", status="old", iostat=ios)
  if (ios /= 0) stop "OPENING input_parameters.dat ERROR"

  ! Read parameters
  read(11,fmt=*) line, npoints

  close(11)

  mass = 1._xp / npoints

  open(unit=12, file="initial_conditions.dat", action="write", iostat=ios)
  if (ios /= 0) stop "OPENING initial_conditions.dat ERROR"

  write(12, '(I8)') npoints

  do i = 1, npoints
     radius = 2
     do while(radius > 1)
        !       call init_random_seed()
        call random_number(x)
        call random_number(y)
        call random_number(z)
        radius = x**2 + y**2 + z**2
     end do

     write(12, '(4(e16.8e2))') x, y, z, mass
  end do

  close(12)

end program gen

module initial_conditions
  use constants

  implicit none

  private

  public :: read_initial_conditions
contains
  
  subroutine read_initial_conditions (filename, un)
    character, intent(in) :: filename
    integer, intent(in)   :: un
    integer               :: stat
    integer               :: npoints
    real(kind = xp), dimension(:), allocatable :: r, mu, phi
    integer :: i = 0

    ! FIXME: check read status
    open(unit=un, file=filename)

    read (un, *) npoints
    ! FIXME: check allocation status
    allocate (r(npoints), stat = stat)
    allocate (mu(npoints), stat = stat)
    allocate (phi(npoints), stat = stat)
    
    do i = 1, npoints
       read (un, *) r(i), mu(i), phi(i)
    end do 
    
  end subroutine read_initial_conditions

  
end module initial_conditions

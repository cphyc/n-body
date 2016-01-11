module initial_conditions
  use constants

  implicit none

  private

  public :: read_npoints, read_xyz
contains
  
  subroutine read_npoints (un, npoints)
    integer, intent(in) :: un
    integer, intent(out)                                    :: npoints

    read (un, *) npoints
  end subroutine read_npoints
  
  subroutine read_xyz(un, x, y, z, npoints)
    integer, intent(in) :: un
    integer, intent(in) :: npoints
    real(kind = xp), dimension(:), intent(out) :: x, y, z

    integer                                                 :: i = 0
    
    do i = 1, npoints
       read (un, *) x(i), y(i), z(i)
    end do 
  
  end subroutine read_xyz

  
end module initial_conditions

module initial_conditions
  use constants

  implicit none

  private

  public :: read_npoints, read_mxyz
contains
  
  subroutine read_npoints (un, npoints)
    integer, intent(in)  :: un
    integer, intent(out) :: npoints

    read (un, *) npoints
  end subroutine read_npoints
  
  subroutine read_mxyz(un, m, r, npoints)
    integer, intent(in)                        :: un
    integer, intent(in)                        :: npoints
    real(kind = xp), dimension(:), intent(out) :: m
    real(kind = xp), dimension(:, :), intent(out) :: r

    integer                                                 :: i = 0
    
    do i = 1, npoints
       read (un, *) r(i, 1), r(i, 2), r(i, 3), m(i)
    end do 
  
  end subroutine read_mxyz

  
end module initial_conditions

module physics
  use constants
  private

  public :: initial_speeds
contains
  subroutine initial_speeds (x, y, z, npoints, vx, vy, vz)
    real(kind = xp), intent(in), dimension(:) :: x, y, z
    integer, intent(in) :: npoints

    real(kind = xp), intent(out), dimension(:) :: vx, vy, vz

    integer :: i
    real(kind = xp) :: d ! distance to the z axis
    real(kind = xp) :: theta ! angle with x axis
    real(kind = xp) :: v ! speed

    do i = 1, npoints
       d = sqrt(x(i)**2 + y(i)**2)
       theta = atan2(y(i), x(i))
       v = 1._xp/(6._xp*pi)*d

       vx(i) = -v*sin(theta)
       vy(i) = v*cos(theta)
       vz(i) = 0._xp
       
    end do
    
  end subroutine initial_speeds

end module physics

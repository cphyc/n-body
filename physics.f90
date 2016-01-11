module physics
  use constants
  private

  public :: initial_speeds, compute_force
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

  subroutine compute_force (m, r, npoints, a)
    real(kind=xp), dimension(:), intent(in)     :: m
    real(kind=xp), dimension(:, :), intent(in)  :: r
    integer, intent(in) :: npoints
    
    real(kind=xp), dimension(:, :), intent(out) :: a
    
    real(kind=xp), dimension(3) :: tmp

    integer :: i, j

    ! FIXME: do optimisation using force is antisymmetric
    a = 0._xp
    do i = 1, npoints
       do j = 1, npoints
          if (j /= i) then
             tmp = G*m(j)*((r(i) - r(j))**2)**(1.5_xp)
             a(i) = a(i) + tmp*(r(j) - r(i))
          end if
       end do
    end do
    
  end subroutine compute_force


end module physics

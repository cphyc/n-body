module physics
  use constants
  private

  public :: initial_speeds, compute_force
contains
  subroutine initial_speeds (r, npoints, v)
    real(kind = xp), intent(in), dimension(:, :) :: r
    integer, intent(in) :: npoints

    real(kind = xp), intent(out), dimension(:, :) :: v

    integer :: i
    real(kind = xp) :: d ! distance to the z axis
    real(kind = xp) :: theta ! angle with x axis
    real(kind = xp) :: tmp

    do i = 1, npoints
       d = sqrt(r(i, 1)**2 + r(i, 2)**2)
       theta = atan2(r(i, 1), r(i, 2))
       tmp = 1._xp/(6._xp*pi)*d

       v(i, 1) = -tmp*sin(theta)
       v(i, 2) = tmp*cos(theta)
       v(i, 3) = 0._xp
       
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
             tmp = G*m(j)*((r(i, :) - r(j, :))**2)**(1.5_xp)
             a(i, :) = a(i, :) + tmp*(r(j, :) - r(i, :))
          end if
       end do
    end do
    
  end subroutine compute_force


end module physics

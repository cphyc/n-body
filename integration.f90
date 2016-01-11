module integration

use constants

   subroutine integrate(f, df, dt)
      implicit none

      real(xp), dimension(npoints), intent(inout) :: f, df
      real(xp),                     intent(in)    :: dt

      f = f + df * dt

   end subroutine

end module integration

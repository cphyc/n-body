module io_tools
   use constants

   implicit none

   private

   public :: read_mpos, write_dump, write_dump_headers
contains

   subroutine read_mpos(u, m, r)
      implicit none

      integer, intent(in) :: u

      real(kind = xp), intent(out) :: m(:)
      real(kind = xp), intent(out) :: r(:,:)

      integer :: i

      do i = 1, npoints
         read (u, *) r(1, i), r(2, i), r(3, i), m(i)
      end do

   end subroutine read_mpos

   subroutine write_dump_headers (u, ua)
      implicit none
      integer, intent(in) :: u, ua

      write (u, '(10(a16))') 'iter', 't', 'Ep', 'Ec', 'E'
      write (ua, '(10(a16))') 'x', 'y', 'z', 'vx', 'vy', 'vz'

   end subroutine write_dump_headers

   subroutine write_dump (ua, u, iter, Ec, Ep, E, t, r, v)
     implicit none

     real(kind=xp), intent(in) :: Ec, Ep, E, t
     real(kind=xp), intent(in) :: r(:,:)
     real(kind=xp), intent(in) :: v(:,:)
     integer,       intent(in) :: ua, u, iter

     integer :: i

     write (u, '(i16, 10(e16.8e2))') iter, t, Ep, Ec, E

     !write (ua, *) iter
     !do i = 1, npoints
     !   write (ua, '(10(e16.8e2))') r(1, i), r(2, i), r(3, i), v(1, i), v(2, i), v(3, i)
     !end do

   end subroutine write_dump

end module io_tools

module io_tools
   use constants

   implicit none

   private

   public :: read_init, write_dump, write_dump_headers

contains

   subroutine read_init(u, start, end, m, r, v) !FIXME: Need to be modified again to old code when flag_mpi=1 will be fixed
      implicit none

      integer, intent(in) :: u, start, end

      real(xp), intent(out) :: m(:)
      real(xp), intent(out) :: r(:,:)
      real(xp), intent(out) :: v(:,:)

      integer :: i

      ! Read until start - 1
      do i = 1, start - 1
         read (u, *)
      end do

      ! Read end-start
      do i = 1, end - start
         read (u, *) r(1, i), r(2, i), r(3, i), m(i), v(1, i), v(2, i), v(3, i)
      end do

   end subroutine read_init

   subroutine write_dump_headers (u, ua) ! FIXME: We will probably need to split this one when flag_mpi=1 will be fixed
      implicit none
      integer, intent(in) :: u, ua

      write (u, '(10(a16))') 'iter', 't', 'Ep', 'Ec'
      write (ua, '(10(a16))') 'x', 'y', 'z', 'vx', 'vy', 'vz'

   end subroutine write_dump_headers

   subroutine write_dump (ua, u, iter, Ec, Ep, t, r, v) ! FIXME: We will probably need to split this one when flag_mpi=1 will be fixed
      implicit none

      real(xp), intent(in) :: Ec, Ep, t
      real(xp), intent(in) :: r(:,:)
      real(xp), intent(in) :: v(:,:)
      integer,  intent(in) :: ua, u, iter

      integer :: i

      write (u, '(i16, 10(e16.8e2))') iter, t, Ep, Ec

      !write (ua, *) iter
      !do i = 1, npoints
      !   write (ua, '(10(e16.8e2))') r(1, i), r(2, i), r(3, i), v(1, i), v(2, i), v(3, i)
      !end do

   end subroutine write_dump

end module io_tools

module io_tools
   use constants

   implicit none

   private

   public :: read_npoints, read_mpos, write_dump, write_dump_headers
contains

   subroutine read_npoints (un)
      implicit none
      integer, intent(in)  :: un

      read (un, *) npoints

   end subroutine read_npoints

   subroutine read_mpos(un, m, r)
      implicit none

      integer, intent(in) :: un

      real(kind = xp), dimension(:),    intent(out) :: m
      real(kind = xp), dimension(:, :), intent(out) :: r

      integer :: i = 0

      do i = 1, npoints
         read (un, *) r(i, 1), r(i, 2), r(i, 3), m(i)
      end do

   end subroutine read_mpos

   subroutine write_dump_headers (un_int)
      implicit none
      integer, intent(in) :: un_int

      write (un_int, '(10(a16))') 'iter', 't', 'Ep', 'Ec', 'E'

   end subroutine write_dump_headers

   subroutine write_dump (un_all, un_int, iter, Ec, Ep, E, t, r, v)
      implicit none

      real(kind=xp),                 intent(in) :: Ec, Ep, E, t
      real(kind=xp), dimension(:,:), intent(in) :: r
      real(kind=xp), dimension(:,:), intent(in) :: v
      integer,                       intent(in) :: un_all, un_int, iter

      integer :: i

      write (un_int, '(i16, 10(e16.8e2))') iter, t, Ep, Ec, E

!      write (un_all, *) iter
!      write (un_all, '(10(a16))') 'x', 'y', 'z', 'vx', 'vy', 'vz'
!      do i = 1, npoints
!         write (un_all, '(10(e16.8e2))') r(i, 1), r(i, 2), r(i, 3), v(i, 1), v(i, 2), v(i, 3)
!      end do

  end subroutine write_dump

end module io_tools

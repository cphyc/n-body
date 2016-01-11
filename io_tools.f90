module io_tools
  use constants

  implicit none

  private

  public :: read_npoints, read_mpos, write_dump
contains
  
  subroutine read_npoints (un)
    integer, intent(in)  :: un

    read (un, *) npoints
  end subroutine read_npoints
  
  subroutine read_mpos(un, m, r)
    integer, intent(in)                           :: un
    real(kind = xp), dimension(:), intent(out)    :: m
    real(kind = xp), dimension(:, :), intent(out) :: r

    integer                                       :: i = 0
    
    do i = 1, npoints
       read (un, *) r(i, 1), r(i, 2), r(i, 3), m(i)
    end do 
  
  end subroutine read_mpos

  subroutine write_dump (un, r, v)
    real(kind=xp), intent(in), dimension(:,:) :: r
    real(kind=xp), intent(in), dimension(:,:) :: v
    integer, intent(in) :: un
    
    integer :: i

    do i = 1, npoints
       write (un,*) r(i, 1), r(i, 2), r(i, 3), v(i, 1), v(i, 2), v(i, 3)
    end do
    
  end subroutine write_dump

  
end module io_tools

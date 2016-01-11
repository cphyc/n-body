module io_tools
  use constants

  implicit none

  private

  public :: read_npoints, read_mpos, write_pos
contains
  
  subroutine read_npoints (un, npoints)
    integer, intent(in)  :: un
    integer, intent(out) :: npoints

    read (un, *) npoints
  end subroutine read_npoints
  
  subroutine read_mpos(un, m, r, npoints)
    integer, intent(in)                        :: un
    integer, intent(in)                        :: npoints
    real(kind = xp), dimension(:), intent(out) :: m
    real(kind = xp), dimension(:, :), intent(out) :: r

    integer                                                 :: i = 0
    
    do i = 1, npoints
       read (un, *) r(i, 1), r(i, 2), r(i, 3), m(i)
    end do 
  
  end subroutine read_mpos

  subroutine write_pos (un, r, npoints)
    real(kind=xp), intent(in), dimension(:,:) :: r
    integer, intent(in) :: un, npoints
    
    integer :: i

    do i = 1, npoints
       write (un,*) r(i, 1), r(i, 2), r(i, 3)
    end do
    
  end subroutine write_pos

  
end module io_tools

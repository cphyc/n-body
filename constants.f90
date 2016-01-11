module constants
  private
  integer, parameter :: xp = selected_real_kind(8)
  real(kind = xp)    :: pi = 3.14159265359_xp
  real(kind = xp)    :: G = 1._xp

  integer :: npoints

  public :: xp, pi, G, npoints
  
contains
  
end module constants

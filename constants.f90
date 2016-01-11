module constants
  private
  integer, parameter :: xp = selected_real_kind(8)
  real(kind = xp), parameter    :: pi = 3.14159265359_xp
  real(kind = xp), parameter    :: G = 1._xp

  real(kind = xp)    :: epsilon
  integer            :: npoints

  public :: xp, pi, G, epsilon, npoints
  
contains
  
end module constants

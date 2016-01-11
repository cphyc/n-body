module constants
  private
  integer, parameter :: xp = selected_real_kind(8)
  real(kind = xp)    :: pi = 3.14159265359_xp
  real(kind = xp)    :: G = 1._xp

  public :: xp, pi, G
  
contains
  
end module constants

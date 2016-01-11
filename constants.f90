module constants
  private
  integer, parameter :: xp = selected_real_kind(8)
  real(kind = xp)    :: pi = 3.14159265359_xp

  public :: xp, pi
  
contains
  
end module constants

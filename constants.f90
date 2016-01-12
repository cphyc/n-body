module constants
   implicit none

   private

   integer,         parameter    :: xp = selected_real_kind(8)
   real(kind = xp), parameter    :: pi = 4.0_xp*atan(1.0_xp)
   real(kind = xp), parameter    :: G  = 1._xp
   integer,         parameter    :: npoints = 1000

   integer         :: un, una  ! Unit numbers for file operation
   real(kind = xp) :: epsilon2

   public :: xp, pi, G, npoints, epsilon2, un, una

contains

end module constants

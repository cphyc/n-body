module constants
   implicit none

   private

   integer,         parameter    :: xp = selected_real_kind(8)
   real(kind = xp), parameter    :: pi = 4.0_xp*atan(1.0_xp)
   real(kind = xp), parameter    :: G  = 1._xp

   real(kind = xp) :: epsilon2
   integer         :: npoints

   public :: xp, pi, G, epsilon2, npoints

contains

end module constants

module constants
   implicit none

   private

   integer,         parameter :: xp = selected_real_kind(8)
   real(kind = xp), parameter :: pi = 4.0_xp*atan(1.0_xp)
   real(kind = xp), parameter :: G  = 1._xp
   integer,         parameter :: npoints = 1000
   real(kind = xp), parameter :: epsilon2 = npoints**(2._xp/3._xp) / 400._xp ! We take 1/20th of the initial caracteristic distance

   integer :: un, una  ! Unit numbers for file operation

   public :: xp, pi, G, npoints, epsilon2, un, una

contains

end module constants

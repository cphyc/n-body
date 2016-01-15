program gen

   use constants

   integer  :: ios     ! I/O test variable
   integer  :: i       ! Points iteration variable

   real(xp) :: x, y, z
   real(xp) :: radius
   real(xp) :: mass

   mass = 1._xp / npoints

   open(newunit=un, file="initial_conditions.dat", action="write", status="replace", iostat=ios)
   if (ios /= 0) stop "OPENING initial_conditions.dat ERROR"

   do i = 1, npoints

      radius = 2

      do while(radius > 1)
         call random_number(x)
         call random_number(y)
         call random_number(z)
         x = x*2._xp - 1_xp
         y = y*2._xp - 1_xp
         z = z*2._xp - 1_xp

         radius = x**2 + y**2 + z**2
      end do

      write(un, '(4(e16.8e2))') x, y, z, mass

   end do

   close(un)

end program gen

program n_body
   use io_tools
   use constants
   use physics

   implicit none

   real(kind = xp) :: m(npoints)        ! Masses of the particles
   real(kind = xp) :: r(npoints,3)      ! Positions of the particles (3-dim vectors)
   real(kind = xp) :: v(npoints,3)      ! Speeds of the particles (3-dim vectors)
   real(kind = xp) :: a(npoints,3)      ! Acceleration of the particles (3-dim vectors)
   real(kind = xp) :: t = 0._xp         ! Total time elapsed in the simulation
   real(kind = xp) :: dt                ! Timestep
   real(kind = xp) :: Ec                ! Total kinetic energy
   real(kind = xp) :: Ep                ! Total potential energy
   real(kind = xp) :: E                 ! Total energy
   integer         :: iter      = 0     ! Number of iterations ran
   integer         :: dump_freq = 10    ! Frequency at which the system is sampled
   integer         :: maxiter   = 500   ! Maximum number of iterations

   !---------------------------------------------
   ! Read initial positions
   !---------------------------------------------
   open(newunit=un, file='initial_conditions.dat', status="old")
   call read_mpos(un, m, r)
   close(un)

   !---------------------------------------------
   ! Compute initial speeds, accelerations, variables
   !---------------------------------------------
   call initial_speeds(r, v)
   call compute_force(m, r, a)
   call compute_initial_variables()

   !---------------------------------------------
   ! Compute time step as a fraction of dynamic time
   !---------------------------------------------
   dt = 1.e-3_xp
   print *, '# Simulation parameters'
   print *, '# dt', dt
   print *, '# npoints', npoints
   print *, '# maxiter', maxiter

   !---------------------------------------------
   ! Open files for output, add headers
   !---------------------------------------------
   open(newunit=una, file='output.dat', status="replace")
   open(newunit=un, file='output_int.dat', status="replace")

   call write_dump_headers(un)

   !---------------------------------------------
   ! Loop over time
   !---------------------------------------------
   do while (iter < maxiter)

      call integrate(v, a, dt/2)  ! Compute v(t+dt/2)
      call integrate(r, v, dt)    ! Compute r(t+dt)
      call compute_force_omp(m, r, a) ! Compute a(t+dt)
      call integrate(v, a, dt/2)  ! Compute v(t+dt)

      t = t + dt
      iter = iter + 1

      if (mod(iter, dump_freq) == 0) then
         print *, 'Dump', iter
         call compute_energy(m, r, v, Ec, Ep, E) ! Compute Ec, Ep, E at t+dt
         call write_dump(una, un, iter, Ec, Ep, E, t, r, v)
      end if

   end do

   !---------------------------------------------
   ! Close files
   !---------------------------------------------
   close(un)
   close(una)

end program n_body

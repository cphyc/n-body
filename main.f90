program n_body
   use io_tools
   use constants
   use physics

   implicit none

   real(kind = xp), dimension(:),    allocatable :: m ! Masses of the particles
   real(kind = xp), dimension(:, :), allocatable :: r ! Positions of the particles (3-dim vectors)
   real(kind = xp), dimension(:, :), allocatable :: v ! Speeds of the particles (3-dim vectors)
   real(kind = xp), dimension(:, :), allocatable :: a ! Acceleration of the particles (3-dim vectors)
   real(kind = xp) :: t = 0._xp ! Total time elapsed in the simulation
   real(kind = xp) :: dt        ! Timestep
   real(kind = xp) :: Ec        ! Total kinetic energy
   real(kind = xp) :: Ep        ! Total potential energy
   real(kind = xp) :: E         ! Total energy
   integer :: iter = 0          ! Number of iterations ran
   integer :: dump_freq = 10    ! Frequency at which the system is sampled
   integer :: maxiter = 10000   ! Maximum number of iterations


   !---------------------------------------------
   ! Read number of points
   !---------------------------------------------
   open(unit=10, file='initial_conditions.dat')
   call read_npoints(10)

   !---------------------------------------------
   ! Allocations
   !---------------------------------------------
   allocate(m(npoints))
   allocate(r(npoints, 3))
   allocate(v(npoints, 3))
   allocate(a(npoints, 3))

   !---------------------------------------------
   ! Read initial positions
   !---------------------------------------------
   call read_mpos(10, m, r)
   close(unit=10)

   !---------------------------------------------
   ! Compute initial speeds, accelerations, variables
   !---------------------------------------------
   call initial_speeds(r, v)
   call compute_force(m, r, a)
   call compute_initial_variables()

   !---------------------------------------------
   ! Compute time step
   !---------------------------------------------
   ! the timestep here is experimental and corresponds to 1% of the
   ! dynamic time 
   dt = 1.e-2_xp
   print *, dt

   !---------------------------------------------
   ! Loop over time
   !---------------------------------------------
   open(unit=11, file='output.dat')
   open(unit=12, file='output_int.dat')

   call write_dump_headers(12)

   do while (iter < maxiter)

      call integrate(v, a, dt/2)  ! Compute v(t+dt/2)
      call integrate(r, v, dt)    ! Compute r(t+dt)
      call compute_force(m, r, a) ! Compute a(t+dt)
      call integrate(v, a, dt/2)  ! Compute v(t+dt)

      t = t + dt
      iter = iter + 1

      if (mod(iter, dump_freq) == 0) then
         print*, 'Dump!'
         call compute_energy(m, r, v, Ec, Ep, E) ! Compute Ec, Ep, E at t+dt
         call write_dump(11, 12, iter, Ec, Ep, E, t, r, v)
      end if

   end do

   close(unit=11)
   close(unit=12)

   !---------------------------------------------
   ! Deallocations
   !---------------------------------------------
   deallocate(r)
   deallocate(v)
   deallocate(a)

end program n_body

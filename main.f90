program n_body
   use io_tools
   use constants
   use physics

   use mpi

   implicit none

   real(kind = xp) :: m(npoints)            ! Masses of the particles
   real(kind = xp) :: r(3,npoints)          ! Positions of the particles (3-dim vectors)
   real(kind = xp) :: v(3,npoints)          ! Speeds of the particles (3-dim vectors)
   real(kind = xp) :: a(3,npoints)          ! Acceleration of the particles (3-dim vectors)
   real(kind = xp) :: a_reduced(3,npoints)
   real(kind = xp) :: t = 0._xp             ! Total time elapsed in the simulation
   real(kind = xp) :: dt                    ! Timestep
   real(kind = xp) :: Ec                    ! Total kinetic energy
   real(kind = xp) :: Ep                    ! Total potential energy
   real(kind = xp) :: E                     ! Total energy
   integer         :: iter      = 0         ! Number of iterations ran
   integer         :: dump_freq = 10        ! Frequency at which the system is sampled
   real            :: maxtime = real(npoints)/10000 ! Maximum time (ad hoc)
   integer         :: maxiter

   integer         :: err, nprocs, rank = 0
   integer         :: MASTER = 0

   integer         :: istart, iend, domain_size

   !---------------------------------------------
   ! Initialize MPI
   !---------------------------------------------
   call mpi_init(err)
   call mpi_comm_size(mpi_comm_world, nprocs, err)
   call mpi_comm_rank(mpi_comm_world, rank, err)

   !---------------------------------------------
   ! Read initial positions
   !---------------------------------------------
   open(newunit=un, file='initial_conditions.dat', status="old")
   call read_mpos(un, m, r)
   close(un)

   !---------------------------------------------
   ! Compute number of domains
   !---------------------------------------------
   domain_size = ceiling(npoints * 1.0 / nprocs)

   !---------------------------------------------
   ! Compute initial speeds and accelerations
   !---------------------------------------------
   call initial_speeds(r, v)
   call compute_force(m, r, a)

   !---------------------------------------------
   ! Set time step and output parameters
   !---------------------------------------------
   dt = 1.e-2_xp
   maxiter = nint(maxtime / dt)
   if (rank == MASTER) then
      print *, '# Simulation parameters'
      print *, '# dt', dt
      print *, '# npoints', npoints
      print *, '# maxtime', maxtime
      print *, '# nprocs', nprocs

      !---------------------------------------------
      ! Open files for output, add headers
      !---------------------------------------------
      open(newunit=una, file='output.dat', status="replace")
      open(newunit=un, file='output_int.dat', status="replace")

      call write_dump_headers(un)
   end if

   !---------------------------------------------
   ! Loop over time
   !---------------------------------------------
   do while (iter < maxiter)

      call integrate(v, a, dt/2)  ! Compute v(t+dt/2)
      call integrate(r, v, dt)    ! Compute r(t+dt)
      !call compute_force(m, r, a) ! Compute a(t+dt)

      ! get the domain for integration from the number of nodes
      istart = domain_size*rank + 1
      iend   = min(istart + domain_size - 1, npoints)
      call compute_force_omp(m, r, istart, iend, a) ! Compute a(t+dt)

      !--------------------------------
      ! reduce accelerations
      !--------------------------------
      call mpi_allreduce(a, a_reduced, npoints, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, err)
      a = a_reduced

      call integrate(v, a, dt/2)  ! Compute v(t+dt)

      t = t + dt
      iter = iter + 1

      if (mod(iter, dump_freq) == 0 .and. rank == MASTER) then
         print *, 'Dump', iter, t
         !call compute_energy(m, r, v, Ec, Ep, E) ! Compute Ec, Ep, E at t+dt
         !call compute_energy_omp(m, r, v, Ec, Ep, E) ! Compute Ec, Ep, E at t+dt
         call compute_energy_omp_nn_1(m, r, v, Ec, Ep, E) ! Compute Ec, Ep, E at t+dt
         call write_dump(una, un, iter, Ec, Ep, E, t, r, v)
      end if

   end do

   !---------------------------------------------
   ! Close files
   !---------------------------------------------
   if (rank == MASTER) then
      close(un)
      close(una)
   end if

   !---------------------------------------------
   ! Stop MPI
   !---------------------------------------------
   call mpi_finalize(err)
end program n_body

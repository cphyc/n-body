program n_body
   use io_tools
   use constants
   use physics

   use mpi_tools

   implicit none

   ! Physics variables
   real(xp), allocatable :: m(:)    ! Masses of the particles
   real(xp), allocatable :: r(:, :) ! Positions of the particles (3-dim vectors)
   real(xp), allocatable :: v(:, :) ! Speeds of the particles (3-dim vectors)
   real(xp), allocatable :: a(:, :) ! Acceleration of the particles (3-dim vectors)
   real(xp)              :: Ec      ! Total kinetic energy
   real(xp)              :: Ep      ! Total potential energy

   ! Counters for the main loop
   real(xp):: t  = 0._xp ! Total time elapsed in the simulation
   integer :: iter = 0   ! Number of iterations ran

   ! Variables for MPI
   integer :: nprocs = 1 ! Must be a power of 2, default to 1 when we don’t run with MPI
   integer :: rank = 0   ! Rank of the tasks among the nprocs, start at 0 and default to this for running without MPI
   integer :: N          ! Length of each block/domain for the MPI nodes

   !---------------------------------------------
   ! Initialize MPI
   !---------------------------------------------
   call initialize_mpi_groups(nprocs, rank)

   !---------------------------------------------
   ! Perform some sanity checks
   !---------------------------------------------
   if (IAND(nprocs,nprocs-1) /= 0) stop 'E: The number of MPI_PROC is not a power of 2'
   if (npoints < nprocs) stop 'E: You need at least the same number of particles than MPI_PROC'
   if (flag_memory .and. flag_diag .and. npoints < 2*nprocs) stop 'E: You need at least twice more particles than MPI_PROC'

   !---------------------------------------------
   ! Compute size of domains, allocate space
   !---------------------------------------------
   if (flag_memory) then
      N = npoints / nprocs

      allocate(m(N))
      allocate(r(3, N))
      allocate(v(3, N))
      allocate(a(3, N))
   else
      if (flag_diag) then
         N = npoints / (2 * nprocs) ! We subdivide the domains to make couples !TODO: Check whether we could avoid this
      else
         N = npoints / nprocs
      end if

      allocate(m(npoints))
      allocate(r(3, npoints))
      allocate(v(3, npoints))
      allocate(a(3, npoints))
   end if

   !---------------------------------------------
   ! Read initial positions                      ! FIXME: Not compatible with all code modes, and not really good in MPI_low
   !---------------------------------------------
   if (flag_memory) then
      open(newunit=un, file='initial_conditions.dat', status="old")
      call read_init(un, rank*N + 1, (rank+1)*N, m, r, v)
      close(un)
   else
      open(newunit=un, file='initial_conditions.dat', status="old")
      call read_init(un, 1, npoints, m, r, v)
      close(un)
   end if

   !---------------------------------------------
   ! Print parameters
   !---------------------------------------------
   if (rank == MASTER) then
      print *, '# Simulation parameters'
      print *, '# dt            :', dt
      print *, '# npoints       :', npoints
      print *, '# maxtime       :', maxtime
      print *, '# maxiter       :', maxiter
      print *, '# nprocs        :', nprocs
      print *, '# N             :', N
      print *, '# diag          :', flag_diag
      print *, '# memory        :', flag_memory
      print *, '# memory_factor :', memory_factor

      !---------------------------------------------
      ! Open files for output, add headers
      !---------------------------------------------
      open(newunit=una, file='output.dat', status="replace")
      open(newunit=un, file='output_int.dat', status="replace")

      call write_dump_headers(un, una)
   end if

   !---------------------------------------------
   ! Compute initial energy and accelerations
   !---------------------------------------------
   call compute_force_wrap(N, rank, nprocs, m, r, a)
   call compute_energy_wrap(N, t, 0, rank, nprocs, m, r, v, Ec, Ep)

   !---------------------------------------------
   ! Loop over time
   !---------------------------------------------
   do while (iter < maxiter)

      call integrate(v, a, dt/2)                        ! Compute v(t+dt/2)
      call integrate(r, v, dt)                          ! Compute r(t+dt)

      call compute_force_wrap(N, rank, nprocs, m, r, a) ! Compute a(t+dt)

      call integrate(v, a, dt/2)                        ! Compute v(t+dt)

      t = t + dt
      iter = iter + 1

      if (mod(iter, dump_freq) == 0) then
         if (rank == MASTER) then
            print *, 'Dump', iter, t
         end if

         call compute_energy_wrap(N, t, iter, rank, nprocs, m, r, v, Ec, Ep)

      end if

   end do

   !---------------------------------------------
   ! Close files and free allocated
   !---------------------------------------------
   if (rank == MASTER) then
      close(un)
      close(una)
   end if

   deallocate(m)
   deallocate(r)
   deallocate(v)
   deallocate(a)

   !---------------------------------------------
   ! Stop MPI
   !---------------------------------------------
   call finalize_mpi_groups()

end program n_body

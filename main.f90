program n_body
   use io_tools
   use constants
   use physics

   use mpi

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
   integer :: err = 0    ! Store MPI error code resulting from MPI calls
   integer :: nprocs = 1 ! Must be a power of 2, default to 1 when we don’t run with MPI
   integer :: rank = 0   ! Rank of the tasks among the nprocs, start at 0 and default to this for running without MPI
   integer :: N          ! Length of each block/domain for the MPI nodes

   !---------------------------------------------
   ! Initialize MPI
   !---------------------------------------------
   call mpi_init(err)
   call mpi_comm_size(mpi_comm_world, nprocs, err)
   call mpi_comm_rank(mpi_comm_world, rank, err)

   !---------------------------------------------
   ! Compute size of domains
   !---------------------------------------------
   if (IAND(nprocs,nprocs-1) /= 0) stop 'E: The number of MPI_PROC is not a power of 2'

   if (flag_diag) then
      N = npoints / (2 * nprocs) ! We subdivide the domains to make couples
   else
      N = npoints / nprocs
   end if

   !---------------------------------------------
   ! Allocate space
   !---------------------------------------------
   select case (flag_mpi)
      case(0)
         allocate(m(npoints))
         allocate(r(3, npoints))
         allocate(v(3, npoints))
         allocate(a(3, npoints))
      case(1)
         if (nprocs == 1) stop 'E: you are trying to parallelize a task that requires 2 or more threads'
         N = npoints / nprocs ! FIXME: Currently, the flag_mpi=1 code doesn’t use subdomains
         allocate(m(N))
         allocate(r(3, N))
         allocate(v(3, N))
         allocate(a(3, N))
      case default
         stop "Unknown value of flag_mpi"
   end select

   !---------------------------------------------
   ! Read initial positions                      ! FIXME: Not compatible with all code modes, and not really good in MPI_low
   !---------------------------------------------
   select case (flag_mpi)
      case(0)
         open(newunit=un, file='initial_conditions.dat', status="old")
         call read_init(un, 1, npoints, m, r, v)
         close(un)
      case(1)
         open(newunit=un, file='initial_conditions.dat', status="old")
         call read_init(un, rank*N + 1, (rank+1)*N + 1, m, r, v)
         close(un)
      case default
         stop "Unknown value of flag_mpi"
   end select

   !---------------------------------------------
   ! Print parameters
   !---------------------------------------------
   if (rank == MASTER) then
      print *, '# Simulation parameters'
      print *, '# dt     :', dt
      print *, '# npoints:', npoints
      print *, '# maxtime:', maxtime
      print *, '# maxiter:', maxiter
      print *, '# nprocs :', nprocs
      print *, '# N      :', N
      print *, '# diag   :', flag_diag
      print *, '# mpi    :', flag_mpi

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
   call compute_energy_wrap(N, rank, nprocs, m, r, v, Ec, Ep)

   if (rank == MASTER) then !FIXME: In flag_mpi=1 mode, everyone should write its own data
      print *, 'Dump', iter, t
      call write_dump(una, un, iter, Ec, Ep, t, r, v)
   end if

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

         call compute_energy_wrap(N, rank, nprocs, m, r, v, Ec, Ep)

         if (rank == MASTER) then !FIXME: In flag_mpi=1 mode, everyone should write its own data
            print *, 'Dump', iter, t
            call write_dump(una, un, iter, Ec, Ep, t, r, v)
         end if

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
   call mpi_finalize(err)

end program n_body

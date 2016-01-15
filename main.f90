program n_body
   use io_tools
   use constants
   use physics

   use mpi

   implicit none

   real(kind = xp), allocatable :: m(:)                 ! Masses of the particles
   real(kind = xp), allocatable :: r(:, :)               ! Positions of the particles (3-dim vectors)
   real(kind = xp), allocatable :: v(:, :)               ! Speeds of the particles (3-dim vectors)
   real(kind = xp), allocatable :: a(:, : )               ! Acceleration of the particles (3-dim vectors)
   real(kind = xp) :: Ec                         ! Total kinetic energy
   real(kind = xp) :: Ep                         ! Total potential energy
   real(kind = xp) :: E                          ! Total energy

   real(kind = xp) :: t  = 0._xp                 ! Total time elapsed in the simulation
   integer         :: iter = 0                   ! Number of iterations ran

   integer         :: err = 0
   integer         :: rank = 0
   integer         :: nprocs = 1                 ! Assumed to be either 1 or pair
   integer         :: MASTER = 0
   integer         :: N                          ! Length of each block/domain

   !---------------------------------------------
   ! Initialize MPI
   !---------------------------------------------
   call mpi_init(err)
   call mpi_comm_size(mpi_comm_world, nprocs, err)
   call mpi_comm_rank(mpi_comm_world, rank, err)

   !---------------------------------------------
   ! Compute size of domains                     ! FIXME: Should not depend on flag, rest of code must adapt
   !---------------------------------------------
   select case (flag_compute_mpi)
      case(0)
         N = ceiling(npoints * 0.5 / nprocs)
      case(1)
         N = ceiling(npoints * 1.0 / nprocs)
      case default
         stop "Unknown value of flag_compute_mpi"
   end select

   !---------------------------------------------
   ! Allocate space                              ! FIXME: See above
   !---------------------------------------------
   select case (flag_compute_mpi)
      case(0)
         allocate(m(npoints))
         allocate(r(3, npoints))
         allocate(v(3, npoints))
         allocate(a(3, npoints))
      case(1)
         allocate(m(N))
         allocate(r(3, N))
         allocate(v(3, N))
         allocate(a(3, N))
      case default
         stop "Unknown value of flag_compute_mpi"
   end select

   !---------------------------------------------
   ! Read initial positions                      ! FIXME: Not compatible with all code modes, and not really good in MPI_low
   !---------------------------------------------
   select case (flag_compute_mpi)
      case(0)
         open(newunit=un, file='initial_conditions.dat', status="old")
         call read_mpos(un, 1, npoints, m, r)
         close(un)
      case(1)
         open(newunit=un, file='initial_conditions.dat', status="old")
         call read_mpos(un, rank*N + 1, (rank+1)*N + 1, m, r)
         close(un)
      case default
         stop "Unknown value of flag_compute_mpi"
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

      !---------------------------------------------
      ! Open files for output, add headers
      !---------------------------------------------
      open(newunit=una, file='output.dat', status="replace")
      open(newunit=un, file='output_int.dat', status="replace")

      call write_dump_headers(un, una)
   end if

   !---------------------------------------------
   ! Compute initial speeds and accelerations
   !---------------------------------------------
   call initial_speeds(r, v)
   call compute_force_wrap(N, rank, nprocs, m, r, a)

   !---------------------------------------------
   ! Loop over time
   !---------------------------------------------
   do while (iter < maxiter)

      call integrate(v, a, dt/2)  ! Compute v(t+dt/2)
      call integrate(r, v, dt)    ! Compute r(t+dt)

      call compute_force_wrap(N, rank, nprocs, m, r, a)

      call integrate(v, a, dt/2)  ! Compute v(t+dt)

      t = t + dt
      iter = iter + 1

      if (mod(iter, dump_freq) == 0 .and. rank == MASTER) then

         print *, 'Dump', iter, t

         select case (flag_compute_energy)
            case(0)
               call compute_energy(m, r, v, Ec, Ep, E)          ! Compute Ec, Ep, E at t+dt with sequential version
            case(1)
               call compute_energy_omp(m, r, v, Ec, Ep, E)      ! Compute Ec, Ep, E at t+dt with naive OpenMP version
            case(2)
               call compute_energy_omp_diag(m, r, v, Ec, Ep, E) ! Compute Ec, Ep, E at t+dt with fast OpenMP version
            case default
               stop "Unknown value of flag_compute_energy"
         end select

         call write_dump(una, un, iter, Ec, Ep, E, t, r, v)

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

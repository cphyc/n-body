program n_body
   use io_tools
   use constants
   use physics

   use mpi

   implicit none

   real(kind = xp) :: m(npoints)                ! Masses of the particles
   real(kind = xp) :: r(3,npoints)              ! Positions of the particles (3-dim vectors)
   real(kind = xp) :: v(3,npoints)              ! Speeds of the particles (3-dim vectors)
   real(kind = xp) :: a(3,npoints)              ! Acceleration of the particles (3-dim vectors)
   real(kind = xp) :: a_reduced(3,npoints)
   real(kind = xp) :: t = 0._xp                 ! Total time elapsed in the simulation
   real(kind = xp) :: dt                        ! Timestep
   real(kind = xp) :: Ec                        ! Total kinetic energy
   real(kind = xp) :: Ep                        ! Total potential energy
   real(kind = xp) :: E                         ! Total energy
   integer         :: iter      = 0             ! Number of iterations ran
   integer         :: dump_freq = 10            ! Frequency at which the system is sampled
   real            :: maxtime = npoints/128._xp ! Maximum time (ad hoc)
   integer         :: maxiter
   integer         :: flag_compute_force = 1    ! Choose force computation subroutine, 0=sequential, 1=omp, 2=omp_nn_1
   integer         :: flag_compute_energy = 2   ! Choose energy computation subroutine, 0=sequential, 1=omp, 2=omp_nn_1

   integer         :: err = 0
   integer         :: rank = 0
   integer         :: nprocs = 1
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
   istart = 1
   iend = npoints
   select case (flag_compute_force)
      case(0)
         call compute_force(m, r, istart, iend, a)
      case(1)
         call compute_force_omp(m, r, istart, iend, a)
      case(2)
         call compute_force_omp_nn_1(m, r, istart, iend, a)
      case default
         stop "Unknown value of flag_compute_force"
   end select

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

      call write_dump_headers(un,una)
   end if

   !---------------------------------------------
   ! Loop over time
   !---------------------------------------------
   do while (iter < maxiter)

      call integrate(v, a, dt/2)  ! Compute v(t+dt/2)
      call integrate(r, v, dt)    ! Compute r(t+dt)

      !--------------------------------
      ! Get the domain for integration from the number of nodes
      !--------------------------------
      istart = domain_size*rank + 1
      iend   = min(istart + domain_size - 1, npoints)
      select case (flag_compute_force)
         case(0)
            call compute_force(m, r, istart, iend, a)          ! Compute a(t+dt) with sequential version
         case(1)
            call compute_force_omp(m, r, istart, iend, a)      ! Compute a(t+dt) with naive OpenMP version
         case(2)
            call compute_force_omp_nn_1(m, r, istart, iend, a) ! Compute a(t+dt) with fast OpenMP version
         case default
            stop "Unknown value of flag_compute_force"
      end select

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

         select case (flag_compute_energy)
            case(0)
               call compute_energy(m, r, v, Ec, Ep, E)          ! Compute Ec, Ep, E at t+dt with sequential version
            case(1)
               call compute_energy_omp(m, r, v, Ec, Ep, E)      ! Compute Ec, Ep, E at t+dt with naive OpenMP version
            case(2)
               call compute_energy_omp_nn_1(m, r, v, Ec, Ep, E) ! Compute Ec, Ep, E at t+dt with fast OpenMP version
            case default
               stop "Unknown value of flag_compute_energy"
         end select

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

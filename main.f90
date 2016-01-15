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
   real(kind = xp), allocatable :: a_comm(:, :), a_right(:, :)
   real(kind = xp), allocatable :: r_i(:, :), r_np_i(:, :), r_right(:, :)
   real(kind = xp), allocatable :: a_reduced(:, :)
   real(kind = xp) :: Ec                         ! Total kinetic energy
   real(kind = xp) :: Ep                         ! Total potential energy
   real(kind = xp) :: E                          ! Total energy

   real(kind = xp) :: t  = 0._xp                 ! Total time elapsed in the simulation
   integer         :: iter = 0                   ! Number of iterations ran
   integer         :: i                          ! Counter for elements
   integer         :: N                          ! Length of each block

   integer         :: err = 0
   integer         :: rank = 0
   integer         :: nprocs = 1
   integer         :: MASTER = 0
   integer         :: stat(MPI_STATUS_SIZE)
   integer         :: domain_size

   integer         :: istart, iend, jstart, jend, istart2, iend2

   !---------------------------------------------
   ! Initialize MPI
   !---------------------------------------------
   call mpi_init(err)
   call mpi_comm_size(mpi_comm_world, nprocs, err)
   call mpi_comm_rank(mpi_comm_world, rank, err)

   !---------------------------------------------
   ! Compute number of domains
   !---------------------------------------------
   domain_size = ceiling(npoints * 1.0 / nprocs)
   N = ceiling(1._xp * npoints / nprocs)

   !---------------------------------------------
   ! Allocate space
   !---------------------------------------------
   allocate(m(N))
   allocate(r(3, N))
   allocate(v(3, N))
   allocate(a(3, N))
   allocate(a_comm(3, N))
   allocate(a_right(3, N))
   allocate(a_reduced(3, N))
   allocate(r_i(3, N))
   allocate(r_np_i(3, N))
   allocate(r_right(3, N))

   !---------------------------------------------
   ! Read initial positions
   !---------------------------------------------
   open(newunit=un, file='initial_conditions.dat', status="old")
   call read_mpos(un, rank*N + 1, (rank+1)*N + 1, m, r)
   close(un)

   !---------------------------------------------
   ! Compute initial speeds and accelerations
   !---------------------------------------------
   call initial_speeds(r, v)
   istart = 1
   iend = npoints
   select case (flag_compute_force)
      case(0)
         call compute_force(m, r, istart, iend, r, istart, iend, a)
      case(1)
         call compute_force_omp(m, r, istart, iend, r, istart, iend, a)
      case(2)
         call compute_force_omp_nn_1(m, r, istart, iend, r, istart, iend, a)
      case default
         stop "Unknown value of flag_compute_force"
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

      !---------------------------------------------
      ! Open files for output, add headers
      !---------------------------------------------
      open(newunit=una, file='output.dat', status="replace")
      open(newunit=un, file='output_int.dat', status="replace")

      call write_dump_headers(un,una)
   end if

   !---------------------------------------------
   ! Compute number of domains
   !---------------------------------------------
   domain_size = ceiling(npoints * 0.5 / nprocs)

   !---------------------------------------------
   ! Compute initial speeds and accelerations
   !---------------------------------------------
   call initial_speeds(r, v)
   istart = 1
   iend = npoints
   select case (flag_compute_force)
      case(0)
         call compute_force(m, r, istart, iend, r, istart, iend, a)
      case(1)
         call compute_force_diag(m, r, istart, iend/2, npoints, a)
      case(2)
         call compute_force_omp(m, r, istart, iend, r, istart, iend, a)
      case(3)
         call compute_force_omp_diag(m, r, istart, iend/2, npoints, a)
      case default
         stop "Unknown value of flag_compute_force"
   end select

   !---------------------------------------------
   ! Loop over time
   !---------------------------------------------
   do while (iter < maxiter)

      call integrate(v, a, dt/2)  ! Compute v(t+dt/2)
      call integrate(r, v, dt)    ! Compute r(t+dt)

      select case(flag_compute_mpi)
         case(0)
            !--------------------------------
            ! Get the domain for integration from the number of nodes
            !--------------------------------
            istart = domain_size*rank + 1
            iend   = min(istart + domain_size - 1, npoints/2)
            istart2 = istart*2-1
            iend2   = iend*2
            select case (flag_compute_force)
               case(0)
                  call compute_force(m, r, istart2, iend2, r, istart2, iend2, a)          ! Compute a(t+dt) with sequential version
               case(1)
                  call compute_force_diag(m, r, istart, iend, npoints, a) ! Compute a(t+dt) with fast OpenMP version
               case(2)
                  call compute_force_omp(m, r, istart2, iend2, r, istart2, iend2, a)      ! Compute a(t+dt) with naive OpenMP version
               case(3)
                  call compute_force_omp_diag(m, r, istart, iend, npoints, a) ! Compute a(t+dt) with fast OpenMP version
               case default
                  stop "Unknown value of flag_compute_force"
            end select
         case(1)
            !--------------------------------
            ! Scatter positions across all processes
            ! and compute interactions
            !--------------------------------
            a = 0._xp
            N = ceiling(1._xp * npoints / nprocs)
            do i = 1, nprocs / 2
               ! Get data from nprocs-i in i
               ! FIXME send r from right to left
               call mpi_sendrecv(r, 3*N, MPI_REAL_XP, nprocs-i+1, 0, &
                    r_right, 3*N, MPI_REAL_XP, i, 0, &
                    MPI_COMM_WORLD, stat, err)

               if (rank == i) then
                  r_i = r
               else if (rank == nprocs - i) then
                  r_np_i = r
               end if

               ! Broadcast i-th data
               call mpi_bcast(r_i, 3*N, MPI_REAL_XP, i, MPI_COMM_WORLD, err)
               ! Broadcast n-i-th data
               call mpi_bcast(r_np_i, 3*N, MPI_REAL_XP, nprocs - i, MPI_COMM_WORLD, err)

               ! If i == rank, compute both diagonals
               if (rank == i) then
                  r_right = r_i
                  call compute_force_diag(m, r, r, a)
                  call compute_force_diag(m, r_right, r_right, a_right)
                  a_comm = a

                  ! compute on right side but communicate nothing
               else if (rank < i) then
                  call compute_force(m, r_np_i, r_right, a_right)
                  a_comm = 0._xp

                  ! compute on left side and communicate interaction
               else
                  call compute_force(m, r, r_i, a_comm)
                  a = a - a_comm
               end if

               ! Receive interaction in i
               call mpi_reduce(a_comm, a, 3*N, MPI_REAL_XP, MPI_SUM, i, MPI_COMM_WORLD, err)
            end do

            !--------------------------------
            ! Gather missing data
            !--------------------------------
            do i = 1, nprocs / 2
               if (rank == i .or. rank == nprocs - i + 1) then
                  call mpi_sendrecv(a_right, 3*N, MPI_REAL_XP, i, 0, &
                       a_comm, 3*N, MPI_REAL_XP, nprocs - i + 1, 0, &
                       MPI_COMM_WORLD, stat, err)
                  if (rank == nprocs - i) then
                     a = a + a_comm
                  end if
               end if
            end do

            !--------------------------------
            ! Update positions
            !--------------------------------
            call integrate(v, a, dt/2)

         case default
            stop "Unknown value of flag_compute_mpi"
      end select
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIX BELOW
      !--------------------------------
      ! reduce accelerations
      !--------------------------------
      call mpi_allreduce(a, a_reduced, npoints*3, MPI_REAL_XP, MPI_SUM, MPI_COMM_WORLD, err)
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
               call compute_energy_omp_diag(m, r, v, Ec, Ep, E) ! Compute Ec, Ep, E at t+dt with fast OpenMP version
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

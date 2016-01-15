module physics
use constants

use mpi

implicit none

private

public :: initial_speeds, integrate, compute_force_wrap, compute_energy_wrap

contains

   subroutine initial_speeds (r, v)
      implicit none

      real(kind = xp), intent(in)  :: r(:,:)
      real(kind = xp), intent(out) :: v(:,:)

      v(1, :) =   r(2, :)
      v(2, :) = - r(1, :)
      v(3, :) = 0._xp

   end subroutine initial_speeds

   subroutine integrate(f, df, dt)
      implicit none

      real(xp), intent(in) :: dt
      real(xp), intent(inout) :: f(:,:), df(:,:)

      !$OMP PARALLEL
      !$OMP WORKSHARE
      f = f + df * dt
      !$OMP END WORKSHARE
      !$OMP END PARALLEL

   end subroutine integrate

   subroutine compute_force (m, r1, istart, iend, r2, jstart, jend, a)
      implicit none

      real(kind=xp), intent(in) :: m(:)
      real(kind=xp), intent(in) :: r1(:, :), r2(:, :)
      integer,       intent(in) :: istart, iend, jstart, jend

      real(kind=xp), intent(out) :: a(:, :)

      real(kind=xp) :: vec(3), tmp(3)
      integer :: i, j

      a = 0._xp

      ! No REDUCTION(+:a) because only one thread is modifying a(:,i)
      !$OMP PARALLEL PRIVATE(j, vec, tmp)
      !$OMP DO SCHEDULE(RUNTIME)
      do i = istart, iend

         do j = jstart, jend

            if (i /= j) then
               vec = r1(:, i) - r2(:, j)
               tmp = G / (norm2(vec)**2 + epsilon2)**1.5_xp

               a(:, i) = a(:, i) - tmp*m(j)*vec
            end if

         end do

      end do
      !$OMP END DO
      !$OMP END PARALLEL

   end subroutine compute_force

   subroutine compute_force_diag (m, r, istart, iend, length, a)
      implicit none

      real(kind=xp), intent(in) :: m(:)
      real(kind=xp), intent(in) :: r(:,:)
      integer,       intent(in) :: istart, iend
      integer,       intent(in) :: length

      real(kind=xp), intent(out) :: a(:,:)

      real(kind=xp) :: vec(3), tmp(3)
      integer :: i, j, k

      a = 0._xp

      !$OMP PARALLEL PRIVATE(j, k, vec, tmp) REDUCTION(+:a)
      !$OMP DO SCHEDULE(RUNTIME)
      do i = istart, iend

         do j = 1, i-1

            vec = r(:, i) - r(:, j)
            tmp = G / (norm2(vec)**2 + epsilon2)**1.5_xp

            a(:, i) = a(:, i) - tmp*m(j)*vec
            a(:, j) = a(:, j) + tmp*m(i)*vec

         end do

         k = length + 1 - i

         do j = 1, k-1

            vec = r(:, k) - r(:, j)
            tmp = G / (norm2(vec)**2 + epsilon2)**1.5_xp

            a(:, k) = a(:, k) - tmp*m(j)*vec
            a(:, j) = a(:, j) + tmp*m(k)*vec

         end do

      end do
      !$OMP END DO
      !$OMP END PARALLEL

   end subroutine compute_force_diag

   subroutine compute_force_wrap (N, rank, nprocs, m, r, a)
      implicit none

      integer,  intent(in) :: N
      integer,  intent(in) :: nprocs
      integer,  intent(in) :: rank
      real(xp), intent(in) :: m(:)
      real(xp), intent(in) :: r(:,:)

      real(xp), intent(out) :: a(:,:)

      real(xp), allocatable :: a_comm(:, :), a_right(:, :)
      real(xp), allocatable :: r_i(:, :), r_np_i(:, :), r_right(:, :)
      real(xp), allocatable :: a_reduced(:, :)

      integer         :: err = 0
      integer         :: stat(MPI_STATUS_SIZE)
      integer         :: i                          ! Counter for elements

      integer         :: istart, iend

      logical :: debug = .true.
!      integer         :: jstart, jend

      select case(flag_mpi)
         case(0)
            allocate(a_reduced(3, npoints))
            !--------------------------------
            ! Get the domain for integration from the number of nodes
            !--------------------------------
            istart = N*rank + 1
            iend   = min(istart + N - 1, npoints/2)
            if (flag_diag) then
               call compute_force_diag(m, r, istart, iend, npoints, a)                ! Compute a(t+dt) with fast version
            else
               call compute_force(m, r, istart*2-1, iend*2, r, 1, npoints, a)         ! Compute a(t+dt) with slow version
            end if
            !--------------------------------
            ! reduce accelerations
            !--------------------------------
            call mpi_allreduce(a, a_reduced, npoints*3, MPI_REAL_XP, MPI_SUM, MPI_COMM_WORLD, err)
            a = a_reduced
            deallocate(a_reduced)
         case(1)
            allocate(a_comm(3, N))
            allocate(a_right(3, N))
            allocate(r_i(3, N))
            allocate(r_np_i(3, N))
            allocate(r_right(3, N))

            !--------------------------------
            ! Scatter positions across all processes
            ! and compute interactions
            !--------------------------------
            a = 0._xp
            if (nprocs == 1) then
               stop 'E: you are trying to parallelize a task that requires 2 or more threads'
            else if (mod(nprocs, 2) /= 0) then
               print*, 'W: you are trying to parallelize a task that runs better with an even number of threads'
            end if

            do i = 0, nprocs / 2 - 1
               ! Get data from nprocs-i in i
               ! Note: use nprocs-i, because nprocs start at 0

               if (debug) print*, 'R:', rank, 'i:', i
               if (rank == nprocs-i-1) then
                  if (debug) print*, 'R:', rank, '(1/2) sends to:', i
                  call mpi_send(r, 3*N, MPI_REAL_XP, i, 0, MPI_COMM_WORLD, err)
               else if (rank == i) then
                  if (debug) print*, 'R:', rank, '(1/2) get from:', nprocs-i-1
                  call mpi_recv(r_right, 3*N, MPI_REAL_XP, nprocs-i-1, 0, MPI_COMM_WORLD, stat, err)
               end if

               if (debug) call mpi_barrier(MPI_COMM_WORLD, err)

               if (rank == i) then
                  r_i = r
               else if (rank == nprocs - i) then
                  r_np_i = r
               end if

               if (rank == i) then
                  if (debug) print*, 'R:', rank, '(1) broadcasts'
               end if
               ! Broadcast i-th data
               call mpi_bcast(r_i, 3*N, MPI_REAL_XP, i, MPI_COMM_WORLD, err)
               if (debug) call mpi_barrier(MPI_COMM_WORLD, err)

               if (rank == nprocs-i-1) then
                  if (debug) print*, 'R:', rank, '(2) broadcasts'
               end if
               ! Broadcast n-i-th data
               call mpi_bcast(r_np_i, 3*N, MPI_REAL_XP, nprocs-i-1, MPI_COMM_WORLD, err)
               if (debug) call mpi_barrier(MPI_COMM_WORLD, err)

               ! If i == rank, compute both diagonals
               if (rank == i) then
                  r_right = r_i
                  call compute_force_diag(m, r, 1, N, N, a)
                  call compute_force_diag(m, r_right, 1, N, N, a_right)
                  a_comm = a

                  ! compute on right side but communicate nothing
               else if (rank < i) then
                  call compute_force(m, r_np_i, 1, N, r_right, 1, N, a_right)
                  a_comm = 0._xp

                  ! compute on left side and communicate interaction
               else
                  call compute_force(m, r, 1, N, r_i, 1, N, a_comm)
                  a = a - a_comm
               end if

               if (rank == i) then
                  if (debug) print*, 'R:', rank, 'getting reduce'
               else
                  if (debug) print*, 'R:', rank, 'sending reduce'
               end if
               ! Receive interaction in i
               call mpi_reduce(a_comm, a, 3*N, MPI_REAL_XP, MPI_SUM, i, MPI_COMM_WORLD, err)
               if (debug) call mpi_barrier(MPI_COMM_WORLD, err)

               if (rank == 0) then
                  if (debug) print*, ''
               end if
            end do

            !--------------------------------
            ! Gather missing data
            !--------------------------------
            if (rank == 0 .and. debug) then
               print*, 'Exchanging data'
            end if
            do i = 0, nprocs / 2 - 1
               if (rank == i) then
                  if (debug) print*, 'R:', i, '(2/2) sends to:', nprocs-i-1
                  call mpi_send(a_right, 3*N, MPI_REAL_XP, nprocs-i-1, 0, MPI_COMM_WORLD, err)
               else if (rank == nprocs-i-1) then
                  if (debug) print*, 'R:', rank, '(2/2) gets from:', i
                  call mpi_recv(a_comm, 3*N, MPI_REAL_XP, i, 0, MPI_COMM_WORLD, stat, err)
                  a = a + a_comm
               end if
               if (debug) call mpi_barrier(MPI_COMM_WORLD, err)
            end do
            if (debug) call mpi_barrier(MPI_COMM_WORLD, err)
            if (rank == 0 .and. debug) then
               print*, 'Exchanged'
               print*, ''
            end if
            deallocate(a_comm)
            deallocate(a_right)
            deallocate(r_i)
            deallocate(r_np_i)
            deallocate(r_right)
         case default
            stop "Unknown value of flag_mpi"
      end select

   end subroutine compute_force_wrap

   subroutine compute_energy (m, r, v, Ec, Ep, E)
      implicit none

      real(kind=xp), intent(in) :: m(:)
      real(kind=xp), intent(in) :: r(:,:), v(:,:)

      real(kind=xp), intent(out) :: Ec, Ep, E

      integer :: i, j

      Ec = 0._xp
      Ep = 0._xp

      !$OMP PARALLEL REDUCTION(+:Ec,Ep) PRIVATE(j)
      !$OMP DO SCHEDULE(RUNTIME)
      do i = 1, npoints

         Ec = Ec + 0.5_xp * m(i) * norm2(v(:,i))**2

         do j = 1, npoints

            if (i /= j) then
               Ep = Ep - 0.5_xp * G * m(j) * m(i) / sqrt(norm2(r(:, i) - r(:, j))**2 + epsilon2)
            end if

         end do

      end do
      !$OMP END DO
      !$OMP END PARALLEL

      E = Ec + Ep

   end subroutine compute_energy

   subroutine compute_energy_diag (m, r, v, Ec, Ep, E)
      implicit none

      real(kind=xp), intent(in) :: m(:)
      real(kind=xp), intent(in) :: r(:,:), v(:,:)

      real(kind=xp), intent(out) :: Ec, Ep, E

      integer :: i, j, k

      Ec = 0._xp
      Ep = 0._xp

      !$OMP PARALLEL REDUCTION(+:Ec,Ep) PRIVATE(j,k)
      !$OMP DO SCHEDULE(RUNTIME)
      do i = 1, npoints/2

         Ec = Ec + 0.5_xp * m(i) * norm2(v(:,i))**2

         do j = 1, i-1

            Ep = Ep -  G * m(j) * m(i) / sqrt(norm2(r(:, i) - r(:, j))**2 + epsilon2)

         end do

         k = npoints + 1 - i
         Ec = Ec + 0.5_xp * m(k) * norm2(v(:,k))**2

         do j = 1, k-1

            Ep = Ep -  G * m(j) * m(k) / sqrt(norm2(r(:, k) - r(:, j))**2 + epsilon2)

         end do

      end do
      !$OMP END DO
      !$OMP END PARALLEL

      E = Ec + Ep

   end subroutine compute_energy_diag

   subroutine compute_energy_wrap(N, rank, nprocs, m, r, v, Ec, Ep, E)
      implicit none

      integer,  intent(in) :: N, rank, nprocs
      real(xp), intent(in) :: m(:)
      real(xp), intent(in) :: r(:,:), v(:,:)

      real(xp), intent(out) :: Ec, Ep, E

      real(xp), allocatable :: r_gathered(:, :), v_gathered(:, :)
      integer :: err


      select case(flag_mpi)
      case(0)
         if (rank == MASTER) then
            if (flag_diag) then
               call compute_energy_diag(m, r, v, Ec, Ep, E) ! Compute Ec, Ep, E at t+dt with fast version
            else
               call compute_energy(m, r, v, Ec, Ep, E)      ! Compute Ec, Ep, E at t+dt with slow version
            end if
         end if
      case(1)
         ! gather r and v on master node
         if (rank == MASTER) then
            allocate(r_gathered(3, npoints))
            allocate(v_gathered(3, npoints))
         end if

         call mpi_gather(r, 3*N, MPI_REAL_XP, r_gathered, 3*npoints, MPI_REAL_XP, 0, MPI_COMM_WORLD, err)
         call mpi_gather(v, 3*N, MPI_REAL_XP, v_gathered, 3*npoints, MPI_REAL_XP, 0, MPI_COMM_WORLD, err)

         if (rank == MASTER) then

            if (flag_diag) then
               call compute_energy_diag(m, r_gathered, v_gathered, Ec, Ep, E) ! Compute Ec, Ep, E at t+dt with fast version
            else
               call compute_energy(m, r_gathered, v_gathered, Ec, Ep, E)      ! Compute Ec, Ep, E at t+dt with slow version
            end if

            deallocate(r_gathered)
            deallocate(v_gathered)

         end if
      case default
         stop "Unknown value of flag_mpi"
      end select

   end subroutine compute_energy_wrap

end module physics

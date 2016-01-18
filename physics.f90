module physics
   use constants

   use mpi
   use mpi_tools
   use io_tools

   implicit none

   private

   public :: integrate, compute_force_wrap, compute_energy_wrap

contains

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

      real(xp), intent(in) :: m(:)
      real(xp), intent(in) :: r1(:, :), r2(:, :)
      integer,  intent(in) :: istart, iend, jstart, jend

      real(xp), intent(inout) :: a(:, :)

      real(xp) :: vec(3), tmp(3)
      integer  :: i, j

      ! No REDUCTION(+:a) because only one thread is modifying a(:,i)
      !$OMP PARALLEL PRIVATE(j, vec, tmp)
      !$OMP DO SCHEDULE(RUNTIME)
      do i = istart, iend

         do j = jstart, jend

            ! When i = j, this code returns zero if r1=r2, but is important otherwise
            vec = r1(:, i) - r2(:, j)
            tmp = vec * G / sqrt(sum(vec**2) + epsilon2)**3

            a(:, i) = a(:, i) - m(j)*tmp

         end do

      end do
      !$OMP END DO
      !$OMP END PARALLEL

   end subroutine compute_force

   subroutine compute_force_mpi (m1, r1, istart, iend, m2, r2, jstart, jend, a1, a2)
      implicit none

      real(xp), intent(in) :: m1(:), m2(:)
      real(xp), intent(in) :: r1(:, :), r2(:, :)
      integer,  intent(in) :: istart, iend, jstart, jend

      real(xp), intent(inout) :: a1(:, :), a2(:, :)

      real(xp) :: vec(3), tmp(3)
      integer  :: i, j

      ! No REDUCTION(+:a) because only one thread is modifying a(:,i)
      !$OMP PARALLEL PRIVATE(j, vec, tmp)
      !$OMP DO SCHEDULE(RUNTIME)
      do i = istart, iend

         do j = jstart, jend

            ! When i = j, this code returns zero if r1=r2, but is important otherwise
            vec = r1(:, i) - r2(:, j)
            tmp = vec * G / sqrt(sum(vec**2) + epsilon2)**3

            a1(:, i) = a1(:, i) - m2(j)*tmp
            a2(:, j) = a2(:, j) + m1(i)*tmp

         end do

      end do
      !$OMP END DO
      !$OMP END PARALLEL

   end subroutine compute_force_mpi

   subroutine compute_force_diag (m, r, istart, iend, length, a)
      implicit none

      real(xp), intent(in) :: m(:)
      real(xp), intent(in) :: r(:,:)
      integer,  intent(in) :: istart, iend
      integer,  intent(in) :: length

      real(xp), intent(inout) :: a(:,:)

      real(xp) :: vec(3), tmp(3)
      integer  :: i, j, k

      !$OMP PARALLEL PRIVATE(j, k, vec, tmp) REDUCTION(+:a)
      !$OMP DO SCHEDULE(RUNTIME)
      do i = istart, iend

         do j = 1, i-1

            vec = r(:, i) - r(:, j)
            tmp = vec * G / sqrt(sum(vec**2) + epsilon2)**3

            a(:, i) = a(:, i) - m(j)*tmp
            a(:, j) = a(:, j) + m(i)*tmp

         end do

         k = length + 1 - i

         do j = 1, k-1

            vec = r(:, k) - r(:, j)
            tmp = vec * G / sqrt(sum(vec**2) + epsilon2)**3

            a(:, k) = a(:, k) - m(j)*tmp
            a(:, j) = a(:, j) + m(k)*tmp

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

      real(xp), allocatable :: a_comm_i(:, :), a_comm_np_i(:, :), a_right(:, :)
      real(xp), allocatable :: r_i(:, :), r_np_i(:, :), r_right(:, :)
      real(xp), allocatable :: m_i(:), m_np_i(:), m_right(:)
      real(xp), allocatable :: a_reduced(:, :)

      integer :: err = 0
      integer :: stat(MPI_STATUS_SIZE)
      integer :: i, np_i               ! Counters for elements
      integer :: j, k, s               ! Tmp variables

      integer :: i_to_translate(1), i_translated(1)

      integer :: istart, iend

      a = 0._xp

      select case(flag_mpi)
         case(0)
            allocate(a_reduced(3, npoints))
            !--------------------------------
            ! Get the domain for integration from the number of nodes
            !--------------------------------
            istart = N*rank + 1
            iend   = istart + N - 1
            if (flag_diag) then
               call compute_force_diag(m, r, istart, iend, npoints, a)  ! Compute a(t+dt) with fast version
            else
               call compute_force(m, r, istart, iend, r, 1, npoints, a) ! Compute a(t+dt) with slow version
            end if
            !--------------------------------
            ! reduce accelerations
            !--------------------------------
            call mpi_allreduce(a, a_reduced, npoints*3, MPI_REAL_XP, MPI_SUM, MPI_COMM_WORLD, err)
            a = a_reduced
            deallocate(a_reduced)
         case(1)
            !--------------------------------
            ! Allocate and initialize
            !--------------------------------
            allocate(a_comm_i(3, N))
            allocate(a_comm_np_i(3, N))
            allocate(a_right(3, N))
            allocate(r_i(3, N))
            allocate(r_np_i(3, N))
            allocate(r_right(3, N))
            allocate(m_i(N))
            allocate(m_np_i(N))
            allocate(m_right(N))

            a_right = 0._xp

            !--------------------------------
            ! Scatter positions across all processes
            ! and compute interactions
            !--------------------------------

            do i = 0, nprocs / 2 - 1
               ! Get data from nprocs-i in i
               ! Note: use nprocs-i-1, because nprocs start at 0
               np_i = nprocs - i - 1

               if (rank == np_i) then !TODO : Check wether doing two sends this way is safe, MPI_TAG?
                  call mpi_send(r,       3*N, MPI_REAL_XP,    i, 0, MPI_COMM_WORLD, err)
                  call mpi_send(m,         N, MPI_REAL_XP,    i, 0, MPI_COMM_WORLD, err)
               else if (rank == i) then
                  call mpi_recv(r_right, 3*N, MPI_REAL_XP, np_i, 0, MPI_COMM_WORLD, stat, err)
                  call mpi_recv(m_right,   N, MPI_REAL_XP, np_i, 0, MPI_COMM_WORLD, stat, err)
               end if

               if (rank == i) then
                  r_i = r
                  m_i = m
               else if (rank == np_i) then
                  r_np_i = r
                  m_np_i = m
               end if

               if (communicate_right(rank, i)) then
                  i_to_translate(1) = i
                  call mpi_group_translate_ranks(wgroup, 1, i_to_translate, mpi_group_to_right(i), i_translated, err)

                  ! Broadcast i-th data to the right of i
                  call mpi_bcast(r_i, 3*N, MPI_REAL_XP, i_translated(1), mpi_comm_to_right(i), err)
                  call mpi_bcast(m_i,   N, MPI_REAL_XP, i_translated(1), mpi_comm_to_right(i), err)
               end if

               if (communicate_left(rank, i)) then
                  i_to_translate(1) = np_i
                  call mpi_group_translate_ranks(wgroup, 1, i_to_translate, mpi_group_to_left(i), i_translated, err)

                  ! Broadcast n-i-th data to the left of i (included)
                  call mpi_bcast(r_np_i, 3*N, MPI_REAL_XP, i_translated(1), mpi_comm_to_left(i), err)
                  call mpi_bcast(m_np_i,   N, MPI_REAL_XP, i_translated(1), mpi_comm_to_left(i), err)
               end if


               ! compute own interactions
               a_comm_i = 0._xp
               a_comm_np_i = 0._xp
               if (rank == i) then
                  a_right = 0._xp
                  call compute_force_diag(m,       r,       1, N/2, N, a)
                  call compute_force_diag(m_right, r_right, 1, N/2, N, a_right)
                  a_comm_i = a

                  ! compute interaction on right side
               else if (rank < i) then
                  call compute_force_mpi(m_np_i, r_np_i, 1, N, m_right, r_right, 1, N, a_comm_np_i, a_right)

                  ! compute interaction on left side
               else
                  call compute_force_mpi(m_i,    r_i,    1, N, m,       r,       1, N, a_comm_i,    a)

                  if (rank == np_i) then
                     a_comm_np_i = a
                  end if

               end if

               if (communicate_right(rank, i)) then
                  i_to_translate(1) = i
                  call mpi_group_translate_ranks(wgroup, 1, i_to_translate, mpi_group_to_right(i), i_translated, err)
                  ! Receive interaction in i
                  call mpi_reduce(a_comm_i, a, 3*N, MPI_REAL_XP, MPI_SUM, i_translated(1), mpi_comm_to_right(i), err)
               end if

               if (communicate_left(rank, i)) then
                  i_to_translate(1) = np_i
                  call mpi_group_translate_ranks(wgroup, 1, i_to_translate, mpi_group_to_left(i), i_translated, err)
                  ! Receive interaction in np-i-1
                  call mpi_reduce(a_comm_np_i, a, 3*N, MPI_REAL_XP, MPI_SUM, i_translated(1), mpi_comm_to_left(i), err)
               end if

            end do

            !--------------------------------
            ! Gather missing data
            !--------------------------------
            do i = 0, nprocs / 2 - 1
               np_i = nprocs - i - 1
               if (rank == i) then
                  call mpi_send(a_right,     3*N, MPI_REAL_XP, np_i, 0, MPI_COMM_WORLD, err)
               else if (rank == np_i) then
                  call mpi_recv(a_comm_np_i, 3*N, MPI_REAL_XP,    i, 0, MPI_COMM_WORLD, stat, err)
                  a = a + a_comm_np_i
               end if
            end do

            deallocate(a_comm_i)
            deallocate(a_comm_np_i)
            deallocate(a_right)
            deallocate(r_i)
            deallocate(r_np_i)
            deallocate(r_right)
            deallocate(m_i)
            deallocate(m_np_i)
            deallocate(m_right)

         case(2) !FIXME: Does not work with different masses
            allocate(a_comm_i(3, N))
            allocate(r_i(3, N))
            ! Iterate over each proc. sending each time the local position

            do i = 0, nprocs - 1

               if (i == rank) then
                  r_i = r
               end if

               call mpi_bcast(r_i, 3*N, MPI_REAL_XP, i, MPI_COMM_WORLD, err)

               a_comm_i = 0._xp
               call compute_force(m, r_i, 1, N, r, 1, N, a_comm_i)

               call mpi_reduce(a_comm_i, a, 3*N, MPI_REAL_XP, MPI_SUM, i, MPI_COMM_WORLD, err)

            end do

            deallocate(a_comm_i)
            deallocate(r_i)
         case(3)
            s = N / memory_factor
            if (s*memory_factor /= N) then
               print*, 'E: Memory_factor', memory_factor
               print*, 'E: Number of particles per process', N
               stop 'E: The memory factor should be a divider of the number of particles per process'
            end if

            allocate(a_comm_i(3, s))
            allocate(r_i(3, s))
            allocate(m_i(s))

            do i = 0, nprocs-1

               do j = 1, memory_factor

                  if (i == rank) then
                     r_i = r(:, (j-1)*s+1:j*s)
                  end if

                  call mpi_bcast(r_i, 3*s, MPI_REAL_XP, i, MPI_COMM_WORLD, err)
                  a_comm_i = 0._xp

                  call compute_force(m, r_i, 1, s, r, 1, N, a_comm_i)

                  call mpi_reduce(a_comm_i, a(:, (j-1)*s+1:j*s), 3*s, MPI_REAL_XP, MPI_SUM, i, MPI_COMM_WORLD, err)

               end do
            end do

            deallocate(r_i)
            deallocate(m_i)
            deallocate(a_comm_i)
         case(4)
            !--------------------------------
            ! Allocate and initialize
            !--------------------------------
            allocate(a_comm_i(3, N/2))
            allocate(a_comm_np_i(3, N))
            allocate(r_i(3, N/2))
            allocate(r_np_i(3, N))
            allocate(m_i(N/2))
            allocate(m_np_i(N))

            !--------------------------------
            ! Scatter positions across all processes
            ! and compute interactions
            !--------------------------------

            do i = 0, nprocs - 1
               ! Get data from nprocs-i in i

               if (rank == i) then
                  r_i = r(:, 1:N/2)
                  m_i = m(1:N/2)
                  r_np_i = r
                  m_np_i = m
               end if

               if (communicate_right(rank, i)) then
                  i_to_translate(1) = i
                  call mpi_group_translate_ranks(wgroup, 1, i_to_translate, mpi_group_to_right(i), i_translated, err)

                  ! Broadcast i-th data to the right of i
                  call mpi_bcast(r_i, 3*N/2, MPI_REAL_XP, i_translated(1), mpi_comm_to_right(i), err)
           print *, rank, "I was here"
                  call mpi_bcast(m_i,   N/2, MPI_REAL_XP, i_translated(1), mpi_comm_to_right(i), err)
               end if

               if (communicate_left(rank, i)) then
                  i_to_translate(1) = i
                  call mpi_group_translate_ranks(wgroup, 1, i_to_translate, mpi_group_to_left(i), i_translated, err)

                  ! Broadcast n-i-th data to the left of i (included)
                  call mpi_bcast(r_np_i, 3*N, MPI_REAL_XP, i_translated(1), mpi_comm_to_left(i), err)
                  call mpi_bcast(m_np_i,   N, MPI_REAL_XP, i_translated(1), mpi_comm_to_left(i), err)
               end if

               ! compute own interactions
               a_comm_i = 0._xp
               a_comm_np_i = 0._xp
               if (rank == i) then
                  call compute_force_diag(m,       r,       1, N/2, N, a)
                  a_comm_i = a(:,1:N/2)

                  ! compute interaction on right side
               else if (rank < i) then
                  call compute_force_mpi(m_np_i, r_np_i, 1, N, m(N/2+1:N), r(:,N/2+1:N), 1, N/2, a_comm_np_i, a(:,N/2+1:N))

                  ! compute interaction on left side
               else
                  call compute_force_mpi(m_i,    r_i,    1, N/2, m,       r,       1, N, a_comm_i,    a)

               end if

               if (communicate_right(rank, i)) then
                  i_to_translate(1) = i
                  call mpi_group_translate_ranks(wgroup, 1, i_to_translate, mpi_group_to_right(i), i_translated, err)
                  ! Receive interaction in i
                  call mpi_reduce(a_comm_i, a, 3*N/2, MPI_REAL_XP, MPI_SUM, i_translated(1), mpi_comm_to_right(i), err)
               end if

               if (rank == i) a_comm_np_i = a

               if (communicate_left(rank, i)) then
                  i_to_translate(1) = i
                  call mpi_group_translate_ranks(wgroup, 1, i_to_translate, mpi_group_to_left(i), i_translated, err)
                  ! Receive interaction in np-i-1
                  call mpi_reduce(a_comm_np_i, a, 3*N, MPI_REAL_XP, MPI_SUM, i_translated(1), mpi_comm_to_left(i), err)
               end if

            end do

            deallocate(a_comm_i)
            deallocate(a_comm_np_i)
            deallocate(r_i)
            deallocate(r_np_i)
            deallocate(m_i)
            deallocate(m_np_i)

         case default
            stop "Unknown value of flag_mpi"
      end select

   end subroutine compute_force_wrap

   subroutine compute_energy (m, r, v, istart, iend, Ec, Ep)
      implicit none

      real(xp), intent(in) :: m(:)
      real(xp), intent(in) :: r(:,:), v(:,:)
      integer,  intent(in) :: istart, iend

      real(xp), intent(out) :: Ec, Ep

      integer :: i, j

      !$OMP PARALLEL REDUCTION(+:Ec,Ep) PRIVATE(j)
      !$OMP DO SCHEDULE(RUNTIME)
      do i = istart, iend

         Ec = Ec + 0.5_xp * m(i) * sum(v(:,i)**2)

         do j = 1, npoints

            ! TODO: Check this condition when using mpi_flag=1 and distributed energy computation
            if (i /= j) then
               Ep = Ep - 0.5_xp * G * m(j) * m(i) / sqrt(sum((r(:, i) - r(:, j))**2) + epsilon2)
            end if

         end do

      end do
      !$OMP END DO
      !$OMP END PARALLEL

   end subroutine compute_energy

   subroutine compute_energy_diag (m, r, v, istart, iend, length, Ec, Ep)
      implicit none

      real(xp), intent(in) :: m(:)
      real(xp), intent(in) :: r(:,:), v(:,:)
      integer,  intent(in) :: istart, iend
      integer,  intent(in) :: length

      real(xp), intent(out) :: Ec, Ep

      integer :: i, j, k

      !$OMP PARALLEL REDUCTION(+:Ec,Ep) PRIVATE(j,k)
      !$OMP DO SCHEDULE(RUNTIME)
      do i = istart, iend

         Ec = Ec + 0.5_xp * m(i) * sum(v(:,i)**2)

         do j = 1, i-1

            Ep = Ep - G * m(j) * m(i) / sqrt(sum((r(:, i) - r(:, j))**2) + epsilon2)

         end do

         k = length + 1 - i
         Ec = Ec + 0.5_xp * m(k) * sum(v(:,k)**2)

         do j = 1, k-1

            Ep = Ep - G * m(j) * m(k) / sqrt(sum((r(:, k) - r(:, j))**2) + epsilon2)

         end do

      end do
      !$OMP END DO
      !$OMP END PARALLEL

   end subroutine compute_energy_diag

   subroutine compute_energy_wrap(N, t, iter, rank, nprocs, m, r, v, Ec, Ep)
      implicit none

      integer,  intent(in) :: N
      real(xp), intent(in) :: t
      integer,  intent(in) :: iter, rank, nprocs
      real(xp), intent(in) :: m(npoints)
      real(xp), intent(in) :: r(:, :), v(:, :)

      real(xp), intent(out) :: Ec, Ep

      integer  :: istart, iend
      real(xp) :: Ec_loc, Ep_loc

      real(xp) :: r_gathered(3, npoints), v_gathered(3, npoints), m_gathered(npoints)

      integer :: err

      Ec = 0._xp
      Ep = 0._xp

      select case(flag_mpi)
      case(0)
         !--------------------------------
         ! Get the domain for integration from the number of nodes
         !--------------------------------
         istart = N*rank + 1
         iend   = istart + N - 1
         Ec_loc = 0._xp
         Ep_loc = 0._xp
         if (flag_diag) then
            call compute_energy_diag(m, r, v, istart, iend, npoints, Ec_loc, Ep_loc) ! Compute Ec & Ep at t+dt with fast version
         else
            call compute_energy(m, r, v, istart, iend, Ec_loc, Ep_loc)      ! Compute Ec & Ep at t+dt with slow version
         end if
         !--------------------------------
         ! reduce accelerations
         !--------------------------------
         call mpi_reduce(Ec_loc, Ec, 1, MPI_REAL_XP, MPI_SUM, MASTER, MPI_COMM_WORLD, err)
         call mpi_reduce(Ep_loc, Ep, 1, MPI_REAL_XP, MPI_SUM, MASTER, MPI_COMM_WORLD, err)

         r_gathered = r
         v_gathered = v
      case(1, 2, 3, 4)
         ! gather r and v on master node

         call mpi_gather(r, 3*N, MPI_REAL_XP, r_gathered, 3*N, MPI_REAL_XP, 0, MPI_COMM_WORLD, err)
         call mpi_gather(v, 3*N, MPI_REAL_XP, v_gathered, 3*N, MPI_REAL_XP, 0, MPI_COMM_WORLD, err)
         call mpi_gather(m,   N, MPI_REAL_XP, m_gathered,   N, MPI_REAL_XP, 0, MPI_COMM_WORLD, err)

         if (rank == MASTER) then

            if (flag_diag) then
               call compute_energy_diag(m_gathered, r_gathered, v_gathered, 1, npoints/2, npoints, Ec, Ep) ! Compute Ec & Ep at t+dt with fast version
            else
               ! TODO: This would be broken when computing using subdomains, because Ec must not be calculated elsewhere than the
               ! real diag
               call compute_energy(m_gathered, r_gathered, v_gathered, 1, npoints, Ec, Ep)      ! Compute Ec & Ep at t+dt with slow version
            end if

         end if
      case default
         stop "Unknown value of flag_mpi"
      end select

      if (rank == MASTER) then !FIXME: In flag_mpi=1 mode, everyone should write its own data
         call write_dump(iter, Ec, Ep, t, r_gathered, v_gathered)
      end if


   end subroutine compute_energy_wrap

end module physics

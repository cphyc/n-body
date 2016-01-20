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

      ! No REDUCTION(+:a1) because only one thread is modifying a(:,i)
      !$OMP PARALLEL PRIVATE(j, vec, tmp) REDUCTION(+:a2)
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

      real(xp), allocatable :: a_comm_i(:, :), a_comm_np_i(:, :)
      real(xp), allocatable :: r_i(:, :), r_np_i(:, :)
      real(xp), allocatable :: m_i(:), m_np_i(:)
      real(xp), allocatable :: a_reduced(:, :)

      integer :: err = 0

      integer :: i      ! Counters for elements
      integer :: j, s   ! Tmp variables

      integer :: i_to_translate(1), i_translated(1)

      logical :: communicate_right, communicate_left

      a = 0._xp

      if (flag_memory) then

         !---------------------------------------------------------------------------
         ! Compute the force using only own particles and communicating each other
         !---------------------------------------------------------------------------
         if (flag_diag) then
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
            ! Loop over all processes
            !--------------------------------

            do i = 0, nprocs - 1
               ! Get data from nprocs-i in i

               communicate_right = rank >= i
               communicate_left  = rank <= i

               if (rank == i) then
                  r_i = r(:, 1:N/2)
                  m_i = m(1:N/2)
                  r_np_i = r
                  m_np_i = m
               end if

               !--------------------------------
               ! Scatter positions and masses
               !--------------------------------
               if (communicate_right) then
                  i_to_translate(1) = i
                  call mpi_group_translate_ranks(wgroup, 1, i_to_translate, mpi_group_to_right(i), i_translated, err)

                  ! Broadcast i-th data to the right of i
                  call mpi_bcast(r_i, 3*N/2, MPI_REAL_XP, i_translated(1), mpi_comm_to_right(i), err)
                  call mpi_bcast(m_i,   N/2, MPI_REAL_XP, i_translated(1), mpi_comm_to_right(i), err)
               end if

               if (communicate_left) then
                  i_to_translate(1) = i
                  call mpi_group_translate_ranks(wgroup, 1, i_to_translate, mpi_group_to_left(i), i_translated, err)

                  ! Broadcast n-i-th data to the left of i (included)
                  call mpi_bcast(r_np_i, 3*N, MPI_REAL_XP, i_translated(1), mpi_comm_to_left(i), err)
                  call mpi_bcast(m_np_i,   N, MPI_REAL_XP, i_translated(1), mpi_comm_to_left(i), err)
               end if

               !--------------------------------
               ! Compute interactions
               !--------------------------------
               a_comm_i = 0._xp
               a_comm_np_i = 0._xp
               if (rank == i) then
                  call compute_force_diag(m, r, 1, N/2, N, a)
                  a_comm_i = a(:,1:N/2)

                  ! compute interaction on right side
               else if (rank < i) then
                  call compute_force_mpi(m_np_i, r_np_i, 1, N, m(N/2+1:N), r(:,N/2+1:N), 1, N/2, a_comm_np_i, a(:,N/2+1:N))

                  ! compute interaction on left side
               else
                  call compute_force_mpi(m_i, r_i, 1, N/2, m, r, 1, N, a_comm_i, a)

               end if


               !--------------------------------
               ! Reduce accelerations
               !--------------------------------
               ! here, we get the accelerations from nodes on the right
               if (communicate_right) then
                  i_to_translate(1) = i
                  call mpi_group_translate_ranks(wgroup, 1, i_to_translate, mpi_group_to_right(i), i_translated, err)
                  ! Receive interaction in i
                  call mpi_reduce(a_comm_i, a, 3*N/2, MPI_REAL_XP, MPI_SUM, i_translated(1), mpi_comm_to_right(i), err)
               end if

               if (rank == i) a_comm_np_i = a

               ! here, we get the accelerations from nodes on the left
               if (communicate_left) then
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

         !---------------------------------------------------------------------------
         ! Compute the force using only own particles and communicating a subset to
         ! each other to spare memory
         !---------------------------------------------------------------------------
         else
            s = N / memory_factor

            allocate(a_comm_i(3, s))
            allocate(r_i(3, s))

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
            deallocate(a_comm_i)
         end if

      else
         !---------------------------------------------------------------------------
         ! Compute the force knowing all particles. Compute all the accelerations with
         ! a local subset of all the particles, then reduce all accelerations
         !---------------------------------------------------------------------------
         allocate(a_reduced(3, npoints))

         !--------------------------------
         ! Get the domain for integration from the number of nodes
         !--------------------------------
         if (flag_diag) then
            ! Compute a(t+dt) with fast version
            call compute_force_diag(m, r, N*rank/2 + 1, N*(rank + 1)/2, npoints, a)
         else
            ! Compute a(t+dt) with slow version
            call compute_force(m, r, N*rank + 1, N*(rank + 1), r, 1, npoints, a)
         end if

         !--------------------------------
         ! Reduce accelerations
         !--------------------------------
         call mpi_allreduce(a, a_reduced, npoints*3, MPI_REAL_XP, MPI_SUM, MPI_COMM_WORLD, err)
         a = a_reduced

         deallocate(a_reduced)

      end if

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

      real(xp) :: Ec_loc, Ep_loc

      real(xp) :: r_gathered(3, npoints), v_gathered(3, npoints), m_gathered(npoints)

      integer :: err

      Ec = 0._xp
      Ep = 0._xp

      if (flag_memory) then
         ! Gather r and v on master node !FIXME: Everyone should calculate its own energy

         call mpi_gather(r, 3*N, MPI_REAL_XP, r_gathered, 3*N, MPI_REAL_XP, 0, MPI_COMM_WORLD, err)
         call mpi_gather(v, 3*N, MPI_REAL_XP, v_gathered, 3*N, MPI_REAL_XP, 0, MPI_COMM_WORLD, err)
         call mpi_gather(m,   N, MPI_REAL_XP, m_gathered,   N, MPI_REAL_XP, 0, MPI_COMM_WORLD, err)

         if (rank == MASTER) then

            if (flag_diag) then
               ! Compute Ec & Ep at t+dt with fast version
               call compute_energy_diag(m_gathered, r_gathered, v_gathered, 1, npoints/2, npoints, Ec, Ep)
            else
               ! Compute Ec & Ep at t+dt with slow version
               ! TODO: This would be broken when computing using subdomains, because Ec must not be calculated elsewhere than the
               ! real diag
               call compute_energy(m_gathered, r_gathered, v_gathered, 1, npoints, Ec, Ep)
            end if

         end if
      else
         Ec_loc = 0._xp
         Ep_loc = 0._xp
         if (flag_diag) then
            ! Compute Ec & Ep at t+dt with fast version
            call compute_energy_diag(m, r, v, N*rank/2 + 1, N*(rank + 1)/2, npoints, Ec_loc, Ep_loc)
         else
            ! Compute Ec & Ep at t+dt with slow version
            call compute_energy(m, r, v, N*rank + 1, N*(rank + 1), Ec_loc, Ep_loc)
         end if
         !--------------------------------
         ! Reduce energies
         !--------------------------------
         call mpi_reduce(Ec_loc, Ec, 1, MPI_REAL_XP, MPI_SUM, MASTER, MPI_COMM_WORLD, err)
         call mpi_reduce(Ep_loc, Ep, 1, MPI_REAL_XP, MPI_SUM, MASTER, MPI_COMM_WORLD, err)

         r_gathered = r
         v_gathered = v
      end if

      if (rank == MASTER) then !FIXME: In flag_memory=.true. mode, everyone should write its own data
         call write_dump(iter, Ec, Ep, t, r_gathered, v_gathered)
      end if


   end subroutine compute_energy_wrap

end module physics

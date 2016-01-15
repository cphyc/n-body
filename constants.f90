module constants
   use mpi

   implicit none

   private

   integer,         parameter :: xp = selected_real_kind(8)
   integer,         parameter :: MPI_REAL_XP = MPI_DOUBLE_PRECISION

   real(kind = xp), parameter :: pi = 4.0_xp*atan(1.0_xp)
   real(kind = xp), parameter :: G  = 1._xp

   integer,         parameter :: npoints = 2**10
   real(kind = xp), parameter :: epsilon2 = npoints**(-2._xp/3._xp) / 400._xp ! We take 1/20th of the initial caracteristic distance

   real(kind = xp), parameter :: dt = 1.e-3_xp             ! Timestep
   real(kind = xp), parameter :: maxtime = npoints/128._xp ! Maximum time (ad hoc)
   integer,         parameter :: maxiter = maxtime/dt      ! Number of iteration

   integer,         parameter :: dump_freq = 10            ! Frequency at which the system is sampled
   integer,         parameter :: flag_compute_force = 3    ! Choose force computation subroutine, 0=sequential, 1=omp, 2=omp_nn_1
   integer,         parameter :: flag_compute_energy = 2   ! Choose energy computation subroutine, 0=sequential, 1=omp, 2=omp_nn_1
   integer,         parameter :: flag_compute_mpi = 0   ! Choose energy computation subroutine, 0=sequential, 1=omp, 2=omp_nn_1

   integer :: un, una ! Unit numbers for file operation

   public :: xp, MPI_REAL_XP, &
             pi, G, &
             npoints, epsilon2, dt, maxtime, maxiter, &
             dump_freq, flag_compute_force, flag_compute_energy, flag_compute_mpi, &
             un, una

contains

end module constants

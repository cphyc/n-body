module constants
   use mpi

   implicit none

   private

   integer,  parameter :: xp = selected_real_kind(8)
   integer,  parameter :: MPI_REAL_XP = MPI_DOUBLE_PRECISION

   real(xp), parameter :: pi = 4.0_xp*atan(1.0_xp)
   real(xp), parameter :: G  = 1._xp

   integer,  parameter :: npoints = 2**10
   real(xp), parameter :: epsilon2 = npoints**(-2._xp/3._xp) / 400._xp ! We take 1/20th of the initial caracteristic distance

   real(xp), parameter :: dt = 1.e-3_xp             ! Timestep
   real(xp), parameter :: maxtime = 8.e0_xp         ! Maximum time (ad hoc)
   integer,  parameter :: maxiter = int(maxtime/dt) ! Number of iteration

   integer,  parameter :: dump_freq = 10            ! Frequency at which the system is sampled

   logical,  parameter :: flag_diag = .true.        ! Choose force/energy computation subroutine, 0=dummy/1=diag
   integer,  parameter :: flag_mpi = 1              ! Choose mpi logic, 0=full copy of pos. arrays/1=low mem. impact
   real,     parameter :: maxmemory = 1.5           ! Only usefull when using flag_mpi=3, this is the max mem you can use
   integer,  parameter :: MASTER = 0                ! Choose which rank should be chosen as master node

   integer :: un, una ! Unit numbers for file operation

   public :: xp, MPI_REAL_XP, &
             pi, G, &
             npoints, epsilon2, dt, maxtime, maxiter, &
             dump_freq, flag_diag, flag_mpi, &
             maxmemory, MASTER, &
             un, una

contains

end module constants

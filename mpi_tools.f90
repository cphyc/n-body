module mpi_tools

  use constants
  use mpi

  implicit none

  integer, allocatable :: mpi_group_to_left(:), mpi_group_to_right(:), &
       mpi_comm_to_left(:), mpi_comm_to_right(:)
  integer :: wgroup

  private

  public :: initialize_mpi_groups, finalize_mpi_groups, &
       mpi_group_to_left, mpi_group_to_right, wgroup, &
       mpi_comm_to_left, mpi_comm_to_right

contains

  subroutine initialize_mpi_groups(nprocs, rank)
    integer, intent(out) :: nprocs
    integer, intent(out) :: rank

    integer :: ranges(3, 2)
    integer :: i
    integer :: err

    !---------------------------------------------
    ! Initialize MPI
    !---------------------------------------------
    call mpi_init(err)
    call mpi_comm_size(mpi_comm_world, nprocs, err)
    call mpi_comm_rank(mpi_comm_world, rank, err)
    !---------------------------------------------
    ! Create subgroups depending on flag
    !---------------------------------------------
    select case (flag_mpi)
    case(0)
    case(1)
    case(2)
       allocate(mpi_group_to_left(1:nprocs/2-1))
       allocate(mpi_group_to_right(1:nprocs/2-1))
       allocate(mpi_comm_to_left(1:nprocs/2-1))
       allocate(mpi_comm_to_right(1:nprocs/2-1))

       call mpi_comm_group(MPI_COMM_WORLD, wgroup, err)

       do i = 0, nprocs / 2 - 1
          !--------------------------------------
          ! Create two groups
          ! - mpi_groups_to_right all ranks ≥ i
          ! - mpi_groups_to_left all ranks ≤ i and nprocs-1-i
          !--------------------------------------
          if (rank == MASTER) then
             print*, 'Creating group right', i, 'from', i, 'to', nprocs - 1
             print*, '      and group left', i, 'from', 0, 'to', i, 'with', nprocs-i-1
          end if
          ranges(:, 1) = (/i, nprocs - 1, 1/)
          ranges(:, 2) = 0
          call mpi_group_range_incl(wgroup, 1, ranges, mpi_group_to_right(i), err)

          ranges(:, 1) = (/0, i, 1/)
          ranges(:, 2) = (/nprocs-1-i, nprocs-1-i, 1/)
          call mpi_group_range_incl(wgroup, 1, ranges, mpi_group_to_left(i), err)

          !--------------------------------------
          ! Create communicators along
          !--------------------------------------
          call mpi_comm_create(MPI_COMM_WORLD, mpi_group_to_right(i), mpi_comm_to_right(i), err)
          call mpi_comm_create(MPI_COMM_WORLD, mpi_group_to_left(i), mpi_comm_to_left(i), err)
       end do

    end select
  end subroutine initialize_mpi_groups

  subroutine finalize_mpi_groups()
    integer :: i, err
    integer :: i0, in

    i0 = lbound(mpi_group_to_left, 1)
    in = ubound(mpi_group_to_left, 1)

    do i = i0, in
       call mpi_comm_free(mpi_comm_to_left(i), err)
       call mpi_comm_free(mpi_comm_to_right(i), err)

       call mpi_group_free(mpi_group_to_left(i), err)
       call mpi_group_free(mpi_group_to_right(i), err)
    end do

    call mpi_finalize(err)
  end subroutine finalize_mpi_groups
end module mpi_tools

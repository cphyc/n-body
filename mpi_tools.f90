module mpi_tools

  use constants
  use mpi

  implicit none

  integer, allocatable :: mpi_group_to_left(:), mpi_group_to_right(:), &
       mpi_comm_to_left(:), mpi_comm_to_right(:)
  integer :: wgroup

  integer :: nprocs_internal, rank_internal

  private

  public :: initialize_mpi_groups, finalize_mpi_groups, &
       mpi_group_to_left, mpi_group_to_right, wgroup, &
       mpi_comm_to_left, mpi_comm_to_right, &
       communicate_right, communicate_left

contains

  subroutine initialize_mpi_groups(nprocs, rank)

    ! The subroutine creates multiple groups to enable efficient broadcasting. It creates two kind of groups
    ! 'to_right' groups, that link together the i-th process to the j≥i elements
    ! 'to_left' groups, that link together the i-th process to the j≤i elements, plus
    ! the nprocs-1-i element.
    !
    ! when a group is defined, one has also to create a communicator to enable
    ! communication among them


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

    nprocs_internal = nprocs
    rank_internal = rank
    !---------------------------------------------
    ! Create subgroups depending on flag
    !---------------------------------------------
    select case (flag_mpi)
    case(0, 2, 3)
    case(1)
       ! we initialize starting at 0, so that
       ! groups linking i are labelled by i
       allocate(mpi_group_to_left(0:nprocs/2-1))
       allocate(mpi_group_to_right(0:nprocs/2-1))
       allocate(mpi_comm_to_left(0:nprocs/2-1))
       allocate(mpi_comm_to_right(0:nprocs/2-1))

       call mpi_comm_group(MPI_COMM_WORLD, wgroup, err)

       do i = 0, nprocs / 2 - 1
          !--------------------------------------
          ! Create two groups
          ! - mpi_groups_to_right all ranks ≥ i
          ! - mpi_groups_to_left all ranks ≤ i and nprocs-1-i
          !--------------------------------------
          if (rank == MASTER) then
             write(*, '(a,i3,1x,a,i3,1x,a,i3)') &
                  'Creating group right n°', i, 'from', i, 'to', nprocs - 1
             write(*, '(a,i3,1x,a,i3,1x,a,i3,1x,a,i3)') &
                  '      and group left n°', i, 'from', 0, 'to', i, 'with', nprocs-i-1
          end if
          ranges(:, 1) = (/i, nprocs - 1, 1/)
          ranges(:, 2) = 0
          call mpi_group_range_incl(wgroup, 1, ranges, mpi_group_to_right(i), err)

          ranges(:, 1) = (/0, i, 1/)
          ranges(:, 2) = (/nprocs-1-i, nprocs-1-i, 1/)
          call mpi_group_range_incl(wgroup, 2, ranges, mpi_group_to_left(i), err)

          !--------------------------------------
          ! Create communicators along
          !--------------------------------------
          call mpi_comm_create(MPI_COMM_WORLD, mpi_group_to_right(i), mpi_comm_to_right(i), err)
          call mpi_comm_create(MPI_COMM_WORLD, mpi_group_to_left(i), mpi_comm_to_left(i), err)
       end do
    case(4)
       ! we initialize starting at 0, so that
       ! groups linking i are labelled by i
       allocate(mpi_group_to_left(0:nprocs-1))
       allocate(mpi_group_to_right(0:nprocs-1))
       allocate(mpi_comm_to_left(0:nprocs-1))
       allocate(mpi_comm_to_right(0:nprocs-1))

       call mpi_comm_group(MPI_COMM_WORLD, wgroup, err)

       do i = 0, nprocs - 1
          !--------------------------------------
          ! Create two groups
          ! - mpi_groups_to_right all ranks ≥ i
          ! - mpi_groups_to_left all ranks ≤ i and nprocs-1-i
          !--------------------------------------
          if (rank == MASTER) then
             write(*, '(a,i3,1x,a,i3,1x,a,i3)') &
                  'Creating group right n°', i, 'from', i, 'to', nprocs - 1
             write(*, '(a,i3,1x,a,i3,1x,a,i3,1x,a,i3)') &
                  '      and group left n°', i, 'from', 0, 'to', i
          end if
          ranges(:, 1) = (/i, nprocs - 1, 1/)
          ranges(:, 2) = 0
          call mpi_group_range_incl(wgroup, 1, ranges, mpi_group_to_right(i), err)

          ranges(:, 1) = (/0, i, 1/)
          ranges(:, 2) = 0
          call mpi_group_range_incl(wgroup, 1, ranges, mpi_group_to_left(i), err)

          !--------------------------------------
          ! Create communicators along
          !--------------------------------------
          call mpi_comm_create(MPI_COMM_WORLD, mpi_group_to_right(i), mpi_comm_to_right(i), err)
          call mpi_comm_create(MPI_COMM_WORLD, mpi_group_to_left(i), mpi_comm_to_left(i), err)
       end do
    end select
  end subroutine initialize_mpi_groups

  function communicate_right (rank, i)
    logical :: communicate_right
    integer, intent(in) :: rank, i

    communicate_right = rank >= i

  end function communicate_right

  function communicate_left (rank, i)
    logical :: communicate_left
    integer, intent(in) :: rank, i

    communicate_left = ((rank <= i) .or. (rank == nprocs_internal - i - 1))

  end function communicate_left


  subroutine finalize_mpi_groups()
    integer :: err

    call mpi_finalize(err)
  end subroutine finalize_mpi_groups
end module mpi_tools

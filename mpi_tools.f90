module mpi_tools

  use constants
  use mpi

  implicit none

  integer, allocatable :: mpi_group_to_left(:), mpi_group_to_right(:), &
                          mpi_comm_to_left(:), mpi_comm_to_right(:)
  integer :: wgroup

  private

  public :: initialize_mpi_groups, finalize_mpi_groups, wgroup, &
            mpi_group_to_left, mpi_group_to_right, &
            mpi_comm_to_left, mpi_comm_to_right

contains

  subroutine initialize_mpi_groups(nprocs, rank)

    ! The subroutine creates multiple groups to enable efficient broadcasting. It creates two kind of groups
    ! 'to_right' groups, that link together the i-th process to the j≥i elements
    ! 'to_left' groups, that link together the i-th process to the j≤i elements
    !
    ! when a group is defined, one has also to create a communicator to enable
    ! communication among them


    integer, intent(out) :: nprocs
    integer, intent(out) :: rank

    integer :: ranges(3,1)
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
    if (flag_memory) then
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
          ! - mpi_groups_to_left all ranks ≤ i
          !--------------------------------------
          if (rank == MASTER) then
             write(*, '(a,i3,1x,a,i3,1x,a,i3)')         'Creating group right of n°', i, 'from', i, 'to', nprocs - 1
             write(*, '(a,i3,1x,a,i3,1x,a,i3,1x,a,i3)') '      and group left of n°', i, 'from', 0, 'to', i
          end if

          ranges(:,1) = (/i, nprocs - 1, 1/)
          call mpi_group_range_incl(wgroup, 1, ranges, mpi_group_to_right(i), err)

          ranges(:,1) = (/0, i, 1/)
          call mpi_group_range_incl(wgroup, 1, ranges, mpi_group_to_left(i), err)

          !--------------------------------------
          ! Create communicators along
          !--------------------------------------
          call mpi_comm_create(MPI_COMM_WORLD, mpi_group_to_right(i), mpi_comm_to_right(i), err)
          call mpi_comm_create(MPI_COMM_WORLD, mpi_group_to_left(i), mpi_comm_to_left(i), err)
       end do

    end if

  end subroutine initialize_mpi_groups

  subroutine finalize_mpi_groups()
     implicit none

     integer :: err

     call mpi_finalize(err)

  end subroutine finalize_mpi_groups

end module mpi_tools

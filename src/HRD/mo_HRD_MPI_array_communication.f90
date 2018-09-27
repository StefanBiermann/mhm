!> \file mo_HRD_MPI_array_communication.f90

!> \authors Maren Kaluza
!> \date September 2018

MODULE mo_HRD_MPI_array_communication
  use mo_kind, only : i4, dp
  use mo_HRD_types, only: processSchedule, subtreeMeta
  !$ use omp_lib,      only: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
  use mpi

  ! Written Maren Kaluza, June 2018

  IMPLICIT NONE

  public :: distribute_array, collect_array, get_array, send_array
            
  private

CONTAINS
  ! two routines for the master process to distribute and collect an array
  ! cut into subarrays defined via the tree decomposition
  subroutine distribute_array(iBasin,nSubtrees,STmeta,permNodes,schedule,array)
    implicit none
    integer(i4),                     intent(in)  :: iBasin
    integer(i4),                     intent(in)  :: nSubtrees
    type(subtreeMeta), dimension(:), intent(in)  :: STmeta
    integer(i4),       dimension(:), intent(in)  :: permNodes
    type(processSchedule), dimension(:), intent(in) :: schedule
    integer(i4),       dimension(:), intent(in)  :: array
    ! local variables
    integer(i4) :: kk,jj,ii,iPerm,iproc,sizST,iST
    integer(i4) :: nproc,rank,ierror
    integer(i4), dimension(:), allocatable :: sendarray

    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierror)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

    do kk=1,nproc-1
       do jj=1,schedule(kk)%nTrees
          iST=schedule(kk)%trees(jj)
          sizST=STmeta(iST)%iEnd+1-STmeta(iST)%iStart
          allocate(sendarray(sizST))
          do ii=STmeta(iST)%iStart,STmeta(iST)%iEnd
             iPerm=permNodes(ii)
             sendarray(ii-STmeta(iST)%iStart+1)=array(iPerm)
          end do
          call MPI_Send(sendarray(1:sizST),sizST,MPI_INTEGER,kk,3,MPI_COMM_WORLD,ierror)
          deallocate(sendarray)
       end do
    end do
  end subroutine distribute_array

  ! collects data from the other processes into one array in
  ! original order
  subroutine collect_array(iBasin,nSubtrees,STmeta,permNodes,schedule,array)
    implicit none
    integer(i4),                     intent(in)    :: iBasin
    integer(i4),                     intent(in)    :: nSubtrees
    type(subtreeMeta), dimension(:), intent(in)    :: STmeta
    integer(i4),       dimension(:), intent(in)    :: permNodes
    type(processSchedule), dimension(:), intent(in) :: schedule
    integer(i4),       dimension(:), intent(inout) :: array
    ! local variables
    integer(i4) :: kk,jj,ii,iPerm,iproc,sizST,iST
    integer(i4) :: nproc,rank,ierror
    integer status(MPI_STATUS_SIZE)
    integer(i4), dimension(:), allocatable :: recvarray

    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierror)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

    do kk=1,nproc-1
       do jj=1,schedule(kk)%nTrees
          iST=schedule(kk)%trees(jj)
          sizST=STmeta(iST)%iEnd+1-STmeta(iST)%iStart
          allocate(recvarray(sizST))
          call MPI_Recv(recvarray(1:sizST),sizST,MPI_INTEGER,kk,4,MPI_COMM_WORLD,status,ierror)
          do ii=STmeta(iST)%iStart,STmeta(iST)%iEnd
             iPerm=permNodes(ii)
             array(iPerm)=recvarray(ii-STmeta(iST)%iStart+1)
          end do
          deallocate(recvarray)
       end do
    end do
  end subroutine collect_array

  ! two corresponding processes for the other processes to receive and send
  ! arrays assigned to that process back to the master process
  subroutine get_array(iBasin,STmeta,array)
    implicit none
    integer(i4),                                  intent(in)    :: iBasin
    type(subtreeMeta), dimension(:),              intent(in)    :: STmeta
    integer(i4),       dimension(:), allocatable, intent(inout) :: array
    ! local variables
    integer(i4) :: kk
    integer(i4) :: sizST,nST
    integer(i4) :: nproc,rank,ierror
    integer status(MPI_STATUS_SIZE)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierror)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

    nST=size(STmeta)
    allocate(array(STmeta(nST)%iEnd))
    ! ToDo: case: less subtrees than processes

    do kk=1,size(STmeta)
       sizST=STmeta(kk)%iEnd+1-STmeta(kk)%iStart
       call MPI_Recv(array(STmeta(kk)%iStart:STmeta(kk)%iEnd),sizST,MPI_INTEGER,0,3,MPI_COMM_WORLD,status,ierror)
    end do
  end subroutine get_array

  subroutine send_array(iBasin,STmeta,array)
    implicit none
    integer(i4),                     intent(in) :: iBasin
    type(subtreeMeta), dimension(:), intent(in) :: STmeta
    integer(i4),       dimension(:), intent(in) :: array
    ! local variables
    integer(i4) :: kk
    integer(i4) :: sizST
    integer(i4) :: nproc,rank,ierror
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierror)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

    ! ToDo: case: less subtrees than processes

    do kk=1,size(STmeta)
       sizST=STmeta(kk)%iEnd+1-STmeta(kk)%iStart
       call MPI_Send(array(STmeta(kk)%iStart:STmeta(kk)%iEnd),sizST,MPI_INTEGER,0,4,MPI_COMM_WORLD,ierror)
    end do

  end subroutine send_array

END MODULE mo_HRD_MPI_array_communication

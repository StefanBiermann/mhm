!> \file mo_HRD_MPI_array_communication.f90

!> \authors Maren Kaluza
!> \date September 2018

MODULE mo_HRD_MPI_array_communication
  use mo_kind, only : i4, dp
  use mo_HRD_types, only: processSchedule, subtreeMeta, MPI_Parameter
  !$ use omp_lib,      only: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
  use mpi_f08

  ! Written Maren Kaluza, June 2018

  IMPLICIT NONE

  public :: distribute_array, collect_array, get_array, send_array, &
            distribute_array_dp,             get_array_dp, &
            collect_array_dp,                send_array_dp, &
            distribute_array_log,            get_array_log, &
            distribute_full_array_dp,        get_full_array_dp, &
            collect_full_array_dp,           send_full_array_dp
            
  private

CONTAINS
  ! two routines for the master process to distribute and collect an array
  ! cut into subarrays defined via the tree decomposition
  subroutine distribute_array(iBasin,nproc,rank,comm,nSubtrees,STmeta,permNodes,schedule,array)
    implicit none
    integer(i4),                         intent(in)  :: iBasin
    integer(i4),                         intent(in)  :: nproc
    integer(i4),                         intent(in)  :: rank
    type(MPI_Comm)                                   :: comm
    integer(i4),                         intent(in)  :: nSubtrees
    type(subtreeMeta),     dimension(:), intent(in)  :: STmeta
    integer(i4),           dimension(:), intent(in)  :: permNodes
    type(processSchedule), dimension(:), intent(in)  :: schedule
    integer(i4),           dimension(:), intent(in)  :: array
    ! local variables
    integer(i4) :: kk, jj, ii, iPerm, iproc, sizST, iST, indST
    integer(i4) :: ierror
    integer(i4), dimension(:), allocatable :: sendarray

    do kk = 1, nproc-1
       do jj = 1, schedule(kk)%nTrees
          iST   = schedule(kk)%trees(jj)
          sizST = STmeta(iST)%sizST
          indST = STmeta(iST)%indST
          allocate(sendarray(sizST))
          do ii = STmeta(iST)%iStart, STmeta(iST)%iEnd
             iPerm = permNodes(ii)
             sendarray(ii - STmeta(iST)%iStart + 1) = array(iPerm)
          end do
          call MPI_Send(sendarray(1:sizST),sizST,MPI_INTEGER,kk,indST,comm,ierror)
          deallocate(sendarray)
       end do
    end do
  end subroutine distribute_array

  subroutine distribute_array_dp(iBasin,nproc,rank,comm,nSubtrees,STmeta,permNodes,schedule,array)
    implicit none
    integer(i4),                         intent(in)  :: iBasin
    integer(i4),                         intent(in)  :: nproc
    integer(i4),                         intent(in)  :: rank
    type(MPI_Comm)                                   :: comm
    integer(i4),                         intent(in)  :: nSubtrees
    type(subtreeMeta),     dimension(:), intent(in)  :: STmeta
    integer(i4),           dimension(:), intent(in)  :: permNodes
    type(processSchedule), dimension(:), intent(in)  :: schedule
    real(dp),              dimension(:), intent(in)  :: array
    ! local variables
    integer(i4) :: kk, jj, ii, iPerm, iproc, sizST, iST, indST
    integer(i4) :: ierror
    real(dp), dimension(:), allocatable :: sendarray

    do kk = 1, nproc-1
       do jj = 1, schedule(kk)%nTrees
          iST   = schedule(kk)%trees(jj)
          sizST = STmeta(iST)%sizST
          indST = STmeta(iST)%indST
          allocate(sendarray(sizST))
          do ii = STmeta(iST)%iStart, STmeta(iST)%iEnd
             iPerm = permNodes(ii)
             sendarray(ii - STmeta(iST)%iStart + 1) = array(iPerm)
          end do
          call MPI_Send(sendarray(1:sizST),sizST,MPI_DOUBLE_PRECISION,kk,indST,comm,ierror)
          deallocate(sendarray)
       end do
    end do
  end subroutine distribute_array_dp

  subroutine distribute_full_array_dp(iBasin,MPIparam,nSubtrees,STmeta,permNodes,schedule,array)
    implicit none
    integer(i4),                            intent(in)  :: iBasin
    type(MPI_parameter),                    intent(in)  :: MPIparam
    integer(i4),                            intent(in)  :: nSubtrees
    type(subtreeMeta),     dimension(:),    intent(in)  :: STmeta
    integer(i4),           dimension(:),    intent(in)  :: permNodes
    type(processSchedule), dimension(:),    intent(in)  :: schedule
    real(dp),              dimension(:, :), intent(in)  :: array
    ! local variables
    integer(i4) :: kk, jj, ii, hh, iPerm, iproc, sizST, iST, indST
    integer(i4) :: ierror
    real(dp), dimension(:, :), allocatable :: sendarray

    do kk = 1, MPIParam%nproc-1
      do jj = 1, schedule(kk)%nTrees
        iST   = schedule(kk)%trees(jj)
        sizST = STmeta(iST)%sizST
        indST = STmeta(iST)%indST
        allocate(sendarray(sizST,MPIparam%bufferLength))
        do ii = STmeta(iST)%iStart, STmeta(iST)%iEnd
          do hh = 1, MPIparam%bufferLength
            iPerm = permNodes(ii)
            sendarray(ii - STmeta(iST)%iStart + 1, hh) = array(iPerm, hh)
          end do
        end do
        call MPI_Send(sendarray(:, :),sizST*MPIparam%bufferLength, &
                      MPI_DOUBLE_PRECISION,kk,indST,MPIparam%comm,ierror)
        deallocate(sendarray)
      end do
    end do
  end subroutine distribute_full_array_dp

  subroutine distribute_array_log(iBasin,nproc,rank,comm,nSubtrees,STmeta,permNodes,schedule,array)
    implicit none
    integer(i4),                         intent(in)  :: iBasin
    integer(i4),                         intent(in)  :: nproc
    integer(i4),                         intent(in)  :: rank
    type(MPI_Comm)                                   :: comm
    integer(i4),                         intent(in)  :: nSubtrees
    type(subtreeMeta),     dimension(:), intent(in)  :: STmeta
    integer(i4),           dimension(:), intent(in)  :: permNodes
    type(processSchedule), dimension(:), intent(in)  :: schedule
    logical,               dimension(:), intent(in)  :: array
    ! local variables
    integer(i4) :: kk, jj, ii, iPerm, iproc, sizST, iST, indST
    integer(i4) :: ierror
    logical, dimension(:), allocatable :: sendarray

    do kk = 1, nproc-1
       do jj = 1, schedule(kk)%nTrees
          iST   = schedule(kk)%trees(jj)
          sizST = STmeta(iST)%sizST
          indST = STmeta(iST)%indST
          allocate(sendarray(sizST))
          do ii = STmeta(iST)%iStart, STmeta(iST)%iEnd
             iPerm = permNodes(ii)
             sendarray(ii - STmeta(iST)%iStart + 1) = array(iPerm)
          end do
          call MPI_Send(sendarray(1:sizST),sizST,MPI_LOGICAL,kk,indST,comm,ierror)
          deallocate(sendarray)
       end do
    end do
  end subroutine distribute_array_log

  ! collects data from the other processes into one array in
  ! original order
  subroutine collect_array(iBasin,nproc,rank,comm,nSubtrees,STmeta,permNodes,schedule,array)
    implicit none
    integer(i4),                         intent(in)    :: iBasin
    integer(i4),                         intent(in)    :: nproc
    integer(i4),                         intent(in)    :: rank
    type(MPI_Comm)                                     :: comm
    integer(i4),                         intent(in)    :: nSubtrees
    type(subtreeMeta),     dimension(:), intent(in)    :: STmeta
    integer(i4),           dimension(:), intent(in)    :: permNodes
    type(processSchedule), dimension(:), intent(in)    :: schedule
    integer(i4),           dimension(:), intent(inout) :: array
    ! local variables
    integer(i4)      :: kk, jj, ii, iPerm, iproc, sizST, iST, indST
    integer(i4)      :: ierror
    type(MPI_Status) :: status
    integer(i4), dimension(:), allocatable :: recvarray

    do kk = 1, nproc-1
       do jj = 1, schedule(kk)%nTrees
          iST   = schedule(kk)%trees(jj)
          sizST = STmeta(iST)%sizST
          indST = STmeta(iST)%indST
          allocate(recvarray(sizST))
          call MPI_Recv(recvarray(1:sizST),sizST,MPI_INTEGER,kk,indST,comm,status,ierror)
          do ii = STmeta(iST)%iStart, STmeta(iST)%iEnd
             iPerm = permNodes(ii)
             array(iPerm) = recvarray(ii - STmeta(iST)%iStart + 1)
          end do
          deallocate(recvarray)
       end do
    end do
  end subroutine collect_array

  subroutine collect_array_dp(iBasin,nproc,rank,comm,nSubtrees,STmeta,permNodes,schedule,array)
    implicit none
    integer(i4),                         intent(in)    :: iBasin
    integer(i4),                         intent(in)    :: nproc
    integer(i4),                         intent(in)    :: rank
    type(MPI_Comm)                                     :: comm
    integer(i4),                         intent(in)    :: nSubtrees
    type(subtreeMeta),     dimension(:), intent(in)    :: STmeta
    integer(i4),           dimension(:), intent(in)    :: permNodes
    type(processSchedule), dimension(:), intent(in)    :: schedule
    real(dp),              dimension(:), intent(inout) :: array
    ! local variables
    integer(i4)      :: kk, jj, ii, iPerm, iproc, sizST, iST, indST
    integer(i4)      :: ierror
    type(MPI_Status) :: status
    real(dp), dimension(:), allocatable :: recvarray

    do kk = 1, nproc-1
       do jj = 1, schedule(kk)%nTrees
          iST   = schedule(kk)%trees(jj)
          sizST = STmeta(iST)%sizST
          indST = STmeta(iST)%indST
          allocate(recvarray(sizST))
          call MPI_Recv(recvarray(1:sizST),sizST,MPI_DOUBLE_PRECISION,kk,indST,comm,status,ierror)
          do ii = STmeta(iST)%iStart, STmeta(iST)%iEnd
             iPerm = permNodes(ii)
             array(iPerm) = recvarray(ii - STmeta(iST)%iStart + 1)
          end do
          deallocate(recvarray)
       end do
    end do
  end subroutine collect_array_dp

  subroutine collect_full_array_dp(iBasin,MPIparam,nSubtrees,STmeta,permNodes,schedule,array)
    implicit none
    integer(i4),                            intent(in)    :: iBasin
    type(MPI_parameter),                    intent(in)    :: MPIparam
    integer(i4),                            intent(in)    :: nSubtrees
    type(subtreeMeta),     dimension(:),    intent(in)    :: STmeta
    integer(i4),           dimension(:),    intent(in)    :: permNodes
    type(processSchedule), dimension(:),    intent(in)    :: schedule
    real(dp),              dimension(:, :), intent(inout) :: array
    ! local variables
    integer(i4)      :: kk, jj, ii, hh, iPerm, iproc, sizST, iST, indST
    integer(i4)      :: ierror
    type(MPI_Status) :: status
    real(dp), dimension(:, :), allocatable :: recvarray

    do kk = 1, MPIparam%nproc-1
      do jj = 1, schedule(kk)%nTrees
        iST   = schedule(kk)%trees(jj)
        sizST = STmeta(iST)%sizST
        indST = STmeta(iST)%indST
        allocate(recvarray(sizST, MPIparam%bufferLength))
        call MPI_Recv(recvarray(1:sizST,:),sizST*MPIparam%bufferLength,&
                        MPI_DOUBLE_PRECISION,kk,indST,MPIparam%comm,status,ierror)
        do ii = STmeta(iST)%iStart, STmeta(iST)%iEnd
          do hh = 1, MPIparam%bufferLength
            iPerm = permNodes(ii)
            array(iPerm, hh) = recvarray(ii - STmeta(iST)%iStart + 1, hh)
          end do
        end do
        deallocate(recvarray)
      end do
    end do
  end subroutine collect_full_array_dp

  ! two corresponding processes for the other processes to receive and send
  ! arrays assigned to that process back to the master process
  subroutine get_array(iBasin,nproc,rank,comm,STmeta,array)
    implicit none
    integer(i4),                                  intent(in)    :: iBasin
    integer(i4),                                  intent(in)    :: nproc
    integer(i4),                                  intent(in)    :: rank
    type(MPI_Comm)                                              :: comm
    type(subtreeMeta), dimension(:),              intent(in)    :: STmeta
    integer(i4),       dimension(:), allocatable, intent(inout) :: array
    ! local variables
    integer(i4)      :: kk
    integer(i4)      :: sizST,nST
    integer(i4)      :: ierror
    type(MPI_Status) :: status

    nST = size(STmeta)
    if (nST > 0) then
      allocate(array(STmeta(nST)%iEnd))
    else
      allocate(array(0))
    end if
    ! ToDo: case: less subtrees than processes

    do kk = 1, size(STmeta)
       sizST = STmeta(kk)%sizST
       call MPI_Recv(array(STmeta(kk)%iStart:STmeta(kk)%iEnd),sizST,MPI_INTEGER, &
                                            0,STmeta(kk)%indST,comm,status,ierror)
    end do
  end subroutine get_array

  subroutine get_array_dp(iBasin,nproc,rank,comm,STmeta,array)
    implicit none
    integer(i4),                                  intent(in)    :: iBasin
    integer(i4),                                  intent(in)    :: nproc
    integer(i4),                                  intent(in)    :: rank
    type(MPI_Comm)                                              :: comm
    type(subtreeMeta), dimension(:),              intent(in)    :: STmeta
    real(dp),       dimension(:), allocatable, intent(inout)    :: array
    ! local variables
    integer(i4)      :: kk
    integer(i4)      :: sizST,nST
    integer(i4)      :: ierror
    type(MPI_Status) :: status

    nST = size(STmeta)
    if (nST > 0) then
      allocate(array(STmeta(nST)%iEnd))
    else
      allocate(array(0))
    end if
    ! ToDo: case: less subtrees than processes

    do kk = 1, size(STmeta)
       sizST = STmeta(kk)%sizST
       call MPI_Recv(array(STmeta(kk)%iStart:STmeta(kk)%iEnd),sizST,MPI_DOUBLE_PRECISION, &
                                            0,STmeta(kk)%indST,comm,status,ierror)
    end do
  end subroutine get_array_dp

  subroutine get_full_array_dp(iBasin,MPIparam,STmeta,array)
    implicit none
    integer(i4),                                  intent(in)    :: iBasin
    type(MPI_parameter),                          intent(in)    :: MPIparam
    type(subtreeMeta), dimension(:),              intent(in)    :: STmeta
    real(dp),       dimension(:, :), allocatable, intent(inout) :: array
    ! local variables
    real(dp), dimension(:, :), allocatable :: recvArray
    integer(i4)      :: kk
    integer(i4)      :: sizST,nST
    integer(i4)      :: ierror
    type(MPI_Status) :: status

    nST = size(STmeta)
    if (nST > 0) then
      allocate(array(STmeta(nST)%iEnd, MPIparam%bufferLength))
    else
      allocate(array(0, 0))
    end if
    ! ToDo: case: less subtrees than processes

    do kk = 1, size(STmeta)
       sizST = STmeta(kk)%sizST
       allocate(recvArray(sizST, MPIparam%bufferLength))
       call MPI_Recv(recvArray,&
                                   sizST*MPIparam%bufferLength,MPI_DOUBLE_PRECISION, &
                                            0,STmeta(kk)%indST,MPIparam%comm,status,ierror)
       array(STmeta(kk)%iStart:STmeta(kk)%iEnd, :) = recvArray(:, :)
       deallocate(recvArray)
    end do
  end subroutine get_full_array_dp

  subroutine get_array_log(iBasin,nproc,rank,comm,STmeta,array)
    implicit none
    integer(i4),                                  intent(in)    :: iBasin
    integer(i4),                                  intent(in)    :: nproc
    integer(i4),                                  intent(in)    :: rank
    type(MPI_Comm)                                              :: comm
    type(subtreeMeta), dimension(:),              intent(in)    :: STmeta
    logical,        dimension(:), allocatable, intent(inout)    :: array
    ! local variables
    integer(i4)      :: kk
    integer(i4)      :: sizST,nST
    integer(i4)      :: ierror
    type(MPI_Status) :: status

    nST = size(STmeta)
    if (nST > 0) then
      allocate(array(STmeta(nST)%iEnd))
    else
      allocate(array(0))
    end if
    ! ToDo: case: less subtrees than processes

    do kk = 1, size(STmeta)
       sizST = STmeta(kk)%sizST
       call MPI_Recv(array(STmeta(kk)%iStart:STmeta(kk)%iEnd),sizST,MPI_LOGICAL, &
                                            0,STmeta(kk)%indST,comm,status,ierror)
    end do
  end subroutine get_array_log

  subroutine send_array(iBasin,nproc,rank,comm,STmeta,array)
    implicit none
    integer(i4),                     intent(in) :: iBasin
    integer(i4),                     intent(in) :: nproc
    integer(i4),                     intent(in) :: rank
    type(MPI_Comm)                              :: comm
    type(subtreeMeta), dimension(:), intent(in) :: STmeta
    integer(i4),       dimension(:), intent(in) :: array
    ! local variables
    integer(i4) :: kk
    integer(i4) :: sizST
    integer(i4) :: ierror

    ! ToDo: case: less subtrees than processes

    do kk = 1, size(STmeta)
       sizST = STmeta(kk)%sizST
       call MPI_Send(array(STmeta(kk)%iStart:STmeta(kk)%iEnd),sizST,MPI_INTEGER, &
                                                   0,STmeta(kk)%indST,comm,ierror)
    end do

  end subroutine send_array

  subroutine send_array_dp(iBasin,nproc,rank,comm,STmeta,array)
    implicit none
    integer(i4),                     intent(in) :: iBasin
    integer(i4),                     intent(in) :: nproc
    integer(i4),                     intent(in) :: rank
    type(MPI_Comm)                              :: comm
    type(subtreeMeta), dimension(:), intent(in) :: STmeta
    real(dp),          dimension(:), intent(in) :: array
    ! local variables
    integer(i4) :: kk
    integer(i4) :: sizST
    integer(i4) :: ierror

    ! ToDo: case: less subtrees than processes

    do kk = 1, size(STmeta)
       sizST = STmeta(kk)%sizST
       call MPI_Send(array(STmeta(kk)%iStart:STmeta(kk)%iEnd),sizST,MPI_DOUBLE_PRECISION, &
                                                   0,STmeta(kk)%indST,comm,ierror)
    end do

  end subroutine send_array_dp

  subroutine send_full_array_dp(iBasin,MPIparam,STmeta,array)
    implicit none
    integer(i4),                        intent(in) :: iBasin
    type(MPI_parameter),                intent(in) :: MPIparam
    type(subtreeMeta), dimension(:),    intent(in) :: STmeta
    real(dp),          dimension(:, :), intent(in) :: array
    ! local variables
    real(dp), dimension(:, :), allocatable :: sendArray
    integer(i4) :: kk
    integer(i4) :: sizST
    integer(i4) :: ierror

    ! ToDo: case: less subtrees than processes

    do kk = 1, size(STmeta)
       sizST = STmeta(kk)%sizST
       allocate(sendArray(sizST, MPIparam%bufferLength))
       sendArray(:, :) = array(STmeta(kk)%iStart:STmeta(kk)%iEnd, :)
       call MPI_Send(sendArray ,sizST*MPIparam%bufferLength,MPI_DOUBLE_PRECISION, &
                                                   0,STmeta(kk)%indST,MPIparam%comm,ierror)
       deallocate(sendArray)
    end do

  end subroutine send_full_array_dp

END MODULE mo_HRD_MPI_array_communication

!> \file mo_HRD_subtree_meta.f90

!> \authors Maren Kaluza
!> \date June 2018

MODULE mo_HRD_subtree_meta
  use mo_kind, only : i4, dp
  use mo_HRD_types, only : ptrTreeNode, subtreeMeta, processSchedule, &
          subtreeBuffer
  use mo_HRD_tree_tools, only : write_tree_to_array
  !$ use omp_lib,      only: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
  use mpi_f08

  ! Written Maren Kaluza, June 2018

  IMPLICIT NONE

  public :: init_subtree_metadata, distribute_subtree_meta, get_subtree_meta, &
            init_subtree_buffer, destroy_subtree_buffer
            

  private

CONTAINS

  ! A subtree data structure makes communication between the subtrees
  ! much easier for the master. Processing the data is more efficient
  ! with array, so everything gets written into a nice array in
  ! routing order. Therefore we need a permutation array permNodes
  ! which is at the same time fromNodes, and a toNode array
  subroutine init_subtree_metadata(iBasin,nNodes,subtrees,STmeta,permNodes,toNodes,toInNodes)
    implicit none
    integer(i4),                                  intent(in)     :: iBasin
    integer(i4),                                  intent(in)     :: nNodes
    type(ptrTreeNode), dimension(:),              intent(in)     :: subtrees ! the array of
    type(subtreeMeta), dimension(:),              intent(inout)  :: STmeta
    integer(i4),       dimension(:), allocatable, intent(inout)  :: permNodes
    integer(i4),       dimension(:), allocatable, intent(inout)  :: toNodes
    integer(i4),       dimension(:), allocatable, intent(inout)  :: toInNodes
    ! local
    integer(i4) :: nSubtrees
    integer(i4) :: kk, jj, ind
    integer(i4) :: nInNodes

    nSubtrees=size(STmeta)
    nInNodes = 0
    do kk = 1, nSubtrees
      nInNodes = nInNodes + subtrees(kk)%tN%ST%NpraeST
    end do
    allocate(permNodes(nNodes), toNodes(nNodes), toInNodes(nInNodes))

    STmeta(1)%indST=1
    STmeta(1)%sizST=subtrees(1)%tN%ST%sizST
    STmeta(1)%iStart=1
    STmeta(1)%iEnd=STmeta(1)%iStart+STmeta(1)%sizST-1
    STmeta(1)%nIn=subtrees(1)%tN%ST%NpraeST
    STmeta(1)%iInStart=1
    STmeta(1)%iInEnd=STmeta(1)%iInStart+STmeta(1)%nIn-1
    do kk=2,nSubtrees
       STmeta(kk)%indST=kk
       STmeta(kk)%sizST=subtrees(kk)%tN%ST%sizST
       STmeta(kk)%iStart=Stmeta(kk-1)%iEnd+1
       STmeta(kk)%iEnd=STmeta(kk)%iStart+STmeta(kk)%sizST-1
       STmeta(kk)%nIn=subtrees(kk)%tN%ST%NpraeST
       STmeta(kk)%iInStart=STmeta(kk-1)%iInEnd+1
       STmeta(kk)%iInEnd=STmeta(kk)%iInStart+STmeta(kk)%nIn-1
    end do
    toNodes(:)=0
    toInNodes(:)=0
    do kk=nSubtrees,1,-1
       ind = STmeta(kk)%sizST
       call write_tree_to_array(subtrees(kk),STmeta(kk)%iStart,ind,&
            permNodes(STmeta(kk)%iStart:STmeta(kk)%iEnd), &
              toNodes(STmeta(kk)%iStart:STmeta(kk)%iEnd))
       do jj = 1, STmeta(kk)%nIn
         ind = subtrees(kk)%tN%ST%praeST(jj)%tN%post%tN%ind
         toInNodes(STmeta(kk)%iInStart+jj-1) = ind
       end do
    end do
    ! ToDo: repair this part
    ! correct last toNode
    do kk=1,nSubtrees
       if (associated(subtrees(kk)%tN%post%tN)) then
          toNodes(STmeta(kk)%iEnd)=subtrees(kk)%tN%post%tN%ind
       end if
    end do
  end subroutine init_subtree_metadata

  subroutine distribute_subtree_meta(iBasin,nproc,comm,nSubtrees,STmeta,toNodes,toInNodes,schedule,subtrees)
    implicit none
    integer(i4),                         intent(in)    :: iBasin
    integer(i4),                         intent(in)    :: nproc
    type(MPI_Comm)                                     :: comm
    integer(i4),                         intent(in)    :: nSubtrees
    type(subtreeMeta),     dimension(:), intent(in)    :: STmeta
    integer(i4),           dimension(:), intent(in)    :: toNodes
    integer(i4),           dimension(:), intent(in)    :: toInNodes
    type(processSchedule), dimension(:), intent(inout) :: schedule
    type(ptrTreeNode),     dimension(:), intent(inout) :: subtrees ! the array of
    ! local variables
    integer(i4) :: kk,ii,jj,i,iPerm,iproc,sizST,iST,nIn
    integer(i4) :: rank,ierror
    integer(i4), dimension(:,:), allocatable :: iSends ! the number and over all
                                                       ! length of arrays to be send to a process
    integer(i4), dimension(:), allocatable :: sendarray

    ! send number of subtrees and total number of tree nodes assigned
    ! to the processes to the corresponding process, so arrays can be
    ! allocated in the receiving subroutines
    allocate(iSends(3,nproc-1))
    do kk=1,nproc-1
       iSends(1,kk)=schedule(kk)%nTrees
       iSends(2,kk)=schedule(kk)%overallSize
       iSends(3,kk)=schedule(kk)%overallInSize
       call MPI_Send(iSends(:,kk),3,MPI_INTEGER,kk,0,comm,ierror)
    end do
    ! send metadata of the subtrees to the nodes where they are assigned to
    ! size, identifying index corresponding to the subtree array and number of
    ! nodes, where data is flowing into the tree
    do kk=1,nproc-1
       do jj=1,schedule(kk)%nTrees
          iST=schedule(kk)%trees(jj)
          sizST=STmeta(iST)%sizST
          call MPI_Send(iST,1,MPI_INTEGER,kk,1,comm,ierror)
          call MPI_Send(sizST,1,MPI_INTEGER,kk,iST,comm,ierror)
          call MPI_Send(STmeta(iST)%nIn,1,MPI_INTEGER,kk,iST,comm,ierror)
       end do
    end do
    ! for each subtree send corresponding toNodes to the node. Move indices from
    ! toNodes, so they start from 0
    ! They can then be moved on the receiving process corresponding to
    ! its own offset
    do kk=1,nproc-1
       do jj=1,schedule(kk)%nTrees
          iST=schedule(kk)%trees(jj)
          sizST=STmeta(iST)%sizST
          nIn=STmeta(iST)%nIn
          allocate(sendarray(sizST))
          do ii=STmeta(iST)%iStart,STmeta(iST)%iEnd
             !ToDo: is the minus term still correct?
             sendarray(ii-STmeta(iST)%iStart+1)=toNodes(ii)-STmeta(iST)%iStart+1
          end do
          call MPI_Send(sendarray(1:sizST),sizST,MPI_INTEGER,kk,iST,comm,ierror)
          ! ToDo: should do with the old sendarray
          deallocate(sendarray)
          if (nIn > 0) then
             allocate(sendarray(nIn))
             do ii=STmeta(iST)%iInStart,STmeta(iST)%iInEnd
                !ToDo: is the minus term still correct?
                sendarray(ii-STmeta(iST)%iInStart+1)=toInNodes(ii)-STmeta(iST)%iStart+1
             end do
             call MPI_Send(sendarray(1:nIn),nIn,MPI_INTEGER,kk,iST,comm,ierror)
             ! send indices of subtrees of halo leaves
             do ii=STmeta(iST)%iInStart,STmeta(iST)%iInEnd
                i = ii-STmeta(iST)%iInStart+1
                sendarray(i)=subtrees(iST)%tN%ST%praeST(i)%tN%ST%indST
             end do
             call MPI_Send(sendarray(1:nIn),nIn,MPI_INTEGER,kk,iST,comm,ierror)
             deallocate(sendarray)
          end if
       end do
    end do

    deallocate(iSends)
  end subroutine distribute_subtree_meta

  subroutine get_subtree_meta(iBasin,comm,STmeta,toNodes,toInNodes,inInds)
    implicit none
    integer(i4),                                  intent(in)    :: iBasin
    type(MPI_Comm)                                              :: comm
    type(subtreeMeta), dimension(:), allocatable, intent(inout) :: STmeta
    integer(i4),       dimension(:), allocatable, intent(inout) :: toNodes
    integer(i4),       dimension(:), allocatable, intent(inout) :: toInNodes
    integer(i4),       dimension(:), allocatable, intent(inout) :: inInds
    ! local variables
    integer(i4) :: kk
    integer(i4), dimension(3) :: nDatasets ! number of incoming data sets
                                           ! total size of datasets
    integer(i4)      :: sizST, indST, nIn
    integer(i4)      :: nSubtrees, totSizeOfSubtrees, totNIn
    integer(i4)      :: ierror
    type(MPI_Status) :: status

    ! recieves number of subtrees and total number of tree nodes assigned to
    ! this process
    call MPI_Recv(nDatasets(:),3,MPI_INTEGER,0,0,comm,status,ierror)
    nSubtrees=nDatasets(1)
    totSizeOfSubtrees=nDatasets(2)
    totNIn=nDatasets(3)
    allocate(STmeta(nSubtrees))
    allocate(toNodes(totSizeOfSubtrees))
    allocate(toInNodes(totNIn),inInds(totNIn))
    if (nSubtrees > 0) then
      call MPI_Recv(STmeta(1)%indST,1,MPI_INTEGER,0,1,comm,status,ierror)
      call MPI_Recv(STmeta(1)%sizST,1,MPI_INTEGER,0,STmeta(1)%indST,comm,status,ierror)
      call MPI_Recv(STmeta(1)%nIn,1,MPI_INTEGER,0,STmeta(1)%indST,comm,status,ierror)
      STmeta(1)%iStart=1
      STmeta(1)%iEnd=STmeta(1)%iStart+STmeta(1)%sizST-1
      STmeta(1)%iInStart=1
      STmeta(1)%iInEnd=STmeta(1)%iInStart+STmeta(1)%nIn-1
    end if
    do kk=2,nSubtrees
      call MPI_Recv(STmeta(kk)%indST,1,MPI_INTEGER,0,1,comm,status,ierror)
      call MPI_Recv(STmeta(kk)%sizST,1,MPI_INTEGER,0,STmeta(kk)%indST,comm,status,ierror)
      call MPI_Recv(STmeta(kk)%nIn,1,MPI_INTEGER,0,STmeta(kk)%indST,comm,status,ierror)
      STmeta(kk)%iStart=Stmeta(kk-1)%iEnd+1
      STmeta(kk)%iEnd=STmeta(kk)%iStart+STmeta(kk)%sizST-1
      STmeta(kk)%iInStart=STmeta(kk-1)%iInEnd+1
      STmeta(kk)%iInEnd=STmeta(kk)%iInStart+STmeta(kk)%nIn-1
    end do

    do kk=1,nSubtrees
      sizST=STmeta(kk)%sizST
      nIn=STmeta(kk)%nIn
      call MPI_Recv(toNodes(STmeta(kk)%iStart:STmeta(kk)%iEnd),sizST,MPI_INTEGER,0,STmeta(kk)%indST,comm,status,ierror)
      if (nIn > 0) then
         call MPI_Recv(toInNodes(STmeta(kk)%iInStart:STmeta(kk)%iInEnd),nIn,MPI_INTEGER,0,STmeta(kk)%indST,comm,status,ierror)
         call MPI_Recv(inInds(STmeta(kk)%iInStart:STmeta(kk)%iInEnd),nIn,MPI_INTEGER,0,STmeta(kk)%indST,comm,status,ierror)
      end if
      ! ToDo: why not -1
      ! the toNodes have been moved, so they start from index 1, before sending
      ! now the last entry gets moved to the starting point of the subtree in the array
      ! toNodes(STmeta(kk)%iStart:STmeta(kk)%iEnd)=toNodes(STmeta(kk)%iStart:STmeta(kk)%iEnd)+STmeta(kk)%iStart
      toNodes(STmeta(kk)%iEnd)=toNodes(STmeta(kk)%iEnd)+STmeta(kk)%iStart
    end do

  end subroutine get_subtree_meta

  subroutine init_subtree_buffer(STmeta, bufferLength, STbuffer)
    type(subtreeMeta),   dimension(:),              intent(in)  :: STmeta
    integer(i4)                                                 :: bufferLength
    type(subtreeBuffer), dimension(:), allocatable, intent(out) :: STbuffer
    ! local variables
    integer(i4) :: nSubtrees
    integer(i4) :: i

    nSubtrees = size(STmeta)
    allocate(STbuffer(nSubtrees))

    ! initiate buffering space for all data send from the subtree roots to other
    ! subtree leaves. nIn is the number of leaves, so we need nIn+1 space for
    ! every message collective
    do i = 1, nSubtrees
      allocate(STbuffer(i)%statuses(STmeta(i)%nIn+1))
      allocate(STbuffer(i)%requests(STmeta(i)%nIn+1))
      allocate(STbuffer(i)%buffer(bufferLength+1, STmeta(i)%nIn+1))
    end do

  end subroutine init_subtree_buffer

  subroutine destroy_subtree_buffer(STbuffer)
    type(subtreeBuffer), dimension(:), allocatable, intent(inout) :: STbuffer
    ! local variables
    integer(i4) :: nSubtrees
    integer(i4) :: i

    nSubtrees = size(STbuffer)

    do i = 1, nSubtrees
      deallocate(STbuffer(i)%statuses)
      deallocate(STbuffer(i)%requests)
      deallocate(STbuffer(i)%buffer)
    end do

    deallocate(STbuffer)

  end subroutine destroy_subtree_buffer

END MODULE mo_HRD_subtree_meta

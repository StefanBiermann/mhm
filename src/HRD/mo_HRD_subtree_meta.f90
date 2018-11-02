!> \file mo_HRD_subtree_meta.f90

!> \authors Maren Kaluza
!> \date June 2018

MODULE mo_HRD_subtree_meta
  use mo_kind, only : i4, dp
  use mo_HRD_types, only : ptrTreeNode, subtreeMeta, processSchedule 
  use mo_HRD_tree_tools, only : write_tree_to_array
  !$ use omp_lib,      only: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
  use mpi

  ! Written Maren Kaluza, June 2018

  IMPLICIT NONE

  public :: init_subtree_metadata, distribute_subtree_meta, get_subtree_meta
            

  private

CONTAINS

  ! A subtree data structure makes communication between the subtrees
  ! much easier for the master. Processing the data is more efficient
  ! with array, so everything gets written into a nice array in
  ! routing order. Therefore we need a permutation array permNodes
  ! which is at the same time fromNodes, and a toNode array
  subroutine init_subtree_metadata(iBasin,subtrees,STmeta,permNodes,toNodes)
    implicit none
    integer(i4),                     intent(in)     :: iBasin
    type(ptrTreeNode), dimension(:), intent(in)     :: subtrees ! the array of
    type(subtreeMeta), dimension(:), intent(inout)  :: STmeta
    integer(i4),       dimension(:), intent(inout)  :: permNodes
    integer(i4),       dimension(:), intent(inout)  :: toNodes
    ! local
    integer(i4) :: nNodes, nSubtrees
    integer(i4) :: kk,ind

    nSubtrees=size(STmeta)
    nNodes=size(permNodes)

    STmeta(1)%iStart=1
    STmeta(1)%iEnd=subtrees(1)%tN%ST%sizST
    STmeta(1)%indST=1
    STmeta(1)%nIn=subtrees(1)%tN%ST%NpraeST
    do kk=2,nSubtrees
       STmeta(kk)%iStart=Stmeta(kk-1)%iEnd+1
       STmeta(kk)%iEnd=STmeta(kk)%iStart+subtrees(kk)%tN%ST%sizST-1
       STmeta(kk)%indST=kk
       STmeta(kk)%nIn=subtrees(kk)%tN%ST%NpraeST
    end do
    toNodes(:)=0
    do kk=nSubtrees,1,-1
       ind = subtrees(kk)%tN%ST%sizST
       call write_tree_to_array(subtrees(kk),STmeta(kk)%iStart,ind,&
            permNodes(STmeta(kk)%iStart:STmeta(kk)%iEnd), &
              toNodes(STmeta(kk)%iStart:STmeta(kk)%iEnd))
    end do
    ! ToDo: repair this part
    ! correct last toNode
    do kk=1,nSubtrees
       if (associated(subtrees(kk)%tN%post%tN)) then
          toNodes(STmeta(kk)%iEnd)=subtrees(kk)%tN%post%tN%ind
       end if
    end do
  end subroutine init_subtree_metadata

  subroutine distribute_subtree_meta(iBasin,nproc,nSubtrees,STmeta,toNodes,schedule,subtrees)
    implicit none
    integer(i4),                     intent(in)  :: iBasin
    integer(i4),                     intent(in)  :: nproc
    integer(i4),                     intent(in)  :: nSubtrees
    type(subtreeMeta), dimension(:), intent(in)  :: STmeta
    integer(i4),       dimension(:), intent(in)  :: toNodes
    type(processSchedule), dimension(:), intent(inout) :: schedule
    type(ptrTreeNode), dimension(:), intent(inout)     :: subtrees ! the array of
    ! local variables
    integer(i4) :: kk,ii,jj,iPerm,iproc,sizST,iST
    integer(i4) :: rank,ierror
    integer(i4), dimension(:,:), allocatable :: iSends ! the number and over all
                                                       ! length of arrays to be send to a process
    integer(i4), dimension(:), allocatable :: sendarray

    ! send number of subtrees and total number of tree nodes assigned
    ! to the processes to the corresponding process, so arrays can be
    ! allocated in the receiving subroutines
    allocate(iSends(nproc-1,2))
    do kk=1,nproc-1
       iSends(kk,1)=schedule(kk)%nTrees
       iSends(kk,2)=schedule(kk)%overallSize
       call MPI_Send(iSends(kk,:),2,MPI_INTEGER,kk,0,MPI_COMM_WORLD,ierror)
    end do
    ! send metadata of the subtrees to the nodes where they are assigned to
    ! size, identifying index corresponding to the subtree array and number of
    ! nodes, where data is flowing into the tree
    do kk=1,nproc-1
       do jj=1,schedule(kk)%nTrees
          iST=schedule(kk)%trees(jj)
          sizST=STmeta(iST)%iEnd+1-STmeta(iST)%iStart
          call MPI_Send(sizST,1,MPI_INTEGER,kk,1,MPI_COMM_WORLD,ierror)
          call MPI_Send(iST,1,MPI_INTEGER,kk,1,MPI_COMM_WORLD,ierror)
          call MPI_Send(STmeta(iST)%nIn,1,MPI_INTEGER,kk,1,MPI_COMM_WORLD,ierror)
       end do
    end do
    ! for each subtree send corresponding toNodes to the node. Move indices from
    ! toNodes, so they start from 0
    ! They can then be moved on the receiving process corresponding to
    ! its own offset
    do kk=1,nproc-1
       do jj=1,schedule(kk)%nTrees
          iST=schedule(kk)%trees(jj)
          sizST=STmeta(iST)%iEnd+1-STmeta(iST)%iStart
          allocate(sendarray(sizST))
          do ii=STmeta(iST)%iStart,STmeta(iST)%iEnd
             sendarray(ii-STmeta(iST)%iStart+1)=toNodes(ii)-STmeta(iST)%iStart+1
          end do
          call MPI_Send(sendarray(1:sizST),sizST,MPI_INTEGER,kk,2,MPI_COMM_WORLD,ierror)
          deallocate(sendarray)
       end do
    end do

    deallocate(iSends)
  end subroutine distribute_subtree_meta

  subroutine get_subtree_meta(iBasin,nproc,rank,STmeta,toNodes)
    implicit none
    integer(i4),               intent(in)                       :: iBasin
    integer(i4),               intent(in)                       :: nproc
    integer(i4),               intent(in)                       :: rank
    type(subtreeMeta), dimension(:), allocatable, intent(inout) :: STmeta
    integer(i4),       dimension(:), allocatable, intent(inout) :: toNodes
    ! local variables
    integer(i4) :: kk
    integer(i4), dimension(2) :: nDatasets ! number of incoming data sets
                                           ! total size of datasets
    integer(i4) :: sizST,indST
    integer(i4) :: nSubtrees, totSizeOfSubtrees
    integer(i4) :: ierror
    integer status(MPI_STATUS_SIZE)

    ! recieves number of subtrees and total number of tree nodes assigned to
    ! this process
    call MPI_Recv(nDatasets(:),2,MPI_INTEGER,0,0,MPI_COMM_WORLD,status,ierror)
    nSubtrees=nDatasets(1)
    totSizeOfSubtrees=nDatasets(2)
    allocate(STmeta(nSubtrees))
    allocate(toNodes(totSizeOfSubtrees))
    ! ToDo: case: less subtrees than processes
    call MPI_Recv(sizST,1,MPI_INTEGER,0,1,MPI_COMM_WORLD,status,ierror)
    STmeta(1)%iStart=1
    STmeta(1)%iEnd=sizST
    call MPI_Recv(STmeta(1)%indST,1,MPI_INTEGER,0,1,MPI_COMM_WORLD,status,ierror)
    call MPI_Recv(STmeta(1)%nIn,1,MPI_INTEGER,0,1,MPI_COMM_WORLD,status,ierror)
    do kk=2,nSubtrees
       call MPI_Recv(sizST,1,MPI_INTEGER,0,1,MPI_COMM_WORLD,status,ierror)
       STmeta(kk)%iStart=Stmeta(kk-1)%iEnd+1
       STmeta(kk)%iEnd=STmeta(kk)%iStart+sizST-1
       call MPI_Recv(STmeta(kk)%indST,1,MPI_INTEGER,0,1,MPI_COMM_WORLD,status,ierror)
       call MPI_Recv(STmeta(kk)%nIn,1,MPI_INTEGER,0,1,MPI_COMM_WORLD,status,ierror)
    end do

    do kk=1,nSubtrees
       sizST=STmeta(kk)%iEnd+1-STmeta(kk)%iStart
       call MPI_Recv(toNodes(STmeta(kk)%iStart:STmeta(kk)%iEnd),sizST,MPI_INTEGER,0,2,MPI_COMM_WORLD,status,ierror)
       ! ToDo: why not -1
       ! the toNodes have been moved, so they start from index 1, before sending
       ! now the last entry gets moved to the starting point of the subtree in the array
       ! toNodes(STmeta(kk)%iStart:STmeta(kk)%iEnd)=toNodes(STmeta(kk)%iStart:STmeta(kk)%iEnd)+STmeta(kk)%iStart
       toNodes(STmeta(kk)%iEnd)=toNodes(STmeta(kk)%iEnd)+STmeta(kk)%iStart
    end do

  end subroutine get_subtree_meta

END MODULE mo_HRD_subtree_meta

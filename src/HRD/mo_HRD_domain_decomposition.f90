!> \file mo_HRD_domain_decomposition.f91

!> \brief decomposes the basins into subdomains for parallel computing

!> \details The data processing with mRM is dependend on
!>      the (geological) river network structure. A river network
!>      in mHM is simplified to a (mathematical) tree.
!>      This model creates a tree structure out of the
!>      input data, cuts it into subtrees. The subtrees
!>      again can be interpreted as vertices in another tree,
!>      referred to as the _subtreetree_ throughout the code. Data
!>      corresponding to a subtree then get send to
!>      different (computing) nodes. The
!>      data gets processed there in routingorder, and a
!>      master process handles the communication. Afterwards
!>      it collects the data corresponding to the subtrees
!>      from the nodes.

!> \authors Maren Kaluza
!> \date June 2018

MODULE mo_HRD_domain_decomposition
  use mo_kind, only : i4, dp

  use mo_HRD_write
  use mo_HRD_types, only: ptrTreeNode, subtreeMeta, processSchedule, &
                  subtreeBuffer

  use mo_HRD_tree_init_and_destroy, only: tree_init_global, tree_init, tree_destroy, &
                                          forest_init, forest_destroy

  use mo_HRD_decompose, only: decompose

  use mo_HRD_schedule, only: create_schedule, create_schedule_hu, schedule_destroy

  use mo_HRD_subtree_meta, only: init_subtree_metadata, distribute_subtree_meta, &
                                 get_subtree_meta, init_subtree_buffer, destroy_subtree_buffer
  use mo_HRD_MPI_array_communication, only: distribute_array, collect_array, &
                                 get_array, send_array

  !$ use omp_lib,      only: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
  use mpi_f08

  ! Written Maren Kaluza, June 2018

  IMPLICIT NONE

  public :: domain_decomposition

  private

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         domain_decomposition

  !     PURPOSE
  !>        \brief decomposes the basins into subdomains for parallel computing
  !
  !>        \details decomposes basins into subdomains depending
  !>             on wether mRM is active or not, variable will then not be read in again.
  !
  !     INTENT(IN)
  !         None
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !         None
  !
  !     RESTRICTIONS
  !>        None
  !
  !     EXAMPLE
  !       None
  !
  !     LITERATURE
  !       None

  !     HISTORY
  !>        \author Maren Kaluza
  !>        \date June 2018
  !         Modified, June 2018 - Maren Kaluza, start of implementation

  subroutine domain_decomposition()
 ! USE mo_timer, ONLY : &
 !         timers_init, timer_start, timer_stop, timer_get, timer_clear  ! Timing of processes
    implicit none
    ! input variables

    ! local variables
    ! ToDo: Later to be moved to globally known variables
    ! ToDo: handle case without MPI, or only one process
    integer(i4)                                      :: iBasin
    integer(i4)                                      :: lowBound            ! a subdomain should include at least
                                                                            ! lowBound tree nodes
    integer(i4)                                      :: ind                 ! index of link/edge for subdomain
    type(ptrTreeNode),     dimension(:), allocatable :: roots               ! the root node of the tree structure
    type(ptrTreeNode)                                :: root                ! the root node of the tree structure
    type(ptrTreeNode),     dimension(:), allocatable :: subtrees            ! the array of
                                                                            ! subtrees in routing order
    integer(i4)                                      :: nNodes              ! the number of forest nodes in the original forest
    integer(i4)                                      :: nLinks              ! the number of edges in the forest
    integer(i4)                                      :: nBasins             ! nBasins
    integer(i4)                                      :: nSubtrees           ! number of subtrees in treedecomposition

    type(subtreeMeta),     dimension(:), allocatable :: STmeta
    type(subtreeBuffer),   dimension(:), allocatable :: STbuffer
    integer(i4),           dimension(:), allocatable :: permNodes           ! in-Nodes in routing order corresponding
    integer(i4),           dimension(:), allocatable :: toNodes             ! out-Nodes of corresponding (in-)Nodes
                                                                            ! to decomposition
    integer(i4),           dimension(:), allocatable :: fromNodes           ! in-Nodes of corresponding Nodes
                                                                            ! to decomposition
    type(processSchedule), dimension(:), allocatable :: schedule            ! knows for each process the number
                                                                            ! of subtrees assigned to it, the
                                                                            ! indices of the subtrees and the overall
                                                                            ! size of the subtrees
    integer(i4)                                      :: bufferLength        ! length of arrays buffered and send via MPI
    integer(i4)                                      :: nproc, rank, ierror
    type(MPI_Comm)                                   :: comm                ! MPI communicator

    ! for testing purposes
    integer(i4),           dimension(:), allocatable :: testarray
    integer(i4)                                      :: kk, mes, ll
    ! ToDo: remove later
    integer(i4)                                      :: iTimer              ! Current timer number

#ifdef MRM2MHM
    iBasin = 1
    ! create a dublicate communicator of MPI_COMM_WORLD
    ! ToDo: Later, have an optional argument with another communicator
    call MPI_Comm_dup(MPI_COMM_WORLD, comm, ierror)
    ! find number of processes nproc
    call MPI_Comm_size(comm, nproc, ierror)
    ! find the number the process is referred to, called rank
    call MPI_Comm_rank(comm, rank, ierror)
  !  do lowBound = 10,210,20
  !  do nproc = 2,6,2!96,2
    bufferLength = 2
    if (rank .eq. 0) then
      write(*,*) 'the domain decomposition with mRM gets implemented now...'
      ! this subroutine is called by all processes, but only the
      ! (perhaps) master process knows global variables
      ! so these should not be loaded from other modules inside this subroutine
      call get_number_of_basins(iBasin, nBasins)
      call get_L11_information(iBasin, nLinks, nNodes, toNodes, fromNodes, permNodes)
      call init_testarray(nNodes-1, testarray)

      lowBound = 3
      ! In this subroutine the tree structure gets initialized for
      ! the flownetwork of the iBasin-th basin.
      ! In each tree node the size of the smallest subtree
      ! larger than lowBound is saved, so cutting off
      ! subtrees can be done efficiently.
      ! Root is the root node of the tree. All other tree nodes
      ! can be accessed through it.
      ! call tree_init_global(iBasin, lowBound, root)
      call forest_init(nLinks, nNodes, toNodes, lowBound, roots, fromNodes=fromNodes, perm=permNodes)
     ! call write_tree(roots(1), lowBound)
      deallocate(toNodes, permNodes, fromNodes)

      ! ToDo: that's possibly too much, but maybe more efficient than reallocating?
      allocate(subtrees(nNodes))
      ! this subroutine cuts down the tree into a subtreetree and
      ! writes each subtreetree node into an array (subtrees) in routing order
      call decompose(lowBound, roots, subtrees(:), nSubtrees)

      ! call write_domain_decomposition(root)
      allocate(STmeta(nSubtrees), schedule(nproc-1))
      ! create schedule:
      ! to each process in the array schedule the number of trees, the
      ! indices of the trees and the over all size is assigned
      call create_schedule_hu(nSubtrees, nproc, subtrees(:), schedule)
      !call schedule_destroy(iBasin,schedule)
      !allocate(schedule(nproc-1))
      ! call create_schedule(nSubtrees,subtrees,schedule)
      ! call write_graphviz_output_forest(roots)
      ! A subtree data structrure makes communication between the subtrees
      ! much easier for the master. Processing the data is more efficient
      ! with array, so everything gets written into a nice array in
      ! routing order. Therefore we need a permutation array permNodes
      ! which is at the same time fromNodes, and a toNode array
      call init_subtree_metadata(iBasin, nNodes, subtrees(:), STmeta, permNodes, toNodes)

      ! sends the meta data from master process to all the others
      call distribute_subtree_meta(iBasin, nproc, comm, nSubtrees, STmeta, toNodes, schedule, subtrees(:))
      ! - sends data (testarray) corresponding to subtrees to nodes
      ! - collects processed data from roots from subtrees and sends this
      !   data to corresponding leaves in connected subtrees
      ! - collects the data in the end
      iTimer = 26
     ! do ll = 1, 10
     ! call timer_start(itimer)
     ! do kk = 1, 1
      call routing(iBasin, nproc, rank, comm, bufferLength, subtrees(:), nSubtrees, STmeta, permNodes, schedule, testarray)
     ! end do
     ! call timer_stop(itimer)
     ! write(*,*) timer_get(itimer), 'seconds.', nSubtrees, 'subtrees, lowBound:', lowBound
     ! write(*,*) nproc, timer_get(itimer), nSubtrees, lowBound
     ! call timer_clear(itimer)
     ! end do
    !  call write_forest_with_array(roots, lowBound,testarray)

      call schedule_destroy(schedule)
      deallocate(STmeta)

      call forest_destroy(roots)
      deallocate(subtrees, permNodes, toNodes)

      call destroy_testarray(testarray)
    else if (rank < nproc) then
      ! lowBound=3
      ! all other processes receive meta data for their
      ! individual subtrees from the master process
      call get_subtree_meta(iBasin, comm, STmeta, toNodes)
      call subtree_init(lowBound, STmeta, toNodes, subtrees)
      call init_subtree_buffer(STmeta, bufferLength, STbuffer)
      ! - receives data corresponding to an array and assigned
      !   subtrees
      ! - receive input data from connected subtrees
      ! - processes data
      ! - sends root data to master
      ! - send data to master
    !  do ll=1,10
    !  do kk=1,1
      call subtree_routing(iBasin, nproc, rank, comm, bufferLength, toNodes, subtrees, STmeta, STbuffer, testarray)
    !  end do
    !  end do
      call destroy_subtree_buffer(STbuffer)
      do kk = 1, size(subtrees)
        call tree_destroy(subtrees(kk))
      end do
      deallocate(subtrees)
      ! ToDo: deallocating might be nicer on the same level, so
      ! either outside or an extra subroutine?
      deallocate(STmeta, toNodes)
    endif
  !  end do
  !  end do
#else
    if (rank .eq. 0) then
      write(*,*) 'the domain decomposition without mRM is not implemented yet'
    else
      write(*,*) 'need to have something to eat'
    endif
#endif
  end subroutine domain_decomposition

  ! - sends data (array) corresponding to subtrees to nodes
  ! - collects processed data from roots from subtrees and sends this
  !   data to corresponding leaves in connected subtrees
  ! - collects the data in the end
  subroutine routing(iBasin, nproc, rank, comm, bufferLength, subtrees, nSubtrees, STmeta, permNodes, schedule, array)
    implicit none
    integer(i4),                         intent(in)    :: iBasin
    integer(i4),                         intent(in)    :: nproc
    integer(i4),                         intent(in)    :: rank
    type(MPI_Comm),                      intent(in)    :: comm
    integer(i4),                         intent(in)    :: bufferLength
    type(ptrTreeNode),     dimension(:), intent(in)    :: subtrees ! the array of
    integer(i4),                         intent(in)    :: nSubtrees
    type(subtreeMeta),     dimension(:), intent(in)    :: STmeta
    integer(i4),           dimension(:), intent(in)    :: permNodes
    type(processSchedule), dimension(:), intent(in)    :: schedule
    integer(i4),           dimension(:), intent(inout) :: array
    ! local
    integer(i4)                            :: kk,iproc,next,ind, indST
    integer(i4), dimension(bufferLength+1) :: buffer
    integer(i4)                            :: ierror
    type(MPI_Status)                       :: status

    call distribute_array(iBasin, nproc, rank, comm, nSubtrees, STmeta, permNodes, schedule, array)

    do kk = 1, nSubtrees
      ! for each subtree the master process 0 gets the data of
      ! the root tree node
      call MPI_Recv(buffer, bufferLength+1, MPI_INTEGER, MPI_ANY_SOURCE, 0, comm, status, ierror)
      indST = buffer(bufferLength+1)
      ! if the root node of the subtree has a parent
      if (associated(subtrees(indST)%tN%post%tN)) then
        ! next becomes the index of the parent subtree root node (the subtree
        ! that gets data)
        next = subtrees(indST)%tN%ST%postST%tN%ST%indST

        ! ind becomes the index of the parent of the root node (the tree node
        ! in the subtree that gets data, note: there is an offset)
        ind = subtrees(indST)%tN%post%tN%ind
        ! iproc is the process the parent subtree is assigned to
        iproc = subtrees(next)%tN%ST%sched(1)
        ! the index, where the value will be added to (we remove the offset)
        buffer(bufferLength+1) = ind - STmeta(next)%iStart
        call MPI_Send(buffer, bufferLength+1, MPI_INTEGER, iproc, next, comm, ierror)
        ! write(*,*) 'master sent', buffer(1), 'to process',iproc,'tree',next !, 'ind_diff', ind-STmeta(next)%iStart
      end if
    end do

    array(:) = 0
    call collect_array(iBasin, nproc, rank, comm, nSubtrees, STmeta, permNodes, schedule, array)
  end subroutine

  subroutine subtree_init(lowBound, STmeta, toNodes, subtrees)
    implicit none
    integer(i4),                                  intent(in)    :: lowBound
    type(subtreeMeta), dimension(:),              intent(in)    :: STmeta
    integer(i4),       dimension(:),              intent(in)    :: toNodes
    type(ptrTreeNode), dimension(:), allocatable, intent(inout) :: subtrees ! the array of
    ! local
    integer(i4) :: nST ! number of subtrees scheduled on this computational node
    integer(i4) :: kk
    integer(i4) :: sizST,iStart,iEnd

    nST = size(STmeta)

    allocate(subtrees(nST))

    do kk = 1, nST
      iStart = STmeta(kk)%iStart
      iEnd   = STmeta(kk)%iEnd
      sizST  = iEnd - iStart + 1
   !   write(*,*) '---------------------------',sizST, iStart, iEnd, STmeta(kk)%nIn
      call tree_init(sizST-1, sizST, toNodes(iStart:iEnd), lowBound, subtrees(kk))
     ! call write_subtree(subtrees(kk), lowBound)
    end do
    
  end subroutine subtree_init

  ! - receives data corresponding to an array and assigned
  !   subtrees
  ! - receive input data from connected subtrees
  ! - processes data
  ! - sends root data to master
  ! - send data to master
  subroutine subtree_routing(iBasin, nproc, rank, comm, bufferLength, toNodes, subtrees, STmeta, STbuffer, array)
    implicit none
    integer(i4),                                    intent(in)    :: iBasin
    integer(i4),                                    intent(in)    :: nproc
    integer(i4),                                    intent(in)    :: rank
    type(MPI_Comm),                                 intent(in)    :: comm
    integer(i4),                                    intent(in)    :: bufferLength
    integer(i4),         dimension(:),              intent(in)    :: toNodes
    type(ptrTreeNode),   dimension(:),              intent(in)    :: subtrees ! the array of
    type(subtreeMeta),   dimension(:),              intent(in)    :: STmeta
    type(subtreeBuffer), dimension(:),              intent(inout) :: STbuffer
    integer(i4),         dimension(:), allocatable, intent(inout) :: array
    ! local
    integer(i4)                                    :: nST ! number of subtrees scheduled
                                                          ! on this computational node
    integer(i4)                                    :: kk,jj,next,iStart,iEnd, ii, nIn
    integer(i4)                                    :: ierror

    call get_array(iBasin, nproc, rank, comm, STmeta, array)

    ! number of subtrees assigned to that process
    nST = size(STmeta)
    do kk = 1, nST
      do jj = 1, STmeta(kk)%nIn
        ! for every subtree for every leave, where data comes in, receive this
        ! data. Via STmeta(kk)%indST we write the right incoming data to the
        ! matching buffer. In the last entry of the buffer the index of the
        ! leaves is trasmitted, where the data goes into
        call MPI_IRecv(STbuffer(kk)%buffer(:,jj), bufferLength+1, MPI_INTEGER, 0, STmeta(kk)%indST, comm, STbuffer(kk)%requests(jj), ierror)
      end do
    end do
    ! now we run over the subtrees on the computing node in routing order
    do kk = 1, nST
      nIn = STmeta(kk)%nIn
      do jj = 1, nIn
        ! we can start calculations in that subtree, when als its inflows are
        ! reveived
        call MPI_Wait(STbuffer(kk)%requests(jj), STbuffer(kk)%statuses(jj), ierror)
      !   write(*,*) 'process', rank, 'tree', STmeta(kk)%indST
      end do
      iStart = STmeta(kk)%iStart
      iEnd   = STmeta(kk)%iEnd
      do ii = 1, bufferLength
        do jj = 1, STmeta(kk)%nIn
          ! shift the index with respect to the start of the subtree
          next = STbuffer(kk)%buffer(bufferLength+1, jj) + STmeta(kk)%iStart
          array(next) = array(next) + STbuffer(kk)%buffer(ii,jj)
        end do
        !call nodeinternal_routing_array(toNodes(iStart:iEnd),array(iStart:iEnd))
        !!$OMP parallel num_threads(jj) private(rank) shared(testarray)
        !$OMP parallel private(rank) shared(array,subtrees)
        !$OMP single
        call nodeinternal_routing(subtrees(kk), array(iStart:iEnd))
     !   call nodeinternal_routing_array(toNodes(iStart:iEnd),array(iStart:iEnd))
        !$OMP end single
        !$OMP barrier
        !$OMP end parallel
        STbuffer(kk)%buffer(ii, STmeta(kk)%nIn+1) = array(STmeta(kk)%iEnd)
      end do
      STbuffer(kk)%buffer(bufferLength+1, STmeta(kk)%nIn+1)=STmeta(kk)%indST
      ! send the outflow to the master process
      call MPI_Send(STbuffer(kk)%buffer(:, STmeta(kk)%nIn+1), bufferLength+1, MPI_INTEGER, 0, 0, comm, ierror)
    end do

    call send_array(iBasin, nproc, rank, comm, STmeta, array)
    deallocate(array)
  end subroutine subtree_routing

  subroutine nodeinternal_routing_array(toNodes, array)
    implicit none
    integer(i4),       dimension(:),              intent(in)    :: toNodes
    integer(i4),       dimension(:),              intent(inout) :: array
    ! local
    integer(i4) :: jj
    integer(i4) :: ii

    do jj = 1, size(toNodes)-1
      array(toNodes(jj)) = array(toNodes(jj)) + array(jj)
    end do
  end subroutine nodeinternal_routing_array

  recursive subroutine nodeinternal_routing(root, array)
    implicit none
    type(ptrTreeNode),                            intent(in)    :: root ! the array of
    integer(i4),       dimension(:),              intent(inout) :: array
    ! local
    integer(i4) :: jj
    integer(i4) :: ii, ithread
    integer(i4) :: tNode

    do jj = 1, root%tN%Nprae
      !$OMP task shared(root,array)
      call nodeinternal_routing(root%tN%prae(jj), array)
      !$OMP end task
    end do
    !$OMP taskwait
    if (associated(root%tN%post%tN)) then
      tNode = root%tN%post%tN%origind
      !$OMP critical
     ! !$ ithread=omp_get_thread_num()
     ! !$ write(*,*) ithread
      array(tNode) = array(tNode) + array(root%tN%origind)
      !$OMP end critical
    end if
  end subroutine nodeinternal_routing

  subroutine get_number_of_basins(iBasin, numBasins)
    use mo_common_variables, only : &
            nBasins
    implicit none
    integer, intent(in)    :: iBasin
    integer, intent(inout) :: numBasins
    numBasins = nBasins
  end subroutine get_number_of_basins

  subroutine get_L11_information(iBasin, nLinks, nNodes, toNodes, fromNodes, permNodes)
    use mo_mrm_global_variables, only : &
            level11,      & ! IN: for number of nCells
            L11_fromN,    & ! IN: for an edge this is the incoming tree node
            L11_toN,      & ! IN: for an edge this is the outgoing tree node
            L11_netPerm,  & ! IN: network routing order
            L11_nOutlets    ! IN: number of nodes minus number of outlets is the
                            !        number of links
    implicit none
    integer(i4),                            intent(in)    :: iBasin
    integer(i4),                            intent(out)   :: nLinks
    integer(i4),                            intent(out)   :: nNodes
    integer(i4), dimension(:), allocatable, intent(out)   :: toNodes
    integer(i4), dimension(:), allocatable, intent(out)   :: fromNodes
    integer(i4), dimension(:), allocatable, intent(out)   :: permNodes

    nNodes = level11(iBasin)%nCells
    nLinks = nNodes - L11_nOutlets(iBasin)
    allocate(toNodes(nLinks))
    toNodes(:)   = L11_toN(1:nLinks)
    allocate(fromNodes(nLinks))
    fromNodes(:) = L11_fromN(1:nLinks)
    allocate(permNodes(nLinks))
    permNodes(:) = L11_netPerm(1:nLinks)

  end subroutine get_L11_information

  subroutine init_testarray(nLinks, testarray)
    implicit none
    integer(i4),               intent(in)                  :: nLinks
    integer(i4), dimension(:), allocatable, intent(inout)  :: testarray
    integer(i4) :: kk ! loop variable to run over all edges/links

    allocate(testarray(nLinks+1))
    do kk = 1, nLinks+1
      !testarray(kk)=kk
      testarray(kk) = 1
    end do
  end subroutine init_testarray

  subroutine destroy_testarray(testarray)
    implicit none
    integer(i4), dimension(:), allocatable, intent(inout)  :: testarray

    deallocate(testarray)
  end subroutine destroy_testarray

  ! Not needed, can be deleted soon
  subroutine find_L11_fAcc(iBasin,L11_fAcc)
    use mo_mrm_global_variables, only : &
            level11,      & ! IN: for number of nCells
            L11_fromN,    & ! IN: for an edge this is the incoming tree node
            L11_toN,      & ! IN: for an edge this is the outgoing tree node
            L11_label,    & ! IN: label on edge, Id [0='', 1=HeadWater, 2=Sink]
            L11_rOrder,   & ! IN: network routing order
            L11_netPerm     ! IN: network routing order
    implicit none
    integer(i4),               intent(in)    :: iBasin
    integer(i4), dimension(:), intent(inout) :: L11_fAcc

    ! local variables
    integer(i4)                                  :: kk           ! loop variable to run over all edges/links
    integer(i4)                                  :: i            ! index of edge in routing order
    integer(i4)                                  :: iNode, tNode ! in and to tree node of corresponding edge
    integer(i4)                                  :: nLinks       ! number of edges
    type(ptrTreeNode), dimension(:), allocatable :: tree
    integer(i4),       dimension(:), allocatable :: Nprae

    nLinks = level11(iBasin)%ncells - 1
    ! ToDo: Why is nLinks != size(L11_netPerm(:))? For trees it is m=n-1
    do kk = 1, nLinks+1
      L11_fAcc(kk) = 1
    enddo
    do kk = 1, nLinks
      ! get LINK routing order -> i
      i     = L11_netPerm(kk)
      iNode = L11_fromN(i)
      tNode = L11_toN(i)
      L11_fAcc(tNode) = L11_fAcc(tNode) + L11_fAcc(iNode)
    end do
  end subroutine find_L11_fAcc
  
END MODULE mo_HRD_domain_decomposition

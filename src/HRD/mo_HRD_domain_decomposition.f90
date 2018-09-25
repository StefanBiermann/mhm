!> \file mo_HRD_domain_decomposition.f90

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
  use mo_HRD_types, only: ptrTreeNode, subtreeMeta, processSchedule

  use mo_HRD_tree_init_and_destroy

  use mo_HRD_decompose, only: decompose

  use mo_HRD_schedule, only: create_schedule, create_schedule_hu

  use mo_HRD_subtree_meta_init_and_destroy
  !$ use omp_lib,      only: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
  use mpi

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
    implicit none
    ! input variables

    ! local variables
    ! ToDo: Later to be moved to globally known variables
    ! ToDo: handle case without MPI, or only one process
    ! ToDo: handle more outlets
    integer(i4) :: iBasin
    integer(i4) :: lowBound,uppBound ! a subdomain should include at least
                                     ! lowBound tree nodes and at most uppBound
    integer(i4) :: ind               ! index of link/edge for subdomain
    type(ptrTreeNode) :: root        ! the root node of the tree structure
    type(ptrTreeNode), dimension(:), allocatable :: subtrees ! the array of
                                     ! subtrees in routing order
    integer(i4)       :: nNodes      ! the number of tree nodes in the original tree
    integer(i4)       :: nBasins     ! nBasins
    integer(i4)       :: nSubtrees   ! number of subtrees in treedecomposition

    type(subtreeMeta), dimension(:), allocatable :: STmeta
    integer(i4) , dimension(:), allocatable      :: permNodes   ! in-Nodes in routing order corresponding
    integer(i4) , dimension(:), allocatable      :: toNodes     ! out-Nodes of corresponding in-Nodes
                                                                ! to decomposition
    type(processSchedule), dimension(:), allocatable :: schedule! knows for each process the number
                                                                ! of subtrees assigned to it, the
                                                                ! indices of the subtrees and the overall
                                                                ! size of the subtrees
    integer(i4) :: nproc,rank,ierror

    ! for testing purposes
    integer(i4), dimension(:), allocatable :: testarray
    integer(i4) :: kk,mes

#ifdef MRM2MHM
    iBasin=1
    ! find number of processes nproc
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierror)
    ! find the number the process is referred to, called rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)
    if (rank .eq. 0) then
       write(*,*) 'the domain decomposition with mRM gets implemented now...'
       ! this subroutine is called by all processes, but only the
       ! (perhaps) master process knows global variables
       ! so these should not be loaded from other modules inside this subroutine
       call get_number_of_basins_and_nodes(iBasin,nNodes,nBasins)
       call init_testarray(iBasin,testarray)

       lowBound=60
       uppBound=5
       ! In this subroutine the tree structure gets initialized for
       ! the flownetwork of the iBasin-th basin.
       ! In each tree node the size of the smallest subtree
       ! larger than lowBound is saved, so cutting of
       ! subtrees can be done efficiently.
       ! Root is the root node of the tree. All other tree nodes
       ! can be accessed through it.
       call init_tree(iBasin, lowBound, root)

       ! ToDo: that's possibly a bit too much, but maybe more efficient than reallocating?
       allocate(subtrees(nNodes/lowBound+1))
       ! this subroutine cuts down the tree into a subtreetree and
       ! writes each subtreetree node into an array (subtrees) in routing order
       call decompose(iBasin,lowBound,root,subtrees,nSubtrees)

       ! call write_domain_decomposition(root)
       allocate(STmeta(nSubtrees),permNodes(nNodes),toNodes(nNodes),schedule(nproc-1))
       ! create schedule:
       ! to each process in the array schedule the number of trees, the
       ! indices of the trees and the over all size is assigned
       call create_schedule_hu(iBasin,nSubtrees,subtrees,schedule)
       ! call write_graphviz_output(root)
       !call schedule_destroy(iBasin,schedule)
       !allocate(schedule(nproc-1))
       ! call create_schedule(iBasin,nSubtrees,subtrees,schedule)
       ! A subtree data structrure makes communication between the subtrees
       ! much easier for the master. Processing the data is more efficient
       ! with array, so everything gets written into a nice array in
       ! routing order. Therefore we need a permutation array permNodes
       ! which is at the same time fromNodes, and a toNode array
       call init_subtree_metadata(iBasin,subtrees,STmeta,permNodes,toNodes)

       ! sends the meta data from master process to all the others
       call distribute_subtree_meta(iBasin,nSubtrees,STmeta,toNodes,schedule,subtrees)
       ! - sends data (testarray) corresponding to subtrees to nodes
       ! - collects processed data from roots from subtrees and sends this
       !   data to corresponding leaves in connected subtrees
       ! - collects the data in the end
       call routing(iBasin,subtrees,nSubtrees,STmeta,permNodes,schedule,testarray)
      ! call write_tree_with_array(root, lowBound,testarray)

       call schedule_destroy(iBasin,schedule)
       deallocate(STmeta)

       call tree_destroy(iBasin,root)
       deallocate(subtrees,permNodes,toNodes)

       call destroy_testarray(testarray)
    else
       ! all other processes receive meta data for their
       ! individual subtrees from the master process
       call get_subtree_meta(iBasin,STmeta,toNodes)
       ! - receives data corresponding to an array and assigned
       !   subtrees
       ! - receive input data from connected subtrees
       ! - processes data
       ! - sends root data to master
       ! - send data to master
       call subtree_routing(iBasin,toNodes,STmeta,testarray)
       ! ToDo: deallocating might be nicer on the same level, so
       ! either outside or an extra subroutine?
       deallocate(STmeta,toNodes)
    endif
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
  subroutine routing(iBasin,subtrees,nSubtrees,STmeta,permNodes,schedule,array)
    implicit none
    integer(i4),                     intent(in)     :: iBasin
    type(ptrTreeNode), dimension(:), intent(in)     :: subtrees ! the array of
    integer(i4),                     intent(in)     :: nSubtrees
    type(subtreeMeta), dimension(:), intent(in)     :: STmeta
    integer(i4),       dimension(:), intent(in)     :: permNodes
    type(processSchedule), dimension(:), intent(in) :: schedule
    integer(i4),       dimension(:), intent(inout)  :: array
    ! local
    integer(i4) :: kk,iproc,next,ind
    integer(i4), dimension(2) :: value_ind
    integer(i4) :: nproc,rank,ierror
    integer status(MPI_STATUS_SIZE)

    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierror)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)
    call distribute_array(iBasin,nSubtrees,STmeta,permNodes,schedule,array)
    
    do kk=1,nSubtrees
       iproc=subtrees(kk)%tN%ST%sched(1)
       call MPI_Recv(value_ind,2,MPI_INTEGER,iproc,7,MPI_COMM_WORLD,status,ierror)
       if (associated(subtrees(value_ind(2))%tN%post%tN)) then
          next=subtrees(value_ind(2))%tN%ST%postST%tN%ST%indST

          ind=subtrees(value_ind(2))%tN%post%tN%ind
          iproc=subtrees(next)%tN%ST%sched(1)
          call MPI_Send([value_ind(1),ind-STmeta(next)%iStart],2,MPI_INTEGER,iproc,next,MPI_COMM_WORLD,ierror)
       !   write(*,*) 'master sent', value_ind(1), 'to process',iproc,'tree',next, 'ind_diff', ind-STmeta(next)%iStart
       end if
    end do

    array(:)=0
    call collect_array(iBasin,nSubtrees,STmeta,permNodes,schedule,array)
  end subroutine

  ! - receives data corresponding to an array and assigned
  !   subtrees
  ! - receive input data from connected subtrees
  ! - processes data
  ! - sends root data to master
  ! - send data to master
  subroutine subtree_routing(iBasin,toNodes,STmeta,array)
    implicit none
    integer(i4),                                  intent(in)    :: iBasin
    integer(i4),       dimension(:),              intent(in)    :: toNodes
    type(subtreeMeta), dimension(:),              intent(in)    :: STmeta
    integer(i4),       dimension(:), allocatable, intent(inout) :: array
    ! local
    integer(i4) :: nST ! number of subtrees scheduled on this computational node
    integer(i4) :: kk,jj,next
    integer(i4), dimension(2) :: value_ind
    integer(i4) :: nproc,rank,ierror
    integer status(MPI_STATUS_SIZE)

    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierror)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)
    call get_array(iBasin,STmeta,array)

    nST=size(STmeta)
    do kk=1,nST
       do jj=1,STmeta(kk)%nIn
          call MPI_Recv(value_ind,2,MPI_INTEGER,0,STmeta(kk)%indST,MPI_COMM_WORLD,status,ierror)
          next=value_ind(2)+STmeta(kk)%iStart
          array(next)=array(next)+value_ind(1)
   !       write(*,*) 'process',rank, 'tree', STmeta(kk)%indST, 'with indices', &
   !               STmeta(kk)%iStart,'-',STmeta(kk)%iEnd ,&
   !               'gets', value_ind(1), 'for ind', next, 'ind_diff', value_ind(2)
       end do
       call nodeinternal_routing(kk,toNodes,STmeta,array)
       call MPI_Send([array(STmeta(kk)%iEnd),STmeta(kk)%indST],2,MPI_INTEGER,0,7,MPI_COMM_WORLD,ierror)
    end do

    call send_array(iBasin,STmeta,array)
    deallocate(array)
  end subroutine subtree_routing

  subroutine nodeinternal_routing(kk,toNodes,STmeta,array)
    implicit none
    integer(i4),                                  intent(in)    :: kk
    integer(i4),       dimension(:),              intent(in)    :: toNodes
    type(subtreeMeta), dimension(:),              intent(in)    :: STmeta
    integer(i4),       dimension(:),              intent(inout) :: array
    ! local
    integer(i4) :: jj
    do jj=STmeta(kk)%iStart,STmeta(kk)%iEnd-1
       array(toNodes(jj))=array(toNodes(jj))+array(jj)
    end do
  end subroutine nodeinternal_routing

  subroutine distribute_subtree_meta(iBasin,nSubtrees,STmeta,toNodes,schedule,subtrees)
    implicit none
    integer(i4),                     intent(in)  :: iBasin
    integer(i4),                     intent(in)  :: nSubtrees
    type(subtreeMeta), dimension(:), intent(in)  :: STmeta
    integer(i4),       dimension(:), intent(in)  :: toNodes
    type(processSchedule), dimension(:), intent(inout) :: schedule
    type(ptrTreeNode), dimension(:), intent(inout)     :: subtrees ! the array of
    ! local variables
    integer(i4) :: kk,ii,jj,iPerm,iproc,sizST,iST
    integer(i4) :: nproc,rank,ierror
    integer(i4), dimension(:,:), allocatable :: iSends ! the number and over all
                                                       ! length of arrays to be send to a process
    integer(i4), dimension(:), allocatable :: sendarray

    ! find number of processes nproc
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierror)
    ! find the number the process is referred to, called rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)
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
             sendarray(ii-STmeta(iST)%iStart+1)=toNodes(ii)-STmeta(iST)%iStart
          end do
          call MPI_Send(sendarray(1:sizST),sizST,MPI_INTEGER,kk,2,MPI_COMM_WORLD,ierror)
          deallocate(sendarray)
       end do
    end do

    deallocate(iSends)
  end subroutine distribute_subtree_meta

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

  subroutine get_subtree_meta(iBasin,STmeta,toNodes)
    implicit none
    integer(i4),               intent(in)                       :: iBasin
    type(subtreeMeta), dimension(:), allocatable, intent(inout) :: STmeta
    integer(i4),       dimension(:), allocatable, intent(inout) :: toNodes
    ! local variables
    integer(i4) :: kk
    integer(i4), dimension(2) :: nDatasets ! number of incoming data sets
                                           ! total size of datasets
    integer(i4) :: sizST,indST
    integer(i4) :: nproc,rank,ierror
    integer status(MPI_STATUS_SIZE)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierror)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

    ! recieves number of subtrees and total number of tree nodes assigned to
    ! this process
    call MPI_Recv(nDatasets(:),2,MPI_INTEGER,0,0,MPI_COMM_WORLD,status,ierror)
    allocate(STmeta(nDatasets(1)))
    allocate(toNodes(nDatasets(2)))
    ! ToDo: case: less subtrees than processes
    call MPI_Recv(sizST,1,MPI_INTEGER,0,1,MPI_COMM_WORLD,status,ierror)
    STmeta(1)%iStart=1
    STmeta(1)%iEnd=sizST
    call MPI_Recv(STmeta(1)%indST,1,MPI_INTEGER,0,1,MPI_COMM_WORLD,status,ierror)
    call MPI_Recv(STmeta(1)%nIn,1,MPI_INTEGER,0,1,MPI_COMM_WORLD,status,ierror)
    do kk=2,nDatasets(1)
       call MPI_Recv(sizST,1,MPI_INTEGER,0,1,MPI_COMM_WORLD,status,ierror)
       STmeta(kk)%iStart=Stmeta(kk-1)%iEnd+1
       STmeta(kk)%iEnd=STmeta(kk)%iStart+sizST-1
       call MPI_Recv(STmeta(kk)%indST,1,MPI_INTEGER,0,1,MPI_COMM_WORLD,status,ierror)
       call MPI_Recv(STmeta(kk)%nIn,1,MPI_INTEGER,0,1,MPI_COMM_WORLD,status,ierror)
    end do

    do kk=1,nDatasets(1)
       sizST=STmeta(kk)%iEnd+1-STmeta(kk)%iStart
       call MPI_Recv(toNodes(STmeta(kk)%iStart:STmeta(kk)%iEnd),sizST,MPI_INTEGER,0,2,MPI_COMM_WORLD,status,ierror)
       ! ToDo: why not -1
       ! the toNodes have been moved, so they start from index 1, before sending
       ! now they get moved to the starting point of the subtree in the array
       toNodes(STmeta(kk)%iStart:STmeta(kk)%iEnd)=toNodes(STmeta(kk)%iStart:STmeta(kk)%iEnd)+STmeta(kk)%iStart
    end do

  end subroutine get_subtree_meta

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

  ! destroy schedule variable
  subroutine schedule_destroy(iBasin,schedule)
    implicit none
    integer(i4),                     intent(in) :: iBasin
    type(processSchedule), dimension(:), allocatable, intent(inout) :: schedule
    ! local variables
    integer(i4) :: kk

    do kk=1,size(schedule)
       deallocate(schedule(kk)%trees)
    end do
    deallocate(schedule)
  end subroutine schedule_destroy

  subroutine get_number_of_basins_and_nodes(iBasin,nNodes,numBasins)
    use mo_mrm_global_variables, only : &
            level11        ! IN: for number of nCells
    use mo_common_variables, only : &
            nBasins
    implicit none
    integer, intent(in)    :: iBasin
    integer, intent(inout) :: nNodes
    integer, intent(inout) :: numBasins
    nNodes=level11(iBasin)%ncells
    numBasins=nBasins
  end subroutine get_number_of_basins_and_nodes

  subroutine init_testarray(iBasin,testarray)
    use mo_mrm_global_variables, only : &
            level11       ! IN: for number of nCells
    implicit none
    integer(i4),               intent(in)                  :: iBasin
    integer(i4), dimension(:), allocatable, intent(inout)  :: testarray
    integer(i4) :: kk ! loop variable to run over all edges/links
    integer(i4) :: nLinks ! number of edges

   nLinks=level11(iBasin)%ncells - 1

   allocate(testarray(nLinks+1))
   do kk=1,nLinks+1
      !testarray(kk)=kk
      testarray(kk)=1
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
    integer(i4) :: kk ! loop variable to run over all edges/links
    integer(i4) :: i ! index of edge in routing order
    integer(i4) :: iNode, tNode ! in and to tree node of corresponding edge
    integer(i4) :: nLinks ! number of edges
    type(ptrTreeNode), dimension(:), allocatable :: tree
    integer(i4), dimension(:), allocatable :: Nprae

   nLinks=level11(iBasin)%ncells - 1
   ! ToDo: Why is nLinks != size(L11_netPerm(:))? For trees it is m=n-1
   do kk = 1, nLinks+1
      L11_fAcc(kk)=1
   enddo
   do kk = 1, nLinks
      ! get LINK routing order -> i
      i = L11_netPerm(kk)
      iNode = L11_fromN(i)
      tNode = L11_toN(i)
      L11_fAcc(tNode)=L11_fAcc(tNode)+L11_fAcc(iNode)
   end do
  end subroutine find_L11_fAcc

END MODULE mo_HRD_domain_decomposition

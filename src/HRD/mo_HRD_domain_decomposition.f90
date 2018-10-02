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

  use mo_HRD_tree_init_and_destroy, only: tree_init_global, tree_init, tree_destroy

  use mo_HRD_decompose, only: decompose

  use mo_HRD_schedule, only: create_schedule, create_schedule_hu, schedule_destroy

  use mo_HRD_subtree_meta, only: init_subtree_metadata, distribute_subtree_meta, &
                                 get_subtree_meta
  use mo_HRD_MPI_array_communication, only: distribute_array, collect_array, &
                                 get_array, send_array

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
    integer(i4) , dimension(:), allocatable      :: toNodes     ! out-Nodes of corresponding (in-)Nodes
                                                                ! to decomposition
    integer(i4) , dimension(:), allocatable      :: fromNodes   ! in-Nodes of corresponding Nodes
                                                                ! to decomposition
    type(processSchedule), dimension(:), allocatable :: schedule! knows for each process the number
                                                                ! of subtrees assigned to it, the
                                                                ! indices of the subtrees and the overall
                                                                ! size of the subtrees
    integer(i4) :: bufferLength                                 ! length of arrays buffered and send via MPI
    integer(i4) :: nproc,rank,ierror

    ! for testing purposes
    integer(i4), dimension(:), allocatable :: testarray
    integer(i4) :: kk,mes,ll

#ifdef MRM2MHM
    iBasin=1
    ! find number of processes nproc
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierror)
    ! find the number the process is referred to, called rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)
    bufferLength=2
    if (rank .eq. 0) then
       write(*,*) 'the domain decomposition with mRM gets implemented now...'
       ! this subroutine is called by all processes, but only the
       ! (perhaps) master process knows global variables
       ! so these should not be loaded from other modules inside this subroutine
       call get_number_of_basins_and_nodes(iBasin,nNodes,nBasins)
       call get_L11_information(iBasin, toNodes, fromNodes, permNodes)
       call init_testarray(iBasin,testarray)

       lowBound=50
       uppBound=5
       ! In this subroutine the tree structure gets initialized for
       ! the flownetwork of the iBasin-th basin.
       ! In each tree node the size of the smallest subtree
       ! larger than lowBound is saved, so cutting off
       ! subtrees can be done efficiently.
       ! Root is the root node of the tree. All other tree nodes
       ! can be accessed through it.
       ! call tree_init_global(iBasin, lowBound, root)
       call tree_init(nNodes,toNodes,lowBound,root,fromNodes=fromNodes,perm=permNodes)
       deallocate(toNodes,permNodes,fromNodes)

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
       !call schedule_destroy(iBasin,schedule)
       !allocate(schedule(nproc-1))
       ! call create_schedule(iBasin,nSubtrees,subtrees,schedule)
      !  call write_graphviz_output(root)
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
       call routing(iBasin,bufferLength,subtrees,nSubtrees,STmeta,permNodes,schedule,testarray)
        call write_tree_with_array(root, lowBound,testarray)

       call schedule_destroy(iBasin,schedule)
       deallocate(STmeta)

       call tree_destroy(root)
       deallocate(subtrees,permNodes,toNodes)

       call destroy_testarray(testarray)
    else
       lowBound=3
       ! all other processes receive meta data for their
       ! individual subtrees from the master process
       call get_subtree_meta(iBasin,STmeta,toNodes)
       call subtree_init(lowBound,STmeta,toNodes,subtrees)
       ! - receives data corresponding to an array and assigned
       !   subtrees
       ! - receive input data from connected subtrees
       ! - processes data
       ! - sends root data to master
       ! - send data to master
       call subtree_routing(iBasin,bufferLength,toNodes,subtrees,STmeta,testarray)
       do kk=1, size(subtrees)
          call tree_destroy(subtrees(kk))
       end do
       deallocate(subtrees)
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
  subroutine routing(iBasin,bufferLength,subtrees,nSubtrees,STmeta,permNodes,schedule,array)
    implicit none
    integer(i4),                     intent(in)     :: iBasin
    integer(i4),                     intent(in)     :: bufferLength
    type(ptrTreeNode), dimension(:), intent(in)     :: subtrees ! the array of
    integer(i4),                     intent(in)     :: nSubtrees
    type(subtreeMeta), dimension(:), intent(in)     :: STmeta
    integer(i4),       dimension(:), intent(in)     :: permNodes
    type(processSchedule), dimension(:), intent(in) :: schedule
    integer(i4),       dimension(:), intent(inout)  :: array
    ! local
    integer(i4) :: kk,iproc,next,ind, indST
    integer(i4), dimension(bufferLength+1) :: buffer
    integer(i4) :: nproc,rank,ierror
    integer status(MPI_STATUS_SIZE)

    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierror)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)
    call distribute_array(iBasin,nSubtrees,STmeta,permNodes,schedule,array)

    do kk=1,nSubtrees
       ! the process where subtree kk is assigned to
       iproc=subtrees(kk)%tN%ST%sched(1)
       call MPI_Recv(buffer,bufferLength+1,MPI_INTEGER,iproc,7,MPI_COMM_WORLD,status,ierror)
       indST=buffer(bufferLength+1)
       ! if the root node of the subtree has a parent
       if (associated(subtrees(indST)%tN%post%tN)) then
          ! next becomes the index of the parent subtree root node
          next=subtrees(indST)%tN%ST%postST%tN%ST%indST

          ! ind becomes the index of the parent of the root node
          ind=subtrees(indST)%tN%post%tN%ind
          ! iproc is the process the parent subtree is assigned to
          iproc=subtrees(next)%tN%ST%sched(1)
          ! the index, where the value will be added to
          buffer(bufferLength+1)=ind-STmeta(next)%iStart
          call MPI_Send(buffer,bufferLength+1,MPI_INTEGER,iproc,next,MPI_COMM_WORLD,ierror)
          ! write(*,*) 'master sent', buffer(1), 'to process',iproc,'tree',next, 'ind_diff', ind-STmeta(next)%iStart
       end if
    end do

    array(:)=0
    call collect_array(iBasin,nSubtrees,STmeta,permNodes,schedule,array)
  end subroutine

  subroutine subtree_init(lowBound,STmeta,toNodes,subtrees)
    implicit none
    integer(i4),                                  intent(in)    :: lowBound
    type(subtreeMeta), dimension(:),              intent(in)    :: STmeta
    integer(i4),       dimension(:),              intent(in)    :: toNodes
    type(ptrTreeNode), dimension(:), allocatable, intent(inout) :: subtrees ! the array of
    ! local
    integer(i4) :: nST ! number of subtrees scheduled on this computational node
    integer(i4) :: kk
    integer(i4) :: sizST,iStart,iEnd

    nST=size(STmeta)

    allocate(subtrees(nST))

    do kk=1,nST
      iStart=STmeta(kk)%iStart
      iEnd=STmeta(kk)%iEnd
      sizST=iEnd-iStart+1
     ! write(*,*) sizST
      call tree_init(sizST,toNodes(iStart:iEnd),lowBound,subtrees(kk))
     ! call write_subtree(subtrees(kk), lowBound)
    end do
    
  end subroutine subtree_init

  ! - receives data corresponding to an array and assigned
  !   subtrees
  ! - receive input data from connected subtrees
  ! - processes data
  ! - sends root data to master
  ! - send data to master
  subroutine subtree_routing(iBasin,bufferLength,toNodes,subtrees,STmeta,array)
    implicit none
    integer(i4),                                  intent(in)    :: iBasin
    integer(i4),                                  intent(in)    :: bufferLength
    integer(i4),       dimension(:),              intent(in)    :: toNodes
    type(ptrTreeNode), dimension(:),              intent(in)    :: subtrees ! the array of
    type(subtreeMeta), dimension(:),              intent(in)    :: STmeta
    integer(i4),       dimension(:), allocatable, intent(inout) :: array
    ! local
    integer(i4) :: nST ! number of subtrees scheduled on this computational node
    integer(i4) :: kk,jj,next,iStart,iEnd, ii
    integer(i4), dimension(:,:), allocatable :: buffer
    integer(i4) :: nproc,rank,ierror
    integer status(MPI_STATUS_SIZE)

    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierror)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)
    call get_array(iBasin,STmeta,array)

    nST=size(STmeta)
    do kk=1,nST
       allocate(buffer(STmeta(kk)%nIn+1,bufferLength+1))
       buffer(:,:)=0
       do jj=1,STmeta(kk)%nIn
          ! receive the value from the root of the child subtree and the index
          ! the value will be added to
          call MPI_Recv(buffer(jj,:),bufferLength+1,MPI_INTEGER,0,STmeta(kk)%indST,MPI_COMM_WORLD,status,ierror)
       !   write(*,*) 'process',rank, 'tree', STmeta(kk)%indST, 'with indices', &
       !           STmeta(kk)%iStart,'-',STmeta(kk)%iEnd ,&
       !           'gets', buffer(jj,1), 'for ind', next, 'ind_diff', buffer(jj,bufferLength+1)
       end do
       iStart=STmeta(kk)%iStart
       iEnd=STmeta(kk)%iEnd
       do ii=1,bufferLength
          do jj=1,STmeta(kk)%nIn
             ! shift the index with respect to the start of the subtree
             next=buffer(jj,bufferLength+1)+STmeta(kk)%iStart
             array(next)=array(next)+buffer(jj,ii)
          end do
          !call nodeinternal_routing_array(toNodes(iStart:iEnd),array(iStart:iEnd))
          !!$OMP parallel num_threads(jj) private(rank) shared(testarray)
          !$OMP parallel private(rank) shared(array,subtrees)
          !$OMP single
          call nodeinternal_routing(subtrees(kk),array(iStart:iEnd))
          !$OMP end single
          !$OMP barrier
          !$OMP end parallel
          buffer(STmeta(kk)%nIn+1,ii)=array(STmeta(kk)%iEnd)
       end do
       buffer(STmeta(kk)%nIn+1,bufferLength+1)=STmeta(kk)%indST
       call MPI_Send(buffer(STmeta(kk)%nIn+1,:),bufferLength+1,MPI_INTEGER,0,7,MPI_COMM_WORLD,ierror)
       deallocate(buffer)
    end do

    call send_array(iBasin,STmeta,array)
    deallocate(array)
  end subroutine subtree_routing

  subroutine nodeinternal_routing_array(toNodes,array)
    implicit none
    integer(i4),       dimension(:),              intent(in)    :: toNodes
    integer(i4),       dimension(:),              intent(inout) :: array
    ! local
    integer(i4) :: jj
    integer(i4) :: ii

    do jj=1,size(toNodes)-1
       array(toNodes(jj))=array(toNodes(jj))+array(jj)
    end do
  end subroutine nodeinternal_routing_array

  subroutine nodeinternal_routing(root,array)
    implicit none
    type(ptrTreeNode),                            intent(in)    :: root ! the array of
    integer(i4),       dimension(:),              intent(inout) :: array
    ! local
    integer(i4) :: jj
    integer(i4) :: ii,ithread
    integer(i4) :: tNode

    do jj=1,root%tN%Nprae
       !$OMP task shared(root,array)
       call nodeinternal_routing(root%tN%prae(jj),array)
       !$OMP end task
    end do
    !$OMP taskwait
    if (associated(root%tN%post%tN)) then
       tNode=root%tN%post%tN%origind
       !$OMP critical
      ! !$ ithread=omp_get_thread_num()
      ! !$ write(*,*) ithread
       array(tNode)=array(tNode)+array(root%tN%origind)
       !$OMP end critical
    end if
  end subroutine nodeinternal_routing

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

  subroutine get_L11_information(iBasin, toNodes, fromNodes, permNodes)
    use mo_mrm_global_variables, only : &
            level11,      & ! IN: for number of nCells
            L11_fromN,    & ! IN: for an edge this is the incoming tree node
            L11_toN,      & ! IN: for an edge this is the outgoing tree node
            L11_netPerm     ! IN: network routing order
    implicit none
    integer(i4),                            intent(in)    :: iBasin
    integer(i4), dimension(:), allocatable, intent(out)   :: toNodes
    integer(i4), dimension(:), allocatable, intent(out)   :: fromNodes
    integer(i4), dimension(:), allocatable, intent(out)   :: permNodes

   allocate(toNodes(size(L11_toN)))
   toNodes(:)=L11_toN(:)
   allocate(fromNodes(size(L11_fromN)))
   fromNodes(:)=L11_fromN(:)
   allocate(permNodes(size(L11_netPerm)))
   permNodes(:)=L11_netPerm(:)

  end subroutine get_L11_information

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

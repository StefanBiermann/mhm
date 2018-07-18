!> \file mo_common_mHM_mRM_domain_decomposition.f90

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

MODULE mo_common_mHM_mRM_domain_decomposition
  use mo_kind, only : i4, dp
  use mpi

  ! Written Maren Kaluza, June 2018

  IMPLICIT NONE

  public :: domain_decomposition

  private

  ! in Fortran one cannot have a allocatable array of
  ! pointers directly due to syntax inconsistency otherwise.
  ! The "fortran 95/2003 explained" book recommends to
  ! have a type with just one attribute, which is a pointer.
  ! Sadly this leads to an extra layer of attribute calls
  type ptrTreeNode
    type(treeNode), pointer :: tN
  end type ptrTreeNode

  ! treeNode is a tree and a node in a tree at the
  ! same time. Most important attributs concerning the
  ! tree strucuture are post and prae.
  ! post is a pointer to the parent node and prae an array
  ! of pointers to the child vertices.
  ! Nprae is the number of children, which is in the
  ! beginning the same as size(prae). During the
  ! cut of process of subtrees Nprae gets reduced and
  ! the children get switched around. When a
  ! subtree is cut of, the values of Nprae and the
  ! order of children in prae is not changed anymore.
  ! In the end, this
  ! data structure can hold the original tree and the
  ! reduced treestructure, where for each subtreetree node
  ! Nprae is set to the number of children in that subtree
  ! represented by the node
  type treeNode
    integer(i4)                                :: origind   ! index in node in original array
    integer(i4)                                :: ind       ! index in node in routing ordered array

    type(ptrTreeNode)                          :: post      ! node downstream,
                                                            ! parent
                                                            ! null, if root
    integer(i4)                                :: Nprae     ! number of children
    type(ptrTreeNode),dimension(:),allocatable :: prae      ! array of children
    integer(i4)                                :: siz       ! size of subtree
    integer(i4)                                :: sizUp     ! size of smallest
                                                            ! subtree > lowBound
                                                            ! upstream
    type(subtreeNode), pointer                 :: ST        ! data of a node, that is also a subtreetree node

    integer(i4)                                :: NSTinBranch ! number of cut of subtrees in branch

  end type treeNode

  type subtreeNode
    integer(i4)                                :: indST     ! index in node in Subtree array
    type(ptrTreeNode)                          :: postST    ! next subtree downstream,
                                                            ! parent
    integer(i4)                                :: NpraeST   ! number of cut of subtrees in node
    type(ptrTreeNode),dimension(:),allocatable :: praeST    ! array of subtree children
    integer(i4)                                :: sizST     ! size of cut of subtree
  end type subtreeNode

  type subtreeMeta
     integer(i4)        :: iStart
     integer(i4)        :: iEnd
     integer(i4)        :: indST
     integer(i4)        :: nIn
  end type subtreeMeta
CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         domain_decomposition

  !     PURPOSE
  !>        \brief decomposes the basins into subdomains for parallel computing
  !
  !>        \details decomposes basins into subdomains depending
  !>             on wether mRM is active or notse variable will then not be read in again.
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

       lowBound=3
       uppBound=5
       ! In this subroutine the tree structure gets initialized for
       ! the flownetwork of the iBasin-th basin.
       ! In each tree node the size of the smallest subtree
       ! larger than lowBound is saved, so cutting of
       ! subtrees can be done efficiently.
       ! Root is the root node of the tree. All other tree nodes
       ! can be accessed throug it.
       call init_tree(iBasin, lowBound, root)

       ! ToDo: that's possibly a bit too much, but maybe more efficient than reallocating?
       allocate(subtrees(nNodes/lowBound+1))
       ! this subroutine cuts down the tree into a subtreetree and
       ! writes each subtreetree node into an array (subtrees) in routing order
       call decompose(iBasin,lowBound,root,subtrees,nSubtrees)
       call write_domain_decomposition(root)
       allocate(STmeta(nSubtrees),permNodes(nNodes),toNodes(nNodes))
       ! A subtree data structrure makes communication between the subtrees
       ! much easier for the master. Processing the data is more efficient
       ! with array, so everything gets written into a nice array in
       ! routing order. Therefore we need a permutation array permNodes
       ! which is at the same time fromNodes, and a toNode array
       call init_subtree_metadata(iBasin,subtrees,STmeta,permNodes,toNodes)

       ! sends the meta data from master process to all the others
       call distribute_subtree_meta(iBasin,nSubtrees,STmeta,toNodes)
       ! - sends data (testarray) corresponding to subtrees to nodes
       ! - collects processed data from roots from subtrees and sends this
       !   data to corresponding leaves in connected subtrees
       ! - collects the data in the end
       call routing(iBasin,subtrees,nSubtrees,STmeta,permNodes,testarray)
    !   call write_tree_with_array(root, lowBound,testarray)

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
  subroutine routing(iBasin,subtrees,nSubtrees,STmeta,permNodes,array)
    implicit none
    integer(i4),                     intent(in)     :: iBasin
    type(ptrTreeNode), dimension(:), intent(in)     :: subtrees ! the array of
    integer(i4),                     intent(in)     :: nSubtrees
    type(subtreeMeta), dimension(:), intent(in)     :: STmeta
    integer(i4),       dimension(:), intent(in)     :: permNodes
    integer(i4),       dimension(:), intent(inout)  :: array
    ! local
    integer(i4) :: kk,iproc,next,ind
    integer(i4), dimension(2) :: value_ind
    integer(i4) :: nproc,rank,ierror
    integer status(MPI_STATUS_SIZE)

    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierror)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)
    call distribute_array(iBasin,nSubtrees,STmeta,permNodes,array)
    
    do kk=1,nSubtrees
       call iST_to_iproc(kk,nproc,iproc)
       call MPI_Recv(value_ind,2,MPI_INTEGER,iproc,7,MPI_COMM_WORLD,status,ierror)
       if (associated(subtrees(value_ind(2))%tN%post%tN)) then
          next=subtrees(value_ind(2))%tN%ST%postST%tN%ST%indST

          ind=subtrees(value_ind(2))%tN%post%tN%ind
          call iST_to_iproc(next,nproc,iproc)
          call MPI_Send([value_ind(1),ind-STmeta(next)%iStart],2,MPI_INTEGER,iproc,next,MPI_COMM_WORLD,ierror)
   !       write(*,*) 'master sent', value_ind(1), 'to process',iproc,'tree',next, 'ind_diff', ind-STmeta(next)%iStart
       end if
    end do

    array(:)=0
    call collect_array(iBasin,nSubtrees,STmeta,permNodes,array)
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
       enddo
       call nodeinternal_routing(kk,toNodes,STmeta,array)
       call MPI_Send([array(STmeta(kk)%iEnd),STmeta(kk)%indST],2,MPI_INTEGER,0,7,MPI_COMM_WORLD,ierror)
    enddo

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

  ! A subtree data structrure makes communication between the subtrees
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
  end subroutine init_subtree_metadata

  ! gets the root node of a tree, after tree decomposition.
  ! Nprae then corresponds to the cut of subtree.
  ! This subroutine starts from the root node, which is the
  ! last one in routing order.
  ! The value for ind, when the function gets called, should therefore
  ! be the size of the subtree.
  ! Then the root node gets written into the array, ind gets reduced
  ! by one, and the subroutine gets recurcively called.
  ! The last tree node in array will be the root node, and all
  ! other nodes have post nodes. So size-1 entries of the toArray
  ! can be set as the toNodes. The last entry could be usefull
  ! aswell, but is not needed.
  recursive subroutine write_tree_to_array(tree,start,ind,array,toArray)
    implicit none
    type(ptrTreeNode),         intent(in)    :: tree
    integer(i4),               intent(in)    :: start
    integer(i4),               intent(inout) :: ind
    integer(i4), dimension(:), intent(inout) :: array
    integer(i4), dimension(:), intent(inout) :: toArray
    ! local
    integer(i4) :: kk

    array(ind)=tree%tN%origind
    tree%tN%ind=start+ind-1
    if (associated(tree%tN%post%tN)) then
       toArray(ind)=tree%tN%post%tN%ind
    else
       toArray(ind)=0
    end if
    ind=ind-1
    do kk=1,tree%tN%Nprae
       call write_tree_to_array(tree%tN%prae(kk),start,ind,array,toArray)
    end do
  end subroutine write_tree_to_array

  subroutine distribute_subtree_meta(iBasin,nSubtrees,STmeta,toNodes)
    implicit none
    integer(i4),                     intent(in)  :: iBasin
    integer(i4),                     intent(in)  :: nSubtrees
    type(subtreeMeta), dimension(:), intent(in)  :: STmeta
    integer(i4),       dimension(:), intent(in)  :: toNodes
    ! local variables
    integer(i4) :: kk,ii,iPerm,iproc,sizST
    integer(i4) :: nproc,rank,ierror
    integer(i4), dimension(:,:), allocatable :: iSends ! the number and over all
                                                       ! length of arrays to be send to a process
    integer(i4), dimension(:), allocatable :: sendarray

    ! find number of processes nproc
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierror)
    ! find the number the process is referred to, called rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)
    allocate(iSends(nproc-1,2))
    iSends(:,:)=0
    ! Each process gets a subset of the subtrees in a round robin fashion.
    ! In each process an array will be allocated with a length, so that
    ! every subtree fits into it.
    do kk=1,nSubtrees
       ! for each subtree estimate the process iproc, where it gets send to
       call iST_to_iproc(kk,nproc,iproc)
       ! count the number of subtrees that will be send to process iproc
       iSends(iproc,1)=iSends(iproc,1)+1
       ! find the size of the subtree
       sizST=STmeta(kk)%iEnd+1-STmeta(kk)%iStart
       ! sum up the over all number of tree nodes send to process iproc
       iSends(iproc,2)=iSends(iproc,2)+sizST
    end do
    ! send number of subtrees and total number of tree nodes assigned
    ! to the processes to the corresponding process, so arrays can be
    ! allocated in the receiving subroutines
    do kk=1,nproc-1
       call MPI_Send(iSends(kk,:),2,MPI_INTEGER,kk,0,MPI_COMM_WORLD,ierror)
    end do
    ! send metadata of the subtrees to the noders where theres are assigned to
    ! size, identifiing index corresponding to the subtree array and number of
    ! nodes, where data is flowing into the tree
    do kk=1,nSubtrees
       call iST_to_iproc(kk,nproc,iproc)
       sizST=STmeta(kk)%iEnd+1-STmeta(kk)%iStart
       call MPI_Send(sizST,1,MPI_INTEGER,iproc,1,MPI_COMM_WORLD,ierror)
       call MPI_Send(kk,1,MPI_INTEGER,iproc,1,MPI_COMM_WORLD,ierror)
       call MPI_Send(STmeta(kk)%nIn,1,MPI_INTEGER,iproc,1,MPI_COMM_WORLD,ierror)
    end do
    ! for each subtree send corresponding toNodes to the node. Move indices from
    ! toNodes, so they start from 0
    ! They can then be moved on the receiving process corresponding to
    ! its own offset
    do kk=1,nSubtrees
       call iST_to_iproc(kk,nproc,iproc)
       sizST=STmeta(kk)%iEnd+1-STmeta(kk)%iStart
       allocate(sendarray(SizST))
       do ii=STmeta(kk)%iStart,STmeta(kk)%iEnd
          sendarray(ii-STmeta(kk)%iStart+1)=toNodes(ii)-STmeta(kk)%iStart
       end do
       call MPI_Send(sendarray(1:sizST),sizST,MPI_INTEGER,iproc,2,MPI_COMM_WORLD,ierror)
       deallocate(sendarray)
    end do
    deallocate(iSends)
  end subroutine distribute_subtree_meta

  ! calculates the rank of the process (a number, a
  ! process can be referred by) a subtree gets associated
  ! and send to by the master process
  subroutine iST_to_iproc(kk,nproc,iproc)
    implicit none
    integer(i4),                     intent(in)  :: kk
    integer(i4),                     intent(in)  :: nproc
    integer(i4),                     intent(out) :: iproc

    iproc=mod(kk-1,nproc-1)+1
  end subroutine iST_to_iproc

  subroutine distribute_array(iBasin,nSubtrees,STmeta,permNodes,array)
    implicit none
    integer(i4),                     intent(in)  :: iBasin
    integer(i4),                     intent(in)  :: nSubtrees
    type(subtreeMeta), dimension(:), intent(in)  :: STmeta
    integer(i4),       dimension(:), intent(in)  :: permNodes
    integer(i4),       dimension(:), intent(in)  :: array
    ! local variables
    integer(i4) :: kk,ii,iPerm,iproc,sizST
    integer(i4) :: nproc,rank,ierror
    integer(i4), dimension(:), allocatable :: sendarray
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierror)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

    do kk=1,nSubtrees
       call iST_to_iproc(kk,nproc,iproc)
       sizST=STmeta(kk)%iEnd+1-STmeta(kk)%iStart
       allocate(sendarray(sizST))
       do ii=STmeta(kk)%iStart,STmeta(kk)%iEnd
          iPerm=permNodes(ii)
          sendarray(ii-STmeta(kk)%iStart+1)=array(iPerm)
       end do
       call MPI_Send(sendarray(1:sizST),sizST,MPI_INTEGER,iproc,3,MPI_COMM_WORLD,ierror)
       deallocate(sendarray)
    end do
  end subroutine distribute_array

  ! collects data from the other processes into one array in
  ! original order
  subroutine collect_array(iBasin,nSubtrees,STmeta,permNodes,array)
    implicit none
    integer(i4),                     intent(in)    :: iBasin
    integer(i4),                     intent(in)    :: nSubtrees
    type(subtreeMeta), dimension(:), intent(in)    :: STmeta
    integer(i4),       dimension(:), intent(in)    :: permNodes
    integer(i4),       dimension(:), intent(inout) :: array
    ! local variables
    integer(i4) :: kk,ii,iPerm,iproc,sizST,maxSizST
    integer(i4) :: nproc,rank,ierror
    integer status(MPI_STATUS_SIZE)
    integer(i4), dimension(:), allocatable :: recvarray

    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierror)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)
    do kk=1,nSubtrees
       call iST_to_iproc(kk,nproc,iproc)
       sizST=STmeta(kk)%iEnd+1-STmeta(kk)%iStart
       allocate(recvarray(sizST))
       call MPI_Recv(recvarray(1:sizST),sizST,MPI_INTEGER,iproc,4,MPI_COMM_WORLD,status,ierror)
       do ii=STmeta(kk)%iStart,STmeta(kk)%iEnd
          iPerm=permNodes(ii)
          array(iPerm)=recvarray(ii-STmeta(kk)%iStart+1)
       end do
       deallocate(recvarray)
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

    ! receices number of subtrees and total number of tree nodes assigned to
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

  subroutine init_tree(iBasin,lowBound,root)
    use mo_mrm_global_variables, only : &
            level11,      & ! IN: for number of nCells
            L11_fromN,    & ! IN: for an edge this is the incoming tree node
            L11_toN,      & ! IN: for an edge this is the outgoing tree node
            L11_label,    & ! IN: label on edge, Id [0='', 1=HeadWater, 2=Sink]
            L11_netPerm     ! IN: network routing order
    implicit none
    integer(i4),               intent(in)    :: iBasin
    integer(i4),               intent(in)    :: lowBound
    type(ptrTreeNode),         intent(inout) :: root

    ! local variables
    integer(i4) :: kk ! loop variable to run over all edges/links
    integer(i4) :: i ! index of edge in routing order
    integer(i4) :: iNode, tNode ! in and to tree node of corresponding edge
    integer(i4) :: nLinks ! number of edges
    type(ptrTreeNode), dimension(:), allocatable :: tree
    integer(i4), dimension(:), allocatable :: Nprae

   nLinks=level11(iBasin)%ncells - 1
   allocate(tree(nLinks+1),Nprae(nLinks+1))
   do kk = 1, nLinks+1
      Nprae(kk)=0
   end do
   ! for every edge the tree node where it points to
   ! gets another child count
   do kk = 1, nLinks
      i = L11_netPerm(kk)
      tNode = L11_toN(i)
      iNode = L11_fromN(i)
      Nprae(tNode)=Nprae(tNode)+1
   end do
   ! allocate array for children, initialize tree node
   do kk = 1, nLinks+1
      allocate(tree(kk)%tN)
      tree(kk)%tN%Nprae=Nprae(kk)
      allocate(tree(kk)%tN%prae(Nprae(kk)))
      do i = 1, Nprae(kk)
         tree(kk)%tN%prae(i)%tN=>null()
      end do
      tree(kk)%tN%siz=1
      tree(kk)%tN%sizUp=1
      tree(kk)%tN%origind=kk
      tree(kk)%tN%ind=0
      tree(kk)%tN%post%tN=>null()
      tree(kk)%tN%NSTinBranch=1
      tree(kk)%tN%ST=>null()
   end do
   ! assign the pointers of the children and the parent
   ! use the array of number of children:
   ! if a child is assigned, reduce the number by one
   do kk = 1, nLinks
      i = L11_netPerm(kk)
      iNode = L11_fromN(i)
      tNode = L11_toN(i)
      ! if more children were to be assigned than Nprae
      ! then we counted wrong beforehand or we have
      ! too many children
      if (Nprae(tNode)<1) then
         write(*,*) 'something went wrong while preparing tree'
      else
         tree(tNode)%tN%prae(Nprae(tNode))%tN=>tree(iNode)%tN
         Nprae(tNode)=Nprae(tNode)-1
      end if

      ! we return only the root, and now set the parent nodes
      if (L11_label(i).eq.2) then
         root%tN => tree(tNode)%tN
      end if
      ! set the parent node
      tree(iNode)%tN%post%tN=>tree(tNode)%tN
      ! while assigning, we can calculate the size of the subtree
      tree(tNode)%tN%siz=tree(tNode)%tN%siz+tree(iNode)%tN%siz
   end do
   ! assign size of smallest subtree greater lowBound to each tree node
   do kk = 1, nLinks
      i = L11_netPerm(kk)
      tNode = L11_toN(i)
      ! call this subroutine only in routing order. Otherwise
      ! it fails
      ! this subroutine also considers sizUp of each child to
      ! be initialized with 1
      call find_sizUp_of_node(tree(tNode),lowBound)
   end do

    ! call write_tree(root,lowBound)

   deallocate(tree,Nprae)
  end subroutine init_tree

  ! can only be used in routing order, otherwise it fails
  ! this subroutine also considers sizUp of each child to
  ! be initialized with 1
  subroutine find_sizUp_of_node(node,lowBound)
    implicit none
    type(ptrTreeNode),         intent(in) :: node
    integer(i4),               intent(in) :: lowBound
    ! local variables
    integer(i4) :: kk ! loop variable to run over all tree nodes
    integer(i4) :: sizUpOfBestChild

    ! the size of the smallest subtree greater than lowBound
    sizUpOfBestChild=-1
    do kk = 1, node%tN%Nprae
       ! if one of the children already has an assigned value for
       ! the size of the smallest subtree greater than lowBound
       ! then we take it
       if (node%tN%prae(kk)%tN%sizUp .gt. 1) then
          if (sizUpOfBestChild .eq. -1) then
             sizUpOfBestChild = node%tN%prae(kk)%tN%sizUp
           ! we take the lowest of the values greater lowBound
          else if (sizUpOfBestChild .gt. node%tN%prae(kk)%tN%sizUp) then
             sizUpOfBestChild = node%tN%prae(kk)%tN%sizUp
          end if
       end if
    end do
    ! if there was a child with a value for the size of
    ! the smallest subtree greater than lowBound, we
    ! assign it to our tree node
    if (.not. sizUpOfBestChild .eq. -1) then
       node%tN%sizUp = sizUpOfBestChild
    ! else we assign its own size to that value, if it
    ! is large enough
    else if (node%tN%siz .ge. lowBound) then
       node%tN%sizUp = node%tN%siz
    else
       node%tN%sizUp = 1
    end if
    ! else it stays 1
    ! we need to set it again in case of update after tree change

  end subroutine find_sizUp_of_node

  recursive subroutine write_tree(root, lowBound)
    implicit none
    type(ptrTreeNode),         intent(in) :: root
    integer(i4),               intent(in) :: lowBound
    ! local variables
    integer(i4) :: kk ! loop variable to run over all tree nodes
    integer(i4) :: NChildren

    NChildren=size(root%tN%prae)
    write(*,*) '**********************************************************************'
    write(*,*) '* node:', root%tN%origind, 'new:',root%tN%ind, '                     *'
    write(*,*) '**********************************************************************'
    write(*,*) 'has size: ', root%tN%siz
    write(*,*) 'size of smallest subtree larger than', lowBound, 'is: ', root%tN%sizUp
    if (.not. associated(root%tN%post%tN)) then
       write(*,*) 'it is the root node'
    else
       write(*,*) 'its parent is', root%tN%post%tN%origind, 'new:', root%tN%post%tN%ind
    end if
    write(*,*) 'has', root%tN%Nprae ,'children: '
    do kk = 1, NChildren
       write(*,*) '  old:', root%tN%prae(kk)%tN%origind, 'new:', root%tN%prae(kk)%tN%ind
    end do
    do kk = 1, NChildren
       call write_tree(root%tN%prae(kk),lowBound)
    end do

  end subroutine write_tree

  recursive subroutine write_tree_with_array(root, lowBound,array)
    implicit none
    type(ptrTreeNode),         intent(in) :: root
    integer(i4),               intent(in) :: lowBound
    integer(i4), dimension(:), intent(in) :: array
    ! local variables
    integer(i4) :: kk ! loop variable to run over all tree nodes
    integer(i4) :: NChildren

    NChildren=size(root%tN%prae)
    write(*,*) '**********************************************************************'
    write(*,*) '* node:', root%tN%origind, 'new:',root%tN%ind, '                     *'
    write(*,*) '* value:', array(root%tN%origind), '                                 *'
    write(*,*) '**********************************************************************'
    write(*,*) 'has size: ', root%tN%siz
    write(*,*) 'size of smallest subtree larger than', lowBound, 'is: ', root%tN%sizUp
    if (.not. associated(root%tN%post%tN)) then
       write(*,*) 'it is the root node'
    else
       write(*,*) 'its parent is', root%tN%post%tN%origind, 'new:', root%tN%post%tN%ind
    end if
    write(*,*) 'has', root%tN%Nprae ,'children: '
    do kk = 1, NChildren
       write(*,*) '  old:', root%tN%prae(kk)%tN%origind, 'new:', root%tN%prae(kk)%tN%ind
    end do
    do kk = 1, NChildren
       call write_tree_with_array(root%tN%prae(kk),lowBound,array)
    end do

  end subroutine write_tree_with_array

  recursive subroutine write_domain_decomposition(root)
    implicit none
    type(ptrTreeNode),         intent(in) :: root
    ! local variables
    integer(i4) :: kk ! loop variable to run over all tree nodes
    integer(i4) :: NChildren
    NChildren=size(root%tN%ST%praeST)
    write(*,*) '**********************************************************************'
    write(*,*) '* node:', root%tN%origind, '                                             *'
    write(*,*) '**********************************************************************'
    write(*,*) 'has size: ', root%tN%ST%sizST
    if (.not. associated(root%tN%post%tN)) then
       write(*,*) 'it is the root node'
    else
       write(*,*) 'its parent is', root%tN%ST%postST%tN%origind
    end if
    write(*,*) 'has', root%tN%ST%NpraeST ,'children: '
    do kk = 1, NChildren
       write(*,*) '   ', root%tN%ST%praeST(kk)%tN%origind
    end do
    do kk = 1, NChildren
       call write_domain_decomposition(root%tN%ST%praeST(kk))
    end do

  end subroutine write_domain_decomposition

  ! subtrees is an array of pointers to subtrees, where we
  ! write the subtrees in routing order to.
  ! It is at least as long as the numnber we will get
  ! nSubtrees counts the subtrees.
  !
  ! From root we crawl along a fitting branch to the
  ! subtree, we want to cut of. We know which
  ! path to go because of metadata in the tree nodes
  subroutine decompose(iBasin,lowBound,root,subtrees,nSubtrees)
    implicit none
    integer(i4),                     intent(in)    :: iBasin
    integer(i4),                     intent(in)    :: lowBound
    type(ptrTreeNode),               intent(inout) :: root
    type(ptrTreeNode), dimension(:), intent(inout) :: subtrees
    integer(i4),                     intent(out)   :: nSubtrees ! number of subtrees

    ! local variables
    type(ptrTreeNode) :: subtree
    integer(i4)       :: kk

    nSubtrees=0

    ! if root has no children, the tree is probably too small
    if (root%tN%Nprae .eq. 0) then
       write(*,*) 'warning: There came a tree in with only one tree node'
       subtree%tN => root%tN
       nSubtrees=1
       subtrees(nSubtrees)%tN => root%tN
       subtree%tN%ST%indST = nSubtrees
    else
       ! set subtree to a subtree with a parent
       subtree%tN => root%tN%prae(1)%tN
       ! cut of subtrees while the root subtree
       ! is larger than one tree node
       do while (associated(subtree%tN%post%tN))
          ! cut of one subtree following a rule defined
          ! in find_branch
          call cut_of_subtree(lowBound,0,root,subtree)
          ! update the number of cut of subtrees in that branch
          ! update sizes of the smallest subtree larger
          ! than lowBound in that branch
          call update_tree(lowBound,root,subtree)
          nSubtrees=nSubtrees+1
          subtrees(nSubtrees)%tN => subtree%tN
          subtree%tN%ST%indST = nSubtrees
       end do
    end if

    ! subtrees is now the array of all subtrees and
    ! therefore the tree decomposition. But they
    ! are not connected until now. This subroutine
    ! takes care of pointers to parents and
    ! children in the subtreetree
    call init_subtreetree(nSubtrees,root,subtrees)

  end subroutine decompose

  recursive subroutine cut_of_subtree(lowBound,childInd,root,subtree)
    implicit none
    integer(i4),               intent(in)    :: lowBound
    integer(i4),               intent(in)    :: childInd
    type(ptrTreeNode),         intent(inout) :: root
    type(ptrTreeNode),         intent(inout) :: subtree

    ! local variables
    integer(i4)          :: kk,ll
    integer(i4)          :: minST,minsize,indOfST
    type(ptrTreeNode)    :: lastSibling
    logical              :: found

    ! ToDo: Check all cases
    ! There are some special cases, if we
    ! cut of the root tree:
    ! if the tree node, we handle, is the root node:
    if (.not. associated(root%tN%post%tN)) then
       found = .true.
       ! if root has no children then this is our subtree
       if (root%tN%Nprae .eq. 0) then
          subtree%tN => root%tN
          ! initialize the tree node as one of the subtreetree
          call initiate_subtreetreenode(subtree)
       ! if root has children, we have a look into
       ! the meta data of the children and decide, if
       ! we crawl along a branch
       else
          call find_branch(root,found,indOfST)
       endif
       ! if there was a branch matching the criteria, we call this
       ! function from there (and will proceed with the else, not root, case)
       if (.not. found) then
          call cut_of_subtree(lowBound,indOfST,root%tN%prae(indOfST),subtree)
       else
       ! if there was no branch matching the critera, we cut of root
          subtree%tN => root%tN
          ! initialize the tree node as one of the subtreetree
          call initiate_subtreetreenode(subtree)
       end if
    ! if root is not the tree node, we handle
    else
       ! we have a look to the branches
       call find_branch(root,found,indOfST)
       ! if there was a branch matching the criteria, we call this
       ! function from there (and will proceed with the else, not root, case)
       if (.not. found) then
          call cut_of_subtree(lowBound,indOfST,root%tN%prae(indOfST),subtree)
       ! if there was no branch matching the critera, we cut of the subtree
       else
          subtree%tN => root%tN
          ! If we cut of a subtree, all tree nodes downstream get
          ! their size reduced by the size of that subtree.
          ! Not necessary if that subtree is root.
          call update_sizes(subtree%tN%siz,subtree)
          ! initialize the tree node as one of the subtreetree
          call initiate_subtreetreenode(subtree)

          ! the parent gets one child removed
          ! it is not removed from the array
          ! it gets switched with the last child, and Nprae reduced by 1
          ! with this, we have an updated datastructure where the cut
          ! of subtree is missing
          ! but still we have the original one, because we know the original
          ! number of children was size(prae)
          lastSibling%tN => root%tN%post%tN%prae(root%tN%post%tN%Nprae)%tN
          ! root%tN%post%tN%prae(root%tN%post%tN%Nprae)%tN => root%tN%post%tN%prae(childInd)%tN
          root%tN%post%tN%prae(root%tN%post%tN%Nprae)%tN => subtree%tN
          root%tN%post%tN%prae(childInd)%tN => lastSibling%tN
          root%tN%post%tN%Nprae = root%tN%post%tN%Nprae - 1
       endif
    end if
  end subroutine cut_of_subtree

  subroutine initiate_subtreetreenode(subtree)
    implicit none
    type(ptrTreeNode),         intent(in)    :: subtree

    ! initialize the tree node as one of the subtreetree, so
    ! we can later derive this tree
    allocate(subtree%tN%ST)
    ! the size of the subtree at the moment we cut it of has
    ! exactly the size of the subtree. Each time, we cut of
    ! a subtree, we reduce the sizes of all tree nodes downstrem
    subtree%tN%ST%sizST = subtree%tN%siz
    subtree%tN%ST%NpraeST = 0
  end subroutine initiate_subtreetreenode

  subroutine find_branch(root,found,indOfST)
    implicit none
    type(ptrTreeNode),         intent(in)    :: root
    integer(i4),               intent(inout) :: indOfST
    logical,                   intent(inout) :: found

    ! local variables
    integer(i4)          :: ll,ii
    integer(i4)          :: minST,minsize
    type(ptrTreeNode)    :: lastSibling

    found = .true.
    ! find it
    !*******************************************************************************
    ! variant, where the smallest subtree is removed, regrardless of anything else *
    !*******************************************************************************
!    do kk=1,root%tN%Nprae
!       if (root%tN%prae(kk)%tN%sizUp .eq. root%tN%sizUp) then
!          found = .false.
!        !  root%tN%siz=root%tN%siz-root%tN%sizUp
!          call find_subtree(lowBound,kk,root%tN%prae(kk),subtree)
!          exit
!       end if
!    end do
    !*******************************************************************************
    ! variant, where the smallest subtree larger than lowbound is removed          *
    ! on the branch with the fewest already cut of children                        *
    !*******************************************************************************
    ! of all children find the lowest number of subtrees in its branch
    minST=root%tN%NSTinBranch
    do ll=1,root%tN%Nprae
       if (root%tN%prae(ll)%tN%NSTinBranch .lt. minST) then
          minST = root%tN%prae(ll)%tN%NSTinBranch
       end if
    end do
    ! in order of the number of subtrees in the branch of each child, find
    ! the child with the smallest subtree > k
    do ll=minST,root%tN%NSTinBranch
       minsize=1
       indOfST=1
       do ii=1,root%tN%Nprae
          if ((root%tN%prae(ii)%tN%NSTinBranch .eq. ll) .and. (root%tN%prae(ii)%tN%sizUp .gt. 1)) then
             found = .false.
             if (minsize .eq. 1) then
                minsize = root%tN%prae(ii)%tN%sizUp
                indOfST=ii
             else if (root%tN%prae(ii)%tN%sizUp .le. minsize) then
                minsize = root%tN%prae(ii)%tN%sizUp
                indOfST=ii
             end if
          end if
       end do
       if (.not. found) then
          exit
       end if
    end do
  end subroutine find_branch

  ! if a subtree is cut of, the sizes in all tree nodes
  ! downstream get reduced by the size of that subtree
  recursive subroutine update_sizes(redSize,subtree)
    implicit none
    integer(i4),               intent(in)    :: redSize
    type(ptrTreeNode),         intent(inout) :: subtree

    if (associated(subtree%tN%post%tN)) then
       subtree%tN%post%tN%siz = subtree%tN%post%tN%siz - redSize
       call update_sizes(redSize,subtree%tN%post)
    end if
  end subroutine update_sizes

  ! Quite similar to init_tree, but only
  ! for the subtreetree nodes.
  subroutine init_subtreetree(nSubtrees,root,subtrees)
    implicit none
    integer(i4),                     intent(in)    :: nSubtrees
    type(ptrTreeNode),               intent(inout) :: root
    type(ptrTreeNode), dimension(:), intent(inout) :: subtrees
    ! local variables
    integer(i4)         :: kk
    integer(i4)         :: NpraeST
    type(ptrTreeNode)   :: next


    ! find the link downstream between two tree nodes, if it is not root
    ! root lies in subtrees at position nSubtrees
    do kk=1,nSubtrees-1
       ! set next to the next tree node downstream
       next%tN => subtrees(kk)%tN%post%tN
       ! while the next tree node is not a subtreetree node
       ! set next to the next tree node downstream
       ! ToDo: check, if not root, and if, return with error?
       do while (.not. associated(next%tN%ST))
          next%tN => next%tN%post%tN
       end do
       ! ToDo: check, if next%tN is a subtreetree node?
       ! we now have next set to the parent subtreetree node
       
       ! count number of subtrees
       next%tN%ST%NpraeST=next%tN%ST%NpraeST+1

       ! set post tree node
       subtrees(kk)%tN%ST%postST%tN => next%tN
    end do
    ! allocate array of childsubtrees
    do kk=1,nSubtrees
       NpraeST=subtrees(kk)%tN%ST%NpraeST
       allocate(subtrees(kk)%tN%ST%praeST(NpraeST))
    end do
    ! again find the link downstream between two subtreetree nodes
    do kk=1,nSubtrees-1
       next%tN => subtrees(kk)%tN%post%tN
       do while ((.not. associated(next%tN%ST)) .and. (associated(next%tN%post%tN)))
          next%tN => next%tN%post%tN
       end do
       ! set the links to the children
       next%tN%ST%praeST(next%tN%ST%NpraeST)%tN => subtrees(kk)%tN
       ! sadly reduce the number of children, so we can write the
       ! children to the right place in the tree structure
       next%tN%ST%NpraeST=next%tN%ST%NpraeST-1
    end do
    ! ToDo: This is awfull... repair the subtree sizes
    ! we also could live with NpraeST being 0 now, because
    ! NpraeST = size(subtrees(kk)%tN%ST%praeST)
    do kk=1,nSubtrees
       subtrees(kk)%tN%ST%NpraeST=size(subtrees(kk)%tN%ST%praeST)
    end do
    
  end subroutine init_subtreetree

  recursive subroutine update_tree(lowBound,root,subtree)
    implicit none
    integer(i4),               intent(in)    :: lowBound
    type(ptrTreeNode),         intent(inout) :: root
    type(ptrTreeNode),         intent(inout) :: subtree

    ! local variables
    integer(i4)          :: kk
    type(ptrTreeNode)    :: parent

    ! if the parent node is not root
    if (associated(subtree%tN%post%tN)) then
       ! downstream of the cut of tree the number of cut of trees from
       ! that branch increases by 1
       subtree%tN%post%tN%NSTinBranch=subtree%tN%post%tN%NSTinBranch+1
       ! have a look to all children and update the smallest
       ! subtree in that branch larger than lowBound
       call find_sizUp_of_node(subtree%tN%post,lowBound)
       call update_tree(lowBound,root,subtree%tN%post)
    endif
  end subroutine update_tree

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

  ! better destroy virtual trees than real ones
  recursive subroutine tree_destroy(iBasin,root)
    implicit none
    integer(i4),               intent(in)    :: iBasin
    type(ptrTreeNode),         intent(inout) :: root

    ! local variables
    integer(i4)          :: kk
    
    do kk=1,size(root%tN%prae)
       call tree_destroy(iBasin,root%tN%prae(kk))
    enddo
    if (associated(root%tN%ST)) then
       deallocate(root%tN%ST%praeST)
       deallocate(root%tN%ST)
    end if
    deallocate(root%tN%prae)
    deallocate(root%tN)
  end subroutine tree_destroy

END MODULE mo_common_mHM_mRM_domain_decomposition

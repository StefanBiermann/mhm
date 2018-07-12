!> \file mo_common_mHM_mRM_domain_decomposition.f90

!> \brief decomposes the basins into subdomains for parallel computing

!> \details decomposes basins into subdomains depending
!>      on wether mRM is active or not, also hopefully in future lots of
!>      small helpfull data structure and helping subroutins

!> \authors Maren Kaluza
!> \date June 2018

MODULE mo_common_mHM_mRM_domain_decomposition
  use mo_kind, only : i4, dp
  use mpi

  ! This module sets decomposes the basins into subdomains for parallel
  ! computing

  ! Written  Maren Kaluza, June 2018

  IMPLICIT NONE

  public :: domain_decomposition

  private

  type ptrTreeNode
    type(treeNode), pointer :: tN
  end type ptrTreeNode

  type treeNode
    ! general basin information
    integer(i4)                                :: origind       ! index in node
                                                            ! array
    logical                                    :: root      ! true if the node is root
    type(ptrTreeNode)                          :: post      ! node downstream,
                                                            ! parent
    integer(i4)                                :: Nprae     ! number of children
    type(ptrTreeNode),dimension(:),allocatable :: prae      ! array of children
    integer(i4)                                :: siz       ! size of subtree
    integer(i4)                                :: sizUp     ! size of smallest
                                                            ! subtree > lowBound
                                                            ! downstream
    type(ptrTreeNode)                          :: postST    ! next subtree downstream,
                                                            ! parent
    integer(i4)                                :: NpraeST   ! number of cut of subtrees in node
    type(ptrTreeNode),dimension(:),allocatable :: praeST    ! array of subtree children
    integer(i4)                                :: sizST     ! size of cut of subtree

    integer(i4)                                :: NSTinBranch ! number of cut of subtrees in branch

  end type treeNode

  type subtreeMeta
     integer(i4)        :: iStart
     integer(i4)        :: iEnd
     integer(i4)        :: iIn
     integer(i4)        :: iOut
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
    integer(i4) :: iBasin
    integer(i4) :: lowBound,uppBound ! a subdomain should include at least
                                     ! lowBound nodes and at most uppBound
    integer(i4) :: ind               ! index of link/edge for subdomain
    type(ptrTreeNode) :: root        ! the root node of the tree structure
    type(ptrTreeNode), dimension(:), allocatable :: subtrees ! the array of
                                     ! subtrees in routing order
    integer(i4)       :: nNodes      ! the number of nodes in the original tree
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
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierror)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)
    if (rank .eq. 0) then
       write(*,*) 'the domain decomposition with mRM gets implemented now...'
       call get_number_of_basins_and_nodes(iBasin,nNodes,nBasins)
       call init_testarray(iBasin,testarray)

       lowBound=3
       uppBound=5
       call init_tree(iBasin, lowBound, root)

       ! ToDo: that's possibly a bit too much, but maybe more efficient than reallocating?
       allocate(subtrees(nNodes/lowBound+1))
       call decompose(iBasin,lowBound,root,subtrees,nSubtrees)
       ! call write_domain_decomposition(root)
       allocate(STmeta(nSubtrees),permNodes(nNodes),toNodes(nNodes))
       call init_subtree_metadata(iBasin,subtrees,STmeta,permNodes,toNodes)

       call distribute_subtrees(iBasin,nSubtrees,toNodes,STmeta,permNodes,testarray)

       deallocate(STmeta)

       call tree_destroy(iBasin,root)
       deallocate(subtrees,permNodes,toNodes)

       call destroy_testarray(testarray)
    else
       call get_subtree(iBasin)
    endif
#else
    if (rank .eq. 0) then
       write(*,*) 'the domain decomposition without mRM is not implemented yet'
    else
       write(*,*) 'need to have something to eat'
    endif
#endif
  end subroutine domain_decomposition

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
    STmeta(1)%iEnd=subtrees(1)%tN%sizST
    STmeta(1)%iIn=0
    STmeta(1)%iOut=subtrees(1)%tN%origind
    do kk=2,nSubtrees
       STmeta(kk)%iStart=Stmeta(kk-1)%iEnd+1
       STmeta(kk)%iEnd=STmeta(kk)%iStart+subtrees(kk)%tN%sizST-1
       STmeta(kk)%iIn=0
       STmeta(kk)%iOut=subtrees(kk)%tN%origind
    end do
    toNodes(:)=0
    do kk=1,nSubtrees
       ind = subtrees(kk)%tN%sizST
       call write_tree_to_array(subtrees(kk),ind,permNodes(STmeta(kk)%iStart:STmeta(kk)%iEnd), &
                                                   toNodes(STmeta(kk)%iStart:STmeta(kk)%iEnd))
    end do
    write(*,*) toNodes
  end subroutine init_subtree_metadata

  recursive subroutine write_tree_to_array(tree,ind,array,toArray)
    implicit none
    type(ptrTreeNode),         intent(in)    :: tree
    integer(i4),               intent(inout) :: ind
    integer(i4), dimension(:), intent(inout) :: array
    integer(i4), dimension(:), intent(inout) :: toArray
    ! local
    integer(i4) :: kk

    array(ind)=tree%tN%origind
    if (.not. tree%tN%root) then
       toArray(ind)=tree%tN%post%tN%origind
    else
       toArray(ind)=0
    end if
    ind=ind-1
    do kk=1,tree%tN%Nprae
       call write_tree_to_array(tree%tN%prae(kk),ind,array,toArray)
    end do
  end subroutine write_tree_to_array

  subroutine distribute_subtrees(iBasin,nSubtrees,toNodes,STmeta,permNodes,testarray)
    implicit none
    integer(i4),                     intent(in)     :: iBasin
    integer(i4),                     intent(in)     :: nSubtrees
    integer(i4),       dimension(:), intent(in)     :: toNodes
    type(subtreeMeta), dimension(:), intent(inout)  :: STmeta
    integer(i4),       dimension(:), intent(inout)  :: permNodes
    integer(i4),       dimension(:), intent(inout)  :: testarray
    ! local variables
    integer(i4) :: kk,ii,iPerm,iproc,sizST,maxSizST
    integer(i4) :: nproc,rank,ierror
    integer(i4), dimension(:,:), allocatable :: iSends ! the number and over all
                                                       ! length of arrays to be send to a process
    integer(i4), dimension(:), allocatable :: sendarray

    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierror)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)
    allocate(iSends(nproc-1,2))
    write(*,*) 'I am root'
    iSends(:,:)=0
    maxSizST=0
    do kk=1,nSubtrees
       iproc=mod(kk-1,nproc-1)+1 
       iSends(iproc,1)=iSends(iproc,1)+1
       sizST=STmeta(kk)%iEnd+1-STmeta(kk)%iStart
       iSends(iproc,2)=iSends(iproc,2)+sizST
       if (maxSizST .lt. sizST) maxSizST=sizST
    end do
    allocate(sendarray(maxSizST))
    sendarray(:)=0
    do kk=1,nproc-1
       write(*,*) 'I send',iSends(kk,1),'to',kk
       call MPI_Send(iSends(kk,:),2,MPI_INTEGER,kk,0,MPI_COMM_WORLD,ierror)
    end do
    do kk=1,nSubtrees
       iproc=mod(kk-1,nproc-1)+1 
       sizST=STmeta(kk)%iEnd+1-STmeta(kk)%iStart
       call MPI_Send(sizST,1,MPI_INTEGER,iproc,1,MPI_COMM_WORLD,ierror)
    end do
    do kk=1,nSubtrees
       iproc=mod(kk-1,nproc-1)+1 
       sizST=STmeta(kk)%iEnd+1-STmeta(kk)%iStart
       do ii=STmeta(kk)%iStart,STmeta(kk)%iEnd
          sendarray(ii-STmeta(kk)%iStart+1)=toNodes(ii)
       end do
       call MPI_Send(sendarray(1:sizST),sizST,MPI_INTEGER,iproc,2,MPI_COMM_WORLD,ierror)

       do ii=STmeta(kk)%iStart,STmeta(kk)%iEnd
          iPerm=permNodes(ii)
          sendarray(ii-STmeta(kk)%iStart+1)=testarray(iPerm)
       end do
       call MPI_Send(sendarray(1:sizST),sizST,MPI_INTEGER,iproc,3,MPI_COMM_WORLD,ierror)
    end do
    write(*,*) 'I sent some messages to everyone else'
    deallocate(sendarray)
    deallocate(iSends)
  end subroutine distribute_subtrees

  subroutine get_subtree(iBasin)
    implicit none
    integer(i4),               intent(in)                  :: iBasin
    ! local variables
    integer(i4) :: kk
    integer(i4), dimension(2) :: nDatasets ! number of incoming data sets
    integer(i4), dimension(:), allocatable :: toNodes
    integer(i4), dimension(:), allocatable :: subtreeArray
    type(subtreeMeta), dimension(:), allocatable :: STmeta
    integer(i4) :: mes,sizST
    integer(i4) :: nproc,rank,ierror
    integer status(MPI_STATUS_SIZE)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierror)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

    mes=2
    call MPI_Recv(nDatasets(:),2,MPI_INTEGER,0,0,MPI_COMM_WORLD,status,ierror)
    allocate(STmeta(nDatasets(1)))
    allocate(subtreeArray(nDatasets(2)),toNodes(nDatasets(2)))
    ! ToDo: case: less subtrees than processes
    call MPI_Recv(sizST,1,MPI_INTEGER,0,1,MPI_COMM_WORLD,status,ierror)
    write(*,*) 'process', rank, 'gets', nDatasets(1), 'data sets, with', nDatasets(2), 'length'
    STmeta(1)%iStart=1
    STmeta(1)%iEnd=sizST
    STmeta(1)%iIn=0
    STmeta(1)%iOut=0
    do kk=2,nDatasets(1)
       call MPI_Recv(sizST,1,MPI_INTEGER,0,1,MPI_COMM_WORLD,status,ierror)
       STmeta(kk)%iStart=Stmeta(kk-1)%iEnd+1
       STmeta(kk)%iEnd=STmeta(kk)%iStart+sizST-1
       STmeta(kk)%iIn=0
       STmeta(kk)%iOut=0
    end do

    do kk=1,nDatasets(1)
       sizST=STmeta(kk)%iEnd+1-STmeta(kk)%iStart
       call MPI_Recv(toNodes(STmeta(kk)%iStart:STmeta(kk)%iEnd),sizST,MPI_INTEGER,0,2,MPI_COMM_WORLD,status,ierror)
       write(*,*) 'process', rank, 'ate', toNodes(STmeta(kk)%iStart:STmeta(kk)%iEnd)
       call MPI_Recv(subtreeArray(STmeta(kk)%iStart:STmeta(kk)%iEnd),sizST,MPI_INTEGER,0,3,MPI_COMM_WORLD,status,ierror)
       write(*,*) 'process', rank, 'ate', subtreeArray(STmeta(kk)%iStart:STmeta(kk)%iEnd)
    end do
    deallocate(subtreeArray,toNodes)
    deallocate(STmeta)

  end subroutine get_subtree

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
      testarray(kk)=kk
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
            L11_fromN,    & ! IN: from node
            L11_toN,      & ! IN: to node
            L11_label,    & ! IN: label on edge, Id [0='', 1=HeadWater, 2=Sink]
            L11_rOrder,   & ! IN: network routing order
            L11_netPerm     ! IN: network routing order
    implicit none
    integer(i4),               intent(in)    :: iBasin
    integer(i4),               intent(in)    :: lowBound
    type(ptrTreeNode),         intent(inout) :: root

    ! local variables
    integer(i4) :: kk ! loop variable to run over all edges/links
    integer(i4) :: i ! index of edge in routing order
    integer(i4) :: iNode, tNode ! in and to node of corresponding edge
    integer(i4) :: nLinks ! number of edges
    type(ptrTreeNode), dimension(:), allocatable :: tree
    integer(i4), dimension(:), allocatable :: Nprae

   nLinks=level11(iBasin)%ncells - 1
   allocate(tree(nLinks+1),Nprae(nLinks+1))
   do kk = 1, nLinks+1
      Nprae(kk)=0
   end do
   ! for every edge the node where it points to
   ! gets another child count
   do kk = 1, nLinks
      i = L11_netPerm(kk)
      tNode = L11_toN(i)
      iNode = L11_fromN(i)
      Nprae(tNode)=Nprae(tNode)+1
   end do
   ! allocate array for children
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
      tree(kk)%tN%root=.false.
      tree(kk)%tN%NpraeST=-1
      tree(kk)%tN%NSTinBranch=1
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
         tree(tNode)%tN%root=.true.
         nullify(tree(tNode)%tN%post%tN)
      end if
      ! set the parent node
      tree(iNode)%tN%post%tN=>tree(tNode)%tN
      ! while assigning, we can calculate the size of the subtree
      tree(tNode)%tN%siz=tree(tNode)%tN%siz+tree(iNode)%tN%siz
   end do
   ! assign size of smallest subtree greater lowBound to each node
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
    type(ptrTreeNode),         intent(in) :: node
    integer(i4),               intent(in) :: lowBound
    ! local variables
    integer(i4) :: kk ! loop variable to run over all nodes
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
    ! assign it to our node
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
    type(ptrTreeNode),         intent(in) :: root
    integer(i4),               intent(in) :: lowBound
    ! local variables
    integer(i4) :: kk ! loop variable to run over all nodes
    integer(i4) :: NChildren

    NChildren=size(root%tN%prae)
    write(*,*) '**********************************************************************'
    write(*,*) '* node:', root%tN%origind, '                                             *'
    write(*,*) '**********************************************************************'
    write(*,*) 'has size: ', root%tN%siz
    write(*,*) 'size of smallest subtree larger than', lowBound, 'is: ', root%tN%sizUp
    if (root%tN%root) then
       write(*,*) 'it is the root node'
    else
       write(*,*) 'its parent is', root%tN%post%tN%origind
    end if
    write(*,*) 'has', root%tN%Nprae ,'children: '
    do kk = 1, NChildren
       write(*,*) '   ', root%tN%prae(kk)%tN%origind
    end do
    do kk = 1, NChildren
       call write_tree(root%tN%prae(kk),lowBound)
    end do

  end subroutine write_tree

  recursive subroutine write_domain_decomposition(root)
    type(ptrTreeNode),         intent(in) :: root
    ! local variables
    integer(i4) :: kk ! loop variable to run over all nodes
    integer(i4) :: NChildren

    NChildren=size(root%tN%praeST)
    write(*,*) '**********************************************************************'
    write(*,*) '* node:', root%tN%origind, '                                             *'
    write(*,*) '**********************************************************************'
    write(*,*) 'has size: ', root%tN%sizST
    if (root%tN%root) then
       write(*,*) 'it is the root node'
    else
       write(*,*) 'its parent is', root%tN%postST%tN%origind
    end if
    write(*,*) 'has', root%tN%NpraeST ,'children: '
    do kk = 1, NChildren
       write(*,*) '   ', root%tN%praeST(kk)%tN%origind
    end do
    do kk = 1, NChildren
       call write_domain_decomposition(root%tN%praeST(kk))
    end do

  end subroutine write_domain_decomposition

  subroutine decompose(iBasin,lowBound,root,subtrees,nSubtrees)
    use mo_mrm_global_variables, only : &
            level11        ! IN: for number of nCells
    integer(i4),                     intent(in)    :: iBasin
    integer(i4),                     intent(in)    :: lowBound
    type(ptrTreeNode),               intent(inout) :: root
    type(ptrTreeNode), dimension(:), intent(inout) :: subtrees
    integer(i4),                     intent(out)   :: nSubtrees ! number of subtrees

    ! local variables
    type(ptrTreeNode) :: subtree
    integer(i4)       :: kk

    nSubtrees=0

    do while (root%tN%sizUp .gt. 1)
       call cut_of_subtree(lowBound,0,root,subtree)
       call update_tree(lowBound,root,subtree)
       !call write_tree(root, lowBound)
       nSubtrees=nSubtrees+1
       subtrees(nSubtrees)%tN => subtree%tN
    end do

    ! cut of root
    call cut_of_subtree(lowBound,0,root,subtree)
    call update_tree(lowBound,root,subtree)
    nSubtrees=nSubtrees+1
    subtrees(nSubtrees)%tN => subtree%tN
    
    call init_subtreetree(nSubtrees,root,subtrees)

  end subroutine decompose

  recursive subroutine cut_of_subtree(lowBound,childInd,root,subtree)
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
    if (root%tN%root) then
       root%tN%NpraeST = 0
       found = .true.
       ! if root has no children then this is our subtree
       if (root%tN%Nprae .eq. 0) then
          subtree%tN => root%tN
          root%tN%sizST = root%tN%siz
          ! initialize the node as one of the subtreetree, so
          ! we can later derive this tree
          subtree%tN%NpraeST = 0
       else
          call find_branch(root,found,indOfST)
       endif
       if (.not. found) then
          call cut_of_subtree(lowBound,indOfST,root%tN%prae(indOfST),subtree)
       else
          subtree%tN => root%tN
          subtree%tN%sizST = subtree%tN%siz
          ! initialize the node as one of the subtreetree, so
          ! we can later derive this tree
          subtree%tN%NpraeST = 0
       end if
    else
       call find_branch(root,found,indOfST)
       if (.not. found) then
          call cut_of_subtree(lowBound,indOfST,root%tN%prae(indOfST),subtree)
       else
          ! if found, cut it of
          subtree%tN => root%tN
          subtree%tN%sizST = subtree%tN%sizUp
          call update_sizes(subtree%tN%siz,subtree)
          ! initialize the node as one of the subtreetree, so
          ! we can later derive this tree
          subtree%tN%NpraeST = 0
          ! the parent gets one child removed
          ! it is not removed from the array
          ! it gets switched with the last child, and Nprae reduced by 1
          lastSibling%tN => root%tN%post%tN%prae(root%tN%post%tN%Nprae)%tN
          root%tN%post%tN%prae(root%tN%post%tN%Nprae)%tN => root%tN%post%tN%prae(childInd)%tN
          root%tN%post%tN%prae(childInd)%tN => lastSibling%tN
          root%tN%post%tN%Nprae = root%tN%post%tN%Nprae - 1
       endif
    end if
  end subroutine cut_of_subtree

  subroutine find_branch(root,found,indOfST)
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
    ! variant, where the smallest subtree is removed, on the branch with the fewest*
    ! already cut of children                                                      *
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

  recursive subroutine update_sizes(redSize,subtree)
    integer(i4),               intent(in)    :: redSize
    type(ptrTreeNode),         intent(inout) :: subtree

    if (.not. subtree%tN%root) then
       subtree%tN%post%tN%siz = subtree%tN%post%tN%siz - redSize
       call update_sizes(redSize,subtree%tN%post)
    end if
  end subroutine

  ! ToDo: awful subroutine, to be exchanged later
  subroutine init_subtreetree(nSubtrees,root,subtrees)
    integer(i4),                     intent(in)    :: nSubtrees
    type(ptrTreeNode),               intent(inout) :: root
    type(ptrTreeNode), dimension(:), intent(inout) :: subtrees
    ! local variables
    integer(i4)                            :: kk
    integer(i4), dimension(:), allocatable :: NpraeST
    type(ptrTreeNode)                      :: next

    allocate(NpraeST(nSubtrees))

    ! find the link downstream between two nodes, if it is not root
    do kk=1,nSubtrees-1
       next%tN => subtrees(kk)%tN%post%tN
       do while ((next%tN%NpraeST .eq. -1) .and. (.not. next%tN%root))
          next%tN => next%tN%post%tN
       end do
       ! count number of subtrees
       next%tN%NpraeST=next%tN%NpraeST+1
       ! set post node
       subtrees(kk)%tN%postST%tN => next%tN
    end do
    ! allocate array of childsubtrees
    do kk=1,nSubtrees
       NpraeST(kk)=subtrees(kk)%tN%NpraeST
       allocate(subtrees(kk)%tN%praeST(NpraeST(kk)))
    end do
 !   allocate(root%tN%praeST(root%tN%NpraeST))
    ! again find the link downstream between two nodes
    do kk=1,nSubtrees-1
       next%tN => subtrees(kk)%tN%post%tN
       do while ((next%tN%NpraeST .eq. -1) .and. (.not. next%tN%root))
          next%tN => next%tN%post%tN
       end do
       ! set the links to the children
       next%tN%praeST(next%tN%NpraeST)%tN => subtrees(kk)%tN
       next%tN%NpraeST=next%tN%NpraeST-1
    end do
    ! ToDo: This is awfull... repair the subtree sizes
    do kk=1,nSubtrees
       subtrees(kk)%tN%NpraeST=NpraeSt(kk)
    end do

    !ToDo: why is that?
    next%tN => subtrees(1)%tN

    deallocate(NpraeST)
    
  end subroutine init_subtreetree

  recursive subroutine update_tree(lowBound,root,subtree)
    integer(i4),               intent(in)    :: lowBound
    type(ptrTreeNode),         intent(inout) :: root
    type(ptrTreeNode),         intent(inout) :: subtree

    ! local variables
    integer(i4)          :: kk
    type(ptrTreeNode)    :: parent

    if (.not. subtree%tN%root) then
       subtree%tN%post%tN%NSTinBranch=subtree%tN%post%tN%NSTinBranch+1
       call find_sizUp_of_node(subtree%tN%post,lowBound)
       call update_tree(lowBound,root,subtree%tN%post)
    endif
  end subroutine update_tree

  subroutine find_L11_fAcc(iBasin,L11_fAcc)
    use mo_mrm_global_variables, only : &
            level11,      & ! IN: for number of nCells
            L11_fromN,    & ! IN: from node
            L11_toN,      & ! IN: to node
            L11_label,    & ! IN: label on edge, Id [0='', 1=HeadWater, 2=Sink]
            L11_rOrder,   & ! IN: network routing order
            L11_netPerm     ! IN: network routing order
    implicit none
    integer(i4),               intent(in)    :: iBasin
    integer(i4), dimension(:), intent(inout) :: L11_fAcc

    ! local variables
    integer(i4) :: kk ! loop variable to run over all edges/links
    integer(i4) :: i ! index of edge in routing order
    integer(i4) :: iNode, tNode ! in and to node of corresponding edge
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
!   do kk=1, nLinks+1
!      write(*,*) kk, L11_fAcc(kk)
!   end do
  end subroutine find_L11_fAcc

  recursive subroutine tree_destroy(iBasin,root)
    integer(i4),               intent(in)    :: iBasin
    type(ptrTreeNode),         intent(inout) :: root

    ! local variables
    integer(i4)          :: kk
    
    do kk=1,size(root%tN%prae)
       call tree_destroy(iBasin,root%tN%prae(kk))
    enddo
    deallocate(root%tN%prae)
    deallocate(root%tN)
  end subroutine tree_destroy

END MODULE mo_common_mHM_mRM_domain_decomposition

!> \file mo_HRD_tree_init_and_destroy.f90

!> \authors Maren Kaluza
!> \date June 2018

MODULE mo_HRD_tree_init_and_destroy
  use mo_kind, only : i4, dp
  use mo_HRD_types, only : ptrTreeNode
  use mo_HRD_tree_tools, only : find_sizUp_of_node,find_revLevel_of_node
  !$ use omp_lib,      only: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
  use mpi

  ! Written Maren Kaluza, June 2018

  IMPLICIT NONE

  public :: tree_init, tree_init_global, tree_destroy
            

  private

CONTAINS

  subroutine tree_init_global(iBasin,lowBound,root)
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
      tree(kk)%tN%level=0
      tree(kk)%tN%revLevel=0
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
   ! set level
   do kk = nLinks,1,-1
      i = L11_netPerm(kk)
      iNode = L11_fromN(i)
      tree(iNode)%tN%level=tree(iNode)%tN%post%tN%level+1
   end do
   ! assign size of smallest subtree greater lowBound to each tree node
   ! assign farthest distant to a leave
   do kk = 1, nLinks
      i = L11_netPerm(kk)
      tNode = L11_toN(i)
      ! call this subroutine only in routing order. Otherwise
      ! it fails
      ! this subroutine also considers sizUp of each child to
      ! be initialized with 1
      call find_sizUp_of_node(tree(tNode),lowBound)
      call find_revLevel_of_node(tree(tNode))
   end do

    ! call write_tree(root,lowBound)

   deallocate(tree,Nprae)
  end subroutine tree_init_global

  subroutine tree_init(iBasin,nLinks,toNodes,fromNodes,perm,lowBound,root)
    implicit none
    integer(i4),                         intent(in)    :: iBasin
    integer(i4),                         intent(in)    :: nLinks
    integer(i4), dimension(:),           intent(in)    :: toNodes
    integer(i4), dimension(:), optional, intent(in)    :: fromNodes
    integer(i4), dimension(:), optional, intent(in)    :: perm
    integer(i4),                         intent(in)    :: lowBound
    type(ptrTreeNode),                   intent(inout) :: root

    ! local variables
    integer(i4) :: kk ! loop variable to run over all edges/links
    integer(i4) :: i ! index of edge in routing order
    integer(i4) :: iNode, tNode ! in and to tree node of corresponding edge
    type(ptrTreeNode), dimension(:), allocatable :: tree
    integer(i4), dimension(:), allocatable :: Nprae

    allocate(tree(nLinks),Nprae(nLinks))
    do kk = 1, nLinks
       Nprae(kk)=0
    end do
    ! for every edge the tree node where it points to
    ! gets another child count
    call count_children(nLinks-1,toNodes,perm,Nprae)
    ! allocate array for children, initialize tree node
    do kk = 1, nLinks
       call treenode_init(Nprae(kk),kk,tree(kk))
    end do
    ! assign the pointers of the children and the parent
    ! use the array of number of children:
    ! if a child is assigned, reduce the number by one
    call set_edges(nLinks-1,toNodes,fromNodes,perm,tree,Nprae,root)
    ! set level
    call set_level(nLinks-1,fromNodes,perm,tree)
    ! assign size of smallest subtree greater lowBound to each tree node
    ! assign farthest distant to a leave
    do kk = 1, nLinks-1
       if (present(perm)) then
          i = perm(kk)
       else
          i = kk
       end if
       tNode = toNodes(i)
       ! call this subroutine only in routing order. Otherwise
       ! it fails
       ! this subroutine also considers sizUp of each child to
       ! be initialized with 1
       call find_sizUp_of_node(tree(tNode),lowBound)
       call find_revLevel_of_node(tree(tNode))
    end do
    ! call write_tree(root,lowBound)
    deallocate(tree,Nprae)
  end subroutine tree_init

  subroutine treenode_init(Nprae,kk,treenode)
    implicit none
    integer(i4),       intent(in)    :: Nprae
    integer(i4),       intent(in)    :: kk ! original index of the node
    type(ptrTreeNode), intent(inout) :: treenode
    !local
    integer(i4) :: i

    allocate(treenode%tN)
    treenode%tN%Nprae=Nprae
    allocate(treenode%tN%prae(Nprae))
    do i = 1, Nprae
       treenode%tN%prae(i)%tN=>null()
    end do
    treenode%tN%siz=1
    treenode%tN%sizUp=1
    treenode%tN%origind=kk
    treenode%tN%ind=0
    treenode%tN%post%tN=>null()
    treenode%tN%NSTinBranch=1
    treenode%tN%level=0
    treenode%tN%revLevel=0
    treenode%tN%ST=>null()
  end subroutine treenode_init

  subroutine set_level(nLinks,fromNodes,perm,tree)
    implicit none
    integer(i4),                         intent(in)    :: nLinks
    integer(i4), dimension(:), optional, intent(in)    :: fromNodes
    integer(i4), dimension(:), optional, intent(in)    :: perm
    type(ptrTreeNode), dimension(:),     intent(in)    :: tree
    ! local
    integer(i4) :: kk, i, iNode

    do kk = nLinks,1,-1
       if (present(perm)) then
          i = perm(kk)
       else
          i = kk
       end if
       if (present(fromNodes)) then
          iNode = fromNodes(i)
       else
          iNode = i
       end if
       tree(iNode)%tN%level=tree(iNode)%tN%post%tN%level+1
    end do
  end subroutine set_level

  subroutine set_edges(nLinks,toNodes,fromNodes,perm,tree,Nprae,root)
    implicit none
    integer(i4),                         intent(in)    :: nLinks
    integer(i4), dimension(:),           intent(in)    :: toNodes
    integer(i4), dimension(:), optional, intent(in)    :: fromNodes
    integer(i4), dimension(:), optional, intent(in)    :: perm
    type(ptrTreeNode), dimension(:),     intent(in)    :: tree
    integer(i4), dimension(:),           intent(inout) :: Nprae
    type(ptrTreeNode),                   intent(inout) :: root
    ! local
    integer(i4) :: kk, i, iNode, tNode
    do kk = 1, nLinks
       if (present(perm)) then
          i = perm(kk)
       else
          i = kk
       end if
       if (present(fromNodes)) then
          iNode = fromNodes(i)
       else
          iNode = i
       end if
       tNode=toNodes(i)
       ! if more children were to be assigned than Nprae
       ! then we counted wrong beforehand or we have
       ! too many children
       if (Nprae(tNode)<1) then
          write(*,*) 'something went wrong while preparing tree'
       else
          tree(tNode)%tN%prae(Nprae(tNode))%tN=>tree(iNode)%tN
          Nprae(tNode)=Nprae(tNode)-1
       end if

    !   ! we return only the root, and now set the parent nodes
    !   if (L11_label(i).eq.2) then
    !      root%tN => tree(tNode)%tN
    !   end if
       ! set the parent node
       tree(iNode)%tN%post%tN=>tree(tNode)%tN
       ! while assigning, we can calculate the size of the subtree
       tree(tNode)%tN%siz=tree(tNode)%tN%siz+tree(iNode)%tN%siz
    end do
    ! find root node
    do kk = 1, nLinks
       if (.not. associated(tree(kk)%tN%post%tN)) then
          root%tN => tree(kk)%tN
       end if
    end do
  end subroutine set_edges

  subroutine count_children(nLinks,toNodes,perm,Nprae)
    implicit none
    integer(i4),                         intent(in)    :: nLinks
    integer(i4), dimension(:),           intent(in)    :: toNodes
    integer(i4), dimension(:), optional, intent(in)    :: perm
    integer(i4), dimension(:),           intent(inout) :: Nprae
    ! local
    integer(i4) :: i,iNode,tNode, kk

    do kk = 1, nLinks
      if (present(perm)) then
         i = perm(kk)
      else
         i = kk
      end if
      tNode=toNodes(i)
      Nprae(tNode)=Nprae(tNode)+1
    end do
  end subroutine count_children

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


END MODULE mo_HRD_tree_init_and_destroy

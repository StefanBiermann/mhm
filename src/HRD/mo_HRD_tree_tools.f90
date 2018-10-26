!> \file mo_HRD_tree_tools.f90

!> \authors Maren Kaluza
!> \date June 2018

MODULE mo_HRD_tree_tools
  use mo_kind, only : i4, dp
  use mo_HRD_types, only: ptrTreeNode
  !$ use omp_lib,      only: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS

  ! Written Maren Kaluza, June 2018

  IMPLICIT NONE

  public :: find_sizUp_of_node, find_revLevel_of_node, update_sizes, write_tree_to_array
            

  private

CONTAINS

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

  subroutine find_revLevel_of_node(node)
    implicit none
    type(ptrTreeNode),         intent(in) :: node
    ! local variables
    integer(i4) :: kk ! loop variable to run over all tree nodes
    integer(i4) :: maxDist

    maxDist=-1
    ! find maxDist of all children
    do kk = 1, node%tN%Nprae
       if (maxDist .lt. node%tN%prae(kk)%tN%revLevel) then
          maxDist=node%tN%prae(kk)%tN%revLevel
       end if
    end do
    node%tN%revLevel=maxDist+1
  end subroutine find_revLevel_of_node

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

  ! gets the root node of a tree, after tree decomposition.
  ! Nprae then corresponds to the cut of subtree.
  ! This subroutine starts from the root node, which is the
  ! last one in routing order.
  ! The value for ind, when the function gets called first, should therefore
  ! be set to the size of the subtree.
  ! Then the root node gets written into the array, ind gets reduced
  ! by one, and the subroutine gets recursively called.
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

END MODULE mo_HRD_tree_tools

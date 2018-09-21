!> \file mo_HRD_decompose.f90

!> \authors Maren Kaluza
!> \date June 2018

MODULE mo_HRD_decompose
  use mo_kind, only : i4, dp
  use mo_HRD_types, only: ptrTreeNode
  use mo_HRD_tree_tools, only : update_sizes,find_sizUp_of_node
  !$ use omp_lib,      only: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
  use mpi

  ! Written Maren Kaluza, June 2018

  IMPLICIT NONE

  public :: decompose
            

  private

CONTAINS

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
    integer(i4)          :: indOfST
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
          ! ToDo: moved update_tree from decompose here. Check if this made
          ! new errors
          call update_sizes(subtree%tN%siz,subtree)
          ! update the number of cut of subtrees in that branch
          ! update sizes of the smallest subtree larger
          ! than lowBound in that branch
       !   call update_tree(lowBound,root,subtree)
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
    subtree%tN%ST%postST%tN => null()
    subtree%tN%ST%NpraeST = 0
    subtree%tN%ST%levelST(1) = 0
    subtree%tN%ST%levelST(2) = 0
    !ToDo: remove this later perhabs, debugging information
    subtree%tN%ST%sched(1) = 0
    subtree%tN%ST%sched(2) = 0
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
    !*******************************************************************************
    ! variant, where the smallest subtree is removed, regrardless of anything else *
    !*******************************************************************************
    ! of all children find the lowest number of subtrees in its branch
    ! find the child with the smallest subtree > k
    minsize=1
    indOfST=1
    do ii=1,root%tN%Nprae
       if (root%tN%prae(ii)%tN%sizUp .gt. 1) then
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
  end subroutine find_branch

  subroutine find_branch_least_cut_first(root,found,indOfST)
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
    sizloop : do ll=minST,root%tN%NSTinBranch
       minsize=1
       indOfST=1
       do ii=1,root%tN%Nprae
          if ((root%tN%prae(ii)%tN%NSTinBranch .eq. ll) &
                  .and. (root%tN%prae(ii)%tN%sizUp .gt. 1)) then
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
          exit sizloop
       end if
    end do sizloop
  end subroutine find_branch_least_cut_first


  ! Quite similar to init_tree, but only
  ! for the subtreetree nodes.
  subroutine init_subtreetree(nSubtrees,root,subtrees)
    implicit none
    integer(i4),                     intent(in)    :: nSubtrees
    type(ptrTreeNode),               intent(inout) :: root
    type(ptrTreeNode), dimension(:), intent(inout) :: subtrees
    ! local variables
    integer(i4)         :: kk,jj
    integer(i4)         :: NpraeST,maxDist
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
    ! set farthest distance to a leave
    do kk=1,nSubtrees
       if (associated(subtrees(kk)%tN%ST%postST%tN)) then
          maxDist=0
          next%tN=>subtrees(kk)%tN%ST%postST%tN
          do jj=1,next%tN%ST%NpraeST
             if (maxDist .lt. next%tN%ST%praeST(jj)%tN%ST%levelST(2)) then
                maxDist = next%tN%ST%praeST(jj)%tN%ST%levelST(2)
             end if
          end do
          next%tN%ST%levelST(2)=maxDist+1
       end if
    end do
    ! set level of nodes
    do kk=nSubtrees-1,1,-1
       subtrees(kk)%tN%ST%levelST(1)=subtrees(kk)%tN%ST%postST%tN%ST%levelST(1)+1
    end do
    
  end subroutine init_subtreetree

  recursive subroutine update_tree(lowBound,root,subtree)
    implicit none
    integer(i4),               intent(in)    :: lowBound
    type(ptrTreeNode),         intent(inout) :: root
    type(ptrTreeNode),         intent(inout) :: subtree

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

END MODULE mo_HRD_decompose

!> \file mo_HRD_write.f90

!> \authors Maren Kaluza
!> \date June 2018

MODULE mo_HRD_write
  use mo_kind, only : i4, dp
  use mo_HRD_types
  !$ use omp_lib,      only: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS

  ! Written Maren Kaluza, June 2018

  IMPLICIT NONE

  public :: write_tree, write_subtree, write_tree_with_array, &
            write_domain_decomposition, write_graphviz_output

  private

CONTAINS

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

  recursive subroutine write_subtree(root, lowBound)
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
       call write_subtree(root%tN%prae(kk),lowBound)
    end do

  end subroutine write_subtree

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
    write(*,*) 'distance to root: ',root%tN%level
    write(*,*) 'distance to farthest leave: ',root%tN%revLevel
    write(*,*) 'size of smallest subtree larger than', lowBound, 'is: ', root%tN%sizUp
    if (.not. associated(root%tN%post%tN)) then
       write(*,*) 'it is the root node'
    else
       write(*,*) 'its parent is', root%tN%post%tN%origind, 'new:', root%tN%post%tN%ind
    end if
    if (associated(root%tN%ST)) then
       write(*,*) 'it is a tree node in the subtreetree with index:', root%tN%ST%indST
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
    write(*,*) '* node:', root%tN%origind,'new ind:', root%tN%ind, 'subtreetree ind:', root%tN%ST%indST, '*'
    write(*,*) '**********************************************************************'
    write(*,*) 'farthest distant leave:', root%tN%ST%levelST(2)
    write(*,*) 'level:', root%tN%ST%levelST(1)
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

  recursive subroutine write_graphviz_output(root)
    implicit none
    type(ptrTreeNode),         intent(in) :: root
    write(*,*) 'graph ""'
    write(*,*) '{'
    call write_graphviz_output_nodes(root)
    write(*,*) '}'

  end subroutine write_graphviz_output

  recursive subroutine write_graphviz_output_nodes(root)
    implicit none
    type(ptrTreeNode),         intent(in) :: root
    ! local variables
    integer(i4) :: kk ! loop variable to run over all tree nodes
    integer(i4) :: NChildren
    NChildren=size(root%tN%ST%praeST)
    if (.not. associated(root%tN%post%tN)) then
       write(*,*) root%tN%origind,';'
    end if
    write(*,*) root%tN%origind, '[label="',root%tN%ST%indST
    write(*,*) 'proc:',root%tN%ST%sched(1),'slot:',root%tN%ST%sched(2)
    write(*,*) 'siz:', root%tN%ST%sizST
    write(*,*) '"]', ';'
    do kk = 1, NChildren
       write(*,*) root%tN%origind,'--', root%tN%ST%praeST(kk)%tN%origind
    end do
    do kk = 1, NChildren
       call write_graphviz_output_nodes(root%tN%ST%praeST(kk))
    end do
  end subroutine write_graphviz_output_nodes

END MODULE mo_HRD_write

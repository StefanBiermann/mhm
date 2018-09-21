!> \file mo_HRD_subtree_meta_init_and_destroy.f90

!> \authors Maren Kaluza
!> \date June 2018

MODULE mo_HRD_subtree_meta_init_and_destroy
  use mo_kind, only : i4, dp
  use mo_HRD_types, only : ptrTreeNode, subtreeMeta
  use mo_HRD_tree_tools, only : write_tree_to_array
  !$ use omp_lib,      only: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
  use mpi

  ! Written Maren Kaluza, June 2018

  IMPLICIT NONE

  public :: init_subtree_metadata
            

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

END MODULE mo_HRD_subtree_meta_init_and_destroy

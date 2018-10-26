!> \file mo_HRD_dLinkedList.f90

!> \authors Maren Kaluza
!> \date September 2018

MODULE mo_HRD_dLinkedList
  use mo_kind, only : i4, dp
  use mo_HRD_types, only: ptrTreeNode
  !$ use omp_lib,      only: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS

  ! Written Maren Kaluza, June 2018

  IMPLICIT NONE

  public :: dLinkedList, list_init, list_destroy, list_insert, list_del,&
                                      list_write, list_append, list_prepend
            

  private

  type dLinkedList
     type(ptrTreeNode)            :: content
     type(dLinkedList), pointer   :: next
     type(dLinkedList), pointer   :: prev
  end type dLinkedList

CONTAINS
! Initialize a head node self
  subroutine list_init(tN, self)
    implicit none
    type(ptrTreeNode),          intent(in)    :: tN
    type(dLinkedList), pointer, intent(inout) :: self

    allocate(self)
    nullify(self%next)
    nullify(self%prev)

    self%content = tN
  end subroutine list_init

  ! destroy list at self
  subroutine list_destroy(self)
    implicit none
    type(dLinkedList), pointer, intent(inout) :: self
    ! local
    type(dLinkedList), pointer :: current
    type(dLinkedList), pointer :: next

    current => self
    do while (associated(current))
       next => current%next
       deallocate(current)
       nullify(current)
       current => next
    end do
  end subroutine list_destroy

  ! Insert a list node after self
  subroutine list_insert(tN, self)
    implicit none
    type(ptrTreeNode),          intent(in)    :: tN
    type(dLinkedList), pointer, intent(inout) :: self
    ! local
    type(dLinkedList), pointer :: next

    allocate(next)

    next%content = tN

    next%next => self%next
    self%next => next

    next%prev => self
    if (associated(next%next)) then
       next%next%prev => next
    end if
  end subroutine list_insert

  subroutine list_del(self)
    implicit none
    type(dLinkedList), pointer, intent(inout) :: self

    if (associated(self%next)) then
       self%next%prev=>self%prev
    end if
    if (associated(self%prev)) then
       self%prev%next=>self%next
    end if
    deallocate(self)
  end subroutine list_del

  subroutine list_append(self,list)
    implicit none
    type(dLinkedList), pointer, intent(inout) :: self
    type(dLinkedList), pointer, intent(inout) :: list

    if (.not. associated(self)) then
       self=>list
    else if (associated(list)) then
       self%next=>list
       list%prev=>self
    end if
  end subroutine list_append

  subroutine list_prepend(list,self)
    implicit none
    type(dLinkedList), pointer, intent(inout) :: self
    type(dLinkedList), pointer, intent(inout) :: list

    if (.not. associated(self)) then
       self=>list
    else if (associated(list)) then
       self%prev=>list
       list%next=>self
    end if
  end subroutine list_prepend

  subroutine list_write(self)
    implicit none
    type(dLinkedList), pointer, intent(in) :: self
    !local
    type(dLinkedList), pointer :: element

    element=>self
    
    write(*,*) '*******************'
    if (associated(element)) then
       write(*,*) element%content%tN%ST%indST,element%content%tN%ST%levelST(1)
       do while(associated(element%next))
          element=>element%next
          write(*,*) element%content%tN%ST%indST,element%content%tN%ST%levelST(1)
       end do
    else
       write(*,*) 'the list has no elements'
    end if
  end subroutine list_write

END MODULE mo_HRD_dLinkedList

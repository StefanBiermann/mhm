!> \file mo_HRD_schedule.f90

!> \authors Maren Kaluza
!> \date September 2018

MODULE mo_HRD_schedule
  use mo_kind, only : i4, dp
  use mo_HRD_types, only: ptrTreeNode, processSchedule
  use mo_HRD_tree_tools, only : update_sizes,find_sizUp_of_node
  !$ use omp_lib,      only: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
  use mpi

  ! Written Maren Kaluza, June 2018

  IMPLICIT NONE

  public :: create_schedule
            

  private

CONTAINS

  subroutine create_schedule(iBasin,nSubtrees,subtrees,schedule)
    implicit none
    integer(i4),                     intent(in)  :: iBasin
    integer(i4),                     intent(in)  :: nSubtrees
    type(ptrTreeNode), dimension(:), intent(inout)     :: subtrees ! the array of
    type(processSchedule), dimension(:), allocatable, intent(inout) :: schedule
    ! local variables
    integer(i4) :: kk,iproc,sizST,place
    integer(i4) :: nproc,ierror
    integer(i4), dimension(:), allocatable :: sendarray

    ! find number of processes nproc
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierror)

    ! sort the array of subtrees, so that distant leaves come first
    call sort_subtrees(nSubtrees,subtrees)

    ! initialize
    do kk=1,nproc
       schedule(kk)%nTrees=0
       schedule(kk)%overallSize=0
    end do
    ! Each process gets a subset of the subtrees in a round robin fashion.
    ! In each process an array will be allocated with a length, so that
    ! every subtree fits into it.
    do kk=1,nSubtrees
       ! for each subtree estimate the process iproc, where it gets send to
       call iST_to_iproc(kk,nproc,iproc)
       ! count the number of subtrees that will be send to process iproc
       schedule(iproc)%nTrees=schedule(iproc)%nTrees+1
       ! find the size of the subtree
       sizST=subtrees(kk)%tN%ST%sizST
       ! sum up the over all number of tree nodes send to process iproc
       schedule(iproc)%overallSize=schedule(iproc)%overallSize+sizST
    end do
    ! allocate array of subtree indices
    do kk=1,nproc-1
       allocate(schedule(kk)%trees(schedule(kk)%nTrees))
       schedule(kk)%trees(:)=0
    end do
    ! write subtree indices to index array
    do kk=1,nSubtrees
       ! for each subtree estimate the process iproc, where it gets send to
       call iST_to_iproc(kk,nproc,iproc)
       place=size(schedule(iproc)%trees)-schedule(iproc)%nTrees+1
       schedule(iproc)%trees(place)=kk
       schedule(iproc)%nTrees=schedule(iproc)%nTrees-1
    end do
    ! repair sizes
    do kk=1,nproc-1
       schedule(kk)%Ntrees=size(schedule(kk)%trees)
    end do
    call schedule_to_subtree_nodes(iBasin,schedule,nSubtrees,subtrees)
  end subroutine create_schedule

 !ToDo: remove parts of this, if necessary, debugging information
  subroutine schedule_to_subtree_nodes(iBasin,schedule,nSubtrees,subtrees)
    implicit none
    integer(i4),                     intent(in)  :: iBasin
    type(processSchedule), dimension(:), allocatable, intent(in) :: schedule
    integer(i4),                     intent(in)  :: nSubtrees
    type(ptrTreeNode), dimension(:), intent(inout)     :: subtrees ! the array of
    ! local variables
    integer(i4) :: kk,jj,iTree,iproc
    integer(i4) :: nproc,ierror
    type(ptrTreeNode) :: childsubtree

    do kk=1,size(schedule)
       do jj=1,schedule(kk)%nTrees
          iTree=schedule(kk)%trees(jj)
          subtrees(iTree)%tN%ST%sched(1)=kk
          subtrees(iTree)%tN%ST%sched(2)=jj
       end do
    end do
    do kk=1,nSubtrees
       do jj=1,subtrees(kk)%tN%ST%NpraeST
          childsubtree%tN=>subtrees(kk)%tN%ST%praeST(jj)%tN
          if (childsubtree%tN%ST%sched(2).ge.subtrees(kk)%tN%ST%sched(2)) then
             subtrees(kk)%tN%ST%sched(2)=childsubtree%tN%ST%sched(2)+1
             exit
          end if
       end do
    end do
  end subroutine schedule_to_subtree_nodes

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

  ! sort the tree nodes first by their farthest distant leave and
  ! if this is equal, then by their distance to root
  subroutine sort_subtrees(nSubtrees,subtrees)
    implicit none
    integer(i4),                     intent(in)    :: nSubtrees
    type(ptrTreeNode), dimension(:), intent(inout) :: subtrees
    ! local variables
    integer(i4) :: kk,maxRevLevel
    integer(i4), dimension(:), allocatable :: posRevLevel
    integer(i4), dimension(:), allocatable :: kLevel

    call sort_by_component(2,subtrees(1:nSubtrees),posRevLevel)
    ! sort equal revRevel parts by level
    maxRevlevel=size(posRevLevel)-2
    do kk=1,maxRevlevel+1
       call sort_by_component(1,subtrees(posRevLevel(kk):posRevLevel(kk+1)-1),kLevel)
       deallocate(kLevel)
    end do
    ! update indices
    do kk=1,nSubtrees
       subtrees(kk)%tN%ST%indST=kk
    end do
    deallocate(posRevLevel)
  end subroutine sort_subtrees

  ! mix the sorted array so that successive arrays are distant
  ! cache line problems with openMP
  subroutine unsort_subtrees(nSubtrees,subtrees)
    implicit none
    integer(i4),                     intent(in)    :: nSubtrees
    type(ptrTreeNode), dimension(:), intent(inout) :: subtrees
    ! local variables
    type(ptrTreeNode), dimension(:), allocatable :: newSubtrees
    real(dp) :: goldenRatio
    integer(i4) :: step
    integer(i4) :: ii

    ! set golden ratio
    goldenRatio=(1.0+sqrt(5.0))/2.0
    ! find a nice stepsize
    step=int(real(nSubtrees)/goldenRatio)
    if (mod(nSubtrees,step) .eq. 0) then
       step=step+1
    end if
    ! ToDo: catch error
    if (step .ge. nSubtrees) then
      write(*,*) 'module decomposition, unsort_subtrees, strange things happend'
    end if
    allocate(newSubtrees(nSubtrees))
    do ii=1,nSubtrees
       newSubtrees(ii)%tN=>subtrees(mod(step*ii,nSubtrees)+1)%tN
    end do
    do ii=1,nSubtrees
      subtrees(ii)%tN=>newSubtrees(ii)%tN
      subtrees(ii)%tN%ST%indST=ii
    end do
    deallocate(newSubtrees)

  end subroutine unsort_subtrees

  subroutine sort_by_component(ind,subtrees,posLevel)
    implicit none
    integer(i4),                            intent(in)    :: ind
    type(ptrTreeNode), dimension(:),        intent(inout) :: subtrees
    integer(i4), dimension(:), allocatable, intent(inout) :: posLevel
    ! local variables
    integer(i4) :: nSubtrees
    integer(i4) :: kk,pos
    integer(i4) :: maxLevel
    integer(i4), dimension(:), allocatable :: kLevel
    type(ptrTreeNode), dimension(:), allocatable :: newSubtrees
    ! find the maximum level
    maxLevel=0
    nSubtrees=size(subtrees)
    do kk=1,nSubtrees
       if (maxLevel .lt. subtrees(kk)%tN%ST%levelST(ind)) then
          maxLevel = subtrees(kk)%tN%ST%levelST(ind)
       end if
    end do
    allocate(kLevel(maxLevel+2),posLevel(maxLevel+2),newSubtrees(nSubtrees))
    do kk=1,nSubtrees
       newSubtrees(kk)%tN=>null()
    end do
    kLevel(:)=0
    posLevel(:)=0
    ! lowest level is 0. We want to count the occurency of this level in component 2
    do kk=1,nSubtrees
       kLevel(subtrees(kk)%tN%ST%levelST(ind)+2) &
               =kLevel(subtrees(kk)%tN%ST%levelST(ind)+2)+1
    end do
    kLevel(1)=1
    ! find position, where the parts with equal level in array start and end
    do kk=2,maxLevel+1
       kLevel(kk)=kLevel(kk)+kLevel(kk-1)
    end do
    ! the last component should be the number of subtrees +1, so
    ! kLevel(kk)-kLevel(kk-1) gives us the number of occurencies of level kk-1
    kLevel(maxLevel+2)=nSubtrees+1
    ! save this positions in posLevel to give it back
    posLevel(:)=kLevel(:)
    ! do the actual sorting
    do kk=1,nSubtrees
       pos=kLevel(subtrees(kk)%tN%ST%levelST(ind)+1)
       kLevel(subtrees(kk)%tN%ST%levelST(ind)+1) &
               =kLevel(subtrees(kk)%tN%ST%levelST(ind)+1)+1
       newSubtrees(pos)%tN=>subtrees(kk)%tN
    end do
    ! write to array
    ! ToDo: rather reverse tag?
    do kk=1,nSubtrees
       if (ind .eq. 2) then
          subtrees(kk)%tN=>newSubtrees(kk)%tN
       else
          subtrees(nSubtrees-kk+1)%tN=>newSubtrees(kk)%tN
       end if
    end do
    deallocate(kLevel,newSubtrees)

  end subroutine sort_by_component

END MODULE mo_HRD_schedule

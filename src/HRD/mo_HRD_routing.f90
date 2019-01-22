!> \file mo_HRD_routing.f90

!> \authors Maren Kaluza
!> \date January 2019

MODULE mo_HRD_routing
  use mo_kind, only : i4, dp

  use mo_HRD_types, only: subtreeMeta, MPI_parameter, ptrTreeNode, processSchedule
!
!  use mo_HRD_tree_init_and_destroy, only: tree_init_global, tree_init, tree_destroy, &
!                                          forest_init, forest_destroy, tree_init_buffer
!
  use mo_HRD_tree_tools, only: tree_init_qTIN_with_array, tree_init_qTR_with_array, &
                   tree_init_C1_with_array, tree_init_C2_with_array, tree_init_qOut_with_array, &
                   tree_extract_qTIN_in_array, tree_extract_qTR_in_array
!
!  use mo_HRD_decompose, only: decompose
!
!  use mo_HRD_schedule, only: create_schedule, create_schedule_hu, schedule_destroy
!
!  use mo_HRD_subtree_meta, only: init_subtree_metadata, distribute_subtree_meta, &
!                                 get_subtree_meta, destroy_subtree_meta
  use mo_HRD_MPI_array_communication, only: get_array_dp, send_array_dp, get_full_array_dp, send_full_array_dp
!
  !$ use omp_lib,      only: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
  use mpi_f08

  ! Written Maren Kaluza, January 2019

  IMPLICIT NONE

  public :: muskignum_subtree_routing_process, muskignum_master_routing

  private

CONTAINS
  
  subroutine muskignum_master_routing(iBasin, MPIparam, subtrees, nSubtrees, STmeta, schedule)
    implicit none
    integer(i4),                         intent(in)    :: iBasin
    type(MPI_parameter),                 intent(in)    :: MPIparam
    type(ptrTreeNode),     dimension(:), intent(in)    :: subtrees
    integer(i4),                         intent(in)    :: nSubtrees
    type(subtreeMeta),     dimension(:), intent(inout) :: STmeta
    type(processSchedule), dimension(:), intent(in)    :: schedule
    ! local
    integer(i4)                            :: kk, iproc, next, indST
    integer(i4), dimension(MPIparam%bufferLength+1) :: buffer
    real(dp), dimension(MPIparam%bufferLength+1) :: qTIN, qTR
    integer(i4)                            :: ierror
    type(MPI_Status)                       :: status

    do kk = 1, nSubtrees
      ! for each subtree the master process 0 gets the data of
      ! the root tree node
      call MPI_Recv(buffer, MPIparam%bufferLength+1, MPI_INTEGER, MPI_ANY_SOURCE, 0, MPIparam%comm, status, ierror)
      indST = buffer(1)
      call MPI_Recv(qTIN,   MPIparam%bufferLength+1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, indST, MPIparam%comm, status, ierror)
      call MPI_Recv(qTR,    MPIparam%bufferLength+1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, indST, MPIparam%comm, status, ierror)
      ! if the root node of the subtree has a parent
      if (associated(subtrees(indST)%tN%post%tN)) then
        ! next becomes the index of the parent subtree root node (the subtree
        ! that gets data)
        next = subtrees(indST)%tN%ST%postST%tN%ST%indST
        ! iproc is the process the parent subtree is assigned to
        iproc = subtrees(next)%tN%ST%sched(1)
        call MPI_Send(buffer, MPIparam%bufferLength+1, MPI_INTEGER, iproc, indST, MPIparam%comm, ierror)
        call MPI_Send(qTIN,   MPIparam%bufferLength+1, MPI_DOUBLE_PRECISION, iproc, indST, MPIparam%comm, ierror)
        call MPI_Send(qTR,    MPIparam%bufferLength+1, MPI_DOUBLE_PRECISION, iproc, indST, MPIparam%comm, ierror)
      end if
    end do
  end subroutine muskignum_master_routing

  subroutine muskignum_subtree_routing_process(MPIParam, routLoop, STmeta, subtrees, inInds, trees, inTrees)
    type(MPI_parameter),                          intent(in)    :: MPIparam
    integer(i4),                                  intent(in)    :: routLoop
    type(subtreeMeta), dimension(:),              intent(inout) :: STmeta
    type(ptrTreeNode), dimension(:), allocatable, intent(in)    :: subtrees            ! the array of
    integer(i4),       dimension(:), allocatable, intent(in)    :: inInds
    type(ptrTreeNode), dimension(:), allocatable, intent(inout) :: trees               ! the array of all nodes of the trees
    type(ptrTreeNode), dimension(:), allocatable, intent(in)    :: inTrees             ! the array of halo leaves
    ! local
    integer(i4) :: ierror
    integer(i4)                            :: nInflowGauges
    type(MPI_Status)                       :: status
    integer(i4) :: iBasin
    
    real(dp), dimension(:, :), allocatable :: L11_C1
    real(dp), dimension(:, :), allocatable :: L11_C2
    real(dp), dimension(:, :), allocatable :: L11_qOut
    real(dp), dimension(:),    allocatable :: L11_qTIN
    real(dp), dimension(:),    allocatable :: L11_qTR
    real(dp), dimension(:, :), allocatable :: L11_buf_qTIN
    real(dp), dimension(:, :), allocatable :: L11_buf_qTR
    logical,  dimension(:),    allocatable :: InflowGaugeHeadwater
    integer(i4), dimension(:), allocatable :: InflowGaugeNodeList

    integer(i4) :: tt
    ! ToDo: change
    iBasin = 1
    do tt = 1, routLoop
      call MPIparam%increment()
      if (MPIparam%buffered) then
        call MPI_Recv(nInflowGauges, 1, MPI_INTEGER, 0, 2, MPIparam%comm, status, ierror)
        allocate(InflowGaugeHeadwater(nInflowGauges), InflowGaugeNodeList(nInflowGauges))
        call MPI_Recv(InflowGaugeHeadwater, nInflowGauges, MPI_LOGICAL, 0, 2, MPIparam%comm, status, ierror)
        call MPI_Recv(InflowGaugeNodeList, nInflowGauges, MPI_INTEGER, 0, 2, MPIparam%comm, status, ierror)
        call get_full_array_dp(iBasin, MPIParam, STmeta, L11_C1)
        call tree_init_C1_with_array(L11_C1, trees)
        call get_full_array_dp(iBasin, MPIParam, STmeta, L11_C2)
        call tree_init_C2_with_array(L11_C2, trees)
        call get_full_array_dp(iBasin, MPIParam, STmeta, L11_qOut)
        call tree_init_qOut_with_array(L11_qOut, trees)
        call get_array_dp(iBasin, MPIParam%nproc, MPIParam%rank, MPIParam%comm, STmeta, L11_qTIN)
        call tree_init_qTIN_with_array(L11_qTIN, trees)
        call get_array_dp(iBasin, MPIParam%nproc, MPIParam%rank, MPIParam%comm, STmeta, L11_qTR)
        call tree_init_qTR_with_array(L11_qTR, trees)
        call muskignum_subprocess_routing(MPIParam, InflowGaugeHeadwater, InflowGaugeNodeList,&
                                       STmeta, subtrees, inInds, inTrees, trees)
        call MPI_Barrier(MPIparam%comm)
        allocate(L11_buf_qTIN(size(L11_qTIN, dim=1), MPIparam%bufferLength))
        allocate(L11_buf_qTR( size(L11_qTIN, dim=1), MPIparam%bufferLength))
        call tree_extract_qTIN_in_array(trees, L11_buf_qTIN)
        call send_full_array_dp(iBasin, MPIparam, STmeta, L11_buf_qTIN)
        call tree_extract_qTR_in_array(trees, L11_buf_qTR)
        call send_full_array_dp(iBasin, MPIparam, STmeta, L11_buf_qTR)
        call MPI_Barrier(MPIparam%comm)
        deallocate(L11_buf_qTIN, L11_buf_qTR)
        deallocate(InflowGaugeHeadwater, InflowGaugeNodeList)
      end if
    end do
    
  end subroutine muskignum_subtree_routing_process

  subroutine muskignum_subprocess_routing(MPIParam, InflowGaugeHeadwater, InflowGaugeNodeList, &
                                       STmeta, subtrees, inInds, inTrees, trees)
    type(MPI_parameter), intent(in) :: MPIparam
    logical,           dimension(:), allocatable, intent(in)    :: InflowGaugeHeadwater
    integer(i4),       dimension(:), allocatable, intent(in)    :: InflowGaugeNodeList
    type(subtreeMeta), dimension(:),              intent(inout) :: STmeta
    type(ptrTreeNode), dimension(:), allocatable, intent(in)    :: subtrees            ! the array of
    integer(i4),       dimension(:), allocatable, intent(in)    :: inInds
    type(ptrTreeNode), dimension(:), allocatable, intent(in)    :: inTrees             ! the array of halo leaves
    type(ptrTreeNode), dimension(:), allocatable, intent(inout) :: trees               ! the array of all nodes of the trees
    ! local
    integer(i4) :: ierror
    integer(i4)                            :: routLoop
    integer(i4)                            :: nInflowGauges
    type(MPI_Status)                       :: status
    integer(i4) :: iBasin

    integer(i4) :: nST, kk, jj, ii, nIn

    nST = size(STmeta)

    do kk = 1, nST
      do jj = 1, STmeta(kk)%nIn
        ii = jj + STmeta(kk)%iInStart - 1
        ! for every subtree for every leaf, where data comes in, receive this
        ! data. Via inInds(ii) we write the right incoming data to the
        ! matching buffer.
        call MPI_IRecv(inTrees(ii)%tN%values%buffer, MPIparam%bufferLength+1, &
                                         MPI_INTEGER, 0, inInds(ii), MPIParam%comm, inTrees(ii)%tN%values%request, ierror)
        call MPI_IRecv(inTrees(ii)%tN%qTIN%buffer,   MPIparam%bufferLength+1, &
                                         MPI_DOUBLE_PRECISION, 0, inInds(ii), MPIParam%comm, inTrees(ii)%tN%qTIN%request,   ierror)
        call MPI_IRecv(inTrees(ii)%tN%qTR%buffer,    MPIparam%bufferLength+1, &
                                         MPI_DOUBLE_PRECISION, 0, inInds(ii), MPIParam%comm, inTrees(ii)%tN%qTR%request,    ierror)
      end do
    end do
    ! now we run over the subtrees on the computing node in routing order
    do kk = 1, nST
      nIn = STmeta(kk)%nIn
      do jj = 1, nIn
        ii = jj + STmeta(kk)%iInStart - 1
        ! we can start calculations in that subtree, when all its inflows are
        ! reveived
        call MPI_Wait(inTrees(ii)%tN%values%request, inTrees(ii)%tN%values%status, ierror)
        call MPI_Wait(inTrees(ii)%tN%qTIN%request,   inTrees(ii)%tN%qTIN%status,   ierror)
        call MPI_Wait(inTrees(ii)%tN%qTR%request,    inTrees(ii)%tN%qTR%status,    ierror)
      !  if (.not. MPI_ASYNC_PROTECTS_NONBLOCKING) call MPI_F_SYNC_REG(inTrees(ii)%tN%values%buffer)
      !  if (.not. MPI_ASYNC_PROTECTS_NONBLOCKING) call MPI_F_SYNC_REG(inTrees(ii)%tN%qTIN%buffer)
      !  if (.not. MPI_ASYNC_PROTECTS_NONBLOCKING) call MPI_F_SYNC_REG(inTrees(ii)%tN%qTR%buffer)
      end do
      !!$OMP parallel num_threads(jj) private(rank)
      !$OMP parallel private(rank) shared(subtrees)
      !$OMP single
      call muskignum_subtree_routing(MPIParam, InflowGaugeHeadwater, InflowGaugeNodeList, &
                                       STmeta, subtrees(kk), inInds, inTrees)
      !$OMP end single
      !$OMP barrier
      !$OMP end parallel
      subtrees(kk)%tN%values%buffer(1) = STmeta(kk)%indST
      ! send the outflow to the master process
      call MPI_Send(subtrees(kk)%tN%values%buffer, MPIparam%bufferLength+1, MPI_INTEGER, 0, 0, MPIparam%comm, ierror)
      call MPI_Send(subtrees(kk)%tN%qTIN%buffer,   MPIparam%bufferLength+1, MPI_DOUBLE_PRECISION, 0, STmeta(kk)%indST, MPIparam%comm, ierror)
      call MPI_Send(subtrees(kk)%tN%qTR%buffer,    MPIparam%bufferLength+1, MPI_DOUBLE_PRECISION, 0, STmeta(kk)%indST, MPIparam%comm, ierror)
    end do
  end subroutine muskignum_subprocess_routing

  recursive subroutine muskignum_subtree_routing(MPIParam, InflowGaugeHeadwater, InflowGaugeNodeList, &
                                       STmeta, root, inInds, inTrees)
    type(MPI_parameter), intent(in) :: MPIparam
    logical,           dimension(:), allocatable, intent(in)    :: InflowGaugeHeadwater
    integer(i4),       dimension(:), allocatable, intent(in)    :: InflowGaugeNodeList
    type(subtreeMeta), dimension(:),              intent(in)    :: STmeta
    type(ptrTreeNode),                            intent(in)    :: root             ! the array of
    integer(i4),       dimension(:), allocatable, intent(in)    :: inInds
    type(ptrTreeNode), dimension(:), allocatable, intent(in)    :: inTrees             ! the array of halo leaves
    ! local
    integer(i4) :: kk

    do kk = 1, root%tN%ST%NpraeST
      !$OMP task shared(root)
      call muskignum_subtree_routing(MPIParam, InflowGaugeHeadwater, InflowGaugeNodeList, &
                                       STmeta, root%tN%ST%praeST(kk), inInds, inTrees)
      !$OMP end task
    end do
    !$OMP taskwait
    call muskignum_subtree_routing_serial(MPIParam, InflowGaugeHeadwater, InflowGaugeNodeList, &
                                       STmeta, root, inInds, inTrees)
  end subroutine muskignum_subtree_routing

  recursive subroutine muskignum_subtree_routing_serial(MPIParam, InflowGaugeHeadwater, InflowGaugeNodeList, &
                                       STmeta, root, inInds, inTrees)
    type(MPI_parameter), intent(in) :: MPIparam
    logical,           dimension(:), allocatable, intent(in)    :: InflowGaugeHeadwater
    integer(i4),       dimension(:), allocatable, intent(in)    :: InflowGaugeNodeList
    type(subtreeMeta), dimension(:),              intent(in)    :: STmeta
    type(ptrTreeNode),                            intent(in)    :: root             ! the array of
    integer(i4),       dimension(:), allocatable, intent(in)    :: inInds
    type(ptrTreeNode), dimension(:), allocatable, intent(in)    :: inTrees             ! the array of halo leaves
    ! local
    integer(i4) :: ierror
    integer(i4)                            :: routLoop
    integer(i4)                            :: nInflowGauges
    type(MPI_Status)                       :: status
    integer(i4) :: iBasin

    integer(i4) :: kk, jj, bufferLength

    do kk = 1, root%tN%Nprae
      call muskignum_subtree_routing_serial(MPIParam, InflowGaugeHeadwater, InflowGaugeNodeList, &
                                     STmeta, root%tN%prae(kk), inInds, inTrees)
    end do
    bufferLength = MPIparam%bufferLength
    ! add routed water to downstream node
   !   netNode_qTIN(tNode, IT) = netNode_qTIN(tNode, IT) + netNode_qTR(iNode, IT)
    do jj = 2, bufferLength+1
      if (size(root%tN%prae) > 0 .and. .not. root%tN%isIn) then
       ! write(*,*) root%tN%qTIN%buffer(jj), jj
        do kk = 1, size(root%tN%prae)
         ! write(*,*) root%tN%prae(kk)%tN%qTR%buffer(jj-1), &
         ! root%tN%prae(kk)%tN%qTR%buffer(jj), root%tN%permind, root%tN%prae(kk)%tN%permind, root%tN%prae(kk)%tN%isIn
          root%tN%qTIN%buffer(jj) = root%tN%qTIN%buffer(jj) + &
                                      root%tN%prae(kk)%tN%qTR%buffer(jj)
        end do
      end if

    !  if (size(root%tN%prae) == 0 .and. .not. root%tN%isIn) then
    !    root%tN%qTR%buffer(jj) = root%tN%qTR%buffer(jj-1)
    !    root%tN%qTIN%buffer(jj) = root%tN%qTIN%buffer(jj-1)
    !  else if (size(root%tN%prae) > 0) then
    !    do kk = 1, size(root%tN%prae)
    !      root%tN%qTIN%buffer(jj) = root%tN%qTIN%buffer(jj) + &
    !                                  root%tN%prae(kk)%tN%qTR%buffer(jj)
    !    end do
    !  end if
      if (.not. root%tN%isIn) then
      ! accumulate all inputs in iNode
     ! netNode_qTIN(iNode, IT) = netNode_qTIN(iNode, IT) + netNode_qOUT(iNode)
      root%tN%qTIN%buffer(jj) = root%tN%qTIN%buffer(jj) + root%tN%qOut%buffer(jj-1)
             !  write(*,*) root%tN%qTIN%buffer(jj), root%tN%qOut%buffer(jj-1), root%tN%permInd

      ! routing iNode
      ! Here is a difference compared to the old code. This is also done for the
      ! root node, but wasn't before. But qTR is a temporary result and this last
      ! value has no impact anyhow. So we leave it as it is.
     !   netNode_qTR(iNode, IT) = netNode_qTR(iNode, IT1)                               &
     !           + netLink_C1(i) * (netNode_qTIN(iNode, IT1) - netNode_qTR (iNode, IT1)) &
     !           + netLink_C2(i) * (netNode_qTIN(iNode, IT) - netNode_qTIN(iNode, IT1))
      root%tN%qTR%buffer(jj) = root%tN%qTR%buffer(jj-1) &
                 + root%tN%C1%buffer(jj-1) * ( root%tN%qTIN%buffer(jj-1) -  root%tN%qTR%buffer(jj-1) ) &
                 + root%tN%C2%buffer(jj-1) * ( root%tN%qTIN%buffer(jj)   - root%tN%qTIN%buffer(jj-1) )
               !  write(*,*) root%tN%C1%buffer(jj-1),root%tN%qTIN%buffer(jj-1), root%tN%qTR%buffer(jj-1),&
               !             root%tN%C2%buffer(jj-1),root%tN%qTIN%buffer(jj), &
               !             root%tN%permind
               !  write(*,*) '*****************'
               !   write(*,*) root%tN%qTR%buffer(jj), root%tN%permind

      !ToDo: is only necessery in incomplete catchments
      ! check if the inflow from upstream cells should be deactivated
     !   if (nInflowGauges .GT. 0) then
     !     do i = 1, nInflowGauges
     !       ! check if downstream Node (tNode) is inflow gauge and headwaters should be ignored
     !       if ((tNode == InflowNodeList(i)) .AND. (.NOT. InflowHeadwater(i))) netNode_qTR(iNode, IT) = 0.0_dp
     !     end do
     !   end if

      end if
   end do
  end subroutine muskignum_subtree_routing_serial

END MODULE mo_HRD_routing

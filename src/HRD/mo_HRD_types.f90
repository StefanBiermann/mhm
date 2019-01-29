!> \file mo_HRD_types.f90

!> \authors Maren Kaluza
!> \date June 2018

MODULE mo_HRD_types
  use mo_kind, only : i4, dp
  use mpi_f08

  ! Written Maren Kaluza, June 2018

  IMPLICIT NONE

  public :: ptrTreeNode, treeNode, subtreeMeta, subtreeBuffer, processSchedule, treeNodeBuffer, &
            MPI_parameter, MPIParam_increment

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
  ! same time. Most important attributes concerning the
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
    integer(i4)                                :: permind   ! index in node in original array, 
                                                            ! ToDo: because there is a mess at the moment
    type(treeNodeBuffer), pointer              :: values    ! a buffer that holds the values associated to that node
    type(treeNodeBuffer_dp), pointer           :: qTIN      ! a buffer that holds the values associated to that node
    type(treeNodeBuffer_dp), pointer           :: qTR       ! a buffer that holds the values associated to that node
    type(treeNodeBuffer_dp), pointer           :: C1        ! a buffer that holds the values associated to that node
    type(treeNodeBuffer_dp), pointer           :: C2        ! a buffer that holds the values associated to that node
    type(treeNodeBuffer_dp), pointer           :: qOut      ! a buffer that holds the values associated to that node

    type(ptrTreeNode)                          :: post      ! node downstream,
                                                            ! parent
                                                            ! null, if root
    integer(i4)                                :: Nprae     ! number of children
    type(ptrTreeNode),dimension(:),allocatable :: prae      ! array of children
    integer(i4)                                :: siz       ! size of subtree
    integer(i4)                                :: sizOrig   ! size of subtree before cut down
    integer(i4)                                :: sizUp     ! size of smallest
                                                            ! subtree > lowBound
                                                            ! upstream
    type(subtreeNode), pointer                 :: ST        ! data of a node, that is also a subtreetree node

    !ToDo: Do we need this three components? If so, array of three?
    integer(i4)                                :: NSTinBranch ! number of cut of subtrees in branch
    integer(i4)                                :: level     ! distance to root
    integer(i4)                                :: revLevel  ! distance from farthest leave

    logical                                    :: isIn      ! true if this is a halo leaf

  end type treeNode

  type subtreeNode
    integer(i4)                                :: indST     ! index in node in Subtree array
    type(ptrTreeNode)                          :: postST    ! next subtree downstream,
                                                            ! parent
    integer(i4)                                :: NpraeST   ! number of cut of subtrees in node
    type(ptrTreeNode),dimension(:),allocatable :: praeST    ! array of subtree children
    integer(i4)                                :: sizST     ! size of cut of subtree
    integer(i4), dimension(2)                  :: levelST   ! component 1: distance to root
                                                            ! component 2: distance from farthest leave
                                                            
    integer(i4), dimension(2)                  :: sched     ! component 1: process
                                                            ! component 2: time slot
  end type subtreeNode

  type treeNodeBuffer
    integer(i4),   dimension(:), allocatable  :: buffer
    type(MPI_Status)                          :: status
    type(MPI_Request)                         :: request
  end type treeNodeBuffer

  type treeNodeBuffer_dp
    real(dp),      dimension(:), allocatable  :: buffer
    type(MPI_Status)                          :: status
    type(MPI_Request)                         :: request
  end type treeNodeBuffer_dp

  type subtreeMeta
    integer(i4)        :: indST
    integer(i4)        :: iStart
    integer(i4)        :: iEnd
    integer(i4)        :: sizST
    integer(i4)        :: iInStart
    integer(i4)        :: iInEnd
    integer(i4)        :: nIn
    type(MPI_Status),  dimension(:),   allocatable :: statuses
    type(MPI_Request), dimension(:),   allocatable :: requests
  end type subtreeMeta

  type subtreeBuffer
    type(MPI_Status),  dimension(:),   allocatable :: statuses
    type(MPI_Request), dimension(:),   allocatable :: requests
  end type subtreeBuffer

  type processSchedule
     integer(i4)                                :: nTrees
     integer(i4), dimension(:), allocatable     :: trees
     integer(i4)                                :: overallSize
     integer(i4)                                :: overallInSize ! number of all halo nodes
  end type processSchedule

  type MPI_parameter
    type(MPI_Comm)                   :: comm                ! MPI communicator
    integer(i4)                      :: nproc
    integer(i4)                      :: rank
    integer(i4)                      :: bufferIndex
    logical                          :: buffered
    logical                          :: bufferWrite
    integer(i4)                      :: bufferLength
    integer(i4)                      :: lowBound
    integer(i4)                      :: lowBoundOMP
    contains
    procedure :: increment => MPIParam_increment
  end type MPI_parameter
  
  contains

  subroutine MPIParam_increment(this)
    class(MPI_parameter), intent(inout) :: this

    if (this%bufferIndex < this%bufferLength) then
      this%bufferIndex = this%bufferIndex + 1
      this%buffered = .false.
    else
      this%bufferIndex = 1
      this%buffered = .true.
      this%bufferWrite = .true.
    end if
  end subroutine MPIParam_increment

END MODULE mo_HRD_types

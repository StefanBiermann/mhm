!>       \file mo_mrm_routing.f90

!>       \brief Performs runoff routing for mHM at level L11.

!>       \details This module performs flood routing at a given time step
!>       through the stream network at level L11 to the sink cell.
!>       The Muskingum flood routing algorithm is used.

!>       \authors Luis Samaniego

!>       \date Dec 2012

! Modifications:
! Stephan Thober Aug 2015 - adapted to mRM

MODULE mo_mrm_routing

  ! This module performs runoff flood routing for mHM.

  ! Written Luis Samaniego, Dec 2012

  USE mo_kind, ONLY : i4, dp
  use mpi_f08

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mRM_routing, mRM_routing_par

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        mRM_routing

  !    PURPOSE
  !>       \brief route water given runoff

  !>       \details This routine first performs mpr for the routing variables
  !>       if required, then accumulates the runoff to the routing resolution
  !>       and eventually routes the water in a third step. The last step is
  !>       repeated multiple times if the routing timestep is smaller than
  !>       the timestep of the hydrological timestep

  !    INTENT(IN)
  !>       \param[in] "logical :: read_states"                            whether states are derived from restart file
  !>       \param[in] "integer(i4) :: processCase"                        Process switch for routing
  !>       \param[in] "real(dp), dimension(:) :: global_routing_param"    routing parameters
  !>       \param[in] "real(dp), dimension(:) :: L1_total_runoff"         total runoff from L1 grid cells
  !>       \param[in] "real(dp), dimension(:) :: L1_areaCell"             L1 cell area
  !>       \param[in] "integer(i4), dimension(:) :: L1_L11_Id"            L1 cell ids on L11
  !>       \param[in] "real(dp), dimension(:) :: L11_areaCell"            L11 cell area
  !>       \param[in] "integer(i4), dimension(:) :: L11_L1_Id"            L11 cell ids on L1
  !>       \param[in] "integer(i4), dimension(:) :: L11_netPerm"          L11 routing order
  !>       \param[in] "integer(i4), dimension(:) :: L11_fromN"            L11 source grid cell order
  !>       \param[in] "integer(i4), dimension(:) :: L11_toN"              L11 target grid cell order
  !>       \param[in] "integer(i4) :: L11_nOutlets"                       L11 number of outlets/sinks
  !>       \param[in] "integer(i4) :: timestep"                           simulation timestep in [h]
  !>       \param[in] "real(dp) :: tsRoutFactor"                          factor between routing timestep and
  !>       hydrological timestep
  !>       \param[in] "integer(i4) :: nNodes"                             number of nodes
  !>       \param[in] "integer(i4) :: nInflowGauges"                      number of inflow gauges
  !>       \param[in] "integer(i4), dimension(:) :: InflowGaugeIndexList" index list of inflow gauges
  !>       \param[in] "logical, dimension(:) :: InflowGaugeHeadwater"     flag for headwater cell of inflow gauge
  !>       \param[in] "integer(i4), dimension(:) :: InflowGaugeNodeList"  gauge node list at L11
  !>       \param[in] "real(dp), dimension(:) :: InflowDischarge"         inflowing discharge at discharge gauge at
  !>       current day
  !>       \param[in] "integer(i4) :: nGauges"                            number of recording gauges
  !>       \param[in] "integer(i4), dimension(:) :: gaugeIndexList"       index list for outflow gauges
  !>       \param[in] "integer(i4), dimension(:) :: gaugeNodeList"        gauge node list at L11
  !>       \param[in] "logical :: map_flag"                               flag indicating whether routing resolution
  !>       iscoarser than hydrologic resolution
  !>       \param[in] "real(dp), dimension(:) :: L11_length"              L11 link length
  !>       \param[in] "real(dp), dimension(:) :: L11_slope"               L11 slope
  !>       \param[in] "real(dp), dimension(:) :: L11_FracFPimp"           L11 fraction of flood plain with impervios
  !>       cover

  !    INTENT(INOUT)
  !>       \param[inout] "real(dp), dimension(:) :: L11_C1"         L11 muskingum parameter 1
  !>       \param[inout] "real(dp), dimension(:) :: L11_C2"         L11 muskingum parameter 2
  !>       \param[inout] "real(dp), dimension(:) :: L11_qOut"       total runoff from L11 grid cells
  !>       \param[inout] "real(dp), dimension(:, :) :: L11_qTIN"    L11 inflow to the reach
  !>       \param[inout] "real(dp), dimension(:, :) :: L11_qTR"     L11 routed outflow
  !>       \param[inout] "real(dp), dimension(:) :: L11_qMod"       modelled discharge at each grid cell
  !>       \param[inout] "real(dp), dimension(:) :: GaugeDischarge" modelled discharge at each gauge

  !    HISTORY
  !>       \authors Stephan Thober

  !>       \date Aug 2015

  ! Modifications:
  ! Stephan Thober Sep 2015 - using arguments instead of global variables
  ! Stephan Thober Sep 2015 - added variables for routing resolution higher than hydrologic resolution
  ! Stephan Thober May 2016 - added check whether gauge is actually inside modelling domain before copying simulated runoff
  ! Stephan Thober Nov 2016 - implemented second routing process i.e. adaptive timestep
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine mRM_routing(read_states, processCase, global_routing_param, L1_total_runoff, L1_areaCell, L1_L11_Id, &
                        L11_areaCell, L11_L1_Id, L11_netPerm, L11_fromN, L11_toN, L11_nOutlets, timestep, tsRoutFactor, &
                        nNodes, nInflowGauges, InflowGaugeIndexList, InflowGaugeHeadwater, InflowGaugeNodeList, &
                        InflowDischarge, nGauges, gaugeIndexList, gaugeNodeList, map_flag, L11_length, L11_slope, &
                        L11_FracFPimp, L11_C1, L11_C2, L11_qOut, L11_qTIN, L11_qTR, L11_qMod, GaugeDischarge)

    use mo_mrm_global_variables, only : is_start
    use mo_mrm_mpr, only : reg_rout

    implicit none

    ! whether states are derived from restart file
    logical, intent(in) :: read_states

    ! Process switch for routing
    integer(i4), intent(in) :: processCase

    ! routing parameters
    real(dp), dimension(:), intent(in) :: global_routing_param

    ! total runoff from L1 grid cells
    real(dp), dimension(:), intent(in) :: L1_total_runoff

    ! L1 cell area
    real(dp), dimension(:), intent(in) :: L1_areaCell

    ! L1 cell ids on L11
    integer(i4), dimension(:), intent(in) :: L1_L11_Id

    ! L11 cell area
    real(dp), dimension(:), intent(in) :: L11_areaCell

    ! L11 cell ids on L1
    integer(i4), dimension(:), intent(in) :: L11_L1_Id

    ! L11 routing order
    integer(i4), dimension(:), intent(in) :: L11_netPerm

    ! L11 source grid cell order
    integer(i4), dimension(:), intent(in) :: L11_fromN

    ! L11 target grid cell order
    integer(i4), dimension(:), intent(in) :: L11_toN

    ! L11 number of outlets/sinks
    integer(i4), intent(in) :: L11_nOutlets

    ! simulation timestep in [h]
    integer(i4), intent(in) :: timestep

    ! factor between routing timestep and hydrological timestep
    real(dp), intent(in) :: tsRoutFactor

    ! number of nodes
    integer(i4), intent(in) :: nNodes

    ! number of inflow gauges
    integer(i4), intent(in) :: nInflowGauges

    ! index list of inflow gauges
    integer(i4), dimension(:), intent(in) :: InflowGaugeIndexList

    ! flag for headwater cell of inflow gauge
    logical, dimension(:), intent(in) :: InflowGaugeHeadwater

    ! gauge node list at L11
    integer(i4), dimension(:), intent(in) :: InflowGaugeNodeList

    ! inflowing discharge at discharge gauge at current day
    real(dp), dimension(:), intent(in) :: InflowDischarge

    ! number of recording gauges
    integer(i4), intent(in) :: nGauges

    ! index list for outflow gauges
    integer(i4), dimension(:), intent(in) :: gaugeIndexList

    ! gauge node list at L11
    integer(i4), dimension(:), intent(in) :: gaugeNodeList

    ! flag indicating whether routing resolution iscoarser than hydrologic resolution
    logical, intent(in) :: map_flag

    ! L11 link length
    real(dp), dimension(:), intent(in) :: L11_length

    ! L11 slope
    real(dp), dimension(:), intent(in) :: L11_slope

    ! L11 fraction of flood plain with impervios cover
    real(dp), dimension(:), intent(in) :: L11_FracFPimp

    ! L11 muskingum parameter 1
    real(dp), dimension(:), intent(inout) :: L11_C1

    ! L11 muskingum parameter 2
    real(dp), dimension(:), intent(inout) :: L11_C2

    ! total runoff from L11 grid cells
    real(dp), dimension(:), intent(inout) :: L11_qOut

    ! L11 inflow to the reach
    real(dp), dimension(:, :), intent(inout) :: L11_qTIN

    ! L11 routed outflow
    real(dp), dimension(:, :), intent(inout) :: L11_qTR

    ! modelled discharge at each grid cell
    real(dp), dimension(:), intent(inout) :: L11_qMod

    ! modelled discharge at each gauge
    real(dp), dimension(:), intent(inout) :: GaugeDischarge

    integer(i4) :: gg

    integer(i4) :: tt

    ! number of routing loops
    integer(i4) :: rout_loop

    ! variable for accumulation
    real(dp), dimension(size(L11_qMod, dim = 1)) :: L11_qAcc


    if (is_start) then
      is_start = .false.
    end if

    ! this is using the sealed fraction for determining the routing parameters
    ! MPR has already been done
    if (processCase .eq. 1_i4 .AND. (.not. read_states)) then
      ! for a single node model run
      if (nNodes .GT. 1) then
        call reg_rout(global_routing_param, &
                L11_length, L11_slope, L11_FracFPimp(: nNodes - L11_nOutlets), &
                real(timeStep, dp), L11_C1(: nNodes - L11_nOutlets), L11_C2(: nNodes - L11_nOutlets))
      end if
    end if

    ! =====================================================================
    ! NOW, EXECUTE ROUTING
    ! ====================================================================
    ! calculate number of routing loops
    rout_loop = max(1_i4, nint(1._dp / tsRoutFactor))

    ! runoff accumulation from L1 to L11 level
    call L11_runoff_acc(L1_total_runoff, L1_areaCell, L1_L11_Id, &
            L11_areaCell, L11_L1_Id, timeStep, & ! Intent IN
            map_flag, & ! Intent IN
            L11_qOut) ! Intent OUT

    ! add inflow
    call add_inflow(nInflowGauges, &
            InflowGaugeIndexList, &
            InflowGaugeHeadwater, &
            InflowGaugeNodeList, &
            InflowDischarge, & ! Intent IN
            L11_qOUT) ! Intent INOUT

    ! for a single node model run
    if(nNodes .GT. 1) then
      ! routing multiple times if timestep is smaller than 1
      !
      L11_qAcc = 0._dp
      do tt = 1, rout_loop
        ! routing of water within river reaches
        call L11_routing(nNodes, nNodes - L11_nOutlets, &
                L11_netPerm, &
                L11_fromN, & ! Intent IN
                L11_toN, & ! Intent IN
                L11_C1, & ! Intent IN
                L11_C2, & ! Intent IN
                L11_qOut, & ! Intent IN
                nInflowGauges, & ! Intent IN
                InflowGaugeHeadwater, & ! Intent IN
                InflowGaugeNodeList, & ! Intent IN
                L11_qTIN, & ! Intent INOUT
                L11_qTR, & ! Intent INOUT
                L11_Qmod) ! Intent OUT
        ! accumulate values of individual subtimesteps
        L11_qAcc = L11_qAcc + L11_qMod
      end do
      ! calculate mean over routing period (timestep)
      L11_qMod = L11_qAcc / real(rout_loop, dp)
    else
      L11_Qmod = L11_qOUT
    end if

    !----------------------------------------------------------------------
    ! FOR STORING the optional arguments
    ! 
    ! FOR RUNOFF
    ! NOTE:: Node ID for a given gauging station is stored at gaugeindex's
    !        index in runoff. In consequence the gauges in runoff are 
    !        ordered corresponing to gauge%Q(:,:)
    !----------------------------------------------------------------------
    do gg = 1, nGauges
      GaugeDischarge(gaugeIndexList(gg)) = L11_Qmod(gaugeNodeList(gg))
    end do

  end subroutine mRM_routing

  subroutine mRM_routing_par(read_states, processCase, global_routing_param, L1_total_runoff, L1_areaCell, L1_L11_Id, &
                        L11_areaCell, L11_L1_Id, L11_netPerm, L11_fromN, L11_toN, L11_nOutlets, timestep, tsRoutFactor, &
                        nNodes, nInflowGauges, InflowGaugeIndexList, InflowGaugeHeadwater, InflowGaugeNodeList, &
                        InflowDischarge, nGauges, gaugeIndexList, gaugeNodeList, map_flag, L11_length, L11_slope, &
                        L11_FracFPimp, L11_buf_C1, L11_buf_C2, L11_buf_qOut, L11_qTIN, L11_qTR, L11_qMod, GaugeDischarge, &
                        iBasin, MPIparam, subtrees, nSubtrees, STmeta, permNodes, schedule)

    use mo_mrm_global_variables, only : is_start
    use mo_mrm_mpr, only : reg_rout
    use mo_HRD_types, only : MPI_parameter, ptrTreeNode, subtreeMeta, processSchedule, MPIParam_increment
    use mo_HRD_domain_decomposition, only : routing
    use mo_HRD_routing, only : muskignum_master_routing
    use mo_HRD_MPI_array_communication, only : distribute_array_dp, collect_array_dp, distribute_full_array_dp, collect_full_array_dp

    implicit none

    ! whether states are derived from restart file
    logical, intent(in) :: read_states

    ! Process switch for routing
    integer(i4), intent(in) :: processCase

    ! routing parameters
    real(dp), dimension(:), intent(in) :: global_routing_param

    ! total runoff from L1 grid cells
    real(dp), dimension(:), intent(in) :: L1_total_runoff

    ! L1 cell area
    real(dp), dimension(:), intent(in) :: L1_areaCell

    ! L1 cell ids on L11
    integer(i4), dimension(:), intent(in) :: L1_L11_Id

    ! L11 cell area
    real(dp), dimension(:), intent(in) :: L11_areaCell

    ! L11 cell ids on L1
    integer(i4), dimension(:), intent(in) :: L11_L1_Id

    ! L11 routing order
    integer(i4), dimension(:), intent(in) :: L11_netPerm

    ! L11 source grid cell order
    integer(i4), dimension(:), intent(in) :: L11_fromN

    ! L11 target grid cell order
    integer(i4), dimension(:), intent(in) :: L11_toN

    ! L11 number of outlets/sinks
    integer(i4), intent(in) :: L11_nOutlets

    ! simulation timestep in [h]
    integer(i4), intent(in) :: timestep

    ! factor between routing timestep and hydrological timestep
    real(dp), intent(in) :: tsRoutFactor

    ! number of nodes
    integer(i4), intent(in) :: nNodes

    ! number of inflow gauges
    integer(i4), intent(in) :: nInflowGauges

    ! index list of inflow gauges
    integer(i4), dimension(:), intent(in) :: InflowGaugeIndexList

    ! flag for headwater cell of inflow gauge
    logical, dimension(:), intent(in) :: InflowGaugeHeadwater

    ! gauge node list at L11
    integer(i4), dimension(:), intent(in) :: InflowGaugeNodeList

    ! inflowing discharge at discharge gauge at current day
    real(dp), dimension(:), intent(in) :: InflowDischarge

    ! number of recording gauges
    integer(i4), intent(in) :: nGauges

    ! index list for outflow gauges
    integer(i4), dimension(:), intent(in) :: gaugeIndexList

    ! gauge node list at L11
    integer(i4), dimension(:), intent(in) :: gaugeNodeList

    ! flag indicating whether routing resolution iscoarser than hydrologic resolution
    logical, intent(in) :: map_flag

    ! L11 link length
    real(dp), dimension(:), intent(in) :: L11_length

    ! L11 slope
    real(dp), dimension(:), intent(in) :: L11_slope

    ! L11 fraction of flood plain with impervios cover
    real(dp), dimension(:), intent(in) :: L11_FracFPimp

    ! L11 muskingum parameter 1
    real(dp), dimension(:, :), intent(inout) :: L11_buf_C1

    ! L11 muskingum parameter 2
    real(dp), dimension(:, :), intent(inout) :: L11_buf_C2

    ! total runoff from L11 grid cells
    real(dp), dimension(:, :), intent(inout) :: L11_buf_qOut

    ! L11 inflow to the reach
    real(dp), dimension(:, :), intent(inout) :: L11_qTIN

    ! L11 routed outflow
    real(dp), dimension(:, :), intent(inout) :: L11_qTR

    ! modelled discharge at each grid cell
    real(dp), dimension(:, :), intent(inout) :: L11_qMod

    ! modelled discharge at each gauge
    real(dp), dimension(:), intent(inout) :: GaugeDischarge

    integer(i4), intent(in) :: iBasin
    type(MPI_parameter), intent(in) :: MPIparam
    type(ptrTreeNode),     dimension(:), intent(in)    :: subtrees ! the array of
    integer(i4),                         intent(in)    :: nSubtrees
    type(subtreeMeta),     dimension(:), intent(inout) :: STmeta
    integer(i4),           dimension(:), intent(in)    :: permNodes
    type(processSchedule), dimension(:), intent(in)    :: schedule


    integer(i4) :: gg

    integer(i4) :: tt

    ! number of routing loops
    integer(i4) :: rout_loop

    ! variable for accumulation
    real(dp), dimension(size(L11_qMod, dim = 1), size(L11_qMod, dim=2)) :: L11_qAcc

    integer(i4)                                  :: ierror
    integer(i4)                                  :: iproc, ind

    real(dp), dimension(size(L11_buf_C1, dim=1)) :: L11_C1

    real(dp), dimension(size(L11_buf_C2, dim=1)) :: L11_C2

    real(dp), dimension(size(L11_buf_qOut, dim=1)) :: L11_qOut

    if (is_start) then
      is_start = .false.
    end if

    ! this is using the sealed fraction for determining the routing parameters
    ! MPR has already been done
    if (processCase .eq. 1_i4 .AND. (.not. read_states)) then
      ! for a single node model run
      if (nNodes .GT. 1) then
        call reg_rout(global_routing_param, &
                L11_length, L11_slope, L11_FracFPimp(: nNodes - L11_nOutlets), &
                real(timeStep, dp), L11_buf_C1(: nNodes - L11_nOutlets, MPIparam%bufferIndex), &
                                    L11_buf_C2(: nNodes - L11_nOutlets, MPIparam%bufferIndex))
        ! ToDo: Test this
        L11_buf_C1(: nNodes - L11_nOutlets,mod(MPIparam%bufferIndex+1,MPIparam%bufferLength)+1) = &
                                          L11_buf_C1(: nNodes - L11_nOutlets, MPIparam%bufferIndex)
        L11_buf_C2(: nNodes - L11_nOutlets,mod(MPIparam%bufferIndex+1,MPIparam%bufferLength)+1) = &
                                          L11_buf_C2(: nNodes - L11_nOutlets, MPIparam%bufferIndex)
      end if
    end if

    ! =====================================================================
    ! NOW, EXECUTE ROUTING
    ! ====================================================================
    ! calculate number of routing loops
    rout_loop = max(1_i4, nint(1._dp / tsRoutFactor))

    ! runoff accumulation from L1 to L11 level
    call L11_runoff_acc(L1_total_runoff, L1_areaCell, L1_L11_Id, &
            L11_areaCell, L11_L1_Id, timeStep, & ! Intent IN
            map_flag, & ! Intent IN
            L11_buf_qOut(:, MPIparam%bufferIndex)) ! Intent OUT

    ! add inflow
    call add_inflow(nInflowGauges, &
            InflowGaugeIndexList, &
            InflowGaugeHeadwater, &
            InflowGaugeNodeList, &
            InflowDischarge, & ! Intent IN
            L11_buf_qOUT(:, MPIparam%bufferIndex)) ! Intent INOUT

    ! for a single node model run
    if(nNodes .GT. 1) then
      ! routing multiple times if timestep is smaller than 1
      !
      do tt = 1, rout_loop
       ! L11_buf_C1(:, MPIparam%bufferIndex)   = L11_C1(:)
       ! L11_buf_C2(:, MPIparam%bufferIndex)   = L11_C2(:)
       ! L11_buf_qOut(:, MPIparam%bufferIndex) = L11_qOut(:)
        call MPIparam%increment()
        ! sending intent ins of L11_routing to all nodes ToDo: send arrays and buffer
        if (MPIparam%buffered) then
          do iproc = 1, MPIparam%nproc-1
            call MPI_Send(nInflowGauges, 1, MPI_INTEGER, iproc, 2, MPIparam%comm, ierror)
            call MPI_Send(InflowGaugeHeadwater, nInflowGauges, MPI_LOGICAL, iproc, 2, MPIparam%comm, ierror)
            call MPI_Send(InflowGaugeNodeList, nInflowGauges, MPI_INTEGER, iproc, 2, MPIparam%comm, ierror)
          end do
          ! ToDo: If C1, C2 were not constant maybe the sorting would be wrong
          call distribute_full_array_dp(iBasin, MPIparam, &
                                nSubtrees, STmeta,permNodes, schedule, L11_buf_C1(:, :))
          call distribute_full_array_dp(iBasin, MPIparam, &
                                nSubtrees, STmeta,permNodes, schedule, L11_buf_C2(:, :))
          call distribute_full_array_dp(iBasin, MPIparam, &
                                nSubtrees, STmeta,permNodes, schedule, L11_buf_qOut(:, :))
          call distribute_array_dp(iBasin, MPIparam%nproc, MPIparam%rank, MPIparam%comm, &
                                nSubtrees, STmeta,permNodes, schedule, L11_qTIN(:, 1))
          call distribute_array_dp(iBasin, MPIparam%nproc, MPIparam%rank, MPIparam%comm, &
                                nSubtrees, STmeta,permNodes, schedule, L11_qTR(:, 1))
          ! routing of water within river reaches
          call muskignum_master_routing(iBasin, MPIparam, subtrees, nSubtrees, STmeta, schedule)
          call MPI_Barrier(MPIparam%comm)
          ! backflow t-> t-1
          call collect_full_array_dp(iBasin, MPIParam, &
                                nSubtrees, STmeta, permNodes, schedule, L11_qTIN(:, :))
          call collect_full_array_dp(iBasin, MPIParam, &
                                nSubtrees, STmeta, permNodes, schedule, L11_qTR(:, :))
          call MPI_Barrier(MPIparam%comm)
          ! store generated discharge
          do iproc = 1, MPIparam%bufferLength
            write(*,*) L11_qTIN(1 : nNodes, iproc), iproc, MPIparam%bufferLength
            read(*,*)
          end do
          L11_qMod(1 : nNodes, :) = L11_qTIN(1 : nNodes, :)
          L11_qTIN(:, 1) = L11_qTIN(:, MPIparam%bufferLength)
          L11_qTR(:, 1)  = L11_qTR(:, MPIparam%bufferLength)
        end if
      end do
      ! ToDo: also when was buffered
      if (MPIparam%buffered) then
        tt = 0
        ind = 1
        do gg = 1, MPIparam%bufferLength
          L11_qAcc = 0._dp
          ! accumulate values of individual subtimesteps
          L11_qAcc(:, ind) = L11_qAcc(:, ind) + L11_qMod(:, gg)
          tt = tt + 1
          ! calculate mean over routing period (timestep) ToDo: handle
          ! bufferLength % routloop /= 0
          if (tt == rout_loop) then
            do tt = ind, gg
              L11_qMod(:, tt) = L11_qAcc(:, gg) / real(rout_loop, dp)
            end do
            ind = ind + 1
            tt = 0
          end if
        end do
      end if
    else
      ! ToDo: Check this case
      L11_Qmod(:, :) = L11_buf_qOUT(:, :)
    end if

    !----------------------------------------------------------------------
    ! FOR STORING the optional arguments
    ! 
    ! FOR RUNOFF
    ! NOTE:: Node ID for a given gauging station is stored at gaugeindex's
    !        index in runoff. In consequence the gauges in runoff are 
    !        ordered corresponing to gauge%Q(:,:)
    !----------------------------------------------------------------------
    do gg = 1, nGauges
      GaugeDischarge(gaugeIndexList(gg)) = L11_Qmod(gaugeNodeList(gg), 1)
    end do

  end subroutine mRM_routing_par

  ! ------------------------------------------------------------------

  !    NAME
  !        L11_runoff_acc

  !    PURPOSE
  !>       \brief total runoff accumulation at L11.

  !>       \details Upscales runoff in space from L1 to L11 if routing resolution
  !>       is higher than hydrology resolution (map_flag equals .true.) or
  !>       downscales runoff from L1 to L11 if routing resolution is lower
  !>       than hydrology resolution.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: qall"         total runoff L1 [mm tst-1]
  !>       \param[in] "real(dp), dimension(:) :: efecarea"     effective area in [km2] at Level 1
  !>       \param[in] "integer(i4), dimension(:) :: L1_L11_Id" L11 Ids mapped on L1
  !>       \param[in] "real(dp), dimension(:) :: L11_areacell" effective area in [km2] at Level 11
  !>       \param[in] "integer(i4), dimension(:) :: L11_L1_Id" L1 Ids mapped on L11
  !>       \param[in] "integer(i4) :: TS"                      time step in [s]
  !>       \param[in] "logical :: map_flag"                    Flag indicating whether routing resolution is higher than
  !>       hydrologic one

  !    INTENT(OUT)
  !>       \param[out] "real(dp), dimension(:) :: qAcc" aggregated runoff at L11 [m3 s-1]

  !    HISTORY
  !>       \authors Luis Samaniego

  !>       \date Jan 2013

  ! Modifications:
  ! Matthias Zink  Mar 2014 - added inflow from upstream areas
  ! Matthias Zink  Dec 2014 - adopted inflow gauges to ignore headwater cells
  ! Stephan Thober Sep 2015 - included downscaling of runoff
  ! Stephan Thober Feb 2016 - refactored upscaling of discharge from L1 to L11
  ! Stephan Thober Feb 2016 - refactored downscaling of discharge from L1 to L11
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  SUBROUTINE L11_runoff_acc(qAll, efecArea, L1_L11_Id, L11_areaCell, L11_L1_Id, TS, map_flag, qAcc)

    use mo_common_constants, only : HourSecs, nodata_dp

    implicit none

    ! total runoff L1 [mm tst-1]
    real(dp), intent(in), dimension(:) :: qall

    ! effective area in [km2] at Level 1
    real(dp), intent(in), dimension(:) :: efecarea

    ! L11 Ids mapped on L1
    integer(i4), intent(in), dimension(:) :: L1_L11_Id

    ! effective area in [km2] at Level 11
    real(dp), intent(in), dimension(:) :: L11_areacell

    ! L1 Ids mapped on L11
    integer(i4), intent(in), dimension(:) :: L11_L1_Id

    ! time step in [s]
    integer(i4), intent(in) :: TS

    ! Flag indicating whether routing resolution is higher than hydrologic one
    logical, intent(in) :: map_flag

    ! aggregated runoff at L11 [m3 s-1]
    real(dp), intent(out), dimension(:) :: qAcc

    integer(i4) :: k

    ! [s] time step
    real(dp) :: TST


    ! ------------------------------------------------------------------
    ! ACCUMULATION OF DISCHARGE TO A ROUTING CELL
    ! ------------------------------------------------------------------
    ! Hydrologic timestep in seconds
    TST = HourSecs * TS

    if (map_flag) then
      ! Estimate specific runoff at  L11
      ! NOTE:
      ! 1) Total discharge depth aggregated at L11 level [mm/TST]
      ! 2) Transform  depth [mm/TST] to discharge [m3/s]
      ! Total runoff should be divided by total_area to get
      ! specific discharge at L11. Then, to transform specific
      ! discharge from [mm/TST] to [m3/s], it should be multiplied by
      ! total_area [km2]*10^3 and divided by TST.
      ! Therefore, in this operation total_area cancels out.
      qAcc = 0._dp
      ! loop over high-resolution cells (L1) and add discharge to
      ! corresponding low-resolution cells (L11)
      do k = 1, size(qAll, 1)
        qAcc(L1_L11_Id(k)) = qAcc(L1_L11_Id(k)) + qAll(k) * efecArea(k)
      end do
      qAcc = qAcc * 1000.0_dp / TST
      !
    else
      ! initialize qout
      qAcc = nodata_dp
      do k = 1, size(qAcc, 1)
        ! map flux from coarse L1 resolution to fine L11 resolution
        qAcc(k) = qAll(L11_L1_Id(k))
      end do
      ! adjust flux by area cell
      qAcc(:) = qAcc(:) * L11_areaCell(:) * 1000.0_dp / TST
    end if

  END SUBROUTINE L11_runoff_acc

  ! ------------------------------------------------------------------

  !    NAME
  !        add_inflow

  !    PURPOSE
  !>       \brief
  !>          Adds inflow discharge to the runoff produced at the
  !>          cell where the inflow is occurring.
  !>       \details
  !>          If a inflow gauge is given, then this routine is adding the
  !>          values to the runoff produced at the grid cell where the
  !>          inflow is happening. The values are not directly added to the
  !>          river network. If this cell is not a headwater then the streamflow
  !>          produced upstream will be neglected.

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: nInflowGauges"                 [-] number of inflow points
  !>       \param[in] "integer(i4), dimension(:) :: InflowIndexList" [-] index of inflow points
  !>       \param[in] "logical, dimension(:) :: InflowHeadwater"     [-] if to consider headwater cells of inflow gauge
  !>       \param[in] "integer(i4), dimension(:) :: InflowNodeList"  [-]        L11 ID of inflow points
  !>       \param[in] "real(dp), dimension(:) :: QInflow"            [m3 s-1]   inflowing water

  !    INTENT(INOUT)
  !>       \param[inout] "real(dp), dimension(:) :: qOut" [m3 s-1] Series of attenuated runoff

  !    HISTORY
  !>       \authors Stephan Thober & Matthias Zink

  !>       \date Jul 2016

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine add_inflow(nInflowGauges, InflowIndexList, InflowHeadwater, InflowNodeList, QInflow, qOut)

    use mo_kind, only : dp, i4

    implicit none

    ! [-] number of inflow points
    integer(i4), intent(in) :: nInflowGauges

    ! [-] index of inflow points
    integer(i4), intent(in), dimension(:) :: InflowIndexList

    ! [-] if to consider headwater cells of inflow gauge
    logical, intent(in), dimension(:) :: InflowHeadwater

    ! [-]        L11 ID of inflow points
    integer(i4), intent(in), dimension(:) :: InflowNodeList

    ! [m3 s-1]   inflowing water
    real(dp), intent(in), dimension(:) :: QInflow

    ! [m3 s-1] Series of attenuated runoff
    real(dp), intent(inout), dimension(:) :: qOut

    integer(i4) :: ii


    ! discharge for inflow gauges (e.g. for missing upstream catchments) is added here
    ! should be put after UH attenuation because it is measured runoff at this cell 
    if (nInflowGauges .gt. 0) then
      do ii = 1, nInflowGauges
        if (InflowHeadwater(ii)) then
          ! add inflowing water to water produced by upstream/headwater cells
          qOut(InflowNodeList(ii)) = qOut(InflowNodeList(ii)) + QInflow(InflowIndexList(ii))
        else
          ! put only timeseries and cut upstream/headwater cells produced water for routing
          qOut(InflowNodeList(ii)) = QInflow(InflowIndexList(ii))
        end if
      end do
    end if
  end subroutine add_inflow

  ! ------------------------------------------------------------------

  !    NAME
  !        L11_routing

  !    PURPOSE
  !>       \brief Performs runoff routing for mHM at L11 upscaled network
  !>       (\ref fig_routing "Routing Network").
  !>       \details
  !>       Hydrograph routing is carried out with the Muskingum algorithm
  !>       \cite CMM1988.  This simplification of the St. Venant
  !>       equations is justified in mHM because the potential areas of
  !>       application of this model would hardly exhibit abruptly
  !>       changing hydrographs with supercritical flows.  The discharge
  !>       leaving the river reach located on cell \f$ i \f$ \f$
  !>       Q_{i}^{1}(t) \f$ at time step \f$ t \f$ can be determined by
  !>       \f[ Q_{i}^{1}(t) =  Q_{i}^{1}(t-1)
  !>       + c_{1} \left( Q_{i}^{0}(t-1) - Q_{i}^{1}(t-1) \right)
  !>       + c_{2} \left( Q_{i}^{0}(t)   - Q_{i}^{0}(t-1) \right) \f]
  !>       with
  !>       \f[  Q_{i}^{0}(t) = Q_{i'}(t) + Q_{i'}^{1}(t) \f]
  !>       \f[ c_{1}= \frac{\Delta t} { \kappa (1- \xi ) + \frac{\Delta t}{2} } \f]
  !>       \f[ c_{2}= \frac{ \frac{\Delta t}{2} - \kappa \xi} { \kappa (1- \xi)
  !>       + \frac{\Delta t}{2} } \f]
  !>       where
  !>       \f$ Q_{i}^{0} \f$ and \f$ Q_{i}^{1} \f$ denote the discharge
  !>       entering and leaving the river reach located on cell \f$ i \f$
  !>       respectively.
  !>       \f$ Q_{i'} \f$ is the contribution from the upstream cell \f$
  !>       i'\f$.
  !>       \f$ \kappa \f$ Muskingum travel time parameter.
  !>       \f$ \xi \f$ Muskingum attenuation parameter.
  !>       \f$ \Delta t \f$ time interval in hours.
  !>       \f$ t \f$ Time index for each \f$ \Delta t \f$ interval.
  !>       To improve performance, a routing sequence "netPerm" is
  !>       required. This permutation is determined in the mo_init_mrm
  !>       routine.

  !>       \details TODO: add description

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: nNodes"                       number of network nodes = nCells1
  !>       \param[in] "integer(i4) :: nLinks"                       number of stream segment (reaches)
  !>       \param[in] "integer(i4), dimension(:) :: netPerm"        routing order of a given basin (permutation)
  !>       \param[in] "integer(i4), dimension(:) :: netLink_fromN"  from node
  !>       \param[in] "integer(i4), dimension(:) :: netLink_toN"    to node
  !>       \param[in] "real(dp), dimension(:) :: netLink_C1"        routing parameter  C1 (\cite CMM1988 p. 25-41)
  !>       \param[in] "real(dp), dimension(:) :: netLink_C2"        routing parameters C2 (id)
  !>       \param[in] "real(dp), dimension(:) :: netNode_qOUT"      Total outflow from cells (given basin) L11 at time
  !>       tt in [m3 s-1]
  !>       \param[in] "integer(i4) :: nInflowGauges"                [-]      number of inflow points
  !>       \param[in] "logical, dimension(:) :: InflowHeadwater"    [-]      if to consider headwater cells of inflow
  !>       gauge
  !>       \param[in] "integer(i4), dimension(:) :: InflowNodeList" [-]      L11 ID of inflow points

  !    INTENT(INOUT)
  !>       \param[inout] "real(dp), dimension(:, :) :: netNode_qTIN" [m3 s-1] Total inputs at t-1 and t
  !>       \param[inout] "real(dp), dimension(:, :) :: netNode_qTR"  [m3 s-1] Transformed outflow leaving
  !>       node I (Muskingum)

  !    INTENT(OUT)
  !>       \param[out] "real(dp), dimension(nNodes) :: netNode_Qmod" [m3 s-1] Simulated routed discharge

  !    HISTORY
  !>       \authors Luis Samaniego

  !>       \date Dec 2005

  ! Modifications:
  ! Luis Samaniego Feb 2008 - routing module (cells)
  ! Rohini Kumar   Aug 2011 - vector version of mHM-UFZ
  !                Nov 2011 - parallel version
  ! Luis Samaniego Jan 2013 - modularization, documentation
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine L11_routing(nNodes, nLinks, netPerm, netLink_fromN, netLink_toN, netLink_C1, netLink_C2, netNode_qOUT, &
                        nInflowGauges, InflowHeadwater, InflowNodeList, netNode_qTIN, netNode_qTR, netNode_Qmod)
    implicit none

    ! number of network nodes = nCells1
    integer(i4), intent(in) :: nNodes

    ! number of stream segment (reaches)
    integer(i4), intent(in) :: nLinks

    ! routing order of a given basin (permutation)
    integer(i4), dimension(:), intent(in) :: netPerm

    ! from node
    integer(i4), dimension(:), intent(in) :: netLink_fromN

    ! to node
    integer(i4), dimension(:), intent(in) :: netLink_toN

    ! routing parameter  C1 (\cite CMM1988 p. 25-41)
    real(dp), dimension(:), intent(in) :: netLink_C1

    ! routing parameters C2 (id)
    real(dp), dimension(:), intent(in) :: netLink_C2

    ! Total outflow from cells (given basin) L11 at time tt in [m3 s-1]
    real(dp), dimension(:), intent(in) :: netNode_qOUT

    ! [-]      number of inflow points
    integer(i4), intent(in) :: nInflowGauges

    ! [-]      if to consider headwater cells of inflow gauge
    logical, dimension(:), intent(in) :: InflowHeadwater

    ! [-]      L11 ID of inflow points
    integer(i4), dimension(:), intent(in) :: InflowNodeList

    ! [m3 s-1] Total inputs at t-1 and t
    real(dp), dimension(:, :), intent(inout) :: netNode_qTIN

    ! [m3 s-1] Transformed outflow leaving
    ! node I (Muskingum)
    real(dp), dimension(:, :), intent(inout) :: netNode_qTR

    ! [m3 s-1] Simulated routed discharge
    real(dp), dimension(nNodes), intent(out) :: netNode_Qmod

    integer(i4) :: i, k, iNode, tNode

    ! current routing state (2)
    integer(i4), parameter :: IT = 2

    ! past routing state (1)
    integer(i4), parameter :: IT1 = 1


    ! Entry value for the auxiliary vectors
    !   netNode_qTIN(iNode,:)
    !   netNode_qTR(iNode,:)
    ! which store current and past states of
    ! incoming and outgoing of discharge at iNode
    !--------------------------------------------------------------------------
    !                             Muskingum Flood Routing
    !--------------------------------------------------------------------------
    ! initialize total input at point time IT in all nodes
    netNode_qTIN(:, IT) = 0.0_dp
    !--------------------------------------------------------------------------
    ! Links in sequential mode .... with single node
    !--------------------------------------------------------------------------
    ! ST - decent parallelization has to be done!!!
    !!$OMP parallel
    !!$OMP do private( i, inode, tnode)
    do k = 1, nLinks
      ! get LINK routing order -> i
      i = netPerm(k)
      iNode = netLink_fromN(i)
      tNode = netLink_toN(i)

      ! accumulate all inputs in iNode
      netNode_qTIN(iNode, IT) = netNode_qTIN(iNode, IT) + netNode_qOUT(iNode)

     ! write(0,*) '----', netNode_qTR(iNode, IT1), netNode_qTIN(iNode, IT1), netNode_qTIN(iNode, IT), netLink_C1(i), netLink_C2(i)
     ! write(0,*) '----', netNode_qTIN(iNode, IT), netLink_C1(i), netLink_C2(i)
      ! routing iNode
      netNode_qTR(iNode, IT) = netNode_qTR(iNode, IT1)                               &
              + netLink_C1(i) * (netNode_qTIN(iNode, IT1) - netNode_qTR (iNode, IT1)) &
              + netLink_C2(i) * (netNode_qTIN(iNode, IT) - netNode_qTIN(iNode, IT1))

      ! check if the inflow from upstream cells should be deactivated
      if (nInflowGauges .GT. 0) then
        do i = 1, nInflowGauges
          ! check if downstream Node (tNode) is inflow gauge and headwaters should be ignored
          if ((tNode == InflowNodeList(i)) .AND. (.NOT. InflowHeadwater(i))) netNode_qTR(iNode, IT) = 0.0_dp
        end do
      end if

      ! add routed water to downstream node
      netNode_qTIN(tNode, IT) = netNode_qTIN(tNode, IT) + netNode_qTR(iNode, IT)
    !  write(0,*) '----', netNode_qTR(iNode, IT), netNode_qTIN(iNode, IT), tNode
    !  write(0,*) iNode
    end do
    !!$OMP end do
    !!$OMP end parallel

    !--------------------------------------------------------------------------
    ! Accumulate all inputs in tNode (netNode_qOUT) ONLY for last link
    !--------------------------------------------------------------------------
    tNode = netLink_toN(netPerm(nLinks))
    netNode_qTIN(tNode, IT) = netNode_qTIN(tNode, IT) + netNode_qOUT(tNode)

    !--------------------------------------------------------------------------
    ! save modeled discharge at time step tt then shift flow storages
    ! (NOTE aggregation to daily values to be done outside)
    !--------------------------------------------------------------------------
    ! !!$OMP parallel
    ! store generated discharge
    netNode_Qmod(1 : nNodes) = netNode_qTIN(1 : nNodes, IT)
    ! backflow t-> t-1
    netNode_qTR(1 : nNodes, IT1) = netNode_qTR(1 : nNodes, IT)
    netNode_qTIN(1 : nNodes, IT1) = netNode_qTIN(1 : nNodes, IT)
    ! !!$OMP end parallel

  end subroutine L11_routing


END MODULE mo_mrm_routing

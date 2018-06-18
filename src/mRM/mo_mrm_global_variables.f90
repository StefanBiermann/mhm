!>       \file mo_mrm_global_variables.f90

!>       \brief Global variables for mRM only

!>       \details 

!>       \authors s Luis Samaniego, Stephan Thober

!>       \date Aug 2015

! Modifications:
! Robert Schweppe Dec 2017 - merged duplicated variables with mhm into common variables

module mo_mrm_global_variables

  use mo_kind, only : i4, i8, dp
  use mo_mrm_constants, only : nOutFlxState
  use mo_common_variables, only : Grid, GridRemapper

  implicit none

  PUBLIC :: gaugingStation

  ! -------------------------------------------------------------------
  ! General variables
  ! -------------------------------------------------------------------
  logical :: is_start              ! flag for first timestep for mpr

  ! -------------------------------------------------------------------
  ! DEFINE OUTPUTS
  ! -------------------------------------------------------------------
  !
  integer(i4) :: timeStep_model_outputs_mrm ! timestep for writing model outputs

  logical, dimension(nOutFlxState) :: outputFlxState_mrm         ! Define model outputs see "mhm_outputs.nml"
  !                                                            dim1 = number of output variables to be written

  ! ------------------------------------------------------------------
  ! DIRECTORIES
  ! ------------------------------------------------------------------
  ! has the dimension of nBasins
  character(256), dimension(:), allocatable, public :: dirGauges ! Directory where discharge files are located
  character(256), dimension(:), allocatable, public :: dirTotalRunoff ! Directory where simulated total runoff files are located
  character(256), public :: filenameTotalRunoff ! Filename of simulated total runoff file
  character(256), public :: varnameTotalRunoff ! variable name of total runoff


  ! -------------------------------------------------------------------
  ! GRID description
  ! -------------------------------------------------------------------
  type(Grid), dimension(:), allocatable, target, public :: level11 ! Reference of the routing variables
  type(GridRemapper), dimension(:), allocatable, public :: l0_l11_remap ! grid information at runoff level
  type(GridRemapper), dimension(:), allocatable, public :: l1_l11_remap ! grid information at runoff level


  ! -----------------------------------------------------------------
  ! RUNOFF variable
  ! -----------------------------------------------------------------
  real(dp), dimension(:, :), allocatable, public :: mRM_runoff ! variable containing runoff for each basin and gauge

  ! -----------------------------------------------------------------
  ! GAUGED station data
  ! -----------------------------------------------------------------
  integer(i4), public :: nGaugesTotal ! Number of evaluation gauges for all basins 
  integer(i4), public :: nInflowGaugesTotal ! Number of evaluation gauges for all basins 
  integer(i4), public :: nMeasPerDay ! Number of observations per day,
  !                                  ! e.g. 24 -> hourly discharge, 1 -> daily discharge
  type gaugingStation
    integer(i4), dimension(:), allocatable :: basinId ! Basin Id
    integer(i4), dimension(:), allocatable :: gaugeId ! Gauge Id (e.g. 0000444)
    character(256), dimension(:), allocatable :: fname ! Name runoff file
    real(dp), dimension(:, :), allocatable :: Q ! [m3 s-1] observed daily mean discharge (simPer)
    !                                          ! dim1=number observations, dim2=number of gauges
  end type gaugingStation
  type(gaugingStation), public :: gauge ! Gauging station information
  type(gaugingStation), public :: InflowGauge ! inflow gauge information

  ! -------------------------------------------------------------------
  ! BASIN general description
  ! -------------------------------------------------------------------
  type basinInfo_mRM
    ! dim1 = maximum number of gauges in a given basin
    ! discharge measurement gauges
    integer(i4) :: nGauges        ! Number of gauges within a basin
    integer(i4), dimension(:), allocatable :: gaugeIdList    ! Gauge Id list (e.g. 0000444 0000445)
    integer(i4), dimension(:), allocatable :: gaugeIndexList ! Gauge index list (e.g. 1 for 00444, 2 for 00445)
    integer(i4), dimension(:), allocatable :: gaugeNodeList  ! Gauge node list at L11

    ! discharge inflow gauges (e.g if headwar bsins are missing)
    integer(i4) :: nInflowGauges        ! Number of gauges within a basin
    integer(i4), dimension(:), allocatable :: InflowGaugeIdList    ! Gauge Id list (e.g. 0000444 0000445)
    integer(i4), dimension(:), allocatable :: InflowGaugeIndexList ! Gauge index list (e.g. 1 for 00444, 2 for 00445)
    integer(i4), dimension(:), allocatable :: InflowGaugeNodeList  ! Gauge node list at L11
    logical, dimension(:), allocatable :: InflowGaugeHeadwater ! if headwater cells of inflow gauge will be considered

    ! basin outlet
    ! TODO: move this out of here since it is mrm_net_startup relevant only for domain0
    integer(i4) :: L0_Noutlet
    integer(i4), dimension(:), allocatable :: L0_rowOutlet   ! Outlet locations in L0
    integer(i4), dimension(:), allocatable :: L0_colOutlet   ! Outlet locations in L0
  end type basinInfo_mRM

  ! dim1 = basinId
  type(basinInfo_mRM), dimension(:), allocatable, public, target :: basin_mrm ! Basin structure
  ! -------------------------------------------------------------------
  ! L0 DOMAIN description -> those are needed within the mrm_net_startup routines only
  ! TODO: deallocate when net_startup is done
  ! -------------------------------------------------------------------
  ! dim1 = number grid cells
  ! input data - morphological variables
  integer(i4), public, dimension(:), allocatable :: L0_gaugeLoc ! Location of gauges within the catchment
  integer(i4), public, dimension(:), allocatable :: L0_InflowGaugeLoc ! Location of inflow gauges within catchment
  integer(i4), public, dimension(:), allocatable :: L0_fAcc ! Flow accumulation
  integer(i4), public, dimension(:), allocatable :: L0_fDir ! Flow direction (standard ArcGIS)
  !
  ! mRM derived variables
  ! dim1 = number grid cells L0
  integer(i4), public, dimension(:), allocatable :: L0_draSC      ! Index of draining cell of each sub catchment
  !                                                               ! i.e. a routing cell L11
  integer(i4), public, dimension(:), allocatable :: L0_draCell    ! Draining cell id at L11 of ith cell of L0
  integer(i4), public, dimension(:), allocatable :: L0_streamNet  ! Stream network
  integer(i4), public, dimension(:), allocatable :: L0_floodPlain ! Floodplains of stream i

  ! -------------------------------------------------------------------
  ! L1 DOMAIN description
  ! -------------------------------------------------------------------
  ! dim1 = number grid cells L1
  integer(i4), public, dimension(:), allocatable :: L11_L1_ID  ! Mapping of L11 Id on L1
  ! -------------------------------------------------------------------
  ! L1 variables
  ! -------------------------------------------------------------------
  ! dim1 = number grid cells L1
  ! dim2 = number of timesteps
  real(dp), public, dimension(:, :), allocatable :: L1_total_runoff_in

  ! -------------------------------------------------------------------
  ! L11 DOMAIN description
  ! -------------------------------------------------------------------
  ! dim1 = number grid cells L11
  ! dim2 = 2
  integer(i4), public, dimension(:), allocatable :: L1_L11_ID  ! Mapping of L1 Id on L11
  integer(i4), public, dimension(:), allocatable :: L11_fDir ! Flow direction (standard notation)
  integer(i4), public, dimension(:), allocatable :: L11_nOutlets

  ! Constants
  ! dim1 = number grid cells L11
  integer(i4), public, dimension(:), allocatable :: L11_rowOut ! Grid vertical location of the Outlet
  integer(i4), public, dimension(:), allocatable :: L11_colOut ! Grid horizontal location  of the Outlet

  ! -------------------------------------------------------------------
  ! L11 NETWORK description
  ! -------------------------------------------------------------------
  ! Fluxes
  ! dim1 = number grid cells L11
  ! dim2 = 2
  real(dp), public, dimension(:), allocatable :: L11_Qmod        ! [m3 s-1] Simulated discharge
  real(dp), public, dimension(:), allocatable :: L11_qOUT        ! [m3 s-1] Total outflow from cells L11 at time tt
  real(dp), public, dimension(:, :), allocatable :: L11_qTIN        !          Total discharge inputs at t-1 and t
  real(dp), public, dimension(:, :), allocatable :: L11_qTR         !          Routed outflow leaving a node

  integer(i4), public, dimension(:), allocatable :: L11_fromN       !         From node (sinks are at the end)
  integer(i4), public, dimension(:), allocatable :: L11_toN         !         To node (sinks are at the end)
  integer(i4), public, dimension(:), allocatable :: L11_netPerm     !         Routing sequence (permutation of L11_rOrder)
  integer(i4), public, dimension(:), allocatable :: L11_fRow        !         From row in L0 grid
  integer(i4), public, dimension(:), allocatable :: L11_fCol        !         From col in L0 grid
  integer(i4), public, dimension(:), allocatable :: L11_tRow        !         To row in L0 grid
  integer(i4), public, dimension(:), allocatable :: L11_tCol        !         To col in L0 grid
  integer(i4), public, dimension(:), allocatable :: L11_rOrder      !         Network routing order
  integer(i4), public, dimension(:), allocatable :: L11_label       !         Label Id [0='', 1=HeadWater, 2=Sink]
  logical, public, dimension(:), allocatable :: L11_sink        !         .true. if sink node reached
  real(dp), public, dimension(:), allocatable :: L11_length      ! [m]     Total length of river link
  real(dp), public, dimension(:), allocatable, target :: L11_aFloodPlain ! [m2]    Area of the flood plain
  !                                                                  !         impervious cover
  real(dp), public, dimension(:), allocatable :: L11_slope       ! [1]     Average slope of river link

  ! Parameters
  ! dim1 = number grid cells L11
  real(dp), public, dimension(:, :), allocatable :: L11_nLinkFracFPimp  !     fraction of impervious lcover at the floodplain
  !         dims: nLink, LCYearID

  real(dp), public, dimension(:), allocatable :: L11_K           ! [d]     kappa: Muskingum travel time parameter.
  real(dp), public, dimension(:), allocatable :: L11_xi          ! [1]     xi:    Muskingum diffusion parameter
  !                                                                  !                (attenuation).
  real(dp), public, dimension(:), allocatable :: L11_tsRout      ! [s]     Routing timestep
  real(dp), public, dimension(:), allocatable :: L11_C1          ! [-]     Routing parameter C1=f(K,xi, DT) (Chow, 25-41)
  real(dp), public, dimension(:), allocatable :: L11_C2          ! [-]     Routing parameter C2 (")

end module mo_mrm_global_variables

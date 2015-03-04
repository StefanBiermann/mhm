!> \file mo_read_optional_data.f90

!> \brief Read optional data for mHM calibration.

!> \details Data have to be provided in resolution of the hydrology.

!> \authors Matthias Zink
!> \date Mar 2015

MODULE mo_read_optional_data

  USE mo_kind, ONLY: i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_soil_moisture ! MZMZMZ
 
  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         read_soil_moisture ! MZMZMZMZ edit docu
  
  !     PURPOSE
  !>        \brief read soil moisture data for calibration

  !>        \details Prepare meteorological forcings data for a given variable.
  !>                 Internally this subroutine calls another routine meteo_wrapper   
  !>                 for different meterological variables

  !     CALLING SEQUENCE

  !     INTENT(IN)
  !>        \param[in] "integer(i4)              :: iBasin"        Basin Id

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None


  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS

  !     EXAMPLE

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Matthias Zink
  !>        \date Mar 2015

  subroutine read_soil_moisture(iBasin)
    use mo_message,          only: message
    use mo_string_utils,     only: num2str
    use mo_init_states,      only: get_basin_info
    use mo_read_meteo,       only: read_meteo_nc
    use mo_timer,            only:                         &
         timer_start, timer_stop, timer_get, timer_clear     ! Timing of processes
    use mo_append,           only: append                    ! append data
    use mo_mhm_constants,    only: nodata_dp
    use mo_global_variables, only:                         &
         dirSoil_moisture,                                 & ! directory of meteo input
         simPer,                                           & ! chunk read in config                           
         L1_sm, L1_sm_mask                                     ! soil mositure data and mask

    implicit none

    integer(i4), intent(in)  :: iBasin                         ! Basin Id ! MZMZMZMZ data packing

    ! local variables
    integer(i4)                             :: nTimeSteps,t    ! loop  vars packing L1_data to L1_data_packed
    integer(i4)                             :: nrows1, ncols1  ! level 1 number of culomns and rows
    logical, dimension(:,:), allocatable    :: mask1           ! mask of level 1 for packing
    integer(i4)                             :: ncells1         ! ncells1 of level 1
    real(dp), dimension(:,:,:), allocatable :: L1_data         ! data at level-1
    real(dp), dimension(:,:), allocatable   :: L1_data_packed  ! packed data at level-1 from 3D to 2D
    logical,  dimension(:,:,:), allocatable :: L1_mask         ! mask at level-1
    logical,  dimension(:,:), allocatable   :: L1_mask_packed  ! packed mask at level-1 from 3D to 2D
    
    ! get basic basin information at level-1
    call get_basin_info( iBasin, 1, nrows1, ncols1, nCells=nCells1, mask=mask1 ) 
       
    !  basin characteristics and read meteo header
    call message( '  Reading optional data for basin:           ', trim(adjustl(num2str(iBasin))),' ...')
    call timer_start(1)

    call read_meteo_nc( dirSoil_moisture(iBasin), nRows1, nCols1, simPer(iBasin), trim('sm'), L1_data, mask1, &
         nctimestep=-2, nocheck=.TRUE., maskout=L1_mask) !MZMZMZMZ

    ! pack variables
    nTimeSteps = size(L1_data, 3)
    allocate( L1_data_packed(nCells1, nTimeSteps))
    allocate( L1_mask_packed(nCells1, nTimeSteps))
    do t = 1, nTimeSteps
       L1_data_packed(:,t) = pack( L1_data(:,:,t), MASK=mask1(:,:) ) 
       L1_mask_packed(:,t) = pack( L1_mask(:,:,t), MASK=mask1(:,:) ) 
    end do
    
    ! append
    call append( L1_sm,      L1_data_packed(:,:), fill_value=nodata_dp )
    call append( L1_sm_mask, L1_mask_packed(:,:), fill_value=.FALSE. )

    !free space
    deallocate(L1_data, L1_data_packed) 

    call timer_stop(1)
    call message('    in ', trim(num2str(timer_get(1),'(F9.3)')), ' seconds.')
    call timer_clear(1)
    
  end subroutine read_soil_moisture
  
END MODULE mo_read_optional_data

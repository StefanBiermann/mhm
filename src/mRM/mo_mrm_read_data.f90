!>       \file mo_mrm_read_data.f90

!>       \brief This module contains all routines to read mRM data from file.
!>       \details

!>       \details TODO: add description

!>       \authors Stephan Thober

!>       \date Aug 2015

! Modifications:

module mo_mrm_read_data
  use mo_kind, only : i4, dp
  implicit none
  public :: mrm_read_L0_data
  public :: mrm_read_discharge
  public :: mrm_read_total_runoff
  public :: mrm_read_bankfull_runoff
  private
contains
  ! ------------------------------------------------------------------

  !    NAME
  !        mrm_read_L0_data

  !    PURPOSE
  !>       \brief read L0 data from file

  !>       \details With the exception of L0_mask, L0_elev, and L0_LCover, all
  !>       L0 variables are read from file. The former three are only read if they
  !>       are not provided as variables.

  !    INTENT(IN)
  !>       \param[in] "logical :: do_reinit"
  !>       \param[in] "logical :: do_readlatlon"
  !>       \param[in] "logical :: do_readlcover"

  !    HISTORY
  !>       \authors Juliane Mai, Matthias Zink, and Stephan Thober

  !>       \date Aug 2015

  ! Modifications:
  ! Stephan Thober Sep 2015 - added L0_mask, L0_elev, and L0_LCover
  ! Robert Schweppe Jun 2018 - refactoring and reformatting
  ! Stephan Thober Jun 2018 - including varying celerity functionality

  subroutine mrm_read_L0_data(do_reinit, do_readlatlon, do_readlcover)

    use mo_append, only : append
    use mo_common_constants, only : nodata_i4, nodata_dp
    use mo_common_read_data, only : read_dem, read_lcover
    use mo_common_variables, only : Grid, L0_Basin, L0_LCover, dirMorpho, level0, nBasins, processMatrix
    use mo_mpr_file, only: file_slope, uslope
    use mo_mpr_global_variables, only: L0_slope
    use mo_message, only : message
    use mo_mrm_file, only : file_facc, file_fdir, &
                            file_gaugeloc, ufacc, ufdir, ugaugeloc
    use mo_mrm_global_variables, only : L0_InflowGaugeLoc, L0_fAcc, L0_fDir, L0_gaugeLoc, basin_mrm
    use mo_read_latlon, only : read_latlon
    use mo_read_spatial_data, only : read_spatial_data_ascii
    use mo_string_utils, only : num2str

    implicit none

    logical, intent(in) :: do_reinit

    logical, intent(in) :: do_readlatlon

    logical, intent(in) :: do_readlcover

    integer(i4) :: iBasin

    integer(i4) :: iVar

    integer(i4) :: iGauge

    character(256) :: fname

    integer(i4) :: nunit

    integer(i4), dimension(:, :), allocatable :: data_i4_2d

    real(dp), dimension(:, :), allocatable :: data_dp_2d

    integer(i4), dimension(:, :), allocatable :: dataMatrix_i4

    logical, dimension(:, :), allocatable :: mask_2d

    logical, dimension(:, :), allocatable :: mask_global

    type(Grid), pointer :: level0_iBasin => null()


    ! ************************************************
    ! READ SPATIAL DATA FOR EACH BASIN
    ! ************************************************
    if (do_reinit) then
      call read_dem()
    end if

    if (do_readlcover .and. processMatrix(8, 1) .eq. 1) then
      call read_lcover()
    else if (do_readlcover .and. ((processMatrix(8, 1) .eq. 2) .or. (processMatrix(8, 1) .eq. 3))) then
      allocate(dataMatrix_i4(count(mask_global), 1))
      dataMatrix_i4 = nodata_i4
      call append(L0_LCover, dataMatrix_i4)
      ! free memory
      deallocate(dataMatrix_i4)
    end if

    do iBasin = 1, nBasins

      level0_iBasin => level0(L0_Basin(iBasin))

      ! check whether L0 data is shared
      if (iBasin .gt. 1) then
        if (L0_Basin(iBasin) .eq. L0_Basin(iBasin - 1)) then
          !
          call message('      Using data of basin ', &
                  trim(adjustl(num2str(L0_Basin(iBasin)))), ' for basin: ',&
                  trim(adjustl(num2str(iBasin))), '...')
          cycle
          !
        end if
      end if
      !
      call message('      Reading data for basin: ', trim(adjustl(num2str(iBasin))), ' ...')

      if (do_readlatlon) then
        ! read lat lon coordinates of each basin
        call read_latlon(iBasin, "lon_l0", "lat_l0", "level0", level0_iBasin)
      end if

      ! read fAcc, fDir, gaugeLoc
      do iVar = 1, 4
        select case (iVar)
        case(1) ! flow accumulation
          fName = trim(adjustl(dirMorpho(iBasin))) // trim(adjustl(file_facc))
          nunit = ufacc
        case(2) ! flow direction
          fName = trim(adjustl(dirMorpho(iBasin))) // trim(adjustl(file_fdir))
          nunit = ufdir
        case(3) ! location of gauging stations
          fName = trim(adjustl(dirMorpho(iBasin))) // trim(adjustl(file_gaugeloc))
          nunit = ugaugeloc
       case(4)
          fName = trim(adjustl(dirMorpho(iBasin))) // trim(adjustl(file_slope))
          nunit = uslope
       end select

       if (iVar .le. 3) then
          !
          ! reading and transposing
          call read_spatial_data_ascii(trim(fName), nunit, &
               level0_iBasin%nrows, level0_iBasin%ncols, level0_iBasin%xllcorner, &
               level0_iBasin%yllcorner, level0_iBasin%cellsize, data_i4_2d, mask_2d)

          ! put global nodata value into array (probably not all grid cells have values)
          data_i4_2d = merge(data_i4_2d, nodata_i4, mask_2d)
       end if

#ifndef MRM2MHM
       ! only read slope if not coupled to mHM and using third process Matrix
       if ((iVar .eq. 4) .and. (processMatrix(8, 1) .eq. 3)) then
          ! reading
          call read_spatial_data_ascii(trim(fName), nunit, &
               level0_iBasin%nrows, level0_iBasin%ncols, level0_iBasin%xllcorner, &
               level0_iBasin%yllcorner, level0_iBasin%cellsize, data_dp_2d, mask_2d)
          ! put global nodata value into array (probably not all grid cells have values)
          data_dp_2d = merge(data_dp_2d, nodata_dp, mask_2d)
          ! put data in variable
          call append(L0_slope, pack(data_dp_2d, level0_iBasin%mask))

          ! deallocate arrays
          deallocate(data_dp_2d, mask_2d)
       end if
#endif

        ! put data into global L0 variable
        select case (iVar)
        case(1) ! flow accumulation
          call append(L0_fAcc, pack(data_i4_2d, level0_iBasin%mask))
        case(2) ! flow direction
          ! rotate flow direction and any other variable with directions
          ! NOTE: ONLY when ASCII files are read
          call rotate_fdir_variable(data_i4_2d)
          ! append
          call append(L0_fDir, pack(data_i4_2d, level0_iBasin%mask))
        case(3) ! location of evaluation and inflow gauging stations
          ! evaluation gauges
          ! Input data check
          do iGauge = 1, basin_mrm(iBasin)%nGauges
            ! If gaugeId is found in gauging location file?
            if (.not. any(data_i4_2d .EQ. basin_mrm(iBasin)%gaugeIdList(iGauge))) then
              call message()
              call message('***ERROR: Gauge ID "', trim(adjustl(num2str(basin_mrm(iBasin)%gaugeIdList(iGauge)))), &
                      '" not found in ')
              call message('          Gauge location input file: ', &
                      trim(adjustl(dirMorpho(iBasin))) // trim(adjustl(file_gaugeloc)))
              stop
            end if
          end do

          call append(L0_gaugeLoc, pack(data_i4_2d, level0_iBasin%mask))

          ! inflow gauges
          ! if no inflow gauge for this subbasin exists still matirx with dim of subbasin has to be paded
          if (basin_mrm(iBasin)%nInflowGauges .GT. 0_i4) then
            ! Input data check
            do iGauge = 1, basin_mrm(iBasin)%nInflowGauges
              ! If InflowGaugeId is found in gauging location file?
              if (.not. any(data_i4_2d .EQ. basin_mrm(iBasin)%InflowGaugeIdList(iGauge))) then
                call message()
                call message('***ERROR: Inflow Gauge ID "', &
                        trim(adjustl(num2str(basin_mrm(iBasin)%InflowGaugeIdList(iGauge)))), &
                        '" not found in ')
                call message('          Gauge location input file: ', &
                        trim(adjustl(dirMorpho(iBasin))) // trim(adjustl(file_gaugeloc)))
                stop 1
              end if
            end do
          end if

          call append(L0_InflowGaugeLoc, pack(data_i4_2d, level0_iBasin%mask))

        end select
        !
        ! deallocate arrays
        if (allocated(data_i4_2d)) deallocate(data_i4_2d)
        if (allocated(mask_2d)) deallocate(mask_2d)
        !
      end do
    end do

  end subroutine mrm_read_L0_data
  ! ---------------------------------------------------------------------------

  !    NAME
  !        mrm_read_discharge

  !    PURPOSE
  !>       \brief Read discharge timeseries from file

  !>       \details Read Observed discharge at the outlet of a catchment
  !>       and at the inflow of a catchment. Allocate global runoff
  !>       variable that contains the simulated runoff after the simulation.

  !    HISTORY
  !>       \authors Matthias Zink & Stephan Thober

  !>       \date Aug 2015

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine mrm_read_discharge

    use mo_append, only : paste
    use mo_common_constants, only : nodata_dp
    use mo_common_mHM_mRM_variables, only : evalPer, nTstepDay, opti_function, optimize, simPer
    use mo_common_variables, only : nBasins
    use mo_message, only : message
    use mo_mrm_file, only : udischarge
    use mo_mrm_global_variables, only : InflowGauge, gauge, mRM_runoff, nGaugesTotal, &
                                        nInflowGaugesTotal, nMeasPerDay
    use mo_read_timeseries, only : read_timeseries
    use mo_string_utils, only : num2str

    implicit none

    integer(i4) :: iGauge

    integer(i4) :: iBasin

    integer(i4) :: maxTimeSteps

    ! file name of file to read
    character(256) :: fName

    integer(i4), dimension(3) :: start_tmp, end_tmp

    real(dp), dimension(:), allocatable :: data_dp_1d

    logical, dimension(:), allocatable :: mask_1d


    !----------------------------------------------------------
    ! INITIALIZE RUNOFF
    !----------------------------------------------------------
    maxTimeSteps = maxval(simPer(1 : nBasins)%julEnd - simPer(1 : nBasins)%julStart + 1) * nTstepDay
    allocate(mRM_runoff(maxTimeSteps, nGaugesTotal))
    mRM_runoff = nodata_dp

    ! READ GAUGE DATA
    do iGauge = 1, nGaugesTotal
      ! get basin id
      iBasin = gauge%basinId(iGauge)
      ! get start and end dates
      start_tmp = (/evalPer(iBasin)%yStart, evalPer(iBasin)%mStart, evalPer(iBasin)%dStart/)
      end_tmp = (/evalPer(iBasin)%yEnd, evalPer(iBasin)%mEnd, evalPer(iBasin)%dEnd  /)
      ! evaluation gauge
      fName = trim(adjustl(gauge%fname(iGauge)))
      call read_timeseries(trim(fName), udischarge, &
              start_tmp, end_tmp, optimize, opti_function, &
              data_dp_1d, mask = mask_1d, nMeasPerDay = nMeasPerDay)
      data_dp_1d = merge(data_dp_1d, nodata_dp, mask_1d)
      call paste(gauge%Q, data_dp_1d, nodata_dp)
      deallocate (data_dp_1d)
    end do
    !
    ! inflow gauge
    !
    ! in mhm call InflowGauge%Q has to be initialized -- dummy allocation with period of basin 1 and initialization
    if (nInflowGaugesTotal .EQ. 0) then
      allocate(data_dp_1d(maxval(simPer(:)%julEnd - simPer(:)%julStart + 1)))
      data_dp_1d = nodata_dp
      call paste(InflowGauge%Q, data_dp_1d, nodata_dp)
    else
      do iGauge = 1, nInflowGaugesTotal
        ! get basin id
        iBasin = InflowGauge%basinId(iGauge)
        ! get start and end dates
        start_tmp = (/simPer(iBasin)%yStart, simPer(iBasin)%mStart, simPer(iBasin)%dStart/)
        end_tmp = (/simPer(iBasin)%yEnd, simPer(iBasin)%mEnd, simPer(iBasin)%dEnd  /)
        ! inflow gauge
        fName = trim(adjustl(InflowGauge%fname(iGauge)))
        call read_timeseries(trim(fName), udischarge, &
                start_tmp, end_tmp, optimize, opti_function, &
                data_dp_1d, mask = mask_1d, nMeasPerDay = nMeasPerDay)
        if (.NOT. (all(mask_1d))) then
          call message()
          call message('***ERROR: Nodata values in inflow gauge time series. File: ', trim(fName))
          call message('          During simulation period from ', num2str(simPer(iBasin)%yStart) &
                  , ' to ', num2str(simPer(iBasin)%yEnd))
          stop
        end if
        data_dp_1d = merge(data_dp_1d, nodata_dp, mask_1d)
        call paste(InflowGauge%Q, data_dp_1d, nodata_dp)
        deallocate (data_dp_1d)
      end do
    end if

  end subroutine mrm_read_discharge

  ! ---------------------------------------------------------------------------

  !    NAME
  !        mrm_read_total_runoff

  !    PURPOSE
  !>       \brief read simulated runoff that is to be routed

  !>       \details read spatio-temporal field of total runoff that has been
  !>       simulated by a hydrologic model or land surface model. This
  !>       total runoff will then be aggregated to the level 11 resolution
  !>       and then routed through the stream network.

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iBasin" basin id

  !    HISTORY
  !>       \authors Stephan Thober

  !>       \date Sep 2015

  ! Modifications:
  ! Stephan Thober  Feb 2016 - refactored deallocate statements
  ! Stephan Thober  Sep 2016 - added ALMA convention
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine mrm_read_total_runoff(iBasin)

    use mo_append, only : append
    use mo_common_constants, only : HourSecs, nodata_dp
    use mo_common_mHM_mRM_variables, only : simPer, timestep
    use mo_common_variables, only : ALMA_convention, level1
    use mo_mrm_global_variables, only : L1_total_runoff_in, dirTotalRunoff, filenameTotalRunoff, &
                                        varnameTotalRunoff
    use mo_read_forcing_nc, only : read_forcing_nc

    implicit none

    ! basin id
    integer(i4), intent(in) :: iBasin

    integer(i4) :: tt

    integer(i4) :: nTimeSteps

    ! tell nc file to read daily or hourly values
    integer(i4) :: nctimestep

    ! read data from file
    real(dp), dimension(:, :, :), allocatable :: L1_data

    real(dp), dimension(:, :), allocatable :: L1_data_packed


    if (timestep .eq. 1) nctimestep = -4 ! hourly input
    if (timestep .eq. 24) nctimestep = -1 ! daily input
    call read_forcing_nc(trim(dirTotalRunoff(iBasin)), level1(iBasin)%nrows, level1(iBasin)%ncols, &
            varnameTotalRunoff, level1(iBasin)%mask, L1_data, target_period = simPer(iBasin), &
            nctimestep = nctimestep, filename = filenameTotalRunoff)
    ! pack variables
    nTimeSteps = size(L1_data, 3)
    allocate(L1_data_packed(level1(iBasin)%nCells, nTimeSteps))
    do tt = 1, nTimeSteps
      L1_data_packed(:, tt) = pack(L1_data(:, :, tt), mask = level1(iBasin)%mask)
    end do
    ! free space immediately
    deallocate(L1_data)

    ! convert if ALMA conventions have been given
    if (ALMA_convention) then
      ! convert from kg m-2 s-1 to mm TST-1
      ! 1 kg m-2 -> 1 mm depth
      ! multiply with time to obtain per timestep
      L1_data_packed = L1_data_packed * timestep * HourSecs
    end if

    ! append
    call append(L1_total_runoff_in, L1_data_packed(:, :), nodata_dp)

    !free space
    deallocate(L1_data_packed)

  end subroutine mrm_read_total_runoff

  subroutine mrm_read_bankfull_runoff(iBasin)
  ! ---------------------------------------------------------------------------

  !      NAME
  !          mrm_read_bankfull_runoff

  !>         \brief reads the bankfull runoff for approximating the channel widths

  !>         \details reads the bankfull runoff, which can be calculated with
  !>         the script in mhm/post_proc/bankfull_discharge.py
  
  !     INTENT(IN)
  !>        \param[in] "integer(i4)               :: iBasin"  basin id

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
  !>        \note The file read in must contain a double precision float variable with the name
  !>        "Q_bkfl".

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !         \author Lennart Schueler
  !         \date    May 2018

    use mo_mrm_global_variables, only: level11
    use mo_common_constants, only: HourSecs
    use mo_common_mHM_mRM_variables, only: timestep
    use mo_read_forcing_nc, only: read_const_forcing_nc
    use mo_mrm_global_variables, only: &
         dirBankfullRunoff, &   ! directory of bankfull_runoff file for each basin
         L11_bankfull_runoff_in ! bankfull runoff at L1
    use mo_common_variables, only: ALMA_convention

    implicit none

    ! input variables
    integer(i4), intent(in) :: iBasin

    ! local variables
    logical, dimension(:,:), allocatable :: mask
    real(dp), dimension(:,:), allocatable :: L11_data ! read data from file
    real(dp), dimension(:), allocatable :: L11_data_packed

    call read_const_forcing_nc(trim(dirBankfullRunoff(iBasin)), &
                               level11(iBasin)%nrows, &
                               level11(iBasin)%ncols, &
                               "Q_bkfl", mask, L11_data)

    allocate(L11_data_packed(level11(iBasin)%nCells))
    L11_data_packed(:) = pack(L11_data(:,:), mask=level11(iBasin)%mask)

    ! append
    if (allocated(L11_bankfull_runoff_in)) then
        L11_bankfull_runoff_in = [L11_bankfull_runoff_in, L11_data_packed]
    else
        allocate(L11_bankfull_runoff_in(size(L11_data_packed)))
        L11_bankfull_runoff_in = L11_data_packed
    end if

    deallocate(L11_data)
    deallocate(L11_data_packed)

  end subroutine mrm_read_bankfull_runoff

  !    NAME
  !        rotate_fdir_variable

  !    PURPOSE
  !>       \brief TODO: add description

  !>       \details TODO: add description

  !    INTENT(INOUT)
  !>       \param[inout] "integer(i4), dimension(:, :) :: x"

  !    HISTORY
  !>       \authors L. Samaniego & R. Kumar

  !>       \date Jun 2018

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine rotate_fdir_variable(x)

    use mo_common_constants, only : nodata_i4

    implicit none

    integer(i4), dimension(:, :), intent(INOUT) :: x


    ! 0  -1   0 |
    integer(i4) :: i, j


    !-------------------------------------------------------------------
    ! NOTE:
    !
    !     Since the DEM was transposed from (lat,lon) to (lon,lat), i.e.
    !     new DEM = transpose(DEM_old), then
    !
    !     the flow direction X (which was read) for every i, j cell needs
    !     to be rotated as follows
    !
    !                 X(i,j) = [R] * {uVector}
    !
    !     with
    !     {uVector} denoting the fDir_old (read) in vector form, and
    !               e.g. dir 8 => {-1, -1, 0 }
    !     [R]       denting a full rotation matrix to transform the flow
    !               direction into the new coordinate system (lon,lat).
    !
    !     [R] = [rx][rz]
    !
    !     with
    !           |      1       0      0 |
    !     [rx] =|      0   cos T  sin T | = elemental rotation along x axis
    !           |      0  -sin T  cos T |
    !
    !           |  cos F   sin F      0 |
    !     [rz] =| -sin F   cos F      0 | = elemental rotation along z axis
    !           |      0       0      1 |
    !
    !     and T = pi, F = - pi/2
    !     thus
    !          |  0  -1   0 |
    !     [R] =| -1   0   0 |
    !          |  0   0  -1 |
    !     making all 8 directions the following transformation were
    !     obtained.
    !-------------------------------------------------------------------
    do i = 1, size(x, 1)
      do j = 1, size(x, 2)
        if (x(i, j)  .eq. nodata_i4) cycle
        select case (x(i, j))
        case(1)
          x(i, j) = 4
        case(2)
          x(i, j) = 2
        case(4)
          x(i, j) = 1
        case(8)
          x(i, j) = 128
        case(16)
          x(i, j) = 64
        case(32)
          x(i, j) = 32
        case(64)
          x(i, j) = 16
        case(128)
          x(i, j) = 8
        end select
      end do
    end do

  end subroutine rotate_fdir_variable

end module mo_mrm_read_data

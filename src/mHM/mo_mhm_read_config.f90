!> \file mo_read_config.f90

!> \brief Reading of main model configurations.

!> \details This routine reads the configurations of mHM including, input and
!>          output directories, module usage specification, simulation time periods,
!>          global parameters, ...

!> \authors Matthias Zink
!> \date Dec 2012

MODULE mo_mhm_read_config

  USE mo_kind, ONLY : i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mhm_read_config ! read main directories

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         read_config

  !     PURPOSE
  !>        \brief Read main configurations for mHM

  !>        \details The main configurations in mHM are read from three files:
  !>                 <ol>
  !>                   <li> mhm.nml
  !>                   <li> mhm_parameters.nml
  !>                   <li> mhm_outputs.nml
  !>                 </ol>
  !>                 For details please refer to the above mentioned namelist files.

  !     CALLING SEQUENCE
  !         None

  !     INDENT(IN)
  !         None

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         None

  !     INDENT(IN), OPTIONAL
  !         None

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Matthias Zink
  !>        \date Dec 2012
  !         Modified Luis Samaniego,              Jan 2013 - messages
  !                  Rohini Kumar              
  !                  Matthias Cuntz,              Jan  2013 - namelist consolidation and positioning
  !                  Matthias Zink,               Jan  2013 - bug fix, added gaugeinfo reading
  !                  Rohini Kumar,                Jun  2013 - added restart flags
  !                  R. Kumar &            
  !                  S. Thober,                   Aug. 2013 - code change to incorporate output timestep
  !                                                           during writing of the netcdf file
  !                  Rohini Kumar,                Aug  2013 - name changed from "inputFormat" to inputFormat_meteo_forcings
  !                  Rohini Kumar,                Aug  2013 - added dirSoil_LUT and dirGeology_LUT, and changed in
  !                                                           namelist made accordingly
  !                  Rohini Kumar,                Aug  2013 - added new namelist for LAI related datasets, and changed in
  !                                                           within the code made accordingly
  !                  Matthias Zink,               Aug  2013 - changed read in for land cover period
  !                  Juliane Mai,                 Oct  2013 - adding global_parameters_name
  !                  Matthias Zink,               Nov  2013 - edited documentation and included DEFAULT cases for ptocess Matrix
  !                  Stephan Thober,              Nov  2013 - added read of directories where latitude longitude fields are located
  !                  Matthias Zink,               Feb  2014 - added multiple options for PET process
  !                  Matthias Zink,               Mar  2014 - added inflow from upstream areas and gauge information as namelist
  !                  Rohini Kumar,                May  2014 - added options for the model run coordinate system
  !                  Stephan Thober,              May  2014 - added switch for chunk read in
  !                  Stephan Thober,              Jun  2014 - added option for switching off mpr
  !                  Matthias Cuntz & Juliane Mai Nov  2014 - LAI input from daily, monthly or yearly files
  !                  Matthias Zink,               Dec  2014 - adopted inflow gauges to ignore headwater cells
  !                  Matthias Zink,               Mar  2015 - added optional soil moisture read in for calibration
  !                  Matthias Cuntz,              Jul  2015 - removed adjustl from trim(adjustl()) of Geoparams for PGI compatibilty
  !                  Stephan Thober,              Aug  2015 - added read_config_routing and read_routing_params from mRM
  !                  Oldrich Rakovec,             Oct  2015 - added reading of the basin average TWS data
  !                  Rohini Kumar,                Mar  2016 - options to handle different soil databases
  !                  Stephan Thober,              Nov  2016 - moved nProcesses and processMatrix to common variables
  !                  Rohini Kumar,                Dec  2016 - option to handle monthly mean gridded fields of LAI
  !                  M.Zink & M. Cuneyd Demirel   Mar  2017 - Added Jarvis soil water stress function at SM process(3)  
  !                  M.C. Demirel & Simon Stisen  Apr  2017 - Added FC dependency on root fraction coefficient (ET) at SM process(3)  
  !                  Robert Schweppe              Dec  2017 - switched from fractional julian day to integer


  subroutine mhm_read_config(file_namelist, unamelist)

    use mo_message, only : message
    use mo_string_utils, only : num2str
    use mo_nml, only : open_nml, close_nml, position_nml
    use mo_common_mHM_mRM_read_config, only : common_check_resolution
    use mo_common_constants, only : &
            maxNoBasins, & ! maximum number of allowed basins
            nodata_i4                                 ! nodata values
    use mo_mpr_constants, only : &
            maxNoSoilHorizons ! maximum number of allowed soil layers
    use mo_file, only : &
            file_defOutput, udefOutput ! file specifying which output to write
    use mo_global_variables, only : &
            outputFlxState, & ! definition which output to write
            dirPrecipitation, dirTemperature, & ! directory of meteo input
            dirReferenceET, & ! PET input path  if process 5 is 'PET is input' (case 0)
            dirMinTemperature, dirMaxTemperature, & ! PET input paths if process 5 is Hargreaves-Samani (case 1)
            dirNetRadiation, & ! PET input paths if process 5 is Priestely-Taylor (case 2)
            dirabsVapPressure, dirwindspeed, & ! PET input paths if process 5 is Penman-Monteith (case 3)
            inputFormat_meteo_forcings, & ! input format: nc only
            timestep_model_inputs, & ! read input frequency
            dirSoil_moisture, timeStep_sm_input, & ! directory and time stepping of soil moisture data
            nSoilHorizons_sm_input, & ! No. of mhm soil horizons equivalent to soil moisture input
            basin_avg_TWS_obs, & ! basin avg TWS data
            fileTWS, & ! directory with basin average tws data
            dirNeutrons, timeStep_neutrons_input, & ! directory where neutron data is located
            dirEvapotranspiration, timeStep_et_input, & ! directory and time stepping of evapotranspiration data
            timeStep_model_outputs, & ! timestep for writing model outputs
            evap_coeff, & ! pan evaporation
            read_meteo_weights, & ! flag for read meteo weights
            fday_prec, fnight_prec, fday_pet, & ! day-night fraction
            fnight_pet, fday_temp, fnight_temp ! day-night fraction

    use mo_mpr_global_variables, only : &
            nSoilHorizons_mHM
    use mo_common_mhm_mrm_variables, only : &
            optimize, opti_function
    use mo_common_variables, only : &
            nBasins, & ! number of basins
            processMatrix ! process configuration

    implicit none

    character(*), intent(in) :: file_namelist
    integer, intent(in) :: unamelist

    ! LOCAL variables
    integer(i4) :: iBasin

    integer(i4), dimension(maxNoBasins) :: time_step_model_inputs
    character(256), dimension(maxNoBasins) :: dir_Precipitation
    character(256), dimension(maxNoBasins) :: dir_Temperature
    character(256), dimension(maxNoBasins) :: dir_MinTemperature
    character(256), dimension(maxNoBasins) :: dir_MaxTemperature
    character(256), dimension(maxNoBasins) :: dir_NetRadiation
    character(256), dimension(maxNoBasins) :: dir_windspeed
    character(256), dimension(maxNoBasins) :: dir_absVapPressure
    character(256), dimension(maxNoBasins) :: dir_ReferenceET
    character(256), dimension(maxNoBasins) :: dir_soil_moisture      ! soil moisture input
    character(256), dimension(maxNoBasins) :: file_TWS               ! total water storage input file
    character(256), dimension(maxNoBasins) :: dir_neutrons           ! ground albedo neutron input
    character(256), dimension(maxNoBasins) :: dir_evapotranspiration ! ground albedo neutron input


    ! define namelists
    ! namelist directories
    namelist /directories_mHM/ inputFormat_meteo_forcings, &
            dir_Precipitation, dir_Temperature, dir_ReferenceET, dir_MinTemperature, &
            dir_MaxTemperature, dir_absVapPressure, dir_windspeed, &
            dir_NetRadiation, time_step_model_inputs
    ! optional data used for optimization
    namelist /optional_data/ &
            dir_soil_moisture, &
            nSoilHorizons_sm_input, &
            timeStep_sm_input, &
            file_TWS, &
            dir_neutrons, &
            dir_evapotranspiration, &
            timeStep_et_input
    ! namelist for pan evaporation
    namelist /panEvapo/evap_coeff
    ! namelist for night-day ratio of precipitation, referenceET and temperature
    namelist /nightDayRatio/read_meteo_weights, fnight_prec, fnight_pet, fnight_temp

    ! name list regarding output
    namelist /NLoutputResults/timeStep_model_outputs, outputFlxState

    !===============================================================
    !  Read namelist main directories
    !===============================================================
    call open_nml(file_namelist, unamelist, quiet = .true.)

    allocate(dirPrecipitation(nBasins))
    allocate(dirTemperature(nBasins))
    allocate(dirwindspeed(nBasins))
    allocate(dirabsVapPressure(nBasins))
    allocate(dirReferenceET(nBasins))
    allocate(dirMinTemperature(nBasins))
    allocate(dirMaxTemperature(nBasins))
    allocate(dirNetRadiation(nBasins))
    allocate(dirSoil_Moisture(nBasins))
    allocate(dirNeutrons(nBasins))
    allocate(fileTWS(nBasins))
    ! allocate time periods
    allocate(timestep_model_inputs(nBasins))

    !===============================================================
    !  Read namelist for mainpaths
    !===============================================================
    call position_nml('directories_mHM', unamelist)
    read(unamelist, nml = directories_mHM)

    dirPrecipitation = dir_Precipitation(1 : nBasins)
    dirTemperature = dir_Temperature(1 : nBasins)
    dirReferenceET = dir_ReferenceET(1 : nBasins)
    dirMinTemperature = dir_MinTemperature(1 : nBasins)
    dirMaxTemperature = dir_MaxTemperature(1 : nBasins)
    dirNetRadiation = dir_NetRadiation(1 : nBasins)
    dirwindspeed = dir_windspeed(1 : nBasins)
    dirabsVapPressure = dir_absVapPressure(1 : nBasins)
    timestep_model_inputs = time_step_model_inputs(1 : nBasins)

    ! consistency check for timestep_model_inputs
    if (any(timestep_model_inputs .ne. 0) .and. &
            .not. all(timestep_model_inputs .ne. 0)) then
      call message()
      call message('***ERROR: timestep_model_inputs either have to be all zero or all non-zero')
      stop
    end if
    ! check for optimzation and timestep_model_inputs options
    if (optimize .and. (any(timestep_model_inputs .ne. 0))) then
      call message()
      call message('***ERROR: optimize and chunk read is switched on! (set timestep_model_inputs to zero)')
      stop
    end if

    !===============================================================
    !  Read namelist of optional input data
    !===============================================================
    ! read optional optional data if necessary
    if (optimize) then
      select case (opti_function)
      case(10 : 13, 28)
        ! soil moisture
        call position_nml('optional_data', unamelist)
        read(unamelist, nml = optional_data)
        dirSoil_moisture = dir_Soil_moisture(1 : nBasins)
        if (nSoilHorizons_sm_input .GT. nSoilHorizons_mHM) then
          call message()
          call message('***ERROR: Number of soil horizons representative for input soil moisture exceeded')
          call message('          defined number of soil horizions: ', adjustl(trim(num2str(maxNoSoilHorizons))), '!')
          stop
        end if
      case(17)
        ! neutrons
        call position_nml('optional_data', unamelist)
        read(unamelist, nml = optional_data)
        dirNeutrons = dir_neutrons(1 : nBasins)
        timeStep_neutrons_input = -1 ! TODO: daily, hard-coded, to be flexibilized
      case(27, 29, 30)
        ! evapotranspiration
        call position_nml('optional_data', unamelist)
        read(unamelist, nml = optional_data)
        dirEvapotranspiration = dir_evapotranspiration(1 : nBasins)
      case(15)
        ! basin average TWS data
        call position_nml('optional_data', unamelist)
        read(unamelist, nml = optional_data)
        fileTWS = file_TWS (1 : nBasins)

        allocate(basin_avg_TWS_obs%basinId(nBasins)); basin_avg_TWS_obs%basinId = nodata_i4
        allocate(basin_avg_TWS_obs%fName  (nBasins)); basin_avg_TWS_obs%fName(:) = num2str(nodata_i4)

        do iBasin = 1, nBasins
          if (trim(fileTWS(iBasin)) .EQ. trim(num2str(nodata_i4))) then
            call message()
            call message('***ERROR: mhm.nml: Filename of evaluation TWS data ', &
                    ' for subbasin ', trim(adjustl(num2str(iBasin))), &
                    ' is not defined!')
            call message('          Error occured in namelist: evaluation_tws')
            stop 1
          end if

          basin_avg_TWS_obs%basinId(iBasin) = iBasin
          basin_avg_TWS_obs%fname(iBasin) = trim(file_TWS(iBasin))
        end do
      end select
    end if

    !===============================================================
    ! Read night-day ratios and pan evaporation
    !===============================================================
    ! Evap. coef. for free-water surfaces
    call position_nml('panEvapo', unamelist)
    read(unamelist, nml = panEvapo)
    ! namelist for night-day ratio of precipitation, referenceET and temperature
    call position_nml('nightDayRatio', unamelist)
    read(unamelist, nml = nightDayRatio)
    !
    fday_prec = 1.0_dp - fnight_prec
    fday_pet = 1.0_dp - fnight_pet
    fday_temp = -1.0_dp * fnight_temp

    call common_check_resolution(.true., .false.)

    call close_nml(unamelist)

    !===============================================================
    ! Read output specifications for mHM
    !===============================================================
    call open_nml(file_defOutput, udefOutput, quiet = .true.)
    outputFlxState = .FALSE.
    call position_nml('NLoutputResults', udefOutput)
    read(udefOutput, nml = NLoutputResults)
    call close_nml(udefOutput)

    call message('')
    call message('Following output will be written:')
    call message('  STATES:')
    if (outputFlxState(1)) then
      call message('    interceptional storage                      (L1_inter)     [mm]')
    end if
    if (outputFlxState(2)) then
      call message('    height of snowpack                          (L1_snowpack)  [mm]')
    end if
    if (outputFlxState(3)) then
      call message('    soil water content in the single layers     (L1_soilMoist) [mm]')
    end if
    if (outputFlxState(4)) then
      call message('    volumetric soil moisture in the single layers              [mm/mm]')
    end if
    if (outputFlxState(5)) then
      call message('    mean volum. soil moisture averaged over all soil layers    [mm/mm]')
    end if
    if (outputFlxState(6)) then
      call message('    waterdepth in reservoir of sealed areas     (L1_sealSTW)   [mm]')
    end if
    if (outputFlxState(7)) then
      call message('    waterdepth in reservoir of unsat. soil zone (L1_unsatSTW)  [mm]')
    end if
    if (outputFlxState(8)) then
      call message('    waterdepth in reservoir of sat. soil zone   (L1_satSTW)    [mm]')
    end if
    if (processMatrix(10, 1) .eq. 0) outputFlxState(18) = .false. ! suppress output if process is off
    if (outputFlxState(18)) then
      call message('    ground albedo neutrons                      (L1_neutrons)  [cph]')
    end if

    call message('  FLUXES:')
    if (outputFlxState(9)) then
      call message('    actual evapotranspiration aET      (L1_pet)                [mm/T]')
    end if
    if (outputFlxState(10)) then
      call message('    total discharge generated per cell (L1_total_runoff)       [mm/T]')
    end if
    if (outputFlxState(11)) then
      call message('    direct runoff generated per cell   (L1_runoffSeal)         [mm/T]')
    end if
    if (outputFlxState(12)) then
      call message('    fast interflow generated per cell  (L1_fastRunoff)         [mm/T]')
    end if
    if (outputFlxState(13)) then
      call message('    slow interflow generated per cell  (L1_slowRunoff)         [mm/T]')
    end if
    if (outputFlxState(14)) then
      call message('    baseflow generated per cell        (L1_baseflow)           [mm/T]')
    end if
    if (outputFlxState(15)) then
      call message('    groundwater recharge               (L1_percol)             [mm/T]')
    end if
    if (outputFlxState(16)) then
      call message('    infiltration                       (L1_infilSoil)          [mm/T]')
    end if
    call message('')
    call message('FINISHED reading config')

    ! warning message
    if (any(outputFlxState) .and. optimize) then
      call message('WARNING: FLUXES and STATES netCDF will be not written since optimization flag is TRUE ')
    end if

  end subroutine mhm_read_config

END MODULE mo_mhm_read_config

!>       \file mo_write_ascii.f90

!>       \brief Module to write ascii file output.

!>       \details Module to write ascii file output.
!>       Writing model output to ASCII should be the exception. Therefore, output is written usually as NetCDF
!>       and only:
!>       (1) The configuration file of mHM,
!>       (2) the final parameter set after optimization, and
!>       (3) the simulated vs. observed daily discharge
!>       is written in ASCII file format to allow for a quick assurance of proper model runs.

!>       \authors Christoph Schneider, Juliane Mai, Luis Samaniego

!>       \date May 2013

! Modifications:

MODULE mo_write_ascii

  ! This module is a template for the UFZ CHS mesoscale hydrologic model mHM.

  ! Written  Christoph Schneider, May 2013
  ! Modified, Juliane Mai,        May 2013 - module version and documentation
  ! Modified, Luis Samaniego,     Nov 2013 - improving all formats
  ! Modified, Luis Samaniego,     Mar 2014 - added inflow gauge information write out
  ! Modified, Stephan Thober,     Jun 2014 - bug fixed: in writing network properties
  ! Modified, Rohini Kumar,       Jun 2014 - bug fixed: writing of max and min value of discharge
  ! Modified, Stephan Thober,     Aug 2015 - moved write_daily_obs_sim_discharge to mRM

  USE mo_kind, ONLY : i4, dp
  IMPLICIT NONE

  PUBLIC :: write_configfile                   ! Writes configuration file
  PUBLIC :: write_optifile                     ! Write final OF and best parameter set
  PUBLIC :: write_optinamelist                 ! Write final OF and best parameter set in a namelist format
  ! ------------------------------------------------------------------

  !    NAME
  !        write_configfile

  !    PURPOSE
  !>       \brief This modules writes the results of the configuration into an ASCII-file
  !>       \details

  !>       \details TODO: add description

  !    HISTORY
  !>       \authors Christoph Schneider

  !>       \date May 2013

  ! Modifications:
  ! Juliane Mai    May 2013 - module version and documentation
  ! Stephan Thober Jun 2014 - bug fix in L11 config print out
  ! Stephan Thober Jun 2014 - updated read_restart
  ! Rohini, Luis   Jul 2015 - updated version, L1 level prints
  ! Stephan Thober Nov 2016 - moved processMatrix to common variables
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  PRIVATE

  ! ------------------------------------------------------------------

CONTAINS

  !    NAME
  !        write_configfile

  !    PURPOSE
  !>       \brief TODO: add description

  !>       \details TODO: add description

  !    HISTORY
  !>       \authors Robert Schweppe

  !>       \date Jun 2018

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting
  ! P Shrestha, S Thober Aug 2018 - resolved bug while printing River Network in
  !                                 cases with multiple outlets.

  Subroutine write_configfile

    use mo_common_file, only : file_config, uconfig
    use mo_common_mHM_mRM_variables, only : LCyearId, SimPer, evalPer, read_restart, timeStep, warmPer
    use mo_common_variables, only : L0_Basin, LC_year_end, &
                                    LC_year_start, LCfilename, dirConfigOut, dirLCover, dirMorpho, dirOut, dirRestartOut, &
                                    global_parameters, global_parameters_name, iFlag_cordinate_sys, level0, level1, &
                                    nBasins, nLCoverScene, resolutionHydrology, write_restart
    use mo_file, only : version
    use mo_global_variables, only : dirPrecipitation, dirReferenceET, &
                                    dirTemperature
    use mo_kind, only : i4
    use mo_message, only : message
    use mo_string_utils, only : num2str
#ifdef MRM2MHM
    use mo_common_constants, only : nodata_dp
    use mo_common_mHM_mRM_variables, only : resolutionRouting
    use mo_common_variables, only : processMatrix
    use mo_mrm_global_variables, only : InflowGauge, L11_fromN, L11_label, L11_length, L11_netPerm, L11_rOrder, &
                                        L11_slope, L11_toN, L1_L11_ID, dirGauges, gauge, level11, nGaugesTotal, &
                                        nInflowGaugesTotal, L11_nOutlets
#endif

    implicit none

    character(256) :: fName

    integer(i4) :: i, j, n

    integer(i4) :: err


    fName = trim(adjustl(dirConfigOut)) // trim(adjustl(file_config))
    call message()
    call message('  Log-file written to ', trim(fName))
    open(uconfig, file = fName, status = 'unknown', action = 'write', iostat = err)
    if (err .ne. 0) then
      call message('  Problems while creating File')
      call message('  Error-Code', num2str(err))
      stop
    end if
    write(uconfig, 200)
    write(uconfig, 100) 'mHM-UFZ v-' // trim(version)
    write(uconfig, 100) 'L. Samaniego & R. Kumar, UFZ'
    write(uconfig, 200)
    write(uconfig, 100)
    write(uconfig, 201) '         M A I N  mHM  C O N F I G U R A T I O N  I N F O R M A T I O N         '
    write(uconfig, 100)
    write(uconfig, 103) 'Number of basins            ', nBasins
#ifdef MRM2MHM
    if (processMatrix(8, 1) .ne. 0) then
      write(uconfig, 103) 'Total No. of gauges         ', nGaugesTotal
    end if
#endif
    write(uconfig, 103)    'Time Step [h]               ', timeStep
    do i = 1, nBasins
      write(uconfig, 103) 'Basin  ', i, 'No. of cells L0             ', level0(L0_Basin(i))%nCells
      write(uconfig, 103) 'Basin  ', i, 'No. of cells L1             ', level1(i)%nCells
#ifdef MRM2MHM
    if (processMatrix(8, 1) .ne. 0) then
      write(uconfig, 103) 'Total No. of nodes          ', level11(i)%nCells
      write(uconfig, 103) 'Total No. of reaches        ', level11(i)%nCells - 1
      if (processMatrix(8, 1) .ne. 0) then
        write(uconfig, 103) 'No. of cells L11            ', level11(i)%nCells
        write(uconfig, 103) 'Total No. of gauges         ', nGaugesTotal
      end if
    end if
#endif

      select case (iFlag_cordinate_sys)
      case (0)
        write(uconfig, 301)      'Basin  ', i, '   Hydrology Resolution [m]      ', resolutionHydrology(i)
#ifdef MRM2MHM
        if (processMatrix(8, 1) .ne. 0) then
          write(uconfig, 301)   'Basin  ', i, '   Routing Resolution [m]        ', resolutionRouting(i)
        end if
#endif
       case(1)
        write(uconfig, 302)       'Basin  ', i, '   Hydrology Resolution [o]      ', resolutionHydrology(i)
#ifdef MRM2MHM
        if (processMatrix(8, 1) .ne. 0) then
          write(uconfig, 302)   'Basin  ', i, '   Routing Resolution [o]        ', resolutionRouting(i)
        end if
#endif

      end select
    end do
    write(uconfig, 126)    'Flag READ  restart            ', read_restart
    write(uconfig, 126)    'Flag WRITE restart            ', write_restart
    !
    !******************
    ! Model Run period
    !******************
    do j = 1, nBasins
      write(uconfig, 115) '                      Model Run Periods for Basin ', num2str(j)
      write(uconfig, 116) &
              'From                To', &
              '   Day Month  Year   Day Month  Year'
      write(uconfig, 117)  &
              'Warming Period (1)            ', &
              warmPer(j)%dStart, warmPer(j)%mStart, warmPer(j)%yStart, &
              warmPer(j)%dEnd, warmPer(j)%mEnd, warmPer(j)%yEnd, &
              'Evaluation Period (2)         ', &
              evalPer(j)%dStart, evalPer(j)%mStart, evalPer(j)%yStart, &
              evalPer(j)%dEnd, evalPer(j)%mEnd, evalPer(j)%yEnd, &
              'Simulation Period (1)+(2)     ', &
              SimPer(j)%dStart, SimPer(j)%mStart, SimPer(j)%yStart, &
              SimPer(j)%dEnd, SimPer(j)%mEnd, SimPer(j)%yEnd
    end do

    !*********************************
    ! Model Land Cover Observations
    !*********************************
    do j = 1, nBasins
      write(uconfig, 118) '       Land Cover Observations for Basin ', num2str(i)
      write(uconfig, 119) ' Start Year', ' End Year', '    Land cover scene', 'Land Cover File'
      do i = 1, nLCoverScene
        write(uconfig, 120) LC_year_start(i), LC_year_end(i), &
                LCyearId(max(evalPer(j)%yStart, LC_year_start(i)), j), trim(LCfilename(i))
      end do
    end do
    !*********************************
    ! Initial Parameter Ranges
    !*********************************
    write(uconfig, 121) '  Initial Transfer Function Parameter Ranges (gammas)  '
    !
    ! Transfer functions
    write(uconfig, 122)      &
            '         i', '            min', '            max', '        current', &
            '                               name'
    do i = 1, size(global_parameters, 1)
      write(uconfig, 123) &
              i, global_parameters(i, 1), global_parameters(i, 2), global_parameters(i, 3), &
              trim(adjustl(global_parameters_name(i)))
    end do
#ifdef MRM2MHM
    ! basin runoff data
    if (processMatrix(8, 1) .ne. 0) then
      write(uconfig, 202) '                Basin Runoff Data                '
      write(uconfig, 107) ' Gauge No.', '  Basin Id', '     Qmax[m3/s]', '     Qmin[m3/s]'
      do i = 1, nGaugesTotal
        if(any(gauge%Q(:, i) > nodata_dp)) then
          write(uconfig, 108) i, gauge%basinId(i), maxval(gauge%Q(:, i), gauge%Q(:, i) > nodata_dp), &
                  minval(gauge%Q(:, i), gauge%Q(:, i) > nodata_dp)
        else
          write(uconfig, 108) i, gauge%basinId(i), nodata_dp, nodata_dp
        end if
      end do
    end if
    ! inflow gauge data
    if (nInflowGaugesTotal .GT. 0) then
      write(uconfig, 202) '                Basin Inflow Data                 '
      write(uconfig, 107) ' Gauge No.', '  Basin Id', '     Qmax[m3/s]', '     Qmin[m3/s]'
      do i = 1, nInflowGaugesTotal
        if(all(InflowGauge%Q(:, i) > nodata_dp)) then
          write(uconfig, 108) i, InflowGauge%basinId(i), maxval(InflowGauge%Q(:, i), InflowGauge%Q(:, i) > nodata_dp), &
                  minval(InflowGauge%Q(:, i), InflowGauge%Q(:, i) > nodata_dp)
        else
          write(uconfig, 108) i, InflowGauge%basinId(i), nodata_dp, nodata_dp
        end if
      end do
    end if
#endif
    ! basin config
    write(uconfig, 218) 'Basin-wise Configuration'
    do n = 1, nBasins
      !ST has to be moved to the config write of mRM
      ! if ( processMatrix(8,1) .ne. 0 ) then
      !    write(uconfig,103) 'Basin No.                   ', n, &
      !         'No. of gauges               ', basin%nGauges(n)
      ! end if

      write(uconfig, 222)   'Directory list'

      write(uconfig, 224) 'Directory to morphological input         ', dirMorpho(n)
      write(uconfig, 224) 'Directory to land cover input            ', dirLCover(n)
#ifdef MRM2MHM
      if (processMatrix(8, 1) .ne. 0) then
        write(uconfig, 224) 'Directory to gauging station input       ', dirGauges(n)
      end if
#endif
      write(uconfig, 224) 'Directory to precipitation input         ', dirPrecipitation(n)
      write(uconfig, 224) 'Directory to temperature input           ', dirTemperature(n)
      write(uconfig, 224) 'Directory to reference ET input          ', dirReferenceET(n)
      write(uconfig, 224) 'Directory to write output by default     ', dirOut(n)
      write(uconfig, 224) 'Directory to write output when restarted ', dirRestartOut(n)

#ifdef MRM2MHM
      if (processMatrix(8, 1) .ne. 0) then
        write(uconfig, 102) 'River Network  (Routing level)'
        write(uconfig, 100) 'Label 0 = intermediate draining cell '
        write(uconfig, 100) 'Label 1 = headwater cell             '
        write(uconfig, 100) 'Label 2 = sink cell                  '

        if (processMatrix(8, 1) .eq. 1_i4) then
          write(uconfig, 104) '   Overall', &
                  '      From', &
                  '        To', &
                  '   Routing', &
                  '     Label', &
                  '    Length', &
                  '      Mean', &
                  '      Link', &
                  '   Routing', &
                  '   Routing', &
                  '  Sequence', &
                  '          ', &
                  '          ', &
                  '     Slope'
          !
          write(uconfig, 105) '        Id', &
                  '      Node', &
                  '      Node', &
                  '', &
                  '', &
                  '      [km]', &
                  '    [o/oo]'
          !
          do j = level11(n)%iStart, level11(n)%iEnd -  L11_nOutlets(n)
            i = L11_netPerm(j) + level11(n)%iStart - 1 ! adjust permutation for multi-basin option
            write(uconfig, 106) i, L11_fromN(i), L11_toN(i), L11_rOrder(i), L11_label(i), &
                    L11_length(i) / 1000.0_dp, L11_slope(i) * 1.0e3_dp
          end do

        else if (processMatrix(8, 1) .eq. 2_i4) then
          write(uconfig, 134) '   Overall', &
                  '      From', &
                  '        To', &
                  '   Routing', &
                  '     Label', &
                  '      Link', &
                  '   Routing', &
                  '   Routing', &
                  '  Sequence', &
                  '          '
          !
          write(uconfig, 135) '        Id', &
                  '      Node', &
                  '      Node', &
                  '', &
                  ''
          !
          do j = level11(n)%iStart, level11(n)%iEnd -  L11_nOutlets(n)
            i = L11_netPerm(j) + level11(n)%iStart - 1 ! adjust permutation for multi-basin option
            write(uconfig, 136) i, L11_fromN(i), L11_toN(i), L11_rOrder(i), L11_label(i)
          end do
        end if
        ! draining node at L11
        write(uconfig, 109)  '   Overall', '     Basin', &
                '      Cell', '   Routing', &
                '        Id', '   Node Id'
        do i = 1, level11(n)%nCells
          write(uconfig, 110) i
        end do

        ! L1 level information
        write(uconfig, 111)  '  Modeling', '   Routing', ' Effective', &
                '      Cell', '   Cell Id', '      Area', &
                '        Id', '       [-]', '     [km2]'

        do i = 1, level1(n)%nCells
          write(uconfig, 113) i, L1_L11_Id(i), level1(n)%CellArea(i) *  1.0E-6_dp
        end do
        write(uconfig, 114)  ' Total[km2]', sum(level1(n)%CellArea) *  1.0E-6_dp
      end if
#endif
       !
    end do

    write(uconfig, *)
    close(uconfig)

    !! Formats
    100 format (a80)
    102 format (/ 30('-') / a30 / 30('-'))
    103 format (a20, 10x, i10)
    104 format (/ 75('-') / 5a10, 5x, 2a10 / 5a10, 5x, 2a10)
    105 format (5a10, 5x, 2a10 / 75('-'))
    106 format (5i10, 5x, 2f10.3)
    107 format (2a10, 2a15)
    108 format (2i10, 2f15.3)
    !
    109 format (/ 20('-') / 2a10 / 2a10 / 2a10 / 20('-'))
    110 format (i10)
    !
    111 format (/ 30('-') / 3a10 / 3a10 / 3a10 /  30('-'))
    113 format (2i10, 1f10.3)
    114 format (30('-') / a15, 5x, 1f10.3 /)
    !
    115 format (/61('-')/ a50, a10 /61('-'))
    116 format (39x, a22 / 25x, a36)
    117 format (3(a25, 6(i6)))
    !
    118 format (/50('-')/ a40, a10  /50('-'))
    119 format (a10, a10, a20, a20/)
    120 format (i10, i10, 10x, i10, a20)
    !
    121 format (/55('-')/ a55 /55('-'))
    122 format (a10, 3a15, a35)
    123 format (i10, 3f15.3, a35)
    !
    126 format (a30, 9x, L1)
    !
    134 format (/ 50('-') / 5a10 / 5a10)
    135 format (5a10 / 50('-'))
    136 format (5i10)
    !
    200 format (80('-'))
    201 format (a80)
    202 format (/50('-')/ a50 /50('-'))
    !
    218 format (/ 80('-')/ 26x, a24, 26x, /80('-'))
    222 format (/80('-')/ 26x, a21 /80('-'))
    224 format (a40, 5x, a256)

    301 format (a7, i2, a32, f15.0)
    302 format (a7, i2, a32, es20.8)
  end Subroutine write_configfile


  ! ------------------------------------------------------------------

  !    NAME
  !        write_optifile

  !    PURPOSE
  !>       \brief Write briefly final optimization results.

  !>       \details Write overall best objective function and the best optimized parameter set to a file_opti.

  !    INTENT(IN)
  !>       \param[in] "real(dp) :: best_OF"                             best objective function value as returnedby the
  !>       optimization routine
  !>       \param[in] "real(dp), dimension(:) :: best_paramSet"         best associated global parameter setCalled only
  !>       when optimize is .TRUE.
  !>       \param[in] "character(len = *), dimension(:) :: param_names"

  !    HISTORY
  !>       \authors David Schaefer

  !>       \date July 2013

  ! Modifications:
  ! Rohini Kumar Aug 2013 - change in structure of the code including call statements
  ! Juliane Mai  Oct 2013 - clear parameter names added 
  !                       - double precision written
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine write_optifile(best_OF, best_paramSet, param_names)

    use mo_common_mhm_mrm_file, only : file_opti, uopti
    use mo_common_variables, only : dirConfigOut
    use mo_message, only : message
    use mo_string_utils, only : num2str

    implicit none

    ! best objective function value as returnedby the optimization routine
    real(dp), intent(in) :: best_OF

    ! best associated global parameter setCalled only when optimize is .TRUE.
    real(dp), dimension(:), intent(in) :: best_paramSet

    character(len = *), dimension(:), intent(in) :: param_names

    character(256) :: fName, formHeader, formParams

    integer(i4) :: ii, err, n_params


    ! number of parameters
    n_params = size(best_paramSet)

    ! open file
    fName = trim(adjustl(dirConfigOut)) // trim(adjustl(file_opti))
    open(uopti, file = fName, status = 'unknown', action = 'write', iostat = err, recl = (n_params + 1) * 40)
    if(err .ne. 0) then
      call message ('  IOError while openening ', trim(fName))
      call message ('  Error-Code ', num2str(err))
      stop
    end if

    ! header
    write(formHeader, *) '(a40,', n_params, 'a40)'
    ! len(param_names(1))=256 but only 39 characters taken here
    ! write(uopti, formHeader) 'OF', (trim(adjustl(param_names(ii))), ii=1, n_params)
    write(uopti, formHeader) 'OF', (trim(adjustl(param_names(ii)(1 : 39))), ii = 1, n_params)

    ! output
    write(formParams, *) '( es40.14, ', n_params, '(es40.14) )'
    write(uopti, formParams) best_OF, (best_paramSet(ii), ii = 1, n_params)

    ! close file
    close(uopti)

    ! screen output
    call message()
    call message(' Optimized parameters written to ', trim(fName))

  end subroutine write_optifile

  ! ------------------------------------------------------------------

  !    NAME
  !        write_optinamelist

  !    PURPOSE
  !>       \brief Write final, optimized parameter set in a namelist format.

  !>       \details Write final, optimized parameter set in a namelist format.
  !>       Only parameters of processes which were switched on are written to the namelist.
  !>       All others are discarded.

  !    INTENT(IN)
  !>       \param[in] "integer(i4), dimension(nProcesses, 3) :: processMatrix"                information about which
  !>       process
  !>       case was used
  !>       \param[in] "real(dp), dimension(:, :) :: parameters"                               (min, max, opti)
  !>       \param[in] "logical, dimension(size(parameters, 1)) :: maskpara"                   .true. if parameter was
  !>       calibrated
  !>       \param[in] "character(len = *), dimension(size(parameters, 1)) :: parameters_name" clear names of parameters

  !    HISTORY
  !>       \authors Juliane Mai

  !>       \date Dec 2013

  ! Modifications:
  ! Stephan Thober Nov  2016 - moved nProcesses to common variables
  ! Stephan Thober Nov  2016 - write namelist for routing process 2
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine write_optinamelist(processMatrix, parameters, maskpara, parameters_name)

    use mo_common_mhm_mrm_file, only : file_opti_nml, uopti_nml
    use mo_common_variables, only : dirConfigOut, nProcesses
    use mo_message, only : message
    use mo_string_utils, only : num2str

    implicit none

    ! information about which process
    ! case was used
    integer(i4), dimension(nProcesses, 3), intent(in) :: processMatrix

    ! (min, max, opti)
    real(dp), dimension(:, :), intent(in) :: parameters

    ! .true. if parameter was calibrated
    logical, dimension(size(parameters, 1)), intent(in) :: maskpara

    ! clear names of parameters
    character(len = *), dimension(size(parameters, 1)), intent(in) :: parameters_name

    character(256) :: fName

    character(3) :: flag

    character(len = 28), dimension(nProcesses) :: Process_descr

    integer(i4) :: err

    integer(i4) :: iProc, iPar, iPar_start


    Process_descr(1) = 'interception'
    Process_descr(2) = 'snow'
    Process_descr(3) = 'soilmoisture'
    Process_descr(4) = 'directSealedAreaRunoff'
    Process_descr(5) = 'potential evapotranspiration'
    Process_descr(6) = 'interflow'
    Process_descr(7) = 'percolation'
    Process_descr(8) = 'routing'
    Process_descr(9) = 'geology'
    Process_descr(10) = 'neutrons'

    ! open file
    fName = trim(adjustl(dirConfigOut)) // trim(adjustl(file_opti_nml))
    open(uopti_nml, file = fName, status = 'unknown', action = 'write', iostat = err)
    if(err .ne. 0) then
      call message ('  IOError while openening ', trim(fName))
      call message ('  Error-Code ', num2str(err))
      stop
    end if

    write(uopti_nml, *) '!global_parameters'
    write(uopti_nml, *) '!PARAMETER                       lower_bound  upper_bound          value   FLAG  SCALING'

    iPar_start = 1
    do iProc = 1, nProcesses

      write(uopti_nml, *) '! ', trim(adjustl(process_descr(iProc)))

      select case (iProc)
      case(1)
        if (processMatrix(iProc, 1) .eq. 1) then
          write(uopti_nml, *) '&interception1'
        end if
      case(2)
        if (processMatrix(iProc, 1) .eq. 1) then
          write(uopti_nml, *) '&snow1'
        end if
      case(3)
        select case (processMatrix(iProc, 1))
        case(1)
          write(uopti_nml, *) '&soilmoisture1'
        case(2)
          write(uopti_nml, *) '&soilmoisture2'
        end select
      case(4)
        if (processMatrix(iProc, 1) .eq. 1) then
          write(uopti_nml, *) '&directRunoff1'
        end if
      case(5)
        select case (processMatrix(iProc, 1))
        case(0)
          write(uopti_nml, *) '&PET0'
        case(1)
          write(uopti_nml, *) '&PET1'
        case(2)
          write(uopti_nml, *) '&PET2'
        case(3)
          write(uopti_nml, *) '&PET3'
        end select
      case(6)
        if (processMatrix(iProc, 1) .eq. 1) then
          write(uopti_nml, *) '&interflow1'
        end if
      case(7)
        if (processMatrix(iProc, 1) .eq. 1) then
          write(uopti_nml, *) '&percolation1'
        end if
      case(8)
        if (processMatrix(iProc, 1) .eq. 1) then
          write(uopti_nml, *) '&routing1'
        end if
        if (processMatrix(iProc, 1) .eq. 2) then
          write(uopti_nml, *) '&routing2'
        end if
      case(9)
        if (processMatrix(iProc, 1) .eq. 1) then
          write(uopti_nml, *) '&geoparameter'
        end if
      case(10)
        if (processMatrix(iProc, 1) .ge. 1) then
          write(uopti_nml, *) '&neutrons1'
        end if
      end select

      do iPar = iPar_Start, processMatrix(iProc, 3)

        if (maskpara(iPar)) then
          flag = ' 1 '
        else
          flag = ' 0 '
        end if

        write(uopti_nml, *) trim(adjustl(parameters_name(iPar))), ' = ', &
                parameters(iPar, 1), ' , ', &
                parameters(iPar, 2), ' , ', &
                parameters(iPar, 3), ' , ', &
                flag, ', 1 '
      end do

      iPar_Start = processMatrix(iProc, 3) + 1

      write(uopti_nml, *) '/'
      write(uopti_nml, *) ' '

    end do ! loop over processes

    ! close file
    close(uopti_nml)

    ! screen output
    call message()
    call message(' Optimized parameters written in namelist format to ', trim(fName))

  end subroutine write_optinamelist

END MODULE mo_write_ascii

!>       \file mo_objective_function.f90

!>       \brief Objective Functions for Optimization of mHM.

!>       \details This module provides a wrapper for several objective functions used to optimize mHM against various
!>       variables.
!>       If the objective is only regarding runoff move it to mRM/mo_mrm_objective_function_runoff.f90.
!>       If it contains besides runoff another variable like TWS implement it here.

!>       All the objective functions are supposed to be minimized!
!>       (10) SO: SM:       1.0 - KGE of catchment average soilmoisture
!>       (11) SO: SM:       1.0 - Pattern dissimilarity (PD) of spatially distributed soil moisture
!>       (12) SO: SM:       Sum of squared errors (SSE) of spatially distributed standard score (normalization)
!>       of soil moisture
!>       (13) SO: SM:       1.0 - average temporal correlation of spatially distributed soil moisture
!>       (15) SO: Q + TWS:  [1.0-KGE(Q)]*RMSE(basin_avg_TWS) - objective function using Q and basin average
!>       (standard score) TWS
!>       (17) SO: N:        1.0 - KGE of spatio-temporal neutron data, catchment-average
!>       (27) SO: ET:       1.0 - KGE of catchment average evapotranspiration

!>       \authors Juliane Mai

!>       \date Dec 2012

! Modifications:
! Oldrich Rakovec Oct 2015 - added obj. func. 15 (objective_kge_q_rmse_tws) and extract_basin_avg_tws routine
! Robert Schweppe Jun 2018 - refactoring and reformatting

MODULE mo_objective_function

  ! This module provides objective functions for optimization of the UFZ CHS mesoscale hydrologic model mHM.

  ! Written  Juliane Mai, Dec 2012
  ! Modified Stephan Thober, Oct 2015 moved all runoff only related objectives to mRM

  USE mo_kind, ONLY : i4, dp
  use mo_optimization_utils, only : eval_interface

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: objective ! objective function wrapper for soil moisture only

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        objective

  !    PURPOSE
  !>       \brief Wrapper for objective functions.

  !>       \details The functions selects the objective function case defined in a namelist,
  !>       i.e. the global variable \e opti\_function.
  !>       It return the objective function value for a specific parameter set.

  !    INTENT(IN)
  !>       \param[in] "REAL(dp), DIMENSION(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    INTENT(IN), OPTIONAL
  !>       \param[in] "real(dp), optional :: arg1"

  !    INTENT(OUT), OPTIONAL
  !>       \param[out] "real(dp), optional :: arg2"
  !>       \param[out] "real(dp), optional :: arg3"

  !    RETURN
  !>       \return real(dp) :: objective &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !    HISTORY
  !>       \authors Juliane Mai

  !>       \date Dec 2012

  ! Modifications:
  ! Stephan Thober Oct 2015 - moved all runoff related objective functions to mRM
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION objective(parameterset, eval, arg1, arg2, arg3)

    use mo_common_constants, only : nodata_dp
    use mo_common_mHM_mRM_variables, only : opti_function
    use mo_message, only : message

    implicit none

    REAL(dp), DIMENSION(:), INTENT(IN) :: parameterset

    procedure(eval_interface), INTENT(IN), POINTER :: eval

    real(dp), optional, intent(in) :: arg1

    real(dp), optional, intent(out) :: arg2

    real(dp), optional, intent(out) :: arg3

    REAL(dp) :: objective


    if (present(arg1) .or. present(arg2) .or. present(arg3)) then
      call message("Error mo_objective_function: Received unexpected argument, check optimization settings")
      stop 1
    end if

    ! set these to nan so compiler does not complain
    if (present(arg2)) then
      arg2 = nodata_dp
    end if
    if (present(arg3)) then
      arg3 = nodata_dp
    end if

    select case (opti_function)
    case (10)
      ! KGE of catchment average SM
      objective = objective_sm_kge_catchment_avg(parameterset, eval)
    case (11)
      ! pattern dissimilarity (PD) of SM fields
      objective = objective_sm_pd(parameterset, eval)
    case (12)
      ! sum of squared errors of standard_score SM
      objective = objective_sm_sse_standard_score(parameterset, eval)
    case (13)
      ! soil moisture correlation - temporal
      objective = objective_sm_corr(parameterset, eval)
    case (15)
      ! KGE for Q * RMSE for basin_avg TWS (standarized scored)
      objective = objective_kge_q_rmse_tws(parameterset, eval)
    case (17)
      ! KGE of catchment average SM
      objective = objective_neutrons_kge_catchment_avg(parameterset, eval)
    case (27)
      ! KGE of catchment average ET
      objective = objective_et_kge_catchment_avg(parameterset, eval)
    case (28)
      !  KGE for Q + SSE for SM (standarized scored)
      objective = objective_kge_q_sm_corr(parameterset, eval)
    case (29)
      !  KGE for Q + KGE of catchment average ET
      objective = objective_kge_q_et(parameterset, eval)
    case (30)
      ! KGE for Q * RMSE for basin_avg ET (standarized scored)
      objective = objective_kge_q_rmse_et(parameterset, eval)

    case default
      call message("Error objective: opti_function not implemented yet.")
      stop 1
    end select

  END FUNCTION objective

  ! ------------------------------------------------------------------

  !    NAME
  !        objective_sm_kge_catchment_avg

  !    PURPOSE
  !>       \brief Objective function for soil moisture.

  !>       \details The objective function only depends on a parameter vector.
  !>       The model will be called with that parameter vector and
  !>       the model output is subsequently compared to observed data.

  !>       Therefore, the Kling-Gupta model efficiency \f$ KGE \f$ of the catchment average
  !>       soil mloisture (SM) is calculated
  !>       \f[ KGE = 1.0 - \sqrt{( (1-r)^2 + (1-\alpha)^2 + (1-\beta)^2 )} \f]
  !>       where
  !>       \f$ r \f$ = Pearson product-moment correlation coefficient
  !>       \f$ \alpha \f$ = ratio of simulated mean to observed mean SM
  !>       \f$ \beta  \f$ = ratio of similated standard deviation to observed standard deviation
  !>       is calculated and the objective function for a given basin \f$ i \f$ is
  !>       \f[ \phi_{i} = 1.0 - KGE_{i} \f]
  !>       \f$ \phi_{i} \f$ is the objective since we always apply minimization methods.
  !>       The minimal value of \f$ \phi_{i} \f$ is 0 for the optimal KGE of 1.0.

  !>       Finally, the overall objective function value \f$ OF \f$ is estimated based on the power-6
  !>       norm to combine the \f$ \phi_{i} \f$ from all basins \f$ N \f$.
  !>       \f[ OF = \sqrt[6]{\sum((1.0 - KGE_{i})/N)^6 }.  \f]
  !>       The observed data L1_sm, L1_sm_mask are global in this module.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp) :: objective_sm_kge_catchment_avg &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date May 2015

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION objective_sm_kge_catchment_avg(parameterset, eval)

    use mo_common_constants, only : nodata_dp
    use mo_common_variables, only : level1, nBasins
    use mo_errormeasures, only : KGE
    use mo_global_variables, only : L1_sm, L1_sm_mask
    use mo_message, only : message
    use mo_moment, only : average
    use mo_string_utils, only : num2str

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), POINTER :: eval

    real(dp) :: objective_sm_kge_catchment_avg

    ! basin loop counter
    integer(i4) :: iBasin

    ! time loop counter
    integer(i4) :: iTime

    ! number of time steps in simulated SM
    integer(i4) :: n_time_steps

    ! start and end index for the current basin
    integer(i4) :: s1, e1

    ! ncells1 of level 1
    integer(i4) :: ncells1

    ! number of invalid timesteps
    real(dp) :: invalid_times

    ! for sixth root
    real(dp), parameter :: onesixth = 1.0_dp / 6.0_dp

    ! spatial average of observed soil moisture
    real(dp), dimension(:), allocatable :: sm_catch_avg_basin

    ! spatial avergae of modeled  soil moisture
    real(dp), dimension(:), allocatable :: sm_opti_catch_avg_basin

    ! simulated soil moisture
    ! (dim1=ncells, dim2=time)
    real(dp), dimension(:, :), allocatable :: sm_opti

    ! mask for valid sm catchment avg time steps
    logical, dimension(:), allocatable :: mask_times


    call eval(parameterset, sm_opti = sm_opti)

    ! initialize some variables
    objective_sm_kge_catchment_avg = 0.0_dp

    ! loop over basin - for applying power law later on
    do iBasin = 1, nBasins

      ! get basin information
      ncells1 = level1(iBasin)%ncells
      s1 = level1(iBasin)%iStart
      e1 = level1(iBasin)%iEnd

      ! allocate
      allocate(mask_times             (size(sm_opti, dim = 2)))
      allocate(sm_catch_avg_basin     (size(sm_opti, dim = 2)))
      allocate(sm_opti_catch_avg_basin(size(sm_opti, dim = 2)))

      ! initalize
      mask_times = .TRUE.
      sm_catch_avg_basin = nodata_dp
      sm_opti_catch_avg_basin = nodata_dp

      invalid_times = 0.0_dp
      ! calculate catchment average soil moisture
      n_time_steps = size(sm_opti, dim = 2)
      do iTime = 1, n_time_steps

        ! check for enough data points in timesteps for KGE calculation
        ! more then 10 percent avaiable in current field
        if (count(L1_sm_mask(s1 : e1, iTime)) .LE. (0.10_dp * real(nCells1, dp))) then
          invalid_times = invalid_times + 1.0_dp
          mask_times(iTime) = .FALSE.
          cycle
        end if
        sm_catch_avg_basin(iTime) = average(L1_sm(s1 : e1, iTime), mask = L1_sm_mask(s1 : e1, iTime))
        sm_opti_catch_avg_basin(iTime) = average(sm_opti(s1 : e1, iTime), mask = L1_sm_mask(s1 : e1, iTime))
      end do

      ! user information about invalid times
      if (invalid_times .GT. 0.5_dp) then
        call message('   WARNING: objective_sm: Detected invalid timesteps (.LT. 10 valid data points).')
        call message('                          Fraction of invlauid timestpes: ', &
                num2str(invalid_times / real(n_time_steps, dp), '(F4.2)'))
      end if


      ! calculate average soil moisture KGE over all basins with power law
      ! basins are weighted equally ( 1 / real(nBasin,dp))**6
      objective_sm_kge_catchment_avg = objective_sm_kge_catchment_avg + &
              ((1.0_dp - KGE(sm_catch_avg_basin, sm_opti_catch_avg_basin, mask = mask_times)) / real(nBasins, dp))**6

      ! deallocate
      deallocate(mask_times)
      deallocate(sm_catch_avg_basin)
      deallocate(sm_opti_catch_avg_basin)
    end do

    objective_sm_kge_catchment_avg = objective_sm_kge_catchment_avg**onesixth

    call message('    objective_sm_kge_catchment_avg = ', num2str(objective_sm_kge_catchment_avg, '(F9.5)'))

  END FUNCTION objective_sm_kge_catchment_avg

  ! ------------------------------------------------------------------

  !    NAME
  !        objective_sm_corr

  !    PURPOSE
  !>       \brief Objective function for soil moisture.

  !>       \details The objective function only depends on a parameter vector.
  !>       The model will be called with that parameter vector and
  !>       the model output is subsequently compared to observed data.

  !>       Therefore the Pearson correlation between observed and modeled soil
  !>       moisture on each grid cell \f$ j \f$ is compared
  !>       \f[ r_j = r^2(SM_{obs}^j, SM_{sim}^j) \f]
  !>       where
  !>       \f$ r^2\f$        = Pearson correlation coefficient,
  !>       \f$ SM_{obs} \f$  = observed soil moisture,
  !>       \f$ SM_{sim}  \f$ = simulated soil moisture.
  !>       The observed data \f$ SM_{obs} \f$ are global in this module.

  !>       The the correlation is spatially averaged as
  !>       \f[ \phi_{i} = \frac{1}{K} \cdot \sum_{j=1}^K r_j \f]
  !>       where \f$ K \f$ denotes the number of valid cells in the study domain.
  !>       Finally, the overall objective function value \f$ OF \f$ is estimated based on the power-6
  !>       norm to combine the \f$ \phi_{i} \f$ from all basins \f$ N \f$.
  !>       \f[ OF = \sqrt[6]{\sum((1.0 - \phi_{i})/N)^6 }. \f]
  !>       The observed data L1_sm, L1_sm_mask are global in this module.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp) :: objective_sm_corr &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date March 2015

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION objective_sm_corr(parameterset, eval)

    use mo_common_variables, only : level1, nBasins
    use mo_global_variables, only : L1_sm, L1_sm_mask
    use mo_message, only : message
    use mo_moment, only : correlation
    use mo_string_utils, only : num2str

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), POINTER :: eval

    real(dp) :: objective_sm_corr

    ! basin loop counter
    integer(i4) :: iBasin

    ! cell loop counter
    integer(i4) :: iCell

    ! start and end index for the current basin
    integer(i4) :: s1, e1

    ! ncells1 of level 1
    integer(i4) :: ncells1

    ! number of invalid cells in catchment
    real(dp) :: invalid_cells

    ! basins wise objectives
    real(dp) :: objective_sm_corr_basin

    real(dp), parameter :: onesixth = 1.0_dp / 6.0_dp

    ! simulated soil moisture
    ! (dim1=ncells, dim2=time)
    real(dp), dimension(:, :), allocatable :: sm_opti


    call eval(parameterset, sm_opti = sm_opti)

    ! initialize some variables
    objective_sm_corr = 0.0_dp

    ! loop over basin - for applying power law later on
    do iBasin = 1, nBasins

      ! init
      objective_sm_corr_basin = 0.0_dp
      ! get basin information
      ncells1 = level1(iBasin)%ncells
      s1 = level1(iBasin)%iStart
      e1 = level1(iBasin)%iEnd

      invalid_cells = 0.0_dp
      ! temporal correlation is calculated on individual gridd cells

      do iCell = s1, e1

        ! check for enough data points in time for correlation
        if (count(L1_sm_mask(iCell, :)) .LE. 0.10_dp * real(size(L1_sm, dim = 2), dp)) then
          invalid_cells = invalid_cells + 1.0_dp
          cycle
        end if
        objective_sm_corr_basin = objective_sm_corr_basin + &
                correlation(L1_sm(iCell, :), sm_opti(iCell, :), mask = L1_sm_mask(iCell, :))
      end do

      ! user information about invalid cells
      if (invalid_cells .GT. 0.5_dp) then
        call message('   WARNING: objective_sm: Detected invalid cells in study area (.LT. 10 valid data points).')
        call message('                          Fraction of invlauid cells: ', &
                num2str(invalid_cells / real(nCells1, dp), '(F4.2)'))
      end if


      ! calculate average soil moisture correlation over all basins with power law
      ! basins are weighted equally ( 1 / real(nBasin,dp))**6
      objective_sm_corr = objective_sm_corr + &
              ((1.0_dp - objective_sm_corr_basin / real(nCells1, dp)) / real(nBasins, dp))**6
    end do

    objective_sm_corr = objective_sm_corr**onesixth

    call message('    objective_sm_corr = ', num2str(objective_sm_corr, '(F9.5)'))

  END FUNCTION objective_sm_corr

  ! ------------------------------------------------------------------

  !    NAME
  !        objective_sm_pd

  !    PURPOSE
  !>       \brief Objective function for soil moisture.

  !>       \details The objective function only depends on a parameter vector.
  !>       The model will be called with that parameter vector and
  !>       the model output is subsequently compared to observed data.

  !>       Therefore the Pattern Dissimilarity (PD) of observed and modeled soil
  !>       moisture fields is calculated - aim: matching spatial patters
  !>       \f[ E(t) = PD\left( SM_{obs}(t), SM_{sim}(t) \right) \f]
  !>       where
  !>       \f$ PD \f$        = pattern dissimilarity function,
  !>       \f$ SM_{obs} \f$  = observed soil moisture,
  !>       \f$ SM_{sim}  \f$ = simulated soil moisture.
  !>       \f$ E(t)  \f$     = pattern dissimilarity at timestep \f$ t \f$.
  !>       The the pattern dissimilaity (E) is spatially averaged as
  !>       \f[ \phi_{i} = \frac{1}{T} \cdot \sum_{t=1}^T E_t \f]
  !>       where \f$ T \f$ denotes the number of time steps.
  !>       Finally, the overall objective function value \f$ OF \f$ is estimated based on the power-6
  !>       norm to combine the \f$ \phi_{i} \f$ from all basins \f$ N \f$.
  !>       \f[ OF = \sqrt[6]{\sum((1.0 - \phi_{i})/N)^6 } . \f]
  !>       The observed data L1_sm, L1_sm_mask are global in this module.
  !>       The observed data L1_sm, L1_sm_mask are global in this module.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp) :: objecive_sm_pd &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date May 2015

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION objective_sm_pd(parameterset, eval)

    use mo_common_constants, only : nodata_dp
    use mo_common_variables, only : level1, nBasins
    use mo_global_variables, only : L1_sm, L1_sm_mask
    use mo_message, only : message
    use mo_spatialsimilarity, only : PD
    use mo_string_utils, only : num2str

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), POINTER :: eval

    ! objective function value
    real(dp) :: objective_sm_pd

    ! basin loop counter
    integer(i4) :: iBasin

    ! time loop counter
    integer(i4) :: iTime

    ! level 1 number of culomns and rows
    integer(i4) :: nrows1, ncols1

    ! start and end index for the current basin
    integer(i4) :: s1, e1

    ! for sixth root
    real(dp), parameter :: onesixth = 1.0_dp / 6.0_dp

    ! matrices of SM from vectorized arrays
    real(dp), dimension(:, :), allocatable :: mat1, mat2

    ! pattern dissimilarity (pd) at every time step
    real(dp), dimension(:), allocatable :: pd_time_series

    ! simulated soil moisture
    ! (dim1=ncells, dim2=time)
    real(dp), dimension(:, :), allocatable :: sm_opti

    ! mask of valid cells at level1
    logical, dimension(:, :), allocatable :: mask1

    ! mask of valid sm cells
    logical, dimension(:, :), allocatable :: mask_sm

    ! mask for valid sm catchment avg time steps
    logical, dimension(:), allocatable :: mask_times


    call eval(parameterset, sm_opti = sm_opti)

    ! initialize some variables
    objective_sm_pd = 0.0_dp

    ! loop over basin - for applying power law later on
    do iBasin = 1, nBasins

      ! get basin information
      mask1 = level1(iBasin)%mask
      ncols1 = level1(iBasin)%ncols
      nrows1 = level1(iBasin)%nrows
      s1 = level1(iBasin)%iStart
      e1 = level1(iBasin)%iEnd

      ! allocate
      allocate(mask_times    (size(sm_opti, dim = 2)))
      allocate(pd_time_series(size(sm_opti, dim = 2)))
      allocate(mat1   (nrows1, ncols1))
      allocate(mat2   (nrows1, ncols1))
      allocate(mask_sm(nrows1, ncols1))

      ! initalize
      mask_times = .FALSE.
      pd_time_series = 0.0_dp

      ! calculate pattern similarity criterion
      do iTime = 1, size(sm_opti, dim = 2)
        mat1 = unpack(L1_sm(s1 : e1, iTime), mask1, nodata_dp)
        mat2 = unpack(sm_opti(s1 : e1, iTime), mask1, nodata_dp)
        mask_sm = unpack(L1_sm_mask(s1 : e1, iTime), mask1, .FALSE.)
        pd_time_series = PD(mat1, mat2, mask = mask_sm, valid = mask_times(itime))
      end do

      if (count(mask_times) > 0_i4) then
        ! calculate avergae PD over all basins with power law -basins are weighted equally ( 1 / real(nBasin,dp))**6
        objective_sm_pd = objective_sm_pd + &
                ((1.0_dp - sum(pd_time_series, mask = mask_times) / real(count(mask_times), dp)) / real(nBasins, dp))**6
      else
        call message('***ERROR: mo_objective_funtion: objective_sm_pd: No soil moisture observations available!')
        stop
      end if

      ! deallocate
      deallocate(mask_times)
      deallocate(pd_time_series)
      deallocate(mat1)
      deallocate(mat2)
      deallocate(mask_sm)
    end do

    objective_sm_pd = objective_sm_pd**onesixth

    call message('    objective_sm_pd = ', num2str(objective_sm_pd, '(F9.5)'))

  END FUNCTION objective_sm_pd

  ! ------------------------------------------------------------------

  !    NAME
  !        objective_sm_sse_standard_score

  !    PURPOSE
  !>       \brief Objective function for soil moisture.

  !>       \details The objective function only depends on a parameter vector.
  !>       The model will be called with that parameter vector and
  !>       the model output is subsequently compared to observed data.

  !>       Therefore the sum of squared errors (SSE) of the standard score of observed and
  !>       modeled soil moisture is calculated. The standard score or normalization (anomaly)
  !>       make the objctive function bias insensitive and basically the dynamics of the soil moisture
  !>       is tried to capture by this obejective function.
  !>       \f[ phi_i = \sum_{j=1}^K \{ standard\_score( SM_{obs}(j) )- standard\_score(SM_{sim}(j)) \}^2 \f]
  !>       where
  !>       \f$  standard\_score \f$ = standard score function,
  !>       \f$ SM_{obs} \f$  = observed soil moisture,
  !>       \f$ SM_{sim}  \f$ = simulated soil moisture.
  !>       \f$ K  \f$ = valid elements in study domain.
  !>       Finally, the overall objective function value \f$ OF \f$ is estimated based on the power-6
  !>       norm to combine the \f$ \phi_{i} \f$ from all basins \f$ N \f$.
  !>       \f[ OF = \sqrt[6]{\sum(\phi_{i}/N)^6 }.  \f]
  !>       The observed data L1_sm, L1_sm_mask are global in this module.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp) :: objective_sm_sse_standard_score &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date March 2015

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION objective_sm_sse_standard_score(parameterset, eval)

    use mo_common_variables, only : level1, nBasins
    use mo_errormeasures, only : SSE
    use mo_global_variables, only : L1_sm, L1_sm_mask
    use mo_message, only : message
    use mo_standard_score, only : standard_score
    use mo_string_utils, only : num2str

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), POINTER :: eval

    real(dp) :: objective_sm_sse_standard_score

    ! basin loop counter
    integer(i4) :: iBasin

    ! cell loop counter
    integer(i4) :: iCell

    ! start and end index for the current basin
    integer(i4) :: s1, e1

    ! ncells1 of level 1
    integer(i4) :: ncells1

    ! number of invalid cells in catchment
    real(dp) :: invalid_cells

    ! basins wise objectives
    real(dp) :: objective_sm_sse_standard_score_basin

    real(dp), parameter :: onesixth = 1.0_dp / 6.0_dp

    ! simulated soil moisture
    ! (dim1=ncells, dim2=time)
    real(dp), dimension(:, :), allocatable :: sm_opti


    call eval(parameterset, sm_opti = sm_opti)

    ! initialize some variables
    objective_sm_sse_standard_score = 0.0_dp

    ! loop over basin - for applying power law later on
    do iBasin = 1, nBasins

      ! init
      objective_sm_sse_standard_score_basin = 0.0_dp
      ! get basin information
      nCells1 = level1(iBasin)%nCells
      s1 = level1(iBasin)%iStart
      e1 = level1(iBasin)%iEnd

      invalid_cells = 0.0_dp
      ! standard_score signal is calculated on individual grid cells
      do iCell = s1, e1

        ! check for enough data points in time for statistical calculations (e.g. mean, stddev)
        if (count(L1_sm_mask(iCell, :)) .LE. (0.10_dp * real(size(L1_sm, dim = 2), dp))) then
          invalid_cells = invalid_cells + 1.0_dp
          cycle
        end if
        objective_sm_sse_standard_score_basin = objective_sm_sse_standard_score_basin + &
                SSE(standard_score(L1_sm(iCell, :), mask = L1_sm_mask(iCell, :)), &
                        standard_score(sm_opti(iCell, :), mask = L1_sm_mask(iCell, :)), mask = L1_sm_mask(iCell, :))

      end do

      ! user information about invalid cells
      if (invalid_cells .GT. 0.5_dp) then
        call message('   WARNING: objective_sm: Detected invalid cells in study area (.LT. 10 valid data points).')
        call message('                          Fraction of invlauid cells: ', &
                num2str(invalid_cells / real(nCells1, dp), '(F4.2)'))
      end if

      ! calculate average soil moisture correlation over all basins with power law
      ! basins are weighted equally ( 1 / real(nBasin,dp))**6
      objective_sm_sse_standard_score = objective_sm_sse_standard_score + &
              (objective_sm_sse_standard_score_basin / real(nBasins, dp))**6
    end do

    objective_sm_sse_standard_score = objective_sm_sse_standard_score**onesixth

    call message('    objective_sm_sse_standard_score = ', num2str(objective_sm_sse_standard_score, '(E12.5)'))

  END FUNCTION objective_sm_sse_standard_score


  ! -----------------------------------------------------------------

  !    NAME
  !        objective_kge_q_rmse_tws

  !    PURPOSE
  !>       \brief Objective function of KGE for runoff and RMSE for basin_avg TWS (standarized scores)

  !>       \details Objective function of KGE for runoff and RMSE for basin_avg TWS (standarized scores)

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp) :: objective_kge_q_rmse_tws &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !    HISTORY
  !>       \authors Oldrich Rakovec, Rohini Kumar

  !>       \date Oct. 2015

  ! Modifications:
  ! Stephan Thober Oct 2015 - moved tws optimization from mo_mrm_objective_function_runoff here
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION objective_kge_q_rmse_tws(parameterset, eval)

    use mo_common_constants, only : eps_dp, nodata_dp
    use mo_common_mhm_mrm_variables, only : evalPer
    use mo_common_variables, only : nBasins
    use mo_errormeasures, only : rmse
    use mo_julian, only : caldat
    use mo_message, only : message
    use mo_moment, only : mean
    use mo_standard_score, only : classified_standard_score
    use mo_string_utils, only : num2str
    use mo_temporal_aggregation, only : day2mon_average
#ifdef MRM2MHM
    use mo_errormeasures, only : kge
    use mo_mrm_objective_function_runoff, only : extract_runoff
#endif

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), POINTER :: eval

    real(dp) :: objective_kge_q_rmse_tws

    ! modelled runoff for a given parameter set
    ! dim1=nTimeSteps, dim2=nGauges
    real(dp), allocatable, dimension(:, :) :: runoff

    ! modelled runoff for a given parameter set
    ! dim1=nTimeSteps, dim2=nGauges
    real(dp), allocatable, dimension(:, :) :: tws

    ! basin counter, month counters
    integer(i4) :: ii, pp, mmm

    integer(i4) :: year, month, day

    real(dp), dimension(nBasins) :: initTime

    ! simulated tws
    real(dp), dimension(:), allocatable :: tws_sim

    ! measured tws
    real(dp), dimension(:), allocatable :: tws_obs

    ! mask for measured tws
    logical, dimension(:), allocatable :: tws_obs_mask

    ! total number of months
    integer(i4) :: nMonths

    ! vector with months' classes
    integer(i4), dimension(:), allocatable :: month_classes

    ! monthly values original time series
    real(dp), DIMENSION(:), allocatable :: tws_sim_m, tws_obs_m

    ! monthly values anomaly time series
    real(dp), DIMENSION(:), allocatable :: tws_sim_m_anom, tws_obs_m_anom

    logical, DIMENSION(:), allocatable :: tws_obs_m_mask

    ! rmse_tws(nBasins)
    real(dp), dimension(:), allocatable :: rmse_tws

    ! obj. functions
    real(dp) :: rmse_tws_avg, kge_q_avg

#ifdef MRM2MHM
    integer(i4) :: nGaugesTotal

    ! aggregated simulated runoff
    real(dp), dimension(:), allocatable :: runoff_agg

    ! measured runoff
    real(dp), dimension(:), allocatable :: runoff_obs

    ! mask for measured runoff
    logical, dimension(:), allocatable :: runoff_obs_mask

    ! kge_q(nGaugesTotal)
    real(dp), dimension(:), allocatable :: kge_q

    ! gauges counter
    integer(i4) :: gg
#endif

    ! obtain hourly values of runoff and tws:
    call eval(parameterset, runoff = runoff, basin_avg_tws = tws)

    !--------------------------------------------
    !! TWS
    !--------------------------------------------

    ! allocate
    allocate(rmse_tws(nBasins))
    rmse_tws(:) = nodata_dp

    do ii = 1, nBasins

      ! extract tws the same way as runoff using mrm
      call extract_basin_avg_tws(ii, tws, tws_sim, tws_obs, tws_obs_mask)

      ! check for potentially 2 years of data
      if (count(tws_obs_mask) .lt.  365 * 2) then
        call message('objective_kge_q_rmse_tws: Length of TWS data of Basin ', trim(adjustl(num2str(ii))), &
                ' less than 2 years: this is not recommended')
      end if

      ! get initial time of the evaluation period
      initTime(ii) = real(evalPer(ii)%julStart, dp)

      ! get calendar days, months, year
      call caldat(int(initTime(ii)), yy = year, mm = month, dd = day)

      ! calculate monthly averages from daily values of the model
      call day2mon_average(tws_sim, year, month, day, tws_sim_m, misval = nodata_dp)

      ! remove mean from modelled time series
      tws_sim_m(:) = tws_sim_m(:) - mean(tws_sim_m(:))

      ! calculate monthly averages from daily values of the observations, which already have removed the long-term mean
      call day2mon_average(tws_obs, year, month, day, tws_obs_m, misval = nodata_dp)

      ! get number of months for given basin
      nMonths = size(tws_obs_m)

      allocate (month_classes(nMonths))
      allocate (tws_obs_m_mask(nMonths))
      allocate (tws_obs_m_anom(nMonths))
      allocate (tws_sim_m_anom(nMonths))

      month_classes(:) = 0
      tws_obs_m_mask(:) = .TRUE.
      tws_obs_m_anom(:) = nodata_dp
      tws_sim_m_anom(:) = nodata_dp

      ! define months' classes
      mmm = month
      do pp = 1, nMonths
        month_classes(pp) = mmm
        if (mmm .LT. 12) then
          mmm = mmm + 1
        else
          mmm = 1
        end if
      end do

      ! define mask for missing data in observations (there are always data for simulations)
      where(abs(tws_obs_m - nodata_dp) .lt. eps_dp) tws_obs_m_mask = .FALSE.

      ! calculate standard score
      tws_obs_m_anom = classified_standard_score(tws_obs_m, month_classes, mask = tws_obs_m_mask)
      tws_sim_m_anom = classified_standard_score(tws_sim_m, month_classes, mask = tws_obs_m_mask)
      rmse_tws(ii) = rmse(tws_sim_m_anom, tws_obs_m_anom, mask = tws_obs_m_mask)

      deallocate (month_classes)
      deallocate (tws_obs_m)
      deallocate (tws_sim_m)
      deallocate (tws_obs_m_mask)
      deallocate (tws_sim_m_anom)
      deallocate (tws_obs_m_anom)
      deallocate (tws_sim, tws_obs, tws_obs_mask)

    end do

    rmse_tws_avg = sum(rmse_tws(:), abs(rmse_tws - nodata_dp) .gt. eps_dp) / &
            real(count(abs(rmse_tws - nodata_dp) .gt. eps_dp), dp)
    deallocate(rmse_tws)

    !--------------------------------------------
    !! RUNOFF
    !--------------------------------------------
    kge_q_avg = 0_dp
#ifdef MRM2MHM
    nGaugesTotal = size(runoff, dim = 2)
    allocate(kge_q(nGaugesTotal))
    kge_q(:) = nodata_dp

    do gg = 1, nGaugesTotal

      ! extract runoff
      call extract_runoff(gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask)

      ! check for potentially 2 years of data
      pp = count(runoff_agg .ge. 0.0_dp)
      if (pp .lt.  365 * 2) then
        call message('objective_kge_q_rmse_tws: Length of TWS data of Basin ', trim(adjustl(num2str(ii))), &
        ' less than 2 years: this is not recommended')
      end if
      ! calculate KGE for each basin:
      kge_q(gg) = kge(runoff_obs, runoff_agg, mask = runoff_obs_mask)
      deallocate (runoff_agg, runoff_obs, runoff_obs_mask)

    end do

    ! calculate average KGE value for runoff
    kge_q_avg = sum(kge_q(:), abs(kge_q - nodata_dp) .gt. eps_dp) / &
            real(count(abs(kge_q - nodata_dp) .gt. eps_dp), dp)
    deallocate(kge_q)
#else
    call message('***ERROR: objective_kge_q_rmse_tws: missing routing module for optimization')
    stop 1
#endif

    !
    objective_kge_q_rmse_tws = rmse_tws_avg * (1._dp - kge_q_avg)

    call message('    objective_kge_q_rmse_tws = ', num2str(objective_kge_q_rmse_tws, '(F9.5)'))

  END FUNCTION objective_kge_q_rmse_tws

  ! ------------------------------------------------------------------

  !    NAME
  !        objective_neutrons_kge_catchment_avg

  !    PURPOSE
  !>       \brief Objective function for neutrons.

  !>       \details The objective function only depends on a parameter vector.
  !>       The model will be called with that parameter vector and
  !>       the model output is subsequently compared to observed data.

  !>       Therefore, the Kling-Gupta model efficiency \f$ KGE \f$ of the catchment average
  !>       neutrons (N) is calculated
  !>       \f[ KGE = 1.0 - \sqrt{( (1-r)^2 + (1-\alpha)^2 + (1-\beta)^2 )} \f]
  !>       where
  !>       \f$ r \f$ = Pearson product-moment CORRELATION coefficient
  !>       \f$ \alpha \f$ = ratio of simulated mean to observed mean SM
  !>       \f$ \beta  \f$ = ratio of similated standard deviation to observed standard deviation
  !>       is calculated and the objective function for a given basin \f$ i \f$ is
  !>       \f[ \phi_{i} = 1.0 - KGE_{i} \f]
  !>       \f$ \phi_{i} \f$ is the objective since we always apply minimization methods.
  !>       The minimal value of \f$ \phi_{i} \f$ is 0 for the optimal KGE of 1.0.

  !>       Finally, the overall objective function value \f$ OF \f$ is estimated based on the power-6
  !>       norm to combine the \f$ \phi_{i} \f$ from all basins \f$ N \f$.
  !>       \f[ OF = \sqrt[6]{\sum((1.0 - KGE_{i})/N)^6 }.  \f]
  !>       The observed data L1_neutronsdata, L1_neutronsdata_mask are global in this module.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp) :: objective_neutrons_kge_catchment_avg &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine)

  !    HISTORY
  !>       \authors Martin Schroen

  !>       \date Jun 2015

  ! Modifications:
  ! Maren Kaluza Mar 2018 - changed format string to '(I10)'
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION objective_neutrons_kge_catchment_avg(parameterset, eval)

    use mo_common_constants, only : nodata_dp
    use mo_common_variables, only : level1, nBasins
    use mo_errormeasures, only : KGE
    use mo_global_variables, only : L1_neutronsdata, L1_neutronsdata_mask
    use mo_message, only : message
    use mo_moment, only : average
    use mo_string_utils, only : num2str

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), POINTER :: eval

    real(dp) :: objective_neutrons_kge_catchment_avg

    ! basin loop counter
    integer(i4) :: iBasin

    ! time loop counter
    integer(i4) :: iTime

    ! start and end index for the current basin
    integer(i4) :: s1, e1

    ! for sixth root
    real(dp), parameter :: onesixth = 1.0_dp / 6.0_dp

    ! spatial average of observed neutrons
    real(dp), dimension(:), allocatable :: neutrons_catch_avg_basin

    ! spatial avergae of modeled  neutrons
    real(dp), dimension(:), allocatable :: neutrons_opti_catch_avg_basin

    ! simulated neutrons
    ! (dim1=ncells, dim2=time)
    real(dp), dimension(:, :), allocatable :: neutrons_opti

    ! mask for valid neutrons catchment avg time steps
    logical, dimension(:), allocatable :: mask_times


    call eval(parameterset, neutrons_opti = neutrons_opti)

    ! initialize some variables
    objective_neutrons_kge_catchment_avg = 0.0_dp

    ! loop over basin - for applying power law later on
    do iBasin = 1, nBasins

      ! get basin information
      s1 = level1(iBasin)%iStart
      e1 = level1(iBasin)%iEnd

      ! allocate
      allocate(mask_times             (size(neutrons_opti, dim = 2)))
      allocate(neutrons_catch_avg_basin     (size(neutrons_opti, dim = 2)))
      allocate(neutrons_opti_catch_avg_basin(size(neutrons_opti, dim = 2)))

      ! initalize
      mask_times = .TRUE.
      neutrons_catch_avg_basin = nodata_dp
      neutrons_opti_catch_avg_basin = nodata_dp

      ! calculate catchment average soil moisture
      do iTime = 1, size(neutrons_opti, dim = 2)

        ! check for enough data points in time for correlation
        if (all(.NOT. L1_neutronsdata_mask(s1 : e1, iTime))) then
          call message('WARNING: neutrons data at time ', num2str(iTime, '(I10)'), ' is empty.')
          !call message('WARNING: objective_neutrons_kge_catchment_avg: ignored currrent time step since less than')
          !call message('         10 valid cells available in soil moisture observation')
          mask_times(iTime) = .FALSE.
          cycle
        end if
        neutrons_catch_avg_basin(iTime) = average(L1_neutronsdata(s1 : e1, iTime), mask = L1_neutronsdata_mask(s1 : e1, iTime))
        neutrons_opti_catch_avg_basin(iTime) = average(neutrons_opti(s1 : e1, iTime), mask = L1_neutronsdata_mask(s1 : e1, iTime))
      end do

      ! calculate average neutrons KGE over all basins with power law
      ! basins are weighted equally ( 1 / real(nBasin,dp))**6
      objective_neutrons_kge_catchment_avg = objective_neutrons_kge_catchment_avg + &
              ((1.0_dp - KGE(neutrons_catch_avg_basin, neutrons_opti_catch_avg_basin, mask = mask_times)) / real(nBasins, dp))**6

      ! deallocate
      deallocate(mask_times)
      deallocate(neutrons_catch_avg_basin)
      deallocate(neutrons_opti_catch_avg_basin)

    end do

    objective_neutrons_kge_catchment_avg = objective_neutrons_kge_catchment_avg**onesixth

    call message('    objective_neutrons_kge_catchment_avg = ', num2str(objective_neutrons_kge_catchment_avg, '(F9.5)'))

  END FUNCTION objective_neutrons_kge_catchment_avg

  ! ------------------------------------------------------------------

  !    NAME
  !        objective_et_kge_catchment_avg

  !    PURPOSE
  !>       \brief Objective function for evpotranspirstion (et).

  !>       \details The objective function only depends on a parameter vector.
  !>       The model will be called with that parameter vector and
  !>       the model output is subsequently compared to observed data.

  !>       Therefore, the Kling-Gupta model efficiency \f$ KGE \f$ of the catchment average
  !>       evapotranspiration (et) is calculated
  !>       \f[ KGE = 1.0 - \sqrt{( (1-r)^2 + (1-\alpha)^2 + (1-\beta)^2 )} \f]
  !>       where
  !>       \f$ r \f$ = Pearson product-moment correlation coefficient
  !>       \f$ \alpha \f$ = ratio of simulated mean to observed mean SM
  !>       \f$ \beta  \f$ = ratio of similated standard deviation to observed standard deviation
  !>       is calculated and the objective function for a given basin \f$ i \f$ is
  !>       \f[ \phi_{i} = 1.0 - KGE_{i} \f]
  !>       \f$ \phi_{i} \f$ is the objective since we always apply minimization methods.
  !>       The minimal value of \f$ \phi_{i} \f$ is 0 for the optimal KGE of 1.0.

  !>       Finally, the overall objective function value \f$ OF \f$ is estimated based on the power-6
  !>       norm to combine the \f$ \phi_{i} \f$ from all basins \f$ N \f$.
  !>       \f[ OF = \sqrt[6]{\sum((1.0 - KGE_{i})/N)^6 }.  \f]
  !>       The observed data L1_et, L1_et_mask are global in this module.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp) :: objective_et_kge_catchment_avg &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine)

  !    HISTORY
  !>       \authors Johannes Brenner

  !>       \date Feb 2017

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION objective_et_kge_catchment_avg(parameterset, eval)

    use mo_common_constants, only : nodata_dp
    use mo_common_variables, only : level1, nBasins
    use mo_errormeasures, only : KGE
    use mo_global_variables, only : L1_et, L1_et_mask
    use mo_message, only : message
    use mo_moment, only : average
    use mo_string_utils, only : num2str

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), POINTER :: eval

    real(dp) :: objective_et_kge_catchment_avg

    ! basin loop counter
    integer(i4) :: iBasin

    ! time loop counter
    integer(i4) :: iTime

    ! start and end index for the current basin
    integer(i4) :: s1, e1

    ! for sixth root
    real(dp), parameter :: onesixth = 1.0_dp / 6.0_dp

    ! spatial average of observed et
    real(dp), dimension(:), allocatable :: et_catch_avg_basin

    ! spatial avergae of modeled  et
    real(dp), dimension(:), allocatable :: et_opti_catch_avg_basin

    ! simulated et
    ! (dim1=ncells, dim2=time)
    real(dp), dimension(:, :), allocatable :: et_opti

    ! mask for valid et catchment avg time steps
    logical, dimension(:), allocatable :: mask_times


    call eval(parameterset, et_opti = et_opti)

    ! initialize some variables
    objective_et_kge_catchment_avg = 0.0_dp

    ! loop over basin - for applying power law later on
    do iBasin = 1, nBasins

      ! get basin information
      s1 = level1(iBasin)%iStart
      e1 = level1(iBasin)%iEnd

      ! allocate
      allocate(mask_times             (size(et_opti, dim = 2)))
      allocate(et_catch_avg_basin     (size(et_opti, dim = 2)))
      allocate(et_opti_catch_avg_basin(size(et_opti, dim = 2)))

      ! initalize
      mask_times = .TRUE.
      et_catch_avg_basin = nodata_dp
      et_opti_catch_avg_basin = nodata_dp

      ! calculate catchment average evapotranspiration
      do iTime = 1, size(et_opti, dim = 2)

        ! check for enough data points in time for correlation
        if (all(.NOT. L1_et_mask(s1 : e1, iTime))) then
          !write (*,*) 'WARNING: et data at time ', iTime, ' is empty.'
          !call message('WARNING: objective_et_kge_catchment_avg: ignored currrent time step since less than')
          !call message('         10 valid cells available in evapotranspiration observation')
          mask_times(iTime) = .FALSE.
          cycle
        end if

        et_catch_avg_basin(iTime) = average(L1_et(s1 : e1, iTime), mask = L1_et_mask(s1 : e1, iTime))
        et_opti_catch_avg_basin(iTime) = average(et_opti(s1 : e1, iTime), mask = L1_et_mask(s1 : e1, iTime))
      end do

      ! calculate average ET KGE over all basins with power law
      ! basins are weighted equally ( 1 / real(nBasin,dp))**6

      objective_et_kge_catchment_avg = objective_et_kge_catchment_avg + &
              ((1.0_dp - KGE(et_catch_avg_basin, et_opti_catch_avg_basin, mask = mask_times)) / real(nBasins, dp))**6

      ! deallocate
      deallocate(mask_times)
      deallocate(et_catch_avg_basin)
      deallocate(et_opti_catch_avg_basin)
    end do

    objective_et_kge_catchment_avg = objective_et_kge_catchment_avg**onesixth

    call message('    objective_et_kge_catchment_avg = ', num2str(objective_et_kge_catchment_avg, '(F9.5)'))

  END FUNCTION objective_et_kge_catchment_avg

  ! -----------------------------------------------------------------

  !    NAME
  !        objective_kge_q_sm_corr

  !    PURPOSE
  !>       \brief Objective function of KGE for runoff and correlation for SM

  !>       \details Objective function of KGE for runoff and SSE for soil moisture (standarized scores).
  !>       Further details can be found in the documentation of objective functions
  !>       '14 - objective_multiple_gauges_kge_power6' and '13 - objective_sm_corr'.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp) :: objective_kge_q_sse_sm &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date Mar. 2017

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION objective_kge_q_sm_corr(parameterset, eval)

    use mo_common_variables, only : level1, nBasins
    use mo_global_variables, only : L1_sm, L1_sm_mask
    use mo_message, only : message
    use mo_moment, only : correlation
    use mo_string_utils, only : num2str
#ifdef MRM2MHM
    use mo_errormeasures, only : kge
    use mo_mrm_objective_function_runoff, only : extract_runoff
#endif

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), POINTER :: eval

    real(dp) :: objective_kge_q_sm_corr

    real(dp) :: objective_sm

    real(dp) :: objective_kge

    ! number of invalid cells in catchment
    real(dp) :: invalid_cells

    ! modelled runoff for a given parameter set
    ! dim1=nTimeSteps, dim2=nGauges
    real(dp), allocatable, dimension(:, :) :: runoff

    ! basin loop counter
    integer(i4) :: iBasin

    ! cell loop counter
    integer(i4) :: iCell

    ! start and end index for the current basin
    integer(i4) :: s1, e1

    ! ncells1 of level 1
    integer(i4) :: ncells1

    ! basins wise objectives
    real(dp) :: objective_sm_basin

    ! simulated soil moisture
    ! (dim1=ncells, dim2=time)
    real(dp), dimension(:, :), allocatable :: sm_opti

    real(dp), parameter :: onesixth = 1.0_dp / 6.0_dp

#ifdef MRM2MHM
    ! gauges counter
    integer(i4) :: gg

    integer(i4) :: nGaugesTotal

    ! aggregated simulated runoff
    real(dp), dimension(:), allocatable :: runoff_agg

    ! measured runoff
    real(dp), dimension(:), allocatable :: runoff_obs

    ! mask for measured runoff
    logical, dimension(:), allocatable :: runoff_obs_mask
#endif

    ! run mHM
    call eval(parameterset, runoff = runoff, sm_opti = sm_opti)

    ! -----------------------------
    ! SOIL MOISTURE
    ! -----------------------------

    ! initialize some variables
    objective_sm = 0.0_dp

    ! loop over basin - for applying power law later on
    do iBasin = 1, nBasins

      ! init
      objective_sm_basin = 0.0_dp
      ! get basin information
      nCells1 = level1(iBasin)%nCells
      s1 = level1(iBasin)%iStart
      e1 = level1(iBasin)%iEnd


      ! correlation signal is calculated on individual grid cells
      invalid_cells = 0.0_dp
      do iCell = s1, e1

        ! check for enough data points in time for statistical calculations (e.g. mean, stddev)
        if (count(L1_sm_mask(iCell, :)) .LE. (0.10_dp * real(size(L1_sm, dim = 2), dp))) then
          invalid_cells = invalid_cells + 1.0_dp
          cycle
        end if

        ! calculate ojective function
        objective_sm_basin = objective_sm_basin + &
                correlation(L1_sm(iCell, :), sm_opti(iCell, :), mask = L1_sm_mask(iCell, :))
      end do

      ! user information about invalid cells
      if (invalid_cells .GT. 0.5_dp) then
        call message('   WARNING: objective_sm: Detected invalid cells in study area (.LT. 10 valid data points).')
        call message('                          Fraction of invlauid cells: ', &
                num2str(invalid_cells / real(nCells1, dp), '(F4.2)'))
      end if

      ! calculate average soil moisture objective over all basins with power law
      ! basins are weighted equally ( 1 / real(nBasin,dp))**6
      objective_sm = objective_sm + &
              ((1.0_dp - objective_sm_basin / real(nCells1, dp)) / real(nBasins, dp))**6
    end do

    ! compromise solution - sixth root
    objective_sm = objective_sm**onesixth

    ! -----------------------------
    ! RUNOFF
    ! -----------------------------
    objective_kge = 0.0_dp
#ifdef MRM2MHM
    nGaugesTotal = size(runoff, dim = 2)

    do gg = 1, nGaugesTotal

      ! extract runoff
      call extract_runoff(gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask)

      ! KGE
      objective_kge = objective_kge + &
              ((1.0_dp - kge(runoff_obs, runoff_agg, mask = runoff_obs_mask)) / real(nGaugesTotal, dp))**6

    end do

    deallocate(runoff_agg, runoff_obs, runoff_obs_mask)

    ! compromise solution - sixth root
    objective_kge = objective_kge**onesixth

#else
    call message('***ERROR: objective_kge_q_rmse_tws: missing routing module for optimization')
    stop
#endif

    ! equal weighted compromise objective functions for discharge and soilmoisture
    objective_kge_q_sm_corr = (objective_sm**6 + objective_kge**6)**onesixth
    !    print*, "1-SM 2-Q : ", 1.0_dp-objective_sm, 1.0_dp-objective_kge ! MZMZMZMZ

    call message('    objective_kge_q_sm_corr = ', num2str(objective_kge_q_sm_corr, '(F9.5)'))

  END FUNCTION objective_kge_q_sm_corr


  ! -----------------------------------------------------------------

  !    NAME
  !        objective_kge_q_et

  !    PURPOSE
  !>       \brief Objective function of KGE for runoff and KGE for ET

  !>       \details Objective function of KGE for runoff and KGE for ET.
  !>       Further details can be found in the documentation of objective functions
  !>       '14 - objective_multiple_gauges_kge_power6'.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp) :: objective_kge_q_et &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !    HISTORY
  !>       \authors Johannes Brenner

  !>       \date July 2017

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION objective_kge_q_et(parameterset, eval)

    use mo_common_variables, only : level1, nBasins
    use mo_errormeasures, only : kge
    use mo_global_variables, only : L1_et, L1_et_mask
    use mo_message, only : message
    use mo_string_utils, only : num2str
#ifdef MRM2MHM
    use mo_mrm_objective_function_runoff, only : extract_runoff
#endif

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), POINTER :: eval

    real(dp) :: objective_kge_q_et

    real(dp) :: objective_et

    real(dp) :: objective_q

    ! number of invalid cells in catchment
    real(dp) :: invalid_cells

    ! modelled runoff for a given parameter set
    ! dim1=nTimeSteps, dim2=nGauges
    real(dp), allocatable, dimension(:, :) :: runoff

    ! basin loop counter
    integer(i4) :: iBasin

    ! cell loop counter
    integer(i4) :: iCell

    ! start and end index for the current basin
    integer(i4) :: s1, e1

    ! ncells1 of level 1
    integer(i4) :: nCells1

    ! basins wise objectives
    real(dp) :: objective_et_basin

    ! simulated evapotranspiration
    ! (dim1=ncells, dim2=time)
    real(dp), dimension(:, :), allocatable :: et_opti

    real(dp), parameter :: onesixth = 1.0_dp / 6.0_dp

#ifdef MRM2MHM
    ! gauges counter
    integer(i4) :: gg

    integer(i4) :: nGaugesTotal

    ! aggregated simulated runoff
    real(dp), dimension(:), allocatable :: runoff_agg

    ! measured runoff
    real(dp), dimension(:), allocatable :: runoff_obs

    ! mask for measured runoff
    logical, dimension(:), allocatable :: runoff_obs_mask
#endif

    ! run mHM
    call eval(parameterset, runoff = runoff, et_opti = et_opti)

    ! -----------------------------
    ! EVAPOTRANSPIRATION
    ! -----------------------------

    ! initialize some variables
    objective_et = 0.0_dp

    ! loop over basin - for applying power law later on
    do iBasin = 1, nBasins

      ! init
      objective_et_basin = 0.0_dp
      ! get basin information
      nCells1 = level1(iBasin)%nCells
      s1 = level1(iBasin)%iStart
      e1 = level1(iBasin)%iEnd


      ! correlation signal is calculated on individual grid cells
      invalid_cells = 0.0_dp
      do iCell = s1, e1

        ! check for enough data points in time for statistical calculations (e.g. mean, stddev)
        if (count(L1_et_mask(iCell, :)) .LE. (0.10_dp * real(size(L1_et, dim = 2), dp))) then
          invalid_cells = invalid_cells + 1.0_dp
          cycle
        end if

        ! calculate ojective function
        objective_et_basin = objective_et_basin + &
                kge(L1_et(iCell, :), et_opti(iCell, :), mask = L1_et_mask(iCell, :))
      end do

      ! user information about invalid cells
      if (invalid_cells .GT. 0.5_dp) then
        call message('   WARNING: objective_et: Detected invalid cells in study area (.LT. 10 valid data points).')
        call message('                          Fraction of invlauid cells: ', &
                num2str(invalid_cells / real(nCells1, dp), '(F4.2)'))
      end if

      ! calculate average soil moisture objective over all basins with power law
      ! basins are weighted equally ( 1 / real(nBasin,dp))**6
      objective_et = objective_et + &
              ((1.0_dp - objective_et_basin / real(nCells1, dp)) / real(nBasins, dp))**6
    end do

    ! compromise solution - sixth root
    objective_et = objective_et**onesixth

    ! -----------------------------
    ! RUNOFF
    ! -----------------------------
    objective_q = 0.0_dp
#ifdef MRM2MHM
    nGaugesTotal = size(runoff, dim = 2)

    do gg = 1, nGaugesTotal

      ! extract runoff
      call extract_runoff(gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask)

      ! KGE
      objective_q = objective_q + &
              ((1.0_dp - kge(runoff_obs, runoff_agg, mask = runoff_obs_mask)) / real(nGaugesTotal, dp))**6

    end do

    deallocate(runoff_agg, runoff_obs, runoff_obs_mask)

    ! compromise solution - sixth root
    objective_q = objective_q**onesixth

#else
    call message('***ERROR: objective_kge_q_et: missing routing module for optimization')
    stop
#endif

    ! equal weighted compromise objective functions for discharge and soilmoisture
    objective_kge_q_et = (objective_et**6 + objective_q**6)**onesixth
    !    print*, "1-SM 2-Q : ", 1.0_dp-objective_sm, 1.0_dp-objective_kge ! MZMZMZMZ

    call message('    objective_kge_q_et = ', num2str(objective_kge_q_et, '(F9.5)'))

  END FUNCTION objective_kge_q_et

  ! -----------------------------------------------------------------

  !    NAME
  !        objective_kge_q_rmse_et

  !    PURPOSE
  !>       \brief Objective function of KGE for runoff and RMSE for basin_avg ET (standarized scores)

  !>       \details Objective function of KGE for runoff and RMSE for basin_avg ET (standarized scores)

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp) :: objective_kge_q_rmse_et &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !    HISTORY
  !>       \authors Johannes Brenner

  !>       \date July 2017

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION objective_kge_q_rmse_et(parameterset, eval)

    use mo_common_constants, only : eps_dp, nodata_dp
    use mo_common_mhm_mrm_variables, only : evalPer
    use mo_common_variables, only : level1, nBasins
    use mo_errormeasures, only : rmse
    use mo_global_variables, only : L1_et, L1_et_mask, timeStep_et_input
    use mo_julian, only : caldat
    use mo_message, only : message
    use mo_moment, only : average, mean
    use mo_standard_score, only : classified_standard_score
    use mo_string_utils, only : num2str
    use mo_temporal_aggregation, only : day2mon_average
#ifdef MRM2MHM
    use mo_errormeasures, only : kge
    use mo_mrm_objective_function_runoff, only : extract_runoff
#endif

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), POINTER :: eval

    real(dp) :: objective_kge_q_rmse_et

    ! modelled runoff for a given parameter set
    ! dim1=nTimeSteps, dim2=nGauges
    real(dp), allocatable, dimension(:, :) :: runoff

    ! modelled runoff for a given parameter set
    ! dim1=nTimeSteps, dim2=nGauges
    real(dp), allocatable, dimension(:, :) :: et_opti

    ! time loop counter
    integer(i4) :: iTime

    ! start and end index for the current basin
    integer(i4) :: s1, e1

    ! basin counter, month counters
    integer(i4) :: iBasin, pp, mmm

    integer(i4) :: year, month, day

    real(dp), dimension(nBasins) :: initTime

    ! total number of months
    integer(i4) :: nMonths

    ! vector with months' classes
    integer(i4), dimension(:), allocatable :: month_classes

    ! monthly values original time series
    real(dp), dimension(:), allocatable :: et_sim_m, et_obs_m

    ! monthly values anomaly time series
    real(dp), dimension(:), allocatable :: et_sim_m_anom, et_obs_m_anom

    logical, dimension(:), allocatable :: et_obs_m_mask

    ! kge_q(nGaugesTotal)
    real(dp), dimension(:), allocatable :: rmse_et

    ! obj. functions
    real(dp) :: rmse_et_avg, kge_q_avg

    ! spatial average of observed et
    real(dp), dimension(:), allocatable :: et_catch_avg_basin

    ! spatial avergae of modeled  et
    real(dp), dimension(:), allocatable :: et_opti_catch_avg_basin

    ! simulated et
    ! (dim1=ncells, dim2=time)
    !real(dp), dimension(:,:), allocatable :: et_opti

    ! mask for valid et catchment avg time steps
    logical, dimension(:), allocatable :: mask_times

#ifdef MRM2MHM
    ! rmse_et(nBasins)
    real(dp), dimension(:), allocatable :: kge_q

    ! gauges counter
    integer(i4) :: gg

    integer(i4) :: nGaugesTotal

    ! aggregated simulated runoff
    real(dp), dimension(:), allocatable :: runoff_agg

    ! measured runoff
    real(dp), dimension(:), allocatable :: runoff_obs

    ! mask for measured runoff
    logical, dimension(:), allocatable :: runoff_obs_mask
#endif

    ! obtain simulation values of runoff (hourly) and ET
    ! for ET only valid cells (basins concatenated)
    ! et_opti: aggregate ET to needed time step for optimization (see timeStep_et_input)
    call eval(parameterset, runoff = runoff, et_opti = et_opti)

    !--------------------------------------------
    !! EVAPOTRANSPIRATION
    !--------------------------------------------

    ! allocate
    allocate(rmse_et(nBasins))
    rmse_et(:) = nodata_dp

    do iBasin = 1, nBasins
      ! get basin info
      s1 = level1(iBasin)%iStart
      e1 = level1(iBasin)%iEnd

      ! allocate
      allocate(mask_times             (size(et_opti, dim = 2)))
      allocate(et_catch_avg_basin     (size(et_opti, dim = 2)))
      allocate(et_opti_catch_avg_basin(size(et_opti, dim = 2)))

      ! initalize
      mask_times = .TRUE.
      et_catch_avg_basin = nodata_dp
      et_opti_catch_avg_basin = nodata_dp

      ! calculate catchment average evapotranspiration
      do iTime = 1, size(et_opti, dim = 2)
        ! check for enough data points in time for correlation
        if (all(.NOT. L1_et_mask(s1 : e1, iTime))) then
          !write (*,*) 'WARNING: et data at time ', iTime, ' is empty.'
          !call message('WARNING: objective_et_kge_catchment_avg: ignored currrent time step since less than')
          !call message('         10 valid cells available in evapotranspiration observation')
          mask_times(iTime) = .FALSE.
          cycle
        end if
        ! spatial average of observed ET
        et_catch_avg_basin(iTime) = average(L1_et(s1 : e1, iTime), mask = L1_et_mask(s1 : e1, iTime))
        ! spatial avergae of modeled ET
        et_opti_catch_avg_basin(iTime) = average(et_opti(s1 : e1, iTime), mask = L1_et_mask(s1 : e1, iTime))
      end do

      ! get initial time of the evaluation period
      initTime(iBasin) = real(evalPer(iBasin)%julStart, dp)

      ! get calendar days, months, year
      call caldat(int(initTime(iBasin)), yy = year, mm = month, dd = day)

      ! if evapotranspiration input daily
      select case(timeStep_et_input)
        ! JBJBJB
        ! daily: aggregate to monthly mean
      case(-1)
        ! calculate monthly averages from daily values of the model
        call day2mon_average(et_opti_catch_avg_basin, year, month, day, et_sim_m, misval = nodata_dp)
        ! calculate monthly averages from daily values of the observations
        call day2mon_average(et_catch_avg_basin, year, month, day, et_obs_m, misval = nodata_dp)
        !
        ! monthly: proceed without action
      case(-2)
        ! simulation
        allocate(et_sim_m(size(et_opti, dim = 2)))
        et_sim_m = et_opti_catch_avg_basin
        ! observation
        allocate(et_obs_m(size(et_opti, dim = 2)))
        et_obs_m = et_catch_avg_basin

        ! yearly: ERROR stop program
      case(-3)
        call message('***ERROR: objective_kge_q_rmse_et: time step of evapotranspiration yearly.')
        stop
      end select
      ! remove mean from modelled time series
      et_sim_m(:) = et_sim_m(:) - mean(et_sim_m(:))
      ! remove mean from observed time series
      et_obs_m(:) = et_obs_m(:) - mean(et_obs_m(:))
      ! get number of months for given basin
      nMonths = size(et_obs_m)
      allocate (month_classes(nMonths))
      allocate (et_obs_m_mask(nMonths))
      allocate (et_obs_m_anom(nMonths))
      allocate (et_sim_m_anom(nMonths))

      month_classes(:) = 0
      et_obs_m_mask(:) = .TRUE.
      et_obs_m_anom(:) = nodata_dp
      et_sim_m_anom(:) = nodata_dp

      ! define months' classes
      mmm = month
      do pp = 1, nMonths
        month_classes(pp) = mmm
        if (mmm .LT. 12) then
          mmm = mmm + 1
        else
          mmm = 1
        end if
      end do
      ! double check missing data
      ! define mask for missing data in observations (there are always data for simulations)
      where(abs(et_obs_m - nodata_dp) .lt. eps_dp) et_obs_m_mask = .FALSE.

      ! calculate standard score
      et_obs_m_anom = classified_standard_score(et_obs_m, month_classes, mask = et_obs_m_mask)
      et_sim_m_anom = classified_standard_score(et_sim_m, month_classes, mask = et_obs_m_mask)
      rmse_et(iBasin) = rmse(et_sim_m_anom, et_obs_m_anom, mask = et_obs_m_mask)

      deallocate (month_classes)
      deallocate (et_obs_m)
      deallocate (et_sim_m)
      deallocate (et_obs_m_mask)
      deallocate (et_sim_m_anom)
      deallocate (et_obs_m_anom)
      !end if

    end do

    rmse_et_avg = sum(rmse_et(:), abs(rmse_et - nodata_dp) .gt. eps_dp) / &
            real(count(abs(rmse_et - nodata_dp) .gt. eps_dp), dp)
    deallocate(rmse_et)

    !--------------------------------------------
    !! RUNOFF
    !--------------------------------------------
    kge_q_avg = 0_dp
#ifdef MRM2MHM
    nGaugesTotal = size(runoff, dim = 2)
    allocate(kge_q(nGaugesTotal))
    kge_q(:) = nodata_dp

    do gg = 1, nGaugesTotal

      ! extract runoff
      call extract_runoff(gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask)

      ! check for whether to proceed with this basin or not
      ! potentially 3 years of data
      !pp = count(runoff_agg .ge. 0.0_dp )
      !if (pp .lt.  365*3 ) then
      !    deallocate (runoff_agg, runoff_obs, runoff_obs_mask)
      !    cycle
      ! else
      ! calculate KGE for each basin:
      kge_q(gg) = kge(runoff_obs, runoff_agg, mask = runoff_obs_mask)
      deallocate (runoff_agg, runoff_obs, runoff_obs_mask)
      ! end if

    end do

    ! calculate average KGE value for runoff
    kge_q_avg = sum(kge_q(:), abs(kge_q - nodata_dp) .gt. eps_dp) / &
            real(count(abs(kge_q - nodata_dp) .gt. eps_dp), dp)
    deallocate(kge_q)
#else
    call message('***ERROR: objective_kge_q_rmse_et: missing routing module for optimization')
    stop
#endif

    !
    objective_kge_q_rmse_et = rmse_et_avg * (1._dp - kge_q_avg)

    call message('    objective_kge_q_rmse_et = ', num2str(objective_kge_q_rmse_et, '(F9.5)'))

  END FUNCTION objective_kge_q_rmse_et

  ! ------------------------------------------------------------------

  !    NAME
  !        extract_basin_avg_tws

  !    PURPOSE
  !>       \brief extracts basin average tws data from global variables

  !>       \details extracts simulated and measured basin average tws from global variables,
  !>       such that they overlay exactly. For measured tws, only the tws
  !>       during the evaluation period are cut, not succeeding nodata values.
  !>       For simulated tws, warming days as well as succeeding nodata values
  !>       are neglected.
  !>       see use in this module above

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: basinId"           current basin Id
  !>       \param[in] "real(dp), dimension(:, :) :: tws" simulated basin average tws

  !    INTENT(OUT)
  !>       \param[out] "real(dp), dimension(:) :: tws_sim"     aggregated simulated
  !>       \param[out] "real(dp), dimension(:) :: tws_obs"     extracted measured
  !>       \param[out] "logical, dimension(:) :: tws_obs_mask" mask of no data values

  !    HISTORY
  !>       \authors Stephan Thober

  !>       \date Oct 2015

  ! Modifications:
  ! Stephan Thober Oct 2015 - moved subroutine to objective_function_sm
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine extract_basin_avg_tws(basinId, tws, tws_sim, tws_obs, tws_obs_mask)

    use mo_common_constants, only : eps_dp, nodata_dp
    use mo_common_mhm_mrm_variables, only : evalPer, nTstepDay, warmingDays
    use mo_global_variables, only : basin_avg_TWS_obs, nMeasPerDay_TWS
    use mo_message, only : message

    implicit none

    ! current basin Id
    integer(i4), intent(in) :: basinId

    ! simulated basin average tws
    real(dp), dimension(:, :), intent(in) :: tws

    ! aggregated simulated
    real(dp), dimension(:), allocatable, intent(out) :: tws_sim

    ! extracted measured
    real(dp), dimension(:), allocatable, intent(out) :: tws_obs

    ! mask of no data values
    logical, dimension(:), allocatable, intent(out) :: tws_obs_mask

    ! basin id
    integer(i4) :: iBasin

    ! timestep counter
    integer(i4) :: tt

    ! length of extracted time series
    integer(i4) :: length

    ! between simulated and measured time scale
    integer(i4) :: factor

    ! simulated Timesteps per Day
    integer(i4) :: TPD_sim

    ! observed Timesteps per Day
    integer(i4) :: TPD_obs

    real(dp), dimension(:), allocatable :: dummy


    ! copy time resolution to local variables
    TPD_sim = nTstepDay
    TPD_obs = nMeasPerDay_TWS

    ! check if modelled timestep is an integer multiple of measured timesteps
    if (modulo(TPD_sim, TPD_obs) .eq. 0) then
      factor = TPD_sim / TPD_obs
    else
      call message(' Error: Number of modelled datapoints is no multiple of measured datapoints per day')
      stop
    end if

    ! extract basin Id
    iBasin = basin_avg_TWS_obs%basinId(basinId)

    ! get length of evaluation period times TPD_obs
    length = (evalPer(iBasin)%julEnd - evalPer(iBasin)%julStart + 1) * TPD_obs

    ! extract measurements
    if (allocated(tws_obs)) deallocate(tws_obs)
    allocate(tws_obs(length))
    tws_obs = basin_avg_TWS_obs%TWS(1 : length, basinId)

    ! create mask of observed tws
    if (allocated(tws_obs_mask)) deallocate(tws_obs_mask)
    allocate(tws_obs_mask(length))
    tws_obs_mask = .TRUE.
    where(abs(tws_obs - nodata_dp) .lt. eps_dp) tws_obs_mask = .FALSE.

    ! extract and aggregate simulated tws
    if (allocated(tws_sim)) deallocate(tws_sim)
    allocate(tws_sim(length))
    ! remove warming days
    length = (evalPer(iBasin)%julEnd - evalPer(iBasin)%julStart + 1) * TPD_sim
    allocate(dummy(length))
    dummy = tws(warmingDays(iBasin) * TPD_sim + 1 : warmingDays(iBasin) * TPD_sim + length, basinId)
    ! aggregate tws
    length = (evalPer(iBasin)%julEnd - evalPer(iBasin)%julStart + 1) * TPD_obs
    forall(tt = 1 : length) tws_sim(tt) = sum(dummy((tt - 1) * factor + 1 : tt * factor)) / &
            real(factor, dp)
    ! clean up
    deallocate(dummy)

  end subroutine extract_basin_avg_tws


END MODULE mo_objective_function
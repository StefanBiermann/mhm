MODULE mo_mrm_tools

  ! small subroutines used by mRM and HRD

  ! Written Maren Kaluza January 2019

  USE mo_kind, ONLY : i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: set_input_variables_for_routing

  ! ------------------------------------------------------------------

CONTAINS


  subroutine set_input_variables_for_routing(tt, nTimeSteps, &                                   ! intent in
                processMatrix, timestep, L11_tsRout, HourSecs, &                                 ! intent in, but could "use"
                                               do_rout, tsRoutFactorIn, timestep_rout, &         ! intent out
                                                                L1_total_runoff, InflowGaugeQ, & ! optional in
                                                                RunToRout, InflowDischarge)      ! optional out
  integer(i4), intent(in)  :: tt
  integer(i4), intent(in)  :: nTimeSteps
  integer(i4), intent(in)  :: processMatrix
  integer(i4), intent(in)  :: timestep
  real(dp),    intent(in)  :: L11_tsRout
  real(dp),    intent(in)  :: HourSecs

  logical,     intent(out) :: do_rout
  real(dp),    intent(out) :: tsRoutFactorIn
  integer(i4), intent(out) :: timestep_rout

  real(dp), dimension(:), optional, intent(in)  :: L1_total_runoff ! [m3 TS-1] Generated runoff
  real(dp), dimension(:), optional, intent(in)  :: InflowGaugeQ

  real(dp), dimension(:), optional, intent(out) :: RunToRout
  real(dp), dimension(:), optional, intent(out) :: InflowDischarge
  ! local
  ! factor between routing and hydrological modelling resolution
  real(dp) :: tsRoutFactor

  ! set input variables for routing
  if (processMatrix .eq. 1) then
    ! >>>
    ! >>> original Muskingum routing, executed every time
    ! >>>
    do_rout = .True.
    tsRoutFactorIn = 1._dp
    timestep_rout = timestep
    if (present(RunToRout) .and. present(L1_total_runoff)) then
      RunToRout(:) = L1_total_runoff(:) ! runoff [mm TST-1] mm per timestep
    end if
    if (present(InflowDischarge) .and. present(InflowGaugeQ)) then
      InflowDischarge(:) = InflowGaugeQ(:) ! inflow discharge in [m3 s-1]
    end if
    !
  else if (processMatrix .eq. 2) then
    ! >>>
    ! >>> adaptive timestep
    ! >>>
    do_rout = .False.
    ! calculate factor
    tsRoutFactor = L11_tsRout / (timestep * HourSecs)
    ! print *, 'routing factor: ', tsRoutFactor
    ! prepare routing call
    if (tsRoutFactor .lt. 1._dp) then
      ! ----------------------------------------------------------------
      ! routing timesteps are shorter than hydrologic time steps
      ! ----------------------------------------------------------------
      ! set all input variables
      do_rout = .True.
      timestep_rout = timestep
      tsRoutFactorIn = tsRoutFactor

      if (present(RunToRout) .and. present(L1_total_runoff)) then
        RunToRout(:) = L1_total_runoff(:) ! runoff [mm TST-1] mm per timestep
      end if
      if (present(InflowDischarge) .and. present(InflowGaugeQ)) then
        InflowDischarge(:) = InflowGaugeQ(:) ! inflow discharge in [m3 s-1]
      end if
    else
      ! ----------------------------------------------------------------
      ! routing timesteps are longer than hydrologic time steps
      ! ----------------------------------------------------------------
      ! set all input variables
      tsRoutFactorIn = tsRoutFactor
      if (present(RunToRout) .and. present(L1_total_runoff)) then
        RunToRout(:) = RunToRout(:) + L1_total_runoff(:)
      end if
      if (present(InflowDischarge) .and. present(InflowGaugeQ)) then
        InflowDischarge(:) = InflowDischarge(:) + InflowGaugeQ(:)
      end if
      ! reset tsRoutFactorIn if last period did not cover full period
      if ((tt .eq. nTimeSteps) .and. (mod(tt, nint(tsRoutFactorIn)) .ne. 0_i4)) &
              tsRoutFactorIn = mod(tt, nint(tsRoutFactorIn))
      if ((mod(tt, nint(tsRoutFactorIn)) .eq. 0_i4) .or. (tt .eq. nTimeSteps)) then
        do_rout = .True.
        timestep_rout = timestep * nint(tsRoutFactor, i4)
        if (present(InflowDischarge) .and. present(InflowGaugeQ)) then
          InflowDischarge(:) = InflowDischarge(:) / tsRoutFactorIn
        end if
      end if
    end if
  end if
  end subroutine set_input_variables_for_routing

END MODULE mo_mrm_tools

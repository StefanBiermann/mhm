!>       \file mo_et_cr.f90

!>       \brief Module for calculating complementary relationship-based ET  [mm s-1]

!>       \details This module calculates ET [mm/s] based the method of complementary relationships (CR)
!>       - Szilagyi et al. (2017) formulation of Brutsaert’s nonlinear CR model (2015).

!>       \authors Johannes Brenner

!>       \date Jul 2019

! Modifications:

MODULE mo_et_cr

  ! This module is for the UFZ CHS mesoscale hydrologic model mHM.

  USE mo_kind, ONLY : i4, dp

  IMPLICIT NONE

  PRIVATE :: atm_pres ! atmospheric pressure (Pa) as a function of elevation

  PUBLIC :: et_cr ! nonlinear CR model


  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        et_cr

  !    PURPOSE
  !>       \brief Szilagyi et al. (2017) formulation of Brutsaert’s nonlinear CR model (2015)

  !>       \details Calculates the Actual Evapotranspiration  \f$ET [mm\;d^{-1}] \f$
  !>       Szilagyi et al. (2017) formulation of Brutsaert’s nonlinear CR model (2015)
  !>       model for a given cell by applying the equation
  !>       \f[ ET = (2 - \frac{E^{max}_p - E_p}{E^{max}_p - E_w} \frac{E_w}{E_p}) \left( \frac{E^{max}_p - E_p}{E^{max}_p - E_w}  \frac{E_w}{E_p} \right)^2 E_p \f]
  !>       where \f$ E_p\;[W\;m^{-2}]\f$ is the potential evapotranspiration,
  !>       i.e., the evapotranspiration rate of a small wet patch in a drying
  !>       (i.e., not fully wet) environment, typically expressed by the Penman (1948) equation,
  !>       \f$ E_p^{max} \;[W\;m^{-2}]\f$ is the maximum value that \f$E_p \f$ can,
  !>       in theory, reach during a complete dry-out of the land surface,
  !>       \f$ E_w \;[W\;m^{-2}]\f$ is the wet-environment ET rate,
  !>       observed over a regionally extensive well-watered surface,
  !>       specified by the Priestley-Taylor equation (Priestley & Taylor, 1972).
  !>       \note Ma, N., & Szilagyi, J. (n.d.). The complementary relationship ( CR ) of evaporation :
  !>       a calibration- free diagnostic and benchmarking tool for large-scale
  !>       terrestrial evapotranspiration modeling. Water Resources Research, in Review.

  !    INTENT(IN)
  !>       \param[in] "real(dp) :: net_rad"            ! net radiation \f$[W m^{-2}]\f$
  !>       \param[in] "real(dp) :: tavg"               ! average daily temperature \f$[^{\circ}C]\f$
  !>       \param[in] "real(dp) :: windspeed"          ! average daily wind speed at 2m hight \f$[m s^{−1}]\f$
  !>       \param[in] "real(dp) :: act_vap_pressure"   ! actual vapor pressure \f$[kPa]\f$

  !    RETURN
  !>       \return real(dp) :: et_cr &mdash; daily actual evapotranspiration [mm s-1]

  !    HISTORY
  !>       \authors Johannes Brenner

  !>       \date Jul 2019

  ! Modifications:

  function et_cr(net_rad, tavg, windspeed, act_vap_pressure, elevation, PrieTayParam)

    use mo_pet,        only : pet_penman48, pet_priestly, sat_vap_pressure
    use mo_psychrolib, only : GetRelHumFromVapPres, GetTWetBulbFromRelHum
    use mo_constants,  only : Psychro_dp

    implicit none

    ! net radiation \f$[W m^{-2}]\f$
    real(dp), intent(in) :: net_rad

    ! average daily temperature \f$[^{\circ}C]\f$
    real(dp), intent(in) :: tavg

    ! average daily wind speed at 2m hight \f$[m s^{−1}]\f$
    real(dp), intent(in) :: windspeed

    ! actual vapor pressure \f$[kPa]\f$
    real(dp), intent(in) :: act_vap_pressure

    ! elevation a.s.l \f$[m]\f$
    real(dp), intent(in) :: elevation

    ! Priestley-Taylor coefficient \f$ \alpha [-] \f$
    real(dp), intent(in) :: PrieTayParam

    ! daily potential evapotranspiration \f$[W\;m^{-2}]\f$
    real(dp) :: etp

    ! relative humidity in range [0, 1]
    real(dp) :: relhum

    ! dry-environment air temperature \f$[^{\circ}C]\f$
    real(dp) :: tdry

    ! wet-bulb temperature \f$[^{\circ}C]\f$
    real(dp) :: twetbulb

    ! daily maximum of potential evapotranspiration \f$[W\;m^{-2}]\f$
    real(dp) :: etpmax

    ! daily wet patch evapotranspiration \f$[W\;m^{-2}]\f$
    real(dp) :: ew

    ! dimensionless evapotranspiration term X \f$[-]\f$
    real(dp) :: X

    ! actual evapotranspiration in [mm s-1]
    real(dp) :: et_cr

    ! calculate potential evapotranspiration, Penman (1948) approach
    etp = pet_penman48(net_rad, tavg, windspeed, act_vap_pressure)

    ! calculate relative humidity from actual vapor pressure
    relhum = GetRelHumFromVapPres(TDryBulb = tavg, VapPres = act_vap_pressure)

    ! calculate wet bulb temperature from relative humidity
    twetbulb = GetTWetBulbFromRelHum(TDryBulb = tavg, RelHum = relhum, Pressure = atm_pres(elevation))

    ! calculate dry-environment air temperature
    tdry = twetbulb + sat_vap_pressure(twetbulb) / Psychro_dp

    ! calculate daily maximum of potential evapotranspiration
    etpmax = pet_penman48(net_rad, tdry, windspeed, act_vap_pressure = 0.0_dp)

    ! calculate daily wet patch evapotranspiration
    ew = pet_priestly(PrieTayParam, net_rad, tavg)
    ! \Delta in P.T. formulation should be evaluated at the air temperature, Tw (°C),
    ! observed in a wet environment, instead of the typical, drying environment T
    ! (Szilagyi, 2014; Szilagyi & Jozsa, 2008)

    ! calculate actual evapotranspiration
    X = (etpmax - etp) / (etpmax - ew) * (ew / etp)
    et_cr = (2 - X) * X * X * etp

  end function et_cr

  ! ------------------------------------------------------------------------------
  !    NAME
  !        atm_pres
  !
  !    PURPOSE
  !        calculates the atmospheric pressure as a function of elevation

  !    CALLING SEQUENCE
  !        sm_dtr = atm_pres(elev)

  !    INTENT(IN)
  !        real(dp) :: elev - elevation of site (m a.s.l.)
  !

  !    RESTRICTIONS
  !

  !    EXAMPLE
  !        none

  !    LITERATURE
  !        Iribane, J.V., and W.L. Godson, 1981. Atmospheric Thermodynamics, 2nd
  !		Edition. D. Reidel Publishing Company, Dordrecht, The Netherlands.
  !		(p. 168)
  !    HISTORY
  !        Written,  Johannes Brenner, Sept 2016

function atm_pres(elev)

   use mo_constants, only: Gravity_dp, LRstd_dp, R_dp, MA_dp, Tstd_dp, P0_dp

   implicit none
   !
   ! daily atmospheric pressure (Pa) as a function of elevation (m)
   real(dp), intent(in)   :: elev
   real(dp)               :: atm_pres
   real(dp)               :: t1, t2
   !
   ! (-K m-1) standard temperature lapse rate
   !real(dp)               :: LR_STD = 0.0065
   ! (K) standard temp at 0.0 m elevation
   !real(dp)               :: T_STD = 288.15
   ! (m s-2) standard gravitational accel. - Gravity_dp
   !real(dp)               :: G_STD = 9.80665
   ! (Pa) standard pressure at 0.0 m elevation
   !real(dp)               :: P_STD = 101325.0 - P0_dp
   ! (m3 Pa mol-1 K-1) gas law constant
   !real(dp)               :: R = 8.3143
   ! (kg mol-1) molecular weight of air
   !real(dp)               :: MA = 28.9644e-3

   !
   t1 = 1.0 - (LRstd_dp * elev)/Tstd_dp
   t2 = Gravity_dp / (LRstd_dp * (R_dp / MA_dp))
   atm_pres = P0_dp * (t1 **t2)
   !
end function atm_pres

END MODULE mo_et_cr

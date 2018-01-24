!> \file mo_pet.f90

!> \brief Module for calculating reference/potential evapotranspiration  [mm s-1]

!> \details This module calculates PET [mm/s] based on one of the methods \n
!>          - Hargreaves-Samani (1982) \n
!>          - Priestly-Taylor (1972) \n
!>          - Penman-Monteith FAO (1998) \n

!> \author Matthias Zink, Christoph Schneider, Matthias Cuntz
!> \date   Apr 2014

MODULE mo_pet

  ! This module is for the UFZ CHS mesoscale hydrologic model mHM.

  USE mo_kind,      ONLY: i4, dp

  IMPLICIT NONE

  PRIVATE :: extraterr_rad_approx
  PRIVATE :: slope_satpressure
  PRIVATE :: sat_vap_pressure
  PRIVATE :: nusselt_number
  PRIVATE :: h_c
  PRIVATE :: g_bw
  PRIVATE :: P_a
  PRIVATE :: T_l
  PRIVATE :: E_l
  PRIVATE :: H_l
  PRIVATE :: R_ll

  PUBLIC :: pet_hargreaves ! Hargreaves-Samani
  PUBLIC :: pet_priestly   ! Priestley-Taylor
  PUBLIC :: pet_penman     ! Penman-Monteith


  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         pet_hargreaves

  !>        \brief Reference Evapotranspiration after Hargreaves

  !>        \details Calculates the Reference Evapotranspiration \f$ [mm\;d^{-1}] \f$ based on the Hargreaves-Samani (1982)
  !>                 model for a given cell by applying the equation
  !>                 \f[ PET = HarSamCoeff * R_a * (T_{avg} +  HarSamConst) * \sqrt{ T_{max} - T_{min}} \f]
  !>                 where \f$ R_a\;[W\;m^{-2}]\f$ is the incoming solar radiation and
  !>                 \f$ T_{avg}, T_{max} \f$ and \f$ T_{min}\f$  \f$ [ ^{\circ}C]\f$ are the mean, maximum,
  !>                 and minimum daily temperatures at the given day, respectively.

  !     INTENT(IN)
  !>        \param[in] "real(dp),    intent(in) :: HarSamCoeff" coefficient of Hargreaves-Samani equation [-]
  !>        \param[in] "real(dp),    intent(in) :: HarSamConst" constant    of Hargreaves-Samani equation [-]
  !>        \param[in] "real(dp),    intent(in) :: tavg"        daily men temperature \f$[^{\circ}C]\f$
  !>        \param[in] "real(dp),    intent(in) :: tmax"        daily maximum of temp \f$[^{\circ}C]\f$
  !>        \param[in] "real(dp),    intent(in) :: tmin"        daily minimum of temp \f$[^{\circ}C]\f$
  !>        \param[in] "real(dp),    intent(in) :: latitude"   latitude of the cell for Ra estimation \f$[radians]\f$

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
  !>        \return real(dp) :: pet_hargreaves &mdash; Hargreaves-Samani pot. evapotranspiration [mm s-1]

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         \note Hargreaves, G.H., and Samani, Z.A. (1982). "Estimating potential evapotranspiration."
  !             Tech. Note, J. Irrig. and drain. Engrg., ASCE, 108(3):225-230.

  !     HISTORY
  !>        \author   Matthias Zink
  !>        \date     Dec 2012

  elemental pure FUNCTION pet_hargreaves(HarSamCoeff, HarSamConst, tavg, tmax, tmin, latitude, doy)

    use mo_constants,     only: deg2rad_dp
    use mo_utils,         only: LE

    implicit none

    real(dp),    intent(in) :: HarSamCoeff           !  coefficient of Hargreaves-Samani equation
    real(dp),    intent(in) :: HarSamConst           !  constatnt   of Hargreaves-Samani equation
    real(dp),    intent(in) :: tavg                  !  daily men temperature
    real(dp),    intent(in) :: tmax                  !  daily maximum of temp.
    real(dp),    intent(in) :: tmin                  !  daily minimum of temp.
    real(dp),    intent(in) :: latitude              ! latitude of the cell for Ra estimation
    integer(i4), intent(in) :: doy                   ! day of year for Ra estimation

    real(dp)                :: pet_hargreaves        ! reference evapotranspiration in [mm s-1]

    ! local
    real(dp)                :: delta_temp            ! tmax-Tmin

    ! correction for shity input data (tmax<tmin) and to avoid numerical errors ! MZMZMZMZ
    delta_temp = tmax - tmin
    if( LE(delta_temp, 0.0_dp) .or. LE(tavg, -HarSamConst) ) then
       pet_hargreaves = 0.0_dp
    else
       pet_hargreaves = HarSamCoeff * extraterr_rad_approx(doy, deg2rad_dp * latitude) * (tavg + HarSamConst) * sqrt(delta_temp)
    end if

  END FUNCTION pet_hargreaves


  ! ------------------------------------------------------------------

  !     NAME
  !         pet_priestly

  !>        \brief Reference Evapotranspiration after Priestly-Taylor

  !>        \details Calculates the Reference Evapotranspiration \f$ [mm\;d^{-1}] \f$ based on the
  !>                 Priestly-Taylor (1972) model for every given cell by applying the equation
  !>                 \f[ PET = \alpha * \frac{\Delta}{(\gamma + \Delta)} * R_n \f]
  !>                 where \f$R_n\;[W\;m^{-2}]\f$ is the net solar radiation \f$\Delta =  f(T_{avg})\f$ is the slope
  !>                 of the saturation-vapour pressure curve and \f$\alpha\f$ is a emperical coefficient.

  !     INTENT(IN)
  !>        \param[in] "real(dp) :: PrieTayParam" Priestley-Taylor coefficient \f$ \alpha [-] \f$
  !>        \param[in] "real(dp) :: Rn"           net solar radiation \f$ [W\;m^{-2}] \f$
  !>        \param[in] "real(dp) :: Tavg"         daily mean air temperature \f$ [ ^{\circ}C]\f$

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
  !>         \return real(dp) :: pet_priestly &mdash; Priestley-Taylor pot. evapotranspiration [mm s-1]

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !>         \note Priestley, C.H.B., and R.J. Taylor. 1972. On the assessment of surface heat flux and evaporation using
  !>                   large-scale parameters. Mon. Weather Rev., 100:81-82.
  !>         \note ASAE Standards. 1998. EP406.2: heating, cooling, and ventilating greenhouses. St. Joseph, MI, USA.

  !     HISTORY
  !>        \author  Matthias Zink
  !>        \date    Apr 2014

  elemental pure FUNCTION pet_priestly(PrieTayParam, Rn, tavg)
!  FUNCTION pet_priestly(PrieTayParam, Rn, tavg)

    use mo_mhm_constants, only: DaySecs
    use mo_constants,     only: Psychro_dp, SpecHeatET_dp

    implicit none

    real(dp), intent(in) :: PrieTayParam       ! Priestley-Taylor coefficient
    real(dp), intent(in) :: Rn
    real(dp), intent(in) :: Tavg
    real(dp)             :: pet_priestly       ! reference evapotranspiration in [mm s-1]

    real(dp)             :: delta              ! save slope of saturation vapor pressure curve

    delta        = slope_satpressure(Tavg) ! slope of saturation vapor pressure curve
    ! in [mm d-1]
    pet_priestly = PrieTayParam * delta / (Psychro_dp + delta) * ( Rn * DaySecs / SpecHeatET_dp )

  END FUNCTION pet_priestly



  ! ------------------------------------------------------------------

  !     NAME
  !         pet_penman

  !>        \brief Reference Evapotranspiration after Penman-Monteith

  !>        \details Calculates the reference evapotranspiration \f$ [mm\;d^{-1}] \f$ based on the
  !>                 Penman-Monteith model for every given cell by applying the equation
  !>                 \f[ PET = \frac{1}{\lambda}  \cdot
  !>                           \frac{\Delta \cdot R_n + \rho \cdot c_p \cdot (e_s-e) \cdot \frac{a_sh}{r_a}}
  !>                           {\Delta + \gamma \cdot \frac{a_sh}{a_s} \cdot \left( 1 + \frac{r_s}{r_a} \right) }         \f]
  !>                 where \f$R_n\;[W\;m^{-2}]\f$ is the net solar radiation,
  !>                 \f$\Delta\;[kPa\;K^{-1}]\f$ is the slope of the saturation-vapour pressure curve,
  !>                 \f$ \lambda\;[MJ\;kg^{-1}] \f$ is the latent heat of vaporization,
  !>                 \f$ (e_s-e)\;[kPa] \f$ is the vapour pressure deficit of the air,
  !>                 \f$ \rho\;[kg\;m^{-3}] \f$ is the mean atmospheric density,
  !>                 \f$ c_p=1005.0\;J\;kg^{-1}\;K^{-1} \f$ is the specific heat of the air,
  !>                 \f$ \gamma [kPa\;K^{-1}] \f$ is the psychrometric constant,
  !>                 \f$ r_s [s m^{-1}] \f$ is the bulk canopy resistance,
  !>                 \f$ r_a [s m^{-1}] \f$ is the aerodynamic resistance,
  !>                 \f$ a_s [1] \f$ is the fraction of one-sided leaf area covered by stomata
  !>                 (1 if stomata are on one side only, 2 if they are on both sides) and
  !>                 \f$ a_{sh} [-] \f$ is the fraction of projected area exchanging sensible heat with the air (2)
  !>                 Implementation refers to the so-called Penman-Montheith equation for transpiration.
  !>                 Adjusting the arguments \f$ a_{sh} \f$ and \f$ a_s \f$ we obtain the corrected MU equation (for details
  !>                 see Schymanski and Or, 2017). If \f$ a_{sh} = 1 = a_s \f$ Penman-Montheith equation for transpiration
  !>                 is preserved. For reproducing characteristics of symmetrical amphistomatous leaves use
  !>                 \f$ a_{sh} = 2 = a_s \f$, in which case the classic PM equation is only missing a factor
  !>                 of 2 in the nominator, as pointed out by Jarvis and McNaughton (1986, Eq. A9).
  !>                 These analytical solutions eliminated the non-linearity problem of the saturation vapour pressure curve,
  !>                 but they do not consider the dependency of the long-wave component of the soil surface or leaf energy balance
  !>                 (\f$ R_l \f$) on soil or leaf temperature (\f$ T_l \f$). We assume that net radiation
  !>                 equals the absorbed short-wave radiation, i.e. \f$ R_N = R_s \f$ (p.79 in Monteith and Unsworth, 2013).

  !     INTENT(IN)
  !>        \param[in] "real(dp), intent(in) :: net_rad"                net solar radiation \f$[W m^{-2}]\f$
  !>        \param[in] "real(dp), intent(in) :: lw_rad_up"              long-wave radiation away from soil/leaf \f$[W m^{-2}]\f$
  !>        \param[in] "real(dp), intent(in) :: tavg"                   average daily temperature \f$[^{\circ}C]\f$
  !>        \param[in] "real(dp), intent(in) :: act_vap_pressure"       actual vapur pressure \f$[kPa]\f$
  !>        \param[in] "real(dp), intent(in) :: aerodyn_resistance"     aerodynmaical resistance \f$s\;m^{-1}\f$
  !>        \param[in] "real(dp), intent(in) :: bulksurface_resistance" bulk surface resistance  \f$s\;m^{-1}\f$
  !>        \param[in] "real(dp), intent(in) :: a_s"                    fraction of one-sided leaf area covered by stomata \f$1\f$
  !>        \param[in] "real(dp), intent(in) :: a_sh"                   fraction of projected area exchanging sensible heat with the air \f$1\f$
  !>        \param[in] "real(dp)             :: pet_penman"             reference evapotranspiration \f$[mm\;s^{-1}]\f$

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
  !>         \return real(dp) :: pet_penman &mdash; Reference Evapotranspiration [mm s-1]

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !>        \note Allen, R. G. R., Pereira, L., Raes, D., & Smith, M. (1998). Crop evapotranspiration - Guidelines for
  !>             computing crop water requirements - FAO Irrigation and drainage paper 56. Rome.
  !>        \note Schymanski, S. J., & Or, D. (2017). Leaf-scale experiments reveal an important omission in the
  !>             Penman-Monteith equation. HESS, 21(2), 685-706.
  !>        \note Monteith, J. L. and Unsworth, M. H. (2013) Principles of environmental physics: plants, animals,
  !>             and the atmosphere, 4th Edn., Elsevier/Academic Press, Amsterdam, Boston.

  !     HISTORY
  !>        \author  Matthias Zink
  !>        \date    Apr 2014
  ! Modified,
  ! Johannes Brenner, Nov 2017 - include arguments a_s and a_sh to enable corrected MU approach

  elemental pure FUNCTION pet_penman(net_rad, lw_rad_up, tavg, act_vap_pressure, aerodyn_resistance, bulksurface_resistance, a_s, a_sh)

    use mo_mhm_constants, only: DaySecs
    use mo_constants,     only: Psychro_dp, SpecHeatET_dp, rho0_dp, cp0_dp

    implicit none

    real(dp), intent(in) :: net_rad                ! net solar radiation
    real(dp), intent(in) :: lw_rad_up              ! long-wave radiation away from soil/leaf
    real(dp), intent(in) :: tavg                   ! average daily temperature
    real(dp), intent(in) :: act_vap_pressure       ! actual vapur pressure
    real(dp), intent(in) :: aerodyn_resistance     ! aerodynmaical resistance
    real(dp), intent(in) :: bulksurface_resistance ! bulk surface resistance
    real(dp), intent(in) :: a_s                    ! fraction of one-sided leaf area covered by stomata
    real(dp), intent(in) :: a_sh                   ! fraction of projected area exchanging sensible heat with the air
    real(dp)             :: pet_penman             ! reference evapotranspiration in [mm s-1]

    pet_penman =  DaySecs / SpecHeatET_dp  *           & ! conversion factor [W m-2] to [mm d-1]
                  (slope_satpressure(tavg) * (net_rad - lw_rad_up) + &
                  rho0_dp * cp0_dp * (sat_vap_pressure(tavg) - act_vap_pressure ) * a_sh / aerodyn_resistance) / &
                  (slope_satpressure(tavg) + Psychro_dp * a_sh / a_s * (1.0_dp + bulksurface_resistance/aerodyn_resistance))

  END FUNCTION pet_penman

  ! ------------------------------------------------------------------

  !     NAME
  !         extraterr_rad_approx

  !>        \brief Approximation of extraterrestrial radiation

  !>        \details Approximation of extraterrestrial radiation at the top of the atmosphere \f$ R_a \f$
  !>                 after Duffie and Beckman (1980).
  !>                 \f$ R_a \f$ is converted from \f$ [J\;m^{-2}\;d^{-1}] \f$ in \f$ [mm\;d^{-1}]\f$ .
  !>                 \f[ R_a   = \frac{86400}{ \pi \cdot \lambda} \cdot E_0 \cdot
  !>                  d_r \cdot (\omega \cdot \sin(latitude) \cdot \sin(\delta) + \cos(latitude) \cdot \cos(\delta) \cdot
  !>                 \sin(\omega) \f]
  !>                 where \f$ E_0=1367\;J\;m^{-2}\;s^{-1} \f$ is the solar constant and
  !<                 \f$ \lambda = 2.45 \cdot 10^6\;J\;m^{-2}\;mm^{-1} \f$ is the latent heat of vaporization. \n
  !>                 It is dependent on the following sub equations:\n
  !>                 The relative distance Earth-Sun:
  !>                 \f[ d_r =  1 + 0.033 \cdot \cos \left( \frac{2 \cdot \pi \cdot doy}{365} \right) \f]
  !>                 in which doy is the day of the year.\n
  !>                 The solar declination [radians] defined by
  !>                 \f[ \delta  = 0.4093 \cdot \sin\left( \frac{2 \cdot \pi \cdot doy}{365} - 1.405 \right) \f]
  !>                 The sunset hour angle [radians]:
  !>                 \f[ \omega  = \arccos( - \tan(latitude) * \tan(\delta) )  \f]

  !     INTENT(IN)
  !>        \param[in] "integer(i4), intent(in) :: doy" day of year [-]
  !>        \param[in] "real(dp),    intent(in) :: latitude" latitude [rad]

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
  !>         \return real(dp) :: extraterr_rad_approx &mdash; extraterrestrial radiation approximation \f$[W\;m^{-2}]\f$

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !>        \note Duffie, J.A. and W.A. Beckman. 1980. Solar engineering of thermal processes.
  !>            John Wiley and Sons, New York. pp. 1-109.

  !     HISTORY
  !>        \author   Matthias Zink
  !>        \date     Apr 2014
  !         Modified  R. Kumar and M. Zink,    June 2016 - correction to include NaN in acos(arg)

  elemental pure FUNCTION extraterr_rad_approx(doy, latitude)

    use mo_constants,     only: SolarConst_dp, SpecHeatET_dp, PI_D, TWOPI_D
    use mo_mhm_constants, only: DuffieDr, DuffieDelta1, DuffieDelta2, YearDays, DaySecs

    implicit none

    integer(i4), intent(in)             :: doy
    real(dp),    intent(in)             :: latitude             ! latitude [rad]
    real(dp)                            :: extraterr_rad_approx ! extraterrestrial radiation

    ! local
    real(dp)                            :: dr, delta
    real(dp)                            :: omega
    real(dp)                            :: arg

    ! inverse relative distance Earth-Sun - correction for eccentricity of Earths orbit around the sun
    dr     =  1.0_dp + DuffieDr * cos( TWOPI_D * doy / YearDays )
    ! declination of the sun above the celestial equator in radians
    delta  =       DuffieDelta1 * sin( TWOPI_D * doy / YearDays - DuffieDelta2 )

    ! arccos(x) is only defined between PI and 0 (for x between -1 and 1)
    ! check limits
    arg = - tan(latitude) * tan(delta)
    if( arg .lt. -1.0_dp ) arg = -1.0_dp
    if( arg .gt.  1.0_dp ) arg =  1.0_dp

    ! sunrise hour angle in radians
    omega  = acos( arg )

    ! Ra - converted from [J m-2 d-1] in [mm d-1]
    extraterr_rad_approx   = DaySecs / PI_D / SpecHeatET_dp * SolarConst_dp *  &
         dr * (omega * sin(latitude) * sin(delta) + cos(latitude) * cos(delta) * sin(omega))

  end FUNCTION extraterr_rad_approx


  ! ------------------------------------------------------------------

  !     NAME
  !         slope_satpressure

  !>        \brief slope of saturation vapour pressure curve

  !>        \details slope of saturation vapour pressure curve after Tetens
  !>                 \f[ \Delta = \frac{0.6108 * e_s(T_a)}{e^(2 \cdot \log(T_a + 237.3))} \f]

  !     INTENT(IN)
  !>        \param[in] "real(dp), intent(in) :: tavg" average daily temperature \f$[^{\circ}C]\f$

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
  !>         \return real(dp) :: slope_satpressure &mdash;  slope of saturation vapour pressure curve
  !>                             \f$[kPa\;K{-1}]\f$

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !>         \note Tetens, O., 1930. Ueber einige meteorologische Begriffe. z. Geophys. 6:297-309.
  !>         \note Murray, F.W. 1967. On the computation of saturation vapor pressure. J. Appl. Meteor. 6: 203-204.
  !>         \note Allen, R. G. R., Pereira, L., Raes, D., & Smith, M. (1998). Crop evapotranspiration - Guidelines for
  !>             computing crop water requirements - FAO Irrigation and drainage paper 56. Rome.

  !     HISTORY
  !>        \author  Matthias Zink
  !>        \date    Apr 2014
  !

  elemental pure FUNCTION slope_satpressure(tavg)

    use mo_mhm_constants, only: satpressureslope1, tetens_c3

    implicit none

    real(dp), intent(in) :: tavg                    ! average daily temperature
    real(dp)             :: slope_satpressure       ! slope of saturation vapour pressure curve


    slope_satpressure = satpressureslope1 * sat_vap_pressure(tavg) / exp(2.0_dp*log(Tavg + tetens_c3))

  END FUNCTION slope_satpressure

  ! ------------------------------------------------------------------

  !     NAME
  !         sat_vap_pressure

  !>        \brief calculation of the saturation vapour pressure

  !>        \details Calculation of the saturation vapour pressure
  !>                 \f[ e_s(T_a) = 0.6108 \cdot \exp \left( \frac{17.27 \cdot T_a}{T_a + 237.3} \right)  \f]

  !     INTENT(IN)
  !>        \param[in] "real(dp), intent(in) :: tavg" temperature [degC]

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
  !>        \return real(dp) :: sat_vap_pressure &mdash; saturation vapour pressure [kPa]

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !>         \note Tetens, O., 1930. Ueber einige meteorologische Begriffe. z. Geophys. 6:297-309.
  !>         \note Allen, R. G. R., Pereira, L., Raes, D., & Smith, M. (1998). Crop evapotranspiration - Guidelines for
  !>             computing crop water requirements - FAO Irrigation and drainage paper 56. Rome.

  !     HISTORY
  !>        \author  Matthias Zink
  !>        \date    Apr 2014

  elemental pure FUNCTION sat_vap_pressure(tavg)

    use mo_mhm_constants, only: tetens_c1, tetens_c2, tetens_c3

    implicit none

    real(dp), intent(in) :: tavg                      ! temperature [degC]
    real(dp)             :: sat_vap_pressure          ! saturation vapour pressure [kPa]

    sat_vap_pressure = tetens_c1 * exp(tetens_c2 * tavg / (tavg + tetens_c3))

  END FUNCTION sat_vap_pressure

  ! ------------------------------------------------------------------

  !     NAME
  !         nusselt_number

  !>        \brief calculation of the Nusselt number

  !>        \details Calculation of the Nusselt number
  !>
  !

  !     INTENT(IN)
  !>        \param[in] "real(dp), intent(in) :: tavg" temperature [degC]
  !>        \param[in] "real(dp), intent(in) :: ws" wind speed [m s-1]
  !>        \param[in] "real(dp), intent(in) :: L_l" characteristic length scale of convection (size of leave) [m]

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
  !>        \return real(dp) :: nusselt_number &mdash; Nusselt number [1]

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !>         \note

  !     HISTORY
  !>        \author  Johannes Brenner
  !>        \date    Jan 2018

  elemental pure FUNCTION nusselt_number(tavg, ws, L_l)

    use mo_constants, only: T0_dp
    use mo_mhm_constants, only: prandtl_N, reynolds_Nc

    implicit none

    real(dp), intent(in) :: tavg            ! temperature [degC]
    real(dp), intent(in) :: ws              ! wind speed  [m s-1]
    real(dp), intent(in) :: L_l             ! characteristic length scale of convection (size of leave) [m]
    real(dp)             :: reynolds_N      ! Reynolds number [1]
    real(dp)             :: nusselt_number  ! Nusselt number [1]

    reynolds_N = ws * L_l / ( 9e-8_dp * (tavg + T0_dp) - 1.13e-5_dp )

    nusselt_number = prandtl_N ** (1.0_dp/3.0_dp) * &
         (-0.037_dp * (reynolds_N + reynolds_Nc - 0.5_dp * &
         abs(reynolds_N - reynolds_Nc) ** (4.0_dp/5.0_dp) + &
         0.037_dp * reynolds_N ** (4.0_dp/5.0_dp) + &
         0.664_dp * sqrt(reynolds_N + reynolds_Nc - &
         0.5_dp * abs(reynolds_N - reynolds_Nc))))

  END FUNCTION nusselt_number

  ! ------------------------------------------------------------------

  !     NAME
  !         h_c

  !>        \brief calculation of average one-sided convective transfer coefficient (h_c)

  !>        \details calculation of average one-sided convective transfer coefficient (h_c)
  !>
  !

  !     INTENT(IN)
  !>        \param[in] "real(dp), intent(in) :: tavg" temperature [degC]
  !>        \param[in] "real(dp), intent(in) :: L_l" characteristic length scale of convection (size of leave) [m]

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
  !>        \return real(dp) :: h_c &mdash; average one-sided convective transfer coefficient [J K-1 m-2 s-1]

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !>         \note

  !     HISTORY
  !>        \author  Johannes Brenner
  !>        \date    Jan 2018

  elemental pure FUNCTION h_c(tavg, L_l, nusselt_number)

    use mo_constants, only: T0_dp

    implicit none

    real(dp), intent(in) :: tavg            ! temperature [degC]
    real(dp), intent(in) :: L_l             ! characteristic length scale of convection (size of leave) [m]
    real(dp), intent(in) :: nusselt_number  ! Nusselt number [1]
    real(dp)             :: h_c             ! average one-sided convective transfer coefficient [J K-1 m-2 s-1]

    h_c = 6.84e-5_dp * (tavg + T0_dp) + 5.62e-3_dp * nusselt_number / L_l

  END FUNCTION h_c

  ! ------------------------------------------------------------------

  !     NAME
  !         g_bw

  !>        \brief calculation of boundary layer conductance to water vapour (g_bw)

  !>        \details calculation of boundary layer conductance to water vapour (g_bw)
  !>
  !

  !     INTENT(IN)
  !>        \param[in] "real(dp), intent(in) :: tavg" temperature [degC]
  !>

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
  !>        \return real(dp) :: g_bw &mdash; boundary layer conductance to water vapour [m s-1]

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !>         \note

  !     HISTORY
  !>        \author  Johannes Brenner
  !>        \date    Jan 2018
  !

  elemental pure FUNCTION g_bw(tavg, act_vap_pressure, h_c, a_s)

    use mo_constants, only: T0_dp, MN2_dp, MO2_dp, MH2O_dp, Rmol_dp, TWOTHIRD_dp, cp0_dp

    implicit none

    real(dp), intent(in) :: tavg             ! temperature [degC]
    real(dp), intent(in) :: act_vap_pressure ! actual vapore pressure [kPa]
    real(dp), intent(in) :: h_c              ! average one-sided convective transfer coefficient [J K-1 m-2 s-1]
    real(dp), intent(in) :: a_s              ! fraction of one-sided leaf area covered by stomata [1]
    real(dp)             :: g_bw             ! boundary layer conductance to water vapour [m s-1]
    real(dp)             :: d_wa, alpha_a, rho_a

    d_wa = 1.49e-7_dp * (tavg + T0_dp) + 1.96e-5_dp

    alpha_a = 1.3e-7_dp * (tavg + T0_dp) - 1.73e-5_dp

    rho_a = (79.0_dp * MN2_dp * (sat_vap_pressure(tavg) - act_vap_pressure ) * 100.0_dp + &
             21.0_dp * MO2_dp * (sat_vap_pressure(tavg) - act_vap_pressure ) * 100.0_dp + &
             100.0_dp * MH2O_dp * act_vap_pressure * 100.0_dp) / &
            (100.0_dp * Rmol_dp * (tavg + T0_dp))

    g_bw = a_s * h_c / ((alpha_a / d_wa) ** TWOTHIRD_dp * cp0_dp * rho_a)

  END FUNCTION g_bw

  ! ------------------------------------------------------------------

  !     NAME
  !         P_a

  !>        \brief calculation of air pressure (P_a)

  !>        \details calculation of air pressure (P_a)
  !>
  !

  !     INTENT(IN)
  !>        \param[in] "real(dp), intent(in) :: tavg" temperature [degC]
  !>        \param[in] "real(dp), intent(in) :: h" elevation a.s.l. [m]

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
  !>        \return real(dp) :: P_a &mdash; air pressure [Pa]

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !>         \note

  !     HISTORY
  !>        \author  Johannes Brenner
  !>        \date    Jan 2018
  !
  elemental pure FUNCTION P_a(tavg, h)

    use mo_constants, only: MN2_dp, MO2_dp, MAr_dp, MCO2_dp, P0_dp, Gravity_dp, Rmol_dp, T0_dp

    implicit none

    real(dp), intent(in) :: tavg             ! temperature [degC]
    real(dp), intent(in) :: h                ! height a.s.l. [m]
    real(dp)             :: P_a              ! air pressure [Pa]
    real(dp)             :: M                ! mass of one molecule [kg mol^-1]

    M = 0.7808_dp * MN2_dp + 0.2095_dp * MO2_dp + 0.0093_dp * MAr_dp + 0.0003_dp * MCO2_dp
    P_a = P0_dp * exp(-(M * Gravity_dp) / (Rmol_dp * (tavg + T0_dp)) * h)

  END FUNCTION P_a

  ! ------------------------------------------------------------------

  !     NAME
  !         T_l

  !>        \brief calculation of leaf surace temperature (T_l)

  !>        \details calculation of leaf surace temperature (T_l)
  !>
  !

  !     INTENT(IN)
  !>        \param[in] "real(dp), intent(in) :: Rs" solar short-wave flux [W m^-2]
  !>        \param[in] "real(dp), intent(in) :: tavg" temperature [degC]
  !>        \param[in] "real(dp), intent(in) :: delta_e" slope of saturation vapor pressure curve [kPa K^-1]
  !>        \param[in] "real(dp), intent(in) :: act_vap_pressure"    actual vapur pressure [kPa]
  !>        \param[in] "real(dp), intent(in) :: a_sh"     fraction of projected area exchanging sensible heat with the air \f$1\f$
  !>        \param[in] "real(dp), intent(in) :: g_bw"     boundary layer conductance to water vapour [m s-1]
  !>        \param[in] "real(dp), intent(in) :: g_sw"     stomatal conductance to water vapour [m s-1]
  !>        \param[in] "real(dp), intent(in) :: P_a"      air pressure [Pa]
  !>        \param[in] "real(dp), intent(in) :: h_c"      average one-side convective tranfer coefficient [m s-1]
  !>        \param[in] "real(dp), intent(in) :: eta1"     longwave emissivity of the leaf surface [1]

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
  !>        \return real(dp) :: T_l &mdash; leaf surace temperature [degC]

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !>         \note

  !     HISTORY
  !>        \author  Johannes Brenner
  !>        \date    Jan 2018
  !
  elemental pure FUNCTION T_l(Rs, tavg, delta_e, act_vap_pressure, a_sh, g_bw, g_sw, P_a, h_c, eta1)

    use mo_constants, only: MH2O_dp, sigma_dp, Rmol_dp, SpecHeatET_dp, T0_dp

    implicit none

    real(dp), intent(in) :: Rs               ! solar short-wave flux [W m^-2]
    real(dp), intent(in) :: tavg             ! temperature [degC]
    real(dp), intent(in) :: delta_e          ! slope of saturation vapor pressure curve [kPa K^-1]
    real(dp), intent(in) :: act_vap_pressure ! actual vapur pressure [kPa]
    real(dp), intent(in) :: a_sh             ! fraction of projected area exchanging sensible heat with the air [1]
    real(dp), intent(in) :: g_bw             ! boundary layer conductance to water vapour [m s-1]
    real(dp), intent(in) :: g_sw             ! stomatal conductance to water vapour [m s-1]
    real(dp), intent(in) :: P_a              ! air pressure [Pa]
    real(dp), intent(in) :: h_c              ! average one-sided convective transfer coefficient [J K-1 m-2 s-1]
    real(dp), intent(in) :: eta1             ! long-wave emissivity of the leaf surface [1]
    real(dp)             :: g_tw             ! total leaf layer conductance to water vapour [m s^-1]
    real(dp)             :: g_tw_mol         ! total leaf layer conductance to water vapour [mol m^-2 s^-1]
    real(dp)             :: c_e              ! latent heat tranfer coefficient [J Pa^-1 m^-2 s^-1]
    real(dp)             :: c_h              ! sensible heat tranfer coefficient [J Pa^-1 m^-2 s^-1]
    real(dp)             :: T_l              ! leaf surace temperature [degC]

    ! total leaf layer conductance to water vapour
    g_tw = 1.0_dp / (1.0_dp/g_bw + 1.0_dp/g_sw)
    g_tw_mol = g_tw * P_a / (Rmol_dp * tavg)

    ! latent and sensible heat tranfer coefficient
    c_e = MH2O_dp * SpecHeatET_dp * g_tw_mol / P_a
    c_h = a_sh * h_c

    ! leaf surace temperature
    T_l = (Rs + c_h * (tavg + T0_dp) + c_e * (delta_e * (tavg + T0_dp) +&
     (act_vap_pressure - sat_vap_pressure(tavg)) * 1000.0_dp) +&
     a_sh * eta1 * sigma_dp * (4.0_dp * (tavg + T0_dp)**4)) /&
     (c_h + c_e * delta_e + 4.0_dp * eta1 * sigma_dp * (tavg + T0_dp)**3)

    ! conversion from K in degC
    T_l = T_l - T0_dp

  END FUNCTION T_l

  ! ------------------------------------------------------------------

  !     NAME
  !         E_l

  !>        \brief calculation of latent heat flux from leaf (E_l)

  !>        \details calculation of latent heat flux from leaf (E_l)
  !>
  !

  !     INTENT(IN)
  !>        \param[in] "real(dp), intent(in) :: T_l" leaf surace temperature [degC]
  !>        \param[in] "real(dp), intent(in) :: tavg" temperature [degC]
  !>        \param[in] "real(dp), intent(in) :: delta_e" slope of saturation vapor pressure curve [kPa K^-1]
  !>        \param[in] "real(dp), intent(in) :: act_vap_pressure"    actual vapur pressure [kPa]
  !>        \param[in] "real(dp), intent(in) :: a_sh"     fraction of projected area exchanging sensible heat with the air \f$1\f$
  !>        \param[in] "real(dp), intent(in) :: h_c"     boundary layer conductance to water vapour [m s-1]
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
  !>        \return real(dp) :: E_l &mdash; latent heat flux from leaf [J m^-2 s^-1]

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !>         \note

  !     HISTORY
  !>        \author  Johannes Brenner
  !>        \date    Jan 2018
  !
  elemental pure FUNCTION E_l(T_l, tavg, delta_e, act_vap_pressure, a_sh, h_c)

    use mo_constants, only: Psychro_dp

    implicit none

    real(dp), intent(in) :: T_l              ! leaf surace temperature [degC]
    real(dp), intent(in) :: tavg             ! temperature [degC]
    real(dp), intent(in) :: delta_e          ! slope of saturation vapor pressure curve [kPa K^-1]
    real(dp), intent(in) :: act_vap_pressure ! actual vapur pressure [kPa]
    real(dp), intent(in) :: a_sh             ! fraction of projected area exchanging sensible heat with the air [1]
    real(dp), intent(in) :: h_c              ! average one-sided convective transfer coefficient [J K-1 m-2 s-1]
    real(dp)             :: c_h              ! sensible heat tranfer coefficient [J Pa^-1 m^-2 s^-1]
    real(dp)             :: E_l              ! leaf surace temperature [degC]

    ! sensible heat tranfer coefficient
    c_h = a_sh * h_c

    ! latent heat flux from leaf
    E_l = c_h * (delta_e * (T_l - tavg)) + &
    (sat_vap_pressure(tavg) - act_vap_pressure) * 1000.0_dp / Psychro_dp

  END FUNCTION E_l

  ! ------------------------------------------------------------------

  !     NAME
  !         H_l

  !>        \brief calculation of sensible heat flux from leaf (H_l)

  !>        \details calculation of sensible heat flux from leaf (H_l)
  !>
  !

  !     INTENT(IN)
  !>        \param[in] "real(dp), intent(in) :: T_l"  leaf surace temperature [degC]
  !>        \param[in] "real(dp), intent(in) :: tavg" temperature [degC]
  !>        \param[in] "real(dp), intent(in) :: a_sh" fraction of projected area exchanging sensible heat with the air \f$1\f$
  !>        \param[in] "real(dp), intent(in) :: h_c"  boundary layer conductance to water vapour [m s-1]
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
  !>        \return real(dp) :: H_l &mdash; sensible heat flux from leaf [J m^-2 s^-1]

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !>         \note

  !     HISTORY
  !>        \author  Johannes Brenner
  !>        \date    Jan 2018
  !
  elemental pure FUNCTION H_l(T_l, tavg, a_sh, h_c)

    use mo_constants, only: Psychro_dp

    implicit none

    real(dp), intent(in) :: T_l              ! leaf surace temperature [degC]
    real(dp), intent(in) :: tavg             ! temperature [degC]
    real(dp), intent(in) :: a_sh             ! fraction of projected area exchanging sensible heat with the air [1]
    real(dp), intent(in) :: h_c              ! average one-sided convective transfer coefficient [J K-1 m-2 s-1]
    real(dp)             :: c_h              ! sensible heat tranfer coefficient [J Pa^-1 m^-2 s^-1]
    real(dp)             :: H_l              ! leaf surace temperature [degC]

    ! sensible heat tranfer coefficient
    c_h = a_sh * h_c

    ! latent heat flux from leaf
    H_l = c_h * (T_l - tavg)

  END FUNCTION H_l

  ! ------------------------------------------------------------------

  !     NAME
  !         R_ll

  !>        \brief calculation of long-wave radiation away from leaf (R_ll)

  !>        \details calculation of long-wave radiation away from leaf (R_ll)
  !>
  !

  !     INTENT(IN)
  !>        \param[in] "real(dp), intent(in) :: T_l"  leaf surace temperature [degC]
  !>        \param[in] "real(dp), intent(in) :: tavg" temperature [degC]
  !>        \param[in] "real(dp), intent(in) :: a_sh" fraction of projected area exchanging sensible heat with the air \f$1\f$
  !>        \param[in] "real(dp), intent(in) :: eta1" longwave emissivity of the leaf surface [1]

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
  !>        \return real(dp) :: R_ll &mdash; long-wave radiation away from leaf [W m^-2]

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !>         \note

  !     HISTORY
  !>        \author  Johannes Brenner
  !>        \date    Jan 2018
  !
  elemental pure FUNCTION R_ll(T_l, tavg, a_sh, eta1)

    use mo_constants, only: sigma_dp

    implicit none

    real(dp), intent(in) :: T_l              ! leaf surace temperature [degC]
    real(dp), intent(in) :: tavg             ! temperature [degC]
    real(dp), intent(in) :: a_sh             ! fraction of projected area exchanging sensible heat with the air [1]
    real(dp), intent(in) :: eta1             ! average one-sided convective transfer coefficient [J K-1 m-2 s-1]
    real(dp)             :: R_ll             ! long-wave radiation away from leaf [W m^-2]

    ! long-wave radiation away from leaf
    R_ll = 4 * a_sh * eta1 * sigma_dp * (tavg**4 * T_l - tavg**4)

  END FUNCTION R_ll

END MODULE mo_pet

!> \file mo_mrm_mpr.f90

!> \brief Perform Multiscale Parameter Regionalization on Routing Parameters

!> \details This module contains the subroutine for calculating the regionalized
!> routing parameters (beta-parameters) given the five global routing parameters
!> (gamma) at the level 0 scale.

!> \author Luis Samaniego, Stephan Thober
!> \date Aug 2015
module mo_mrm_mpr
  use mo_kind, only: dp
  implicit none
  public :: reg_rout, L11_calc_celerity
  private
contains
  
  ! ----------------------------------------------------------------------------

  !      NAME
  !         reg_rout

  !>        \brief Regionalized routing

  !>        \details sets up the Regionalized Routing parameters\n
  !>                 Global parameters needed (see mhm_parameter.nml):\n
  !>                    - param(1) = muskingumTravelTime_constant    \n
  !>                    - param(2) = muskingumTravelTime_riverLength \n
  !>                    - param(3) = muskingumTravelTime_riverSlope  \n
  !>                    - param(4) = muskingumTravelTime_impervious  \n
  !>                    - param(5) = muskingumAttenuation_riverSlope \n

  !      INTENT(IN)
  !>        \param[in] "real(dp) :: param(5)"  - five input parameters
  !>        \param[in] "real(dp) :: length(:)" - [m] total length
  !>        \param[in] "real(dp) :: slope(:)"  - average slope
  !>        \param[in] "real(dp) :: fFPimp(:)" - fraction of the flood plain with
  !>                                             impervious layer
  !>        \param[in] "real(dp) :: TS"        - [h] time step in

  !      INTENT(INOUT)
  !          None
  
  !      INTENT(OUT)
  !>        \param[out] "real(dp) :: C1(:)"    - routing parameter C1 (Chow, 25-41)
  !>        \param[out] "real(dp) :: C2(:)"    - routing parameter C2 (")

  !      INTENT(IN), OPTIONAL
  !          None

  !      INTENT(INOUT), OPTIONAL
  !          None

  !      INTENT(OUT), OPTIONAL
  !          None

  !      RETURN
  !          None

  !      RESTRICTIONS
  !          None

  !      EXAMPLE
  !          None

  !      LITERATURE
  !          None

  !      HISTORY
  !>        \author Stephan Thober, Rohini Kumar
  !>        \date Dec 2012
  !         Written Stephan Thober, Dec 2012

  subroutine reg_rout( param, length, slope, fFPimp, TS, &
       C1, C2 )

    implicit none

    ! Input
    real(dp), dimension(5), intent(in)  :: param  ! input parameter
    real(dp), dimension(:), intent(in)  :: length ! [m] total length
    real(dp), dimension(:), intent(in)  :: slope  ! average slope
    real(dp), dimension(:), intent(in)  :: fFPimp ! fraction of the flood plain with
    !                                             ! impervious layer
    real(dp),               intent(in)  :: TS     ! [h] time step in

    ! Output
    real(dp), dimension(:), intent(out) :: C1     ! routing parameter C1 (Chow, 25-41)
    real(dp), dimension(:), intent(out) :: C2     ! routing parameter C2 (")

    ! loval variables
    real(dp)                            :: ssMax  ! stream slope max
    real(dp), dimension(size(fFPimp,1)) :: K      ! [d] Muskingum travel time parameter
    real(dp), dimension(size(fFPimp,1)) :: xi     ! [1] Muskingum diffusion parameter (attenuation)

    ! normalize stream bed slope
    ssMax = maxval( slope(:) )

    ! New regional relationship; K = f(length, slope, & fFPimp)
    K = param(1) + param(2) * (length * 0.001_dp) &
         + param(3) * slope &
         + param(4) * fFPimp

    ! Xi = f(slope)
    xi = param(5)*(1.0_dp + slope / ssMax)

    ! constraints on Xi
    xi = merge( 0.5_dp, xi, xi > 0.5_dp )
    xi = merge( 0.005_dp, xi, xi < 0.005_dp )

    ! constrains on Ki
    K = merge( 0.5_dp * TS / xi,            K, K > 0.5_dp * TS / xi )
    K = merge( 0.5_dp * TS / (1.0_dp - xi), K, K < 0.5_dp * TS / (1.0_dp - xi))

    ! Muskingum parameters
    C1 = TS / ( K * (1.0_dp - xi) + 0.5_dp * TS )
    C2 = 1.0_dp - C1 * K / TS

  end subroutine reg_rout

  ! ----------------------------------------------------------------------------

  !      NAME
  !         L11_calc_celerity

  !>        \brief 

  !>        \details 

  !      INTENT(IN)

  !      INTENT(INOUT)
  !          None
  
  !      INTENT(OUT)

  !      INTENT(IN), OPTIONAL
  !          None

  !      INTENT(INOUT), OPTIONAL
  !          None

  !      INTENT(OUT), OPTIONAL
  !          None

  !      RETURN
  !          None

  !      RESTRICTIONS
  !          None

  !      EXAMPLE
  !          None

  !      LITERATURE
  !          None

  !      HISTORY
  !>        \author Matthias Kelbling
  !>        \date Aug 2017

  subroutine L11_calc_celerity(slope11, LinkIn_fAcc11, meandering11, nNodes)

    use mo_kind,                 only: i4, dp
    use mo_mrm_constants,        only: nodata_dp
    use mo_append,               only: append
    use mo_nml,                  only: open_nml, position_nml, close_nml
    use mo_mrm_file,             only: &
                file_namelist_param_mrm,   & 
                unamelist_param
    use mo_utils,                only: notequal
    use mo_mrm_global_variables, only: L11_celerity

    implicit none

    real(dp)                                 :: g1(5), g2(5), g3(5), g4(5) ! Parameter, CHANGE NAMES
    real(dp), dimension(:)                   :: meandering11
    real(dp), dimension(:)                   :: LinkIn_fAcc11
    real(dp), dimension(:)                   :: slope11
    real(dp), dimension(:), allocatable      :: celerity11
    integer(i4)                              :: nNodes
    integer(i4)                              :: kk

    ! Namelists are likely unnecessary
    namelist /routing3/g1,g2,g3,g4
    call open_nml(file_namelist_param_mrm, unamelist_param, quiet = .true.)
    call position_nml('routing3',unamelist_param)
    read(unamelist_param, nml=routing3)
    call close_nml(unamelist_param)

    ! allocate
    allocate(celerity11    (nNodes))

    ! initilize
    celerity11(:)    = nodata_dp

    ! calculate celerity
    do kk=1, nNodes
      if(  notequal(meandering11  (kk),    nodata_dp) .and.  &
           notequal(slope11       (kk),    nodata_dp) .and.  &
           notequal(LinkIn_fAcc11 (kk),    nodata_dp) ) then
        celerity11(kk) = ((g1(3)*LinkIn_fAcc11(kk)**g2(3))**(2.0/3.0) * & 
                         slope11(kk)**(1.0/2.0) ) / (g3(3)*meandering11(kk)**g4(3))
      end if
    end do

    call append(L11_celerity, celerity11(:))

  end subroutine L11_calc_celerity

end module mo_mrm_mpr

!> \file mo_mpr_neutrons.f90

!> \brief   Multiscale parameter regionalization (MPR) for neutrons

!> \details This module contains all routines required for parametrizing
!>          neutrons processes.

!> \author Maren Kaluza
!> \date Dec 2017

module mo_mpr_neutrons

  use mo_kind, only: i4, dp

  implicit none

  public :: mpr_neutrons

  private

contains
  ! ----------------------------------------------------------------------------

  !      NAME
  !         mpr_neutrons

  !>        \brief multiscale parameter regionalization for neutrons

  !>        \details calculates neutron variables on L0
  !>                 Global parameters needed (see mhm_parameter.nml):\n
  !>               TODO:     - param( 1) = orgMatterContent_forest     \n
  !>                    - param( 2) = orgMatterContent_impervious \n
  !>                    - param( 3) = orgMatterContent_pervious   \n
  !>                    - param( 4) = PTF_lower66_5_constant      \n
  !>                    - param( 5) = PTF_lower66_5_clay          \n
  !>                    - param( 6) = PTF_lower66_5_Db            \n
  !>                    - param( 7) = PTF_higher66_5_constant     \n
  !>                    - param( 8) = PTF_higher66_5_clay         \n
  !>                    - param( 9) = PTF_higher66_5_Db           \n

  !      INTENT(IN)
  !>        \param[in] "real(dp)    :: param(13)"        - global parameters

  !     INTENT(INOUT)
  !         None

  !      INTENT(OUT)
  !>                                                      capacity w.r.t to saturation

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
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
  !>        \author Maren Kaluza
  !>        \date Dec 2017


  subroutine mpr_neutrons( param , & ! IN:  global parameter set
       iFlag_soil          , & ! IN:  flag to handle different soil database
       is_present          , & ! IN:  flag indicating presence of soil
       nHorizons           , & ! IN:  Number of Horizons of Soiltype
       nTillHorizons       , & ! IN:  Number of tillage Horizons
       LCover0             , & ! IN:  land cover ids at level 0
       DbM                 , & ! IN:  mineral Bulk density
       Db                  , & ! IN: Bulk density
       COSMIC_L3_till      , & ! OUT: COSMIC paramter L3 tillage layer
       COSMIC_L3             & ! OUT: COSMIC paramter L3 tillage layer
       !                       !      hydraulic counductivity for Horizantal flow
       !                       !      hydraulic counductivity for Horizantal flow
       !                       !      w.r.t to saturation
       )

    ! lots of lines copy-pasted from mo_mpr_soilmoist.f90
    use mo_message,       only: message
    !$  use omp_lib

    implicit none

    ! Input --------------------------------------------------------------------
    real(dp),    dimension(13),    intent(in)  :: param        ! global parameters
    integer(i4),                   intent(in)  :: iFlag_soil   ! flag to handle different soil database

    integer(i4), dimension(:),     intent(in)  :: is_present   ! indicates whether soiltype is present
    integer(i4), dimension(:),     intent(in)  :: nHorizons    ! Number of Horizons per soiltype
    integer(i4), dimension(:),     intent(in)  :: nTillHorizons! Number of Tillage Horizons
    real(dp),    dimension(:,:),   intent(in)  :: DbM          ! mineral Bulk density
    real(dp),    dimension(:,:,:), intent(in)  :: Db           ! Bulk density
    integer(i4), dimension(:),     intent(in)  :: LCOVER0      ! land cover ids at level 0


    ! Output -------------------------------------------------------------------
    real(dp),    dimension(:,:,:), intent(out) :: COSMIC_L3_till! COSMIC parameter L3 tillage layer
    real(dp),    dimension(:,:),   intent(out) :: COSMIC_L3     ! COSMIC parameter L3
    !                                                           ! field cap. w.r.t to saturation
    ! Local variables
    integer(i4)                               :: i               ! loop index
    integer(i4)                               :: j               ! loop index
    integer(i4)                               :: l               ! loop index
    integer(i4)                               :: tmp_minSoilHorizon

    tmp_minSoilHorizon = minval(nTillHorizons(:))

    COSMIC_L3_till  = 0.0_dp
    COSMIC_L3    = 0.0_dp

    ! select case according to a given soil database flag
    SELECT CASE(iFlag_soil)
       ! classical mHM soil database format
       CASE(0)
          do i = 1, size(is_present)
             if ( is_present(i) .lt. 1 ) cycle
             horizon: do j = 1, nHorizons(i)
                ! calculating other soil hydraulic properties
                ! tillage horizons
                if ( j .le. nTillHorizons(i) ) then
                   ! LC class
                   do L = 1, maxval( LCOVER0 )
                      call calcL3(param(4:9), Db(i,j,L), COSMIC_L3_till(i,j,L))
                   end do
                ! deeper layers
                else
                   call calcL3(param(4:9), DbM(i,j), COSMIC_L3(i,j-tmp_minSoilHorizon))
                end if
             end do horizon
          end do
       ! to handle multiple soil horizons with unique soil class   
       CASE(1)
           do i = 1, size(is_present)
             if ( is_present(i) .lt. 1 ) cycle
             ! **** FOR THE TILLAGE TYPE OF SOIL *****
             ! there is actually no soil horizons/soil type in this case
             ! but we assign of j = 1 to use variables as defined in the classical option (iFlag_soil = 0)
             do j = 1, 1   
                ! tillage horizons properties depending on the LC class
                do L = 1, maxval( LCOVER0 )
                   call calcL3(param(4:9), Db(i,j,L), COSMIC_L3_till(i,j,L))
                end do
                
                ! *** FOR NON-TILLAGE TYPE OF SOILS ***
                ! note j = 1
                call calcL3(param(4:9), DbM(i,j), COSMIC_L3(i,j))

             end do  !>> HORIZON
          end do   !>> SOIL TYPE
       CASE DEFAULT
          call message()
          call message('***ERROR: iFlag_soilDB option given does not exist. Only 0 and 1 is taken at the moment.')
          stop
       END SELECT
    
   end subroutine


  subroutine calcL3(param, bulkDensity, L3)
     implicit none
     real(dp), dimension(4),  intent(in)       :: param
     real(dp),                intent(in)       :: bulkDensity
     real(dp),                intent(inout)    :: L3
      L3 = bulkDensity*106.194175956_dp - 40.987888406_dp
      if (bulkDensity < 0.4) then ! bulkDensity<0.39 yields negative L3, bulkDensity=0.39 yields L3=0
         L3 = 1.0 ! Prevent division by zero later on; added by joost Iwema to COSMIC 1.13, Feb. 2017
      endif
  end subroutine
end module mo_mpr_neutrons

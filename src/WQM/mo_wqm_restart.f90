!> \file mo_wqm_restart.f90

!> \brief reading and writing states, fluxes and configuration for restart of water quality model.

!> \details routines are seperated for reading and writing variables for:\n
!>          - states and fluxes, and \n
!>          - configuration.\n


!> \authors Xiaoqiang Yang, Modified from original mHM code
!> \date  Jul 2017 

MODULE mo_wqm_restart



  IMPLICIT NONE

  PUBLIC :: wqm_read_restart_states     ! 

  PUBLIC :: wqm_write_restart_files     ! 

  PRIVATE

CONTAINS
  ! ------------------------------------------------------------------
  
  !      NAME
  !         wqm_write_restart_files

  !     PURPOSE
  !>        \brief write restart files for each basin

  !>        \details write restart files for each basin. For each basin
  !>        xxx_wqm_states.nc(xxx being the three digit
  !>        basin index) is written. If a variable is added here, it should also be added
  !>        in the read restart routines below.

  !     INTENT(IN)
  !>        \param[in] "character(256), dimension(:) :: OutPath"     Output Path for each basin

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

  !     RESTRICTIONS 
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         Modified from the original mHM and mRM

  !     HISTORY
  !>        \author   Xiaoqiang Yang
  !>        \date     Jul 2017
  !         Modified  

  ! ------------------------------------------------------------------ 
  subroutine wqm_write_restart_files(OutPath )

    use mo_kind,             only: i4, dp
    use mo_message,          only: message
    use mo_init_states,      only: get_basin_info
    use mo_mrm_tools,        only: get_basin_info_mrm
    use mo_string_utils,     only: num2str
    use mo_netcdf,           only: NcDataset, NcDimension, NcVariable
    use mo_mhm_constants,    only: nodata_dp
    use mo_wqm_global_variables, only: &
         L1_humusN, L1_fastN,L1_dissolvedIN,L1_dissolvedON, &
         L1_csoilMoist,L1_cbaseflow,L1_cunsatSTW,L1_cpercolate

    implicit none

    character(256)                           :: Fname
    character(256),dimension(:), intent(in)  :: OutPath ! list of Output paths per Basin
    integer(i4)                              :: iBasin
    integer(i4)                              :: ii, jj
    integer(i4)                              :: s0       ! start index at level 0
    integer(i4)                              :: e0       ! end index at level 0
    integer(i4)                              :: ncols0   ! number of colums at level 0
    integer(i4)                              :: nrows0   ! number of rows at level 0
    logical, dimension(:,:), allocatable     :: mask0    ! mask at level 0
    integer(i4)                              :: s1       ! start index at level 1
    integer(i4)                              :: e1       ! end index at level 1
    integer(i4)                              :: ncols1   ! number of colums at level 1
    integer(i4)                              :: nrows1   ! number of rows at level 1
    logical, dimension(:,:), allocatable     :: mask1    ! mask at level 1
    real(dp), dimension(:,:,:), allocatable  :: dummy_d3 ! dummy variable
    real(dp), dimension(:,:,:,:), allocatable:: dummy_d4
    integer(i4)                              :: s11       ! start index at level 11
    integer(i4)                              :: e11       ! end index at level 11
    integer(i4)                              :: ncols11   ! number of colums at level 11
    integer(i4)                              :: nrows11   ! number of rows at level 11
    logical, dimension(:,:), allocatable     :: mask11    ! mask at level 11

    type(NcDataset)                          :: nc
    type(NcDimension)                        :: rows1, cols1, rows11,cols11, soil1,nsubst
    type(NcVariable)                         :: var
    
    basin_loop: do iBasin = 1, size(OutPath)

       ! get Level0 information about the basin
       call get_basin_info( iBasin, 0, nrows0, ncols0, iStart=s0, iEnd=e0, mask=mask0 )

       ! get Level1 information about the basin
       call get_basin_info( iBasin, 1, nrows1, ncols1, iStart=s1, iEnd=e1, mask=mask1 )
       ! level 11 info
       call get_basin_info_mrm( iBasin, 11, nrows11, ncols11, iStart=s11, iEnd=e11, mask=mask11 )

       ! write restart file for iBasin
       Fname = trim(OutPath(iBasin)) // "WQM_restart_" // trim(num2str(iBasin, "(i3.3)")) // ".nc"
       ! print a message
       call message("    Writing Restart-file: ", trim(adjustl(Fname))," ...")

       nc     = NcDataset(fname, "w")
       rows1  = nc%setDimension("nrows1 ", nrows1)
       cols1  = nc%setDimension("ncols1 ", ncols1)
       rows11  = nc%setDimension("nrows11", nrows11)
       cols11  = nc%setDimension("ncols11", ncols11)
       soil1  = nc%setDimension("soilhorizons", size( L1_csoilMoist, 2))
       nsubst  = nc%setDimension("substances", size( L1_csoilMoist, 3))

       allocate( dummy_d3( nrows1, ncols1, size( L1_csoilMoist, 2) ) )
       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_humusN(s1:e1,ii), mask1, nodata_dp )
       end do
       var = nc%setVariable("L1_humusN","f64",(/rows1,cols1,soil1/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","soil humusN storage at level 1")
	   
       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_fastN(s1:e1,ii), mask1, nodata_dp )
       end do      
       var = nc%setVariable("L1_fastN","f64",(/rows1,cols1,soil1/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","soil fastN storage at level 1")

       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_dissolvedIN(s1:e1,ii), mask1, nodata_dp )
       end do      
       var = nc%setVariable("L1_dissolvedIN","f64",(/rows1,cols1,soil1/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","soil dissolved IN storage at level 1")
      
       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_dissolvedON(s1:e1,ii), mask1, nodata_dp )
       end do      
       var = nc%setVariable("L1_dissolvedON","f64",(/rows1,cols1,soil1/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","soil dissolved ON storage at level 1")

       allocate( dummy_d4( nrows1,ncols1, size(L1_csoilMoist,2), size(L1_csoilMoist,3)))
       do ii = 1, size(dummy_d4,4)
           do jj =1, size(dummy_d4,3)
               dummy_d4(:,:,jj,ii) = unpack( L1_csoilMoist(s1:e1,jj,ii), mask1, nodata_dp )
           end do
       end do
       var = nc%setVariable("L1_csoilMoist","f64",(/rows1,cols1, soil1, nsubst/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d4)
       call var%setAttribute("long_name","soil moisture concentration at level 1")

       deallocate(dummy_d3)
       allocate(dummy_d3( nrows1, ncols1, size( L1_csoilMoist, 3) ))

       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_cbaseflow(s1:e1,ii), mask1, nodata_dp )
       end do      
       var = nc%setVariable("L1_cbaseflow","f64",(/rows1,cols1,nsubst/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","base flow concentration at level 1")

       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_cunsatSTW(s1:e1,ii), mask1, nodata_dp )
       end do      
       var = nc%setVariable("L1_cunsatSTW","f64",(/rows1,cols1,nsubst/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","unsaturated storage concentration at level 1")

       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_cpercolate(s1:e1,ii), mask1, nodata_dp )
       end do      
       var = nc%setVariable("L1_cpercolate","f64",(/rows1,cols1,nsubst/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","percolated water concentration at level 1")



       call nc%close()
       
    end do basin_loop
    
  end subroutine wqm_write_restart_files


  ! ------------------------------------------------------------------

  !      NAME
  !         wqm_read_restart_states

  !     PURPOSE
  !>        \brief reads fluxes and state variables from file

  !>        \details read fluxes and state variables from given 
  !>        restart directory and initialises all state variables
  !>        that are initialized in the subroutine wqm_initialise,
  !>        contained in module mo_water_quality.

  !     INTENT(IN)
  !>        \param[in] "integer(i4)    :: iBasin"        number of basin


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

  !     RESTRICTIONS 
  !>        \note Restart Files must have the format, as if
  !>        it would have been written by subroutine write_restart_files

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Xiaoqiang Yang  -- Modified from the original mHM and mRM
  !>        \date Jun 2017

  subroutine wqm_read_restart_states( iBasin )

    use mo_kind,             only: i4, dp
    ! use mo_message,          only: message
    use mo_netcdf,           only: NcDataset, NcVariable
    use mo_string_utils,     only: num2str
    use mo_init_states,      only: get_basin_info
    !use mo_mhm_constants,    only: YearMonths_i4
    use mo_global_variables, only:  &
         nSoilHorizons_mHM,dirRestartIn
    use mo_wqm_global_variables, only: &
         L1_humusN, L1_fastN,L1_dissolvedIN,L1_dissolvedON, &
         L1_csoilMoist,L1_cbaseflow,L1_cunsatSTW,L1_cpercolate, &
         nsubstances

    implicit none

    integer(i4),    intent(in) :: iBasin


    character(256)                                    :: Fname
    integer(i4)                                       :: ii,jj       ! loop index
    integer(i4)                                       :: s1       ! start index at level 1
    integer(i4)                                       :: e1       ! end index at level 1
    integer(i4)                                       :: ncols1   ! number of colums at level 1
    integer(i4)                                       :: nrows1   ! number of rows at level 1
    integer(i4)                                       :: ncells1  ! number of cells at level 1
    logical, dimension(:,:), allocatable              :: mask1    ! mask at level 1

    !real(dp), dimension(:,:),   allocatable           :: dummyD2  ! dummy, 2 dimension
    real(dp), dimension(:,:,:), allocatable           :: dummyD3  ! dummy, 3 dimension
    real(dp), dimension(:,:,:,:), allocatable         :: dummyD4

    type(NcDataset) :: nc
    type(NcVariable) :: var
    
    Fname = trim(dirRestartIn(iBasin)) // 'WQM_restart_' // trim(num2str(iBasin, '(i3.3)')) // '.nc'
    ! call message('    Reading states from ', trim(adjustl(Fname)),' ...')

    
    ! get basin information at level 1
    call get_basin_info( iBasin, 1, nrows1, ncols1, ncells=ncells1, &
         iStart=s1, iEnd=e1, mask=mask1 )

    nc = NcDataset(fname,"r")
  

    var = nc%getVariable("L1_humusN")
    call var%getData(dummyD3)
    do ii = 1, nSoilHorizons_mHM
       L1_humusN(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do
    var = nc%getVariable("L1_fastN")
    call var%getData(dummyD3)
    do ii = 1, nSoilHorizons_mHM
       L1_fastN(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do
    var = nc%getVariable("L1_dissolvedIN")
    call var%getData(dummyD3)
    do ii = 1, nSoilHorizons_mHM
       L1_dissolvedIN(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do
    var = nc%getVariable("L1_dissolvedON")
    call var%getData(dummyD3)
    do ii = 1, nSoilHorizons_mHM
       L1_dissolvedON(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do
    var = nc%getVariable("L1_csoilMoist")
    call var%getData(dummyD4)
    do ii = 1, nsubstances
       do jj = 1, nSoilHorizons_mHM
       L1_csoilMoist(s1:e1, jj,ii) = pack( dummyD4( :,:,jj, ii), mask1)
       end do
    end do
    deallocate(dummyD3)
    var = nc%getVariable("L1_cbaseflow")
    call var%getData(dummyD3)
    do ii = 1, nsubstances
       L1_cbaseflow(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do
    var = nc%getVariable("L1_cunsatSTW")
    call var%getData(dummyD3)
    do ii = 1, nsubstances
       L1_cunsatSTW(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do
    var = nc%getVariable("L1_cpercolate")
    call var%getData(dummyD3)
    do ii = 1, nsubstances
       L1_cpercolate(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do

    call nc%close()

  end subroutine wqm_read_restart_states

END MODULE mo_wqm_restart

!
! GPM GMI simulator
!
!#define IN_ACME
module GPM_GMI_MOD

   !------------------
   ! Environment setup
   !------------------
   ! Module use
   use GPM_CRTM_sensor_mod, only: gpm_CRTM_sensor_type
   use CRTM_Module!,    only: CRTM_ChannelInfo_type, CRTM_Geometry_type, &
                  !           CRTM_Options_type, CRTM_Atmosphere_type,   &
                  !           CRTM_Surface_type, CRTM_RTSolution_type,   &
                  !           CRTM_Init
   use MOD_COSP_TYPES, only: cosp_gridbox
   use MOD_COSP_CONSTANTS, only: I_LSCLIQ, I_LSCICE, I_LSRAIN, I_LSSNOW, &
                                 I_CVCLIQ, I_CVCICE, I_CVRAIN, I_CVSNOW, &
                                 I_LSGRPL
   implicit none
   private

!-------------------------------
   ! gpm_gmi_toabt is a structure to store calculated TOA BT
   type gpm_gmi_toabt
      ! The sensor associated with this TOA BT
      type(gpm_gmi_sensor) :: sensor
      ! Here we only need the brightness temperature. Instead saving all the
      ! RTSolution structures in CRTM, we only obtain the brightness
      ! temperatures from CRTM RTSolution structures, and then destroy the
      ! RTSolution structures to save memory.
      real, allocatable :: toa_bt(:,:)  ! n_challels x n_profile

   end type gpm_gmi_toabt
!-------------------------------
   ! max_sensors is a private parameter denoting the maximum number of sensors
   ! supported. It determines the sensor array size. If more sensors is needed,
   ! change the value.
   integer, parameter :: max_sensors = 10

   ! To save memory, the RT calculation can be divided into many iterations. The
   ! variable max_profile_iter is the maximum number of profiles processed in
   ! each iteration.
   integer, parameter :: max_profile_iter = 1000

   ! num_sensor is a private variable denoting how many sensors is registered.
   ! This number starts with 0 and incrases by 1 with each calling of gpm_gmi_addsensor
   integer ::   n_sensors=0
   
   ! l_initialized is a logical variable indicating if the gpm-gmi-sensors are
   ! initialized. After initialization, cannot add sensors any more.
   logical :: l_initialized = .FALSE.
   
   ! CRTM coefficient data file location. It is hardcoded for now.
   ! TODO: dynamically set the location of the coefficient data based on
   ! configuration
   character(100), parameter :: dircoef = '/global/homes/y/yxl232/local/CRTM/coeff_data/Big_Endian_ODAS/'
   !character(100), parameter :: dircoef = 'coeff_data/'
   
   ! pre-allocated array of sensors
   type(gpm_gmi_sensor), target :: sensors(max_sensors)

   ! variables associated with sensors. Information for all sensors has to be 
   ! in one array because of the requirements of CRTM library.
   ! Length: n_sensors

   character(20)               , allocatable :: sensor_id(:)
   type(CRTM_ChannelInfo_type) , allocatable :: chinfo(:)



   public :: gpm_gmi_addsensor, &
             gpm_gmi_adddefaultsensor, &
             gpm_gmi_init,      &
             gpm_gmi_run,       &
             gpm_gmi_clean      
   public :: gpm_gmi_register,  &
             gpm_gmi_addfld
   public :: gpm_gmi_toabt
      
!=======================================
CONTAINS
!=======================================
   subroutine gpm_gmi_register
!---------------------------------------
!
! Register gpm gmi fields in the physics buffer
!
!--------------------------------------
   end subroutine gpm_gmi_register
!=======================================
   subroutine gpm_gmi_addfld
!---------------------------------------
!
! Add gpm gmi fields to history files
!
!--------------------------------------

   use cam_history,         only: addfld, add_default, phys_decomp
   use mod_cosp_constants,  only: R_UNDEF
   
   implicit none
   
  ! TODO: The value of cosp_histfile_num is hardcoded to 1 for now. Need to
  ! be modified so that it is consistant with the values in
  ! cospsimulator_intr.F90
   integer, parameter :: cosp_histfile_num = 1
   
   ! looping variable
   integer :: n=0
   type(gpm_gmi_sensor), pointer :: currentsensor
  !----------------
  ! Loop through all sensors
   do n = 1, n_sensors
      currentsensor => sensors(n)
#ifdef IN_ACME
      call addfld(currentsensor%cam_histfld_name, 'K', currentsensor%n_channels, 'A', &
                  'Simulated GPM GMI brightness temperature',phys_decomp, &
                  flag_xyfill = .true., fill_value = R_UNDEF)
#endif
   end do
   end subroutine gpm_gmi_addfld

   
!======================================

! TODO: find a way to read in sensor information, e.g., from an xml file
   subroutine gpm_gmi_addsensor(sensor_id_in, cam_histfld_name_in, &
              sensor_scan_angle_in, sensor_zenith_angle_in,        &
              channel_subset)
!--------------------------
!
!  Add a sensor and update n_sensors
!
!-------------------------
   character(len=*), intent(in) :: sensor_id_in
   character(len=*), intent(in) :: cam_histfld_name_in
   real, intent(in) :: sensor_scan_angle_in
   real, intent(in) :: sensor_zenith_angle_in
   integer, intent(in), optional :: channel_subset(:)
! local variable
   integer :: err_stat
   type(gpm_gmi_sensor), pointer :: currentsensor
   ! check if the number if sensors reach the maximum value
   if (n_sensors == max_sensors) then
      print *, "number of sensors exceeds maximum number. Increase the value of &
      max_sensors in gpm_gmi.F90 and recompile"
      stop
   end if
   ! check if the gpm-gmi structures are already initialized.
   if (l_initialized == .TRUE. ) then
      print *, "GPM GMI structures are already initialized, cannot add more &
      sensors. Call gpm_gmi_addsensor() before calling gpm_gmi_init() "
      stop
   end if

   n_sensors = n_sensors + 1
   currentsensor => sensors(n_sensors)
   ! set properties to current sensor
   currentsensor%sensor_id        = sensor_id_in
   currentsensor%cam_histfld_name = cam_histfld_name_in
   currentsensor%sensor_scan_angle = sensor_scan_angle_in
   currentsensor%sensor_zenith_angle = sensor_zenith_angle_in
   ! allocate space for channel_subset if needed
   if (present(channel_subset) ) then
      allocate(currentsensor%channel_subset(size(channel_subset)))
      currentsensor%channel_subset = channel_subset
   endif
   print *, currentsensor%cam_histfld_name
   end subroutine gpm_gmi_addsensor

!============================================
! TODO: use xml file to store the sensor properties
! The current implementation hard-wired all the prpoerties
   subroutine gpm_gmi_adddefaultsensor(sensorname)
!-----------------------
!
! Add a default sensor
!
!-----------------------
      character(len=*), intent(in) :: sensorname
      !----------
      ! local variable
      integer :: i
      !----------
      select case (sensorname)
         case ('gpm-gmi-lowfreq')
            call gpm_gmi_addsensor('gmi_gpm','GPM_GMI_LOWFREQ',52.8, 48.5,&
            channel_subset = (/(i,i=1,9)/)  )
         case ('gpm-gmi-highfreq')
            call gpm_gmi_addsensor('gmi_gpm','GPM_GMI_LOWFREQ',49.19, 45.36,&
            channel_subset = (/(i,i=10,13)/)  )
         case default
            print *, 'error: sensor not implemented in default list'
            stop
      end select


   end subroutine gpm_gmi_adddefaultsensor
!============================================
   
   subroutine gpm_gmi_init(gbx)
!-----------------------
!
! Initialize variables
!
!-----------------------

      ! input variable
      type(cosp_gridbox), intent(in) :: gbx


      ! local variables
      type(gpm_gmi_sensor), pointer :: currentsensor
      integer :: err_stat
      integer :: alloc_stat
      integer :: n
      integer :: n_layers
      integer :: n_profiles

      ! Check if GPM GMI structures are already initialized
      if ( l_initialized ) then
         print *, "GPM GMI structures are already initialized. Cannot initialize &
            again"
         stop
      end if

      ! Check if there are already sensors added. If no sensors added, stop
      ! running 
      if (n_sensors <= 0) then
         print *, "No sensors added. Add sensor first then run gpm_gmi_init()"
         stop
      end if

      !------------------------
      ! Allocate sensor_id and chinfo arrays
      allocate( sensor_id(n_sensors), chinfo(n_sensors), &
                 STAT = alloc_stat)
      if (alloc_stat /=0) then
         write(*,*) "fail to allocate sensor_id or chinfo"
         stop
      endif
      ! Copy values of sensor_id from each sensor
      do n = 1, n_sensors
         sensor_id(n) = sensors(n)%sensor_id
      end do 
      ! Initialize chinfo
      err_stat = CRTM_Init(sensor_id, chinfo, &
                           File_Path = dircoef, &
                           Quiet = .TRUE.)
      if (err_stat /= SUCCESS) then
         write(*,*) "fail to initialize chinfo"
         STOP
      endif
      n_profiles = gbx%Npoints
      do n=1, n_sensors
         currentsensor => sensors(n)
         
         ! subset channels if necessary
         ! currently, whether a sensor needs channel subsetting is determined by
         ! whether its channel_subset properties is set.
         if (allocated(currentsensor%channel_subset)) then
            err_stat = CRTM_ChannelInfo_Subset(chinfo(n), &
                           Channel_Subset = currentsensor%channel_subset )
            if (err_stat /= SUCCESS) then
               print *, "failed to subset channels"
               stop
            end if
         end if
         call CRTM_ChannelInfo_Inspect(chinfo(n))
         ! allocate arrays
         currentsensor%n_channels = CRTM_ChannelInfo_n_Channels( chinfo(n) )
!         allocate( currentsensor%toa_bt(currentsensor%n_channels, n_profiles))
      end do !n_sensors


      ! finished initialization. Set l_initialized to .TRUE.
      l_initialized = .TRUE.

   end subroutine gpm_gmi_init

!================================================
   subroutine gpm_gmi_run(gbx, toabt)
!-----------------------
!
!  Run GPM GMI calculations
!
!----------------------
      type(COSP_GRIDBOX), intent(in) :: gbx
      type(gpm_gmi_toabt), intent(out), allocatable :: toabt(:) ! size: n_sensors

      ! local variable
      integer :: n
      type(gpm_gmi_sensor), pointer :: currentsensor
      integer :: m
      integer :: n_profiles
      !----------------
      ! Check if GPM GMI structures are already initialized 
      if ( l_initialized == .FALSE. ) then
         print *, "GPM GMI structures are not initialized. Need to initialize &
                  them first, and then run calculations "
         stop
      end if
      ! allocate toabt
      n_profiles = gbx%Npoints
      allocate(toabt(n_sensors))
      do n=1, n_sensors
         toabt(n)%sensor = sensors(n)
         allocate(toabt(n)%toa_bt(sensors(n)%n_channels, n_profiles))
      end do
      
      ! Loop through sensors and run calculation one by one
      print *, "n_sensors=", n_sensors
      do n=1, n_sensors
!         currentsensor => sensors(n)
!         err_stat = CRTM_Forward(atm, sfc, currentsensor%geo, chInfo(n:n),&
!                                 currentsensor%rts)
!      print *, err_stat

         call gpm_gmi_run_iter(n, 1, 1, gbx, toabt(n))
         call gpm_gmi_run_iter(n, 2, 2, gbx, toabt(n))
         do m=1, size(toabt(n)%toa_bt, dim=2)
            print *, "-----results for sensor ", n, "-----"
            print *, toabt(n)%toa_bt(:,m)
            print *, "-----------------"
         end do
      end do
   end subroutine gpm_gmi_run

!==================================================
   subroutine gpm_gmi_run_iter(isensor, pbeg, pend, gbx, toabt)
!--------------------------------
!
! Run GPM GMI for one iterature to save memory
!
!--------------------------------
      integer, intent(in) :: isensor ! index of current sensor
      integer, intent(in) :: pbeg ! starting profile
      integer, intent(in) :: pend ! ending profile
      type(COSP_GRIDBOX), intent(in) :: gbx ! gridbox
      type(gpm_gmi_toabt), intent(out) :: toabt ! output structure
      
      ! local variable
      integer :: err_stat
      ! Geometry and option structure. Length: n_profile
      type(CRTM_Geometry_type)    , allocatable :: geo(:) 
      type(CRTM_Options_type)     , allocatable :: opt(:)       

      ! RTSolution structure used in this iteration. Size: n_channel x np
      type(CRTM_RTSolution_type)  , allocatable :: rts(:,:)  
      type(gpm_gmi_sensor), pointer :: currentsensor
      integer :: np ! number of profiles processed in this iteration
      ! variables associated with atmosphere states. 
      ! They are shared by all the sensors. 
      ! Length: n_profile
      type(CRTM_Atmosphere_type)  , allocatable :: atm(:)       
      type(CRTM_Surface_type)     , allocatable :: sfc(:)       
      integer :: iatm ! index for atm variable
      integer :: igbx ! index for gbx variable
      integer :: n_layers ! number of atmosphere layers
      ! variable associated to the whole profile of gbx
      integer :: n_profiles = 0
      integer ::  n_absorbers = 2
      integer ::  n_clouds    = 9
      integer ::  n_aerosols  = 1
      integer :: tmpint
      !---------------------
      ! Set parameters obtained from gbx
      n_profiles = gbx%Npoints      
      n_layers   = gbx%Nlevels
      ! checking if the pbeg and pend within valid range
      if (pbeg <= 0 .OR. pend > n_profiles) then
         print *, "pbeg or pend out of range"
         stop
      end if
      ! checking if the isensor is within range (isensor < n_sensors)
      if (isensor <=0 .OR. isensor > n_sensors) then
         print *, "isensor out of range"
         stop
      end if
      !-------------------
      np = pend - pbeg +1 ! number of profiles in current iteration

      ! allocate profile-only arrays. atm and sfc are shared by all sensors
      allocate(atm(np), sfc(np) )
print *, "12345678"
      ! create atmosphere profile.
      call CRTM_Atmosphere_Create( atm, n_layers, n_absorbers, n_clouds,&
                                   n_aerosols)
      if ( ANY (.NOT. CRTM_Atmosphere_Associated(atm) ) )  then
         print *, "fail to create atmosphere profile"
         STOP
      endif

      ! Profile and absorber definition
      ! need to loop through profiles
      do iatm=1, np
      igbx = pbeg + iatm -1
      print *, "iatm=",iatm, "  igbx=", igbx
         atm(iatm)%Climatology         = US_STANDARD_ATMOSPHERE
         atm(iatm)%Absorber_Id(1:2)    = (/H2O_ID                 ,  O3_ID/)
         atm(iatm)%Absorber_units(1:2) = (/MASS_MIXING_RATIO_UNITS,  VOLUME_MIXING_RATIO_UNITS/)

         ! Profile data
         ! COSP is from surface to top, but CRTM wants top to surface. Need to flip
         ! the profiles
         atm(iatm)%Level_Pressure = gbx%pint(igbx,(n_layers+1):1:-1)
         atm(iatm)%Pressure       = gbx%p(igbx,n_layers:1:-1)
         atm(iatm)%Temperature    = gbx%T(igbx,n_layers:1:-1)

         ! unit for q in COSP is kg/kg, while CRTM wants g/kg. Need to do unit
         ! conversion
         atm(iatm)%Absorber(:,1)  = gbx%q(igbx,n_layers:1:-1) * 1000
         ! units for O3 mixing ratio in COSP is kg/kg, but CRTM needs ppmv. The
         ! conversion is not implemented here yet.
         
         ! assign cloud properties
         ! For now, the cloud types follow those in COSP. Each type of
         ! hydrometer defined in COSP is set to be a single cloud
         ! TODO: current implementation do not distinguish cloud or clear. May
         ! need to define number of clouds based on clear/cloud
         ! TODO: Make sure the units agree. In CRTM, Effective radius has unit microns,
         ! while water content has unit kg/m^2

         tmpint = I_LSCLIQ
         atm(iatm)%Cloud(1)%Type = WATER_CLOUD
         atm(iatm)%Cloud(1)%Effective_Radius = gbx%Reff(igbx,n_layers:1:-1,tmpint)
         atm(iatm)%Cloud(1)%Water_Content    = gbx%mr_hydro(igbx,n_layers:1:-1,tmpint)
         
         tmpint = I_LSCICE
         atm(iatm)%Cloud(2)%Type = ICE_CLOUD
         atm(iatm)%Cloud(2)%Effective_Radius = gbx%Reff(igbx,n_layers:1:-1,tmpint)
         atm(iatm)%Cloud(2)%Water_Content    = gbx%mr_hydro(igbx,n_layers:1:-1,tmpint)
         
         tmpint = I_LSRAIN
         atm(iatm)%Cloud(3)%Type = RAIN_CLOUD
         atm(iatm)%Cloud(3)%Effective_Radius = gbx%Reff(igbx,n_layers:1:-1,tmpint)
         atm(iatm)%Cloud(3)%Water_Content    = gbx%mr_hydro(igbx,n_layers:1:-1,tmpint)
         
         tmpint = I_LSSNOW
         atm(iatm)%Cloud(4)%Type = SNOW_CLOUD
         atm(iatm)%Cloud(4)%Effective_Radius = gbx%Reff(igbx,n_layers:1:-1,tmpint)
         atm(iatm)%Cloud(4)%Water_Content    = gbx%mr_hydro(igbx,n_layers:1:-1,tmpint)
         
         tmpint = I_CVCLIQ
         atm(iatm)%Cloud(5)%Type = WATER_CLOUD
         atm(iatm)%Cloud(5)%Effective_Radius = gbx%Reff(igbx,n_layers:1:-1,tmpint)
         atm(iatm)%Cloud(5)%Water_Content    = gbx%mr_hydro(igbx,n_layers:1:-1,tmpint)
         
         tmpint = I_CVCICE
         atm(iatm)%Cloud(6)%Type = ICE_CLOUD
         atm(iatm)%Cloud(6)%Effective_Radius = gbx%Reff(igbx,n_layers:1:-1,tmpint)
         atm(iatm)%Cloud(6)%Water_Content    = gbx%mr_hydro(igbx,n_layers:1:-1,tmpint)
         
         tmpint = I_CVRAIN
         atm(iatm)%Cloud(7)%Type = RAIN_CLOUD
         atm(iatm)%Cloud(7)%Effective_Radius = gbx%Reff(igbx,n_layers:1:-1,tmpint)
         atm(iatm)%Cloud(7)%Water_Content    = gbx%mr_hydro(igbx,n_layers:1:-1,tmpint)
         
         tmpint = I_CVSNOW
         atm(iatm)%Cloud(8)%Type = SNOW_CLOUD
         atm(iatm)%Cloud(8)%Effective_Radius = gbx%Reff(igbx,n_layers:1:-1,tmpint)
         atm(iatm)%Cloud(8)%Water_Content    = gbx%mr_hydro(igbx,n_layers:1:-1,tmpint)
         
         tmpint = I_LSGRPL
         atm(iatm)%Cloud(9)%Type = GRAUPEL_CLOUD
         atm(iatm)%Cloud(9)%Effective_Radius = gbx%Reff(igbx,n_layers:1:-1,tmpint)
         atm(iatm)%Cloud(9)%Water_Content    = gbx%mr_hydro(igbx,n_layers:1:-1,tmpint)
         
         
         
         
         ! assign aerosol properties
         ! TODO: aerosol properties and land cover are dummy values for now.
         ! Need to get them from the model
         atm(iatm)%Aerosol(1)%Type = DUST_AEROSOL
         atm(iatm)%Aerosol(1)%Effective_Radius(3:4) = (/7e-16_fp, 3e-3_fp/)
         atm(iatm)%Aerosol(1)%Concentration(3:4) = (/2e-18_fp, 1e-5_fp/)
         !!! assign surface properties
         !sfc(iatm)%Land_Coverage = 0.1_fp
         !sfc(iatm)%Water_Coverage = 0.5_fp
         !sfc(iatm)%Snow_Coverage = 0.25_fp
         !sfc(iatm)%Ice_Coverage = 0.15_fp

         sfc(iatm)%Land_Coverage = 0.0_fp
         sfc(iatm)%Water_Coverage = 1.0_fp
         sfc(iatm)%Snow_Coverage = 0.0_fp
         sfc(iatm)%Ice_Coverage = 0.0_fp
      end do !n_profiles
      ! allocate sensor specific structures
      currentsensor => sensors(isensor)
      allocate(rts(currentsensor%n_channels, np))
      
      allocate( geo(np), opt(np) )
      ! Set values
         
      !--------------
      geo%Latitude  = gbx%latitude
      geo%Longitude = gbx%longitude
      geo%Sensor_Scan_Angle   = currentsensor%sensor_scan_angle
      geo%Sensor_Zenith_Angle = currentsensor%sensor_zenith_angle

      ! set values of Geo structures
      ! TODO: get the model year/month/day from model run
      geo%Year  = 2001
      geo%Month = 01
      geo%Day   = 01
         
print *, "-------12345678--------" 
      err_stat = CRTM_Forward(atm, sfc, geo, chInfo(isensor:isensor),&
                              rts)
      print *, err_stat

      ! obtain brightness temperatures from RTSolution structures
      toabt%toa_bt(:, pbeg:pend) = rts%Brightness_temperature

      ! destroy the rts structure
      call CRTM_Geometry_Destroy(  geo)
      call CRTM_Options_Destroy(   opt)
      call CRTM_RTSolution_Destroy(rts)
      call CRTM_Atmosphere_Destroy(atm)

   end subroutine gpm_gmi_run_iter

!================================================
   subroutine gpm_gmi_clean()
!-----------------------
!
!  Clean GPM GMI structures
!
!----------------------
      integer :: err_stat
      type(gpm_gmi_sensor), pointer :: currentsensor
      integer n
      err_stat = CRTM_Destroy(chinfo)
   end subroutine gpm_gmi_clean





end module GPM_GMI_MOD



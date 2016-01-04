!
! GPM GMI simulator
!

module GPM_GMI_MOD

   !------------------
   ! Environment setup
   !------------------
   ! Module use
   use CRTM_Module!,    only: CRTM_ChannelInfo_type, CRTM_Geometry_type, &
                  !           CRTM_Options_type, CRTM_Atmosphere_type,   &
                  !           CRTM_Surface_type, CRTM_RTSolution_type,   &
                  !           CRTM_Init
   use MOD_COSP_TYPES, only: cosp_gridbox
   implicit none
      ! Dimension variables
      integer :: n_channels  ! l = 1, ... , L 
      integer :: n_profiles  ! m = 1, ... , M
      integer :: n_sensors   ! n = 1, ... , N
      
      ! CRTM coefficient data file location
      character(100) :: dircoef = &
         'coeff_data/'

      ! Processing parameters 
      character(20)               , allocatable :: sensor_id(:)
      type(CRTM_ChannelInfo_type) , allocatable :: chinfo(:)       
      type(CRTM_Geometry_type)    , allocatable :: geo(:)          
      type(CRTM_Options_type)     , allocatable :: opt(:)       

      ! Forward declarations
      type(CRTM_Atmosphere_type)  , allocatable :: atm(:)       
      type(CRTM_Surface_type)     , allocatable :: sfc(:)       
      type(CRTM_RTSolution_type)  , allocatable :: rts(:,:)   
   
CONTAINS
   subroutine GPM_GMI_Init(n_absorbers, n_clouds, n_aerosols, gbx)

      ! input variable
      type(cosp_gridbox), intent(in) :: gbx

      integer, intent(in) ::  n_absorbers
      integer, intent(in) ::  n_clouds
      integer, intent(in) ::  n_aerosols

      ! local variables
      integer :: err_stat
      integer :: alloc_stat
      integer :: n
      integer :: n_layers

      !------------
      n_layers = gbx%Nlevels
      n_profiles = gbx%Npoints      
      !------------------------

      ! allocate arrays
      n_sensors = 1
      allocate( sensor_id(n_sensors), chinfo(n_sensors), &
                 STAT = alloc_stat)
      if (alloc_stat /=0) then
         write(*,*) "fail to allocate sensor_id or chinfo"
         stop
      endif

      ! set variable values
      sensor_id = (/'gmi_gpm'/)
      
      err_stat = CRTM_Init(sensor_id, chinfo, &
                           File_Path = dircoef)
      if (err_stat /= SUCCESS) then
         write(*,*) "fail to initialize chinfo"
         STOP
      endif


      ! allocate profile-only arrays
      allocate(geo(n_profiles), opt(n_profiles), atm(n_profiles),&
               sfc(n_profiles) )

      do n=1, n_sensors

         n_channels = CRTM_ChannelInfo_n_Channels( chinfo(n) )
         allocate( rts(n_channels, n_profiles) )
      end do !n_sensors
      
      ! create atmosphere profile. will be moved out later
      call CRTM_Atmosphere_Create( atm, n_layers, n_absorbers, n_clouds,&
                                   n_aerosols)
      if ( ANY (.NOT. CRTM_Atmosphere_Associated(atm) ) )  then
         print *, "fail to create atmosphere profile"
         STOP
      endif

      
      ! Profile and absorber definition
      ! need to loop through profiles
      do n=1, n_profiles
      print *, n
         atm(n)%Climatology         = US_STANDARD_ATMOSPHERE
         atm(n)%Absorber_Id(1:2)    = (/H2O_ID                 ,  O3_ID/)
         atm(n)%Absorber_units(1:2) = (/MASS_MIXING_RATIO_UNITS,  VOLUME_MIXING_RATIO_UNITS/)

         ! Profile data
         ! COSP is from surface to top, but CRTM wants top to surface. Need to flip
         ! the profiles
         atm(n)%Level_Pressure = gbx%pint(n,(n_layers+1):1:-1)
         atm(n)%Pressure       = gbx%p(n,n_layers:1:-1)
         atm(n)%Temperature    = gbx%T(n,n_layers:1:-1)

         ! unit for q in COSP is kg/kg, while CRTM wants g/kg. Need to do unit
         ! conversion
         atm(n)%Absorber(:,1)  = gbx%q(n,n_layers:1:-1) * 1000
         ! units for O3 mixing ratio in COSP is kg/kg, but CRTM needs ppmv. The
         ! conversion is not implemented here yet.
         
         ! assign cloud properties
         atm(n)%Cloud(1)%Type = WATER_CLOUD
         atm(n)%Cloud(1)%Effective_Radius(2:3) = (/20_fp, 12_fp/)
         atm(n)%Cloud(1)%Water_Content(2:3) = (/5_fp, 2_fp/)
         ! assign aerosol properties
         atm(n)%Aerosol(1)%Type = DUST_AEROSOL
         atm(n)%Aerosol(1)%Effective_Radius(3:4) = (/7e-16_fp, 3e-3_fp/)
         atm(n)%Aerosol(1)%Concentration(3:4) = (/2e-18_fp, 1e-5_fp/)
         !!! assign surface properties
         sfc(n)%Land_Coverage = 0.1_fp
         sfc(n)%Water_Coverage = 0.5_fp
         sfc(n)%Snow_Coverage = 0.25_fp
         sfc(n)%Ice_Coverage = 0.15_fp





      end do !n_profiles

   end subroutine gpm_gmi_init

   !!!!!!!!!!!
   subroutine gpm_gmi_run()
      integer :: err_stat
      err_stat = CRTM_Forward(atm, sfc, geo, chInfo, rts)
      print *, err_stat
      call CRTM_RTSolution_Inspect(rts)
   end subroutine gpm_gmi_run

   !!!!!!!!!!
   subroutine gpm_gmi_clean()
      integer :: err_stat
      err_stat = CRTM_Destroy(chinfo)
      call CRTM_Options_Destroy(opt)
      call CRTM_RTSolution_Destroy(rts)
      call CRTM_Atmosphere_Destroy(atm)

   end subroutine gpm_gmi_clean


end module GPM_GMI_MOD



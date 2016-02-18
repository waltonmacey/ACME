! module related to calling CRTM functions
! This file mimics the cosp_rttov.F90 in COSP package
!
! Jan. 13, 2016  Created by Yinghui Lu
include "../cosp_gpm_debugflag.F90"
module GPM_CRTM_mod
   use CRTM_Module
   use MOD_COSP_CONSTANTS, only: I_LSCLIQ, I_LSCICE, I_LSRAIN, I_LSSNOW, &
                                 I_CVCLIQ, I_CVCICE, I_CVRAIN, I_CVSNOW, &
                                 I_LSGRPL
   implicit none
contains
!--------------------
!
! Calling CRTM functions. Calculate one sensor at a time
!
!--------------------
   subroutine crtm_multiprof( &
      n_profiles,   & ! number of profiles
      n_layers,     & ! number of atmosphere layers
      n_absorbers,  & ! number of absorbers
      n_clouds,     & ! number of cloud profiles
      n_aerosols,   & ! number of aerosols
      pint,         & ! pressure between levels [hPa]
      p,            & ! pressure at middle of levels [hPa]
      T,            & ! temperature [Kelvin]
      mr_vapor,     & ! water vapor mixing ratio [g/kg]
      Reff_hydro,   & ! effective radius of hydrometer particle
      water_content,     & ! mixing ratio of hydrometer
      Reff_aerosol, & ! effective radius of aerosol
      mr_aerosol,   & ! mixing ratio of aerosol
      o3,           & ! ozone mixing ratio [ppmv]
      chinfo,       & ! channel info [type: CRTM_ChannelInfo_type][input][1x1 vector]
      sfc,          & ! CRTM surface properties
      latitude,     & ! latitude
      longitude,    & ! longitude
      scan_angle,   & ! sensor scan angle
      zenith_angle, & ! sensor zenith angle
      year,         & ! model year
      month,        & ! model month of year
      day,          & ! model day of month
      tbs           & ! brightness temperature [output]
      )

      ! input and output variables
      integer, intent(in) :: n_profiles       ! scalar
      integer, intent(in) :: n_layers         ! scalar
      integer, intent(in) :: n_absorbers      ! scalar
      integer, intent(in) :: n_clouds         ! scalar
      integer, intent(in) :: n_aerosols       ! scalar
      real, intent(in) :: pint(:,:)           ! [n_profiles x n_layers]
      real, intent(in) :: p(:,:)              ! [n_profiles x n_layers]
      real, intent(in) :: T(:,:)              ! [n_profiles x n_layers]
      real, intent(in) :: mr_vapor(:,:)       ! [n_profiles x n_layers]
      real, intent(in) :: Reff_hydro(:,:,:)   ! [n_profiles x n_layers x n_clouds]
      real, intent(in) :: water_content(:,:,:)     ! [n_profiles x n_layers x n_clouds]
      real, intent(in) :: Reff_aerosol(:,:,:) ! [n_profiles x n_layers x n_aerosols]
      real, intent(in) :: mr_aerosol(:,:,:)   ! [n_profiles x n_layers x n_aerosols]
      real, intent(in) :: o3(:,:)             ! [n_profiles x n_layers]
      type(CRTM_ChannelInfo_type), intent(in) :: chinfo(:)  ! [1 x 1]
      type(CRTM_Surface_type),     intent(in) :: sfc(:)   ! [n_profiles]
      real, intent(in) :: latitude(:)         ! [n_profiles]
      real, intent(in) :: longitude(:)        ! [n_profiles]         
      real, intent(in) :: scan_angle          ! [1]         
      real, intent(in) :: zenith_angle        ! [1]         
      real, intent(in) :: year                ! [1]         
      real, intent(in) :: month               ! [1]         
      real, intent(in) :: day                 ! [1]         
      
      real, intent(out)   :: tbs(:,:)           ! [n_channels x n_profiles]
      !----- parameters -------
      integer, parameter :: n_surfacetypes = 4
      !---------------
      ! local variables
      ! Geometry and option structure. Length: n_profile
      type(CRTM_Geometry_type),   allocatable :: geo(:) ! geometry structure
      type(CRTM_Options_type),    allocatable :: opt(:) ! option structure
      ! RTSolution structure for results. Size: n_channel x n_profile
      type(CRTM_RTSolution_type), allocatable :: rts(:,:) ! RTsolution
      ! variable associated with atmosphere state. length: n_profile
      type(CRTM_Atmosphere_type), allocatable :: atm(:) ! atmosphere

      integer :: n_channels  ! number of channels
      !---------------
      ! temporary variables
      integer :: i_profile
      integer :: tmpint
      integer :: err_stat
      !---------------------------------
      ! Due to the requirement of CRTM, chinfo must be an array. In this
      ! implementation, the chinfo must be a 1x1 array. Check if the input meets
      ! this requirement first.
      if (size(chinfo) .NE. 1) then
         print *, "chinfo must be a 1 x 1 array of type CRTM_ChannelInfo_type"
         stop
      end if
      
      ! allocate structures used in this calculation
      !FIXME: add check for allocation success
      allocate(atm(n_profiles))

      ! create atmosphere profile
      call CRTM_Atmosphere_Create(atm, n_layers, n_absorbers, n_clouds, &
                                  n_aerosols)
      if ( ANY (.NOT. CRTM_Atmosphere_Associated(atm) ) )  then
         print *, "fail to create atmosphere profile"
         STOP
      endif
      !-------------------------------
      ! set values for atmosphere and surface structures
!print *, "in gpm_crtm.F90"
!print *, maxval(water_content)
!print *, maxval(T)
      do i_profile = 1, n_profiles
         atm(i_profile)%Climatology         = US_STANDARD_ATMOSPHERE
         atm(i_profile)%Absorber_Id(1:2)    = (/H2O_ID                 ,  O3_ID/)
         atm(i_profile)%Absorber_units(1:2) = (/MASS_MIXING_RATIO_UNITS,  VOLUME_MIXING_RATIO_UNITS/)

!#ifdef GMI_DEBUG
  if (.false. .and. i_profile == 11 ) then
    print *, "pint"
    print *, pint(i_profile,:)
    print *, "p"
    print *, p(i_profile, :)
    print *, "T"
    print *, T(i_profile, :)
    print *, "mr_vapor"
    print *, mr_vapor(i_profile, :)

    print *, "o3"
    print *, o3(i_profile, :)

    print *, "Reff_hydro"
    print *, Reff_hydro(i_profile, :,:)
    print *, "water_content"
    print *, water_content(i_profile, :,:)


  end if



!#endif


         ! Profile data
         ! COSP is from surface to top, but CRTM wants top to surface. Need to flip
         ! the profiles
         atm(i_profile)%Level_Pressure(0:n_layers) = pint(i_profile,:)
         atm(i_profile)%Pressure       = p(i_profile, :) 
         atm(i_profile)%Temperature    = T(i_profile, :) 


         ! unit for mr_vapor in COSP is kg/kg, while CRTM wants g/kg. Need to do unit
         ! conversion. This unit conversion is done in gpm_crtm_simulator.F90
         atm(i_profile)%Absorber(:,1)  = mr_vapor(i_profile,:) 
         ! units for O3 mixing ratio in COSP is kg/kg, but CRTM needs ppmv. The
         ! conversion is implemented in gpm_crtm_simulator.F90.
         atm(i_profile)%Absorber(:,2)  = o3(i_profile,:)

#ifdef GMI_DEBUG
         atm(i_profile)%Absorber(:,1) = 0
         atm(i_profile)%Absorber(:,2) = 0
#endif

         ! assign cloud properties
         ! For now, the cloud types follow those in COSP. Each type of
         ! hydrometer defined in COSP is set to be a single cloud
         ! TODO: current implementation do not distinguish cloud or clear. May
         ! need to define number of clouds based on clear/cloud
         ! TODO: Make sure the units agree. In CRTM, Effective radius has unit microns,
         ! while water content has unit kg/m^2

         
#ifndef GMI_DEBUG
         tmpint = I_LSCLIQ
         atm(i_profile)%Cloud(1)%Type = WATER_CLOUD
         atm(i_profile)%Cloud(1)%Effective_Radius = Reff_hydro(i_profile,:,tmpint)
         atm(i_profile)%Cloud(1)%Water_Content    = water_content(i_profile,:,tmpint)
         
         tmpint = I_LSCICE
         atm(i_profile)%Cloud(2)%Type = ICE_CLOUD
         atm(i_profile)%Cloud(2)%Effective_Radius = Reff_hydro(i_profile,:,tmpint)
         atm(i_profile)%Cloud(2)%Water_Content    = water_content(i_profile,:,tmpint)
         
         tmpint = I_LSRAIN
         atm(i_profile)%Cloud(3)%Type = RAIN_CLOUD
         atm(i_profile)%Cloud(3)%Effective_Radius = Reff_hydro(i_profile,:,tmpint)
         atm(i_profile)%Cloud(3)%Water_Content    = water_content(i_profile,:,tmpint)
         
         tmpint = I_LSSNOW
         atm(i_profile)%Cloud(4)%Type = SNOW_CLOUD
         atm(i_profile)%Cloud(4)%Effective_Radius = Reff_hydro(i_profile,:,tmpint)
         atm(i_profile)%Cloud(4)%Water_Content    = water_content(i_profile,:,tmpint)
         
         tmpint = I_CVCLIQ
         atm(i_profile)%Cloud(5)%Type = WATER_CLOUD
         atm(i_profile)%Cloud(5)%Effective_Radius = Reff_hydro(i_profile,:,tmpint)
         atm(i_profile)%Cloud(5)%Water_Content    = water_content(i_profile,:,tmpint)
         
         tmpint = I_CVCICE
         atm(i_profile)%Cloud(6)%Type = ICE_CLOUD
         atm(i_profile)%Cloud(6)%Effective_Radius = Reff_hydro(i_profile,:,tmpint)
         atm(i_profile)%Cloud(6)%Water_Content    = water_content(i_profile,:,tmpint)
         
         tmpint = I_CVRAIN
         atm(i_profile)%Cloud(7)%Type = RAIN_CLOUD
         atm(i_profile)%Cloud(7)%Effective_Radius = Reff_hydro(i_profile,:,tmpint)
         atm(i_profile)%Cloud(7)%Water_Content    = water_content(i_profile,:,tmpint)
         
         tmpint = I_CVSNOW
         atm(i_profile)%Cloud(8)%Type = SNOW_CLOUD
         atm(i_profile)%Cloud(8)%Effective_Radius = Reff_hydro(i_profile,:,tmpint)
         atm(i_profile)%Cloud(8)%Water_Content    = water_content(i_profile,:,tmpint)
         
         tmpint = I_LSGRPL
         atm(i_profile)%Cloud(9)%Type = GRAUPEL_CLOUD
         atm(i_profile)%Cloud(9)%Effective_Radius = Reff_hydro(i_profile,:,tmpint)
         atm(i_profile)%Cloud(9)%Water_Content    = water_content(i_profile,:,tmpint)
         
         ! assign aerosol properties
         ! TODO: aerosol properties and land cover are dummy values for now.
         ! Need to get them from the model

         atm(i_profile)%Aerosol(1)%Type = DUST_AEROSOL
         atm(i_profile)%Aerosol(1)%Effective_Radius = Reff_aerosol(i_profile,:,1)
         atm(i_profile)%Aerosol(1)%Concentration    = mr_aerosol(i_profile,:,1)
#endif         
      end do ! i_profile
      !---------------------
      ! geometry and option structures
      n_channels = CRTM_ChannelInfo_n_Channels(chinfo(1))
      ! FIXME: add check for allocation success
      allocate(geo(n_profiles), opt(n_profiles) )
      geo%Latitude  = latitude
      geo%Longitude = longitude
      geo%Sensor_Scan_angle   = scan_angle
      geo%Sensor_Zenith_Angle = zenith_angle
      geo%Year  = year
      geo%Month = month
      geo%Day   = day
      !-----------------
      ! allocate RTSolution structure
      !FIXME: add check for allocation success
      allocate(rts(n_channels, n_profiles) )
      !-----------------
      ! forward calculation
! For debug use
do i_profile = 1, n_profiles
!  print *, "i_profile", i_profile 
!      call CRTM_Atmosphere_Inspect(atm(i_profile:i_profile))
!      call CRTM_Geometry_Inspect(geo(i_profile))
!      call CRTM_ChannelInfo_Inspect(chinfo)

      err_stat = CRTM_Forward(atm(i_profile:i_profile), sfc(i_profile:i_profile), geo(i_profile:i_profile), chinfo(1:1),rts(:,i_profile:i_profile))
      !err_stat = CRTM_Forward(atm, sfc, geo, chinfo(1:1),rts)

!print *,'in gpm_crtm.F90'
!print *, "BT",rts(6,i_profile)%Brightness_temperature
!print *, "SOD", rts(6,i_profile)%SOD, "surfEm", rts(6, i_profile)%Surface_Emissivity, "a"
end do
!call CRTM_RTSolution_Inspect(rts)
      
      if (err_stat .NE. SUCCESS) then
        print *,'CRTM_Forward not successful'

      end if
      ! obtain brightness temperatures from RTSolution structure


      tbs = rts%Brightness_temperature
! For debug use
!print *, rts(1,1)%Brightness_temperature
!print *, tbs

      ! destroy the CRTM structures
      call CRTM_Geometry_destroy(  geo)
      call CRTM_Options_Destroy(   opt)
      call CRTM_RTSolution_Destroy(rts)
      call CRTM_Atmosphere_Destroy(atm)

   end subroutine crtm_multiprof



end module GPM_CRTM_mod

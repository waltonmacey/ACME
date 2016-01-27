! module converts COSP_GRIDBOX data to the ones used by CRTM
!
! Jan. 14, 2016: Created by Yinghui Lu
! Jan. 22, 2016: Move code related to GPM Results to gpm_gmi_result_mod 
!                (Yinghui Lu)


module gpm_crtm_simulator_mod
   use MOD_COSP_TYPES, only: cosp_gridbox
   use CRTM_Module, only: crtm_channelInfo_type, CRTM_ChannelInfo_n_Channels 
   use GPM_CRTM_mod, only: crtm_multiprof 
   use GPM_CRTM_result_mod, only: GPM_CRTM_result_type, GPM_CRTM_result_init, &
                                  GPM_CRTM_result_destroy, GPM_CRTM_result_inquire, &
                                  GPM_CRTM_result_set
   implicit none
   private
   public:: gpm_crtm_simulator_run


contains

!-----------------------
!
! Translate variables in gbx to the ones used by crtm_multiprof and calculate
! scattering results
!
!-----------------------
   subroutine gpm_crtm_simulator_run(gbx,chinfo, y)
      ! Arguments
      type(cosp_gridbox), intent(in)  :: gbx
      type(CRTM_ChannelInfo_type), intent(in) :: chinfo(:) ! [n_sensors]
      type(GPM_CRTM_result_type), intent(out) :: y(:)  ! [n_sensors]
      
      ! Parameters used for converting model variables to that used by crtm
      real, parameter :: eps    =  0.622
      real, parameter :: Mdry   =  28.966
      real, parameter :: Mo3    =  47.9983
      real, parameter :: Mco2   =  44.0096
      real, parameter :: Mch4   =  16.0426
      real, parameter :: Mn2o   =  44.0129
      real, parameter :: Mco    =  28.0102
      
      ! parameter
      integer, parameter :: n_surfacetypes = 4
      ! variables used as inputs to crtm_multiprof()
      integer :: n_profiles       ! scalar
      integer :: n_layers         ! scalar
      integer, parameter :: n_absorbers = 2      ! scalar
      integer, parameter :: n_clouds = 9         ! scalar
      integer, parameter :: n_aerosols = 1       ! scalar
      real, allocatable :: pint(:,:)           ! [n_profiles x n_layers]
      real, allocatable :: p(:,:)              ! [n_profiles x n_layers]
      real, allocatable :: T(:,:)              ! [n_profiles x n_layers]
      real, allocatable :: q(:,:)              ! [n_profiles x n_layers]
      real, allocatable :: Reff_hydro(:,:,:)   ! [n_profiles x n_layers x n_clouds]
      real, allocatable :: mr_hydro(:,:,:)     ! [n_profiles x n_layers x n_clouds]
      real, allocatable :: Reff_aerosol(:,:,:) ! [n_profiles x n_layers x n_aerosols]
      real, allocatable :: mr_aerosol(:,:,:)   ! [n_profiles x n_layers x n_aerosols]
      real, allocatable :: o3(:,:)             ! [n_profiles x n_layers]
      real, allocatable :: surface_type(:,:)   ! [n_profiles x 4]
      real, allocatable :: latitude(:)         ! [n_profiles]
      real, allocatable :: longitude(:)        ! [n_profiles]         
      real :: scan_angle          ! [1]         
      real :: zenith_angle        ! [1]         
      real :: year                ! [1]         
      real :: month               ! [1]         
      real :: day                 ! [1]         
      real, allocatable :: tbs(:,:)
      !---------- local variables -------
      integer :: n_channels
      integer :: n_sensors
      integer :: n
      !--------------
      ! check inputs
      if (size(chinfo) .NE. size(y)) then
         print *, 'chinfo and y must have the same length'
         stop
      endif
      n_sensors = size(chinfo)
      ! obtain variable dimensions from gbx
      n_profiles = gbx%Npoints
      n_layers   = gbx%Nlevels
      
      ! allocate array
      allocate(pint(n_profiles, n_layers+1), &
                  p(n_profiles, n_layers), &
                  T(n_profiles, n_layers), &
                  q(n_profiles, n_layers), &
                 o3(n_profiles, n_layers) ) 

      allocate(Reff_hydro(n_profiles, n_layers, n_clouds), &
                 mr_hydro(n_profiles, n_layers, n_clouds) )

      allocate(Reff_aerosol(n_profiles, n_layers, n_aerosols), &
                 mr_aerosol(n_profiles, n_layers, n_aerosols) )

      allocate(surface_type(n_profiles, n_surfacetypes) )

      allocate(latitude(n_profiles), &
               longitude(n_profiles) )
      !-------------- construct atmosphere profile -----------
      ! COSP structure layers are from bottom to top, while CRTM wants top to
      ! bottom
      pint = gbx%pint(:, (n_layers+1):1:-1)
      p    = gbx%p(   :, n_layers:1:-1)
      T    = gbx%T(   :, n_layers:1:-1)
      q    = gbx%q(   :, n_layers:1:-1)
      o3   = gbx%mr_ozone(  :, n_layers:1:-1)

      Reff_hydro   = gbx%Reff(    :, n_layers:1:-1, :)
      mr_hydro     = gbx%mr_hydro(:, n_layers:1:-1, :)
      !FIXME: add aerosol properties
      Reff_aerosol = 0
      mr_aerosol   = 0
      !FIXME: add surface type properties
      ! now all water surface
      surface_type(:,1) = 0
      surface_type(:,2) = 1
      surface_type(:,3) = 0
      surface_type(:,4) = 0
      !
      latitude  = gbx%latitude
      longitude = gbx%longitude

      ! FIXME: get model year somehow
      year = 2000
      month = 1
      day = 1

  ! ------- calling crtm_multiprof() for each sensor ------
  ! prepare outputs
print *, '--------- debug gpm_crtm_simulator.F90 -------'
  do n=1, n_sensors
     n_channels = CRTM_ChannelInfo_n_Channels(chinfo(n))
     allocate(tbs(n_channels, n_profiles) )
     call GPM_CRTM_result_destroy(y(n))
     call GPM_CRTM_result_init(n_channels, n_profiles, y(n) )
print *,"channel ", n, " has ", n_channels, " channels"
     call crtm_multiprof( &
      n_profiles,   & ! number of profiles
      n_layers,     & ! number of atmosphere layers
      n_absorbers,  & ! number of absorbers
      n_clouds,     & ! number of cloud profiles
      n_aerosols,   & ! number of aerosols
      pint,         & ! pressure between levels [?]
      p,            & ! pressure at middle of levels [?]
      T,            & ! temperature [?]
      q,            & ! water vapor mixing ratio [g/kg]
      Reff_hydro,   & ! effective radius of hydrometer particle
      mr_hydro,     & ! mixing ratio of hydrometer
      Reff_aerosol, & ! effective radius of aerosol
      mr_aerosol,   & ! mixing ratio of aerosol
      o3,           & ! ozone mixing ratio [ppmv]
      chinfo(n:n),  & ! channel info [type: CRTM_ChannelInfo_type][input][1x1 vector]
      surface_type, & ! surface type [land/water/snow/ice coverage]
      latitude,     & ! latitude
      longitude,    & ! longitude
      scan_angle,   & ! sensor scan angle
      zenith_angle, & ! sensor zenith angle
      year,         & ! model year
      month,        & ! model month of year
      day,          & ! model day of month
      tbs           & ! brightness temperature [output]
      )
      call GPM_CRTM_result_set(y(n), tbs)
      deallocate(tbs)
  end do
   end subroutine gpm_crtm_simulator_run


end module gpm_crtm_simulator_mod

! module converts COSP_GRIDBOX data to the ones used by CRTM
!
! Jan. 14, 2016: Created by Yinghui Lu
! Jan. 22, 2016: Move code related to GPM Results to gpm_gmi_result_mod 
!                (Yinghui Lu)


module gpm_crtm_simulator_mod
   use MOD_COSP_TYPES, only: cosp_gridbox, cosp_sghydro
   use CRTM_Module, only: crtm_channelInfo_type, CRTM_ChannelInfo_n_Channels 
   use GPM_CRTM_mod, only: crtm_multiprof 
!   use GPM_CRTM_result_mod, only: GPM_CRTM_result_type, GPM_CRTM_result_init, &
!                                  GPM_CRTM_result_destroy, GPM_CRTM_result_inquire, &
!                                  GPM_CRTM_result_set
   implicit none
   private
   public:: gpm_crtm_simulator_run, GPM_CRTM_result_type, &
            construct_GPM_CRTM_result, &
            free_GPM_CRTM_result,      &
            GPM_CRTM_result_cpsection

   ! User defined type that used to hold CRTM results (TBs)
   type GPM_CRTM_result_type
     integer :: Npoints
     integer :: Ncolumns
     integer :: Nchannels
     real, allocatable, public :: tbs(:,:,:) ! n_channels x n_profiles
   end type GPM_CRTM_result_type

contains
!===========
   subroutine construct_GPM_CRTM_result(Npoints, Ncolumns, Nchannels, x)
     integer, intent(in) :: Npoints
     integer, intent(in) :: Ncolumns
     integer, intent(in) :: Nchannels
     type(GPM_CRTM_result_type), intent(out) :: x

     x%Npoints = Npoints
     x%Ncolumns = Ncolumns
     x%Nchannels = Nchannels

     allocate(x%tbs(Npoints,Ncolumns,Nchannels))
     x%tbs=0.0
   end subroutine construct_GPM_CRTM_result
!============
   subroutine free_GPM_CRTM_result(x)
     type(GPM_CRTM_result_type), intent(inout) :: x
     x%Npoints   = 0
     x%Ncolumns  = 0
     x%Nchannels = 0
     deallocate(x%tbs)
     
   end subroutine free_GPM_CRTM_result
!=============
   subroutine GPM_CRTM_result_cpsection(ix, iy, x, y)
     integer, intent(in) :: ix(2)
     integer, intent(in) :: iy(2)
     type(GPM_CRTM_result_type), intent(in) :: x
     type(GPM_CRTM_result_type), intent(inout) :: y
     y%tbs(iy(1):iy(2), :,:) = x%tbs(ix(1):ix(2), :, :)
   end subroutine GPM_CRTM_result_cpsection

!-----------------------
!
! Translate variables in gbx to the ones used by crtm_multiprof and calculate
! scattering results
!
!-----------------------
   subroutine gpm_crtm_simulator_run(gbx,sgx,chinfo, scan_angle, zenith_angle, gpm_results)
      ! Arguments
      type(cosp_gridbox), intent(in) :: gbx
      type(cosp_sghydro), intent(in) :: sgx
      type(CRTM_ChannelInfo_type), intent(in) :: chinfo(:) ! [n_sensors]
      real, intent(in) :: scan_angle(:)          ! [n_sensors]         
      real, intent(in) :: zenith_angle(:)        ! [n_sensors]         
      type(GPM_CRTM_result_type), intent(inout) :: gpm_results(:)  ! [n_sensors]
      
      ! Parameters used for converting model variables to that used by crtm
      real, parameter :: eps    =  0.622
      real, parameter :: Mdry   =  28.966
      real, parameter :: Mo3    =  47.9983
      real, parameter :: Mco2   =  44.0096
      real, parameter :: Mch4   =  16.0426
      real, parameter :: Mn2o   =  44.0129
      real, parameter :: Mco    =  28.0102
      real, parameter :: g      =  9.81 ! TODO: possible better values for g
      
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
      real, allocatable :: mr_vapor(:,:)       ! [n_profiles x n_layers]
      real, allocatable :: Reff_hydro(:,:,:)   ! [n_profiles x n_layers x n_clouds]
      real, allocatable :: water_content(:,:,:)     ! [n_profiles x n_layers x n_clouds]
      real, allocatable :: Reff_aerosol(:,:,:) ! [n_profiles x n_layers x n_aerosols]
      real, allocatable :: mr_aerosol(:,:,:)   ! [n_profiles x n_layers x n_aerosols]
      real, allocatable :: o3(:,:)             ! [n_profiles x n_layers]
      real, allocatable :: surface_type(:,:)   ! [n_profiles x 4]
      real, allocatable :: latitude(:)         ! [n_profiles]
      real, allocatable :: longitude(:)        ! [n_profiles]         
      real :: year                ! [1]         
      real :: month               ! [1]         
      real :: day                 ! [1]         
      !---------- local variables -------
      integer :: n_channels
      integer :: n_sensors
      integer :: n
      integer :: i_column
      integer :: i_cloud
      integer :: i_channel
      real, allocatable :: tbs(:,:) ! temporary variable to store Tb result
      real, allocatable :: m_layer(:,:) ! mass per unit area of air in each layer [kg/m^2] 
      !--------------
      n_sensors = size(chinfo)
      ! obtain variable dimensions from gbx
      n_profiles = gbx%Npoints
      n_layers   = gbx%Nlevels
      
      ! allocate array
      allocate(pint(n_profiles, n_layers+1), &
                  p(n_profiles, n_layers), &
                  T(n_profiles, n_layers), &
                  mr_vapor(n_profiles, n_layers), &
                 o3(n_profiles, n_layers) ) 

      allocate(Reff_hydro(n_profiles, n_layers, n_clouds), &
                 water_content(n_profiles, n_layers, n_clouds) )

      allocate(Reff_aerosol(n_profiles, n_layers, n_aerosols), &
                 mr_aerosol(n_profiles, n_layers, n_aerosols) )

      allocate(surface_type(n_profiles, n_surfacetypes) )

      allocate(latitude(n_profiles), &
               longitude(n_profiles) )

      allocate(m_layer(n_profiles, n_layers) )
      !-------------- construct atmosphere profile -----------
      ! COSP structure layers are from bottom to top, while CRTM wants top to
      ! bottom
      pint = gbx%pint(:, (n_layers+1):1:-1) / 100  ! convert Pa to hPa
      p    = gbx%p(   :, n_layers:1:-1) / 100      ! convert Pa to hPa
      T    = gbx%T(   :, n_layers:1:-1)
      mr_vapor  = gbx%sh(   :, n_layers:1:-1)*1000 ! convert kg/kg to g/kg
      o3   = gbx%mr_ozone(  :, n_layers:1:-1)* ( Mdry / Mo3 )  * 1e6 

      m_layer = (gbx%pint(:,n_layers:1:-1) - gbx%pint(:,(n_layers+1):2:-1)) /g
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
!print *, "-----debug mr_hydro values -----"
!print *, maxval(sgx%mr_hydro)
!print *, maxval(m_layer)
!print *, minval(m_layer)
!print *, maxval(gbx%pint), minval(gbx%pint)
!print *, gbx%pint(1,(n_layers+1):2:-1)
!print *, gbx%pint(1,n_layers:1:-1)

  ! ------- calling crtm_multiprof() for each sensor ------
  ! prepare outputs
!print *, '--------- debug gpm_crtm_simulator.F90 -------'
!print *, 'sensor scan angle', scan_angle
!print *, 'sensor zenith angle', zenith_angle
!print *, 'n_sensors:', n_sensors
!call CRTM_ChannelInfo_Inspect(chinfo(1))
!call CRTM_ChannelInfo_Inspect(chinfo(2))
  do n=1, n_sensors
     n_channels = CRTM_ChannelInfo_n_Channels(chinfo(n))

     allocate(tbs(n_channels, n_profiles))
      do i_column = 1, gbx%Ncolumns

!print *, "i_sensor, i_column:", n, i_column
      tbs = 0.0
      ! COSP effect radius is in meters. Convert to micros for CRTM
      Reff_hydro   = sgx%Reff(    :, i_column, n_layers:1:-1, :) * 1e6
      ! COSP use mixing ratio for hydrometers [kg/kg]. Convert to water content [kg/m^2] 
      do i_cloud=1,n_clouds
        water_content(:,:,i_cloud)= sgx%mr_hydro(:, i_column, n_layers:1:-1, i_cloud) * m_layer
      end do
      call crtm_multiprof( &
      n_profiles,   & ! number of profiles
      n_layers,     & ! number of atmosphere layers
      n_absorbers,  & ! number of absorbers
      n_clouds,     & ! number of cloud profiles
      n_aerosols,   & ! number of aerosols
      pint,         & ! pressure between levels [hPa]
      p,            & ! pressure at middle of levels [hPa]
      T,            & ! temperature [Kelvin]
      mr_vapor,     & ! water vapor mixing ratio [g/kg]
      Reff_hydro,   & ! effective radius of hydrometer particle [micron]
      water_content,     & ! mixing ratio of hydrometer [kg/m^2]
      Reff_aerosol, & ! effective radius of aerosol [micron]
      mr_aerosol,   & ! mixing ratio of aerosol [kg/m^2]
      o3,           & ! ozone mixing ratio [ppmv]
      chinfo(n:n),  & ! channel info [type: CRTM_ChannelInfo_type][input][1x1 vector]
      surface_type, & ! surface type [land/water/snow/ice coverage]
      latitude,     & ! latitude
      longitude,    & ! longitude
      scan_angle(n),   & ! sensor scan angle
      zenith_angle(n), & ! sensor zenith angle
      year,         & ! model year
      month,        & ! model month of year
      day,          & ! model day of month
      tbs           & ! brightness temperature [output]
      )
!print *, gpm_results(n)%tbs


!print *, 'in gpm_crtm_simulator'
!print *, tbs
      ! Need to reshape the results
      do i_channel = 1, n_channels
        gpm_results(n)%tbs(:,i_column,i_channel) = tbs(i_channel,:)
      end do
      end do ! i_column
!print *, "Maximum TB in sensor ", n, " is ", maxval(gpm_results(n)%tbs )
!print *, "TB(1,1,1): ",gpm_results(n)%tbs(1,1,1)
    deallocate(tbs)
  end do
   end subroutine gpm_crtm_simulator_run


end module gpm_crtm_simulator_mod

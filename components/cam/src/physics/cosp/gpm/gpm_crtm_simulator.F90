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
     real, allocatable, public :: tbs(:,:,:) ! n_profiles x n_columns x n_channels
     
     real, allocatable :: mr_hydro_sg(:,:,:,:)   ! (nprofiles x n_columns x n_layers x n_hydro)
     real, allocatable :: Reff_hydro_sg(:,:,:,:) ! (nprofiles x n_columns x n_layers x n_hydro)
!     real, allocatable :: mr_hydro_sg_2mo(:,:,:,:)   ! (nprofiles x n_columns x n_layers x n_hydro)
!     real, allocatable :: Reff_hydro_sg_2mo(:,:,:,:) ! (nprofiles x n_columns x n_layers x n_hydro)
   end type GPM_CRTM_result_type

contains
!===========
   subroutine construct_GPM_CRTM_result(Npoints, Ncolumns, Nchannels, x, Nlevels)
     integer, intent(in) :: Npoints
     integer, intent(in) :: Ncolumns
     integer, intent(in) :: Nchannels
     type(GPM_CRTM_result_type), intent(out) :: x
     integer, intent(in), optional :: Nlevels

     x%Npoints = Npoints
     x%Ncolumns = Ncolumns
     x%Nchannels = Nchannels

     allocate(x%tbs(Npoints,Ncolumns,Nchannels))
     x%tbs=0.0
     if (present(Nlevels)) then
     allocate(x%mr_hydro_sg      (Npoints,Ncolumns,Nlevels,9))
     allocate(x%Reff_hydro_sg    (Npoints,Ncolumns,Nlevels,9))
!     allocate(x%mr_hydro_sg_2mo  (Npoints,Ncolumns,Nlevels,9))
!     allocate(x%Reff_hydro_sg_2mo(Npoints,Ncolumns,Nlevels,9))

     x%mr_hydro_sg = 0.0
     x%Reff_hydro_sg = 0.0
!     x%mr_hydro_sg_2mo = 0.0
!     x%Reff_hydro_sg_2mo = 0.0
end if
   end subroutine construct_GPM_CRTM_result
!============
   subroutine free_GPM_CRTM_result(x)
     type(GPM_CRTM_result_type), intent(inout) :: x
     x%Npoints   = 0
     x%Ncolumns  = 0
     x%Nchannels = 0
     deallocate(x%tbs)
     if(allocated(x%mr_hydro_sg)) deallocate(x%mr_hydro_sg)
     if(allocated(x%Reff_hydro_sg)) deallocate(x%Reff_hydro_sg)
!     deallocate(x%mr_hydro_sg_2mo)
!     deallocate(x%Reff_hydro_sg_2mo)
   end subroutine free_GPM_CRTM_result
!=============
   subroutine GPM_CRTM_result_cpsection(ix, iy, x, y)
     integer, intent(in) :: ix(2)
     integer, intent(in) :: iy(2)
     type(GPM_CRTM_result_type), intent(in) :: x
     type(GPM_CRTM_result_type), intent(inout) :: y
     y%tbs(iy(1):iy(2), :,:) = x%tbs(ix(1):ix(2), :, :)
     y%mr_hydro_sg      (iy(1):iy(2), :,:,:) = x%mr_hydro_sg      (ix(1):ix(2), :, :,:)
     y%Reff_hydro_sg    (iy(1):iy(2), :,:,:) = x%Reff_hydro_sg    (ix(1):ix(2), :, :,:)
!     y%mr_hydro_sg_2mo  (iy(1):iy(2), :,:,:) = x%mr_hydro_sg_2mo  (ix(1):ix(2), :, :,:)
!     y%Reff_hydro_sg_2mo(iy(1):iy(2), :,:,:) = x%Reff_hydro_sg_2mo(ix(1):ix(2), :, :,:)
   end subroutine GPM_CRTM_result_cpsection

!-----------------------
!
! Translate variables in gbx to the ones used by crtm_multiprof and calculate
! scattering results
!
!-----------------------
   subroutine gpm_crtm_simulator_run(gbx,sghydro,chinfo, scan_angle, zenith_angle, gpm_results, filter_profile, filter_column)
      ! Arguments
      type(cosp_gridbox), intent(in) :: gbx
      type(cosp_sghydro), intent(in) :: sghydro
      type(CRTM_ChannelInfo_type), intent(in) :: chinfo(:) ! [n_sensors]
      real, intent(in) :: scan_angle(:)          ! [n_sensors]         
      real, intent(in) :: zenith_angle(:)        ! [n_sensors]         
      type(GPM_CRTM_result_type), intent(inout) :: gpm_results(:)  ! [n_sensors]
      logical, intent(in), optional :: filter_profile(:) ! [n_profiles] filter to indicate if one profile is calculated
      logical, intent(in), optional :: filter_column(:)  ! [n_columns]  filter to indicate if one column  is calculated
      
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
      integer :: i_profile
      real, allocatable :: tbs(:,:) ! temporary variable to store Tb result
      real, allocatable :: m_layer(:,:) ! mass per unit area of air in each layer [kg/m^2] 
      logical, allocatable :: filter_profile_local(:)
      logical, allocatable :: filter_column_local(:)
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

      allocate(latitude(n_profiles), &
               longitude(n_profiles) )

      allocate(m_layer(n_profiles, n_layers) )

      allocate(filter_profile_local(n_profiles))
      allocate(filter_column_local(gbx%Ncolumns))

      if(present(filter_profile)) then
        filter_profile_local = filter_profile
      else
        filter_profile_local = .true.
      end if
      if(present(filter_column)) then
        filter_column_local = filter_column
      else
        filter_column_local = .true.
      end if
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
      !
      latitude  = gbx%latitude
      longitude = gbx%longitude

      ! FIXME: get model year somehow
      year = 2000
      month = 1
      day = 1
!print *, "-----debug mr_hydro values -----"
!print *, maxval(sghydro%mr_hydro)
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
        if (filter_column_local(i_column) == .true. ) then
!print *, "i_sensor, i_column:", n, i_column
      tbs = 0.0
      ! COSP effect radius is in meters. Convert to micros for CRTM
      Reff_hydro   = sghydro%Reff(    :, i_column, n_layers:1:-1, :) * 1e6
      ! COSP use mixing ratio for hydrometers [kg/kg]. Convert to water content [kg/m^2] 
      do i_cloud=1,n_clouds
        water_content(:,:,i_cloud)= sghydro%mr_hydro(:, i_column, n_layers:1:-1, i_cloud) * m_layer
      end do

if (.false. ) then      
Reff_hydro(:,:,7) = Reff_hydro(:,:,7) * 10
Reff_hydro(:,:,8) = Reff_hydro(:,:,8) * 10

! turn off the clouds
Reff_hydro(:,:,1) = 0
Reff_hydro(:,:,2) = 0
Reff_hydro(:,:,5) = 0
Reff_hydro(:,:,6) = 0
      
water_content(:,:,1) = 0
water_content(:,:,2) = 0
water_content(:,:,5) = 0
water_content(:,:,6) = 0
end if


if (.false.) then
water_content = 0
Reff_hydro = 0

if (i_column <= 8 ) then
water_content(:,:,i_column)= sghydro%mr_hydro(:, 1, n_layers:1:-1, i_column) * m_layer
Reff_hydro(:,:,i_column)   = sghydro%Reff(    :, 1, n_layers:1:-1, i_column) * 1e6
elseif (i_column == 10) then
water_content(:,:,8)= sghydro%mr_hydro(:, 1, n_layers:1:-1, 8) * m_layer
Reff_hydro(:,:,8)   = sghydro%Reff(    :, 1, n_layers:1:-1, 8) * 1e6 * 10
elseif (i_column == 11) then
water_content(:,:,8)= sghydro%mr_hydro(:, 1, n_layers:1:-1, 8) * m_layer
Reff_hydro(:,:,8)   = 3000
elseif (i_column == 12) then
water_content(:,:,8)= sghydro%mr_hydro(:, 1, n_layers:1:-1, 8) * m_layer
Reff_hydro(:,:,8)   = 10000
elseif (i_column == 13) then
water_content(:,:,8)= sghydro%mr_hydro(:, 1, n_layers:1:-1, 8) * m_layer * 10
Reff_hydro(:,:,8)   = 10000
elseif (i_column == 14) then
water_content(:,:,8)= sghydro%mr_hydro(:, 1, n_layers:1:-1, 8) * m_layer * 100
Reff_hydro(:,:,8)   = 10000
elseif (i_column == 15) then
water_content(:,:,8)= sghydro%mr_hydro(:, 1, n_layers:1:-1, 8) * m_layer * 1000
Reff_hydro(:,:,8)   = 10000
elseif (i_column == 16) then
water_content(:,:,7)= sghydro%mr_hydro(:, 1, n_layers:1:-1, 7) * m_layer
Reff_hydro(:,:,7)   = sghydro%Reff(    :, 1, n_layers:1:-1, 7) * 1e6
water_content(:,:,8)= sghydro%mr_hydro(:, 1, n_layers:1:-1, 8) * m_layer * 10
Reff_hydro(:,:,8)   = 10000
elseif (i_column == 17) then
water_content(:,:,7)= sghydro%mr_hydro(:, 1, n_layers:1:-1, 7) * m_layer
Reff_hydro(:,:,7)   = sghydro%Reff(    :, 1, n_layers:1:-1, 7) * 1e6
water_content(:,:,8)= sghydro%mr_hydro(:, 1, n_layers:1:-1, 8) * m_layer * 100
Reff_hydro(:,:,8)   = 10000
elseif (i_column == 18) then
water_content(:,:,7)= sghydro%mr_hydro(:, 1, n_layers:1:-1, 7) * m_layer
Reff_hydro(:,:,7)   = sghydro%Reff(    :, 1, n_layers:1:-1, 7) * 1e6
water_content(:,:,8)= sghydro%mr_hydro(:, 1, n_layers:1:-1, 8) * m_layer * 1000
Reff_hydro(:,:,8)   = 10000
end if
mr_vapor = 0
end if

print *, 'icolumn=', i_column
!if ( i_column >= 2 .AND. i_column <=10) then
!      do i_cloud=1,n_clouds
!      if (i_cloud == 3 .OR. i_cloud==4 .OR. i_cloud==7 .OR. i_cloud == 8) then
!        water_content(:,:,i_cloud)= sghydro%mr_hydro(:, 1, n_layers:1:-1, i_cloud) * m_layer * 2**(i_column-1)
!      end if
!      end do
!end if 
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
      gbx%gpmsurface,  & ! CRTM surface properties
      latitude,     & ! latitude
      longitude,    & ! longitude
      scan_angle(n),   & ! sensor scan angle
      zenith_angle(n), & ! sensor zenith angle
      year,         & ! model year
      month,        & ! model month of year
      day,          & ! model day of month
      tbs,          & ! brightness temperature [output]
      filter_profile_local  & ! filter indicating if one profile is calculated
      )
!print *, gpm_results(n)%tbs
if (.false.) then
do i_profile = 1, n_profiles
if ( i_column == 1) then
  print *, "latitude ", latitude(i_profile), "longitude", longitude(i_profile)
  print *, "maximum water_content ",max( maxval(water_content(:,:,3:4)),  maxval(water_content(:,:,7:8)) )
  write(*,*), water_content(i_profile, :, 3)
  print *, "---"
  write(*,*), water_content(i_profile, :, 4)
  print *, "---"
  write(*,*), water_content(i_profile, :, 7)
  print *, "---"
  write(*,*), water_content(i_profile, :, 8)

  print *, "Reff "
  write(*,*), Reff_hydro(i_profile, :, 3)
  print *, "---"
  write(*,*), Reff_hydro(i_profile, :, 4)
  print *, "---"
  write(*,*), Reff_hydro(i_profile, :, 7)
  print *, "---"
  write(*,*), Reff_hydro(i_profile, :, 8)

end if
end do
end if
!print *, 'in gpm_crtm_simulator'
!print *, tbs
      ! Need to reshape the results
      do i_channel = 1, n_channels
        gpm_results(n)%tbs(:,i_column,i_channel) = tbs(i_channel,:)
      end do

      
if (.false.) then    
if (n_channels == 9) then
gpm_results(n)%tbs(:,1:n_layers,1) = water_content(:,:,3)
gpm_results(n)%tbs(:,1:n_layers,2) = water_content(:,:,4)
gpm_results(n)%tbs(:,1:n_layers,3) = water_content(:,:,7)
gpm_results(n)%tbs(:,1:n_layers,4) = water_content(:,:,8)
end if 

if (n_channels ==4) then
gpm_results(n)%tbs(:,1:n_layers,1) = Reff_hydro(:,:,3)
gpm_results(n)%tbs(:,1:n_layers,2) = Reff_hydro(:,:,4)
gpm_results(n)%tbs(:,1:n_layers,3) = Reff_hydro(:,:,7)
gpm_results(n)%tbs(:,1:n_layers,4) = Reff_hydro(:,:,8)
endif
end if
        end if ! filter_column
      end do ! i_column
!print *, "Maximum TB in sensor ", n, " is ", maxval(gpm_results(n)%tbs )
!print *, "TB(1,1,1): ",gpm_results(n)%tbs(1,1,1)
    deallocate(tbs)
  end do
   end subroutine gpm_crtm_simulator_run


end module gpm_crtm_simulator_mod

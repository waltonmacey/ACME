#define GPM_KU
#define GPM_KA
#define GPM_GMI2


program gmioffline
  use netcdf
  use shr_kind_mod, only: r8 => shr_kind_r8
  use mod_cosp_types
  use GPM_CRTM_sensor_mod, only: GPM_CRTM_sensor_add, GPM_CRTM_sensor_init, &
                                 chinfo_list, sensor_scan_angle_list, &
                                 sensor_zenith_angle_list, n_sensors
  USE gpm_crtm_simulator_mod, only: gpm_crtm_simulator_run, &
                                    GPM_CRTM_result_type,   &
                                    construct_GPM_CRTM_result, &
                                    free_GPM_CRTM_result
  implicit none

  ! local variables
  character (len=*), parameter :: inputfile = &
  "/global/cscratch1/sd/yxl232/acme_scratch/gpm_capt.v0.1.0.ne30_ne30/run/gpm_capt.v0.1.0.ne30_ne30.cam.h0.2011-01-01-00000.nc" 
  character (len=*), parameter :: outputfile = &
  "/global/cscratch1/sd/yxl232/acme_scratch/gpm_capt.v0.1.0.ne30_ne30/run/gpm_capt.v0.1.0.ne30_ne30.cam.out.2011-01-01-00000.nc" 
  integer                      :: ncol
  integer                      :: cosp_scol
  integer                      :: lev
  integer                      :: time
  integer, parameter :: nch=9
  ! used for loops
  integer                      :: itime , icol, iscol
  integer :: ncid, varid, dimid
  integer :: ncid_out, ncol_id, scol_id, ch_id, time_id
  integer :: lat_varid, lon_varid, tmi_varid, tmi_varid2, scol_varid, ch_varid
  real(r8), allocatable :: tmp2d(:), tmp3d(:,:), tmp4d(:,:,:)
  logical, allocatable  :: filter_profile(:)
  logical, allocatable  :: filter_column(:)
  real(r8)  :: tmp

  type(cosp_gridbox) :: gbx 
  type(cosp_sghydro) :: sghydro
  type(GPM_CRTM_result_type) :: sggpmgmi(1)
  type(GPM_CRTM_result_type) :: sggpmgmi2(1)
  
  ! prepare CRTM sensors
  call GPM_CRTM_sensor_add('trmm-tmi')
  call GPM_CRTM_sensor_init()

  ! open NC file for read
  call check( nf90_open(inputfile, NF90_NOWRITE, ncid) )

  ! inquire dimension
  call check( nf90_inq_dimid(ncid, "ncol", dimid) )
  call check( nf90_inquire_dimension(ncid, dimid, len=ncol) )
  call check( nf90_inq_dimid(ncid, "cosp_scol", dimid) )
  call check( nf90_inquire_dimension(ncid, dimid, len=cosp_scol) )
  call check( nf90_inq_dimid(ncid, "lev", dimid) )
  call check( nf90_inquire_dimension(ncid, dimid, len=lev) )
  call check( nf90_inq_dimid(ncid, "time", dimid) )
  call check( nf90_inquire_dimension(ncid, dimid, len=time) )

  ! allocate arrays
  call construct_cosp_gridbox_lite(ncol, lev, cosp_scol, gbx)
  call construct_cosp_sghydro(ncol, cosp_scol, lev, nch, sghydro)
  call construct_gpm_crtm_result(ncol, cosp_scol, nch, sggpmgmi(1))
  call construct_gpm_crtm_result(ncol, cosp_scol, nch, sggpmgmi2(1))
  allocate(tmp2d(ncol) )
  allocate(tmp3d(ncol, lev) )
  allocate(tmp4d(ncol, cosp_scol, lev ) )
  allocate(filter_profile(ncol) )
  allocate(filter_column(cosp_scol))

  ! load variables
  call check( nf90_inq_varid(ncid, "lat",     varid) )
  call check( nf90_get_var(ncid, varid, gbx%latitude )) 
  call check( nf90_inq_varid(ncid, "lon",     varid) )
  call check( nf90_get_var(ncid, varid, gbx%longitude))

  filter_profile = .false. 
  filter_profile = gbx%latitude   > -40 .AND. gbx%latitude  < 40  .AND. &
           gbx%longitude  > 140 .AND. gbx%longitude < 220
  filter_profile = .true.

  filter_column = .false.
  filter_column(1:2) = .true.
filter_column = .true.
  ! prepare NC file for output
  call check( nf90_create(outputfile, nf90_clobber, ncid_out))
  ! define the dimensions
  call check( nf90_def_dim(ncid_out, "ncol",      ncol,      ncol_id))
  call check( nf90_def_dim(ncid_out, "cosp_scol", cosp_scol, scol_id))
  call check( nf90_def_dim(ncid_out, "tmi_chidx", nch,         ch_id))
  call check( nf90_def_dim(ncid_out, "time",      NF90_UNLIMITED,   time_id))
  ! define dimension variables
  call check( nf90_def_var(ncid_out, "lat", NF90_REAL, ncol_id, lat_varid ))
  call check( nf90_put_att(ncid_out, lat_varid, "units", "degrees_north"))
  call check( nf90_def_var(ncid_out, "lon", NF90_REAL, ncol_id, lon_varid))
  call check( nf90_put_att(ncid_out, lon_varid, "units", "degrees_east"))
  call check( nf90_def_var(ncid_out, "cosp_scol", NF90_REAL, scol_id, scol_varid))
  call check( nf90_put_att(ncid_out, scol_varid, "units", "N/A"))
  call check( nf90_def_var(ncid_out, "tmi_chidx", NF90_REAL, ch_id, ch_varid))
  call check( nf90_put_att(ncid_out, ch_varid, "units", "N/A"))

  ! define variables for TBs
  call check( nf90_def_var(ncid_out, "TRMM_TB_CS", NF90_REAL, (/ncol_id, scol_id, ch_id, time_id/), tmi_varid))
  call check( nf90_put_att(ncid_out, tmi_varid, "units", "K" ) )
  call check( nf90_def_var(ncid_out, "TRMM_TB_CS2", NF90_REAL, (/ncol_id, scol_id, ch_id, time_id/), tmi_varid2))
  call check( nf90_put_att(ncid_out, tmi_varid2, "units", "K" ) )
 

  call check( nf90_enddef(ncid_out))


  call check( nf90_put_var(ncid_out, lon_varid, gbx%longitude))
  call check( nf90_put_var(ncid_out, lat_varid, gbx%latitude))
  call check( nf90_put_var(ncid_out, scol_varid, (/1,2,3,4,5,6,7,8,9,10/) ) )
  call check( nf90_put_var(ncid_out, ch_varid,   (/1,2,3,4,5,6,7,8,9/) ) )
  

  ! loop through time 
  do itime = 3,3! time 
    ! surface variable
    call check( nf90_inq_varid(ncid, "PS_COSP", varid) )
    call check( nf90_get_var(ncid, varid, gbx%psfc,      start=(/1, itime/), count=(/ncol, 1/)  )) 
    call check( nf90_inq_varid(ncid, "TS_COSP", varid) )
    call check( nf90_get_var(ncid, varid, gbx%skt,       start=(/1, itime/), count=(/ncol, 1/)  )) 
    call check( nf90_inq_varid(ncid, "U_COSP",  varid) )
    call check( nf90_get_var(ncid, varid, gbx%u_wind,    start=(/1, itime/), count=(/ncol, 1/)  )) 
    call check( nf90_inq_varid(ncid, "V_COSP",  varid) )
    call check( nf90_get_var(ncid, varid, gbx%v_wind,    start=(/1, itime/), count=(/ncol, 1/)  )) 
    call check( nf90_inq_varid(ncid, "CRTM_LAND_FRAC", varid) )
    call check( nf90_get_var(ncid, varid, gbx%gpmsurface%land_coverage,         start=(/1, itime/), count=(/ncol, 1/)  )) 
    call check( nf90_inq_varid(ncid, "CRTM_WATER_FRAC", varid) )
    call check( nf90_get_var(ncid, varid, gbx%gpmsurface%water_coverage,        start=(/1, itime/), count=(/ncol, 1/)  )) 
    call check( nf90_inq_varid(ncid, "CRTM_SOIL_MOIST", varid) )
    call check( nf90_get_var(ncid, varid, gbx%gpmsurface%Soil_Moisture_Content, start=(/1, itime/), count=(/ncol, 1/)  )) 
    call check( nf90_inq_varid(ncid, "CRTM_VEG_FRAC",   varid) )
    call check( nf90_get_var(ncid, varid, gbx%gpmsurface%Vegetation_Fraction,   start=(/1, itime/), count=(/ncol, 1/)  )) 
    call check( nf90_inq_varid(ncid, "CRTM_SOIL_TEMP",  varid) )
    call check( nf90_get_var(ncid, varid, gbx%gpmsurface%Soil_Temperature,      start=(/1, itime/), count=(/ncol, 1/)  )) 
    call check( nf90_inq_varid(ncid, "CRTM_LAI",        varid) )
    call check( nf90_get_var(ncid, varid, gbx%gpmsurface%LAI,                   start=(/1, itime/), count=(/ncol, 1/)  )) 
    call check( nf90_inq_varid(ncid, "CRTM_SOIL_TYPE",  varid) )
    call check( nf90_get_var(ncid, varid, gbx%gpmsurface%Soil_Type,             start=(/1, itime/), count=(/ncol, 1/)  )) 
    call check( nf90_inq_varid(ncid, "CRTM_VEG_TYPE",   varid) )
    call check( nf90_get_var(ncid, varid, gbx%gpmsurface%Vegetation_Type,       start=(/1, itime/), count=(/ncol, 1/)  )) 
    call check( nf90_inq_varid(ncid, "CRTM_WIND_SPEED", varid) )
    call check( nf90_get_var(ncid, varid, gbx%gpmsurface%Wind_Speed,            start=(/1, itime/), count=(/ncol, 1/)  )) 
    call check( nf90_inq_varid(ncid, "CRTM_WIND_DIR",   varid) )
    call check( nf90_get_var(ncid, varid, gbx%gpmsurface%Wind_Direction,        start=(/1, itime/), count=(/ncol, 1/)  )) 

    call check( nf90_inq_varid(ncid, "PINT0_COSP",   varid) )
    call check( nf90_get_var(ncid, varid, tmp2d,        start=(/1, itime/), count=(/ncol, 1/)  )) 



  do icol = 1, ncol
  tmp = gbx%gpmsurface(icol)%land_coverage + gbx%gpmsurface(icol)%water_coverage
    if ( tmp > 1 .AND. tmp<1.000001 ) then
     gbx%gpmsurface(icol)%land_coverage   = gbx%gpmsurface(icol)%land_coverage   /  tmp 
     gbx%gpmsurface(icol)%water_coverage  = gbx%gpmsurface(icol)%water_coverage  /  tmp
    end if
  end do 
    gbx%gpmsurface%ice_coverage   = 1 - gbx%gpmsurface%land_coverage -  gbx%gpmsurface%water_coverage

    gbx%gpmsurface%land_temperature  = gbx%skt
    gbx%gpmsurface%water_temperature = gbx%skt
    gbx%gpmsurface%snow_temperature  = gbx%skt
    gbx%gpmsurface%ice_temperature   = gbx%skt



    ! level variable
    call check( nf90_inq_varid(ncid, "P_COSP",   varid) )
    call check( nf90_get_var(ncid, varid, tmp3d,        start=(/1, 1, itime/), count=(/ncol,lev, 1/)  )) 
    gbx%p        = tmp3d(:, lev:1:-1) 
    call check( nf90_inq_varid(ncid, "T_COSP",   varid) )
    call check( nf90_get_var(ncid, varid, tmp3d,        start=(/1, 1, itime/), count=(/ncol,lev, 1/)  )) 
    gbx%T        = tmp3d(:, lev:1:-1) 
    call check( nf90_inq_varid(ncid, "Q_COSP",   varid) )
    call check( nf90_get_var(ncid, varid, tmp3d,        start=(/1, 1, itime/), count=(/ncol,lev, 1/)  )) 
    gbx%sh       = tmp3d(:, lev:1:-1) 
    call check( nf90_inq_varid(ncid, "PH_COSP",   varid) )
    call check( nf90_get_var(ncid, varid, tmp3d,        start=(/1, 1, itime/), count=(/ncol,lev, 1/)  )) 
    gbx%ph       = tmp3d(:, lev:1:-1) 
    call check( nf90_inq_varid(ncid, "O3_COSP",   varid) )
    call check( nf90_get_var(ncid, varid, tmp3d,        start=(/1, 1, itime/), count=(/ncol,lev, 1/)  )) 
    gbx%mr_ozone = tmp3d(:, lev:1:-1) 

    gbx%pint(:,1:lev) = gbx%ph(:,1:lev)
    gbx%pint(:,lev+1) = tmp2d(:)

    ! sub-grid hydrometer
    call check( nf90_inq_varid(ncid, "CRTM_MR_HYDRO_1",   varid) )
    call check( nf90_get_var(ncid, varid, tmp4d,        start=(/1, 1, 1, itime/), count=(/ncol, cosp_scol, lev, 1/)  )) 
    sghydro%mr_hydro(:,:,:,1)       = tmp4d(:, :, lev:1:-1) 
    call check( nf90_inq_varid(ncid, "CRTM_MR_HYDRO_2",   varid) )
    call check( nf90_get_var(ncid, varid, tmp4d,        start=(/1, 1, 1, itime/), count=(/ncol, cosp_scol, lev, 1/)  )) 
    sghydro%mr_hydro(:,:,:,2)       = tmp4d(:, :, lev:1:-1) 
    call check( nf90_inq_varid(ncid, "CRTM_MR_HYDRO_3",   varid) )
    call check( nf90_get_var(ncid, varid, tmp4d,        start=(/1, 1, 1, itime/), count=(/ncol, cosp_scol, lev, 1/)  )) 
    sghydro%mr_hydro(:,:,:,3)       = tmp4d(:, :, lev:1:-1) 
    call check( nf90_inq_varid(ncid, "CRTM_MR_HYDRO_4",   varid) )
    call check( nf90_get_var(ncid, varid, tmp4d,        start=(/1, 1, 1, itime/), count=(/ncol, cosp_scol, lev, 1/)  )) 
    sghydro%mr_hydro(:,:,:,4)       = tmp4d(:, :, lev:1:-1) 
    call check( nf90_inq_varid(ncid, "CRTM_MR_HYDRO_5",   varid) )
    call check( nf90_get_var(ncid, varid, tmp4d,        start=(/1, 1, 1, itime/), count=(/ncol, cosp_scol, lev, 1/)  )) 
    sghydro%mr_hydro(:,:,:,5)       = tmp4d(:, :, lev:1:-1) 
    call check( nf90_inq_varid(ncid, "CRTM_MR_HYDRO_6",   varid) )
    call check( nf90_get_var(ncid, varid, tmp4d,        start=(/1, 1, 1, itime/), count=(/ncol, cosp_scol, lev, 1/)  )) 
    sghydro%mr_hydro(:,:,:,6)       = tmp4d(:, :, lev:1:-1) 
    call check( nf90_inq_varid(ncid, "CRTM_MR_HYDRO_7",   varid) )
    call check( nf90_get_var(ncid, varid, tmp4d,        start=(/1, 1, 1, itime/), count=(/ncol, cosp_scol, lev, 1/)  )) 
    sghydro%mr_hydro(:,:,:,7)       = tmp4d(:, :, lev:1:-1) 
    call check( nf90_inq_varid(ncid, "CRTM_MR_HYDRO_8",   varid) )
    call check( nf90_get_var(ncid, varid, tmp4d,        start=(/1, 1, 1, itime/), count=(/ncol, cosp_scol, lev, 1/)  )) 
    sghydro%mr_hydro(:,:,:,8)       = tmp4d(:, :, lev:1:-1) 
    call check( nf90_inq_varid(ncid, "CRTM_MR_HYDRO_9",   varid) )
    call check( nf90_get_var(ncid, varid, tmp4d,        start=(/1, 1, 1, itime/), count=(/ncol, cosp_scol, lev, 1/)  )) 
    sghydro%mr_hydro(:,:,:,9)       = tmp4d(:, :, lev:1:-1) 

    call check( nf90_inq_varid(ncid, "CRTM_REFF_HYDRO_1",   varid) )
    call check( nf90_get_var(ncid, varid, tmp4d,        start=(/1, 1, 1, itime/), count=(/ncol, cosp_scol, lev, 1/)  )) 
    sghydro%Reff(:,:,:,1)       = tmp4d(:, :, lev:1:-1) 
    call check( nf90_inq_varid(ncid, "CRTM_REFF_HYDRO_2",   varid) )
    call check( nf90_get_var(ncid, varid, tmp4d,        start=(/1, 1, 1, itime/), count=(/ncol, cosp_scol, lev, 1/)  )) 
    sghydro%Reff(:,:,:,2)       = tmp4d(:, :, lev:1:-1) 
    call check( nf90_inq_varid(ncid, "CRTM_REFF_HYDRO_3",   varid) )
    call check( nf90_get_var(ncid, varid, tmp4d,        start=(/1, 1, 1, itime/), count=(/ncol, cosp_scol, lev, 1/)  )) 
    sghydro%Reff(:,:,:,3)       = tmp4d(:, :, lev:1:-1) 
    call check( nf90_inq_varid(ncid, "CRTM_REFF_HYDRO_4",   varid) )
    call check( nf90_get_var(ncid, varid, tmp4d,        start=(/1, 1, 1, itime/), count=(/ncol, cosp_scol, lev, 1/)  )) 
    sghydro%Reff(:,:,:,4)       = tmp4d(:, :, lev:1:-1) 
    call check( nf90_inq_varid(ncid, "CRTM_REFF_HYDRO_5",   varid) )
    call check( nf90_get_var(ncid, varid, tmp4d,        start=(/1, 1, 1, itime/), count=(/ncol, cosp_scol, lev, 1/)  )) 
    sghydro%Reff(:,:,:,5)       = tmp4d(:, :, lev:1:-1) 
    call check( nf90_inq_varid(ncid, "CRTM_REFF_HYDRO_6",   varid) )
    call check( nf90_get_var(ncid, varid, tmp4d,        start=(/1, 1, 1, itime/), count=(/ncol, cosp_scol, lev, 1/)  )) 
    sghydro%Reff(:,:,:,6)       = tmp4d(:, :, lev:1:-1) 
    call check( nf90_inq_varid(ncid, "CRTM_REFF_HYDRO_7",   varid) )
    call check( nf90_get_var(ncid, varid, tmp4d,        start=(/1, 1, 1, itime/), count=(/ncol, cosp_scol, lev, 1/)  )) 
    sghydro%Reff(:,:,:,7)       = tmp4d(:, :, lev:1:-1) 
    call check( nf90_inq_varid(ncid, "CRTM_REFF_HYDRO_8",   varid) )
    call check( nf90_get_var(ncid, varid, tmp4d,        start=(/1, 1, 1, itime/), count=(/ncol, cosp_scol, lev, 1/)  )) 
    sghydro%Reff(:,:,:,8)       = tmp4d(:, :, lev:1:-1) 
    call check( nf90_inq_varid(ncid, "CRTM_REFF_HYDRO_9",   varid) )
    call check( nf90_get_var(ncid, varid, tmp4d,        start=(/1, 1, 1, itime/), count=(/ncol, cosp_scol, lev, 1/)  )) 
    sghydro%Reff(:,:,:,9)       = tmp4d(:, :, lev:1:-1) 

  print *, n_sensors, ncol, cosp_scol, lev

! with water vapor on
! first turn off all hydrometers
  sghydro%mr_hydro(:,2:cosp_scol,:,:) = 0
  sghydro%Reff    (:,2:cosp_scol,:,:) = 0
! then turn on one hydrometer for one column
  do iscol = 1,9
  sghydro%mr_hydro(:,iscol+1,:,iscol) = sghydro%mr_hydro(:,1,:,iscol)
  end do


  call gpm_crtm_simulator_run(gbx,sghydro,chinfo_list(1:n_sensors),  &
                              sensor_scan_angle_list(1:n_sensors),   &
                              sensor_zenith_angle_list(1:n_sensors), &
                              sggpmgmi,                              &
                              filter_profile=filter_profile,         &
                              filter_column =filter_column  ) 
! turn off water vapor
  gbx%sh = 0 
  call gpm_crtm_simulator_run(gbx,sghydro,chinfo_list(1:n_sensors),  &
                              sensor_scan_angle_list(1:n_sensors),   &
                              sensor_zenith_angle_list(1:n_sensors), &
                              sggpmgmi2,                             &
                              filter_profile=filter_profile,         &
                              filter_column =filter_column  ) 



  call check( nf90_put_var(ncid_out, tmi_varid,   sggpmgmi(1)%tbs ,  start=(/1, 1, 1, itime/), count=(/ncol, cosp_scol, nch, 1/) ) )
  call check( nf90_put_var(ncid_out, tmi_varid2, sggpmgmi2(1)%tbs ,  start=(/1, 1, 1, itime/), count=(/ncol, cosp_scol, nch, 1/) ) )



  end do
  call free_cosp_gridbox(gbx)
  call free_cosp_sghydro(sghydro)
  call free_gpm_crtm_result(sggpmgmi(1))
  call free_gpm_crtm_result(sggpmgmi2(1))
  call check( nf90_close(ncid))
  call check( nf90_close(ncid_out))

contains
  subroutine check(status)
    integer, intent ( in) :: status
      
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  

  subroutine construct_cosp_gridbox_lite(Npoints, Nlevels, ncolumns, gbx)
  use mod_cosp_types
    integer, intent(in) :: Npoints
    integer, intent(in) :: Nlevels
    integer, intent(in) :: Ncolumns
    type(cosp_gridbox), intent(out) :: gbx
  

   real(r8), parameter :: time = 1.0_r8                 ! time ! Time since start of run [days], set to 1 bc running over single CAM timestep
   real(r8), parameter :: time_bnds(2)=(/0.5_r8,1.5_r8/)        ! time_bnds ! Time boundaries - new in cosp v1.3, set following cosp_test.f90 line 121


   ! namelist variables for COSP input related to radar simulator
   real(r8) :: radar_freq = 94.0_r8                     ! CloudSat radar frequency (GHz) (94.0)
   integer :: surface_radar = 0                         ! surface=1, spaceborne=0 (0)
   integer :: use_mie_tables = 0                        ! use a precomputed lookup table? yes=1,no=0 (0)
   integer :: use_gas_abs = 1                           ! include gaseous absorption? yes=1,no=0 (1)
   integer :: do_ray = 0                                ! calculate/output Rayleigh refl=1, not=0 (0)
   integer :: melt_lay = 0                              ! melting layer model off=0, on=1 (0)
   real(r8) :: k2 = -1                                  ! |K|^2, -1=use frequency dependent default (-1)
#ifdef GPM_KU   
   ! namelist variables for COSP input related to GPM Ku radar simulator                                                !+YLu
   real(r8) :: gpmku_radar_freq = 13.6_r8                     ! GPM Ku band radar frequency (GHz) (13.6)               !+YLu
   integer :: gpmku_surface_radar = 0                         ! surface=1, spaceborne=0 (0)                            !+YLu
   integer :: gpmku_use_mie_tables = 0                        ! use a precomputed lookup table? yes=1,no=0 (0)         !+YLu
   integer :: gpmku_use_gas_abs = 1                           ! include gaseous absorption? yes=1,no=0 (1)             !+YLu
   integer :: gpmku_do_ray = 0                                ! calculate/output Rayleigh refl=1, not=0 (0)            !+YLu
   integer :: gpmku_melt_lay = 0                              ! melting layer model off=0, on=1 (0)                    !+YLu
   real(r8) :: gpmku_k2 = -1                                  ! |K|^2, -1=use frequency dependent default (-1)         !+YLu
#endif
#ifdef GPM_KA
   ! namelist variables for COSP input related to GPM Ka radar simulator                                                !+YLu
   real(r8) :: gpmka_radar_freq = 94.0_r8 !35.5_r8                     ! GPM Ka band radar frequency (GHz) (35.5)               !+YLu
   integer :: gpmka_surface_radar = 0                         ! surface=1, spaceborne=0 (0)                            !+YLu
   integer :: gpmka_use_mie_tables = 0                        ! use a precomputed lookup table? yes=1,no=0 (0)         !+YLu
   integer :: gpmka_use_gas_abs = 1                           ! include gaseous absorption? yes=1,no=0 (1)             !+YLu
   integer :: gpmka_do_ray = 0                                ! calculate/output Rayleigh refl=1, not=0 (0)            !+YLu
   integer :: gpmka_melt_lay = 0                              ! melting layer model off=0, on=1 (0)                    !+YLu
   real(r8) :: gpmka_k2 = -1                                  ! |K|^2, -1=use frequency dependent default (-1)         !+YLu
#endif  
   
   ! namelist variables for COSP input related to lidar simulator
   integer, parameter :: Nprmts_max_hydro = 12          ! Max # params for hydrometeor size distributions (12)
   integer, parameter :: Naero = 1                      ! Number of aerosol species (Not used) (1)
   integer, parameter :: Nprmts_max_aero = 1            ! Max # params for aerosol size distributions (not used) (1)
   integer :: lidar_ice_type = 0                        ! Ice particle shape in lidar calculations 
                                                        ! (0=ice-spheres ; 1=ice-non-spherical) (0)
   integer, parameter :: overlap = 3                    ! overlap type: 1=max, 2=rand, 3=max/rand (3)

   !! namelist variables for COSP input related to ISCCP simulator
   integer :: isccp_topheight = 1                       ! 1 = adjust top height using both a computed infrared 
                                                        ! brightness temperature and the visible
                                                        ! optical depth to adjust cloud top pressure. 
                                                        ! Note that this calculation is most appropriate to compare
                                                        ! to ISCCP data during sunlit hours.
                                                        ! 2 = do not adjust top height, that is cloud top pressure 
                                                        ! is the actual cloud top pressure in the model
                                                        ! 3 = adjust top height using only the computed infrared 
                                                        ! brightness temperature. Note that this calculation is most 
                                                        ! appropriate to compare to ISCCP IR only algortihm (i.e. 
                                                        ! you can compare to nighttime ISCCP data with this option) (1)
   integer :: isccp_topheight_direction = 2             ! direction for finding atmosphere pressure level with 
                                                        ! interpolated temperature equal to the radiance
                                                        ! determined cloud-top temperature
                                                        ! 1 = find the *lowest* altitude (highest pressure) level 
                                                        ! with interpolated temperature 
                                                        ! equal to the radiance determined cloud-top temperature
                                                        ! 2 = find the *highest* altitude (lowest pressure) level 
                                                        ! with interpolated temperature
                                                        ! equal to the radiance determined cloud-top temperature
                                                        ! ONLY APPLICABLE IF top_height EQUALS 1 or 3
                                                        ! 1 = default setting in COSP v1.1, matches all versions of 
                                                        ! ISCCP simulator with versions numbers 3.5.1 and lower
                                                        ! 2 = default setting in COSP v1.3. default since V4.0 of ISCCP simulator
   !per Alejandro: cosp rttov gaseous inputs are also mass mixing ratios
   !values from cosp_test.F90, units (kg/kg)
   !Mixing ratio C02 (5.241e-04), Mixing ratio CH4 (9.139e-07), 
   !Mixing ratio N20 (4.665e-07), Mixing ratio CO (2.098e-07)
   !I get CO2, CH4, N20 from cam radiation interface.

!! Other variables
    integer,parameter :: nhydro = 9                     ! number of COSP hydrometeor classes


   integer :: Npoints_it = 10000

   real(r8), parameter :: emsfc_lw = 0.99_r8            ! longwave emissivity of   surface at 10.5 microns 
                                                        ! set value same as in cloudsimulator.F90
   logical :: use_precipitation_fluxes = .true.

   logical :: use_reff = .true.                                 ! True if effective radius to be used by radar simulator 
                                                        ! (always used by lidar)

   
   integer, parameter :: Platform = 1                   ! satellite platform (1)
   integer, parameter :: Satellite = 15                 ! satellite (15)
   integer, parameter :: Instrument = 0                 ! instrument (0)
   integer, parameter :: Nchannels = 8                  ! Number of channels to be computed (8)
   integer, parameter :: Channels(Nchannels) =  (/1,3,5,6,8,10,11,13/)  
                                                        ! Channel numbers (match and supply Nchannels) 
                                                        ! (1,3,5,6,8,10,11,13,)
   real(r8), parameter :: Surfem(Nchannels) = (/0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8/)               
                                                        ! Surface emissivity (match and supply Nchannels)
                                                        ! (0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,)
   real(r8), parameter :: ZenAng = 50._r8               ! Satellite Zenith Angle (50)
   real(r8), parameter :: co =  2.098e-07_r8            ! Mixing ratio CO (2.098e-07), not used in cospv1.1
   real(r8), dimension(1,1) :: co2, ch4, n2o

!----------------------

   
   co2(1,1) = 1e-07_r8;
   ch4(1,1) = 1e-07_r8;
   n2o(1,1) = 1e-07_r8;
   



   call construct_cosp_gridbox(time, &                          ! 1 double precision = real(r8) X
                                time_bnds, &                    ! 1 double precision = real(r8)
                                radar_freq, &                   ! 2 real(r8) X
                                surface_radar, &                ! 3 integer X
                                use_mie_tables, &               ! 4 integer X
                                use_gas_abs, &                  ! 5 integer X
                                do_ray, &                       ! 6 integer X 
                                melt_lay, &                     ! 7 integer X 
                                k2, &                           ! 8 real(r8) X
#ifdef GPM_KU
                                gpmku_radar_freq, &                   ! 2 real(r8) X
                                gpmku_surface_radar, &                ! 3 integer X
                                gpmku_use_mie_tables, &               ! 4 integer X
                                gpmku_use_gas_abs, &                  ! 5 integer X
                                gpmku_do_ray, &                       ! 6 integer X 
                                gpmku_melt_lay, &                     ! 7 integer X 
                                gpmku_k2, &                           ! 8 real(r8) X
#endif
#ifdef GPM_KA
                                gpmka_radar_freq, &                   ! 2 real(r8) X
                                gpmka_surface_radar, &                ! 3 integer X
                                gpmka_use_mie_tables, &               ! 4 integer X
                                gpmka_use_gas_abs, &                  ! 5 integer X
                                gpmka_do_ray, &                       ! 6 integer X 
                                gpmka_melt_lay, &                     ! 7 integer X 
                                gpmka_k2, &                           ! 8 real(r8) X
#endif

                                Npoints, &                      ! 9 integer X, same as CAM's ncol
                                Nlevels, &                      ! 10 integer X
                                ncolumns,&                      ! 11 integer X
                                nhydro,&                        ! 12 integer X
                                Nprmts_max_hydro,&              ! 13 integer X
                                Naero,&                         ! 14 integer X
                                Nprmts_max_aero,&               ! 15 integer X
                                Npoints_it, &                   ! 16 integer X
                                lidar_ice_type,&                ! 17 integer X
                                isccp_topheight,&               ! 18 integer X
                                isccp_topheight_direction,&     ! 19 integer X
                                overlap,&                       ! 20 integer X
                                emsfc_lw, &                     ! 21 real X
                                use_precipitation_fluxes,&      ! 22 logical X
                                use_reff, &                     ! 23 logical X
                                Platform, &                     ! 24 integer X
                                Satellite, &                    ! 25 integer X
                                Instrument, &                   ! 26 integer X
                                Nchannels, &                    ! 27 integer X
                                ZenAng, &                       ! 28 real(r8) X
                                Channels(1:Nchannels),&         ! 29 integer X
                                Surfem(1:Nchannels),&           ! 30 real(r8) X
                                co2(1,1),&                      ! 31 real(r8) X
                                ch4(1,1),&                      ! 32 real(r8) X
                                n2o(1,1),&                      ! 33 real(r8) X
                                co,&                            ! 34 real(r8) X
                                gbx)                            ! OUT

  end subroutine construct_cosp_gridbox_lite


end program gmioffline

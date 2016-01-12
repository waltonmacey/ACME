#define GPM_KU
#define GPM_KA


program gpmtest
   use gpm_gmi_mod, only:gpm_gmi_addsensor, gpm_gmi_adddefaultsensor, gpm_gmi_init, gpm_gmi_run, gpm_gmi_clean, gpm_gmi_toabt
   use mod_cosp_types, only: cosp_gridbox, construct_cosp_gridbox, free_cosp_gridbox
   use shr_kind_mod, only: r8 =>shr_kind_r8
  ! use netcdf
   implicit none
   
   type(cosp_gridbox) :: gbx
   
  
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

   integer :: Npoints = 2
   integer :: Nlevels = 5
   integer :: ncolumns = 1

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

   character (len = *), parameter :: file_name =   "/project/projectdirs/acme/ylu/CRTM_Test/ERA_Interim/20150701-20150702.nc"
   integer, parameter :: nlat = 241
   integer, parameter :: nlon = 480
   integer, parameter :: nlev = 37
   integer, parameter :: ntime = 2

   integer q_in(nlon, nlat, nlev, ntime)
   integer :: ncid, varid
  
   type(gpm_gmi_toabt), allocatable :: toabt(:)

  ! call check( nf90_open(file_name, NF90_NOWRITE, ncid) )
  ! call check( nf90_inq_varid(ncid, "q", varid) )
  ! call check( nf90_get_var(ncid, varid, q_in) )

  ! call check ( nf90_close(ncid))





   
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

   ! assign some dummy values to gbx fields
   gbx%p(1,:) = (/900, 700, 500, 300, 100/)
   gbx%pint(1,:) = (/1000, 800, 600, 400, 200, 50/)
   gbx%T(1,:) = (/270, 250, 230, 210, 200/)
   gbx%q(1,:) = (/0.003, 0.003, 0.003, 0.003, 0.003/)

   gbx%p(2,:) = (/900, 700, 500, 300, 100/)
   gbx%pint(2,:) = (/1000, 800, 600, 400, 200, 50/)
   gbx%T(2,:) = (/280, 260, 240, 220, 200/)
   gbx%q(2,:) = (/0.005, 0.005, 0.005, 0.005, 0.005/)

   gbx%mr_hydro(1,:,1) = (/0.001, 0.000, 0.000, 0.000, 0.000/)
   gbx%mr_hydro(1,:,2) = (/0.000, 0.001, 0.000, 0.000, 0.000/)
   gbx%mr_hydro(1,:,3) = (/0.000, 0.000, 0.001, 0.000, 0.000/)
   gbx%mr_hydro(1,:,4) = (/0.000, 0.000, 0.000, 0.001, 0.000/)
   gbx%mr_hydro(1,:,5) = (/0.000, 0.000, 0.000, 0.000, 0.001/)
   gbx%mr_hydro(1,:,6) = (/0.002, 0.000, 0.000, 0.000, 0.000/)
   gbx%mr_hydro(1,:,7) = (/0.000, 0.002, 0.000, 0.000, 0.000/)
   gbx%mr_hydro(1,:,8) = (/0.000, 0.000, 0.002, 0.000, 0.000/)
   gbx%mr_hydro(1,:,9) = (/0.000, 0.000, 0.000, 0.002, 0.000/)

   gbx%mr_hydro(2,:,:) = gbx%mr_hydro(1,:,:)*100 

   gbx%Reff = 20 

   ! call gpm_gmi subroutines  
  print *, "debug   1" 
   call gpm_gmi_addsensor('gmi_gpm','GPM_GMI_TB', 52.8, 48.5)
   call gpm_gmi_addsensor('gmi_gpm','GPM_GMI_TB', 0.0, 0.0, channel_subset = (/1/))
!   call gpm_gmi_addsensor('atms_npp', 'ATMS_NPP_TB', 0.0, 0.0)
   call gpm_gmi_adddefaultsensor('gpm-gmi-lowfreq')
   call gpm_gmi_adddefaultsensor('gpm-gmi-highfreq')
   print *, "debug    2"
   call gpm_gmi_init(gbx)   
   call gpm_gmi_run(gbx, toabt)

   print *, "clean starts here"
   call gpm_gmi_clean()
print *, "clean gbx"
   call free_cosp_gridbox(gbx)

   
   print *, 'Run successful'


end program gpmtest


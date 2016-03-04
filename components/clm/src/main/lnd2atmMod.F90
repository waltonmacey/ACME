#include "../../../cam/src/physics/cosp/cosp_gpm_debugflag.F90"
module lnd2atmMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle lnd2atm mapping
  !
  ! !USES:
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use shr_infnan_mod       , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod          , only : errMsg => shr_log_errMsg
  use shr_megan_mod        , only : shr_megan_mechcomps_n
  use clm_varpar           , only : numrad, ndst, nlevgrnd !ndst = number of dust bins.
  use clm_varcon           , only : rair, grav, cpair, hfus, tfrz, spval
  use clm_varctl           , only : iulog, use_c13, use_cn, use_lch4, use_voc
  use tracer_varcon        , only : is_active_betr_bgc
  use seq_drydep_mod       , only : n_drydep, drydep_method, DD_XLND
  use decompMod            , only : bounds_type
  use subgridAveMod        , only : p2g, c2g 
  use lnd2atmType          , only : lnd2atm_type
  use atm2lndType          , only : atm2lnd_type
  use ch4Mod               , only : ch4_type
  use CNCarbonFluxType     , only : carbonflux_type
  use DUSTMod              , only : dust_type
  use DryDepVelocity       , only : drydepvel_type
  use VocEmissionMod       , only : vocemis_type
  use EnergyFluxType       , only : energyflux_type
  use FrictionVelocityType , only : frictionvel_type
  use SolarAbsorbedType    , only : solarabs_type
  use SurfaceAlbedoType    , only : surfalb_type
  use TemperatureType      , only : temperature_type
  use WaterFluxType        , only : waterflux_type
  use WaterstateType       , only : waterstate_type
  use GridcellType         , only : grc                
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: lnd2atm
  public :: lnd2atm_minimal
#ifdef GPM_GMI2
  public :: lnd2atm_GPMGMI
#endif
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine lnd2atm_minimal(bounds, &
      waterstate_vars, surfalb_vars, energyflux_vars, lnd2atm_vars)
    !
    ! !DESCRIPTION:
    ! Compute clm_l2a_vars component of gridcell derived type. This routine computes
    ! the bare minimum of components necessary to get the first step of a
    ! run started.
    !
    ! !USES:
    use clm_varcon, only : sb
    !
    ! !ARGUMENTS:
    type(bounds_type)     , intent(in)    :: bounds  
    type(waterstate_type) , intent(in)    :: waterstate_vars
    type(surfalb_type)    , intent(in)    :: surfalb_vars
    type(energyflux_type) , intent(in)    :: energyflux_vars
    type(lnd2atm_type)    , intent(inout) :: lnd2atm_vars 
    !
    ! !LOCAL VARIABLES:
    integer :: g                                    ! index
    real(r8), parameter :: amC   = 12.0_r8          ! Atomic mass number for Carbon
    real(r8), parameter :: amO   = 16.0_r8          ! Atomic mass number for Oxygen
    real(r8), parameter :: amCO2 = amC + 2.0_r8*amO ! Atomic mass number for CO2
    ! The following converts g of C to kg of CO2
    real(r8), parameter :: convertgC2kgCO2 = 1.0e-3_r8 * (amCO2/amC)
    !------------------------------------------------------------------------

    call c2g(bounds, &
         waterstate_vars%h2osno_col (bounds%begc:bounds%endc), &
         lnd2atm_vars%h2osno_grc    (bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    do g = bounds%begg,bounds%endg
       lnd2atm_vars%h2osno_grc(g) = lnd2atm_vars%h2osno_grc(g)/1000._r8
    end do

    call c2g(bounds, nlevgrnd, &
         waterstate_vars%h2osoi_vol_col (bounds%begc:bounds%endc, :), &
         lnd2atm_vars%h2osoi_vol_grc    (bounds%begg:bounds%endg, :), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    call p2g(bounds, numrad, &
         surfalb_vars%albd_patch (bounds%begp:bounds%endp, :), &
         lnd2atm_vars%albd_grc   (bounds%begg:bounds%endg, :), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    call p2g(bounds, numrad, &
         surfalb_vars%albi_patch (bounds%begp:bounds%endp, :), &
         lnd2atm_vars%albi_grc   (bounds%begg:bounds%endg, :), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    call p2g(bounds, &
         energyflux_vars%eflx_lwrad_out_patch (bounds%begp:bounds%endp), &
         lnd2atm_vars%eflx_lwrad_out_grc      (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    do g = bounds%begg,bounds%endg
       lnd2atm_vars%t_rad_grc(g) = sqrt(sqrt(lnd2atm_vars%eflx_lwrad_out_grc(g)/sb))
    end do

  end subroutine lnd2atm_minimal

  !------------------------------------------------------------------------
  subroutine lnd2atm(bounds, &
       atm2lnd_vars, surfalb_vars, temperature_vars, frictionvel_vars, &
       waterstate_vars, waterflux_vars, energyflux_vars, &
       solarabs_vars, carbonflux_vars, drydepvel_vars, &
       vocemis_vars, dust_vars, ch4_vars, lnd2atm_vars) 
    !
    ! !DESCRIPTION:
    ! Compute lnd2atm_vars component of gridcell derived type
    !
    ! !USES:
    use ch4varcon  , only : ch4offline
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)     :: bounds  
    type(atm2lnd_type)     , intent(in)     :: atm2lnd_vars
    type(surfalb_type)     , intent(in)     :: surfalb_vars
    type(temperature_type) , intent(in)     :: temperature_vars
    type(frictionvel_type) , intent(in)     :: frictionvel_vars
    type(waterstate_type)  , intent(inout)  :: waterstate_vars
    type(waterflux_type)   , intent(in)     :: waterflux_vars
    type(energyflux_type)  , intent(in)     :: energyflux_vars
    type(solarabs_type)    , intent(in)     :: solarabs_vars
    type(carbonflux_type)  , intent(in)     :: carbonflux_vars
    type(drydepvel_type)   , intent(in)     :: drydepvel_vars
    type(vocemis_type)     , intent(in)     :: vocemis_vars
    type(dust_type)        , intent(in)     :: dust_vars
    type(ch4_type)         , intent(in)     :: ch4_vars
    type(lnd2atm_type)     , intent(inout)  :: lnd2atm_vars 
    !
    ! !LOCAL VARIABLES:
    integer :: g             ! index
    real(r8), parameter :: amC   = 12.0_r8 ! Atomic mass number for Carbon
    real(r8), parameter :: amO   = 16.0_r8 ! Atomic mass number for Oxygen
    real(r8), parameter :: amCO2 = amC + 2.0_r8*amO ! Atomic mass number for CO2
    ! The following converts g of C to kg of CO2
    real(r8), parameter :: convertgC2kgCO2 = 1.0e-3_r8 * (amCO2/amC)
    !------------------------------------------------------------------------

    !----------------------------------------------------
    ! lnd -> atm
    !----------------------------------------------------
    
    ! First, compute the "minimal" set of fields.
    call lnd2atm_minimal(bounds, &
         waterstate_vars, surfalb_vars, energyflux_vars, lnd2atm_vars)

    call p2g(bounds, &
         temperature_vars%t_ref2m_patch (bounds%begp:bounds%endp), &
         lnd2atm_vars%t_ref2m_grc       (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         waterstate_vars%q_ref2m_patch (bounds%begp:bounds%endp), &
         lnd2atm_vars%q_ref2m_grc      (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         frictionvel_vars%u10_clm_patch (bounds%begp:bounds%endp), &
         lnd2atm_vars%u_ref10m_grc      (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         energyflux_vars%taux_patch (bounds%begp:bounds%endp), &
         lnd2atm_vars%taux_grc      (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         energyflux_vars%tauy_patch (bounds%begp:bounds%endp), &
         lnd2atm_vars%tauy_grc      (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         waterflux_vars%qflx_evap_tot_patch (bounds%begp:bounds%endp), &
         lnd2atm_vars%qflx_evap_tot_grc     (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    call p2g(bounds, &
         solarabs_vars%fsa_patch (bounds%begp:bounds%endp), &
         lnd2atm_vars%fsa_grc    (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    call p2g(bounds, &
         frictionvel_vars%fv_patch (bounds%begp:bounds%endp), &
         lnd2atm_vars%fv_grc       (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         frictionvel_vars%ram1_patch (bounds%begp:bounds%endp), &
         lnd2atm_vars%ram1_grc       (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g( bounds, &
         energyflux_vars%eflx_sh_tot_patch (bounds%begp:bounds%endp), &
         lnd2atm_vars%eflx_sh_tot_grc      (bounds%begg:bounds%endg), &
         p2c_scale_type='unity',c2l_scale_type='urbanf',l2g_scale_type='unity')
    do g = bounds%begg, bounds%endg
       lnd2atm_vars%eflx_sh_tot_grc(g) =  lnd2atm_vars%eflx_sh_tot_grc(g) - &
            energyflux_vars%eflx_dynbal_grc(g) 
    enddo

    call p2g(bounds, &
         energyflux_vars%eflx_lh_tot_patch (bounds%begp:bounds%endp), &
         lnd2atm_vars%eflx_lh_tot_grc      (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    if (use_cn) then
       call c2g(bounds, &
            carbonflux_vars%nee_col(bounds%begc:bounds%endc), &
            lnd2atm_vars%nee_grc   (bounds%begg:bounds%endg), &
            c2l_scale_type= 'unity', l2g_scale_type='unity')

       if (use_lch4) then
          if (.not. ch4offline) then
             ! Adjust flux of CO2 by the net conversion of mineralizing C to CH4
             do g = bounds%begg,bounds%endg
                ! nem is in g C/m2/sec
                lnd2atm_vars%nee_grc(g) = lnd2atm_vars%nee_grc(g) + lnd2atm_vars%nem_grc(g) 
             end do
          end if
       end if

       ! Convert from gC/m2/s to kgCO2/m2/s
       do g = bounds%begg,bounds%endg
          lnd2atm_vars%nee_grc(g) = lnd2atm_vars%nee_grc(g)*convertgC2kgCO2
       end do
    else
       do g = bounds%begg,bounds%endg
          lnd2atm_vars%nee_grc(g) = 0._r8
       end do
    end if

    ! drydepvel
    if ( n_drydep > 0 .and. drydep_method == DD_XLND ) then
       call p2g(bounds, n_drydep, &
            drydepvel_vars%velocity_patch (bounds%begp:bounds%endp, :), &
            lnd2atm_vars%ddvel_grc        (bounds%begg:bounds%endg, :), &
            p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
    endif

    ! voc emission flux
    if (use_voc .and. shr_megan_mechcomps_n>0) then
       call p2g(bounds, shr_megan_mechcomps_n, &
            vocemis_vars%vocflx_patch(bounds%begp:bounds%endp,:), &
            lnd2atm_vars%flxvoc_grc  (bounds%begg:bounds%endg,:), &
            p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
    end if

    ! dust emission flux
    call p2g(bounds, ndst, &
         dust_vars%flx_mss_vrt_dst_patch(bounds%begp:bounds%endp, :), &
         lnd2atm_vars%flxdst_grc        (bounds%begg:bounds%endg, :), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')


    ! ch4 flux
    if (use_lch4 .and. (.not. is_active_betr_bgc)) then
       call c2g( bounds,     &
            ch4_vars%ch4_surf_flux_tot_col (bounds%begc:bounds%endc), &
            lnd2atm_vars%flux_ch4_grc      (bounds%begg:bounds%endg), &
            c2l_scale_type= 'unity', l2g_scale_type='unity' )
    end if

    !----------------------------------------------------
    ! lnd -> rof
    !----------------------------------------------------

    call c2g( bounds, &
         waterflux_vars%qflx_runoff_col (bounds%begc:bounds%endc), &
         lnd2atm_vars%qflx_rofliq_grc   (bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )
    do g = bounds%begg, bounds%endg
       lnd2atm_vars%qflx_rofliq_grc(g) = lnd2atm_vars%qflx_rofliq_grc(g) - waterflux_vars%qflx_liq_dynbal_grc(g)
    enddo

    call c2g( bounds, &
         waterflux_vars%qflx_surf_col (bounds%begc:bounds%endc), &
         lnd2atm_vars%qflx_rofliq_qsur_grc   (bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

    call c2g( bounds, &
         waterflux_vars%qflx_drain_col (bounds%begc:bounds%endc), &
         lnd2atm_vars%qflx_rofliq_qsub_grc   (bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

    call c2g( bounds, &
         waterflux_vars%qflx_qrgwl_col (bounds%begc:bounds%endc), &
         lnd2atm_vars%qflx_rofliq_qgwl_grc   (bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

    call c2g( bounds, &
         waterflux_vars%qflx_snwcp_ice_col(bounds%begc:bounds%endc),  &
         lnd2atm_vars%qflx_rofice_grc     (bounds%begg:bounds%endg),  & 
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )
    do g = bounds%begg, bounds%endg
       lnd2atm_vars%qflx_rofice_grc(g) = lnd2atm_vars%qflx_rofice_grc(g) - waterflux_vars%qflx_ice_dynbal_grc(g)          
    enddo

    ! calculate total water storage for history files
    ! first set tws to gridcell total endwb
    ! second add river storage as gridcell average depth (1.e-3 converts [m3/km2] to [mm])
    ! TODO - this was in BalanceCheckMod - not sure where it belongs?

    call c2g( bounds, &
         waterstate_vars%endwb_col(bounds%begc:bounds%endc), &
         waterstate_vars%tws_grc  (bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )
    do g = bounds%begg, bounds%endg
       waterstate_vars%tws_grc(g) = waterstate_vars%tws_grc(g) + atm2lnd_vars%volr_grc(g) / grc%area(g) * 1.e-3_r8
    enddo

  end subroutine lnd2atm
#ifdef GPM_GMI2
  !------------------------------------------------------------------------
  subroutine lnd2atm_GPMGMI(bounds, &
       temperature_vars, &
       waterstate_vars,  &
       canopystate_vars, &
       soilstate_vars,   &
       lnd2atm_vars) 
    !
    ! !DESCRIPTION:
    ! Compute lnd2atm_vars component of gridcell derived type for GPM surface
    ! property calculations
    !
    ! !USES:
    use CanopystateType,  only : canopystate_type
    use SoilstateType,    only : soilstate_type
    use PatchType,        only : pft
    use ColumnType,       only : col
    use landunit_varcon,  only : istdlak
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)     :: bounds  
    type(temperature_type) , intent(in)     :: temperature_vars
    type(waterstate_type)  , intent(in)     :: waterstate_vars
    type(canopystate_type) , intent(in)     :: canopystate_vars
    type(soilstate_type)   , intent(in)     :: soilstate_vars
    type(lnd2atm_type)     , intent(inout)  :: lnd2atm_vars 
    !
    ! Parameters used to convert PFT types to vegetation density and water
    ! content. It is copied from CRTM NESDIS_LandEM_Module.f90 
! Specific Density
    REAL(r8), PARAMETER, dimension(0:13) :: veg_rho  = (/ 0.33_r8,     &
                          0.40_r8, 0.40_r8, 0.40_r8, 0.40_r8, 0.40_r8, &
                          0.25_r8, 0.25_r8, 0.40_r8, 0.40_r8, 0.40_r8, &
                          0.40_r8, 0.33_r8, 0.33_r8            /)
! MGE
    REAL(r8), PARAMETER, dimension(0:13) :: veg_mge  = (/ 0.50_r8,     &
                          0.45_r8, 0.45_r8, 0.45_r8, 0.40_r8, 0.40_r8, &
                          0.30_r8, 0.35_r8, 0.30_r8, 0.30_r8, 0.40_r8, &
                          0.30_r8, 0.50_r8, 0.40_r8            /)
    ! parameters used to convert CLM PFT types to GFS vegetation types
    integer, parameter, dimension(0:24) :: pft2gfs = &
  (/ & ! -------------------------------------------------------- 
    11 ,&  !   0  => not vegetated
    4  ,&  !   1  => needleleaf evergreen temperate tree
    4  ,&  !   2  => needleleaf evergreen boreal tree
    5  ,&  !   3  => needleleaf deciduous boreal tree
    1  ,&  !   4  => broadleaf evergreen tropical tree
    1  ,&  !   5  => broadleaf evergreen temperate tree
    2  ,&  !   6  => broadleaf deciduous tropical tree
    2  ,&  !   7  => broadleaf deciduous temperate tree
    2  ,&  !   8  => broadleaf deciduous boreal tree
    8  ,&  !   9  => broadleaf evergreen shrub
    9  ,&  !   10 => broadleaf deciduous temperate shrub
   10  ,&  !   11 => broadleaf deciduous boreal shrub
   10  ,&  !   12 => c3 arctic grass
    7  ,&  !   13 => c3 non-arctic grass
    6  ,&  !   14 => c4 grass
   12  ,&  !   15 => c3_crop
   12  ,&  !   16 => c3_irrigated
   12  ,&  !   17 => corn
   12  ,&  !   18 => irrigated corn
   12  ,&  !   19 => spring temperate cereal
   12  ,&  !   20 => irrigated spring temperate cereal
   12  ,&  !   21 => winter temperate cereal
   12  ,&  !   22 => irrigated winter temperate cereal
   12  ,&  !   23 => soybean
   12   &  !   24 => irrigated soybean
    /)
  ! -------------------------------------------------------- 
    ! !LOCAL VARIABLES:
    integer :: p, c             ! index
    real(r8), allocatable :: local_h2osoi_vol_col (:)
    real(r8), allocatable :: local_vegfrac_patch  (:)
    real(r8), allocatable :: local_soil_t_col     (:)
    real(r8), allocatable :: local_lai_patch      (:)
    real(r8), allocatable :: local_sandfrac_patch (:)
    real(r8), allocatable :: local_clayfrac_patch (:)
    real(r8), allocatable :: local_vegrho_patch   (:)
    real(r8), allocatable :: local_vegmge_patch   (:)
    real(r8), allocatable :: local_landfrac_col   (:)

    !------------------------------------------------------------------------
    ! allocate variables
    allocate(local_h2osoi_vol_col(bounds%begc:bounds%endc))  ;    local_h2osoi_vol_col(bounds%begc:bounds%endc) = spval
    allocate(local_soil_t_col    (bounds%begc:bounds%endc))  ;    local_soil_t_col    (bounds%begc:bounds%endc) = spval
    allocate(local_landfrac_col  (bounds%begc:bounds%endc))  ;    local_landfrac_col  (bounds%begc:bounds%endc) = spval
    
    allocate(local_vegfrac_patch (bounds%begp:bounds%endp))  ;    local_vegfrac_patch (bounds%begp:bounds%endp) = spval
    allocate(local_lai_patch     (bounds%begp:bounds%endp))  ;    local_lai_patch     (bounds%begp:bounds%endp) = spval
    allocate(local_sandfrac_patch(bounds%begp:bounds%endp))  ;    local_sandfrac_patch(bounds%begp:bounds%endp) = spval
    allocate(local_clayfrac_patch(bounds%begp:bounds%endp))  ;    local_clayfrac_patch(bounds%begp:bounds%endp) = spval
    allocate(local_vegrho_patch  (bounds%begp:bounds%endp))  ;    local_vegrho_patch  (bounds%begp:bounds%endp) = spval
    allocate(local_vegmge_patch  (bounds%begp:bounds%endp))  ;    local_vegmge_patch  (bounds%begp:bounds%endp) = spval

    ! assign values to local variables based on column types at patch level
    do p = bounds%begp, bounds%endp
      if  ( col%itype(pft%column(p) ) /= istdlak )  then! not deep lake
        if (pft%itype(p) /= 0)       then  ! not bare soil
          local_vegfrac_patch (p) = 1
          local_lai_patch     (p) = canopystate_vars%tlai_patch(p)
          local_sandfrac_patch(p) = soilstate_vars%sandfrac_patch(p)
          local_clayfrac_patch(p) = soilstate_vars%clayfrac_patch(p)
          local_vegrho_patch  (p) = veg_rho(pft2gfs(pft%itype(p)))
          local_vegmge_patch  (p) = veg_mge(pft2gfs(pft%itype(p)))
        end if 
      end if
    end do !patch

    ! assign values to local variables based on column types at column level
    do c = bounds%begc, bounds%endc
      if ( col%itype(c) /= istdlak ) then
        local_landfrac_col  (c) = 1
        local_h2osoi_vol_col(c) = waterstate_vars%h2osoi_vol_col(c,1) ! use the first layer water content
        local_soil_t_col    (c) = temperature_vars%t_soi10cm_col(c)
      end if
    end do ! column

    !----------------------------------------------------
    ! lnd -> atm
    !----------------------------------------------------
    ! column level variables
    call c2g(bounds, &
         local_h2osoi_vol_col            (bounds%begc:bounds%endc), &
         lnd2atm_vars%unisoilw_grc       (bounds%begg:bounds%endg), &
         c2l_scale_type= 'unity', l2g_scale_type='unity')

    call c2g(bounds, &
         local_soil_t_col                (bounds%begc:bounds%endc), &
         lnd2atm_vars%soiltemp_grc       (bounds%begg:bounds%endg), &
         c2l_scale_type= 'unity', l2g_scale_type='unity')

    call c2g(bounds, &
         local_landfrac_col              (bounds%begc:bounds%endc), &
         lnd2atm_vars%landfrac_grc       (bounds%begg:bounds%endg), &
         c2l_scale_type= 'unity', l2g_scale_type='unity')
    ! patch level variables
    call p2g(bounds, &
         local_vegfrac_patch             (bounds%begp:bounds%endp), &
         lnd2atm_vars%vegfrac_grc        (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         local_lai_patch                 (bounds%begp:bounds%endp), &
         lnd2atm_vars%lai_grc            (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         local_sandfrac_patch            (bounds%begp:bounds%endp), &
         lnd2atm_vars%sandfrac_grc       (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         local_clayfrac_patch            (bounds%begp:bounds%endp), &
         lnd2atm_vars%clayfrac_grc       (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         local_vegrho_patch              (bounds%begp:bounds%endp), &
         lnd2atm_vars%vegrho_grc         (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         local_vegmge_patch              (bounds%begp:bounds%endp), &
         lnd2atm_vars%vegmge_grc         (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
  end subroutine lnd2atm_GPMGMI

#endif  
end module lnd2atmMod

Module GPM_LandEM_Mod
! module use
  use WaterstateType,         only : waterstate_type
  use TemperatureType,        only : temperature_type
  use CanopystateType,        only : canopystate_type
  use SoilstateType,          only : soilstate_type
  use GPM_SurfaceOpticsType,  only : gpm_surfaceoptics_type

  use GridcellType      , only : grc                
  use LandunitType      , only : lun                
  use ColumnType        , only : col                
  use PatchType         , only : pft
  use decompMod         , only : bounds_type
  ! visibility
  private ! everything private by default

  public :: GPM_LandEM_calc



  contains

    subroutine GPM_LandEM_calc( bounds, &
                               waterstate_vars, &
                               temperature_vars,&
                               canopystate_vars,&
                               soilstate_vars,  &
                               gpm_surfaceoptics_vars)
! Calculate land surface microwave emissivity for GPM GMI
! This subroutine is an interface between CLM variables and CRTM functions that
! calculate land surface microwave emissivity
!
!______________________________________________________________________________________________________________________________
! Variable description                       |        CLM variables                     |         CRTM inputs (units)         |
!-----------------------------------------------------------------------------------------------------------------------------|
! volumetric vater content of the soil       |  waterstate_vars%h2osoi_vol_col [%]      | Soil_Moisture_Content    (g/cm^3)   |*
! vegetation fraction of the surface         |  1 or 0 depend on PFT type               | Vegetation_Fraction      (%, 0 to 1)|
! soil temperature                           | temperature_vars%t_soi10cm_col  [K]      | Soil_Temperature         (Kelvin)   |**
! land surface temperature                   | temperature_vars%t_grnd_col     [K]      | t_skin                   (Kelvin)   |
! leaf area index                            | canopystate_vars%tlai_patch     [?]      | LAI                      (m^2/m^2)  |***
! snow depth                                 |  waterstate_vars%snow_depth_col [m]      | Snow_Depth                          |
! fraction of sand                           |   soilstate_vars%sandfrac_patch [%]      | frac_sand                (%, 0 to 1)|
! fraction of clay                           |   soilstate_vars%sandfrac_patch [%]      | frac_clay                (%, 0 to 1)|
! bulk volume density of the soil (1.18-1.12)|   soilstate_vars%bd_col         [kg/m^3] | rhob_soil                (g/cm^3)   |****
!____________________________________________|__________________________________________|_____________________________________|
! Notes:
! *    . Soil moisture constant in CRTM seems have some discrepancy
! **   . There are several soil temperatures in CLM. Use this value for now.
! ***  . The current LAI is "patch canopy one-sided leaf area index, no burying
!        by snow". There is also an corresponding value that is "with burying by
!        snow". Use this value for now.
! **** . This value is for dry soil material. Not clear if the CRTM wants the
!        dry material.

      ! Uses:
      use landunit_varcon, only:  istsoil    ,& ! = 1  soil         landunit type (natural vegetation)
                                  istcrop    ,& ! = 2  crop         landunit type
                                  istice     ,& ! = 3  land ice     landunit type (glacier)
                                  istice_mec ,& ! = 4  land ice (multiple elevation classes) landunit type
                                  istdlak    ,& ! = 5  deep lake    landunit type (now used for all lakes)
                                  istwet     ,& ! = 6  wetland      landunit type (swamp, marsh, etc.)     
                                  isturb_MIN ,& ! = 7  minimum urban type index
                                  isturb_MAX    ! = 9  maximum urban type index
      
      ! Arguments:
      type(bounds_type),            intent(in)    :: bounds    ! bounds    
      type(waterstate_type),        intent(in)    :: waterstate_vars
      type(temperature_type),       intent(in)    :: temperature_vars
      type(canopystate_type),       intent(in)    :: canopystate_vars
      type(soilstate_type),         intent(in)    :: soilstate_vars
      type(gpm_surfaceoptics_type), intent(inout) :: gpm_surfaceoptics_vars
      
      ! local variables
      integer :: p, c, l, g  ! indices
      integer :: begc, endc, begp, endp

  !-----------------------------------------------------------------------
      begc = bounds%begc
      endc = bounds%endc
      begp = bounds%begp
      endp = bounds%endp

      ! associate variables
      associate(&
        soil_moist  =>   waterstate_vars%h2osoi_vol_col  , &
        soil_temp   =>  temperature_vars%t_soi10cm_col   , &
        t_skin      =>  temperature_vars%t_grnd_col      , &
        snow_depth  =>   waterstate_vars%snow_depth_col  , &
        lai         =>  canopystate_vars%tlai_patch      , &
        frac_sand   =>    soilstate_vars%sandfrac_patch  , &
        frac_clay   =>    soilstate_vars%clayfrac_patch  , &
        MWem        =>  gpm_surfaceoptics_vars%emissivityMW_patch &
        )

        ! loop through pfts and calculate surface MW emissivities at patch level
        do p = begp, endp
          c = pft%column(p)
          l = pft%landunit(p)
          g = pft%gridcell(p)

          ! ice
          if (lun%itype(l) == istice  .OR. &
              lun%itype(l) == istice_mec )

          ! water
          elseif (lun%itype(l) == istdlak )

          ! land
          elseif (  lun%itype(l) == istsoil  .OR. &
                    lun%itype(l) == istcrop  .OR. &
                    lun%itype(l) == istwet   .OR. &
                   (lun%itype(l) >= isturb_MIN .AND. lun%itype(l) <= isturb_MAX) )
                   
            
          ! unknown, error
          else


          endif

        end do

      end associate

    end subroutine GPM_LANDEM_Calc

end module GPM_LandEM_Mod

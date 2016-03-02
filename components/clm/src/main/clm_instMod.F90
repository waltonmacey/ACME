module clm_instMod

  !
  !-----------------------------------------
  ! Definition of component types
  !-----------------------------------------
  use AerosolType            , only : aerosol_type
  use CanopyStateType        , only : canopystate_type
  use ch4Mod                 , only : ch4_type
  use CNCarbonFluxType       , only : carbonflux_type
  use CNCarbonStateType      , only : carbonstate_type
  use CNDVType               , only : dgvs_type
  use CNStateType            , only : cnstate_type
  use CNNitrogenFluxType     , only : nitrogenflux_type
  use CNNitrogenStateType    , only : nitrogenstate_type
  use PhosphorusFluxType     , only : phosphorusflux_type
  use PhosphorusStateType    , only : phosphorusstate_type
  use CropType               , only : crop_type
  use DryDepVelocity         , only : drydepvel_type
  use DUSTMod                , only : dust_type
  use EnergyFluxType         , only : energyflux_type
  use FrictionVelocityType   , only : frictionvel_type
  use LakeStateType          , only : lakestate_type
  use PhotosynthesisType     , only : photosyns_type
  use SoilHydrologyType      , only : soilhydrology_type
  use SoilStateType          , only : soilstate_type
  use SolarAbsorbedType      , only : solarabs_type
  use SurfaceRadiationMod    , only : surfrad_type
  use SurfaceAlbedoMod       , only : SurfaceAlbedoInitTimeConst !TODO - can this be merged into the type?
  use SurfaceAlbedoType      , only : surfalb_type
  use TemperatureType        , only : temperature_type
  use WaterfluxType          , only : waterflux_type
  use WaterstateType         , only : waterstate_type
  use UrbanParamsType        , only : urbanparams_type
  use VOCEmissionMod         , only : vocemis_type
  use atm2lndType            , only : atm2lnd_type
  use lnd2atmType            , only : lnd2atm_type
  use lnd2glcMod             , only : lnd2glc_type
  use glc2lndMod             , only : glc2lnd_type
  use glcDiagnosticsMod      , only : glc_diagnostics_type
  use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
  use UrbanParamsType        , only : urbanparams_type   ! Constants
  use CNDecompCascadeConType , only : decomp_cascade_con ! Constants
  use CNDVType               , only : dgv_ecophyscon     ! Constants
  use EcophysConType         , only : ecophyscon         ! Constants
  use SoilorderConType       , only : soilordercon         ! Constants
  use EDEcophysConType       , only : EDecophyscon       ! ED Constants
  use EDBioType              , only : EDbio_type         ! ED type used to interact with CLM variables
  use ChemStateType          , only : chemstate_type     ! structure for chemical indices of the soil, such as pH and Eh

implicit none

  save
  public
  !
  !-----------------------------------------
  ! Instances of component types
  !-----------------------------------------
  !
  type(ch4_type)              :: ch4_vars
  type(carbonstate_type)      :: carbonstate_vars
  type(carbonstate_type)      :: c13_carbonstate_vars
  type(carbonstate_type)      :: c14_carbonstate_vars
  type(carbonflux_type)       :: carbonflux_vars
  type(carbonflux_type)       :: c13_carbonflux_vars
  type(carbonflux_type)       :: c14_carbonflux_vars
  type(nitrogenstate_type)    :: nitrogenstate_vars
  type(nitrogenflux_type)     :: nitrogenflux_vars
  type(dgvs_type)             :: dgvs_vars
  type(crop_type)             :: crop_vars
  type(cnstate_type)          :: cnstate_vars
  type(dust_type)             :: dust_vars
  type(vocemis_type)          :: vocemis_vars
  type(drydepvel_type)        :: drydepvel_vars
  type(aerosol_type)          :: aerosol_vars
  type(canopystate_type)      :: canopystate_vars
  type(energyflux_type)       :: energyflux_vars
  type(frictionvel_type)      :: frictionvel_vars
  type(lakestate_type)        :: lakestate_vars
  type(photosyns_type)        :: photosyns_vars
  type(soilstate_type)        :: soilstate_vars
  type(soilhydrology_type)    :: soilhydrology_vars
  type(solarabs_type)         :: solarabs_vars
  type(surfalb_type)          :: surfalb_vars
  type(surfrad_type)          :: surfrad_vars
  type(temperature_type)      :: temperature_vars
  type(urbanparams_type)      :: urbanparams_vars
  type(waterflux_type)        :: waterflux_vars
  type(waterstate_type)       :: waterstate_vars
  type(atm2lnd_type)          :: atm2lnd_vars
  type(glc2lnd_type)          :: glc2lnd_vars
  type(lnd2atm_type)          :: lnd2atm_vars
  type(lnd2glc_type)          :: lnd2glc_vars
  type(glc_diagnostics_type)  :: glc_diagnostics_vars
  class(soil_water_retention_curve_type), allocatable :: soil_water_retention_curve
  type(EDbio_type)            :: EDbio_vars
  type(phosphorusstate_type)  :: phosphorusstate_vars
  type(phosphorusflux_type)   :: phosphorusflux_vars
  type(chemstate_type)        :: chemstate_vars

end module clm_instMod

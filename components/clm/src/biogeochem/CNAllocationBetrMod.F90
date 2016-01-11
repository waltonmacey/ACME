module CNAllocationBetrMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module holding routines used in allocation model for coupled carbon
  ! nitrogen code.
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use clm_varcon          , only : dzsoi_decomp
  use clm_varctl          , only : use_c13, use_c14, use_nitrif_denitrif
  use abortutils          , only : endrun
  use decompMod           , only : bounds_type
  use subgridAveMod       , only : p2c
  use CanopyStateType     , only : canopystate_type
  use CNCarbonFluxType    , only : carbonflux_type
  use CNCarbonStateType   , only : carbonstate_type
  use CNNitrogenFluxType  , only : nitrogenflux_type
  use CNNitrogenStateType , only : nitrogenstate_type
  use CNStateType         , only : cnstate_type
  use PhotosynthesisType  , only : photosyns_type
  use CropType            , only : crop_type
  use EcophysConType      , only : ecophyscon
  use LandunitType        , only : lun
  use ColumnType          , only : col
  use PatchType           , only : pft
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: readCNAllocBetrParams
  public :: CNAllocationBetrInit         ! Initialization
  public :: calc_plant_nutrient_demand
  public :: plantCNAlloc
  type :: CNAllocParamsType
     real(r8) :: bdnr              ! bulk denitrification rate (1/s)
     real(r8) :: dayscrecover      ! number of days to recover negative cpool
     real(r8) :: compet_plant_no3  ! (unitless) relative compettiveness of plants for NO3
     real(r8) :: compet_plant_nh4  ! (unitless) relative compettiveness of plants for NH4
     real(r8) :: compet_decomp_no3 ! (unitless) relative competitiveness of immobilizers for NO3
     real(r8) :: compet_decomp_nh4 ! (unitless) relative competitiveness of immobilizers for NH4
     real(r8) :: compet_denit      ! (unitless) relative competitiveness of denitrifiers for NO3
     real(r8) :: compet_nit        ! (unitless) relative competitiveness of nitrifiers for NH4
  end type CNAllocParamsType
  !
  ! CNAllocParamsInst is populated in readCNAllocParams which is called in
  type(CNAllocParamsType),protected ::  CNAllocParamsInst
  !
  ! !PUBLIC DATA MEMBERS:
  character(len=*), parameter, public :: suplnAll='ALL'       ! Supplemental Nitrogen for all PFT's
  character(len=*), parameter, public :: suplnNon='NONE'      ! No supplemental Nitrogen
  character(len=15)          , public :: suplnitro = suplnNon ! Supplemental Nitrogen mode
  !
  ! !PRIVATE DATA MEMBERS:
  real(r8)              :: dt                   !decomp timestep (seconds)
  real(r8)              :: bdnr                 !bulk denitrification rate (1/s)
  real(r8)              :: dayscrecover         !number of days to recover negative cpool
  real(r8), allocatable :: arepr(:)             !reproduction allocation coefficient
  real(r8), allocatable :: aroot(:)             !root allocation coefficient
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine readCNAllocBetrParams ( ncid )
    !
    ! !USES:
    use ncdio_pio , only : file_desc_t,ncd_io

    ! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'CNAllocParamsType'                  !
    character(len=100) :: errCode = '-Error reading in parameters file:' !
    logical            :: readv                                          ! has variable been read in or not
    real(r8)           :: tempr                                          ! temporary to read in parameter
    character(len=100) :: tString                                        ! temp. var for reading
    !-----------------------------------------------------------------------

    ! read in parameters

    tString='bdnr'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNAllocParamsInst%bdnr=tempr

    tString='dayscrecover'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNAllocParamsInst%dayscrecover=tempr

    tString='compet_plant_no3'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNAllocParamsInst%compet_plant_no3=tempr

    tString='compet_plant_nh4'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNAllocParamsInst%compet_plant_nh4=tempr

    tString='compet_decomp_no3'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNAllocParamsInst%compet_decomp_no3=tempr

    tString='compet_decomp_nh4'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNAllocParamsInst%compet_decomp_nh4=tempr

    tString='compet_denit'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNAllocParamsInst%compet_denit=tempr

    tString='compet_nit'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNAllocParamsInst%compet_nit=tempr

  end subroutine readCNAllocBetrParams

  !-----------------------------------------------------------------------
  subroutine CNAllocationBetrInit ( bounds)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use clm_varcon      , only: secspday
    use clm_time_manager, only: get_step_size
    use clm_varpar      , only: crop_prog
    use clm_varctl      , only: iulog, cnallocate_carbon_only_set
    use shr_infnan_mod  , only: nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    character(len=32) :: subname = 'CNAllocationInit'
    logical           :: carbon_only
    !-----------------------------------------------------------------------

    if ( crop_prog )then
       allocate(arepr(bounds%begp:bounds%endp)); arepr(bounds%begp : bounds%endp) = nan
       allocate(aroot(bounds%begp:bounds%endp)); aroot(bounds%begp : bounds%endp) = nan
    end if


    ! set time steps
    dt = real( get_step_size(), r8 )

    ! set space-and-time parameters from parameter file
    bdnr         = CNAllocParamsInst%bdnr * (dt/secspday)
    dayscrecover = CNAllocParamsInst%dayscrecover

    ! Change namelist settings into private logical variables
    select case(suplnitro)
    case(suplnNon)
       Carbon_only = .false.
    case(suplnAll)
       Carbon_only = .true.
    case default
       write(iulog,*) 'Supplemental Nitrogen flag (suplnitro) can only be: ', &
            suplnNon, ' or ', suplnAll
       call endrun(msg='ERROR: supplemental Nitrogen flag is not correct'//&
            errMsg(__FILE__, __LINE__))
    end select

  end subroutine CNAllocationBetrInit


!-----------------------------------------------------------------------

  subroutine calc_plant_allometry_force


  do fp=1,num_soilp
    p = filter_soilp(fp)
    c = pft%column(p)


        ! 'ECA' or 'MIC' mode
        ! dynamic allocation based on light limitation (more woody growth) vs nutrient limitations (more fine root growth)
        ! set allocation coefficients
        if (cnallocate_carbon_only()) then
            cn_scalar(p) = 1.0_r8
            cp_scalar(p) = 1.0_r8
        else if (cnallocate_carbonnitrogen_only()) then
            cp_scalar(p) = 1.0_r8
        else if (cnallocate_carbonphosphorus_only()) then
            cn_scalar(p) = 1.0_r8
        end if

        call dynamic_plant_alloc(min(cn_scalar(p),cp_scalar(p)), &
             laisun(p)+laisha(p), allocation_leaf(p), allocation_stem(p), allocation_froot(p))


        f1(p) = allocation_froot(p) / allocation_leaf(p)
        f2(p) = croot_stem(ivt(p))
        f3(p) = allocation_stem(p) / allocation_leaf(p)

        ! modified wood allocation to be 2.2 at npp=800 gC/m2/yr, 0.2 at npp=0,
        ! constrained so that it does not go lower than 0.2 (under negative annsum_npp)
        ! There was an error in this formula in previous version, where the coefficient
        ! was 0.004 instead of 0.0025.
        ! This variable allocation is only for trees. Shrubs have a constant
        ! allocation as specified in the pft-physiology file.  The value is also used
        ! as a trigger here: -1.0 means to use the dynamic allocation (trees).

        f4 = flivewd(ivt(p))
        g1 = grperc(ivt(p))
        g2 = grpnow(ivt(p))

        cnl = leafcn(ivt(p))
        cnfr = frootcn(ivt(p))
        cnlw = livewdcn(ivt(p))
        cndw = deadwdcn(ivt(p))

        cpl =  leafcp(ivt(p))
        cpfr = frootcp(ivt(p))
        cplw = livewdcp(ivt(p))
        cpdw = deadwdcp(ivt(p))

        if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            if (croplive(p)) then
                f1(p) = aroot(p) / aleaf(p)
                f3(p) = astem(p) / aleaf(p)
                f5 = arepr(p) / aleaf(p)
                g1 = 0.25_r8
            else
                f1(p) = 0._r8
                f3(p) = 0._r8
                f5 = 0._r8
                g1 = 0.25_r8
            end if
        end if

        sminn_to_npool(p) = sminn_to_plant_patch(p)
        sminp_to_ppool(p) = sminp_to_plant_patch(p)

        plant_nalloc(p) = sminn_to_npool(p) + avail_retransn(p)
        plant_palloc(p) = sminp_to_ppool(p) + avail_retransp(p)
        plant_calloc(p) = availc(p)

        ! here no down-regulation on allocatable C here, NP limitation is implemented in leaf-level NP control on GPP
        if (woody(ivt(p)) == 1.0_r8) then
            c_allometry(p) = (1._r8+g1)*(1._r8+f1(p)+f3(p)*(1._r8+f2(p)))
            n_allometry(p) = 1._r8/cnl + f1(p)/cnfr + (f3(p)*f4*(1._r8+f2(p)))/cnlw + &
                (f3(p)*(1._r8-f4)*(1._r8+f2(p)))/cndw
            p_allometry(p) = 1._r8/cpl + f1(p)/cpfr + (f3(p)*f4*(1._r8+f2))/cplw + &
                (f3(p)*(1._r8-f4)*(1._r8+f2(p)))/cpdw

        else if (ivt(p) >= npcropmin) then ! skip generic crops
            cng = graincn(ivt(p))
            c_allometry(p) = (1._r8+g1)*(1._r8+f1(p)+f5+f3(p)*(1._r8+f2))
            n_allometry(p) = 1._r8/cnl + f1(p)/cnfr + f5/cng + (f3(p)*f4*(1._r8+f2(p)))/cnlw + &
                (f3(p)*(1._r8-f4)*(1._r8+f2(p)))/cndw
            p_allometry(p) = 1._r8/cpl + f1(p)/cpfr + (f3(p)*f4*(1._r8+f2(p)))/cplw + &
                (f3(p)*(1._r8-f4)*(1._r8+f2(p)))/cpdw

        else
            c_allometry(p) = 1._r8+g1+f1(p)+f1(p)*g1
            n_allometry(p) = 1._r8/cnl + f1(p)/cnfr
            p_allometry(p) = 1._r8/cpl + f1(p)/cpfr
        end if

        nlc_c(p) = plant_calloc(p)/c_allometry(p)
        nlc_n(p) = plant_nalloc(p)/n_allometry(p)
        nlc_p(p) = plant_palloc(p)/p_allometry(p)
    end if
  enddo
  end subroutine calc_plant_allometry_force
end module CNAllocationBetrMod

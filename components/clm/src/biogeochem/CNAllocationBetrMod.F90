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
  public :: CNAllocationBetrInit         ! Initialization
  public :: calc_plant_allometry_force

  ! CNAllocParamsInst is populated in readCNAllocParams which is called in
  type(CNAllocParamsType),protected ::  CNAllocParamsInst
  !
  ! !PUBLIC DATA MEMBERS:
  character(len=*), parameter, public :: suplnAll='ALL'       ! Supplemental Nitrogen for all PFT's
  character(len=*), parameter, public :: suplnNon='NONE'      ! No supplemental Nitrogen
  character(len=15)          , public :: suplnitro = suplnNon ! Supplemental Nitrogen mode
  !


  !-----------------------------------------------------------------------

contains



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


  end subroutine CNAllocationBetrInit


!-----------------------------------------------------------------------

  subroutine calc_plant_allometry_force(bounds, num_soilp, filter_soilp, canopystate_vars,carbonflux_vars, cnstate_vars)

  !
  ! !DESCRIPTION
  ! calculate the plant allometry driving force
  !
  ! USES
  use clm_varctl       , only: iulog,cnallocate_carbon_only,cnallocate_carbonnitrogen_only,&
                               cnallocate_carbonphosphorus_only
  use pftvarcon        , only: grperc, grpnow
  !ARGUMENTS
  implicit none
  type(bounds_type)        , intent(in) :: bounds
  integer                  , intent(in)    :: num_soilp        ! number of soil patches in filter
  integer                  , intent(in)    :: filter_soilp(:)  ! filter for soil patches
  type(cnstate_type)       , intent(inout) :: cnstate_vars
  type(canopystate_type)   , intent(in)    :: canopystate_vars
  type(carbonflux_type)    , intent(inout) :: carbonflux_vars
  real(r8) :: allocation_leaf(bounds%begp : bounds%endp)              ! fraction of NPP allocated into leaf
  real(r8) :: allocation_stem(bounds%begp : bounds%endp)              ! fraction of NPP allocated into stem
  real(r8) :: allocation_froot(bounds%begp : bounds%endp)              ! fraction of NPP allocated into froot
  real(r8) :: g1, g2
  real(r8) :: f4, f5
  real(r8):: cnl,cnfr,cnlw,cndw
  
 associate(                                                                 &
   ivt                          => pft%itype                              , & ! Input:  [integer  (:) ]  pft vegetation type

   f1                           => cnstate_vars%f1_patch                  , &
   f2                           => cnstate_vars%f2_patch                  , &
   f3                           => cnstate_vars%f3_patch                  , &
   c_allometry                  => cnstate_vars%c_allometry_patch         , & ! Output: [real(r8) (:)   ]  C allocation index (DIM)
   n_allometry                  => cnstate_vars%n_allometry_patch         , & ! Output: [real(r8) (:)   ]  N allocation index (DIM)
   p_allometry                  => cnstate_vars%p_allometry_patch         , & ! Output: [real(r8) (:)   ]  N allocation index (DIM)
   nlc_c                        => cnstate_vars%nlc_c_patch               , & !
   nlc_n                        => cnstate_vars%nlc_n_patch               , & !
   nlc_p                        => cnstate_vars%nlc_p_patch               , & !
   cn_scalar                    => cnstate_vars%cn_scalar                 , &
   cp_scalar                    => cnstate_vars%cp_scalar                 , &
   aroot                        => cnstate_vars%aroot_patch               , &
   arepr                        => cnstate_vars%arepr                     , &

   croot_stem                   => ecophyscon%croot_stem                  , & ! Input:  [real(r8) (:)   ]  allocation parameter: new coarse root C per new stem C (gC/gC)
   flivewd                      => ecophyscon%flivewd                     , & ! Input:  [real(r8) (:)   ]  allocation parameter: fraction of new wood that is live (phloem and ray parenchyma) (no units)
   leafcn                       => ecophyscon%leafcn                      , & ! Input:  [real(r8) (:)   ]  leaf C:N (gC/gN)
   frootcn                      => ecophyscon%frootcn                     , & ! Input:  [real(r8) (:)   ]  fine root C:N (gC/gN)
   livewdcn                     => ecophyscon%livewdcn                    , & ! Input:  [real(r8) (:)   ]  live wood (phloem and ray parenchyma) C:N (gC/gN)
   deadwdcn                     => ecophyscon%deadwdcn                    , & ! Input:  [real(r8) (:)   ]  dead wood (xylem and heartwood) C:N (gC/gN)
   livewdcp                     => ecophyscon%livewdcp                                   , & ! Input:  [real(r8) (:)   ]  live wood (phloem and ray parenchyma) C:P (gC/gP)
   deadwdcp                     => ecophyscon%deadwdcp                                   , & ! Input:  [real(r8) (:)   ]  dead wood (xylem and heartwood) C:P (gC/gP)

   leafcp                       => ecophyscon%leafcp                      , & ! Input:  [real(r8) (:)   ]  leaf C:P (gC/gP)
   frootcp                      => ecophyscon%frootcp                     , & ! Input:  [real(r8) (:)   ]  fine root C:P (gC/gP)

   laisun                       => canopystate_vars%laisun_patch          , & ! Input:  [real(r8) (:)   ]  sunlit projected leaf area index
   laisha                       => canopystate_vars%laisha_patch          , & ! Input:  [real(r8) (:)   ]  shaded projected leaf area index
   availc                       => carbonflux_vars%availc_patch           , & ! Output: [real(r8) (:)   ]  C flux available for allocation (gC/m2/s)

 )

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

        plant_nalloc(p) = sminn_to_plant_patch(p) + avail_retransn(p)
        plant_palloc(p) = sminp_to_plant_patch(p) + avail_retransp(p)
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

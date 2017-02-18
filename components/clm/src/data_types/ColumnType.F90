module ColumnType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  !DW  Converted from ColumnType
  !DW  Change the old function into PhysicalPropertiesType
  ! Column data type allocation and initialization
  ! -------------------------------------------------------- 
  ! column types can have values of
  ! -------------------------------------------------------- 
  !   1  => (istsoil)          soil (vegetated or bare soil)
  !   2  => (istcrop)          crop (only for crop configuration)
  !   3  => (istice)           land ice
  !   4  => (istice_mec)       land ice (multiple elevation classes)   
  !   5  => (istdlak)          deep lake
  !   6  => (istwet)           wetland
  !   71 => (icol_roof)        urban roof
  !   72 => (icol_sunwall)     urban sunwall
  !   73 => (icol_shadewall)   urban shadewall
  !   74 => (icol_road_imperv) urban impervious road
  !   75 => (icol_road_perv)   urban pervious road
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use clm_varpar     , only : nlevsno, nlevgrnd, nlevlak
  use clm_varcon     , only : spval, ispval
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  type, public :: column_physical_properties_type
     ! g/l/c/p hierarchy, local g/l/c/p cells only
     integer , pointer :: landunit             (:)   ! index into landunit level quantities
     real(r8), pointer :: wtlunit              (:)   ! weight (relative to landunit)
     integer , pointer :: gridcell             (:)   ! index into gridcell level quantities
     real(r8), pointer :: wtgcell              (:)   ! weight (relative to gridcell)
     integer , pointer :: pfti                 (:)   ! beginning pft index for each column
     integer , pointer :: pftf                 (:)   ! ending pft index for each column
     integer , pointer :: npfts                (:)   ! number of patches for each column

     ! topological mapping functionality
     integer , pointer :: itype                (:)   ! column type
     logical , pointer :: active               (:)   ! true=>do computations on this column 

     ! topography
     real(r8), pointer :: glc_topo             (:)   ! surface elevation (m)
     real(r8), pointer :: micro_sigma          (:)   ! microtopography pdf sigma (m)
     real(r8), pointer :: n_melt               (:)   ! SCA shape parameter
     real(r8), pointer :: topo_slope           (:)   ! gridcell topographic slope
     real(r8), pointer :: topo_std             (:)   ! gridcell elevation standard deviation

     ! vertical levels
     integer , pointer :: snl                  (:)   ! number of snow layers
     real(r8), pointer :: dz                   (:,:) ! layer thickness (m)  (-nlevsno+1:nlevgrnd) 
     real(r8), pointer :: z                    (:,:) ! layer depth (m) (-nlevsno+1:nlevgrnd) 
     real(r8), pointer :: zi                   (:,:) ! interface level below a "z" level (m) (-nlevsno+0:nlevgrnd) 
     real(r8), pointer :: zii                  (:)   ! convective boundary height [m]
     real(r8), pointer :: dz_lake              (:,:) ! lake layer thickness (m)  (1:nlevlak)
     real(r8), pointer :: z_lake               (:,:) ! layer depth for lake (m)
     real(r8), pointer :: lakedepth            (:)   ! variable lake depth (m)                             

     !DW variables from SoilorderconType.F90
     real(r8), allocatable :: smax(:)
     real(r8), allocatable :: ks_sorption(:)
     real(r8), allocatable :: r_weather(:)
     real(r8), allocatable :: r_adsorp(:)
     real(r8), allocatable :: r_desorp(:)
     real(r8), allocatable :: r_occlude(:)
     real(r8), allocatable :: k_s1_biochem(:)
     real(r8), allocatable :: k_s2_biochem(:)
     real(r8), allocatable :: k_s3_biochem(:)
     real(r8), allocatable :: k_s4_biochem(:)

   contains

     procedure, public :: Init => col_pp_init
     procedure, public :: Clean => col_pp_clean

  end type column_physical_properties_type

  type(column_physical_properties_type), public, target :: col_pp !column data structure (soil/snow/canopy columns)
  !------------------------------------------------------------------------

contains
  
  !------------------------------------------------------------------------
  subroutine col_pp_init(this, begc, endc)
    !
    ! !ARGUMENTS:
    !class(column_physical_properties_type)  :: this
   
    !DW  Follow 5 lines from SoilorderConType.F90 !USES:
    use clm_varpar, only : nsoilorder
    use soilorder_varcon, only : smax
    use soilorder_varcon, only : ks_sorption
    use soilorder_varcon, only : r_weather,r_adsorp,r_desorp,r_occlude
    use soilorder_varcon, only : k_s1_biochem,k_s2_biochem,k_s3_biochem,k_s4_biochem


    ! !ARGUMENTS:
    class(column_physical_properties_type)  :: this

    integer, intent(in) :: begc,endc

    !  !LOCAL VARIABLES
    integer :: m
    !------------------------------------------------------------------------

    ! The following is set in initGridCellsMod
    allocate(this%gridcell    (begc:endc))                     ; this%gridcell    (:)   = ispval
    allocate(this%wtgcell     (begc:endc))                     ; this%wtgcell     (:)   = nan
    allocate(this%landunit    (begc:endc))                     ; this%landunit    (:)   = ispval
    allocate(this%wtlunit     (begc:endc))                     ; this%wtlunit     (:)   = nan
    allocate(this%pfti        (begc:endc))                     ; this%pfti        (:)   = ispval
    allocate(this%pftf        (begc:endc))                     ; this%pftf        (:)   = ispval
    allocate(this%npfts       (begc:endc))                     ; this%npfts       (:)   = ispval
    allocate(this%itype       (begc:endc))                     ; this%itype       (:)   = ispval
    allocate(this%active      (begc:endc))                     ; this%active      (:)   = .false.

    ! The following is set in initVerticalMod
    allocate(this%snl         (begc:endc))                     ; this%snl         (:)   = ispval  !* cannot be averaged up
    allocate(this%dz          (begc:endc,-nlevsno+1:nlevgrnd)) ; this%dz          (:,:) = nan
    allocate(this%z           (begc:endc,-nlevsno+1:nlevgrnd)) ; this%z           (:,:) = nan
    allocate(this%zi          (begc:endc,-nlevsno+0:nlevgrnd)) ; this%zi          (:,:) = nan
    allocate(this%zii         (begc:endc))                     ; this%zii         (:)   = nan
    allocate(this%lakedepth   (begc:endc))                     ; this%lakedepth   (:)   = spval  
    allocate(this%dz_lake     (begc:endc,nlevlak))             ; this%dz_lake     (:,:) = nan
    allocate(this%z_lake      (begc:endc,nlevlak))             ; this%z_lake      (:,:) = nan

    allocate(this%glc_topo    (begc:endc))                     ; this%glc_topo    (:)   = nan
    allocate(this%micro_sigma (begc:endc))                     ; this%micro_sigma (:)   = nan
    allocate(this%n_melt      (begc:endc))                     ; this%n_melt      (:)   = nan 
    allocate(this%topo_slope  (begc:endc))                     ; this%topo_slope  (:)   = nan
    allocate(this%topo_std    (begc:endc))                     ; this%topo_std    (:)   = nan


    !DW  following section from SoilorderConType.F90
    allocate(this%smax           (0:nsoilorder))        ; this%smax(:)        =nan
    allocate(this%ks_sorption    (0:nsoilorder))        ; this%ks_sorption(:) =nan
    allocate(this%r_weather      (0:nsoilorder))        ; this%r_weather(:)   =nan
    allocate(this%r_adsorp       (0:nsoilorder))        ; this%r_adsorp(:)    =nan
    allocate(this%r_desorp       (0:nsoilorder))        ; this%r_desorp(:)    =nan
    allocate(this%r_occlude      (0:nsoilorder))        ; this%r_occlude(:)    =nan

    allocate(this%k_s1_biochem      (0:nsoilorder))        ; this%k_s1_biochem(:)    =nan
    allocate(this%k_s2_biochem      (0:nsoilorder))        ; this%k_s2_biochem(:)    =nan
    allocate(this%k_s3_biochem      (0:nsoilorder))        ; this%k_s3_biochem(:)    =nan
    allocate(this%k_s4_biochem      (0:nsoilorder))        ; this%k_s4_biochem(:)    =nan

    do m = 0,nsoilorder

       this%smax(m)         = smax(m)
       this%ks_sorption(m)         = ks_sorption(m)
       this%r_weather(m)         = r_weather(m)
       this%r_adsorp(m)         = r_adsorp(m)
       this%r_desorp(m)         = r_desorp(m)
       this%r_occlude(m)         = r_occlude(m)
       this%k_s1_biochem(m)         = k_s1_biochem(m)
       this%k_s2_biochem(m)         = k_s2_biochem(m)
       this%k_s3_biochem(m)         = k_s3_biochem(m)
       this%k_s4_biochem(m)         = k_s4_biochem(m)

    end do
    !DW   above section from SoilorderConType.F90

  end subroutine col_pp_init

  !------------------------------------------------------------------------
  subroutine col_pp_clean(this)
    !
    ! !ARGUMENTS:
    class(column_physical_properties_type) :: this
    !------------------------------------------------------------------------

    deallocate(this%gridcell   )
    deallocate(this%wtgcell    )
    deallocate(this%landunit   )
    deallocate(this%wtlunit    )
    deallocate(this%pfti       )
    deallocate(this%pftf       )
    deallocate(this%npfts      )
    deallocate(this%itype      )
    deallocate(this%active     )
    deallocate(this%snl        )
    deallocate(this%dz         )
    deallocate(this%z          )
    deallocate(this%zi         )
    deallocate(this%zii        )
    deallocate(this%lakedepth  )
    deallocate(this%dz_lake    )
    deallocate(this%z_lake     )
    deallocate(this%glc_topo   )
    deallocate(this%micro_sigma)
    deallocate(this%n_melt     )
    deallocate(this%topo_slope )
    deallocate(this%topo_std   )

    !DW Moved from SoilorderConType.F90
    deallocate(this%smax)
    deallocate(this%ks_sorption)
    deallocate(this%r_weather)
    deallocate(this%r_adsorp)
    deallocate(this%r_desorp)
    deallocate(this%r_occlude)
    deallocate(this%k_s1_biochem)
    deallocate(this%k_s2_biochem)
    deallocate(this%k_s3_biochem)
    deallocate(this%k_s4_biochem)

  end subroutine col_pp_clean


end module ColumnType

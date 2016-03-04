module GPM_SurfaceOpticsType
  

  !-----------------------------------------------------------------------
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use decompMod       , only : bounds_type
!  use clm_varpar     , only : numrad, nlevcan, nlevsno
  use GPM_clm_util_mod, only : GPM_MW_nchannels, &
                               GPM_DefaultEmissivity_MW,&
                               GPM_npol
  use abortutils      , only : endrun

  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC DATA MEMBERS:
  type, public :: gpm_surfaceoptics_type
    
    real(r8), allocatable :: emissivityMW_patch(:,:,:)   ! microwave emissivity calculated at patch level (0 to 1) [npatch * nchannel * 2]
  contains

    procedure, public  :: Init
    procedure, private :: InitAllocate
    procedure, private :: InitHistory
    procedure, private :: InitCold
    procedure, public  :: Restart

  end type gpm_surfaceoptics_type

  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(gpm_surfaceoptics_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    call this%InitCold(bounds)

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! Allocate module variables and data structures
    !
    ! !USES:
    use shr_infnan_mod, only: nan => shr_infnan_nan, assignment(=)
    use clm_varcon    , only: spval, ispval
    !
    ! !ARGUMENTS:
    class(gpm_surfaceoptics_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc

    allocate(this%emissivityMW_patch   (begp:endp, GPM_MW_nchannel, GPM_npol) )  ;    this%emissivityMW_patch (:,:,:) = nan
  end subroutine InitAllocate


  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! History fields initialization
    !
    ! !USES:
    use shr_infnan_mod, only: nan => shr_infnan_nan, assignment(=)
    use clm_varcon    , only: spval
    use histFileMod   , only: hist_addfld1d, hist_addfld2d
    !
    ! !ARGUMENTS:
    class(gpm_surfaceoptics_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc
  end subroutine InitHistory


  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! Initialize module surface optical properties to reasonable values
    !
    ! !ARGUMENTS:
    class(gpm_surfaceoptics_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc
    integer :: begp, endp
    !-----------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    this%emissivityMW_patch (begp:endp, :, :) = GPM_DefaultEmissivity_MW 
  end subroutine InitCold
  !---------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag, &
       tlai_patch, tsai_patch)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use clm_varctl , only : use_snicar_frc, iulog 
    use spmdMod    , only : masterproc
    use decompMod  , only : bounds_type
    use abortutils , only : endrun
    use ncdio_pio  , only : file_desc_t, ncd_defvar, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
    use restUtilMod

    ! !ARGUMENTS:
    class(gpm_surfaceoptics_type)               :: this
    type(bounds_type) , intent(in)    :: bounds 
    type(file_desc_t) , intent(inout) :: ncid ! netcdf id
    character(len=*)  , intent(in)    :: flag ! 'read' or 'write'
    real(r8)          , intent(in)    :: tlai_patch(bounds%begp:)
    real(r8)          , intent(in)    :: tsai_patch(bounds%begp:)
    !
    ! !LOCAL VARIABLES:
    logical :: readvar      ! determine if variable is on initial file
    integer :: iv
    integer :: begp, endp
    integer :: begc, endc
    !---------------------------------------------------------------------
  end subroutine Restart


end module GPM_SurfaceOpticsType

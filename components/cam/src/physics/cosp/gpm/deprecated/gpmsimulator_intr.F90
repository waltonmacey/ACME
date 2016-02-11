! GPM simulator intr 
!
! Jan. 26, 2016 Created by Yinghui Lu
!

module gpmsimulator_intr_mod
  !------------------
  ! Environment setup
  !------------------
  ! Module use
  use GPM_utility_mod, only: gpm_errorHandler
  use CRTM_Module, only: CRTM_ChannelInfo_Type, CRTM_ChannelInfo_Inspect
  ! Disable implicit typing
  implicit none


  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  private
  ! Private variables with save attribute, indicating which simulator is used
  logical, save, protected :: lgpmgmi_sim = .TRUE. 
  logical, save, protected :: lgpmdpr_sim = .FALSE.

  ! Private variable to indicate if the GPM simulator has been initialized
  logical, save :: linitialized = .FALSE.
  ! Private variable associated with GPM sensors
  integer, save ::  crtm_n_sensors = 0
  type(CRTM_ChannelInfo_type), allocatable, save :: crtm_chinfo(:)
  real, allocatable, save :: sensor_scan_angle_list(:)
  real, allocatable, save :: sensor_zenith_angle_list(:)

  ! Public parameters
  character(len=*), public, parameter :: modulename = &
                  'gpmsimulator_intr_mod'
  integer, public, parameter  :: n_gmi_ch = 13
  integer, public, parameter  :: gmi_chidx(n_gmi_ch) = &
                         (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13/)
  
  ! Public subroutines
  public :: gpmsimulator_intr_readnl, &
            gpmsimulator_intr_init,   &
            gpmsimulator_intr_run,    &
            gpmsimulator_intr_finalize

contains

!################################################################################
!################################################################################
!##                                                                            ##
!##                         ## PUBLIC MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################


!--------------------------------------------------------------------------------
! read namelist (not implemented for now)
!--------------------------------------------------------------------------------
subroutine gpmsimulator_intr_readnl()
  use namelist_utils, only: find_group_name
#ifdef SPMD
  use mpishorthand, only: mpicom, mpilog, mpiint, mpichar
#endif



end subroutine gpmsimulator_intr_readnl


!--------------------------------------------------------------------------------
! initialize 
!--------------------------------------------------------------------------------

subroutine gpmsimulator_intr_init()
  !------------------
  ! Environment setup
  !------------------
  ! Module use
  use GPM_CRTM_sensor_mod, only: GPM_CRTM_sensor_add, GPM_CRTM_sensor_init, &
                                 GPM_CRTM_sensor_inquire
  use cam_history,         only: addfld, add_default, outfld
  use cam_history_support, only: add_hist_coord
  use mod_cosp_constants,  only: R_UNDEF
  use GPM_CRTM_Constants,  only: max_hydro_types 

  ! Disable implicit typing
  implicit none
  character(len=*), parameter :: subroutinename = &
            'gpmsimulator_intr_init'  
  ! temporary cosp_histfile_num
  ! FIXME: implement a way to set this number
  integer, parameter :: tmp_cosp_histfile_num = 1
  !-----------------
  ! check if already initialized
  if (linitialized) then
    call gpm_errorHandler('GPM simulator is already initialized',& 
                          modulename, subroutinename)
  end if
    
  ! if run GPM GMI simulator
  if (lgpmgmi_sim) then
    ! Add GPM GMI sensor to senser list in GPM_CRTM_sensor_mod
    call GPM_CRTM_sensor_add('gpm-gmi-lowfreq')
    !call GPM_CRTM_sensor_add('gpm-gmi-highfreq')
    ! initialize GPM GMI sensors
    call GPM_CRTM_sensor_init()
    ! set n_sensors and chinfo
    call GPM_CRTM_sensor_inquire(crtm_chinfo, crtm_n_sensors,&
            sensor_scan_angle_list, sensor_zenith_angle_list)
    !------------------------
    ! Add to history fields
    ! Add coordinates
    call add_hist_coord('gmi_chidx',n_gmi_ch  ,'GPM GMI channel index', &
                   '1',gmi_chidx)
    ! add field
    call addfld('GPM_GMI_TB',(/'cosp_scol','gmi_chidx'/),'A','K',&
                   'simulated GPM GMI brightness temperature',&
                   flag_xyfill=.true., fill_value=R_UNDEF)
    ! add default calls
    ! FIXME: add to correct tape 
    call add_default('GPM_GMI_TB',tmp_cosp_histfile_num, ' ')
  end if

  ! if run GPM DPR simulator
  if (lgpmdpr_sim) then



  end if


  ! Set initialization state
  linitialized = .TRUE.
end subroutine gpmsimulator_intr_init

!--------------------------------------------------------------------------------
! run
!--------------------------------------------------------------------------------
subroutine gpmsimulator_intr_run(state, pbuf)
  !------------------
  ! Environment setup
  !------------------
  ! Module use
  use physics_types,  only: physics_state
  use physics_buffer, only: physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx
  use GPM_CRTM_sensor_mod, only: GPM_CRTM_sensor_inquire
  use MOD_COSP_TYPES, only: cosp_gridbox, construct_cosp_gridbox, free_cosp_gridbox, &
                            cosp_sghydro, construct_cosp_sghydro, free_cosp_sghydro
  use ppgrid,         only: pcols, pver, pverp, begchunk, endchunk  
  use GPM_CRTM_simulator_mod, only: GPM_CRTM_simulator_run, GPM_CRTM_result_type
  ! Disable implicit typing
  implicit none

  ! inputs
  type(physics_state), intent(in), target :: state
  type(physics_buffer_desc), pointer :: pbuf(:)
  ! parameters
  character(len=*),  parameter :: subroutinename = &
            'gpmsimulator_intr_run'   
  ! local variables
  type(GPM_CRTM_result_type), allocatable :: gpm_results(:)
  integer :: i_sensor
  
  ! gridbox and subgrid information
  type(cosp_gridbox) :: gpmgbx
  type(COSP_SGHYDRO) :: gpmsgx

  integer :: Npoints
  integer :: Nlevels
  ! local variables used to convert CAM states to COSP inputs
  real :: longitude(pcols)
  real :: latitude(pcols)
  real :: zlev(pcols,pver)
  real :: zint(pcols,pverp)
  real :: dlev(pcols,pver)
  real :: p(pcols,pver)
  real :: pint(pcols,pverp)
  real :: T(pcols,pver)
  real :: q(pcols,pver)
  real :: sh(pcols,pver)
  real :: mr_ozone(pcols,pver)
  real :: land(pcols)
  real :: psfc(pcols)
  real :: mr_hydro(pcols,pver,max_hydro_types)
  real :: Reff(pcols,pver,max_hydro_types)
  real :: Np(pcols,pver,max_hydro_types)
   
  integer :: lchnk     ! chunk idendifier
  integer :: ncol      ! number of active atmospheric columns

  ! --------- end of variable declaration ----------!
!  real, allocatable :: tbs(:,:)
  !---------- checks -----------
  ! check if already initialized
  if (.NOT. linitialized) then
    call gpm_errorHandler('GPM simulator is not initialized yet, cannot run simulation',&
                           modulename, subroutinename)
  end if
  ! ------------------
  ! find the chunk and ncol from the state vector
  lchnk = state%lchnk  ! chunk number 
  ncol  = state%ncol   ! number of columns in the chunk

  






  ! run calculations
  call gpm_crtm_simulator_run(gbx, crtm_chinfo,sensor_scan_angle_list, sensor_zenith_angle_list, gpm_results)
  ! inspect the results
  do i_sensor = 1, crtm_n_sensors
    print *, gpm_results(i_sensor)%tbs
  end do


end subroutine gpmsimulator_intr_run


!--------------------------------------------------------------------------------
! finalize
!--------------------------------------------------------------------------------
subroutine gpmsimulator_intr_finalize()
  !------------------
  ! Environment setup
  !------------------
  ! Module use
  use GPM_CRTM_sensor_mod, only: GPM_CRTM_sensor_destroy
  ! Disable implicit typing
  implicit none
  character(len=*),  parameter :: subroutinename = &
            'gpmsimulator_intr_finalize'   
  ! check if already initialized
  if (.NOT. linitialized) then
    call gpm_errorHandler('GPM simulator is not initialized yet, cannot finalize',&
                           modulename, subroutinename)
  end if
  ! destroy the sensor structures
  if (lgpmgmi_sim) then
    call GPM_CRTM_sensor_destroy() 
    lgpmgmi_sim  = .FALSE.
  end if

  if (lgpmdpr_sim) then
    lgpmdpr_sim  = .FALSE.
  end if
  ! deallocate crtm_chinfo
  deallocate(crtm_chinfo)
  ! reset initialization flag and instrument flag
  linitialized = .FALSE.
end subroutine gpmsimulator_intr_finalize


!=============================================
end module gpmsimulator_intr_mod

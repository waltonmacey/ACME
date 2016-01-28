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
  logical, save :: lgpmgmi_sim = .TRUE. 
  logical, save :: lgpmdpr_sim = .FALSE.

  ! Private variable to indicate if the GPM simulator has been initialized
  logical, save :: linitialized = .FALSE.
  ! Private variable associated with GPM sensors
  integer, save ::  crtm_n_sensors = 0
  type(CRTM_ChannelInfo_type), allocatable, save :: crtm_chinfo(:)

  ! Public parameters
  character(len=*), public, parameter :: modulename = &
                  'gpmsimulator_intr_mod'
  
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
  ! Disable implicit typing
  implicit none
  character(len=*), parameter :: subroutinename = &
            'gpmsimulator_intr_init'   
  ! check if already initialized
  if (linitialized) then
    call gpm_errorHandler('GPM simulator is already initialized',& 
                          modulename, subroutinename)
  end if
    
  ! if run GPM GMI simulator
  if (lgpmgmi_sim) then
    ! Add GPM GMI sensor to senser list in GPM_CRTM_sensor_mod
    call GPM_CRTM_sensor_add('gpm-gmi-lowfreq')
    call GPM_CRTM_sensor_add('gpm-gmi-highfreq')
    ! initialize GPM GMI sensors
    call GPM_CRTM_sensor_init()
    ! set n_sensors and chinfo
    call GPM_CRTM_sensor_inquire(crtm_chinfo, crtm_n_sensors)

    ! Add to history fields
    !TODO: implement add to history fields

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
subroutine gpmsimulator_intr_run(gbx)
  !------------------
  ! Environment setup
  !------------------
  ! Module use
  use GPM_CRTM_sensor_mod, only: GPM_CRTM_sensor_inquire
  use MOD_COSP_TYPES, only: cosp_gridbox
  ! Disable implicit typing
  implicit none

  ! inputs
  type(cosp_gridbox), intent(in) :: gbx
  ! parameters
  character(len=*),  parameter :: subroutinename = &
            'gpmsimulator_intr_run'   
  !---------------------
  ! check if already initialized
  if (.NOT. linitialized) then
    call gpm_errorHandler('GPM simulator is not initialized yet, cannot run simulation',&
                           modulename, subroutinename)
  end if
  ! run calculations
  call gpm_crtm_simulator_run(gbx, crtm_chinfo,gpm_result)
  


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
  ! reset initialization flag and instrument flag
  linitialized = .FALSE.
end subroutine gpmsimulator_intr_finalize


!=============================================
end module gpmsimulator_intr_mod

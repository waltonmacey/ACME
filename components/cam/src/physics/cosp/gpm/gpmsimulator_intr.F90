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
  ! Disable implicit typing
  implicit none


  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  private
  ! Private variables with save attribute, indicating which simulator is used
  logical, save :: lgpmgmi_sim = .FALSE. 
  logical, save :: lgpmdpr_sim = .FALSE.

  ! Private variable to indicate if the GPM simulator has been initialized
  logical, save :: linitialized = .FALSE.

  ! Public parameters
  character(len=*), save, public, parameter :: module_name = &
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
subroutine gpmsimulator_intr_readnl(nlfile)
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
  use GPM_CRTM_sensor_mod, only: GPM_CRTM_sensor_add, GPM_CRTM_sensor_init
  ! Disable implicit typing
  implicit none
  character(len=*), save, parameter :: subroutine_name = &
            'gpmsimulator_intr_init'   
  ! check if already initialized
  if (linitialized) then
    call gpm_errorHandler('GPM simulator is already initialized',& 
                          module_name, subroutine_name)
  end if
    
  ! if run GPM GMI simulator
  if (lgpmgmi_sim) then
    ! Add GPM GMI sensor to senser list in GPM_CRTM_sensor_mod
    call GPM_CRTM_sensor_add('gpm-gmi-lowfreq')
    call GPM_CRTM_sensor_add('gpm-gmi-highfreq')
    ! initialize GPM GMI sensors
    call GPM_CRTM_sensor_init()
  
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
subroutine gpmsimulator_intr_run()
  !------------------
  ! Environment setup
  !------------------
  ! Module use
  use GPM_CRTM_sensor_mod, only: GPM_CRTM_sensor_inquire
  ! Disable implicit typing
  implicit none
  character(len=*), save, parameter :: subroutine_name = &
            'gpmsimulator_intr_inquire'   
  ! check if already initialized
  if (.NOT. linitialized) then
    call gpm_errorHandler('GPM simulator is not initialized yet, cannot run simulation',&
                           module_name, subroutine_name)
  end if


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
  character(len=*), save, parameter :: subroutine_name = &
            'gpmsimulator_intr_finalize'   
  ! check if already initialized
  if (.NOT. linitialized) then
    call gpm_errorHandler('GPM simulator is not initialized yet, cannot finalize',&
                           module_name, subroutine_name)
  end if
  ! destroy the sensor structures
  call GPM_CRTM_sensor_destroy() 

  ! reset initialization flag and instrument flag
  lgpmdpr_sim  = .FALSE.
  lgpmgmi_sim  = .FALSE.
  linitialized = .FALSE.
end subroutine gpmsimulator_intr_finalize


!=============================================
end module gpmsimulator_intr_mod

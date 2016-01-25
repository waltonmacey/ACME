!
! GPM_CRTM_CONSTANTS
!
! Module defining the constants used in GPM CRTM
!
!
! CREATION HISTORY:
!       Written by:     Yinghui Lu, Jan.-22-2016
!                       yinghuilu@lbl.gov
!

module GPM_CRTM_Constants

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  use cam_history_support, only: max_fieldname_len
  use CRTM_Parameters, only: STRLEN
  ! Disable implicit typing
  implicit none

  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  private
  ! Public attribute will be set on declaration

  ! -----------------
  ! Module parameters
  ! -----------------
  ! maximum length of the CAM history field name
  integer, parameter, public :: maxlen_camhistfld_name = max_fieldname_len
  ! maximum length of the sensor ID used in CRTM
  integer, parameter, public :: maxlen_sensorid = STRLEN
  ! maximum number of sensors calculated using CRTM
  integer, parameter, public :: max_CRTM_sensors = 10
  ! Directory containing all the CRTM coefficients
  character(100), parameter, public :: CRTM_dircoef = '/global/homes/y/yxl232/local/CRTM/coeff_data/Big_Endian_ODAS/'

end module GPM_CRTM_Constants

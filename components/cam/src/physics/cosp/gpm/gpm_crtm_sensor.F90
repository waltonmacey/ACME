!
!  Code related to GPM GMI sensors
!
module GPM_CRTM_sensor_mod
  !------------------
  ! Environment setup
  !------------------
  ! Module use
  use CRTM_Module, only: crtm_channelInfo_type, CRTM_Init, CRTM_Destroy, &
                         CRTM_ChannelInfo_Inspect, SUCCESS, &
                         CRTM_ChannelInfo_n_Channels, CRTM_ChannelInfo_Subset
  use GPM_CRTM_Constants, only: maxlen_camhistfld_name, maxlen_sensorid, &
                                max_CRTM_sensors, CRTM_dircoef
  
  ! Disable implicit typing
  implicit none
  
  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  private
  ! Public attribute will be set on declaration
!  public :: GPM_CRTM_sensor_type
!  public :: GPM_CRTM_sensor_create, GPM_CRTM_sensor_createdefault
  public :: GPM_CRTM_sensor_add, GPM_CRTM_sensor_init, GPM_CRTM_sensor_destroy, &
            GPM_CRTM_sensor_inquire
  ! ---------------------
  ! Procedure overloading
  ! ---------------------
  interface GPM_CRTM_sensor_add
    module procedure GPM_CRTM_sensor_add_byproperties
    module procedure GPM_CRTM_sensor_add_byname
   end interface GPM_CRTM_sensor_add
  ! -------------------------------
  ! GPM_CRTM_sensor_type is a type that saves the properties of the sensors
  ! -------------------------------
  !:tdoc+:
  type GPM_CRTM_sensor_type
    character(maxlen_sensorid) :: sensor_id = ''                 ! sensor ID string
    character(maxlen_camhistfld_name) :: cam_histfld_name = ''   ! cam history field name of results
    real :: sensor_scan_angle    = 0                       ! sensor zenith angle (degree)
    real :: sensor_zenith_angle  = 0                       ! sensor zenith angle (degree)
    integer, allocatable :: channel_subset(:)             ! channel subset
    logical :: l_channelsubset   = .FALSE.

    !! The following variables are set (and allocated if applicabable) 
    !! in function GPM_CRTM_init()
    ! number of channels of the sensor, for convenience
    integer :: n_channels        = 0

  end type GPM_CRTM_sensor_type
  !:tdoc-:
  ! ---------------
  ! saved variables 
  ! ---------------
  ! Arrays of sensors and sensor informations
  type(GPM_CRTM_sensor_type),  save :: sensor_list(max_CRTM_sensors)
  character(maxlen_sensorid),  save :: sensor_id_list(max_CRTM_sensors)
  type(CRTM_ChannelInfo_type), save :: chinfo_list(max_CRTM_sensors)
  ! logical variable indicating if the sensors are initialized
  logical, save :: l_initialized = .FALSE.
  ! current number of sensors added
  integer, save :: n_sensors = 0
  
contains
!################################################################################
!################################################################################
!##                                                                            ##
!##                         ## PUBLIC MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################
!--------------------------------------------------------------------------------
! Initialize GPM CRTM shared variables
!--------------------------------------------------------------------------------

   subroutine GPM_CRTM_sensor_init()
     ! local variables
    integer :: err_stat
    integer :: n  
     !--------------
     
      ! Check if GPM GMI structures are already initialized
      if ( l_initialized ) then
         print *, "GPM GMI structures are already initialized. Cannot initialize &
            again"
         stop
      end if

      ! Check if there are already sensors added. If no sensors added, stop
      ! running 
      if (n_sensors <= 0) then
         print *, "No sensors added. Add sensor first then run GPM_CRTM_init()"
         stop
      end if

      !------------------------
      ! Copy values of sensor_id from each sensor
      do n = 1, n_sensors
         sensor_id_list(n) = sensor_list(n)%sensor_id
      end do 
      ! Initialize chinfo
      err_stat = CRTM_Init(sensor_id_list(1:n_sensors), chinfo_list(1:n_sensors), &
                           File_Path = CRTM_dircoef, &
                           Quiet = .TRUE.)
      if (err_stat /= SUCCESS) then
         write(*,*) "fail to initialize chinfo"
         STOP
      endif
      l_initialized = .TRUE.
      
         ! subset channels if necessary
         ! currently, whether a sensor needs channel subsetting is determined by
         ! whether its channel_subset properties is set.
      do n=1, n_sensors
         if (allocated(sensor_list(n)%channel_subset)) then
            err_stat = CRTM_ChannelInfo_Subset(chinfo_list(n), &
                           Channel_Subset = sensor_list(n)%channel_subset )
            if (err_stat /= SUCCESS) then
               print *, "failed to subset channels"
               stop
            end if
         end if
         ! allocate arrays
         sensor_list(n)%n_channels = CRTM_ChannelInfo_n_Channels( chinfo_list(n) )
!         allocate( currentsensor%toa_bt(currentsensor%n_channels, n_profiles))
      end do !n_sensors

   end subroutine GPM_CRTM_sensor_init
!--------------------------------------------------------------------------------
! Destroy GPM CRTM shared variables
!--------------------------------------------------------------------------------
   subroutine GPM_CRTM_sensor_destroy()
      ! local variables
      integer :: n
      integer :: err_stat
      !--------------
      if (.NOT. l_initialized ) then
         print *, "GPM CRTM shared variables are not initilized. Cannot destroy."
         stop
      end if
      do n= 1, n_sensors
         sensor_id_list(n) = ''
      end do
      call GPM_CRTM_sensor_clear(sensor_list)
      err_stat =  CRTM_Destroy(chinfo_list)
      if (err_stat .NE. SUCCESS) then
         print *, "Error occurred when destroying CRTM_ChannelInfo arrays"
         stop
      end if
      n_sensors = 0
      l_initialized = .FALSE.


   end subroutine GPM_CRTM_sensor_destroy
!--------------------------------------------------------------------------------
! Inquire GPM CRTM shared variables
!--------------------------------------------------------------------------------
   subroutine GPM_CRTM_sensor_inquire(chinfo_list_out, n_sensors_out)
!  type(GPM_CRTM_sensor_type),        intent(out), optional :: sensor_list(:)
!  character(maxlen_sensorid),  intent(out), optional :: sensor_id_list(:)
  type(CRTM_ChannelInfo_type), intent(out), optional :: chinfo_list_out(:)
  ! current number of sensors added
  integer, intent(out), optional :: n_sensors_out
  ! ---------
  ! check if the GPM CRTM shared variable has been initialized
  if (.NOT. l_initialized) then
     print *, "GPM CRTM shared variables are not initialized yet. Cannot &
              inquire their values"
     stop
  end if
  ! inquire value for chinfo
  if (present(chinfo_list_out) ) then
     ! check if the sizes is consistant
     print *, size(chinfo_list_out), n_sensors
     if(size(chinfo_list_out) .NE. n_sensors ) then
        print *, "output variable size is not consistant with n_sensors"
     end if
     chinfo_list_out = chinfo_list(1:n_sensors)
  end if
  ! inquire n_sensors
  if (present(n_sensors_out) ) then
     n_sensors_out = n_sensors
  end if


   end subroutine GPM_CRTM_sensor_inquire

!################################################################################
!################################################################################
!##                                                                            ##
!##                         ## PRIVATE MODULE ROUTINES ##                      ##
!##                                                                            ##
!################################################################################
!################################################################################
!--------------------------------------------------------------------------------
! Add a sensor by its properties
!--------------------------------------------------------------------------------
  subroutine GPM_CRTM_sensor_add_byproperties(sensor_id, cam_histfld_name,&
                      sensor_scan_angle, sensor_zenith_angle, channel_subset )
      character(len=*), intent(in) :: sensor_id
      character(len=*), intent(in) :: cam_histfld_name
      real, intent(in) :: sensor_scan_angle
      real, intent(in) :: sensor_zenith_angle
      integer, intent(in), optional :: channel_subset(:)
      !------------
    if (n_sensors >= max_CRTM_sensors) then
       print *, "already reach maximum number of sensors. Cannot add any more"
       stop
    end if
    if (l_initialized) then
       print *, "GPM GMI shared variables are already initialized. Cannot &
                  add any more sensors"
       stop
    end if
      n_sensors = n_sensors +1
      call GPM_CRTM_sensor_create(sensor_list(n_sensors),sensor_id, cam_histfld_name,&
                sensor_scan_angle, sensor_zenith_angle, channel_subset) 
    
  end subroutine GPM_CRTM_sensor_add_byproperties
!--------------------------------------------------------------------------------
! Add a sensor by its name
!--------------------------------------------------------------------------------
  subroutine GPM_CRTM_sensor_add_byname(sensor_name)
    character(len=*), intent(in) :: sensor_name
    !----------------
    if (n_sensors >= max_CRTM_sensors) then
       print *, "already reach maximum number of sensors. Cannot add any more"
       stop
    end if
    if (l_initialized) then
       print *, "GPM GMI shared variables are already initialized. Cannot &
                  add any more sensors"
       stop
    end if
      n_sensors = n_sensors +1
     
      call GPM_CRTM_sensor_createdefault(sensor_list(n_sensors),sensor_name)
  end subroutine GPM_CRTM_sensor_add_byname
  

!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       GPM_CRTM_sensor_create
!
! PURPOSE:
!       Create GPM CRTM sensor
!
! CALLING SEQUENCE:
!       call GPM_CRTM_sensor_create( Atm )
!
! OBJECTS:
!       Atm:       Atmosphere structure which is to have its member's
!                  status tested.
!                  UNITS:      N/A
!                  TYPE:       CRTM_Atmosphere_type
!                  DIMENSION:  Scalar or any rank
!                  ATTRIBUTES: INTENT(IN)
!
! FUNCTION RESULT:
!       Status:    The return value is a logical value indicating the
!                  status of the Atmosphere members.
!                    .TRUE.  - if the array components are allocated.
!                    .FALSE. - if the array components are not allocated.
!                  UNITS:      N/A
!                  TYPE:       LOGICAL
!                  DIMENSION:  Same as input
!
!:sdoc-:
!--------------------------------------------------------------------------------
! TODO: hide the constructor of GPM_CRTM_sensor_type
! TODO: finish changing comments and documentation
   subroutine GPM_CRTM_sensor_create(gmi, sensor_id, cam_histfld_name,&
                sensor_scan_angle, sensor_zenith_angle, channel_subset) 
   
   !--------------------------
   !
   !  Add a sensor and update n_sensors
   !
   !-------------------------
      type(GPM_CRTM_sensor_type), intent(out) :: gmi
      character(len=*), intent(in) :: sensor_id
      character(len=*), intent(in) :: cam_histfld_name
      real, intent(in) :: sensor_scan_angle
      real, intent(in) :: sensor_zenith_angle
      integer, intent(in), optional :: channel_subset(:)
      !-------------------
      integer :: i
      !-------------------
      gmi%sensor_id = sensor_id
      gmi%cam_histfld_name = cam_histfld_name
      gmi%sensor_scan_angle = sensor_scan_angle
      gmi%sensor_zenith_angle = sensor_zenith_angle
      if ( present(channel_subset) ) then
         gmi%n_channels = size(channel_subset)
         allocate(gmi%channel_subset(gmi%n_channels) )
         gmi%channel_subset = channel_subset
         gmi%l_channelsubset = .TRUE.

      else
         gmi%l_channelsubset = .FALSE.
         gmi%n_channels = -1
      end if
      
   end subroutine GPM_CRTM_sensor_create
!======================================
! TODO: use xml file to store default sensor properties
! The current implementation hard-wired all the prpoerties
   subroutine GPM_CRTM_sensor_createdefault(gmi, sensorname)
!-----------------------
!
! create a default sensor
!
!-----------------------
      type(GPM_CRTM_sensor_type) :: gmi      
      character(len=*), intent(in) :: sensorname
      !--- local variable ---
      integer :: i
      !----------
      select case (sensorname)
         case ('gpm-gmi-lowfreq')
            call GPM_CRTM_sensor_create(gmi, 'gmi_gpm','GPM_GMI_LOWFREQ',52.8, 48.5,&
                      channel_subset = (/(i,i=1,9)/)  )
         case ('gpm-gmi-highfreq')
            call GPM_CRTM_sensor_create(gmi, 'gmi_gpm','GPM_GMI_LOWFREQ',49.19, 45.36,&
                      channel_subset = (/(i,i=10,13)/)  )
         case default
            print *, 'error: sensor not implemented in default list'
            stop
      end select
   end subroutine GPM_CRTM_sensor_createdefault
!============================================

!----------------------------------
!  clear information in sensor
!----------------------------------
   elemental subroutine GPM_CRTM_sensor_clear(sensor)
      type(GPM_CRTM_sensor_type), intent(inout) :: sensor
      !---------------
      sensor%sensor_id = ''
      sensor%cam_histfld_name = ''
      sensor%sensor_scan_angle = 0
      sensor%sensor_zenith_angle = 0
      if ( sensor%l_channelsubset ) then
         deallocate(sensor%channel_subset )
      end if
         sensor%l_channelsubset = .FALSE.
         sensor%n_channels = -1

   end subroutine GPM_CRTM_sensor_clear

end module GPM_CRTM_sensor_mod

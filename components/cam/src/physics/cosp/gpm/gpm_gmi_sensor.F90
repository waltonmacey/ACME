!
!  Code related to GPM GMI sensors
!
module gpm_gmi_sensor_mod
   !------------------
   ! Environment setup
   !------------------
   ! Module use
   use cam_history_support, only: max_fieldname_len
   implicit none
   private

!--------------------------------
   ! gpm_gmi_sensor is a type that saves the properties of the sensors
   type gpm_gmi_sensor
      character(20) :: sensor_id                          ! sensor ID string
      character(max_fieldname_len) :: cam_histfld_name    ! cam history field name of results
      real :: sensor_scan_angle                           ! sensor zenith angle (degree)
      real :: sensor_zenith_angle                         ! sensor zenith angle (degree)
      integer, allocatable :: channel_subset(:)           ! channel subset
      logical :: l_channelsubset

      !! The following variables are set (and allocated if applicabable) 
      !! in function gpm_gmi_init()
      ! number of channels of the sensor, for convenience
      integer :: n_channels

   end type gpm_gmi_sensor

!----------------------------------
   public :: gpm_gmi_sensor
   public :: gpm_gmi_sensor_create, gpm_gmi_sensor_createdefault


!----------------------------------
contains
!====================================
! TODO: hide the constructor of gpm_gmi_sensor
   function gpm_gmi_sensor_create(sensor_id, cam_histfld_name,&
                sensor_scan_angle, sensor_zenith_angle, channel_subset) &
                result(gmi)
   
   !--------------------------
   !
   !  Add a sensor and update n_sensors
   !
   !-------------------------
      character(len=*), intent(in) :: sensor_id
      character(len=*), intent(in) :: cam_histfld_name
      real, intent(in) :: sensor_scan_angle
      real, intent(in) :: sensor_zenith_angle
      integer, intent(in), optional :: channel_subset(:)
      !-------------------
      integer :: i
      type(gpm_gmi_sensor) :: gmi
      gmi%sensor_id = sensor_id
      gmi%cam_histfld_name = cam_histfld_name
      gmi%sensor_scan_angle = sensor_scan_angle
      gmi%sensor_zenith_angle = sensor_zenith_angle
      if ( present(channel_subset) ) then
         gmi%n_channels = size(channel_subset)
         allocate(gmi%channel_subset(gmi%n_channels) )
         gmi%l_channelsubset = .TRUE.
      else
         gmi%l_channelsubset = .FALSE.
         gmi%n_channels = -1
      end if
      
   end function gpm_gmi_sensor_create
!======================================
! TODO: use xml file to store default sensor properties
! The current implementation hard-wired all the prpoerties
   function gpm_gmi_sensor_createdefault(sensorname) result(gmi)
!-----------------------
!
! create a default sensor
!
!-----------------------
      character(len=*), intent(in) :: sensorname
      !--- return variable --
      type(gpm_gmi_sensor) :: gmi      
      !--- local variable ---
      integer :: i
      !----------
      select case (sensorname)
         case ('gpm-gmi-lowfreq')
            gmi = gpm_gmi_sensor_create('gmi_gpm','GPM_GMI_LOWFREQ',52.8, 48.5,&
                      channel_subset = (/(i,i=1,9)/)  )
         case ('gpm-gmi-highfreq')
            gmi = gpm_gmi_sensor_create('gmi_gpm','GPM_GMI_LOWFREQ',49.19, 45.36,&
                      channel_subset = (/(i,i=10,13)/)  )
         case default
            print *, 'error: sensor not implemented in default list'
            stop
      end select
   end function gpm_gmi_sensor_createdefault
!============================================



end module gpm_gmi_sensor_mod

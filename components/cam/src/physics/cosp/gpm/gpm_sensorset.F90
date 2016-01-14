!
! Code related to all GPM sensors
!
module gpm_sensorset_mod
   ! module use
   use gpm_gmi_sensor_mod
   

   implicit none
   private
!----------------
   integer, parameter :: gpm_gmi_max_sensors = 10 ! maximum number of gpm-gmi sensors 

!---------------
   type gpm_sensorset
      integer gpm_gmi_n_sensors = 0 ! current number of GMI sensors
      type(gpm_gmi_sensor) :: gpm_gmi_sensors(gpm_gmi_max_sensors) ! GMI sensor array


   end type gpm_sensorset
!---------------
   public :: gpm_sensorset
   public :: gpm_addsensor_gpmgmi

contains
   subroutine gpm_addsensor_gpmgmi(sensor)
   !---------------
   !
   ! Add a GPM-GMI sensor
   !
   !---------------
   type(gpm_gmi_sensor), intent(in) :: sensor
   if (gpm_gmi_n_sensors >= gpm_gmi_max_sensors) then
      print *, "number of GPM GMI sensors exceeds maximum number"
      stop
   end if
   gpm_gmi_n_sensors = gpm_gmi_n_sensors +1;
   gpm_gmi_sensors(gpm_gmi_n_sensors) = sensor

   end subroutine gpm_addsensor_gpmgmi


end module gpm_sensorset_mod

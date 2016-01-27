! GPM utilities
!
! Jan. 27, 2016 Created by Yinghui Lu
!
module GPM_utility_mod
  implicit none
  private

  ! Public subroutines
  public :: gpm_errorHandler

contains

!################################################################################
!################################################################################
!##                                                                            ##
!##                         ## PUBLIC MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################


!--------------------------------------------------------------------------------
! handle errors
!--------------------------------------------------------------------------------
  subroutine gpm_errorHandler(errinfo, module_name, subroutine_name )
    ! Inputs and outputs
    character(len=*), intent(in) :: errinfo
    character(len=*), intent(in) :: module_name
    character(len=*), intent(in) :: subroutine_name


    ! modify the code below to handle error
    ! For now, just print out the information associated with the error and stop
    ! running
    print *, "error in: ", module_name, " :: ",&
             subroutine_name, " : ", errinfo
    stop
  end subroutine gpm_errorHandler

end module GPM_utility_mod

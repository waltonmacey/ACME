! GPM GMI result type definition
! 
! Jan. 22, 2016 Created by Yinghui Lu
! Jan. 22, 2016: move code related to GPM Results here from
!                gpm_crtm_simulator_mod
module gpm_crtm_result_mod
   implicit none
   private
   public:: gpm_crtm_result,          &
            gpm_crtm_result_init,     &
            gpm_crtm_result_destroy,  &
            gpm_crtm_result_inquire,  &
            gpm_crtm_result_set 
   
   type gpm_crtm_result
   !----------------------
   !
   !  output for GPM CRTM
   !
   !--------------------
      private ! set all components to private for safety
      integer :: n_channels = -1
      integer :: n_profiles = -1
      real, allocatable :: tbs(:,:) ! n_channels x n_profiles

   end type gpm_crtm_result

contains
   
!-----------------------
!
! Initialize GPM result structure
!
!----------------------
   subroutine gpm_crtm_result_init(n_channels, n_profiles, gpm)
      integer, intent(in) :: n_channels
      integer, intent(in) :: n_profiles
      type(gpm_crtm_result), intent(out) :: gpm
      gpm%n_channels = n_channels
      gpm%n_profiles = n_profiles
      allocate(gpm%tbs(n_channels, n_profiles ) )
   end subroutine gpm_crtm_result_init
!-----------------------
!
! Destroy GPM CRTM result structure 
!
!----------------------
   subroutine gpm_crtm_result_destroy(gpm)
      type(gpm_crtm_result), intent(inout) :: gpm
      gpm%n_channels = -1
      gpm%n_profiles = -1
      if(allocated(gpm%tbs) ) then
         deallocate(gpm%tbs)
      endif
   end subroutine gpm_crtm_result_destroy

!-----------------------
!
! Obtain information in gpm_crtm_result structure.
! It is needed because the variables in gpm_crtm_result are private for safety
!
!----------------------
   subroutine gpm_crtm_result_inquire(gpm, n_channels, n_profiles, tb)
      type(gpm_crtm_result), intent(in) :: gpm
      integer, intent(out), optional :: n_channels
      integer, intent(out), optional :: n_profiles
      real,    intent(out), optional :: tb(:,:)
      if(present(n_channels)) then
         n_channels = gpm%n_channels
      endif
      if(present(n_profiles)) then
         n_profiles = gpm%n_profiles
      endif
      if(present(tb)) then
         tb = gpm%tbs
      endif
   end subroutine gpm_crtm_result_inquire

!-----------------------
!
! Obtain information in gpm_crtm_result structure.
! It is needed because the variables in gpm_crtm_result are private for safety
!
!----------------------
   subroutine gpm_crtm_result_set(gpm, tbs)
      type(gpm_crtm_result), intent(inout) :: gpm
      real, intent(in) :: tbs(:,:)
      ! check if the size is the same
      if ( (size(tbs,1) .NE. gpm%n_channels) .OR. (size(tbs,2) .NE. gpm%n_profiles) ) then
         print *, " the size of TB is not correct "
         stop
      end if

      gpm%tbs = tbs
   end subroutine gpm_crtm_result_set

end module gpm_crtm_result_mod

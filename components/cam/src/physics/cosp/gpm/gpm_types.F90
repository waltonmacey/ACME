module MOD_GPM_types

!=========== type definitions===============!

!----------- GPM Gridbox -------------------!

  ! Input data for GPM simulator.
 type GPM_gridbox
    ! scalars and dimensions
    integer :: Npoints   ! # of gridpoints
    integer :: Nlevels   ! # of levels
    integer :: Ncolumns  ! # of subcolumns
    integer :: Nhydro    ! # of types of hydrometers
    integer :: Naero     ! # of types of aerosols

    ! date, time, and geolocation
    real :: time
    real :: time_bnds(2)
    real, allocatable :: longitude(:) ! longitude [degree East]
    real, allocatable :: latitude(:)  ! latitude  [degree North]

    ! gridbox information. Default dimension: [npoints x nlevels]
    real, allocatable :: zlev(:,:)    ! height of model level [m] 
    real, allocatable :: zint(:,:)    ! height of model level interface [m]
                                      ! [npoints x (nlevels+1)]
    real, allocatable :: dlev(:,:)    ! depth of model level
    real, allocatable :: p(:,:)       ! pressure at model level [Pa]
    real, allocatable :: pint(:,:)    ! pressure at model level interface [Pa]
                                      ! [npoints x (nlevels+1)]
    real, allocatable :: T(:,:)       ! Temperature at model level (K)
    real, allocatable :: q(:,:)       ! Relative humidity of water (%)
    real, allocatable :: sh(:,:)      ! Specific humidity of water (kg/kg)
    real, allocatable :: mr_ozone(:,:)! mixing ratio of ozone (kg/kg)

    ! Point information. dimension: [npoints]
    real, allocatable :: land(:)      ! land mask
    real, allocatable :: psfc(:)      ! surface pressure

    ! Hydrometeor information of each type [npoints x nlevels x Nhydro]
    real, allocatable :: mr_hydro(:,:,:)  ! mixing ratio (kg/kg)
    real, allocatable :: Reff(:,:,:)      ! effective radius (m)
    real, allocatable :: Np(:,:,:)        ! number concentration

  end type GPM_gridbox

contains
!============ subroutines =================!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE CONSTRUCT_GPM_GRIDBOX ------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine construct_GPM_gridbox(time, time_bnds, &
                  Npoints, Nlevels, Ncolumns, Nhydro, Naero, &
                  y)
    !-------- inputs and outputs -----!
    ! date, time, and geolocation
    real, intent(in) :: time
    real, intent(in) :: time_bnds(2)
    ! dimensional variables
    integer, intent(in) :: Npoints   ! # of gridpoints
    integer, intent(in) :: Nlevels   ! # of levels
    integer, intent(in) :: Ncolumns  ! # of subcolumns
    integer, intent(in) :: Nhydro    ! # of types of hydrometers
    integer, intent(in) :: Naero     ! # of types of aerosols

    type(GPM_gridbox), intent(out) :: y


    !--------- assign values -----------!
    y%Npoints  = Npoints
    y%Nlevels  = Nlevels
    y%Ncolumns = Ncolumns
    y%Nhydro   = Nhydro
    y%Naero    = Naero
    y%time     = time
    y%time_bnds = time_bnds
    
    !--------- allocate arrays ---------!
    allocate(y%longitude(Npoints), y%latitude(Npoints))
    allocate(y%zlev(Npoints, Nlevels), y%zint(Npoints, Nlevels+1), &
             y%dlev(Npoints, Nlevels),    y%p(Npoints, Nlevels),   &
             y%pint(Npoints, Nlevels+1),  y%T(Npoints, Nlevels),   &
                y%q(Npoints, Nlevels),   y%sh(Npoints, Nlevels),   &
         y%mr_ozone(Npoints, Nlevels) )
    allocate(y%land(Npoints), y%psfc(Npoints) )
    allocate(y%mr_hydro(Npoints, Nlevels, Nhydro) , &
                 y%Reff(Npoints, Nlevels, Nhydro) , &
                   y%Np(Npoints, Nlevels, Nhydro)  )
    !--------- initialize arrays to zero ---!
    y%longitude  = 0 
    y%latitude   = 0 
                  
    y%zlev       = 0 
    y%zint       = 0 
    
    y%dlev       = 0 
    y%p          = 0 
    y%pint       = 0 
                
    y%T          = 0 
    y%q          = 0
    y%sh         = 0 
    y%mr_ozone   = 0 
                  
    y%land       = 0
    y%psfc       = 0 
                  
    y%mr_hydro   = 0
    y%Reff       = 0  
    y%Np         = 0  

  end subroutine construct_GPM_gridbox

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE CONSTRUCT_GPM_GRIDBOX ------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine free_GPM_gridbox(y)
    type(GPM_gridbox), intent(inout) :: y
    
    y%Npoints    = 0
    y%Nlevels    = 0
    y%Ncolumns   = 0
    y%Nhydro     = 0
    y%Naero      = 0
    y%time       = 0
    y%time_bnds  = 0
    
    deallocate(y%longitude, y%latitude)
    deallocate(y%zlev, y%zint, y%dlev, y%p, y%pint, y%T, y%q, y%sh, y%mr_ozone )
    deallocate(y%land, y%psfc )
    deallocate(y%mr_hydro , y%Reff ,  y%Np  )


  end subroutine free_GPM_gridbox

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE GPM_GRIDBOX_CPSECTION ------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine GPM_gridbox_cpsection(ix, iy, x, y)
    integer, intent(in) :: ix(2)
    integer, intent(in) :: iy(2)
    type(GPM_gridbox), intent(in)    :: x
    type(GPM_gridbox), intent(inout) :: y

    ! FIXME: consider including bound check
    ! copy 1-D data
    y%longitude(iy(1):iy(2))   = x%longitude(ix(1):ix(2))
    y%latitude( iy(1):iy(2))   = x%latitude( ix(1):ix(2))
    y%land(     iy(1):iy(2))   = x%land(     ix(1):ix(2))
    y%psfc(     iy(1):iy(2))   = x%psfc(     ix(1):ix(2))

    ! copy 2-D data
    y%zlev(    iy(1):iy(2),:)   = x%zlev(    ix(1):ix(2),:)
    y%zint(    iy(1):iy(2),:)   = x%zint(    ix(1):ix(2),:)
    y%dlev(    iy(1):iy(2),:)   = x%dlev(    ix(1):ix(2),:)
    y%p(       iy(1):iy(2),:)   = x%p(       ix(1):ix(2),:)
    y%pint(    iy(1):iy(2),:)   = x%pint(    ix(1):ix(2),:)
    y%T(       iy(1):iy(2),:)   = x%T(       ix(1):ix(2),:)
    y%q(       iy(1):iy(2),:)   = x%q(       ix(1):ix(2),:)
    y%sh(      iy(1):iy(2),:)   = x%sh(      ix(1):ix(2),:)
    y%mr_ozone(iy(1):iy(2),:)   = x%mr_ozone(ix(1):ix(2),:)

    ! copy 3-D data
    y%mr_hydro(iy(1):iy(2),:,:) = x%mr_hydro(ix(1):ix(2),:,:)
    y%Reff(    iy(1):iy(2),:,:) = x%Reff(    ix(1):ix(2),:,:)
    y%Np(      iy(1):iy(2),:,:) = x%Np(      ix(1):ix(2),:,:)

  end subroutine GPM_gridbox_cpsection
end module MOD_GPM_TYPES

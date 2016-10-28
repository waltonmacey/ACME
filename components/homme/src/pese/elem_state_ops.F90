!
! routines performing simple operations on all prognostic variables
!_______________________________________________________________________
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module elem_state_ops

  ! subroutines for performing operations on the element state structure

  use control_mod,      only: statefreq
  use derivative_mod,   only: derivative_t ! divergence_sphere, gradient_sphere
  use dimensions_mod,   only: np, nlev, nlevp, qsize, qsize_d
  use edge_mod,         only: initEdgeBuffer, edgevpack, edgevunpack
  use edgetype_mod,     only: edgebuffer_t
  use element_mod,      only: element_t, elem_state_t, derived_state_t
  use hybvcoord_mod,    only: hvcoord_t
  use hybrid_mod,       only: hybrid_t
  use kinds,            only: rl => real_kind
  use physical_constants, only: g,p0
  use shr_const_mod,    only: unset => shr_const_spval
  use vertical_se,      only: eta_derivative, elem_height, evenly_spaced_eta_coords, vertical_interp_matrix, vertical_dss

  implicit none

  real(rl), allocatable:: s_interp(:), M_interp(:,:)

contains

  !_____________________________________________________________________
  subroutine apply_vertical_dss(elem, nt)

    ! apply vertical dss to all prognostics and tracers

    type (element_t), intent(inout), target :: elem
    integer, intent(in) :: nt
    integer::l

    call vertical_dss( elem%state%T (:,:,  :,nt) )
    call vertical_dss( elem%state%v (:,:,1,:,nt) )
    call vertical_dss( elem%state%v (:,:,2,:,nt) )

    do l=1,qsize_d
      call vertical_dss(elem%state%Qdp(:,:,:,l,nt))
    enddo

  end subroutine

  !_____________________________________________________________________
  subroutine pack_edge_data(buffer, elem, nt, ie)

    ! pack prognostics into edge buffer

    type (EdgeBuffer_t),  intent(inout)       :: buffer                 ! buffer for edge-data exchange
    type (element_t),     intent(in), target  :: elem                   ! array of element_t structures
    integer,              intent(in)          :: nt, ie                 ! time level and element index

    integer :: i
    i = 0
    call edgeVpack(buffer, elem%state%ps_v(:,:,    nt), 1     , i, ie); i=i+1
    call edgeVpack(buffer, elem%state%T   (:,:,:,  nt), nlev  , i, ie); i=i+nlev
    call edgeVpack(buffer, elem%state%v   (:,:,:,:,nt), 2*nlev, i, ie); i=i+2*nlev

  end subroutine

  !_____________________________________________________________________
  subroutine unpack_edge_data(buffer, elem, nt, ie)

    ! unpack prognostics from edge buffer

    type (EdgeBuffer_t),  intent(inout)       :: buffer                 ! buffer for edge-data exchange
    type (element_t),     intent(inout), target  :: elem                ! array of element_t structures
    integer,              intent(in)          :: nt, ie                 ! time level and element index

    integer :: i
    i = 0
    call edgeVunpack(buffer, elem%state%ps_v(:,:,    nt), 1     , i, ie); i=i+1
    call edgeVunpack(buffer, elem%state%T   (:,:,:,  nt), nlev  , i, ie); i=i+nlev
    call edgeVunpack(buffer, elem%state%v   (:,:,:,:,nt), 2*nlev, i, ie); i=i+2*nlev

  end subroutine

  !_____________________________________________________________________
  subroutine apply_map(M, elem, nt)

    ! apply pointwise 2d map to each prognostic

    real(rl),         intent(in)            :: M(np,np)
    type (element_t), intent(inout), target :: elem                     ! array of element_t structures
    integer,          intent(in)            :: nt                       ! time level

    integer :: k
    elem%state%ps_v(:,:, nt) = elem%state%ps_v(:,:,nt)*M
    do k=1,nlev
      elem%state%T(:,:,k,  nt) = elem%state%T(:,:,  k,nt)*M
      elem%state%v(:,:,1,k,nt) = elem%state%v(:,:,1,k,nt)*M
      elem%state%v(:,:,2,k,nt) = elem%state%v(:,:,2,k,nt)*M
    enddo

  end subroutine

  !_____________________________________________________________________
  subroutine display_max_and_min(elem, hybrid, ie, nt, count)

      ! display verbose diagnostics

      type (element_t),   intent(inout), target :: elem                 ! element
      type (hybrid_t),    intent(in)		:: hybrid												! mpi/omp data struct
      integer,            intent(in)    :: ie                           ! element number
      integer,            intent(in)    :: nt                           ! time level to display
      integer,            intent(in)    :: count                        ! iteration counter

      real(rl):: maxu,minu,maxT,minT,maxv,minv,maxps,minps

      if( mod(count, statefreq)==0 ) then

        maxu = maxval(elem%state%v   (:,:,1,:,nt))
        minu = minval(elem%state%v   (:,:,1,:,nt))
        maxv = maxval(elem%state%v   (:,:,2,:,nt))
        minv = minval(elem%state%v   (:,:,2,:,nt))
        maxps= maxval(elem%state%ps_v(:,:,    nt))
        minps= minval(elem%state%ps_v(:,:,    nt))
        maxT = maxval(elem%state%T   (:,:,:,  nt))
        minT = minval(elem%state%T   (:,:,:,  nt))

        if (hybrid%masterthread) &
        print *,"nt =",nt, "ie =",ie, " max u =",maxu," max v =",maxv," max T =",maxT," max ps=",maxps

        if (any(isnan(elem%state%v(:,:,1,:,nt))))  stop 'detected NaN in u'
        if (any(isnan(elem%state%v(:,:,2,:,nt))))  stop 'detected NaN in v'
        if (any(isnan(elem%state%T(:,:,:,nt))))    stop 'detected NaN in T'
        if (any(isnan(elem%state%ps_v(:,:,nt))))   stop 'detected NaN in ps_v'

      endif

  end subroutine

 !_____________________________________________________________________________
  function v_interpolate(f,ni) result(f_i)

    real(rl), intent(in) :: f(nlev)
    integer,  intent(in) :: ni
    real(rl) :: f_i(ni)
    f_i = matmul(M_interp,f)

  end function

  !_____________________________________________________________________________
  function get_scalar_field(element, short_name, n0, ni) result(sfield)

    ! Get vertically interpolated 3d scalar-field data by name

    character*(*),    intent(in)  :: short_name
    type(element_t),  intent(in)  :: element
    integer,          intent(in)  :: n0       ! time level index
    integer,          intent(in)  :: ni       ! num vertical interp pts

    real(rl) :: sfield(np,np,ni)
    real(rl) :: var(np,np,nlev),var2d(np,np)
    integer i,j

    if(ni==1) then

      ! get 2d scalar field by name
      select case(short_name)
        case('ps'  ); var2d = element%state%ps_v(:,:,n0)
        case default; var2d = unset
      end select
      sfield(:,:,1) = var2d
    endif

    if(ni>1) then

      ! get 3d scalar field by name
      select case(short_name)

        ! prognostic variables
        case('T'  ); var = element%state%T(:,:,:,n0)
        case('u'  ); var = element%state%v(:,:,1,:,n0)
        case('v'  ); var = element%state%v(:,:,2,:,n0)
        case('p'  ); var = element%state%dp3d(:,:,:,n0)
        case('Q'  ); var = element%state%Q(:,:,:,1)
        case('Q2' ); var = element%state%Q(:,:,:,2)
        case('Q3' ); var = element%state%Q(:,:,:,3)
        case('Q4' ); var = element%state%Q(:,:,:,4)
        case('Q5' ); var = element%state%Q(:,:,:,5)

        case('omega');var = element%state%omega

        ! diagnostic variables
        case('geo'); var = element%derived%phi(:,:,:)
        case default; var = unset                                      ! assign special "missing" value
      endselect

      ! interpolate each column
      do i=1,np
        do j=1,np
          sfield(i,j,:) = v_interpolate(var(i,j,:),ni)
        enddo
      enddo
    endif

  end function

  !_____________________________________________________________________________
  function get_vector_field(element, short_name, n0, ni) result(vfield)

    ! Get vertically interpolated vector-field data by name

    character*(*),    intent(in)  :: short_name
    type(element_t),  intent(in)  :: element
    integer,          intent(in)  :: n0       ! time level index
    integer,          intent(in)  :: ni       ! num vertical interp pts

    real(rl) :: vfield(np,np,2,ni)
    real(rl) :: var(np,np,2,nlev)
    integer  :: i,j
    select case(short_name)                                             ! switch on short variable name

      case('v');
        var = element%state%v(:,:,:,:,n0)

      case default; var = unset                                      ! assign special "missing" value
    endselect

    ! interpolate each column
    do i=1,np
      do j=1,np
        vfield(i,j,1,:) = v_interpolate(var(i,j,1,:),ni)
        vfield(i,j,2,:) = v_interpolate(var(i,j,2,:),ni)
      enddo
    enddo

  end function

  !_____________________________________________________________________________
  function get_vertical_levels(ni) result(levels)

    ! Get vertically interpolated levels

    integer, intent(in) :: ni   ! number of interpolated levels
    real(rl) :: levels(ni)

    ! allocate vertical interpolation points and interpolation matrix

    if(.not. allocated(s_interp)) then
      allocate(s_interp(ni))
      allocate(M_interp(ni,nlev))
      s_interp = evenly_spaced_eta_coords(ni)
      M_interp = vertical_interp_matrix(s_interp,ni)
    endif

    levels = s_interp

  end function

end module
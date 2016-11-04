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

  real(rl), allocatable:: s_interp(:), M_interp(:,:) !  pts, and projection matrix for vertical interpolation

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

end module
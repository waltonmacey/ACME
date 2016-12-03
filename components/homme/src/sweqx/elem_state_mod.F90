#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module elem_state_mod

  ! inherited variables
  use elem_state_base, only: &
    statecomponents,        &
    timelevels

  use kinds,          only: real_kind, long_kind, int_kind
  use dimensions_mod, only: np, nc, npsq, nlev, nlevp, qsize_d, max_neigh_edges

  implicit none
  public

  type, public :: elem_state_t

    ! prognostic variables for shallow-water solver
     real (kind=real_kind) :: p(np,np,nlev,timelevels)
     real (kind=real_kind) :: ps(np,np)                               ! surface geopotential
     real (kind=real_kind) :: gradps(np,np,2)                         ! gradient of surface geopotential
     real (kind=real_kind) :: v(np,np,2,nlev,timelevels)              ! contravarient comp

  end type elem_state_t

  !___________________________________________________________________
  type, public :: derived_state_t
     real (kind=real_kind) :: vstar(np,np,2,nlev)                     ! velocity on Lagrangian surfaces
  end type derived_state_t

end module elem_state_mod


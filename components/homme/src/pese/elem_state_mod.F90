!
! Element Module for the PESE dynamics solver
!_______________________________________________________________________
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module elem_state_mod

  use elem_state_base, only: &
    derived_state_t,         &
    timelevels

  use dimensions_mod, only: np, nlev, nlevp, qsize_d
  use kinds,          only: real_kind, long_kind, int_kind

  implicit none
  public

  !_____________________________________________________________________
  type, public :: elem_state_t

    real (kind=real_kind) :: v   (np,np,2,nlev,timelevels)                ! horizontal velocity
    real (kind=real_kind) :: T   (np,np,nlev,timelevels)                  ! temperature
    real (kind=real_kind) :: lnps(np,np,timelevels)                       ! log surface pressure
    real (kind=real_kind) :: ps_v(np,np,timelevels)                       ! surface pressure
    real (kind=real_kind) :: phis(np,np)                                  ! surface geopotential (prescribed)
    real (kind=real_kind) :: Q   (np,np,nlev,qsize_d)                     ! tracer concentration
    real (kind=real_kind) :: Qdp (np,np,nlev,qsize_d,2)                   ! tracer mass
    real (kind=real_kind) :: dp3d(np,np,nlev,timelevels)                  ! delta p on levels

    real (kind=real_kind) :: O1  (np,np,nlev,timelevels)                  ! diagnostic output fields
    real (kind=real_kind) :: O2  (np,np,nlev,timelevels)                  ! diagnostic output fields
    real (kind=real_kind) :: O3  (np,np,nlev,timelevels)                  ! diagnostic output fields
    real (kind=real_kind) :: O4  (np,np,nlev,timelevels)                  ! diagnostic output fields

    real (kind=real_kind) :: omega(np,np,nlev)                            ! dp/dt

  end type elem_state_t

  ! num vars in state struct used for standalone HOMME restart file
  integer(kind=int_kind),public,parameter :: StateComponents=8

end module elem_state_mod


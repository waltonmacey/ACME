! element data structures for primitive equatios solver PREQx
!_______________________________________________________________________
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module elem_state_base

  use kinds,                  only: real_kind, long_kind, int_kind
  use coordinate_systems_mod, only: spherical_polar_t, cartesian2D_t, cartesian3D_t, distance
  use dimensions_mod,         only: np, nc, npsq, nlev, nlevp, qsize_d, max_neigh_edges
  use edgetype_mod,           only: edgedescriptor_t, rotation_t
  use gridgraph_mod,          only: gridvertex_t

  implicit none
  public

  integer, public, parameter :: timelevels = 3

  type, public :: elem_state_t

    ! prognostic variables for preqx solver

    ! prognostics must match those in prim_restart_mod.F90
    ! vertically-lagrangian code advects dp3d instead of ps_v
    ! tracers Q, Qdp always use 2 level time scheme

    real (kind=real_kind) :: v   (np,np,2,nlev,timelevels)            ! velocity                           1
    real (kind=real_kind) :: T   (np,np,nlev,timelevels)              ! temperature                        2
    real (kind=real_kind) :: dp3d(np,np,nlev,timelevels)              ! delta p on levels                  8
    real (kind=real_kind) :: lnps(np,np,timelevels)                   ! log surface pressure               3
    real (kind=real_kind) :: ps_v(np,np,timelevels)                   ! surface pressure                   4
    real (kind=real_kind) :: phis(np,np)                              ! surface geopotential (prescribed)  5
    real (kind=real_kind) :: Q   (np,np,nlev,qsize_d)                 ! Tracer concentration               6
    real (kind=real_kind) :: Qdp (np,np,nlev,qsize_d,2)               ! Tracer mass                        7

  end type elem_state_t

  integer(kind=int_kind),public,parameter::StateComponents=8! num prognistics variables (for prim_restart_mod.F90)

  !___________________________________________________________________
  type, public :: derived_state_t

    ! diagnostic variables for preqx solver

    ! storage for subcycling tracers/dynamics

    real (kind=real_kind) :: vn0  (np,np,2,nlev)                      ! velocity for SE tracer advection
    real (kind=real_kind) :: vstar(np,np,2,nlev)                      ! velocity on Lagrangian surfaces
    real (kind=real_kind) :: dpdiss_biharmonic(np,np,nlev)            ! mean dp dissipation tendency, if nu_p>0
    real (kind=real_kind) :: dpdiss_ave(np,np,nlev)                   ! mean dp used to compute psdiss_tens

    ! diagnostics for explicit timestep
    real (kind=real_kind) :: phi(np,np,nlev)                          ! geopotential
    real (kind=real_kind) :: omega_p(np,np,nlev)                      ! vertical tendency (derived)
    real (kind=real_kind) :: eta_dot_dpdn(np,np,nlevp)                ! mean vertical flux from dynamics

    ! semi-implicit diagnostics: computed in explict-component, reused in Helmholtz-component.
    real (kind=real_kind) :: grad_lnps(np,np,2)                       ! gradient of log surface pressure
    real (kind=real_kind) :: zeta(np,np,nlev)                         ! relative vorticity
    real (kind=real_kind) :: div(np,np,nlev,timelevels)               ! divergence

    ! tracer advection fields used for consistency and limiters
    real (kind=real_kind) :: dp(np,np,nlev)                           ! for dp_tracers at physics timestep
    real (kind=real_kind) :: divdp(np,np,nlev)                        ! divergence of dp
    real (kind=real_kind) :: divdp_proj(np,np,nlev)                   ! DSSed divdp

#ifdef CAM
    ! forcing terms for CAM
    real (kind=real_kind) :: FQ(np,np,nlev,qsize_d, 1)                ! tracer forcing
    real (kind=real_kind) :: FM(np,np,2,nlev, 1)                      ! momentum forcing
    real (kind=real_kind) :: FT(np,np,nlev, 1)                        ! temperature forcing
    real (kind=real_kind) :: etadot_prescribed(np,np,nlevp)           ! prescribed vertical tendency
    real (kind=real_kind) :: u_met(np,np,nlev)                        ! zonal component of prescribed meteorology winds
    real (kind=real_kind) :: dudt_met(np,np,nlev)                     ! rate of change of zonal component of prescribed meteorology winds
    real (kind=real_kind) :: v_met(np,np,nlev)                        ! meridional component of prescribed meteorology winds
    real (kind=real_kind) :: dvdt_met(np,np,nlev)                     ! rate of change of meridional component of prescribed meteorology winds
    real (kind=real_kind) :: T_met(np,np,nlev)                        ! prescribed meteorology temperature
    real (kind=real_kind) :: dTdt_met(np,np,nlev)                     ! rate of change of prescribed meteorology temperature
    real (kind=real_kind) :: ps_met(np,np)                            ! surface pressure of prescribed meteorology
    real (kind=real_kind) :: dpsdt_met(np,np)                         ! rate of change of surface pressure of prescribed meteorology
    real (kind=real_kind) :: nudge_factor(np,np,nlev)                 ! nudging factor (prescribed)
    real (kind=real_kind) :: Utnd(npsq,nlev)                          ! accumulated U tendency due to nudging towards prescribed met
    real (kind=real_kind) :: Vtnd(npsq,nlev)                          ! accumulated V tendency due to nudging towards prescribed met
    real (kind=real_kind) :: Ttnd(npsq,nlev)                          ! accumulated T tendency due to nudging towards prescribed met
#else
    ! forcing terms for HOMME
    real (kind=real_kind) :: FQ(np,np,nlev,qsize_d, timelevels)       ! tracer forcing
    real (kind=real_kind) :: FM(np,np,2,nlev, timelevels)             ! momentum forcing
    real (kind=real_kind) :: FT(np,np,nlev, timelevels)               ! temperature forcing
#endif

    ! forcing terms for both CAM and HOMME
    ! FQps for conserving dry mass in the presence of precipitation

    real (kind=real_kind) :: pecnd(np,np,nlev)                        ! pressure perturbation from condensate
    real (kind=real_kind) :: FQps(np,np,timelevels)                   ! forcing of FQ on ps_v

  end type derived_state_t

end module elem_state_base


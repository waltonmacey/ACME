#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module element_mod

  use kinds,                  only: real_kind, long_kind, int_kind
  use coordinate_systems_mod, only: spherical_polar_t, cartesian2D_t, cartesian3D_t, distance
  use dimensions_mod,         only: np, npsq, nlev, nlevp, qsize_d, max_neigh_edges
  use edge_mod,               only: edgedescriptor_t, rotation_t
  use gridgraph_mod,          only: gridvertex_t
#ifdef _FVM
  use dimensions_mod,         only: nc, nhc, ntrac_d
#endif

  implicit none
  private
  integer, public, parameter :: timelevels = 3

#ifdef _PRIM

! =========== PRIMITIVE-EQUATION DATA-STRUCTURES =====================

#if USE_OPENACC

  public :: setup_element_pointers

  real (kind=real_kind), allocatable, target, public :: state_v   (:,:,:,:,:,:) !(np,np,2,nlev,timelevels,nelemd)
  real (kind=real_kind), allocatable, target, public :: state_T   (:,:,:,:,:)   !(np,np,nlev,timelevels,nelemd)  
  real (kind=real_kind), allocatable, target, public :: state_dp3d(:,:,:,:,:)   !(np,np,nlev,timelevels,nelemd)  
  real (kind=real_kind), allocatable, target, public :: state_lnps(:,:,:,:)     !(np,np,timelevels,nelemd)       
  real (kind=real_kind), allocatable, target, public :: state_ps_v(:,:,:,:)     !(np,np,timelevels,nelemd)       
  real (kind=real_kind), allocatable, target, public :: state_phis(:,:,:)       !(np,np,nelemd)                  
  real (kind=real_kind), allocatable, target, public :: state_Q   (:,:,:,:,:)   !(np,np,nlev,qsize_d,nelemd)     
  real (kind=real_kind), allocatable, target, public :: state_Qdp (:,:,:,:,:,:) !(np,np,nlev,qsize_d,2,nelemd)   
  type, public :: elem_state_t
    ! prognostic variables for preqx solver
    ! prognostics must match those in prim_restart_mod.F90
    ! vertically-lagrangian code advects dp3d instead of ps_v
    ! tracers Q, Qdp always use 2 level time scheme
    real (kind=real_kind), pointer :: v   (:,:,:,:,:)  ! velocity                           1  (np,np,2,nlev,timelevels)
    real (kind=real_kind), pointer :: T   (:,:,:,:)    ! temperature                        2  (np,np,nlev,timelevels)  
    real (kind=real_kind), pointer :: dp3d(:,:,:,:)    ! delta p on levels                  8  (np,np,nlev,timelevels)  
    real (kind=real_kind), pointer :: lnps(:,:,:)      ! log surface pressure               3  (np,np,timelevels)       
    real (kind=real_kind), pointer :: ps_v(:,:,:)      ! surface pressure                   4  (np,np,timelevels)       
    real (kind=real_kind), pointer :: phis(:,:)        ! surface geopotential (prescribed)  5  (np,np)                  
    real (kind=real_kind), pointer :: Q   (:,:,:,:)    ! Tracer concentration               6  (np,np,nlev,qsize_d)     
    real (kind=real_kind), pointer :: Qdp (:,:,:,:,:)  ! Tracer mass                        7  (np,np,nlev,qsize_d,2)   
  end type elem_state_t

  integer(kind=int_kind),public,parameter::StateComponents=8  ! num prognistics variables (for prim_restart_mod.F90)

  real (kind=real_kind), allocatable, target, public :: derived_vn0              (:,:,:,:,:)      ! (np,np,2,nlev,nelemd)                   velocity for SE tracer advection
  real (kind=real_kind), allocatable, target, public :: derived_vstar            (:,:,:,:,:)      ! (np,np,2,nlev,nelemd)                   velocity on Lagrangian surfaces
  real (kind=real_kind), allocatable, target, public :: derived_dpdiss_biharmonic(:,:,:,:)        ! (np,np,nlev,nelemd)                     mean dp dissipation tendency, if nu_p>0
  real (kind=real_kind), allocatable, target, public :: derived_dpdiss_ave       (:,:,:,:)        ! (np,np,nlev,nelemd)                     mean dp used to compute psdiss_tens
  real (kind=real_kind), allocatable, target, public :: derived_phi              (:,:,:,:)        ! (np,np,nlev,nelemd)                     geopotential
  real (kind=real_kind), allocatable, target, public :: derived_omega_p          (:,:,:,:)        ! (np,np,nlev,nelemd)                     vertical tendency (derived)       
  real (kind=real_kind), allocatable, target, public :: derived_eta_dot_dpdn     (:,:,:,:)        ! (np,np,nlevp,nelemd)                    mean vertical flux from dynamics
  real (kind=real_kind), allocatable, target, public :: derived_grad_lnps        (:,:,:,:)        ! (np,np,2,nelemd)                        gradient of log surface pressure               
  real (kind=real_kind), allocatable, target, public :: derived_zeta             (:,:,:,:)        ! (np,np,nlev,nelemd)                     relative vorticity                             
  real (kind=real_kind), allocatable, target, public :: derived_div              (:,:,:,:,:)      ! (np,np,nlev,timelevels,nelemd)          divergence                          
  real (kind=real_kind), allocatable, target, public :: derived_dp               (:,:,:,:)        ! (np,np,nlev,nelemd)                     for dp_tracers at physics timestep
  real (kind=real_kind), allocatable, target, public :: derived_divdp            (:,:,:,:)        ! (np,np,nlev,nelemd)                     divergence of dp
  real (kind=real_kind), allocatable, target, public :: derived_divdp_proj       (:,:,:,:)        ! (np,np,nlev,nelemd)                     DSSed divdp
#ifdef CAM
  real (kind=real_kind), allocatable, target, public :: derived_FQ               (:,:,:,:,:,:)    ! (np,np,nlev,qsize_d,1,nelemd)           tracer forcing  
  real (kind=real_kind), allocatable, target, public :: derived_FM               (:,:,:,:,:,:)    ! (np,np,2,nlev,1,nelemd)                 momentum forcing
  real (kind=real_kind), allocatable, target, public :: derived_FT               (:,:,:,:,:)      ! (np,np,nlev,1,nelemd)                   temperature forcing
  real (kind=real_kind), allocatable, target, public :: derived_omega_prescribed (:,:,:,:)        ! (np,np,nlev,nelemd)                     prescribed vertical tendency
#else
  real (kind=real_kind), allocatable, target, public :: derived_FQ               (:,:,:,:,:,:)    ! (np,np,nlev,qsize_d,timelevels,nelemd)  tracer forcing 
  real (kind=real_kind), allocatable, target, public :: derived_FM               (:,:,:,:,:,:)    ! (np,np,2,nlev,timelevels,nelemd)        momentum forcing
  real (kind=real_kind), allocatable, target, public :: derived_FT               (:,:,:,:,:)      ! (np,np,nlev,timelevels,nelemd)          temperature forcing 
#endif
  real (kind=real_kind), allocatable, target, public :: derived_pecnd            (:,:,:,:)        ! (np,np,nlev,nelemd)                     pressure perturbation from condensate
  real (kind=real_kind), allocatable, target, public :: derived_FQps             (:,:,:,:)        ! (np,np,timelevels,nelemd)               forcing of FQ on ps_v 
  type, public :: derived_state_t
    ! diagnostic variables for preqx solver
    ! storage for subcycling tracers/dynamics
    ! if (compute_mean_flux==1) vn0=time_avg(U*dp) else vn0=U at tracer-time t
  real (kind=real_kind), pointer :: vn0              (:,:,:,:)       ! (np,np,2,nlev)                  velocity for SE tracer advection
  real (kind=real_kind), pointer :: vstar            (:,:,:,:)       ! (np,np,2,nlev)                  velocity on Lagrangian surfaces
  real (kind=real_kind), pointer :: dpdiss_biharmonic(:,:,:)         ! (np,np,nlev)                    mean dp dissipation tendency, if nu_p>0
  real (kind=real_kind), pointer :: dpdiss_ave       (:,:,:)         ! (np,np,nlev)                    mean dp used to compute psdiss_tens
  real (kind=real_kind), pointer :: phi              (:,:,:)         ! (np,np,nlev)                    geopotential
  real (kind=real_kind), pointer :: omega_p          (:,:,:)         ! (np,np,nlev)                    vertical tendency (derived)       
  real (kind=real_kind), pointer :: eta_dot_dpdn     (:,:,:)         ! (np,np,nlevp)                   mean vertical flux from dynamics
  real (kind=real_kind), pointer :: grad_lnps        (:,:,:)         ! (np,np,2)                       gradient of log surface pressure               
  real (kind=real_kind), pointer :: zeta             (:,:,:)         ! (np,np,nlev)                    relative vorticity                             
  real (kind=real_kind), pointer :: div              (:,:,:,:)       ! (np,np,nlev,timelevels)         divergence                          
  real (kind=real_kind), pointer :: dp               (:,:,:)         ! (np,np,nlev)                    for dp_tracers at physics timestep
  real (kind=real_kind), pointer :: divdp            (:,:,:)         ! (np,np,nlev)                    divergence of dp
  real (kind=real_kind), pointer :: divdp_proj       (:,:,:)         ! (np,np,nlev)                    DSSed divdp
#ifdef CAM
  real (kind=real_kind), pointer :: FQ               (:,:,:,:,:)     ! (np,np,nlev,qsize_d,1)          tracer forcing  
  real (kind=real_kind), pointer :: FM               (:,:,:,:,:)     ! (np,np,2,nlev,1)                momentum forcing
  real (kind=real_kind), pointer :: FT               (:,:,:,:)       ! (np,np,nlev,1)                  temperature forcing
  real (kind=real_kind), pointer :: omega_prescribed (:,:,:)         ! (np,np,nlev)                    prescribed vertical tendency
#else
  real (kind=real_kind), pointer :: FQ               (:,:,:,:,:)     ! (np,np,nlev,qsize_d,timelevels) tracer forcing 
  real (kind=real_kind), pointer :: FM               (:,:,:,:,:)     ! (np,np,2,nlev,timelevels)       momentum forcing
  real (kind=real_kind), pointer :: FT               (:,:,:,:)       ! (np,np,nlev,timelevels)         temperature forcing 
#endif
  real (kind=real_kind), pointer :: pecnd            (:,:,:)         ! (np,np,nlev)                    pressure perturbation from condensate
  real (kind=real_kind), pointer :: FQps             (:,:,:)         ! (np,np,timelevels)              forcing of FQ on ps_v 
  end type derived_state_t

!else for USE_OPENACC if
#else

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
    ! if (compute_mean_flux==1) vn0=time_avg(U*dp) else vn0=U at tracer-time t

    real (kind=real_kind) :: vn0(np,np,2,nlev)                        ! velocity for SE tracer advection
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
    real (kind=real_kind) :: omega_prescribed(np,np,nlev)             ! prescribed vertical tendency
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

!ending USE_OPENACC if
#endif

  !___________________________________________________________________
  type, public :: elem_accum_t

#ifdef ENERGY_DIAGNOSTICS

    ! Energy equation:
    ! KE_t  = T1 + T2  + D1   + Err   +  vertical & horizontal advection terms
    ! IE_t  = S1 + D2                 +  vertical & horizontal advection terms
    ! PE_t  = S2        
    !
    ! KEvert*  =  KE net vertical advection    (should be zero)
    ! KEhoriz* =  KE net horizonatl advection  (should be zero)
    ! IEvert*  =  IE net vertical advection    (should be zero)
    ! IEhoriz* =  IE net horizonatl advection  (should be zero)
    !
    ! With leapfrog, energy equations are all exact except KE 
    ! (has an Err term that goes to zero as dt**2)
    !
    ! Transfer terms:
    ! T1   = -< dp/dn u, RT_v/p grad_p >     KE<->IE:   T1 + T2-T2_s = S1
    ! T2   = -< dp/dn u, grad_phi >          KE<->PE:   T2_s         = S2
    ! T2_s = -< dp/dn u, grad_phis >
    ! S1   = < Cp_star dp/dn , RT omega_p/Cp_star >  
    ! S2   = -< div (u dp/dn), phis >                

    real (kind=real_kind) :: KEvert1(np,np)                           ! term from continuity equ
    real (kind=real_kind) :: KEvert2(np,np)                           ! term from momentum equ
    real (kind=real_kind) :: IEvert1(np,np)                           ! term from continuity equ
    real (kind=real_kind) :: IEvert2(np,np)                           ! term from T equ
    real (kind=real_kind) :: IEvert1_wet(np,np)                       ! wet term from continuity equ
    real (kind=real_kind) :: IEvert2_wet(np,np)                       ! wet term from T equ

    real (kind=real_kind) :: KEhorz1(np,np)                           ! at time t
    real (kind=real_kind) :: KEhorz2(np,np)                           ! after calling time_advance, these will be at time t-1
    real (kind=real_kind) :: IEhorz1(np,np)    
    real (kind=real_kind) :: IEhorz2(np,np)    
    real (kind=real_kind) :: IEhorz1_wet(np,np)    
    real (kind=real_kind) :: IEhorz2_wet(np,np)    

    real (kind=real_kind) :: T1(np,np)    
    real (kind=real_kind) :: T2(np,np) 
    real (kind=real_kind) :: T2_s(np,np)    
    real (kind=real_kind) :: S1(np,np)    
    real (kind=real_kind) :: S1_wet(np,np)    
    real (kind=real_kind) :: S2(np,np)    

    ! the KE conversion term and diffusion term
    real (kind=real_kind) :: DIFF(np,np,2,nlev)                       ! net hypervis term
    real (kind=real_kind) :: DIFFT(np,np,nlev)                        ! net hypervis term
    real (kind=real_kind) :: CONV(np,np,2,nlev)                       ! dpdn u dot CONV = T1 + T2
#endif

    ! the "4" timelevels represents data computed at:
    !  1  t-.5   
    !  2  t+.5   after dynamics
    !  3  t+.5   after forcing
    !  4  t+.5   after Robert
    ! after calling TimeLevelUpdate, all times above decrease by 1.0

    real (kind=real_kind) :: KEner(np,np,4)       
    real (kind=real_kind) :: PEner(np,np,4)       
    real (kind=real_kind) :: IEner(np,np,4)    
    real (kind=real_kind) :: IEner_wet(np,np,4)    
    real (kind=real_kind) :: Qvar(np,np,qsize_d,4)                    ! Q variance at half time levels   
    real (kind=real_kind) :: Qmass(np,np,qsize_d,4)                   ! Q mass at half time levels
    real (kind=real_kind) :: Q1mass(np,np,qsize_d)                    ! Q mass at full time levels

  end type elem_accum_t

#elif defined _PRIMDG

! ============ DISCONTINUOUS-GALERKIN DATA-STRUCTURES ================

  type, public :: elem_state_t

    ! prognostic variables for DG primitive-eqn solver

    real (kind=real_kind) :: p(np,np,nlev,timelevels)
    real (kind=real_kind) :: phis(np,np)                              ! surface geopotential
    real (kind=real_kind) :: gradps(np,np,2)                          ! gradient of surface geopotential
    real (kind=real_kind) :: v(np,np,2,nlev,timelevels)               ! contravarient comp

    real (kind=real_kind) :: couv(np,np,2,nlev)                       ! covariant velocities
    real (kind=real_kind) :: uv(np,np,2,nlev)                         ! redundant copy of v (eliminate?)
    real (kind=real_kind) :: uv0(np,np,2,nlev)                        ! temp variable for velocities (eliminate?)
    real (kind=real_kind) :: pgrads(np,np,2,nlev)

    real (kind=real_kind) :: psi(np,np,nlev)
    real (kind=real_kind) :: phi(np,np,nlev)
    real (kind=real_kind) :: ht(np,np,nlev)
    real (kind=real_kind) :: T(np,np,nlev,timelevels)                 ! temperature
    real (kind=real_kind) :: Q(np,np,nlev,timelevels)                 ! moisture
    real (kind=real_kind) :: pt3d(np,np,nlev)                         ! potential temperature
    real (kind=real_kind) :: qt3d(np,np,nlev)
    real (kind=real_kind) :: peta(np,np,nlev)
    real (kind=real_kind) :: dp3d(np,np,nlev)

    real (kind=real_kind) :: zeta(np,np,nlev)
    real (kind=real_kind) :: pr3d(np,np,nlev+1)
    real (kind=real_kind) :: pr3d_ref(np,np,nlev+1)
    real (kind=real_kind) :: gp3d(np,np,nlev+1)

    real (kind=real_kind) :: ptop(np,np)
    real (kind=real_kind) :: sgp(np,np)
    real (kind=real_kind) :: tbar(nlev)

    ! note: non-state variables should be moved into derived_state_t

  end type elem_state_t

  !___________________________________________________________________
  type, public :: derived_state_t
     real (kind=real_kind) :: dummmy
     real (kind=real_kind) :: vstar(np,np,2,nlev)                     ! velocity on Lagrangian surfaces
  end type derived_state_t

  !___________________________________________________________________
  type, public :: elem_accum_t
    real (kind=real_kind) :: u(np,np,nlev)                            ! zonal velocity on sphere
    real (kind=real_kind) :: T(np,np,nlev)                            ! temperature
    real (kind=real_kind) :: ke(np,np,nlev)                           ! kinetic energy
  end type

#else
! ================== SHALLOW-WATER DATA-STRUCTURES ===================

  type, public :: elem_state_t

    ! prognostic variables for shallow-water solver
     real (kind=real_kind) :: p(np,np,nlev,timelevels)
     real (kind=real_kind) :: ps(np,np)                               ! surface geopotential
     real (kind=real_kind) :: gradps(np,np,2)                         ! gradient of surface geopotential     
     real (kind=real_kind) :: v(np,np,2,nlev,timelevels)              ! contravarient comp

#ifdef _SWDG
     ! variables for discontinuous-Galerkin shallow-water
     real (kind=real_kind) :: couv(np,np,2,nlev)                      ! covarient velocity
     real (kind=real_kind) :: uv(np,np,2,nlev)                        ! spherical (u,v) velocity  vector        
     real (kind=real_kind) :: psi(np,np,nlev)                         ! (ht+hs)*metdet(G)
     real (kind=real_kind) :: ht(np,np,nlev)                          ! height variable     
     real (kind=real_kind) :: hs(np,np,nlev)                          ! mountain height
#endif

  end type elem_state_t

  !___________________________________________________________________
  type, public :: derived_state_t
     real (kind=real_kind) :: vstar(np,np,2,nlev)                     ! velocity on Lagrangian surfaces
  end type derived_state_t

#endif

! ============= DATA-STRUCTURES COMMON TO ALL SOLVERS ================

  type, public :: index_t
     integer(kind=int_kind) :: ia(npsq),ja(npsq)
     integer(kind=int_kind) :: is,ie
     integer(kind=int_kind) :: NumUniquePts
     integer(kind=int_kind) :: UniquePtOffset
  end type index_t

  !___________________________________________________________________
  type, public :: element_t
     integer(kind=int_kind) :: LocalId
     integer(kind=int_kind) :: GlobalId

     ! Coordinate values of element points
     type (spherical_polar_t) :: spherep(np,np)                       ! Spherical coords of GLL points

     ! Equ-angular gnomonic projection coordinates
     type (cartesian2D_t)     :: cartp(np,np)                         ! gnomonic coords of GLL points 
     type (cartesian2D_t)     :: corners(4)                           ! gnomonic coords of element corners
     real (kind=real_kind)    :: u2qmap(4,2)                          ! bilinear map from ref element to quad in cubedsphere coordinates
                                                                      ! SHOULD BE REMOVED
     ! 3D cartesian coordinates
     type (cartesian3D_t)     :: corners3D(4)  

     ! Element diagnostics
     real (kind=real_kind)    :: area                                 ! Area of element
     real (kind=real_kind)    :: max_eig                              ! max singular value of metinv
     real (kind=real_kind)    :: min_eig                              ! min singular value of metinv
     real (kind=real_kind)    :: max_eig_ratio                        ! max ratio of singular values
     real (kind=real_kind)    :: dx_short                             ! short length scale
     real (kind=real_kind)    :: dx_long                              ! long length scale

     real (kind=real_kind)    :: variable_hyperviscosity(np,np)       ! hyperviscosity based on above
     real (kind=real_kind)    :: hv_courant                           ! hyperviscosity courant number
     real (kind=real_kind)    :: tensorVisc(2,2,np,np)                !og, matrix V for tensor viscosity

     ! Edge connectivity information
     integer(kind=int_kind)   :: node_numbers(4)
     integer(kind=int_kind)   :: node_multiplicity(4)                 ! number of elements sharing corner node

     type (GridVertex_t)      :: vertex                               ! element grid vertex information
     type (EdgeDescriptor_t)  :: desc

     type (elem_state_t)      :: state

     type (derived_state_t)   :: derived
#if defined _PRIM || defined _PRIMDG
     type (elem_accum_t)       :: accum
#endif
     ! Metric terms 
     real (kind=real_kind)    :: met(2,2,np,np)                       ! metric tensor on velocity and pressure grid
     real (kind=real_kind)    :: metinv(2,2,np,np)                    ! metric tensor on velocity and pressure grid
     real (kind=real_kind)    :: metdet(np,np)                        ! g = SQRT(det(g_ij)) on velocity and pressure grid
     real (kind=real_kind)    :: rmetdet(np,np)                       ! 1/metdet on velocity pressure grid
     real (kind=real_kind)    :: D(2,2,np,np)                         ! Map covariant field on cube to vector field on the sphere
     real (kind=real_kind)    :: Dinv(2,2,np,np)                      ! Map vector field on the sphere to covariant v on cube

     ! Convert vector fields from spherical to rectangular components
     ! The transpose of this operation is its pseudoinverse.
     real (kind=real_kind)    :: vec_sphere2cart(np,np,3,2)

     ! Mass matrix terms for an element on a cube face
     real (kind=real_kind)    :: mp(np,np)                            ! mass matrix on v and p grid
     real (kind=real_kind)    :: rmp(np,np)                           ! inverse mass matrix on v and p grid

     ! Mass matrix terms for an element on the sphere
     ! This mass matrix is used when solving the equations in weak form
     ! with the natural (surface area of the sphere) inner product
     real (kind=real_kind)    :: spheremp(np,np)                      ! mass matrix on v and p grid
     real (kind=real_kind)    :: rspheremp(np,np)                     ! inverse mass matrix on v and p grid

     integer(kind=long_kind)  :: gdofP(np,np)                         ! global degree of freedom (P-grid)

     real (kind=real_kind)    :: fcor(np,np)                          ! Coreolis term

     type (index_t) :: idxP
     type (index_t),pointer :: idxV
     integer :: FaceNum

     ! force element_t to be a multiple of 8 bytes.  
     ! on BGP, code will crash (signal 7, or signal 15) if 8 byte alignment is off
     ! check core file for:
     ! core.63:Generated by interrupt..(Alignment Exception DEAR=0xa1ef671c ESR=0x01800000 CCR0=0x4800a002)
     integer :: dummy
  end type element_t

  !___________________________________________________________________
  public :: element_coordinates
  public :: element_var_coordinates
  public :: element_var_coordinates3D
  public :: GetColumnIdP,GetColumnIdV
  public :: allocate_element_desc

contains

! ===================== ELEMENT_MOD METHODS ==========================

  function GetColumnIdP(elem,i,j) result(col_id)

    ! Get unique identifier for a Physics column on the P-grid

    type(element_t), intent(in) :: elem
    integer, intent(in) :: i,j
    integer :: col_id
    col_id = elem%gdofP(i,j)
  end function GetColumnIdP

  !___________________________________________________________________
  function GetColumnIdV(elem,i,j) result(col_id)

    !  Get unique identifier for a Physics column on the V-grid

    type(element_t), intent(in) :: elem
    integer, intent(in) :: i,j
    integer :: col_id
    col_id = elem%gdofP(i,j)
  end function GetColumnIdV

  !___________________________________________________________________
  function element_coordinates(start,end,points) result(cart)

    ! Initialize 2D rectilinear element colocation points

    use kinds, only : longdouble_kind
    type (cartesian2D_t), intent(in) :: start
    type (cartesian2D_t), intent(in) :: end
    real (kind=longdouble_kind), intent(in) :: points(:)

    type (cartesian2D_t) :: cart(SIZE(points),SIZE(points))
    type (cartesian2D_t) :: length, centroid
    real (kind=longdouble_kind) :: y 
    integer i,j

    length%x   = 0.50D0*(end%x-start%x)
    length%y   = 0.50D0*(end%y-start%y)
    centroid%x = 0.50D0*(end%x+start%x) 
    centroid%y = 0.50D0*(end%y+start%y) 
    do j=1,SIZE(points)
       y = centroid%y + length%y*points(j)
       do i=1,SIZE(points)
          cart(i,j)%x = centroid%x + length%x*points(i)
          cart(i,j)%y = y
       end do
    end do
  end function element_coordinates

  !___________________________________________________________________
  function element_var_coordinates(c,points) result(cart)

    use kinds, only : longdouble_kind
    type (cartesian2D_t), intent(in) :: c(4)
    real (kind=longdouble_kind), intent(in) :: points(:)
    type (cartesian2D_t) :: cart(SIZE(points),SIZE(points))

    real (kind=longdouble_kind) :: p(size(points))
    real (kind=longdouble_kind) :: q(size(points))
    integer i,j

    p(:) = (1.0D0-points(:))/2.0D0
    q(:) = (1.0D0+points(:))/2.0D0

    do j=1,SIZE(points)
       do i=1,SIZE(points)
          cart(i,j)%x = p(i)*p(j)*c(1)%x &
                      + q(i)*p(j)*c(2)%x &
                      + q(i)*q(j)*c(3)%x &
                      + p(i)*q(j)*c(4)%x 
          cart(i,j)%y = p(i)*p(j)*c(1)%y &
                      + q(i)*p(j)*c(2)%y &
                      + q(i)*q(j)*c(3)%y &
                      + p(i)*q(j)*c(4)%y 
       end do
    end do
  end function element_var_coordinates

  !___________________________________________________________________
  function element_var_coordinates3d(c,points) result(cart)

    use kinds, only : longdouble_kind
    type (cartesian3D_t), intent(in) :: c(4)
    real (kind=longdouble_kind), intent(in) :: points(:)
    type (cartesian3D_t) :: cart(SIZE(points),SIZE(points))

    real (kind=longdouble_kind) :: p(size(points))
    real (kind=longdouble_kind) :: q(size(points)),r
    integer i,j

    p(:) = (1.0D0-points(:))/2.0D0
    q(:) = (1.0D0+points(:))/2.0D0

    do j=1,SIZE(points)
       do i=1,SIZE(points)
          cart(i,j)%x = p(i)*p(j)*c(1)%x &
                      + q(i)*p(j)*c(2)%x &
                      + q(i)*q(j)*c(3)%x &
                      + p(i)*q(j)*c(4)%x 
          cart(i,j)%y = p(i)*p(j)*c(1)%y &
                      + q(i)*p(j)*c(2)%y &
                      + q(i)*q(j)*c(3)%y &
                      + p(i)*q(j)*c(4)%y 
          cart(i,j)%z = p(i)*p(j)*c(1)%z &
                      + q(i)*p(j)*c(2)%z &
                      + q(i)*q(j)*c(3)%z &
                      + p(i)*q(j)*c(4)%z 

          ! project back to sphere:
          r = distance(cart(i,j))
          cart(i,j)%x = cart(i,j)%x/r
          cart(i,j)%y = cart(i,j)%y/r
          cart(i,j)%z = cart(i,j)%z/r
       end do
    end do
  end function element_var_coordinates3d

  !___________________________________________________________________
  subroutine allocate_element_desc(elem)

    type (element_t), intent(inout)   :: elem(:)
    integer                           :: num, j,i

    num = SIZE(elem)
    
    do j=1,num
       allocate(elem(j)%desc%putmapP(max_neigh_edges))
       allocate(elem(j)%desc%getmapP(max_neigh_edges))
       allocate(elem(j)%desc%putmapP_ghost(max_neigh_edges))
       allocate(elem(j)%desc%getmapP_ghost(max_neigh_edges))
       allocate(elem(j)%desc%reverse(max_neigh_edges))
       allocate(elem(j)%desc%globalID(max_neigh_edges))
       allocate(elem(j)%desc%loc2buf(max_neigh_edges))
       do i=1,max_neigh_edges
          elem(j)%desc%loc2buf(i)=i
          elem(j)%desc%globalID(i)=-1
       enddo
       
    end do
  end subroutine allocate_element_desc


  !___________________________________________________________________
#if USE_OPENACC
  subroutine setup_element_pointers(elem)
    use dimensions_mod, only: nelemd
    implicit none
    type(element_t), intent(inout) :: elem(:)
    integer :: ie
    allocate( state_v                  (np,np,2,nlev,timelevels,nelemd)       )
    allocate( state_T                  (np,np,nlev,timelevels,nelemd)         )
    allocate( state_dp3d               (np,np,nlev,timelevels,nelemd)         )
    allocate( state_lnps               (np,np,timelevels,nelemd)              )
    allocate( state_ps_v               (np,np,timelevels,nelemd)              )
    allocate( state_phis               (np,np,nelemd)                         )
    allocate( state_Q                  (np,np,nlev,qsize_d,nelemd)            )
    allocate( state_Qdp                (np,np,nlev,qsize_d,2,nelemd)          )
    allocate( derived_vn0              (np,np,2,nlev,nelemd)                  )
    allocate( derived_vstar            (np,np,2,nlev,nelemd)                  )
    allocate( derived_dpdiss_biharmonic(np,np,nlev,nelemd)                    )
    allocate( derived_dpdiss_ave       (np,np,nlev,nelemd)                    )
    allocate( derived_phi              (np,np,nlev,nelemd)                    )
    allocate( derived_omega_p          (np,np,nlev,nelemd)                    )
    allocate( derived_eta_dot_dpdn     (np,np,nlevp,nelemd)                   )
    allocate( derived_grad_lnps        (np,np,2,nelemd)                       )
    allocate( derived_zeta             (np,np,nlev,nelemd)                    )
    allocate( derived_div              (np,np,nlev,timelevels,nelemd)         )
    allocate( derived_dp               (np,np,nlev,nelemd)                    )
    allocate( derived_divdp            (np,np,nlev,nelemd)                    )
    allocate( derived_divdp_proj       (np,np,nlev,nelemd)                    )
#ifdef CAM
    allocate( derived_FQ               (np,np,nlev,qsize_d,1,nelemd)          )
    allocate( derived_FM               (np,np,2,nlev,1,nelemd)                )
    allocate( derived_FT               (np,np,nlev,1,nelemd)                  )
    allocate( derived_omega_prescribed (np,np,nlev,nelemd)                    )
#else
    allocate( derived_FQ               (np,np,nlev,qsize_d,timelevels,nelemd) )
    allocate( derived_FM               (np,np,2,nlev,timelevels,nelemd)       )
    allocate( derived_FT               (np,np,nlev,timelevels,nelemd)         )
#endif
    allocate( derived_pecnd            (np,np,nlev,nelemd)                    )
    allocate( derived_FQps             (np,np,timelevels,nelemd)              )
    do ie = 1 , nelemd
      elem(ie)%state%v                   => state_v                  (:,:,:,:,:,ie)
      elem(ie)%state%T                   => state_T                  (:,:,:,:,ie)  
      elem(ie)%state%dp3d                => state_dp3d               (:,:,:,:,ie)  
      elem(ie)%state%lnps                => state_lnps               (:,:,:,ie)    
      elem(ie)%state%ps_v                => state_ps_v               (:,:,:,ie)    
      elem(ie)%state%phis                => state_phis               (:,:,ie)      
      elem(ie)%state%Q                   => state_Q                  (:,:,:,:,ie)  
      elem(ie)%state%Qdp                 => state_Qdp                (:,:,:,:,:,ie)
      elem(ie)%derived%vn0               => derived_vn0              (:,:,:,:,ie)  
      elem(ie)%derived%vstar             => derived_vstar            (:,:,:,:,ie)  
      elem(ie)%derived%dpdiss_biharmonic => derived_dpdiss_biharmonic(:,:,:,ie)    
      elem(ie)%derived%dpdiss_ave        => derived_dpdiss_ave       (:,:,:,ie)    
      elem(ie)%derived%phi               => derived_phi              (:,:,:,ie)    
      elem(ie)%derived%omega_p           => derived_omega_p          (:,:,:,ie)    
      elem(ie)%derived%eta_dot_dpdn      => derived_eta_dot_dpdn     (:,:,:,ie)    
      elem(ie)%derived%grad_lnps         => derived_grad_lnps        (:,:,:,ie)    
      elem(ie)%derived%zeta              => derived_zeta             (:,:,:,ie)    
      elem(ie)%derived%div               => derived_div              (:,:,:,:,ie)  
      elem(ie)%derived%dp                => derived_dp               (:,:,:,ie)    
      elem(ie)%derived%divdp             => derived_divdp            (:,:,:,ie)    
      elem(ie)%derived%divdp_proj        => derived_divdp_proj       (:,:,:,ie)    
#ifdef CAM                                                                        
      elem(ie)%derived%FQ                => derived_FQ               (:,:,:,:,:,ie)
      elem(ie)%derived%FM                => derived_FM               (:,:,:,:,:,ie)
      elem(ie)%derived%FT                => derived_FT               (:,:,:,:,ie)  
      elem(ie)%derived%omega_prescribed  => derived_omega_prescribed (:,:,:,ie)    
#else                                                                             
      elem(ie)%derived%FQ                => derived_FQ               (:,:,:,:,:,ie)
      elem(ie)%derived%FM                => derived_FM               (:,:,:,:,:,ie)
      elem(ie)%derived%FT                => derived_FT               (:,:,:,:,ie)  
#endif                                                                            
      elem(ie)%derived%pecnd             => derived_pecnd            (:,:,:,ie)    
      elem(ie)%derived%FQps              => derived_FQps             (:,:,:,ie)    
    enddo
  end subroutine setup_element_pointers
#endif

  


end module element_mod

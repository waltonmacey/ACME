#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_advance_mod
  !OVERWRITING: applycamforcing, applycamforcing_dynamics, prim_advance_init, prim_advance_exp, init_dp3d, prim_step_init, advance_hypervis_dp
  use prim_advance_mod_base, only: prim_advance_si, preq_robert3, smooth_phis, overwrite_SEdensity, advance_hypervis, advance_hypervis_lf
#if USE_OPENACC
  use openacc, only: acc_async_sync
  use dimensions_mod, only: np,nelemd
  use edgetype_mod, only : EdgeDescriptor_t, EdgeBuffer_t
  use kinds, only : real_kind, iulog
  use perf_mod, only: t_startf, t_stopf, t_barrierf, t_adj_detailf ! _EXTERNAL
  use parallel_mod, only : abortmp, parallel_t, iam
  use control_mod, only : se_prescribed_wind_2d
  implicit none
  private

  integer, parameter :: asyncid = acc_async_sync

  public :: prim_advance_si, preq_robert3, smooth_phis, overwrite_SEdensity

  public :: prim_advance_init
  public :: prim_advance_exp
  public :: init_dp3d, prim_step_init
  public :: copy_dynamics_d2h
  public :: copy_dynamics_h2d
  public :: applycamforcing
  public :: applycamforcing_dynamics

  real (kind=real_kind) :: initialized_for_dt   = 0
  real (kind=real_kind), allocatable :: ur_weights(:)
  type (EdgeBuffer_t) :: edge3p1


  real (kind=real_kind), allocatable :: grad_ps      (:,:,:,:)          ! lat-lon coord version
  real (kind=real_kind), allocatable :: grad_tmp1    (:,:,:)            !
  real (kind=real_kind), allocatable :: grad_tmp2    (:,:,:,:)          !
  real (kind=real_kind), allocatable :: dp           (:,:,:,:)          ! delta pressure
  real (kind=real_kind), allocatable :: rdp          (:,:,:,:)          ! inverse of delta pressure
  real (kind=real_kind), allocatable :: p            (:,:,:,:)          ! pressure
  real (kind=real_kind), allocatable :: grad_p       (:,:,:,:,:)        !
  real (kind=real_kind), allocatable :: vdp          (:,:,:,:,:)        !                            
  real (kind=real_kind), allocatable :: vgrad_p      (:,:,:,:)          ! v.grad(p)
  real (kind=real_kind), allocatable :: grad_p_m_pmet(:,:,:,:,:)        ! gradient(p - p_met)
  real (kind=real_kind), allocatable :: divdp        (:,:,:,:)          !
  real (kind=real_kind), allocatable :: vort         (:,:,:,:)          ! vorticity
  real (kind=real_kind), allocatable :: T_v          (:,:,:,:)          !
  real (kind=real_kind), allocatable :: kappa_star   (:,:,:,:)
  real (kind=real_kind), allocatable :: omega_p      (:,:,:,:)          !
  real (kind=real_kind), allocatable :: eta_dot_dpdn (:,:,:,:)          ! half level vertical velocity on p-grid
  real (kind=real_kind), allocatable :: T_vadv       (:,:,:,:)          ! temperature vertical advection
  real (kind=real_kind), allocatable :: v_vadv       (:,:,:,:,:)        ! velocity vertical advection
  real (kind=real_kind), allocatable :: sdot_sum     (:,:,:)            ! temporary field
  real (kind=real_kind), allocatable :: vtemp1       (:,:,:,:,:)        ! generic gradient storage
  real (kind=real_kind), allocatable :: vtemp2       (:,:,:,:,:)        ! generic gradient storage
  real (kind=real_kind), allocatable :: Ephi         (:,:,:,:)          ! kinetic energy + PHI term
  real (kind=real_kind), allocatable :: grads_tmp    (:,:,:,:,:)
  real (kind=real_kind), allocatable :: vtens        (:,:,:,:,:)
  real (kind=real_kind), allocatable :: ttens        (:,:,:,:)  
  real (kind=real_kind), allocatable :: dptens       (:,:,:,:)  
  real (kind=real_kind), allocatable :: div_tmp      (:,:,:,:)  
  real (kind=real_kind), allocatable :: vort_tmp     (:,:,:,:)  
  real (kind=real_kind), allocatable :: lap_t        (:,:,:,:)
  real (kind=real_kind), allocatable :: lap_dp       (:,:,:,:)
  real (kind=real_kind), allocatable :: lap_v        (:,:,:,:,:)



contains

  subroutine copy_dynamics_h2d( elem , nt )
    use element_mod, only: element_t, state_v, state_ps_v, state_t, state_qdp, state_dp3d, derived_vn0, derived_pecnd
    implicit none
    type(element_t), intent(inout) :: elem(:)
    integer        , intent(in   ) :: nt
    integer :: ie
    !$omp barrier
    !$omp master
    call t_startf('update device')
    !$acc update device(state_t,state_v,state_ps_v,derived_pecnd)
#if ( defined CAM )
    do ie = 1 , nelemd
      !$acc update device(elem(ie)%derived%u_met,elem(ie)%derived%v_met,elem(ie)%derived%dudt_met,elem(ie)%derived%dvdt_met,elem(ie)%derived%nudge_factor,elem(ie)%derived%Utnd,&
      !$acc&              elem(ie)%derived%Vtnd,elem(ie)%derived%t_met,elem(ie)%derived%dtdt_met,elem(ie)%derived%ttnd,elem(ie)%derived%ps_met,elem(ie)%derived%dpsdt_met)
    enddo
#endif
    call t_stopf('update device')
    !$omp end master
    !$omp barrier

  end subroutine copy_dynamics_h2d

  subroutine copy_dynamics_d2h( elem , nt )
    use element_mod, only: element_t, state_v, state_ps_v, state_t, state_dp3d, derived_vn0, derived_eta_dot_dpdn, derived_omega_p
    implicit none
    type(element_t), intent(inout) :: elem(:)
    integer        , intent(in   ) :: nt
    integer :: ie
    !$omp barrier
    !$omp master
    call t_startf('update host')
    do ie = 1 , nelemd
#if ( defined CAM )
      !$acc update host(elem(ie)%derived%utnd,elem(ie)%derived%vtnd,elem(ie)%derived%ttnd)
      !$acc update host(grad_p_m_pmet)
#endif
    enddo
    !$acc update host(state_t,state_v,state_ps_v,state_dp3d) 
    call t_stopf('update host')
    !$omp end master
    !$omp barrier
  end subroutine copy_dynamics_d2h

  subroutine prim_advance_init(par, elem,integration)
    use edge_mod, only : initEdgeBuffer
    use element_mod, only : element_t
    use dimensions_mod, only : nlev, nelemd
    use control_mod, only : qsplit,rsplit
    implicit none
    type (parallel_t) :: par
    type (element_t), intent(inout), target   :: elem(:)
    character(len=*)    , intent(in) :: integration
    integer :: i
    integer :: ie
    if (rsplit==0) then
      call initEdgeBuffer(par,edge3p1,elem,3*nlev+1)
    else
      ! need extra buffer space for dp3d
      call initEdgeBuffer(par,edge3p1,elem,4*nlev+1)
    endif
    !$acc enter data pcopyin(edge3p1         )
    !$acc enter data pcopyin(edge3p1%buf     )
    !$acc enter data pcopyin(edge3p1%receive )
    !$acc enter data pcopyin(edge3p1%putmap  )
    !$acc enter data pcopyin(edge3p1%getmap  )
    !$acc enter data pcopyin(edge3p1%reverse )
    !$acc enter data pcopyin(edge3p1%tag     )
    !$acc enter data pcopyin(edge3p1%srequest)
    !$acc enter data pcopyin(edge3p1%rrequest)

    ! compute averaging weights for RK+LF (tstep_type=1) timestepping:
    allocate(ur_weights(qsplit))
    ur_weights(:)=0.0d0
    if(mod(qsplit,2).NE.0)then
      ur_weights(1)=1.0d0/qsplit
      do i=3,qsplit,2
        ur_weights(i)=2.0d0/qsplit
      enddo
    else
      do i=2,qsplit,2
        ur_weights(i)=2.0d0/qsplit
      enddo
    endif

    allocate(grad_ps      (np,np,2       ,nelemd))
    allocate(grad_tmp1    (np,np         ,nelemd))
    allocate(grad_tmp2    (np,np,2       ,nelemd))
    allocate(dp           (np,np  ,nlev  ,nelemd))
    allocate(rdp          (np,np  ,nlev  ,nelemd))
    allocate(p            (np,np  ,nlev  ,nelemd))
    allocate(grad_p       (np,np,2,nlev  ,nelemd))
    allocate(vdp          (np,np,2,nlev  ,nelemd))
    allocate(vgrad_p      (np,np  ,nlev  ,nelemd))
    allocate(grad_p_m_pmet(np,np,2,nlev  ,nelemd))
    allocate(divdp        (np,np  ,nlev  ,nelemd))
    allocate(vort         (np,np  ,nlev  ,nelemd))
    allocate(t_v          (np,np  ,nlev  ,nelemd))
    allocate(kappa_star   (np,np  ,nlev  ,nelemd))
    allocate(omega_p      (np,np  ,nlev  ,nelemd))
    allocate(eta_dot_dpdn (np,np  ,nlev+1,nelemd))
    allocate(t_vadv       (np,np  ,nlev  ,nelemd))
    allocate(v_vadv       (np,np,2,nlev  ,nelemd))
    allocate(sdot_sum     (np,np         ,nelemd))
    allocate(vtemp1       (np,np,2,nlev  ,nelemd))
    allocate(vtemp2       (np,np,2,nlev  ,nelemd))
    allocate(ephi         (np,np  ,nlev  ,nelemd))
    allocate(grads_tmp    (np,np,2,nlev  ,nelemd))
    allocate(vtens        (np,np,2,nlev  ,nelemd))
    allocate(ttens        (np,np  ,nlev  ,nelemd))
    allocate(dptens       (np,np  ,nlev  ,nelemd))
    allocate(div_tmp      (np,np  ,nlev  ,nelemd))
    allocate(vort_tmp     (np,np  ,nlev  ,nelemd))
    allocate(lap_t        (np,np  ,nlev  ,nelemd))
    allocate(lap_dp       (np,np  ,nlev  ,nelemd))
    allocate(lap_v        (np,np,2,nlev  ,nelemd))
    !$acc enter data pcreate(grad_ps,dp,rdp,p,grad_p,vdp,vgrad_p,grad_tmp1,grad_tmp2,grad_p_m_pmet,divdp,vort,t_v,kappa_star,omega_p,eta_dot_dpdn,t_vadv,v_vadv,sdot_sum,vtemp1,vtemp2,ephi,&
    !$acc&                   grads_tmp,vtens,ttens,dptens,div_tmp,vort_tmp,lap_t,lap_dp,lap_v)
  end subroutine prim_advance_init


  subroutine init_dp3d(elem,hvcoord,nt,nets,nete)
    use dimensions_mod, only: np, nlev
    use element_mod   , only: element_t, state_ps_v, state_dp3d
    use hybvcoord_mod , only : hvcoord_t
    implicit none
    type(element_t), intent(inout) :: elem(:)
    type(hvcoord_t), intent(in   ) :: hvcoord
    integer        , intent(in   ) :: nt,nets,nete
    integer :: i,j,k,ie
    !$omp barrier
    !$omp master
    do ie = 1 , nelemd
      !$acc update device(state_ps_v(:,:,nt,ie))
    enddo
    !$acc parallel loop gang vector collapse(4) present(elem,hvcoord,state_ps_v,state_dp3d)
    do ie = 1 , nelemd
      do k = 1 , nlev
        do j = 1 , np
          do i = 1 , np
            state_dp3d(i,j,k,nt,ie) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                                      ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*state_ps_v(i,j,nt,ie)
          enddo
        enddo
      enddo
    enddo
    do ie = 1 , nelemd
      !$acc update host(state_dp3d(:,:,:,nt,ie))
    enddo
    !$omp end master
    !$omp barrier
  end subroutine init_dp3d


  subroutine prim_step_init(elem,hvcoord,nt,nets,nete)
    use dimensions_mod, only : np, nlev
    use element_mod   , only : element_t, state_ps_v, state_dp3d, derived_vn0, derived_eta_dot_dpdn, derived_omega_p
    use hybvcoord_mod , only : hvcoord_t
    use control_mod   , only : nu_p, rsplit
    implicit none
    type(element_t), intent(inout) :: elem(:)
    type(hvcoord_t), intent(in   ) :: hvcoord
    integer        , intent(in   ) :: nt,nets,nete
    integer :: i,j,k,ie
    !$omp barrier
    !$omp master
    !$acc parallel loop gang vector collapse(4) present(elem,derived_vn0,hvcoord,state_dp3d,state_ps_v,derived_eta_dot_dpdn,derived_omega_p)
    do ie = 1 , nelemd
      do k = 1 , nlev+1
        do j = 1 , np
          do i = 1 , np
            derived_eta_dot_dpdn(i,j,k,ie) = 0     ! mean vertical mass flux
            if (k <= nlev) derived_vn0    (i,j,:,k,ie) = 0     ! mean horizontal mass flux
            if (k <= nlev) derived_omega_p(i,j  ,k,ie) = 0
            if (nu_p>0) then
              if (k <= nlev) elem(ie)%derived%dpdiss_ave       (i,j,k)=0
              if (k <= nlev) elem(ie)%derived%dpdiss_biharmonic(i,j,k)=0
            endif
            if (rsplit==0) then
              ! save dp at time t for use in tracers
              if (k <= nlev) elem(ie)%derived%dp(i,j,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                                                          ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*state_ps_v(i,j,nt,ie)
            else
              ! dp at time t:  use floating lagrangian levels:
              if (k <= nlev) elem(ie)%derived%dp(i,j,k) = state_dp3d(i,j,k,nt,ie)
            endif
          enddo
        enddo
      enddo
    enddo
    !$omp end master
    !$omp barrier
  end subroutine prim_step_init


  subroutine prim_advance_exp(elem, deriv, hvcoord, hybrid, dt, tl, nets, nete, compute_diagnostics)
    use bndry_mod     , only : bndry_exchangev
    use control_mod   , only : prescribed_wind, qsplit, tstep_type, rsplit, qsplit, moisture, integration
    use derivative_mod, only : derivative_t, vorticity, divergence, gradient, gradient_wk
    use dimensions_mod, only : np, nlev, nlevp, nvar, nc, nelemd
    use edge_mod      , only : edgevpack, edgevunpack, initEdgeBuffer
    use edgetype_mod  , only : EdgeBuffer_t
    use element_mod   , only : element_t, state_v, state_t, state_ps_v, state_dp3d, derived_vn0
    use hybvcoord_mod , only : hvcoord_t
    use hybrid_mod    , only : hybrid_t
    use reduction_mod , only : reductionbuffer_ordered_1d_t
    use time_mod      , only : TimeLevel_t,  timelevel_qdp, tevolve
    use diffusion_mod , only : prim_diffusion

#ifndef CAM
    use asp_tests, only : asp_advection_vertical
#else
    use control_mod, only : prescribed_vertwind
#endif

    implicit none
    type (element_t)     , intent(inout), target :: elem(:)
    type (derivative_t)  , intent(in) :: deriv
    type (hvcoord_t)     , intent(in) :: hvcoord
    type (hybrid_t)      , intent(in) :: hybrid
    real (kind=real_kind), intent(in) :: dt
    type (TimeLevel_t)   , intent(in) :: tl
    integer              , intent(in) :: nets
    integer              , intent(in) :: nete
    logical              , intent(in) :: compute_diagnostics
    ! =================
    ! Local
    ! =================
    real (kind=real_kind) ::  dt2, time, dt_vis, x, eta_ave_w
    real (kind=real_kind) ::  eta_dot_dpdn(np,np,nlevp)
    real (kind=real_kind) ::  dp(np,np)
    real (kind=real_kind) ::  tempdp3d(np,np)
    real (kind=real_kind) ::  tempmass(nc,nc)
    real (kind=real_kind) ::  tempflux(nc,nc,4)
    real (kind=real_kind) ::  deta
    integer :: ie,nm1,n0,np1,nstep,method,qsplit_stage,k, qn0
    integer :: n,i,j,lx,lenx

    call t_adj_detailf(+1)
    call t_startf('prim_advance_exp')
    nm1   = tl%nm1
    n0    = tl%n0
    np1   = tl%np1
    nstep = tl%nstep

    ! timelevel to use for accessing Qdp() to compute virtual temperature
    qn0 = -1    ! -1 = disabled (assume dry dynamics)
    if ( moisture /= "dry") then
      call TimeLevel_Qdp(tl, qsplit, qn0)  ! compute current Qdp() timelevel
    endif

    ! integration = "explicit"
    !
    !   tstep_type=0  pure leapfrog except for very first timestep   CFL=1
    !                    typically requires qsplit=4 or 5
    !   tstep_type=1  RK2 followed by qsplit-1 leapfrog steps        CFL=close to qsplit
    !                    typically requires qsplit=4 or 5
    !   tstep_type=2  RK2-SSP 3 stage (as used by tracers)           CFL=.58
    !                    optimal in terms of SSP CFL, but not        CFLSSP=2
    !                    optimal in terms of CFL
    !                    typically requires qsplit=3
    !                    but if windspeed > 340m/s, could use this
    !                    with qsplit=1
    !   tstep_type=3  classic RK3                                    CFL=1.73 (sqrt(3))
    !
    !   tstep_type=4  Kinnmark&Gray RK4 4 stage                      CFL=sqrt(8)=2.8
    !                 should we replace by standard RK4 (CFL=sqrt(8))?
    !                 (K&G 1st order method has CFL=3)
    !   tstep_type=5  Kinnmark&Gray RK3 5 stage 3rd order            CFL=3.87  (sqrt(15))
    !                 From Paul Ullrich.  3rd order for nonlinear terms also
    !                 K&G method is only 3rd order for linear
    !                 optimal: for windspeeds ~120m/s,gravity: 340m/2
    !                 run with qsplit=1
    !                 (K&G 2nd order method has CFL=4. tiny CFL improvement not worth 2nd order)
    !
    ! integration = "full_imp"
    !
    !   tstep_type=1  Backward Euler or BDF2 implicit dynamics
    !
    ! default weights for computing mean dynamics fluxes
    eta_ave_w = 1d0/qsplit

    if(tstep_type==0)then
      method=0                ! pure leapfrog
      if (nstep==0) method=1  ! but use RK2 on first step
    else if (tstep_type==1) then
      method=0                           ! LF
      qsplit_stage = mod(nstep,qsplit)
      if (qsplit_stage==0) method=1      ! RK2 on first of qsplit steps
      ! RK2 + LF scheme has tricky weights:
      eta_ave_w=ur_weights(qsplit_stage+1)
    else
      method = tstep_type                ! other RK variants
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    ! fix dynamical variables, skip dynamics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    if (1==prescribed_wind .and. .not.se_prescribed_wind_2d) then
      time=tl%nstep*dt
      do ie=nets,nete
#ifdef CAM
        if (prescribed_vertwind==1) then
          eta_dot_dpdn(:,:,1) = 0.0d0
          eta_dot_dpdn(:,:,nlev+1) = 0.0d0
          do k = 2,nlev
            dp(:,:) = (hvcoord%hyam(k) - hvcoord%hyam(k-1))*hvcoord%ps0 + &
                      (hvcoord%hybm(k) - hvcoord%hybm(k-1))*elem(ie)%state%ps_v(:,:,tl%n0)
            deta = hvcoord%etam(k)-hvcoord%etam(k-1)
            eta_dot_dpdn(:,:,k) = dp(:,:)*elem(ie)%derived%etadot_prescribed(:,:,k)/deta
          enddo
        endif
#else
        ! assume most fields are constant in time
        elem(ie)%state%ps_v(:,:,np1) = elem(ie)%state%ps_v(:,:,n0)
        elem(ie)%state%lnps(:,:,np1) = elem(ie)%state%lnps(:,:,n0)
        elem(ie)%state%T(:,:,:,np1)   = elem(ie)%state%T(:,:,:,n0)
        elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,n0)
        elem(ie)%derived%div = 0
        !get eta_dot_dpdn at n0, ps_v is not used or calculated in this routine
        call asp_advection_vertical(time,hvcoord,elem(ie)%state%ps_v(:,:,n0),eta_dot_dpdn)
#endif
        ! accumulate mean fluxes for advection
        if (rsplit==0) then
          elem(ie)%derived%eta_dot_dpdn(:,:,:) = elem(ie)%derived%eta_dot_dpdn(:,:,:) + eta_dot_dpdn(:,:,:)*eta_ave_w
        else
          ! lagrangian case.  mean vertical velocity = 0
          ! compute dp3d on floating levels
          elem(ie)%derived%eta_dot_dpdn(:,:,:) = 0
          do k=1,nlev
            elem(ie)%state%dp3d(:,:,k,np1) = elem(ie)%state%dp3d(:,:,k,n0) + dt*(eta_dot_dpdn(:,:,k+1) - eta_dot_dpdn(:,:,k))
          enddo
        endif
        ! subcycling code uses a mean flux to advect tracers
        do k=1,nlev
          if (rsplit==0) then
            dp(:,:) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                      ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,tl%n0) 
          else
            dp(:,:) = elem(ie)%state%dp3d(:,:,k,tl%n0)
          endif
          elem(ie)%derived%vn0(:,:,1,k)=elem(ie)%derived%vn0(:,:,1,k)+eta_ave_w*elem(ie)%state%v(:,:,1,k,n0)*dp(:,:)
          elem(ie)%derived%vn0(:,:,2,k)=elem(ie)%derived%vn0(:,:,2,k)+eta_ave_w*elem(ie)%state%v(:,:,2,k,n0)*dp(:,:)
        enddo
      enddo
      call t_stopf('prim_advance_exp')
      call t_adj_detailf(-1)
      return
    endif

    ! ==================================
    ! Take timestep
    ! ==================================
    dt_vis = dt
    if (method==0) then
      ! regular LF step
      dt2 = 2*dt
      call compute_and_apply_rhs(np1,nm1,n0,qn0,dt2,elem,hvcoord,hybrid,deriv,nets,nete,compute_diagnostics,eta_ave_w)
      dt_vis = dt2  ! dt to use for time-split dissipation
    else if (method==1) then
      ! RK2
      ! forward euler to u(dt/2) = u(0) + (dt/2) RHS(0)  (store in u(np1))
      call compute_and_apply_rhs(np1,n0,n0 ,qn0,dt/2,elem,hvcoord,hybrid,deriv,nets,nete,compute_diagnostics,0d0     )
      ! leapfrog:  u(dt) = u(0) + dt RHS(dt/2)     (store in u(np1))
      call compute_and_apply_rhs(np1,n0,np1,qn0,dt  ,elem,hvcoord,hybrid,deriv,nets,nete,.false.            ,eta_ave_w)
    else if (method==2) then
      ! RK2-SSP 3 stage.  matches tracer scheme. optimal SSP CFL, but
      ! not optimal for regular CFL
      ! u1 = u0 + dt/2 RHS(u0)
      call compute_and_apply_rhs(np1,n0 ,n0 ,qn0,dt/2,elem,hvcoord,hybrid,deriv,nets,nete,compute_diagnostics,eta_ave_w/3)
      ! u2 = u1 + dt/2 RHS(u1)
      call compute_and_apply_rhs(np1,np1,np1,qn0,dt/2,elem,hvcoord,hybrid,deriv,nets,nete,.false.            ,eta_ave_w/3)
      ! u3 = u2 + dt/2 RHS(u2)
      call compute_and_apply_rhs(np1,np1,np1,qn0,dt/2,elem,hvcoord,hybrid,deriv,nets,nete,.false.            ,eta_ave_w/3)
      ! unew = u/3 +2*u3/3  = u + 1/3 (RHS(u) + RHS(u1) + RHS(u2))
      do ie=nets,nete
        elem(ie)%state%v(:,:,:,:,np1)    = elem(ie)%state%v(:,:,:,:,n0) /3 + 2*elem(ie)%state%v(:,:,:,:,np1) /3
        elem(ie)%state%T(:,:,:,np1)      = elem(ie)%state%T(:,:,:,n0)   /3 + 2*elem(ie)%state%T(:,:,:,np1)   /3
        if (rsplit==0) then
          elem(ie)%state%ps_v(:,:,np1)   = elem(ie)%state%ps_v(:,:,n0)  /3 + 2*elem(ie)%state%ps_v(:,:,np1)  /3
        else
          elem(ie)%state%dp3d(:,:,:,np1) = elem(ie)%state%dp3d(:,:,:,n0)/3 + 2*elem(ie)%state%dp3d(:,:,:,np1)/3
        endif
      enddo
    else if (method==3) then
      ! classic RK3  CFL=sqrt(3)
      ! u1 = u0 + dt/3 RHS(u0)
      call compute_and_apply_rhs(np1,n0,n0 ,qn0,dt/3,elem,hvcoord,hybrid,deriv,nets,nete,compute_diagnostics,0d0      )
      ! u2 = u0 + dt/2 RHS(u1)
      call compute_and_apply_rhs(np1,n0,np1,qn0,dt/2,elem,hvcoord,hybrid,deriv,nets,nete,.false.            ,0d0      )
      ! u3 = u0 + dt RHS(u2)
      call compute_and_apply_rhs(np1,n0,np1,qn0,dt  ,elem,hvcoord,hybrid,deriv,nets,nete,.false.            ,eta_ave_w)
    else if (method==4) then
      ! KG 4th order 4 stage:   CFL=sqrt(8)
      ! low storage version of classic RK4
      ! u1 = u0 + dt/4 RHS(u0)
      call compute_and_apply_rhs(np1,n0,n0 ,qn0,dt/4,elem,hvcoord,hybrid,deriv,nets,nete,compute_diagnostics,0d0      )
      ! u2 = u0 + dt/3 RHS(u1)
      call compute_and_apply_rhs(np1,n0,np1,qn0,dt/3,elem,hvcoord,hybrid,deriv,nets,nete,.false.            ,0d0      )
      ! u3 = u0 + dt/2 RHS(u2)
      call compute_and_apply_rhs(np1,n0,np1,qn0,dt/2,elem,hvcoord,hybrid,deriv,nets,nete,.false.            ,0d0      )
      ! u4 = u0 + dt RHS(u3)
      call compute_and_apply_rhs(np1,n0,np1,qn0,dt  ,elem,hvcoord,hybrid,deriv,nets,nete,.false.            ,eta_ave_w)
    else if (method==5) then
      ! Ullrich 3nd order 5 stage:   CFL=sqrt( 4^2 -1) = 3.87
      ! u1 = u0 + dt/5 RHS(u0)  (save u1 in timelevel nm1)
      ! phl: rhs: t=t
      call compute_and_apply_rhs(nm1,n0 ,n0 ,qn0,  dt/5,elem,hvcoord,hybrid,deriv,nets,nete,compute_diagnostics,  eta_ave_w/4)
      ! u2 = u0 + dt/5 RHS(u1)
      ! phl: rhs: t=t+dt/5
      call compute_and_apply_rhs(np1,n0 ,nm1,qn0,  dt/5,elem,hvcoord,hybrid,deriv,nets,nete,.false.            ,0d0          )
      ! u3 = u0 + dt/3 RHS(u2)
      ! phl: rhs: t=t+2*dt/5
      call compute_and_apply_rhs(np1,n0 ,np1,qn0,  dt/3,elem,hvcoord,hybrid,deriv,nets,nete,.false.            ,0d0          )
      ! u4 = u0 + 2dt/3 RHS(u3)
      ! phl: rhs: t=t+2*dt/5+dt/3
      call compute_and_apply_rhs(np1,n0 ,np1,qn0,2*dt/3,elem,hvcoord,hybrid,deriv,nets,nete,.false.            ,0d0          )
      ! compute (5*u1/4 - u0/4) in timelevel nm1:
      !$omp barrier
      !$omp master
      !$acc parallel loop gang vector collapse(4) present(state_v,state_t,state_dp3d)
      do ie=1,nelemd
        do k = 1 , nlev
          do j = 1 , np
            do i = 1 ,np
              state_v(i,j,:,k,nm1,ie) = (5*state_v(i,j,:,k,nm1,ie) - state_v(i,j,:,k,n0,ie) ) / 4
              state_T(i,j  ,k,nm1,ie) = (5*state_T(i,j  ,k,nm1,ie) - state_T(i,j  ,k,n0,ie) ) / 4
              if (rsplit/=0) state_dp3d(i,j,k,nm1,ie) = (5*state_dp3d(i,j,k,nm1,ie) - state_dp3d(i,j,k,n0,ie) ) / 4
            enddo
          enddo
        enddo
      enddo
      if (rsplit==0) then
        !$acc parallel loop gang vector collapse(3) present(state_ps_v)
        do ie = 1 , nelemd
          do j = 1 , np
            do i = 1 , np
              state_ps_v(i,j,nm1,ie)   = (5*state_ps_v(i,j,nm1,ie)   - state_ps_v(i,j,n0,ie)   ) / 4
            enddo
          enddo
        enddo
      endif
      !$omp end master
      !$omp barrier
      ! u5 = (5*u1/4 - u0/4) + 3dt/4 RHS(u4)
      ! phl: rhs: t=t+2*dt/5+dt/3+3*dt/4         -wrong RK times ...
      call compute_and_apply_rhs(np1,nm1,np1,qn0,3*dt/4,elem,hvcoord,hybrid,deriv,nets,nete,.false.            ,3*eta_ave_w/4)
      ! final method is the same as:
      ! u5 = u0 +  dt/4 RHS(u0)) + 3dt/4 RHS(u4)
    else if ((method==11).or.(method==12)) then
      call abortmp('ERROR: tstep_type == 11 or 12 not supported for OpenACC')
    else
      call abortmp('ERROR: bad choice of tstep_type')
    endif

    ! ==============================================
    ! Time-split Horizontal diffusion: nu.del^2 or nu.del^4
    ! U(*) = U(t+1)  + dt2 * HYPER_DIFF_TERM(t+1)
    ! ==============================================
#ifdef ENERGY_DIAGNOSTICS
    if (compute_diagnostics) then
       do ie = nets,nete
          elem(ie)%accum%DIFF(:,:,:,:)=elem(ie)%state%v(:,:,:,:,np1)
          elem(ie)%accum%DIFFT(:,:,:)=elem(ie)%state%T(:,:,:,np1)
       enddo
    endif
#endif

    ! note:time step computes u(t+1)= u(t*) + RHS.
    ! for consistency, dt_vis = t-1 - t*, so this is timestep method dependent
    if (tstep_type==0) then
      ! leapfrog special case
      call advance_hypervis_lf(edge3p1,elem,hvcoord,hybrid,deriv,nm1,n0,np1,nets,nete,dt_vis)
    else if (method<=10) then ! not implicit
      if (rsplit==0) then
        ! forward-in-time, maybe hypervis applied to PS
        call advance_hypervis(edge3p1,elem,hvcoord,hybrid,deriv,np1,nets,nete,dt_vis,eta_ave_w)
      else
        ! forward-in-time, hypervis applied to dp3d
        call advance_hypervis_dp(edge3p1,elem,hvcoord,hybrid,deriv,np1,nets,nete,dt_vis,eta_ave_w)
      endif
    endif

#ifdef ENERGY_DIAGNOSTICS
    if (compute_diagnostics) then
      do ie = nets,nete
        do k=1,nlev  !  Loop index added (AAM)
          elem(ie)%accum%DIFF(:,:,:,k)=( elem(ie)%state%v(:,:,:,k,np1) - elem(ie)%accum%DIFF(:,:,:,k) ) / dt_vis
          elem(ie)%accum%DIFFT(:,:,k) =( elem(ie)%state%T(:,:,k,np1)   - elem(ie)%accum%DIFFT(:,:,k)  ) / dt_vis
        enddo
      enddo
    endif
#endif

    tevolve=tevolve+dt

    call t_stopf('prim_advance_exp')
    call t_adj_detailf(-1)
  end subroutine prim_advance_exp



  ! phl notes: output is stored in first argument. Advances from 2nd argument using tendencies evaluated at 3rd rgument: 
  ! phl: for offline winds use time at 3rd argument (same as rhs currently)
  subroutine compute_and_apply_rhs(np1,nm1,n0,qn0,dt2,elem,hvcoord,hybrid,deriv,nets,nete,compute_diagnostics,eta_ave_w)
    ! ===================================
    ! compute the RHS, accumulate into u(np1) and apply DSS
    !
    !           u(np1) = u(nm1) + dt2*DSS[ RHS(u(n0)) ]
    !
    ! This subroutine is normally called to compute a leapfrog timestep
    ! but by adjusting np1,nm1,n0 and dt2, many other timesteps can be
    ! accomodated.  For example, setting nm1=np1=n0 this routine will
    ! take a forward euler step, overwriting the input with the output.
    !
    !    qn0 = timelevel used to access Qdp() in order to compute virtual Temperature
    !          qn0=-1 for the dry case
    !
    ! if  dt2<0, then the DSS'd RHS is returned in timelevel np1
    !
    ! Combining the RHS and DSS pack operation in one routine
    ! allows us to fuse these two loops for more cache reuse
    !
    ! Combining the dt advance and DSS unpack operation in one routine
    ! allows us to fuse these two loops for more cache reuse
    !
    ! note: for prescribed velocity case, velocity will be computed at
    ! "real_time", which should be the time of timelevel n0.
    ! ===================================
    use kinds                 , only : real_kind
    use dimensions_mod        , only : np, nc, nlev, max_corner_elem
    use hybrid_mod            , only : hybrid_t
    use element_mod           , only : element_t,PrintElem, state_ps_v, timelevels, derived_vn0, state_v, state_t, state_dp3d, derived_eta_dot_dpdn, derived_omega_p, derived_pecnd, state_qdp
    use derivative_mod        , only : derivative_t, divergence_sphere, gradient_sphere, vorticity_sphere
    use derivative_mod        , only : subcell_div_fluxes, subcell_dss_fluxes
    use edge_mod              , only : edgevpack_openacc, edgevunpack_openacc
    use edgetype_mod          , only : edgedescriptor_t
    use bndry_mod             , only : bndry_exchangev => bndry_exchangev_simple_minimize_pcie
    use control_mod           , only : moisture, qsplit, use_cpstar, rsplit, swest
    use hybvcoord_mod         , only : hvcoord_t
    use physical_constants    , only : cp, cpwater_vapor, Rgas, kappa
    use physics_mod           , only : virtual_specific_heat, virtual_temperature1d
    use prim_si_mod           , only : preq_hydrostatic_openacc, preq_omega_ps_openacc, preq_vertadv_openacc
    use derivative_mod        , only : gradient_sphere_openacc, divergence_sphere_openacc, vorticity_sphere_openacc
    use openacc_utils_mod     , only : memset
#if ( defined CAM )     
    use control_mod       , only: se_met_nudge_u, se_met_nudge_p, se_met_nudge_t, se_met_tevolve
#endif
    use time_mod          , only : tevolve
    implicit none
    integer               , intent(in)    :: np1,nm1,n0,qn0,nets,nete
    real*8                , intent(in)    :: dt2
    logical               , intent(in)    :: compute_diagnostics
    type (hvcoord_t)      , intent(in)    :: hvcoord
    type (hybrid_t)       , intent(in)    :: hybrid
    type (element_t)      , intent(inout) :: elem(:)
    type (derivative_t)   , intent(in)    :: deriv
    real (kind=real_kind) , intent(in)    :: eta_ave_w  ! weighting for eta_dot_dpdn mean flux
    ! local

    real(kind=real_kind) :: vtens1,vtens2,ttens,vgrad_t


    type (EdgeDescriptor_t) :: desc
    real (kind=real_kind) ::  cp2,cp_ratio,E,de,Qt,v1,v2,fcor_tmp,vort_tmp
    real (kind=real_kind) ::  glnps1,glnps2,gpterm
    integer :: i,j,k,kptr,ie, ii
    real (kind=real_kind) ::  u_m_umet, v_m_vmet, t_m_tmet 
    real (kind=real_kind) :: tmp1(np,np), tmp2(np,np,nlev)
    !$acc routine(virtual_specific_heat) seq
    !$acc routine(virtual_temperature1d) seq

    call t_adj_detailf(+1)
    call t_startf('compute_and_apply_rhs')
    !$omp barrier
    !$omp master
    ! ==================================================
    ! compute pressure (p) on half levels from ps
    ! using the hybrid coordinates relationship, i.e.
    ! e.g. equation (3.a.92) of the CCM-2 description,
    ! (NCAR/TN-382+STR), June 1993, p. 24.
    ! ==================================================
    ! vertically eulerian only needs grad(ps)
    if (rsplit == 0) call gradient_sphere_openacc(state_ps_v,deriv,elem,grad_ps,1,1,nelemd,timelevels,n0,1,1,asyncid_in=asyncid)
      
    !$acc parallel loop gang vector collapse(4) present(dp,p,grad_p,elem,hvcoord,grad_ps,rdp,state_ps_v,state_dp3d) async(asyncid)
    do ie = 1 , nelemd
      ! ============================
      ! compute p and delta p
      ! ============================
      do k = 1 , nlev
        do j = 1 , np
          do i = 1 , np
            if (rsplit == 0) then
              dp(i,j,k,ie) = (hvcoord%hyai(k+1)*hvcoord%ps0 + hvcoord%hybi(k+1)*state_ps_v(i,j,n0,ie)) &
                           - (hvcoord%hyai(k  )*hvcoord%ps0 + hvcoord%hybi(k  )*state_ps_v(i,j,n0,ie))
              rdp(i,j,k,ie) = 1.D0 / dp(i,j,k,ie)
              p(i,j,k,ie)   = hvcoord%hyam(k  )*hvcoord%ps0 + hvcoord%hybm(k  )*state_ps_v(i,j,n0,ie)
              grad_p(i,j,1:2,k,ie) = hvcoord%hybm(k)*grad_ps(i,j,1:2,ie)
            else
              ! vertically lagrangian code: we advect dp3d instead of ps_v
              ! we also need grad(p) at all levels (not just grad(ps))
              !p(k)= hyam(k)*ps0 + hybm(k)*ps
              !    = .5*(hyai(k+1)+hyai(k))*ps0 + .5*(hybi(k+1)+hybi(k))*ps
              !    = .5*(ph(k+1) + ph(k) )  = ph(k) + dp(k)/2
              !
              ! p(k+1)-p(k) = ph(k+1)-ph(k) + (dp(k+1)-dp(k))/2
              !             = dp(k) + (dp(k+1)-dp(k))/2 = (dp(k+1)+dp(k))/2
              dp(i,j,k,ie) = state_dp3d(i,j,k,n0,ie)
              rdp(i,j,k,ie) = 1.0D0/dp(i,j,k,ie)
            endif
          enddo
        enddo
      enddo
    enddo
    if (rsplit /= 0) then
      !$acc parallel loop gang vector collapse(3) present(dp,hvcoord,p) async(asyncid)
      do ie = 1 , nelemd 
        do j = 1 , np
          do i = 1 , np
            p(i,j,1,ie)=hvcoord%hyai(1)*hvcoord%ps0 + dp(i,j,1,ie)/2
            !$acc loop seq
            do k = 2 , nlev
              p(i,j,k,ie) = p(i,j,k-1,ie) + (dp(i,j,k-1,ie) + dp(i,j,k,ie))/2
            enddo
          enddo
        enddo
      enddo
      call gradient_sphere_openacc(p,deriv,elem,grad_p,nlev,1,nelemd,1,1,1,1,asyncid_in=asyncid)
    endif
#if ( defined CAM )
    if (se_met_nudge_p > 0.D0) then
      !$acc parallel loop gang vector collapse(3) present(grad_tmp1,elem) async(asyncid)
      do ie = 1 , nelemd
        do j = 1 , np
          do i = 1 , np
            grad_tmp1(i,j,ie) = elem(ie)%derived%ps_met(i,j)+tevolve*elem(ie)%derived%dpsdt_met(i,j)
          enddo
        enddo
      enddo
      call gradient_sphere_openacc(grad_tmp1,deriv,elem,grad_tmp2,1,1,nelemd,1,1,1,1,asyncid_in=asyncid)
      !$acc parallel loop gang vector collapse(5) present(grad_tmp2,hvcoord,grad_p,grad_p_m_pmet) async(asyncid)
      do ie = 1 , nelemd
        do k = 1 , nlev
          do ii = 1 , 2
            do j = 1 , np
              do i = 1 , np
                grad_p_m_pmet(i,j,ii,k,ie) = grad_p(i,j,ii,k,ie) - hvcoord%hybm(k) * grad_tmp2(i,j,ii,ie) 
              enddo
            enddo
          enddo
        enddo
      enddo
    endif
#endif
    !$acc parallel loop gang vector collapse(4) present(vgrad_p,vdp,dp,grad_p,state_v,derived_vn0,t_v,state_t,kappa_star,elem,state_qdp) private(v1,v2,qt) async(asyncid)
    do ie = 1 , nelemd
      do k = 1 , nlev
        ! ============================
        ! compute vgrad_lnps
        ! ============================
        do j = 1  ,np
          do i = 1 , np
            v1 = state_v(i,j,1,k,n0,ie)
            v2 = state_v(i,j,2,k,n0,ie)
            vgrad_p(i,j,k,ie) = (v1*grad_p(i,j,1,k,ie) + v2*grad_p(i,j,2,k,ie))
            vdp(i,j,1,k,ie) = v1*dp(i,j,k,ie)
            vdp(i,j,2,k,ie) = v2*dp(i,j,k,ie)
            ! ================================
            ! Accumulate mean Vel_rho flux in vn0
            ! ================================
            derived_vn0(i,j,:,k,ie) = derived_vn0(i,j,:,k,ie) + eta_ave_w * vdp(i,j,:,k,ie)
            if (qn0 == -1) then
              T_v(i,j,k,ie) = state_T(i,j,k,n0,ie)
              kappa_star(i,j,k,ie) = kappa
            else
              Qt = state_Qdp(i,j,k,1,qn0,ie)/dp(i,j,k,ie)
              T_v(i,j,k,ie) = Virtual_Temperature1d(state_T(i,j,k,n0,ie),Qt)
              if (use_cpstar==1) then
                kappa_star(i,j,k,ie) =  Rgas/Virtual_Specific_Heat(Qt)
              else
                kappa_star(i,j,k,ie) = kappa
              endif
            endif
          enddo
        enddo
      enddo
    enddo
    ! =========================================
    ! Compute relative vorticity and divergence
    ! =========================================
    call divergence_sphere_openacc(vdp    ,deriv,elem,divdp,nlev,1,nelemd,1         ,1 ,1,1,asyncid_in=asyncid)
    call vorticity_sphere_openacc (state_v,deriv,elem,vort ,nlev,1,nelemd,timelevels,n0,1,1,asyncid_in=asyncid)
    ! ====================================================
    ! Compute Hydrostatic equation, modeld after CCM-3
    ! ====================================================
    call preq_hydrostatic_openacc(elem,T_v,p,dp,1,nelemd,asyncid_in=asyncid)
    ! ====================================================
    ! Compute omega_p according to CCM-3
    ! ====================================================
    call preq_omega_ps_openacc(omega_p,p,vgrad_p,divdp,1,nelemd,asyncid_in=asyncid)
    ! ==================================================
    ! Compute eta_dot_dpdn
    ! save sdot_sum as this is the -RHS of ps_v equation
    ! ==================================================
    if (rsplit > 0) then
      call memset(product(shape(eta_dot_dpdn)),eta_dot_dpdn,0._real_kind,asyncid_in=asyncid)
      call memset(product(shape(t_vadv      )),t_vadv      ,0._real_kind,asyncid_in=asyncid)
      call memset(product(shape(v_vadv      )),v_vadv      ,0._real_kind,asyncid_in=asyncid)
      call memset(product(shape(sdot_sum    )),sdot_sum    ,0._real_kind,asyncid_in=asyncid)
    else
      !NOT SUPPORTED YET! do ie = 1 , nelemd
      !NOT SUPPORTED YET!   do k = 1 , nlev
      !NOT SUPPORTED YET!     do j = 1 , np
      !NOT SUPPORTED YET!       do i = 1 , np
      !NOT SUPPORTED YET!         tmp2(i,j,k) = divdp(i,j,k,ie)
      !NOT SUPPORTED YET!       enddo
      !NOT SUPPORTED YET!     enddo
      !NOT SUPPORTED YET!   enddo
      !NOT SUPPORTED YET!   do j = 1 , np
      !NOT SUPPORTED YET!     do i = 1 , np
      !NOT SUPPORTED YET!       tmp1(i,j) = 0.
      !NOT SUPPORTED YET!       do k = 1 , nlev
      !NOT SUPPORTED YET!         ! ==================================================
      !NOT SUPPORTED YET!         ! add this term to PS equation so we exactly conserve dry mass
      !NOT SUPPORTED YET!         ! ==================================================
      !NOT SUPPORTED YET!         tmp1(i,j) = tmp1(i,j) + tmp2(i,j,k)
      !NOT SUPPORTED YET!       enddo
      !NOT SUPPORTED YET!       sdot_sum(i,j,ie) = tmp1(i,j)
      !NOT SUPPORTED YET!       eta_dot_dpdn(i,j,     1,ie) = 0.D0
      !NOT SUPPORTED YET!       eta_dot_dpdn(i,j,nlev+1,ie) = 0.D0
      !NOT SUPPORTED YET!     enddo
      !NOT SUPPORTED YET!   enddo
      !NOT SUPPORTED YET!   ! ===========================================================
      !NOT SUPPORTED YET!   ! at this point, eta_dot_dpdn contains integral_etatop^eta[ divdp ]
      !NOT SUPPORTED YET!   ! compute at interfaces:
      !NOT SUPPORTED YET!   !    eta_dot_dpdn = -dp/dt - integral_etatop^eta[ divdp ]
      !NOT SUPPORTED YET!   ! for reference: at mid layers we have:
      !NOT SUPPORTED YET!   !    omega = v grad p  - integral_etatop^eta[ divdp ]
      !NOT SUPPORTED YET!   ! ===========================================================
      !NOT SUPPORTED YET!   do k = 1 , nlev-1
      !NOT SUPPORTED YET!     do j = 1 , np
      !NOT SUPPORTED YET!       do i = 1 , np
      !NOT SUPPORTED YET!         eta_dot_dpdn(i,j,k+1,ie) = ( hvcoord%hybi(k+1) - 1 ) * tmp1(i,j)
      !NOT SUPPORTED YET!       enddo
      !NOT SUPPORTED YET!     enddo
      !NOT SUPPORTED YET!   enddo
      !NOT SUPPORTED YET! enddo
      !NOT SUPPORTED YET! ! ===========================================================
      !NOT SUPPORTED YET! ! Compute vertical advection of T and v from eq. CCM2 (3.b.1)
      !NOT SUPPORTED YET! ! ==============================================
      !NOT SUPPORTED YET! call preq_vertadv(elem, eta_dot_dpdn,rdp,T_vadv,v_vadv,1,nelemd,n0)
    endif
    ! ================================
    ! accumulate mean vertical flux:
    ! ================================
    !$acc parallel loop gang vector collapse(4) present(elem,eta_dot_dpdn,omega_p,ephi,state_v,derived_eta_dot_dpdn) async(asyncid)
    do ie = 1 , nelemd
      do k = 1 , nlev+1  !  Loop index added (AAM)
        do j = 1 , np
          do i = 1 , np
            derived_eta_dot_dpdn(i,j,k,ie) = derived_eta_dot_dpdn(i,j,k,ie) + eta_ave_w*eta_dot_dpdn(i,j,k,ie)
            if (k <= nlev) Ephi(i,j,k,ie) = 0.5D0*( state_v(i,j,1,k,n0,ie)**2 + state_v(i,j,2,k,n0,ie)**2 ) + elem(ie)%derived%phi(i,j,k) + derived_pecnd(i,j,k,ie)
          enddo
        enddo
      enddo
    enddo
    call gradient_sphere_openacc(state_t,deriv,elem,vtemp1,nlev,1,nelemd,timelevels,n0,1,1,asyncid_in=asyncid)
    call gradient_sphere_openacc(ephi,deriv,elem,vtemp2,nlev,1,nelemd,1,1,1,1,asyncid_in=asyncid)
    !$acc parallel loop gang vector collapse(4) present(t_v,p,grad_p,state_v,v_vadv,t_vadv,kappa_star,vort,elem,vtemp1,vtemp2,omega_p,state_t,state_v,state_dp3d,divdp,eta_dot_dpdn,derived_omega_p) &
    !$acc&                                      private(gpterm,glnps1,glnps2,v1,v2,u_m_umet,v_m_vmet,t_m_tmet,vtens1,vtens2,ttens,vgrad_t) async(asyncid)
    do ie=1,nelemd
      do k=1,nlev
        do j=1,np
          do i=1,np
             derived_omega_p(i,j,k,ie) = derived_omega_p(i,j,k,ie) + eta_ave_w*omega_p(i,j,k,ie)
             vgrad_T = state_v(i,j,1,k,n0,ie) * vtemp1(i,j,1,k,ie) + state_v(i,j,2,k,n0,ie) * vtemp1(i,j,2,k,ie)
             gpterm = T_v(i,j,k,ie)/p(i,j,k,ie)
             glnps1 = Rgas*gpterm*grad_p(i,j,1,k,ie)
             glnps2 = Rgas*gpterm*grad_p(i,j,2,k,ie)
             v1     = state_v(i,j,1,k,n0,ie)
             v2     = state_v(i,j,2,k,n0,ie)
             fcor_tmp = elem(ie)%fcor(i,j)
             vort_tmp = vort(i,j,k,ie)
             vtens1 = - v_vadv(i,j,1,k,ie) + v2*(fcor_tmp + vort_tmp) - vtemp2(i,j,1,k,ie) - glnps1
             vtens2 = - v_vadv(i,j,2,k,ie) - v1*(fcor_tmp + vort_tmp) - vtemp2(i,j,2,k,ie) - glnps2
             ttens = - T_vadv(i,j  ,k,ie) - vgrad_T + kappa_star(i,j,k,ie)*T_v(i,j,k,ie)*omega_p(i,j,k,ie)
#if ( defined CAM )
             if (se_prescribed_wind_2d) then
               vtens1 = 0.D0
               vtens2 = 0.D0
               ttens = 0.D0
             else
               if(se_met_nudge_u.gt.0.D0)then
                 u_m_umet = v1 - elem(ie)%derived%u_met(i,j,k) - se_met_tevolve*tevolve*elem(ie)%derived%dudt_met(i,j,k)
                 v_m_vmet = v2 - elem(ie)%derived%v_met(i,j,k) - se_met_tevolve*tevolve*elem(ie)%derived%dvdt_met(i,j,k)
                 vtens1 =   vtens1 - se_met_nudge_u*u_m_umet * elem(ie)%derived%nudge_factor(i,j,k)
                 elem(ie)%derived%Utnd(i+(j-1)*np,k) = elem(ie)%derived%Utnd(i+(j-1)*np,k) + se_met_nudge_u*u_m_umet * elem(ie)%derived%nudge_factor(i,j,k)
                 vtens2 =   vtens2 - se_met_nudge_u*v_m_vmet * elem(ie)%derived%nudge_factor(i,j,k)
                 elem(ie)%derived%Vtnd(i+(j-1)*np,k) = elem(ie)%derived%Vtnd(i+(j-1)*np,k) + se_met_nudge_u*v_m_vmet * elem(ie)%derived%nudge_factor(i,j,k)
               endif
               if(se_met_nudge_p.gt.0.D0)then
                 vtens1 =   vtens1 - se_met_nudge_p*grad_p_m_pmet(i,j,1,k,ie)  * elem(ie)%derived%nudge_factor(i,j,k)
                 vtens2 =   vtens2 - se_met_nudge_p*grad_p_m_pmet(i,j,2,k,ie)  * elem(ie)%derived%nudge_factor(i,j,k)
               endif
               if(se_met_nudge_t.gt.0.D0)then
                 t_m_tmet = state_T(i,j,k,n0,ie) -  elem(ie)%derived%T_met(i,j,k) - se_met_tevolve*tevolve*elem(ie)%derived%dTdt_met(i,j,k)
                 ttens  = ttens - se_met_nudge_t*t_m_tmet * elem(ie)%derived%nudge_factor(i,j,k)
                 elem(ie)%derived%Ttnd(i+(j-1)*np,k) = elem(ie)%derived%Ttnd(i+(j-1)*np,k) + se_met_nudge_t*t_m_tmet * elem(ie)%derived%nudge_factor(i,j,k)
               endif
             endif
#endif
    ! =========================================================
    ! local element timestep, store in np1.
    ! note that we allow np1=n0 or nm1
    ! apply mass matrix
    ! =========================================================
            if (dt2<0) then
              ! calling program just wanted DSS'd RHS, skip time advance
              state_v(i,j,1,k,np1,ie) = elem(ie)%spheremp(i,j)*vtens1
              state_v(i,j,2,k,np1,ie) = elem(ie)%spheremp(i,j)*vtens2
              state_T(i,j  ,k,np1,ie) = elem(ie)%spheremp(i,j)*ttens
              if (rsplit>0) state_dp3d(i,j,k,np1,ie) = -elem(ie)%spheremp(i,j)*(divdp(i,j,k,ie) + eta_dot_dpdn(i,j,k+1,ie)-eta_dot_dpdn(i,j,k,ie))
            else
              state_v(i,j,1,k,np1,ie) = elem(ie)%spheremp(i,j)*( state_v(i,j,1,k,nm1,ie) + dt2*vtens1 )
              state_v(i,j,2,k,np1,ie) = elem(ie)%spheremp(i,j)*( state_v(i,j,2,k,nm1,ie) + dt2*vtens2 )
              state_T(i,j  ,k,np1,ie) = elem(ie)%spheremp(i,j)*( state_T(i,j  ,k,nm1,ie) + dt2*ttens )
              if (rsplit>0) state_dp3d(i,j,k,np1,ie) = elem(ie)%spheremp(i,j) * (state_dp3d(i,j,k,nm1,ie) - dt2 * (divdp(i,j,k,ie) + eta_dot_dpdn(i,j,k+1,ie)-eta_dot_dpdn(i,j,k,ie)))
            endif
          enddo
        enddo
      enddo
    enddo
    !$acc parallel loop gang vector collapse(3) present(elem,sdot_sum,state_ps_v) async(asyncid)
    do ie=1,nelemd
      do j = 1 , np
        do i = 1 , np
          if (dt2<0) then
            state_ps_v(i,j,np1,ie) = -elem(ie)%spheremp(i,j)*sdot_sum(i,j,ie)
          else
            state_ps_v(i,j,np1,ie) = elem(ie)%spheremp(i,j)*( state_ps_v(i,j,nm1,ie) - dt2*sdot_sum(i,j,ie) )
          endif
        enddo
      enddo
    enddo
    call t_startf('c_a_rhs_pack')
    call edgeVpack_openacc(edge3p1,state_ps_v,  1   ,  0     ,elem,1   ,nelemd,timelevels,np1)
    call edgeVpack_openacc(edge3p1,state_t   ,  nlev,  1     ,elem,1   ,nelemd,timelevels,np1)
    call edgeVpack_openacc(edge3p1,state_v   ,2*nlev,  nlev+1,elem,1   ,nelemd,timelevels,np1)
    if (rsplit>0) call edgeVpack_openacc(edge3p1,state_dp3d,  nlev,3*nlev+1,elem,1   ,nelemd,timelevels,np1)
    call t_stopf('c_a_rhs_pack')
    !$omp end master
    !$omp barrier

    call t_startf('c_a_rhs_bndry')
    call bndry_exchangeV(hybrid,edge3p1)
    call t_stopf('c_a_rhs_bndry')

    !$omp barrier
    !$omp master
    call t_startf('c_a_rhs_unpack')
    call edgeVunpack_openacc(edge3p1,state_ps_v,  1   ,  0     ,elem,1   ,nelemd,timelevels,np1)
    call edgeVunpack_openacc(edge3p1,state_t   ,  nlev,  1     ,elem,1   ,nelemd,timelevels,np1)
    call edgeVunpack_openacc(edge3p1,state_v   ,2*nlev,  nlev+1,elem,1   ,nelemd,timelevels,np1)
    if (rsplit>0) call edgeVunpack_openacc(edge3p1,state_dp3d,  nlev,3*nlev+1,elem,1   ,nelemd,timelevels,np1)
    call t_stopf('c_a_rhs_unpack')

    !$acc parallel loop gang vector collapse(4) present(state_t,state_v,state_dp3d,elem)
    do ie=1,nelemd
      ! ====================================================
      ! Scale tendencies by inverse mass matrix
      ! ====================================================
      do k=1,nlev
        do j = 1 , np
          do i = 1 , np
            state_T(i,j  ,k,np1,ie) = elem(ie)%rspheremp(i,j)*state_T(i,j  ,k,np1,ie)
            state_v(i,j,1,k,np1,ie) = elem(ie)%rspheremp(i,j)*state_v(i,j,1,k,np1,ie)
            state_v(i,j,2,k,np1,ie) = elem(ie)%rspheremp(i,j)*state_v(i,j,2,k,np1,ie)
            if (rsplit > 0) state_dp3d(i,j,k,np1,ie)= elem(ie)%rspheremp(i,j)*state_dp3d(i,j,k,np1,ie)
          enddo
        enddo
      enddo
    enddo
    if (rsplit <= 0) then
      !$acc parallel loop gang vector collapse(3) present(state_ps_v,elem)
      do ie = 1 , nelemd
        ! vertically eulerian: complete ps_v timestep:
        do j = 1 , np
          do i = 1 , np
            state_ps_v(i,j,np1,ie) = elem(ie)%rspheremp(i,j)*state_ps_v(i,j,np1,ie)
          enddo
        enddo
      enddo
    endif

    !$omp end master
    !$omp barrier
    call t_stopf('compute_and_apply_rhs')
    call t_adj_detailf(-1)
  end subroutine compute_and_apply_rhs



  subroutine applyCAMforcing(elem,fvm,hvcoord,np1,np1_qdp,dt_q,nets,nete)
    use dimensions_mod        , only : np, nc, nlev, qsize
    use element_mod           , only : element_t
    use hybvcoord_mod         , only : hvcoord_t
    use control_mod           , only : moisture, tracer_grid_type
    use control_mod           , only : TRACER_GRIDTYPE_GLL, TRACER_GRIDTYPE_FVM
    use physical_constants    , only : Cp
    use fvm_control_volume_mod, only : fvm_struct, n0_fvm
    implicit none
    type (element_t)     , intent(inout) :: elem(:)
    type(fvm_struct)     , intent(inout) :: fvm(:)
    real (kind=real_kind), intent(in)    :: dt_q
    type (hvcoord_t)     , intent(in)    :: hvcoord
    integer              , intent(in)    :: np1,nets,nete,np1_qdp
    ! local
    integer :: i,j,k,ie,q
    real (kind=real_kind) :: v1,dp
    real (kind=real_kind) :: beta(np,np),E0(np,np),ED(np,np),dp0m1(np,np),dpsum(np,np)
    logical :: wet
    call t_startf('applyCAMforcing_openacc')
    wet = (moisture /= "dry")
    do ie=nets,nete
      ! apply forcing to Qdp
      elem(ie)%derived%FQps(:,:,1)=0
      ! even when running fvm tracers we need to updates forcing on ps and qv on GLL grid
      ! for fvm tracer qsize is usually 1 (qv)
      do q=1,qsize
        do k=1,nlev
          do j=1,np
            do i=1,np
              v1 = dt_q*elem(ie)%derived%FQ(i,j,k,q,1)
              if (elem(ie)%state%Qdp(i,j,k,q,np1_qdp) + v1 < 0 .and. v1<0) then
                if (elem(ie)%state%Qdp(i,j,k,q,np1_qdp) < 0 ) then
                  v1=0  ! Q already negative, dont make it more so
                else
                  v1 = -elem(ie)%state%Qdp(i,j,k,q,np1_qdp)
                endif
              endif
              elem(ie)%state%Qdp(i,j,k,q,np1_qdp) = elem(ie)%state%Qdp(i,j,k,q,np1_qdp)+v1
              if (q==1) then
                elem(ie)%derived%FQps(i,j,1)=elem(ie)%derived%FQps(i,j,1)+v1/dt_q
              endif
            enddo
          enddo
        enddo
      enddo
      if (wet .and. qsize>0) then
        ! to conserve dry mass in the precese of Q1 forcing:
        elem(ie)%state%ps_v(:,:,np1) = elem(ie)%state%ps_v(:,:,np1) + dt_q*elem(ie)%derived%FQps(:,:,1)
      endif
      ! Qdp(np1) and ps_v(np1) were updated by forcing - update Q(np1)
      do q=1,qsize
        do k=1,nlev
          do j=1,np
            do i=1,np
              dp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                   ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,np1)
              elem(ie)%state%Q(i,j,k,q) = elem(ie)%state%Qdp(i,j,k,q,np1_qdp)/dp
            enddo
          enddo
        enddo
      enddo
      elem(ie)%state%T(:,:,:,np1)   = elem(ie)%state%T(:,:,:,np1)   + dt_q*elem(ie)%derived%FT(:,:,:,1)
      elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,np1) + dt_q*elem(ie)%derived%FM(:,:,:,:,1)
    enddo
    call t_stopf('applyCAMforcing_openacc')
  end subroutine applyCAMforcing



  subroutine applyCAMforcing_dynamics(elem,hvcoord,np1,dt_q,nets,nete)
    use dimensions_mod, only : np, nlev, qsize
    use element_mod   , only : element_t
    use hybvcoord_mod , only : hvcoord_t
    implicit none
    type (element_t)     , intent(inout) :: elem(:)
    real (kind=real_kind), intent(in)    :: dt_q
    type (hvcoord_t)     , intent(in)    :: hvcoord
    integer              , intent(in)    :: np1,nets,nete
    ! local
    integer :: i,j,k,ie,q
    real (kind=real_kind) :: v1,dp
    logical :: wet
    call t_startf('applyCAMforcing_dynamics_openacc')
    do ie=nets,nete
      elem(ie)%state%T(:,:,:,np1)   = elem(ie)%state%T(:,:,:,np1)   + dt_q*elem(ie)%derived%FT(:,:,:,1)
      elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,np1) + dt_q*elem(ie)%derived%FM(:,:,:,:,1)
    enddo
    call t_stopf('applyCAMforcing_dynamics_openacc')
  end subroutine applyCAMforcing_dynamics



  subroutine advance_hypervis_dp(edge3,elem,hvcoord,hybrid,deriv,nt,nets,nete,dt2,eta_ave_w)
    !  take one timestep of:
    !          u(:,:,:,np) = u(:,:,:,np) +  dt2*nu*laplacian**order ( u )
    !          T(:,:,:,np) = T(:,:,:,np) +  dt2*nu_s*laplacian**order ( T )
    !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
    use dimensions_mod       , only : np, np, nlev, nc, max_corner_elem
    use control_mod          , only : nu, nu_div, nu_s, hypervis_order, hypervis_subcycle, nu_p, nu_top, psurf_vis, swest
    use hybrid_mod           , only : hybrid_t
    use hybvcoord_mod        , only : hvcoord_t
    use element_mod          , only : element_t, state_t, state_v, state_dp3d, timelevels
    use derivative_mod       , only : derivative_t
    use derivative_mod       , only : subcell_Laplace_fluxes, subcell_dss_fluxes
    use edge_mod             , only : edgevpack_openacc, edgevunpack_openacc
    use edgetype_mod         , only : EdgeBuffer_t, EdgeDescriptor_t
    use bndry_mod            , only : bndry_exchangev => bndry_exchangeV_simple_minimize_pcie
    use physical_constants   , only : Cp
    use derivative_mod       , only : subcell_Laplace_fluxes
    use viscosity_mod        , only : biharmonic_wk_dp3d_openacc
    use derivative_mod       , only : laplace_sphere_wk_openacc, vlaplace_sphere_wk_openacc
    implicit none
    type (hybrid_t)      , intent(in)    :: hybrid
    type (element_t)     , intent(inout) :: elem(:)
    type (EdgeBuffer_t)  , intent(inout) :: edge3
    type (derivative_t)  , intent(in)    :: deriv
    type (hvcoord_t)     , intent(in)    :: hvcoord
    real (kind=real_kind), intent(in)    :: dt2
    integer              , intent(in)    :: nets,nete
    ! local
    real (kind=real_kind) :: eta_ave_w  ! weighting for mean flux terms
    real (kind=real_kind) :: dpdn,dpdn0, nu_scale_top
    integer :: k,kptr,i,j,ie,ic,nt
    real (kind=real_kind), dimension(np,np,nlev) :: p
    real (kind=real_kind), dimension(np,np) :: dptemp1,dptemp2
    type (EdgeDescriptor_t)                                       :: desc
    real (kind=real_kind) :: v1,v2,dt,heating,utens_tmp,vtens_tmp,ttens_tmp,dptens_tmp
    real (kind=real_kind)                     :: temp      (np,np,nlev)
    real (kind=real_kind)                     :: laplace_fluxes(nc,nc,4)
    if (nu_s == 0 .and. nu == 0 .and. nu_p==0 ) return;

    call t_adj_detailf(+1) 
    call t_startf('advance_hypervis_dp_oacc')

    dt=dt2/hypervis_subcycle
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  regular viscosity
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! NOT YET SUPPORTED ! if (hypervis_order == 1) then
    ! NOT YET SUPPORTED !   if (nu_p>0) call abortmp( 'ERROR: hypervis_order == 1 not coded for nu_p>0')
    ! NOT YET SUPPORTED !   do ic=1,hypervis_subcycle
    ! NOT YET SUPPORTED !     do ie=nets,nete
    ! NOT YET SUPPORTED !       do k=1,nlev
    ! NOT YET SUPPORTED !         lap_t=laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
    ! NOT YET SUPPORTED !         lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)
    ! NOT YET SUPPORTED !         ! advace in time.  (note: DSS commutes with time stepping, so we
    ! NOT YET SUPPORTED !         ! can time advance and then DSS.  this has the advantage of
    ! NOT YET SUPPORTED !         ! not letting any discontinuties accumulate in p,v via roundoff
    ! NOT YET SUPPORTED !         do j=1,np
    ! NOT YET SUPPORTED !           do i=1,np
    ! NOT YET SUPPORTED !             elem(ie)%state%T(i,j,k  ,nt)=elem(ie)%state%T(i,j,k  ,nt)*elem(ie)%spheremp(i,j) + dt*nu_s*lap_t(i,j)
    ! NOT YET SUPPORTED !             elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%state%v(i,j,1,k,nt)*elem(ie)%spheremp(i,j) + dt*nu  *lap_v(i,j,1)
    ! NOT YET SUPPORTED !             elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%state%v(i,j,2,k,nt)*elem(ie)%spheremp(i,j) + dt*nu  *lap_v(i,j,2)
    ! NOT YET SUPPORTED !           enddo
    ! NOT YET SUPPORTED !         enddo
    ! NOT YET SUPPORTED !       enddo
    ! NOT YET SUPPORTED !       kptr=0
    ! NOT YET SUPPORTED !       call edgeVpack(edge3, elem(ie)%state%T(:,:,:,nt),nlev,kptr,ie)
    ! NOT YET SUPPORTED !       kptr=nlev
    ! NOT YET SUPPORTED !       call edgeVpack(edge3,elem(ie)%state%v(:,:,:,:,nt),2*nlev,kptr,ie)
    ! NOT YET SUPPORTED !     enddo
    ! NOT YET SUPPORTED !     call bndry_exchangeV(hybrid,edge3)
    ! NOT YET SUPPORTED !     do ie=nets,nete
    ! NOT YET SUPPORTED !       kptr=0
    ! NOT YET SUPPORTED !       call edgeVunpack(edge3, elem(ie)%state%T(:,:,:,nt), nlev, kptr, ie)
    ! NOT YET SUPPORTED !       kptr=nlev
    ! NOT YET SUPPORTED !       call edgeVunpack(edge3, elem(ie)%state%v(:,:,:,:,nt), 2*nlev, kptr, ie)
    ! NOT YET SUPPORTED !       ! apply inverse mass matrix
    ! NOT YET SUPPORTED !       do k=1,nlev
    ! NOT YET SUPPORTED !         do j=1,np
    ! NOT YET SUPPORTED !           do i=1,np
    ! NOT YET SUPPORTED !             elem(ie)%state%T(i,j,k  ,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%T(i,j,k  ,nt)
    ! NOT YET SUPPORTED !             elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,1,k,nt)
    ! NOT YET SUPPORTED !             elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,2,k,nt)
    ! NOT YET SUPPORTED !           enddo
    ! NOT YET SUPPORTED !         enddo
    ! NOT YET SUPPORTED !       enddo
    ! NOT YET SUPPORTED !     enddo
    ! NOT YET SUPPORTED !     !$OMP BARRIER
    ! NOT YET SUPPORTED !   enddo  ! subcycle
    ! NOT YET SUPPORTED ! endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  hyper viscosity
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! nu_p=0:
    !   scale T dissipaton by dp  (conserve IE, dissipate T^2)
    ! nu_p>0
    !   dont scale:  T equation IE dissipation matches (to truncation error)
    !                IE dissipation from continuity equation
    !                (1 deg: to about 0.1 W/m^2)
    if (hypervis_order == 2) then
      do ic=1,hypervis_subcycle
        call biharmonic_wk_dp3d_openacc(elem,dptens,ttens,vtens,div_tmp,vort_tmp,grads_tmp,deriv,edge3,hybrid,nt,1,nelemd)
        !$omp barrier
        !$omp master
        if (nu_p>0) then
          !$acc parallel loop gang vector collapse(4) present(elem,dptens,state_dp3d)
          do ie=1,nelemd
            do k = 1 , nlev
              do j = 1 , np
                do i = 1 , np
                  ! comptue mean flux
                  elem(ie)%derived%dpdiss_ave       (i,j,k)=elem(ie)%derived%dpdiss_ave       (i,j,k)+eta_ave_w*state_dp3d(i,j,k,nt,ie)/hypervis_subcycle
                  elem(ie)%derived%dpdiss_biharmonic(i,j,k)=elem(ie)%derived%dpdiss_biharmonic(i,j,k)+eta_ave_w*dptens    (i,j,k   ,ie)/hypervis_subcycle
                enddo
              enddo
            enddo
          enddo
        endif
        if (nu_top>0) then
          call  laplace_sphere_wk_openacc(state_t   ,grads_tmp       ,deriv,elem,.false.,lap_t ,nlev,1,nelemd,timelevels,nt,1,1)
          call  laplace_sphere_wk_openacc(state_dp3d,grads_tmp       ,deriv,elem,.false.,lap_dp,nlev,1,nelemd,timelevels,nt,1,1)
          call vlaplace_sphere_wk_openacc(state_v   ,div_tmp,vort_tmp,deriv,elem,.false.,lap_v ,nlev,1,nelemd,timelevels,nt,1,1)
        endif
        !$acc parallel loop gang vector collapse(4) present(lap_v,lap_t,lap_dp,vtens,ttens,dptens,elem,state_dp3d) private(utens_tmp,vtens_tmp,ttens_tmp,dptens_tmp,nu_scale_top)
        do ie=1,nelemd
          do k=1,nlev
            do j=1,np
              do i=1,np
                nu_scale_top = 1
                if (k==1) nu_scale_top=4
                if (k==2) nu_scale_top=2
                ! biharmonic terms need a negative sign:
                if (nu_top>0 .and. k<=3) then
                  utens_tmp =(-nu  *vtens (i,j,1,k,ie) + nu_scale_top*nu_top*lap_v (i,j,1,k,ie))
                  vtens_tmp =(-nu  *vtens (i,j,2,k,ie) + nu_scale_top*nu_top*lap_v (i,j,2,k,ie))
                  ttens_tmp =(-nu_s*ttens (i,j  ,k,ie) + nu_scale_top*nu_top*lap_t (i,j  ,k,ie))
                  dptens_tmp=(-nu_p*dptens(i,j  ,k,ie) + nu_scale_top*nu_top*lap_dp(i,j  ,k,ie))
                else
                  utens_tmp =-nu  *vtens (i,j,1,k,ie)
                  vtens_tmp =-nu  *vtens (i,j,2,k,ie)
                  ttens_tmp =-nu_s*ttens (i,j  ,k,ie)
                  dptens_tmp=-nu_p*dptens(i,j  ,k,ie)
                endif
                ttens (i,j  ,k,ie) = ttens_tmp
                dptens(i,j  ,k,ie) = dptens_tmp
                vtens (i,j,1,k,ie) = utens_tmp
                vtens (i,j,2,k,ie) = vtens_tmp
                ! NOTE: we will DSS all tendicies, EXCEPT for dp3d, where we DSS the new state
                state_dp3d(i,j,k,nt,ie) = state_dp3d(i,j,k,nt,ie)*elem(ie)%spheremp(i,j) + dt*dptens(i,j,k,ie)
              enddo
            enddo
          enddo
        enddo
        call t_startf('a_h_dp_pack')
        call edgeVpack_openacc(edge3, ttens    ,  nlev,     0,elem,1,nelemd,         1, 1)
        call edgeVpack_openacc(edge3, vtens    ,2*nlev,  nlev,elem,1,nelemd,         1, 1)
        call edgeVpack_openacc(edge3,state_dp3d,  nlev,3*nlev,elem,1,nelemd,timelevels,nt)
        call t_stopf('a_h_dp_pack')
        !$omp end master
        !$omp barrier
        call t_startf('a_h_dp_bndry')
        call bndry_exchangeV(hybrid,edge3)
        call t_stopf('a_h_dp_bndry')
        !$omp barrier
        !$omp master
        call t_startf('a_h_dp_unpack')
        call edgeVunpack_openacc(edge3, ttens    ,  nlev,     0,elem,1,nelemd,         1, 1)
        call edgeVunpack_openacc(edge3, vtens    ,2*nlev,  nlev,elem,1,nelemd,         1, 1)
        call edgeVunpack_openacc(edge3,state_dp3d,  nlev,3*nlev,elem,1,nelemd,timelevels,nt)
        call t_stopf('a_h_dp_unpack')
        !$acc parallel loop gang vector collapse(4) present(elem,vtens,ttens,state_v,state_t,state_dp3d) private(heating)
        do ie = 1,nelemd
          ! apply inverse mass matrix, accumulate tendencies
          do k = 1 , nlev
            do j = 1 , np
              do i = 1 , np
                state_dp3d(i,j  ,k,nt,ie) = state_dp3d(i,j  ,k,nt,ie) * elem(ie)%rspheremp(i,j)
                ! update v first (gives better results than updating v after heating)
                state_v   (i,j,:,k,nt,ie) = state_v   (i,j,:,k,nt,ie) + dt*vtens(i,j,:,k,ie) * elem(ie)%rspheremp(i,j)
                heating = dt*( vtens(i,j,1,k,ie)*state_v(i,j,1,k,nt,ie)  + vtens(i,j,2,k,ie)*state_v(i,j,2,k,nt,ie) ) * elem(ie)%rspheremp(i,j)
                state_T   (i,j  ,k,nt,ie) = state_T   (i,j  ,k,nt,ie) + dt*ttens(i,j  ,k,ie) * elem(ie)%rspheremp(i,j) - heating/cp
              enddo
            enddo
          enddo
        enddo
        !$omp end master
        !$omp barrier
      enddo
    endif

    call t_stopf('advance_hypervis_dp_oacc')
    call t_adj_detailf(-1)

  end subroutine advance_hypervis_dp



#endif
end module prim_advance_mod


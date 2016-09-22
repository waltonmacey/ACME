#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_advance_mod
  use prim_advance_mod_base, only: prim_advance_si, preq_robert3,applyCAMforcing_dynamics, applyCAMforcing, smooth_phis, overwrite_SEdensity

  use edgetype_mod  , only: EdgeDescriptor_t, EdgeBuffer_t
  use kinds         , only: real_kind, iulog
  use perf_mod      , only: t_startf, t_stopf, t_barrierf, t_adj_detailf ! _EXTERNAL
  use parallel_mod  , only: abortmp, parallel_t, iam
  use control_mod   , only: se_prescribed_wind_2d
  use dimensions_mod, only: np, nlev, nelemd
  use element_mod    , only: timelevels
  implicit none
  private
  save
  public :: prim_advance_exp, prim_advance_si, prim_advance_init, preq_robert3,&
            applyCAMforcing_dynamics, applyCAMforcing, smooth_phis, overwrite_SEdensity

  type (EdgeBuffer_t) :: edge1
  type (EdgeBuffer_t) :: edge2
  type (EdgeBuffer_t) :: edge3p1
  real (kind=real_kind) :: initialized_for_dt   = 0
  real (kind=real_kind), allocatable :: ur_weights(:)



  real (kind=real_kind),allocatable :: p (:,:,:,:)         ! pressure
  real (kind=real_kind),allocatable :: dp(:,:,:,:)         ! delta pressure
  real (kind=real_kind),allocatable :: grad_p(:,:,:,:,:)
  real (kind=real_kind),allocatable :: vgrad_p(:,:,:,:)    ! v.grad(p)
  real (kind=real_kind),allocatable :: vdp(:,:,:,:,:)       !                            
  real (kind=real_kind),allocatable :: grad_p_m_pmet(:,:,:,:,:)  ! gradient(p - p_met)
  real (kind=real_kind),allocatable :: divdp(:,:,:,:)
  real (kind=real_kind),allocatable :: vort(:,:,:,:)       ! vorticity
  real (kind=real_kind),allocatable :: T_v(:,:,:,:)
  real (kind=real_kind),allocatable :: kappa_star(:,:,:,:)
  real (kind=real_kind),allocatable :: omega_p(:,:,:,:)
  real (kind=real_kind),allocatable :: eta_dot_dpdn(:,:,:,:)  ! half level vertical velocity on p-grid
  real (kind=real_kind),allocatable :: sdot_sum(:,:,:)   ! temporary field
  real (kind=real_kind),allocatable :: T_vadv(:,:,:,:)     ! temperature vertical advection
  real (kind=real_kind),allocatable :: v_vadv(:,:,:,:,:)   ! velocity vertical advection
  real (kind=real_kind),allocatable :: Ephi(:,:,:,:)       ! kinetic energy + PHI term
  real (kind=real_kind),allocatable :: vtemp1(:,:,:,:,:)     ! generic gradient storage
  real (kind=real_kind),allocatable :: vtemp2(:,:,:,:,:)     ! generic gradient storage

contains



  subroutine prim_advance_init(par, elem,integration)
    use edge_mod, only : initEdgeBuffer
    use element_mod, only : element_t
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
    allocate(p (np,np,nlev,nelemd))
    allocate(dp(np,np,nlev,nelemd))
    allocate(grad_p(np,np,2,nlev,nelemd))
    allocate(vgrad_p(np,np,nlev,nelemd))
    allocate(vdp(np,np,2,nlev,nelemd))
    allocate(grad_p_m_pmet(np,np,2,nlev,nelemd))
    allocate(divdp(np,np,nlev,nelemd))
    allocate(vort(np,np,nlev,nelemd))
    allocate(t_v(np,np,nlev,nelemd))
    allocate(kappa_star(np,np,nlev,nelemd))
    allocate(omega_p(np,np,nlev,nelemd))
    allocate(eta_dot_dpdn(np,np,nlev+1,nelemd))
    allocate(sdot_sum(np,np,nelemd))
    allocate(T_vadv(np,np,nlev,nelemd))
    allocate(v_vadv(np,np,2,nlev,nelemd))
    allocate(ephi(np,np,nlev,nelemd))
    allocate(vtemp1(np,np,2,nlev,nelemd))
    allocate(vtemp2(np,np,2,nlev,nelemd))
    !$acc enter data pcreate(p,dp,grad_p,vgrad_p,vdp,grad_p_m_pmet,divdp,vort,t_v,kappa_star,omega_p,eta_dot_dpdn,sdot_sum,t_vadv,v_vadv,ephi,vtemp1,vtemp2)
    !$acc enter data pcopyin(edge3p1         )
    !$acc enter data pcopyin(edge3p1%buf     )
    !$acc enter data pcopyin(edge3p1%receive )
    !$acc enter data pcopyin(edge3p1%putmap  )
    !$acc enter data pcopyin(edge3p1%getmap  )
    !$acc enter data pcopyin(edge3p1%reverse )
    !$acc enter data pcopyin(edge3p1%tag     )
    !$acc enter data pcopyin(edge3p1%srequest)
    !$acc enter data pcopyin(edge3p1%rrequest)
  end subroutine prim_advance_init



  subroutine prim_advance_exp(elem, deriv, hvcoord, hybrid, dt, tl,  nets, nete, compute_diagnostics)
    use bndry_mod, only : bndry_exchangev
    use control_mod, only : prescribed_wind, qsplit, tstep_type, rsplit, qsplit, moisture, integration
    use derivative_mod, only : derivative_t, vorticity, divergence, gradient, gradient_wk
    use dimensions_mod, only : np, nlev, nlevp, nvar, nc, nelemd
    use edge_mod, only : edgevpack, edgevunpack, initEdgeBuffer
    use edgetype_mod, only : EdgeBuffer_t
    use element_mod, only : element_t, state_v, state_t, state_dp3d
    use hybvcoord_mod, only : hvcoord_t
    use hybrid_mod, only : hybrid_t
    use reduction_mod, only : reductionbuffer_ordered_1d_t
    use time_mod, only : TimeLevel_t,  timelevel_qdp, tevolve
    use diffusion_mod, only :  prim_diffusion
#   ifndef CAM
      use asp_tests, only : asp_advection_vertical
#   else
      use control_mod, only : prescribed_vertwind
#   endif
    implicit none
    type (element_t)     , intent(inout), target :: elem(:)
    type (derivative_t)  , intent(in   ) :: deriv
    type (hvcoord_t)     , intent(in   ) :: hvcoord
    type (hybrid_t)      , intent(in   ) :: hybrid
    real (kind=real_kind), intent(in   ) :: dt
    type (TimeLevel_t)   , intent(in   ) :: tl
    integer              , intent(in   ) :: nets
    integer              , intent(in   ) :: nete
    logical              , intent(in   ) :: compute_diagnostics
    ! Local
    real (kind=real_kind) :: dt2, time, dt_vis, x, eta_ave_w
    real (kind=real_kind) :: eta_dot_dpdn(np,np,nlevp)
    real (kind=real_kind) :: dp(np,np)
    real (kind=real_kind) :: tempdp3d(np,np)
    real (kind=real_kind) :: tempmass(nc,nc)
    real (kind=real_kind) :: tempflux(nc,nc,4)
    real (kind=real_kind) :: deta
    integer :: ie,nm1,n0,np1,nstep,method,qsplit_stage,k, qn0
    integer :: n,i,j,lx,lenx
    !JMD    call t_barrierf('sync_prim_advance_exp', hybrid%par%comm)
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
    ! default weights for computing mean dynamics fluxes
    eta_ave_w = 1d0/qsplit
    if(tstep_type==0)then
      call abortmp('ERROR: tstep_type == 0,1 not supported in OpenACC!')
    else if (tstep_type==1) then
      call abortmp('ERROR: tstep_type == 0,1 not supported in OpenACC!')
    else
      method = tstep_type                ! other RK variants
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    ! fix dynamical variables, skip dynamics
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    if (1==prescribed_wind .and. .not.se_prescribed_wind_2d) then
      call abortmp('ERROR: prescribed_wind==1 not supported in OpenACC!')
    endif

    ! ==================================
    ! Take timestep
    ! ==================================
    dt_vis = dt
    if (method==0) then
      call abortmp('ERROR: tstep_type == 0,1 not supported in OpenACC!')
    else if (method==1) then
      call abortmp('ERROR: tstep_type == 0,1 not supported in OpenACC!')
    else if (method==2) then
      call abortmp('ERROR: tstep_type == 2 not supported in OpenACC!')
    else if (method==3) then
      call abortmp('ERROR: tstep_type == 3 not supported in OpenACC!')
    else if (method==4) then
      call abortmp('ERROR: tstep_type == 4 not supported in OpenACC!')
    else if (method==5) then
      ! Ullrich 3nd order 5 stage:   CFL=sqrt( 4^2 -1) = 3.87
      ! u1 = u0 + dt/5 RHS(u0)  (save u1 in timelevel nm1)
      call t_startf("U3-5stage_timestep")
      ! phl: rhs: t=t
      call compute_and_apply_rhs(nm1,n0,n0,qn0,dt/5,elem,hvcoord,hybrid,deriv,nets,nete,compute_diagnostics,eta_ave_w/4)
      ! u2 = u0 + dt/5 RHS(u1)
      ! phl: rhs: t=t+dt/5
      call compute_and_apply_rhs(np1,n0,nm1,qn0,dt/5,elem,hvcoord,hybrid,deriv,nets,nete,.false.,0d0)
      ! u3 = u0 + dt/3 RHS(u2)
      ! phl: rhs: t=t+2*dt/5
      call compute_and_apply_rhs(np1,n0,np1,qn0,dt/3,elem,hvcoord,hybrid,deriv,nets,nete,.false.,0d0)
      ! u4 = u0 + 2dt/3 RHS(u3)
      ! phl: rhs: t=t+2*dt/5+dt/3
      call compute_and_apply_rhs(np1,n0,np1,qn0,2*dt/3,elem,hvcoord,hybrid,deriv,nets,nete,.false.,0d0)
      ! compute (5*u1/4 - u0/4) in timelevel nm1:
      do ie=nets,nete
         state_v   (:,:,:,:,nm1,ie) = (5*state_v   (:,:,:,:,nm1,ie) - state_v   (:,:,:,:,n0,ie) ) / 4
         state_T   (:,:  ,:,nm1,ie) = (5*state_T   (:,:  ,:,nm1,ie) - state_T   (:,:  ,:,n0,ie) ) / 4
         state_dp3d(:,:  ,:,nm1,ie) = (5*state_dp3d(:,:  ,:,nm1,ie) - state_dp3d(:,:  ,:,n0,ie) ) / 4
      enddo
      ! u5 = (5*u1/4 - u0/4) + 3dt/4 RHS(u4)
      ! phl: rhs: t=t+2*dt/5+dt/3+3*dt/4         -wrong RK times ...
      call compute_and_apply_rhs(np1,nm1,np1,qn0,3*dt/4,elem,hvcoord,hybrid,deriv,nets,nete,.false.,3*eta_ave_w/4)
      ! final method is the same as:
      ! u5 = u0 +  dt/4 RHS(u0)) + 3dt/4 RHS(u4)
      call t_stopf("U3-5stage_timestep")
    else if ((method==11).or.(method==12)) then
      ! Fully implicit JFNK method (vertically langragian not active yet)
      call abortmp('ERROR: implicit not supported in OpenACC!')
    else
      call abortmp('ERROR: bad choice of tstep_type')
    endif

    ! ==============================================
    ! Time-split Horizontal diffusion: nu.del^2 or nu.del^4
    ! U(*) = U(t+1)  + dt2 * HYPER_DIFF_TERM(t+1)
    ! ==============================================
#   ifdef ENERGY_DIAGNOSTICS
      !if (compute_diagnostics) then
      !  do ie = nets,nete
      !    elem(ie)%accum%DIFF(:,:,:,:)=elem(ie)%state%v(:,:,:,:,np1)
      !    elem(ie)%accum%DIFFT(:,:,:)=elem(ie)%state%T(:,:,:,np1)
      !  enddo
      !endif
#   endif
    ! note:time step computes u(t+1)= u(t*) + RHS.
    ! for consistency, dt_vis = t-1 - t*, so this is timestep method dependent
    if (tstep_type==0) then
      call abortmp('ERROR: tstep_type == 0 not supported in OpenACC!')
    else if (method<=10) then ! not implicit
      if (rsplit==0) then
        call abortmp('ERROR: rsplit==0 not supported in OpenACC!')
      else
        ! forward-in-time, hypervis applied to dp3d
        call advance_hypervis_dp(edge3p1,elem,hvcoord,hybrid,deriv,np1,nets,nete,dt_vis,eta_ave_w)
      endif
    endif
#   ifdef ENERGY_DIAGNOSTICS
      !if (compute_diagnostics) then
      !  do ie = nets,nete
      !    do k=1,nlev  !  Loop index added (AAM)
      !     elem(ie)%accum%DIFF(:,:,:,k)=( elem(ie)%state%v(:,:,:,k,np1) - elem(ie)%accum%DIFF(:,:,:,k) ) / dt_vis
      !     elem(ie)%accum%DIFFT(:,:,k)=( elem(ie)%state%T(:,:,k,np1) - elem(ie)%accum%DIFFT(:,:,k) ) / dt_vis
      !    enddo
      !  enddo
      !endif
#   endif
    tevolve=tevolve+dt
    call t_stopf('prim_advance_exp')
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
    !
    ! ===================================
    use kinds             , only : real_kind
    use dimensions_mod    , only : np, nc, nlev, ntrac, max_corner_elem
    use hybrid_mod        , only : hybrid_t
    use element_mod       , only : element_t,PrintElem, derived_vn0, state_v, state_qdp, state_t, state_dp3d, state_ps_v
    use derivative_mod    , only : derivative_t, divergence_sphere, gradient_sphere, vorticity_sphere, gradient_sphere_openacc, divergence_sphere_openacc, vorticity_sphere_openacc
    use derivative_mod    , only : subcell_div_fluxes, subcell_dss_fluxes
    use edge_mod          , only : edgevpack, edgevunpack, edgeDGVunpack
    use edgetype_mod      , only : edgedescriptor_t
    use bndry_mod         , only : bndry_exchangev
    use control_mod       , only : moisture, qsplit, use_cpstar, rsplit, swest
    use hybvcoord_mod     , only : hvcoord_t
    use physical_constants, only : cp, cpwater_vapor, Rgas, kappa
    use physics_mod       , only : virtual_specific_heat, virtual_temperature1d
    use prim_si_mod       , only : preq_vertadv, preq_omega_ps_openacc, preq_hydrostatic_openacc
    use time_mod          , only : tevolve
    use openacc_utils_mod , only : memset
#   if ( defined CAM )
      use control_mod     , only: se_met_nudge_u, se_met_nudge_p, se_met_nudge_t, se_met_tevolve
#   endif
    implicit none
    !$acc routine(virtual_specific_heat) seq
    !$acc routine(virtual_temperature1d) seq
    integer              , intent(in   ) :: np1,nm1,n0,qn0,nets,nete
    real*8               , intent(in   ) :: dt2
    logical              , intent(in   ) :: compute_diagnostics
    type (hvcoord_t)     , intent(in   ) :: hvcoord
    type (hybrid_t)      , intent(in   ) :: hybrid
    type (element_t)     , intent(inout), target :: elem(:)
    type (derivative_t)  , intent(in   ) :: deriv
    real (kind=real_kind), intent(in   ) :: eta_ave_w  ! weighting for eta_dot_dpdn mean flux
    ! local
    real (kind=real_kind), dimension(np,np,2)        :: vtemp     ! generic gradient storage
    real (kind=real_kind), dimension(np,np,2     )   :: v         !                            
    real (kind=real_kind), dimension(np,np)          :: vgrad_T    ! v.grad(T)
    real (kind=real_kind), dimension(np,np,2)        :: grad_ps    ! lat-lon coord version
    real (kind=real_kind), dimension(np,np,nlev+1)   :: ph               ! half level pressures on p-grid
    real (kind=real_kind) :: vtens1(np,np,nlev)
    real (kind=real_kind) :: vtens2(np,np,nlev)
    real (kind=real_kind) :: ttens(np,np,nlev)
    real (kind=real_kind) :: stashdp3d (np,np,nlev)
    real (kind=real_kind) :: tempdp3d  (np,np)
    real (kind=real_kind) :: tempflux  (nc,nc,4)
    real (kind=real_kind) :: cp2,cp_ratio,E,de,Qt,v1,v2
    real (kind=real_kind) :: glnps1,glnps2,gpterm
    real (kind=real_kind) :: u_m_umet, v_m_vmet, t_m_tmet 
    type (EdgeDescriptor_t) :: desc
    integer :: i,j,k,kptr,ie

    call t_startf('compute_and_apply_rhs')
    !$omp barrier
    !$omp master
    do ie = 1 , nelemd
      !$acc update device(state_dp3d(:,:,:,n0,ie),state_v(:,:,:,:,n0,ie),state_T(:,:,:,n0,ie),state_Qdp(:,:,:,1,qn0,ie),elem(ie)%state%phis,elem(ie)%derived%omega_p,elem(ie)%derived%eta_dot_dpdn, &
      !$acc&              elem(ie)%derived%pecnd,elem(ie)%derived%phi)
    enddo
    !$acc update device(derived_vn0)
    !$acc parallel loop gang vector collapse(4) present(dp,state_dp3d)
    do ie=1,nelemd
      ! ============================
      ! compute p and delta p
      ! ============================
      do k=1,nlev
        ! vertically lagrangian code: we advect dp3d instead of ps_v
        ! we also need grad(p) at all levels (not just grad(ps))
        !p(k)= hyam(k)*ps0 + hybm(k)*ps
        !    = .5*(hyai(k+1)+hyai(k))*ps0 + .5*(hybi(k+1)+hybi(k))*ps
        !    = .5*(ph(k+1) + ph(k) )  = ph(k) + dp(k)/2
        !
        ! p(k+1)-p(k) = ph(k+1)-ph(k) + (dp(k+1)-dp(k))/2
        !             = dp(k) + (dp(k+1)-dp(k))/2 = (dp(k+1)+dp(k))/2
        do j = 1 , np
          do i = 1 , np
            dp(i,j,k,ie) = state_dp3d(i,j,k,n0,ie)
          enddo
        enddo
      enddo
    enddo
    !$acc parallel loop gang vector collapse(3) present(p,hvcoord,dp)
    do ie=1,nelemd
      do j = 1 , np
        do i = 1 , np
          do k=1,nlev
            if (k==1) then
              p(i,j,k,ie)=hvcoord%hyai(k)*hvcoord%ps0 + dp(i,j,k,ie)/2
            else
              p(i,j,k,ie)=p(i,j,k-1,ie) + dp(i,j,k-1,ie)/2 + dp(i,j,k,ie)/2
            endif
          enddo
        enddo
      enddo
    enddo
    call gradient_sphere_openacc(p,deriv,elem,grad_p,nlev,1,nelemd,1,1)
    ! ============================
    ! compute vgrad_lnps
    ! ============================
    !$acc parallel loop gang vector collapse(4) private(v1,v2) present(elem,vgrad_p,grad_p,vdp,dp,state_v)
    do ie=1,nelemd
      do k=1,nlev
        do j=1,np
          do i=1,np
            v1 = state_v(i,j,1,k,n0,ie)
            v2 = state_v(i,j,2,k,n0,ie)
            vgrad_p(i,j,k,ie) = (v1*grad_p(i,j,1,k,ie) + v2*grad_p(i,j,2,k,ie))
            vdp(i,j,1,k,ie) = v1*dp(i,j,k,ie)
            vdp(i,j,2,k,ie) = v2*dp(i,j,k,ie)
          enddo
        enddo
      enddo
    enddo
#   if ( defined CAM )
      if (se_met_nudge_p.gt.0.D0) then
        !$acc update host(grad_p)
        do ie=1,nelemd
          do k=1,nlev
            ! ============================
            ! compute grad(P-P_met)
            ! ============================
            grad_p_m_pmet(:,:,:,k,ie) = grad_p(:,:,:,k,ie) - hvcoord%hybm(k)* &
                 gradient_sphere( elem(ie)%derived%ps_met(:,:)+tevolve*elem(ie)%derived%dpsdt_met(:,:),deriv,elem(ie)%Dinv)
          enddo
        enddo
      endif
#   endif
    ! ================================
    ! Accumulate mean Vel_rho flux in vn0
    ! ================================
    !$acc parallel loop gang vector collapse(4) present(elem,vdp,derived_vn0)
    do ie=1,nelemd
      do k=1,nlev
        do j=1,np
          do i=1,np
            derived_vn0(i,j,1,k,ie)=derived_vn0(i,j,1,k,ie)+eta_ave_w*vdp(i,j,1,k,ie)
            derived_vn0(i,j,2,k,ie)=derived_vn0(i,j,2,k,ie)+eta_ave_w*vdp(i,j,2,k,ie)
          enddo
        enddo
      enddo
    enddo
    ! =========================================
    ! Compute relative vorticity and divergence
    ! =========================================
    call divergence_sphere_openacc(vdp,deriv,elem,divdp,nlev,1,nelemd,1,1)
    call vorticity_sphere_openacc(state_v,deriv,elem,vort,nlev,1,nelemd,timelevels,n0,1,1)
    if (qn0 == -1 ) then
      ! compute T_v for timelevel n0
      !$acc parallel loop gang vector collapse(4) present(t_v,elem,kappa_star,state_t)
      do ie=1,nelemd
        do k=1,nlev
          do j=1,np
            do i=1,np
              T_v(i,j,k,ie) = state_T(i,j,k,n0,ie)
              kappa_star(i,j,k,ie) = kappa
            enddo
          enddo
        enddo
      enddo
    else
      !$acc parallel loop gang vector collapse(4) private(qt) present(t_v,elem,kappa_star,state_qdp,state_t,dp)
      do ie=1,nelemd
        do k=1,nlev
          do j=1,np
            do i=1,np
              Qt = state_Qdp(i,j,k,1,qn0,ie)/dp(i,j,k,ie)
              T_v(i,j,k,ie) = Virtual_Temperature1d(state_T(i,j,k,n0,ie),Qt)
              if (use_cpstar==1) then
                kappa_star(i,j,k,ie) = Rgas/Virtual_Specific_Heat(Qt)
              else
                kappa_star(i,j,k,ie) = kappa
              endif
            enddo
          enddo
        enddo
      enddo
    endif
    ! ====================================================
    ! Compute Hydrostatic equation, modeld after CCM-3
    ! ====================================================
    call preq_hydrostatic_openacc(elem,T_v,p,dp,1,nelemd)
    ! ====================================================
    ! Compute omega_p according to CCM-3
    ! ====================================================
    call preq_omega_ps_openacc(omega_p,p,vgrad_p,divdp,1,nelemd)
    call memset(product(shape(eta_dot_dpdn)),eta_dot_dpdn,0._real_kind)
    call memset(product(shape(t_vadv      )),t_vadv      ,0._real_kind)
    call memset(product(shape(v_vadv      )),v_vadv      ,0._real_kind)
    call memset(product(shape(sdot_sum    )),sdot_sum    ,0._real_kind)
    ! ================================
    ! accumulate mean vertical flux:
    ! ================================
    !$acc parallel loop gang vector collapse(4) present(elem,eta_dot_dpdn,omega_p)
    do ie=1,nelemd
      do k=1,nlev+1  !  Loop index added (AAM)
        do j = 1 , np
          do i = 1 , np
            elem(ie)%derived%eta_dot_dpdn(i,j,k) = elem(ie)%derived%eta_dot_dpdn(i,j,k) + eta_ave_w*eta_dot_dpdn(i,j,k,ie)
            if (k <= nlev) elem(ie)%derived%omega_p(i,j,k) = elem(ie)%derived%omega_p(i,j,k) + eta_ave_w*omega_p(i,j,k,ie)
          enddo
        enddo
      enddo
    enddo
    !$acc parallel loop gang vector collapse(4) private(v1,v2,E) present(elem,ephi,state_v)
    do ie=1,nelemd
      ! ==============================================
      ! Compute phi + kinetic energy term: 10*nv*nv Flops
      ! ==============================================
      do k=1,nlev
        do j=1,np
          do i=1,np
            v1 = state_v(i,j,1,k,n0,ie)
            v2 = state_v(i,j,2,k,n0,ie)
            E = 0.5D0*( v1*v1 + v2*v2 )
            Ephi(i,j,k,ie)=E+elem(ie)%derived%phi(i,j,k)+elem(ie)%derived%pecnd(i,j,k)
          enddo
        enddo
      enddo
    enddo




    !$acc update host(dp,p,grad_p,vgrad_p,vdp,derived_vn0,divdp,vort,kappa_star,t_v,omega_p,eta_dot_dpdn,t_vadv,v_vadv,sdot_sum,ephi)
    do ie = 1 , nelemd
      !$acc update host(elem(ie)%derived%phi,elem(ie)%derived%omega_p,elem(ie)%derived%eta_dot_dpdn)
    enddo
    !$omp end master
    !$omp barrier


    do ie=nets,nete
      do k=1,nlev
        ! ================================================
        ! compute gradp term (ps/p)*(dp/dps)*T
        ! ================================================
        vtemp(:,:,:) = gradient_sphere(state_T(:,:,k,n0,ie),deriv,elem(ie)%Dinv)
        do j=1,np
          do i=1,np
            v1     = state_v(i,j,1,k,n0,ie)
            v2     = state_v(i,j,2,k,n0,ie)
            vgrad_T(i,j) =  v1*vtemp(i,j,1) + v2*vtemp(i,j,2)
          enddo
        enddo
        ! vtemp = grad ( E + PHI )
        vtemp = gradient_sphere(Ephi(:,:,k,ie),deriv,elem(ie)%Dinv)
        do j=1,np
          do i=1,np
            gpterm = T_v(i,j,k,ie)/p(i,j,k,ie)
            glnps1 = Rgas*gpterm*grad_p(i,j,1,k,ie)
            glnps2 = Rgas*gpterm*grad_p(i,j,2,k,ie)
            v1     = state_v(i,j,1,k,n0,ie)
            v2     = state_v(i,j,2,k,n0,ie)
            vtens1(i,j,k) = - v_vadv(i,j,1,k,ie) + v2*(elem(ie)%fcor(i,j) + vort(i,j,k,ie)) - vtemp(i,j,1) - glnps1
            ! phl: add forcing term to zonal wind u
            vtens2(i,j,k) = - v_vadv(i,j,2,k,ie) - v1*(elem(ie)%fcor(i,j) + vort(i,j,k,ie)) - vtemp(i,j,2) - glnps2
            ! phl: add forcing term to meridional wind v
            ttens(i,j,k)  = - T_vadv(i,j,k,ie) - vgrad_T(i,j) + kappa_star(i,j,k,ie)*T_v(i,j,k,ie)*omega_p(i,j,k,ie)
            ! phl: add forcing term to T
#           if ( defined CAM )
              if (se_prescribed_wind_2d) then
                 vtens1(i,j,k) = 0.D0
                 vtens2(i,j,k) = 0.D0
                 ttens(i,j,k) = 0.D0
              else
                if(se_met_nudge_u.gt.0.D0)then
                   u_m_umet = v1 - elem(ie)%derived%u_met(i,j,k) - se_met_tevolve*tevolve*elem(ie)%derived%dudt_met(i,j,k)
                   v_m_vmet = v2 - elem(ie)%derived%v_met(i,j,k) - se_met_tevolve*tevolve*elem(ie)%derived%dvdt_met(i,j,k)
                   vtens1(i,j,k) =   vtens1(i,j,k) - se_met_nudge_u*u_m_umet * elem(ie)%derived%nudge_factor(i,j,k)
                   elem(ie)%derived%Utnd(i+(j-1)*np,k) = elem(ie)%derived%Utnd(i+(j-1)*np,k) + se_met_nudge_u*u_m_umet * elem(ie)%derived%nudge_factor(i,j,k)
                   vtens2(i,j,k) =   vtens2(i,j,k) - se_met_nudge_u*v_m_vmet * elem(ie)%derived%nudge_factor(i,j,k)
                   elem(ie)%derived%Vtnd(i+(j-1)*np,k) = elem(ie)%derived%Vtnd(i+(j-1)*np,k) + se_met_nudge_u*v_m_vmet * elem(ie)%derived%nudge_factor(i,j,k)
                endif
                if(se_met_nudge_p.gt.0.D0)then
                   vtens1(i,j,k) =   vtens1(i,j,k) - se_met_nudge_p*grad_p_m_pmet(i,j,1,k,ie)  * elem(ie)%derived%nudge_factor(i,j,k)
                   vtens2(i,j,k) =   vtens2(i,j,k) - se_met_nudge_p*grad_p_m_pmet(i,j,2,k,ie)  * elem(ie)%derived%nudge_factor(i,j,k)
                endif
                if(se_met_nudge_t.gt.0.D0)then
                   t_m_tmet = state_T(i,j,k,n0,ie) - elem(ie)%derived%T_met(i,j,k) - se_met_tevolve*tevolve*elem(ie)%derived%dTdt_met(i,j,k)
                   ttens(i,j,k)  = ttens(i,j,k) - se_met_nudge_t*t_m_tmet * elem(ie)%derived%nudge_factor(i,j,k)
                   elem(ie)%derived%Ttnd(i+(j-1)*np,k) = elem(ie)%derived%Ttnd(i+(j-1)*np,k) + se_met_nudge_t*t_m_tmet * elem(ie)%derived%nudge_factor(i,j,k)
                endif
              endif
#           endif
          enddo
        enddo
      enddo
      ! =========================================================
      ! local element timestep, store in np1.
      ! note that we allow np1=n0 or nm1
      ! apply mass matrix
      ! =========================================================
      if (dt2<0) then
        ! calling program just wanted DSS'd RHS, skip time advance
        do k=1,nlev
          state_v   (:,:,1,k,np1,ie) =  elem(ie)%spheremp(:,:)*vtens1(:,:,k)
          state_v   (:,:,2,k,np1,ie) =  elem(ie)%spheremp(:,:)*vtens2(:,:,k)
          state_T   (:,:  ,k,np1,ie) =  elem(ie)%spheremp(:,:)*ttens(:,:,k)
          state_dp3d(:,:  ,k,np1,ie) = -elem(ie)%spheremp(:,:)*(divdp(:,:,k,ie) + eta_dot_dpdn(:,:,k+1,ie)-eta_dot_dpdn(:,:,k,ie))
        enddo
      else
        do k=1,nlev
          state_v   (:,:,1,k,np1,ie) = elem(ie)%spheremp(:,:) * ( state_v   (:,:,1,k,nm1,ie) + dt2*vtens1(:,:,k) )
          state_v   (:,:,2,k,np1,ie) = elem(ie)%spheremp(:,:) * ( state_v   (:,:,2,k,nm1,ie) + dt2*vtens2(:,:,k) )
          state_T   (:,:  ,k,np1,ie) = elem(ie)%spheremp(:,:) * ( state_T   (:,:  ,k,nm1,ie) + dt2*ttens (:,:,k) )
          state_dp3d(:,:  ,k,np1,ie) = elem(ie)%spheremp(:,:) * ( state_dp3d(:,:  ,k,nm1,ie) - dt2*(divdp(:,:,k,ie) + eta_dot_dpdn(:,:,k+1,ie)-eta_dot_dpdn(:,:,k,ie)))
        enddo
      endif
      ! =========================================================
      ! Pack ps(np1), T, and v tendencies into comm buffer
      ! =========================================================
      kptr=0          ;  call edgeVpack(edge3p1, state_ps_v(:,:    ,np1,ie),1     ,kptr,ie)
      kptr=1          ;  call edgeVpack(edge3p1, state_T   (:,:  ,:,np1,ie),nlev  ,kptr,ie)
      kptr=nlev+1     ;  call edgeVpack(edge3p1, state_v   (:,:,:,:,np1,ie),2*nlev,kptr,ie)
      kptr=kptr+2*nlev;  call edgeVpack(edge3p1, state_dp3d(:,:  ,:,np1,ie),nlev  ,kptr,ie)

    enddo

    ! =============================================================
    ! Insert communications here: for shared memory, just a single
    ! sync is required
    ! =============================================================
    call t_startf('caar_bexchV')
    call bndry_exchangeV(hybrid,edge3p1)
    call t_stopf('caar_bexchV')

    do ie=nets,nete
      ! ===========================================================
      ! Unpack the edges for vgrad_T and v tendencies...
      ! ===========================================================
      kptr=0          ;  call edgeVunpack(edge3p1, state_ps_v(:,:    ,np1,ie), 1     , kptr, ie)
      kptr=1          ;  call edgeVunpack(edge3p1, state_T   (:,:  ,:,np1,ie), nlev  , kptr, ie)
      kptr=nlev+1     ;  call edgeVunpack(edge3p1, state_v   (:,:,:,:,np1,ie), 2*nlev, kptr, ie)
      kptr=kptr+2*nlev;  call edgeVunpack(edge3p1, state_dp3d(:,:  ,:,np1,ie), nlev  , kptr, ie)
      ! ====================================================
      ! Scale tendencies by inverse mass matrix
      ! ====================================================
      do k=1,nlev
        state_T   (:,:  ,k,np1,ie) = elem(ie)%rspheremp(:,:)*state_T   (:,:  ,k,np1,ie)
        state_v   (:,:,1,k,np1,ie) = elem(ie)%rspheremp(:,:)*state_v   (:,:,1,k,np1,ie)
        state_v   (:,:,2,k,np1,ie) = elem(ie)%rspheremp(:,:)*state_v   (:,:,2,k,np1,ie)
        state_dp3d(:,:  ,k,np1,ie) = elem(ie)%rspheremp(:,:)*state_dp3d(:,:  ,k,np1,ie)
      enddo
    enddo
    call t_stopf('compute_and_apply_rhs')
  end subroutine compute_and_apply_rhs



  subroutine advance_hypervis(edge3,elem,hvcoord,hybrid,deriv,nt,nets,nete,dt2,eta_ave_w)
  !  take one timestep of:
  !          u(:,:,:,np) = u(:,:,:,np) +  dt2*nu*laplacian**order ( u )
  !          T(:,:,:,np) = T(:,:,:,np) +  dt2*nu_s*laplacian**order ( T )
  !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
  use dimensions_mod    , only: np, np, nlev
  use control_mod       , only: nu, nu_div, nu_s, hypervis_order, hypervis_subcycle, nu_p, nu_top, psurf_vis
  use hybrid_mod        , only: hybrid_t
  use hybvcoord_mod     , only: hvcoord_t
  use element_mod       , only: element_t
  use derivative_mod    , only: derivative_t, laplace_sphere_wk, vlaplace_sphere_wk
  use edge_mod          , only: edgevpack, edgevunpack
  use edgetype_mod      , only: EdgeBuffer_t
  use bndry_mod         , only: bndry_exchangev
  use viscosity_mod     , only: biharmonic_wk
  use physical_constants, only: Cp
  implicit none
  type (hybrid_t)      , intent(in   ) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (EdgeBuffer_t)  , intent(inout) :: edge3
  type (derivative_t)  , intent(in   ) :: deriv
  type (hvcoord_t)     , intent(in   ) :: hvcoord
  real (kind=real_kind), intent(in   ) :: dt2
  integer              , intent(in   ) :: nets,nete,nt
  real (kind=real_kind), intent(in   ) :: eta_ave_w  ! weighting for mean flux terms
  call abortmp('ERROR: rsplit==0 not supported in OpenACC!')
  end subroutine advance_hypervis




  subroutine advance_hypervis_dp(edge3,elem,hvcoord,hybrid,deriv,nt,nets,nete,dt2,eta_ave_w)
  !  take one timestep of:
  !          u(:,:,:,np) = u(:,:,:,np) +  dt2*nu*laplacian**order ( u )
  !          T(:,:,:,np) = T(:,:,:,np) +  dt2*nu_s*laplacian**order ( T )
  !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
  use dimensions_mod    , only: np, np, nlev, nc, ntrac, max_corner_elem
  use control_mod       , only: nu, nu_div, nu_s, hypervis_order, hypervis_subcycle, nu_p, nu_top, psurf_vis, swest
  use hybrid_mod        , only: hybrid_t
  use hybvcoord_mod     , only: hvcoord_t
  use element_mod       , only: element_t
  use derivative_mod    , only: derivative_t, laplace_sphere_wk, vlaplace_sphere_wk
  use derivative_mod    , only: subcell_Laplace_fluxes, subcell_dss_fluxes
  use edge_mod          , only: edgevpack, edgevunpack, edgeDGVunpack
  use edgetype_mod      , only: EdgeBuffer_t, EdgeDescriptor_t
  use bndry_mod         , only: bndry_exchangev
  use viscosity_mod     , only: biharmonic_wk_dp3d
  use derivative_mod    , only: subcell_Laplace_fluxes
  use physical_constants, only: Cp
  implicit none
  type (hybrid_t)      , intent(in   ) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (EdgeBuffer_t)  , intent(inout) :: edge3
  type (derivative_t)  , intent(in   ) :: deriv
  type (hvcoord_t)     , intent(in   ) :: hvcoord
  real (kind=real_kind), intent(in   ) :: dt2
  integer              , intent(in   ) :: nets,nete,nt
  real (kind=real_kind), intent(in   ) :: eta_ave_w  ! weighting for mean flux terms
  ! local
  real (kind=real_kind), dimension(np,np,2,nlev,nets:nete) :: vtens
  real (kind=real_kind), dimension(np,np,nlev,nets:nete)   :: ttens
  real (kind=real_kind), dimension(np,np,nlev,nets:nete)   :: dptens
  real (kind=real_kind), dimension(0:np+1,0:np+1,nlev)     :: corners
  real (kind=real_kind), dimension(2,2,2)                  :: cflux
  real (kind=real_kind), dimension(nc,nc,4,nlev,nets:nete) :: dpflux
  real (kind=real_kind), dimension(np,np,nlev)             :: p
  real (kind=real_kind), dimension(np,np)                  :: dptemp1,dptemp2
  real (kind=real_kind), dimension(np,np)                  :: lap_t,lap_dp
  real (kind=real_kind), dimension(np,np,2)                :: lap_v
  real (kind=real_kind), dimension(np,np,nlev)             :: temp
  real (kind=real_kind), dimension(nc,nc,4)                :: laplace_fluxes
  real (kind=real_kind) :: v1,v2,dt,heating,utens_tmp,vtens_tmp,ttens_tmp,dptens_tmp
  real (kind=real_kind) :: dpdn,dpdn0, nu_scale_top
  integer :: k,kptr,i,j,ie,ic
  type (EdgeDescriptor_t) :: desc
  if (nu_s == 0 .and. nu == 0 .and. nu_p==0 ) return;
  call t_startf('advance_hypervis_dp')
  dt=dt2/hypervis_subcycle
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  regular viscosity
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (hypervis_order == 1) then
    call abortmp('ERROR: regular hyperviscosity not supported in OPENACC!')
  endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  hyper viscosity
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! nu_p=0:
  !   scale T dissipaton by dp  (conserve IE, dissipate T^2)
  ! nu_p>0
  !   dont scale:  T equation IE dissipation matches (to truncation error)
  !                IE dissipation from continuity equation
  !                (1 deg: to about 0.1 W/m^2)
  !
  if (hypervis_order == 2) then
    do ic=1,hypervis_subcycle
      call biharmonic_wk_dp3d(elem,dptens,dpflux,ttens,vtens,deriv,edge3,hybrid,nt,nets,nete)
      do ie=nets,nete
        ! comptue mean flux
        if (nu_p>0) then
          elem(ie)%derived%dpdiss_ave(:,:,:)=elem(ie)%derived%dpdiss_ave(:,:,:)+&
               eta_ave_w*elem(ie)%state%dp3d(:,:,:,nt)/hypervis_subcycle
          elem(ie)%derived%dpdiss_biharmonic(:,:,:)=elem(ie)%derived%dpdiss_biharmonic(:,:,:)+&
               eta_ave_w*dptens(:,:,:,ie)/hypervis_subcycle
        endif
        do k=1,nlev
          ! advace in time.
          ! note: DSS commutes with time stepping, so we can time advance and then DSS.
          ! note: weak operators alreayd have mass matrix "included"
          ! add regular diffusion in top 3 layers:
          if (nu_top>0 .and. k<=3) then
            lap_t=laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
            lap_dp=laplace_sphere_wk(elem(ie)%state%dp3d(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
            lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)
          endif
          nu_scale_top = 1
          if (k==1) nu_scale_top=4
          if (k==2) nu_scale_top=2
          do j=1,np
            do i=1,np
              ! biharmonic terms need a negative sign:
              if (nu_top>0 .and. k<=3) then
                utens_tmp=(-nu*vtens(i,j,1,k,ie) + nu_scale_top*nu_top*lap_v(i,j,1))
                vtens_tmp=(-nu*vtens(i,j,2,k,ie) + nu_scale_top*nu_top*lap_v(i,j,2))
                ttens_tmp=(-nu_s*ttens(i,j,k,ie) + nu_scale_top*nu_top*lap_t(i,j) )
                dptens_tmp=(-nu_p*dptens(i,j,k,ie) + nu_scale_top*nu_top*lap_dp(i,j) )
              else
                utens_tmp=-nu*vtens(i,j,1,k,ie)
                vtens_tmp=-nu*vtens(i,j,2,k,ie)
                ttens_tmp=-nu_s*ttens(i,j,k,ie)
                dptens_tmp=-nu_p*dptens(i,j,k,ie)
              endif
              ttens(i,j,k,ie) = ttens_tmp
              dptens(i,j,k,ie) =dptens_tmp
              vtens(i,j,1,k,ie)=utens_tmp
              vtens(i,j,2,k,ie)=vtens_tmp
            enddo
          enddo
          ! NOTE: we will DSS all tendicies, EXCEPT for dp3d, where we DSS the new state
          elem(ie)%state%dp3d(:,:,k,nt) = elem(ie)%state%dp3d(:,:,k,nt)*elem(ie)%spheremp(:,:) + dt*dptens(:,:,k,ie)
        enddo
        kptr=0     ;  call edgeVpack(edge3, ttens(:,:,:,ie),nlev,kptr,ie)
        kptr=nlev  ;  call edgeVpack(edge3,vtens(:,:,:,:,ie),2*nlev,kptr,ie)
        kptr=3*nlev;  call edgeVpack(edge3,elem(ie)%state%dp3d(:,:,:,nt),nlev,kptr,ie)
      enddo

      call t_startf('ahdp_bexchV2')
      call bndry_exchangeV(hybrid,edge3)
      call t_stopf('ahdp_bexchV2')

      do ie=nets,nete
        kptr=0     ;  call edgeVunpack(edge3, ttens(:,:,:,ie), nlev, kptr, ie)
        kptr=nlev  ;  call edgeVunpack(edge3, vtens(:,:,:,:,ie), 2*nlev, kptr, ie)
        kptr=3*nlev;  call edgeVunpack(edge3, elem(ie)%state%dp3d(:,:,:,nt), nlev, kptr, ie)
        ! apply inverse mass matrix, accumulate tendencies
        do k=1,nlev
          vtens(:,:,1,k,ie)=dt*vtens(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)
          vtens(:,:,2,k,ie)=dt*vtens(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)
          ttens(:,:,k,ie)=dt*ttens(:,:,k,ie)*elem(ie)%rspheremp(:,:)
          elem(ie)%state%dp3d(:,:,k,nt)=elem(ie)%state%dp3d(:,:,k,nt)*elem(ie)%rspheremp(:,:)
        enddo
        ! apply hypervis to u -> u+utens:
        ! E0 = dpdn * .5*u dot u + dpdn * T  + dpdn*PHIS
        ! E1 = dpdn * .5*(u+utens) dot (u+utens) + dpdn * (T-X) + dpdn*PHIS
        ! E1-E0:   dpdn (u dot utens) + dpdn .5 utens dot utens   - dpdn X
        !      X = (u dot utens) + .5 utens dot utens
        !  alt:  (u+utens) dot utens
        do k=1,nlev
          do j=1,np
            do i=1,np
              ! update v first (gives better results than updating v after heating)
              elem(ie)%state%v(i,j,:,k,nt)=elem(ie)%state%v(i,j,:,k,nt) + vtens(i,j,:,k,ie)
              v1=elem(ie)%state%v(i,j,1,k,nt)
              v2=elem(ie)%state%v(i,j,2,k,nt)
              heating = (vtens(i,j,1,k,ie)*v1  + vtens(i,j,2,k,ie)*v2 )
              elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt) + ttens(i,j,k,ie)-heating/cp
            enddo
          enddo
        enddo
      enddo
    enddo
  endif
  call t_stopf('advance_hypervis_dp')
  end subroutine advance_hypervis_dp







end module prim_advance_mod


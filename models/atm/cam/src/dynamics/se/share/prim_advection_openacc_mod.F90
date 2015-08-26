!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! OpenACC implementations of some prim_advection_mod routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_advection_openacc_mod
#if USE_OPENACC
  use kinds          , only: real_kind, int_kind, log_kind
  use dimensions_mod , only: np,nlevp,nlev,qsize,qsize_d,max_corner_elem,max_neigh_edges,nelemd
  use element_mod    , only: timelevels
  use edge_mod       , only: EdgeBuffer_t
  implicit none
  private
  real(kind=real_kind), allocatable :: qmin(:,:,:), qmax(:,:,:)
  real(kind=real_kind), allocatable :: dp0(:)
  real(kind=real_kind), allocatable :: Qtens_biharmonic(:,:,:,:,:)
  real(kind=real_kind), allocatable :: Qtens(:,:,:,:,:)
  real(kind=real_kind), allocatable :: Qmin_pack(:,:,:,:,:)
  real(kind=real_kind), allocatable :: Qmax_pack(:,:,:,:,:)
  real(kind=real_kind), allocatable :: grads_tracer(:,:,:,:,:,:)
  real(kind=real_kind), allocatable :: dp_star(:,:,:,:)
  integer(kind=int_kind), allocatable :: desc_putmapP(:,:)
  integer(kind=int_kind), allocatable :: desc_getmapP(:,:)
  logical(kind=log_kind), allocatable :: desc_reverse(:,:)
  type (EdgeBuffer_t) :: edgeAdv, edgeAdvQ3, edgeAdv_p1, edgeAdvQ2, edgeAdv1, edgeAdv3
  integer,parameter :: DSSeta = 1
  integer,parameter :: DSSomega = 2
  integer,parameter :: DSSdiv_vdp_ave = 3
  integer,parameter :: DSSno_var = -1
  integer :: nbuf

  public :: prim_advection_openacc_init
  public :: precompute_divdp_openacc
  public :: euler_step_openacc
  public :: qdp_time_avg_openacc
  public :: copy_qdp_h2d
  public :: copy_qdp_d2h
  public :: advance_hypervis_scalar_openacc

contains

  subroutine copy_qdp_h2d( elem , tl )
    use element_mod, only: element_t, state_qdp
    implicit none
    type(element_t), intent(in) :: elem(:)
    integer        , intent(in) :: tl
    integer :: ie
    !$omp barrier
    !$omp master
    do ie = 1 , nelemd
      !$acc update device(state_qdp(:,:,:,:,tl,ie))
    enddo
    !$omp end master
    !$omp barrier
  end subroutine copy_qdp_h2d

  subroutine copy_qdp_d2h( elem , tl )
    use element_mod, only: element_t, state_qdp
    implicit none
    type(element_t), intent(in) :: elem(:)
    integer        , intent(in) :: tl
    integer :: ie
    !$omp barrier
    !$omp master
    do ie = 1 , nelemd
      !$acc update host(state_qdp(:,:,:,:,tl,ie))
    enddo
    !$omp end master
    !$omp barrier
  end subroutine copy_qdp_d2h

  subroutine copy_qdp1_h2d( elem , tl )
    use element_mod, only: element_t, state_qdp
    implicit none
    type(element_t), intent(in) :: elem(:)
    integer        , intent(in) :: tl
    integer :: ie
    !$omp barrier
    !$omp master
    do ie = 1 , nelemd
      !$acc update device(state_qdp(:,:,:,1,tl,ie))
    enddo
    !$omp end master
    !$omp barrier
  end subroutine copy_qdp1_h2d

  subroutine copy_qdp1_d2h( elem , tl )
    use element_mod, only: element_t, state_qdp
    implicit none
    type(element_t), intent(in) :: elem(:)
    integer        , intent(in) :: tl
    integer :: ie
    !$omp barrier
    !$omp master
    do ie = 1 , nelemd
      !$acc update host(state_qdp(:,:,:,1,tl,ie))
    enddo
    !$omp end master
    !$omp barrier
  end subroutine copy_qdp1_d2h

  subroutine prim_advection_openacc_init(elem,deriv,hvcoord,hybrid)
    use element_mod
    use derivative_mod, only: derivative_t
    use hybvcoord_mod , only: hvcoord_t
    use hybrid_mod    , only: hybrid_t
    use edge_mod      , only: initEdgeBuffer
    implicit none
    type(element_t)   , intent(in) :: elem(:)
    type(derivative_t), intent(in) :: deriv
    type(hvcoord_t)   , intent(in) :: hvcoord
    type(hybrid_t)    , intent(in) :: hybrid
    integer :: k, ie

    nbuf = 4*(np+max_corner_elem)*nelemd

    call initEdgeBuffer(edgeAdvQ3 ,max(nlev,qsize*nlev*3))
    call initEdgeBuffer(edgeAdv1  ,nlev                  )
    call initEdgeBuffer(edgeAdv   ,qsize*nlev            )
    call initEdgeBuffer(edgeAdv_p1,qsize*nlev + nlev     ) 
    call initEdgeBuffer(edgeAdvQ2 ,qsize*nlev*2          )
    call initEdgeBuffer(edgeAdv3  ,nlev*3                )

    !$OMP BARRIER
    !$OMP MASTER

    allocate(desc_putmapP(max_neigh_edges,nelemd))
    allocate(desc_getmapP(max_neigh_edges,nelemd))
    allocate(desc_reverse(max_neigh_edges,nelemd))
    do ie = 1 , nelemd
      desc_putmapP(:,ie) = elem(ie)%desc%putmapP
      desc_getmapP(:,ie) = elem(ie)%desc%getmapP
      desc_reverse(:,ie) = elem(ie)%desc%reverse
    enddo

    allocate(dp_star(np,np  ,nlev,nelemd))
    allocate(qmin(nlev,qsize,nelemd))
    allocate(qmax(nlev,qsize,nelemd))
    allocate(Qtens_biharmonic(np,np,nlev,qsize,nelemd))
    allocate(qtens(np,np,nlev,qsize,nelemd))
    allocate(qmin_pack(np,np,nlev,qsize,nelemd))
    allocate(qmax_pack(np,np,nlev,qsize,nelemd))
    allocate(grads_tracer(np,np,2,nlev,qsize,nelemd))
    allocate(dp0(nlev))

    do k = 1 , nlev
      dp0(k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0
    enddo

    !$acc enter data pcreate(state_v,state_T,state_dp3d,state_lnps,state_ps_v,state_phis,state_Q,state_Qdp,derived_vn0,derived_vstar,derived_dpdiss_biharmonic,&
    !$acc&                   derived_dpdiss_ave,derived_phi,derived_omega_p,derived_eta_dot_dpdn,derived_grad_lnps,derived_zeta,derived_div,derived_dp,&
    !$acc&                   derived_divdp,derived_divdp_proj)
    !$acc enter data pcreate(qmin,qmax,qmin_pack,qmax_pack,qtens_biharmonic,grads_tracer,dp_star,qtens)
    !$acc enter data pcopyin(elem(1:nelemd),deriv,dp0,desc_putmapP,desc_getmapP,desc_reverse)
    !$acc enter data pcopyin(edgeAdvQ3,edgeAdv1,edgeAdv,edgeAdv_p1,edgeAdvQ2)
    !$acc enter data pcopyin(edgeAdvQ3%buf,edgeAdv1%buf,edgeAdv%buf,edgeAdv_p1%buf,edgeAdvQ2%buf)
    !$acc enter data pcopyin(edgeAdvQ3%receive,edgeAdv1%receive,edgeAdv%receive,edgeAdv_p1%receive,edgeAdvQ2%receive)

    !$OMP END MASTER
    !$OMP BARRIER
  end subroutine prim_advection_openacc_init

  subroutine advance_hypervis_scalar_openacc( edgeAdv_dontuse , elem , hvcoord , hybrid , deriv , nt , nt_qdp , nets , nete , dt2 )
    !  hyperviscsoity operator for foward-in-time scheme
    !  take one timestep of:  
    !          Q(:,:,:,np) = Q(:,:,:,np) +  dt2*nu*laplacian**order ( Q )
    !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
    use kinds          , only : real_kind
    use dimensions_mod , only : np, nlev
    use hybrid_mod     , only : hybrid_t
    use element_mod    , only : element_t, derived_dp, derived_divdp_proj, derived_dpdiss_ave, state_qdp
    use derivative_mod , only : derivative_t
    use edge_mod       , only : EdgeBuffer_t
    use perf_mod       , only : t_startf, t_stopf                          ! _EXTERNAL
    use hybvcoord_mod  , only : hvcoord_t
    use control_mod    , only : nu_q, hypervis_order, hypervis_subcycle_q, nu_p
    implicit none
    type (EdgeBuffer_t)  , intent(inout)         :: edgeAdv_dontuse
    type (element_t)     , intent(inout), target :: elem(:)
    type (hvcoord_t)     , intent(in   )         :: hvcoord
    type (hybrid_t)      , intent(in   )         :: hybrid
    type (derivative_t)  , intent(in   )         :: deriv
    integer              , intent(in   )         :: nt
    integer              , intent(in   )         :: nt_qdp
    integer              , intent(in   )         :: nets
    integer              , intent(in   )         :: nete
    real (kind=real_kind), intent(in   )         :: dt2
    ! local
    real (kind=real_kind), dimension(      nlev,qsize,nets:nete) :: min_neigh
    real (kind=real_kind), dimension(      nlev,qsize,nets:nete) :: max_neigh
    integer :: k , kptr , i , j , ie , ic , q
    real (kind=real_kind), dimension(np,np) :: lap_p
    real (kind=real_kind) :: v1,v2,dt
    real (kind=real_kind) :: dp
    integer :: density_scaling = 0
    if ( nu_q           == 0 ) return
    if ( hypervis_order /= 2 ) return
    call t_startf('advance_hypervis_scalar')
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  hyper viscosity  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dt = dt2 / hypervis_subcycle_q

    do ic = 1 , hypervis_subcycle_q
      !$omp barrier
      !$omp master
      !$acc parallel loop gang vector collapse(4) present(derived_dp,derived_divdp_proj,derived_dpdiss_ave,state_qdp)
      do ie = 1 , nelemd
        ! Qtens = Q/dp   (apply hyperviscsoity to dp0 * Q, not Qdp)
        do k = 1 , nlev
          ! various options:
          !   1)  biharmonic( Qdp )
          !   2)  dp0 * biharmonic( Qdp/dp )    
          !   3)  dpave * biharmonic(Q/dp)
          ! For trace mass / mass consistenciy, we use #2 when nu_p=0
          ! and #e when nu_p>0, where dpave is the mean mass flux from the nu_p
          ! contribution from dynamics.
          do j = 1 , np
            do i = 1 , np
              dp = derived_dp(i,j,k,ie) - dt2 * derived_divdp_proj(i,j,k,ie)
              if (nu_p > 0) then
                do q = 1 , qsize
                  Qtens(i,j,k,q,ie) = derived_dpdiss_ave(i,j,k,ie)*state_Qdp(i,j,k,q,nt_qdp,ie) / dp 
                enddo
              else
                do q = 1 , qsize
                  Qtens(i,j,k,q,ie) = dp0(k)*state_Qdp(i,j,k,q,nt_qdp,ie) / dp
                enddo
              endif
            enddo
          enddo
        enddo
      enddo
      !$omp end master
      !$omp barrier
      ! compute biharmonic operator. Qtens = input and output 
      call biharmonic_wk_scalar( elem , Qtens , grads_tracer , deriv , edgeAdv , hybrid , nets , nete )
      !$omp barrier
      !$omp master
      !$acc parallel loop gang vector collapse(5) present(state_qdp,elem(:),qtens)
      do ie = 1 , nelemd
        do q = 1 , qsize
          do k = 1 , nlev
            do j = 1 , np
              do i = 1 , np
                ! advection Qdp.  For mass advection consistency:
                ! DIFF( Qdp) ~   dp0 DIFF (Q)  =  dp0 DIFF ( Qdp/dp )  
                state_Qdp(i,j,k,q,nt_qdp,ie) = state_Qdp(i,j,k,q,nt_qdp,ie) * elem(ie)%spheremp(i,j) - dt * nu_q * Qtens(i,j,k,q,ie)
              enddo
            enddo
          enddo
        enddo
      enddo
      call limiter2d_zero(state_Qdp,2,nt_qdp)
      call edgeVpack(edgeAdv,state_qdp,qsize*nlev,0,desc_putmapP,desc_reverse,2,nt_qdp)
      !$omp end master
      !$omp barrier

      call t_startf('ah_scalar_bexchV')
      call bndry_exchangeV( hybrid , edgeAdv )
      call t_stopf('ah_scalar_bexchV')
      
      !$omp barrier
      !$omp master
      call edgeVunpack(edgeAdv%buf,edgeAdv%nlyr,state_qdp,qsize*nlev,0,desc_getmapP,2,nt_qdp)
      !$acc parallel loop gang vector collapse(5) present(state_qdp,elem(:))
      do ie = 1 , nelemd
        do q = 1 , qsize    
          ! apply inverse mass matrix
          do k = 1 , nlev
            do j = 1 , np
              do i = 1 , np
                state_Qdp(i,j,k,q,nt_qdp,ie) = elem(ie)%rspheremp(i,j) * state_Qdp(i,j,k,q,nt_qdp,ie)
              enddo
            enddo
          enddo
        enddo
      enddo
      !$omp end master
      !$omp barrier
    enddo
    call copy_qdp1_d2h( elem , nt_qdp )
    call t_stopf('advance_hypervis_scalar')
  end subroutine advance_hypervis_scalar_openacc

  subroutine biharmonic_wk_scalar(elem,qtens,grads,deriv,edgeq,hybrid,nets,nete)
    use hybrid_mod    , only: hybrid_t
    use element_mod   , only: element_t
    use edge_mod      , only: edgeBuffer_t
    use derivative_mod, only: derivative_t
    use control_mod   , only: hypervis_scaling
    use perf_mod      , only: t_startf, t_stopf
    implicit none
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! compute weak biharmonic operator
    !    input:  qtens = Q
    !    output: qtens = weak biharmonic of Q
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type (element_t)     , intent(inout) :: elem(:)
    real (kind=real_kind), intent(inout) :: qtens(np,np,nlev,qsize,nelemd)
    real(kind=real_kind) , intent(inout) :: grads(np,np,2,nlev,qsize,nelemd)
    type (derivative_t)  , intent(in   ) :: deriv
    type (EdgeBuffer_t)  , intent(inout) :: edgeq
    type (hybrid_t)      , intent(in   ) :: hybrid
    integer              , intent(in   ) :: nets,nete
    ! local
    integer :: k,kptr,i,j,ie,ic,q
    logical :: var_coef1
    !if tensor hyperviscosity with tensor V is used, then biharmonic operator is (\grad\cdot V\grad) (\grad \cdot \grad) 
    !so tensor is only used on second call to laplace_sphere_wk
    var_coef1 = .true.
    if(hypervis_scaling > 0) var_coef1 = .false.
    !$omp barrier
    !$omp master
    call laplace_sphere_wk(qtens,grads,deriv,elem,var_coef1,qtens,nlev*qsize,nets,nete)
    call edgeVpack(edgeq,qtens,qsize*nlev,0,desc_putmapP,desc_reverse,1,1)
    !$omp end master
    !$omp barrier

    call t_startf('biwksc_bexchV')
    call bndry_exchangeV(hybrid,edgeq)
    call t_stopf('biwksc_bexchV')
    
    !$omp barrier
    !$omp master
    call edgeVunpack(edgeq%buf,edgeq%nlyr,qtens,qsize*nlev,0,desc_getmapP,1,1)
    !$acc parallel loop gang vector collapse(5) present(qtens,elem(:))
    do ie = 1 , nelemd
      ! apply inverse mass matrix, then apply laplace again
      do q = 1 , qsize      
        do k = 1 , nlev    !  Potential loop inversion (AAM)
          do j = 1 , np
            do i = 1 , np
              qtens(i,j,k,q,ie) = elem(ie)%rspheremp(i,j)*qtens(i,j,k,q,ie)
            enddo
          enddo
        enddo
      enddo
    enddo
    call laplace_sphere_wk(qtens,grads,deriv,elem,.true.,qtens,nlev*qsize,nets,nete)
    !$omp end master
    !$omp barrier
  end subroutine biharmonic_wk_scalar

  subroutine qdp_time_avg_openacc( elem , rkstage , n0_qdp , np1_qdp , limiter_option , nu_p , nets , nete )
    use element_mod, only: element_t, state_qdp
    use control_mod, only: limiter_option
    implicit none
    type(element_t)     , intent(inout) :: elem(:)
    integer             , intent(in   ) :: rkstage , n0_qdp , np1_qdp , nets , nete , limiter_option
    real(kind=real_kind), intent(in   ) :: nu_p
    integer :: ie,q,k,j,i
    !$omp barrier
    !$omp master
    !$acc parallel loop gang vector collapse(5) present(state_qdp)
    do ie = 1 , nelemd
      do q = 1 , qsize
        do k = 1 , nlev
          do j = 1 , np
            do i = 1 , np
              state_Qdp(i,j,k,q,np1_qdp,ie) = ( state_Qdp(i,j,k,q,n0_qdp ,ie) + (rkstage-1)*state_Qdp(i,j,k,q,np1_qdp,ie) ) / rkstage
            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end master
    !$omp barrier
    if (limiter_option == 8) call copy_qdp1_d2h( elem , np1_qdp )
    
  end subroutine qdp_time_avg_openacc

  subroutine euler_step_openacc( np1_qdp , n0_qdp , dt , elem , hvcoord , hybrid , deriv , nets , nete , DSSopt , rhs_multiplier )
  ! ===================================
  ! This routine is the basic foward
  ! euler component used to construct RK SSP methods
  !
  !           u(np1) = u(n0) + dt2*DSS[ RHS(u(n0)) ]
  !
  ! n0 can be the same as np1.  
  !
  ! DSSopt = DSSeta or DSSomega:   also DSS eta_dot_dpdn or omega
  !
  ! ===================================
  use kinds          , only: real_kind
  use dimensions_mod , only: np, npdg, nlev
  use hybrid_mod     , only: hybrid_t
  use element_mod    , only: element_t
  use derivative_mod , only: derivative_t, gradient_sphere, vorticity_sphere
  use hybvcoord_mod  , only: hvcoord_t
  use control_mod    , only: limiter_option, nu_p, nu_q
  use perf_mod       , only: t_startf, t_stopf
  use element_mod    , only: derived_dp, derived_divdp_proj, state_qdp, derived_dpdiss_ave, derived_vn0, derived_divdp, derived_dpdiss_biharmonic
  implicit none
  integer              , intent(in   )         :: np1_qdp, n0_qdp
  real (kind=real_kind), intent(in   )         :: dt
  type (element_t)     , intent(inout), target :: elem(:)
  type (hvcoord_t)     , intent(in   )         :: hvcoord
  type (hybrid_t)      , intent(in   )         :: hybrid
  type (derivative_t)  , intent(in   )         :: deriv
  integer              , intent(in   )         :: nets
  integer              , intent(in   )         :: nete
  integer              , intent(in   )         :: DSSopt
  integer              , intent(in   )         :: rhs_multiplier
  ! local
  real(kind=real_kind) :: vstar(2)
  real(kind=real_kind) :: tmp, mintmp, maxtmp, dp
  integer :: ie,q,i,j,k
  integer :: rhs_viss

  !$omp barrier
  !$omp master
  !$acc update device(derived_dp,derived_divdp_proj,derived_vn0,derived_dpdiss_ave)
  !$omp end master
  !$omp barrier

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   compute Q min/max values for lim8
  !   compute biharmonic mixing term f
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  rhs_viss = 0
  if ( limiter_option == 8  ) then
    ! when running lim8, we also need to limit the biharmonic, so that term needs
    ! to be included in each euler step.  three possible algorithms here:
    ! 1) most expensive:
    !     compute biharmonic (which also computes qmin/qmax) during all 3 stages
    !     be sure to set rhs_viss=1
    !     cost:  3 biharmonic steps with 3 DSS
    !
    ! 2) cheapest:
    !     compute biharmonic (which also computes qmin/qmax) only on first stage
    !     be sure to set rhs_viss=3
    !     reuse qmin/qmax for all following stages (but update based on local qmin/qmax)
    !     cost:  1 biharmonic steps with 1 DSS
    !     main concern:  viscosity 
    !     
    ! 3)  compromise:
    !     compute biharmonic (which also computes qmin/qmax) only on last stage
    !     be sure to set rhs_viss=3
    !     compute qmin/qmax directly on first stage
    !     reuse qmin/qmax for 2nd stage stage (but update based on local qmin/qmax)
    !     cost:  1 biharmonic steps, 2 DSS
    !
    !  NOTE  when nu_p=0 (no dissipation applied in dynamics to dp equation), we should
    !        apply dissipation to Q (not Qdp) to preserve Q=1
    !        i.e.  laplace(Qdp) ~  dp0 laplace(Q)                
    !        for nu_p=nu_q>0, we need to apply dissipation to Q * diffusion_dp
    !
    ! initialize dp, and compute Q from Qdp (and store Q in Qtens_biharmonic)
    !$omp barrier
    !$omp master
    !$acc parallel loop gang vector collapse(4) private(tmp) present(derived_dp,derived_divdp_proj,state_qdp)
    do ie = 1 , nelemd
      ! add hyperviscosity to RHS.  apply to Q at timelevel n0, Qdp(n0)/dp
      do k = 1 , nlev    !  Loop index added with implicit inversion (AAM)
        do j = 1 , np
          do i = 1 , np
            tmp = derived_dp(i,j,k,ie) - rhs_multiplier*dt*derived_divdp_proj(i,j,k,ie) 
            do q = 1 , qsize
              Qtens_biharmonic(i,j,k,q,ie) = state_Qdp(i,j,k,q,n0_qdp,ie) / tmp
            enddo
          enddo
        enddo
      enddo
    enddo
    !$acc parallel loop gang vector collapse(3) private(mintmp,maxtmp) present(qmin,qmax,qtens_biharmonic)
    do ie = 1 , nelemd
      do k = 1 , nlev    
        do q = 1 , qsize
          if (rhs_multiplier == 0 .or. rhs_multiplier == 2) then
            mintmp =  1.D20
            maxtmp = -1.D20
          else
            mintmp = qmin(k,q,ie)
            maxtmp = qmax(k,q,ie)
          endif
          do j = 1 , np
            do i = 1 , np
              mintmp = min(mintmp,qtens_biharmonic(i,j,k,q,ie))
              maxtmp = max(maxtmp,qtens_biharmonic(i,j,k,q,ie))
            enddo
          enddo
          qmin(k,q,ie) = max(mintmp,0d0)
          qmax(k,q,ie) =     maxtmp
        enddo
      enddo
    enddo
    !$omp end master
    !$omp barrier
    if ( rhs_multiplier == 0 ) call neighbor_minmax(elem,hybrid,edgeAdvQ2,nets,nete,qmin,qmax)
    ! compute biharmonic mixing term
    if ( rhs_multiplier == 2 ) then
      rhs_viss = 3
      ! two scalings depending on nu_p:
      ! nu_p=0:    qtens_biharmonic *= dp0                   (apply viscsoity only to q)
      ! nu_p>0):   qtens_biharmonc *= elem()%psdiss_ave      (for consistency, if nu_p=nu_q)
      if ( nu_p > 0 ) then
        !$omp barrier
        !$omp master
        !$acc parallel loop gang vector collapse(5) present(qtens_biharmonic,derived_dpdiss_ave,dp0)
        do ie = 1 , nelemd
          do q = 1 , qsize
            do k = 1 , nlev    
              do j = 1 , np
                do i = 1 , np
                  ! NOTE: divide by dp0 since we multiply by dp0 below
                  Qtens_biharmonic(i,j,k,q,ie)=Qtens_biharmonic(i,j,k,q,ie)*derived_dpdiss_ave(i,j,k,ie)/dp0(k)
                enddo
              enddo
            enddo
          enddo
        enddo
        !$omp end master
        !$omp barrier
      endif
      call biharmonic_wk_scalar_minmax( elem , qtens_biharmonic , grads_tracer , deriv , edgeAdvQ3 , hybrid , nets , nete , qmin , qmax )
      !$omp barrier
      !$omp master
      !$acc parallel loop gang vector collapse(5) present(qtens_biharmonic,dp0,elem(:))
      do ie = 1 , nelemd
        do q = 1 , qsize
          do k = 1 , nlev    !  Loop inversion (AAM)
            do j = 1 , np
              do i = 1 , np
                ! note: biharmonic_wk() output has mass matrix already applied. Un-apply since we apply again below:
                qtens_biharmonic(i,j,k,q,ie) = -rhs_viss*dt*nu_q*dp0(k)*Qtens_biharmonic(i,j,k,q,ie) / elem(ie)%spheremp(i,j)
              enddo
            enddo
          enddo
        enddo
      enddo
      !$omp end master
      !$omp barrier
    endif
  endif  ! compute biharmonic mixing term and qmin/qmax


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   2D Advection step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Compute velocity used to advance Qdp 
  !$omp barrier
  !$omp master
  if (limiter_option == 8) then
    !$acc update device(derived_divdp)
    if ( nu_p > 0 .and. rhs_viss /= 0 ) then
      !$acc update device(derived_dpdiss_biharmonic)
    endif
  endif
  !$acc parallel loop gang vector collapse(4) present(derived_dp,derived_divdp_proj,derived_vn0,dp_star,derived_divdp,derived_dpdiss_biharmonic,grads_tracer,state_qdp) &
  !$acc& private(dp,vstar)
  do ie = 1 , nelemd
    do k = 1 , nlev    !  Loop index added (AAM)
      do j = 1 , np
        do i = 1 , np
          ! derived variable divdp_proj() (DSS'd version of divdp) will only be correct on 2nd and 3rd stage
          ! but that's ok because rhs_multiplier=0 on the first stage:
          dp = derived_dp(i,j,k,ie) - rhs_multiplier * dt * derived_divdp_proj(i,j,k,ie)
          Vstar(1) = derived_vn0(i,j,1,k,ie) / dp
          Vstar(2) = derived_vn0(i,j,2,k,ie) / dp
          if ( limiter_option == 8 ) then
            ! UN-DSS'ed dp at timelevel n0+1:  
            dp_star(i,j,k,ie) = dp - dt * derived_divdp(i,j,k,ie)  
            if ( nu_p > 0 .and. rhs_viss /= 0 ) then
              ! add contribution from UN-DSS'ed PS dissipation
              dp_star(i,j,k,ie) = dp_star(i,j,k,ie) - rhs_viss * dt * nu_q * derived_dpdiss_biharmonic(i,j,k,ie) / elem(ie)%spheremp(i,j)
            endif
          endif
          do q = 1 , qsize
            grads_tracer(i,j,1,k,q,ie) = Vstar(1) * state_Qdp(i,j,k,q,n0_qdp,ie)
            grads_tracer(i,j,2,k,q,ie) = Vstar(2) * state_Qdp(i,j,k,q,n0_qdp,ie)
          enddo
        enddo
      enddo
    enddo
  enddo
  call divergence_sphere( grads_tracer , deriv , elem(:) , qtens , nlev*qsize )
  !$acc parallel loop gang vector collapse(5) present(qtens,state_qdp,qtens_biharmonic)
  do ie = 1 , nelemd
    ! advance Qdp
    do q = 1 , qsize
      do k = 1 , nlev  !  dp_star used as temporary instead of divdp (AAM)
        do j = 1 , np
          do i = 1 , np
            Qtens(i,j,k,q,ie) = state_Qdp(i,j,k,q,n0_qdp,ie) - dt * qtens(i,j,k,q,ie)
            ! optionally add in hyperviscosity computed above:
            if ( rhs_viss /= 0 ) Qtens(i,j,k,q,ie) = Qtens(i,j,k,q,ie) + Qtens_biharmonic(i,j,k,q,ie)
          enddo
        enddo
      enddo
    enddo
  enddo
  if ( limiter_option == 8 ) then
    call limiter_optim_iter_full( Qtens , elem(:) , qmin , qmax , dp_star )   ! apply limiter to Q = Qtens / dp_star 
  endif
  !$acc parallel loop gang vector collapse(5) present(state_Qdp,elem(:),qtens)
  do ie = 1 , nelemd
    ! advance Qdp
    do q = 1 , qsize
      ! apply mass matrix, overwrite np1 with solution:
      ! dont do this earlier, since we allow np1_qdp == n0_qdp 
      ! and we dont want to overwrite n0_qdp until we are done using it
      do k = 1 , nlev
        do j = 1 , np
          do i = 1 , np
            state_Qdp(i,j,k,q,np1_qdp,ie) = elem(ie)%spheremp(i,j) * Qtens(i,j,k,q,ie) 
          enddo
        enddo
      enddo
    enddo
  enddo
  if ( limiter_option == 4 ) then
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    ! sign-preserving limiter, applied after mass matrix
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    call limiter2d_zero( state_Qdp , 2 , np1_qdp )
  endif
  ! note: eta_dot_dpdn is actually dimension nlev+1, but nlev+1 data is
  ! all zero so we only have to DSS 1:nlev
  call edgeVpack(edgeAdv , state_Qdp , nlev*qsize , 0 , desc_putmapP , desc_reverse , 2 , np1_qdp )
  !$omp end master
  !$omp barrier

  call t_startf('eus_bexchV')
  call bndry_exchangeV( hybrid , edgeAdv    )
  call t_stopf('eus_bexchV')

  !$omp barrier
  !$omp master
  call edgeVunpack( edgeAdv%buf , edgeAdv%nlyr , state_Qdp , nlev*qsize , 0 , desc_getmapP , 2 , np1_qdp )
  !$acc parallel loop gang vector collapse(5) present(state_Qdp,elem(:))
  do ie = 1 , nelemd
    do q = 1 , qsize
      do k = 1 , nlev    !  Potential loop inversion (AAM)
        do j = 1 , np
          do i = 1 , np
            state_Qdp(i,j,k,q,np1_qdp,ie) = elem(ie)%rspheremp(i,j) * state_Qdp(i,j,k,q,np1_qdp,ie)
          enddo
        enddo
      enddo
    enddo
  enddo
  !$omp end master
  !$omp barrier
  end subroutine euler_step_openacc

  subroutine biharmonic_wk_scalar_minmax(elem,qtens,grads,deriv,edgeq,hybrid,nets,nete,emin,emax)
    use hybrid_mod, only: hybrid_t
    use element_mod, only: element_t
    use derivative_mod, only: derivative_t
    use control_mod, only: hypervis_scaling, hypervis_power
    use perf_mod, only: t_startf, t_stopf
    implicit none
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! compute weak biharmonic operator
    !    input:  qtens = Q
    !    output: qtens = weak biharmonic of Q and Q element min/max
    !
    !    note: emin/emax must be initialized with Q element min/max.  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type (hybrid_t)      , intent(in   ) :: hybrid
    type (element_t)     , intent(in   ) :: elem(:)
    integer              , intent(in   ) :: nets,nete
    real (kind=real_kind), intent(inout) :: qtens(np,np,nlev,qsize,nelemd)
    real(kind=real_kind) , intent(inout) :: grads(np,np,2,nlev,qsize,nelemd)
    type (EdgeBuffer_t)  , intent(inout) :: edgeq
    type (derivative_t)  , intent(in   ) :: deriv
    real (kind=real_kind), intent(inout) :: emin(nlev,qsize,nelemd)
    real (kind=real_kind), intent(inout) :: emax(nlev,qsize,nelemd)
    ! local
    integer :: k,kptr,i,j,ie,ic,q
    real (kind=real_kind) :: lap_p(np,np)
    logical :: var_coef1
    !if tensor hyperviscosity with tensor V is used, then biharmonic operator is (\grad\cdot V\grad) (\grad \cdot \grad) 
    !so tensor is only used on second call to laplace_sphere_wk
    var_coef1 = .true.
    if(hypervis_scaling > 0)    var_coef1 = .false.
    !$omp barrier
    !$omp master
    call minmax_pack(qmin_pack,qmax_pack,emin,emax)
    call laplace_sphere_wk(qtens,grads,deriv,elem,var_coef1,qtens,nlev*qsize,nets,nete)
    call edgeVpack(edgeq,    qtens,qsize*nlev,0           ,desc_putmapP,desc_reverse,1,1)
    call edgeVpack(edgeq,Qmin_pack,nlev*qsize,nlev*qsize  ,desc_putmapP,desc_reverse,1,1)
    call edgeVpack(edgeq,Qmax_pack,nlev*qsize,2*nlev*qsize,desc_putmapP,desc_reverse,1,1)
    !$omp end master
    !$omp barrier

    call t_startf('biwkscmm_bexchV')
    call bndry_exchangeV(hybrid,edgeq)
    call t_stopf('biwkscmm_bexchV')

    !$omp barrier
    !$omp master
    call edgeVunpack   (edgeq%buf,edgeq%nlyr,    qtens,qsize*nlev,0           ,desc_getmapP,1,1)
    call edgeVunpackMin(edgeq%buf,edgeq%nlyr,Qmin_pack,qsize*nlev,qsize*nlev  ,desc_getmapP,1,1)
    call edgeVunpackMax(edgeq%buf,edgeq%nlyr,Qmax_pack,qsize*nlev,2*qsize*nlev,desc_getmapP,1,1)
    !$acc parallel loop gang vector collapse(5) present(qtens,elem(:))
    do ie = 1 , nelemd
      do q = 1 , qsize      
        do k = 1 , nlev
          do j = 1 , np
            do i = 1 , np
              qtens(i,j,k,q,ie) = elem(ie)%rspheremp(i,j)*qtens(i,j,k,q,ie)  ! apply inverse mass matrix
            enddo
          enddo
        enddo
      enddo
    enddo
    call laplace_sphere_wk(qtens,grads,deriv,elem,.true.,qtens,nlev*qsize,nets,nete)
    call minmax_reduce_corners(qmin_pack,qmax_pack,emin,emax)
    !$omp end master
    !$omp barrier
  end subroutine biharmonic_wk_scalar_minmax

  subroutine laplace_sphere_wk(s,grads,deriv,elem,var_coef,laplace,len,nets,nete)
    use derivative_mod, only: derivative_t
    use element_mod, only: element_t
    use control_mod, only: hypervis_scaling, hypervis_power
    implicit none
    !input:  s = scalar
    !ouput:  -< grad(PHI), grad(s) >   = weak divergence of grad(s)
    !note: for this form of the operator, grad(s) does not need to be made C0
    real(kind=real_kind) , intent(in   ) :: s(np,np,len,nelemd)
    real(kind=real_kind) , intent(inout) :: grads(np,np,2,len,nelemd)
    type (derivative_t)  , intent(in   ) :: deriv
    type (element_t)     , intent(in   ) :: elem(:)
    logical              , intent(in   ) :: var_coef
    real(kind=real_kind) , intent(  out) :: laplace(np,np,len,nelemd)
    integer              , intent(in   ) :: len
    integer              , intent(in   ) :: nets,nete
    integer :: i,j,k,ie
    ! Local
    real(kind=real_kind) :: oldgrads(2)
    grads = gradient_sphere(s,deriv,elem(:),len,nets,nete)
    !$acc parallel loop gang vector collapse(4) present(grads,elem(:))
    do ie = 1 , nelemd
      do k = 1 , len
        do j = 1 , np
          do i = 1 , np
            if (var_coef) then
              if (hypervis_power/=0 ) then
                ! scalar viscosity with variable coefficient
                grads(i,j,1,k,ie) = grads(i,j,1,k,ie)*elem(ie)%variable_hyperviscosity(i,j)
                grads(i,j,2,k,ie) = grads(i,j,2,k,ie)*elem(ie)%variable_hyperviscosity(i,j)
              else if (hypervis_scaling /=0 ) then
                oldgrads = grads(i,j,:,k,ie)
                grads(i,j,1,k,ie) = sum(oldgrads(:)*elem(ie)%tensorVisc(1,:,i,j))
                grads(i,j,2,k,ie) = sum(oldgrads(:)*elem(ie)%tensorVisc(2,:,i,j))
              endif
            endif
          enddo
        enddo
      enddo
    enddo
    ! note: divergnece_sphere and divergence_sphere_wk are identical *after* bndry_exchange
    ! if input is C_0.  Here input is not C_0, so we should use divergence_sphere_wk().  
    laplace = divergence_sphere_wk(grads,deriv,elem(:),len,nets,nete)
  end subroutine laplace_sphere_wk

  function divergence_sphere_wk(v,deriv,elem,len,nets,nete) result(div)
    use element_mod, only: element_t
    use derivative_mod, only: derivative_t
    use physical_constants, only: rrearth
    implicit none
!   input:  v = velocity in lat-lon coordinates
!   ouput:  div(v)  spherical divergence of v, integrated by parts
!   Computes  -< grad(psi) dot v > 
!   (the integrated by parts version of < psi div(v) > )
!   note: after DSS, divergence_sphere () and divergence_sphere_wk() 
!   are identical to roundoff, as theory predicts.
    real(kind=real_kind), intent(in) :: v(np,np,2,len,nelemd)  ! in lat-lon coordinates
    type (derivative_t) , intent(in) :: deriv
    type (element_t)    , intent(in) :: elem(:)
    integer             , intent(in) :: len
    integer             , intent(in) :: nets , nete
    real(kind=real_kind)             :: div(np,np,len,nelemd)
    ! Local
    integer, parameter :: kchunk = 8
    integer :: i,j,l,k,ie,kc,kk
    real(kind=real_kind) :: vtemp(np,np,2,kchunk), tmp
    ! latlon- > contra
    !$acc parallel loop gang collapse(2) present(v,elem(:),div) private(vtemp)
    do ie = 1 , nelemd
      do kc = 1 , len/kchunk+1
        !$acc cache(vtemp)
        !$acc loop vector collapse(3)
        do kk = 1 , kchunk
          do j = 1 , np
            do i = 1 , np
              k = (kc-1)*kchunk+kk
              if (k <= len) then
                vtemp(i,j,1,kk)=elem(ie)%spheremp(i,j)*(elem(ie)%Dinv(1,1,i,j)*v(i,j,1,k,ie) + elem(ie)%Dinv(1,2,i,j)*v(i,j,2,k,ie))
                vtemp(i,j,2,kk)=elem(ie)%spheremp(i,j)*(elem(ie)%Dinv(2,1,i,j)*v(i,j,1,k,ie) + elem(ie)%Dinv(2,2,i,j)*v(i,j,2,k,ie))
              endif
            enddo
          enddo
        enddo
        !$acc loop vector collapse(3) private(tmp)
        do kk = 1 , kchunk
          do j = 1 , np
            do i = 1 , np
              k = (kc-1)*kchunk+kk
              if (k <= len) then
                tmp = 0.
                do l = 1 , np
                  tmp = tmp - ( vtemp(l,j,1,kk)*deriv%Dvv(i,l) + vtemp(i,l,2,kk)*deriv%Dvv(j,l) )
                enddo
                div(i,j,k,ie) = tmp * rrearth
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
  end function divergence_sphere_wk

  function gradient_sphere(s,deriv,elem,len,nets,nete) result(ds)
    use element_mod, only: element_t
    use derivative_mod, only: derivative_t
    use physical_constants, only: rrearth
    implicit none
    !   input s:  scalar
    !   output  ds: spherical gradient of s, lat-lon coordinates
    real(kind=real_kind), intent(in) :: s(np,np,len,nelemd)
    type(derivative_t)  , intent(in) :: deriv
    type(element_t)     , intent(in) :: elem(:)
    integer             , intent(in) :: len
    integer             , intent(in) :: nets,nete
    real(kind=real_kind)             :: ds(np,np,2,len,nelemd)
    integer :: i, j, l, k, ie
    real(kind=real_kind) :: dsdx00, dsdy00
    !$acc parallel loop gang vector collapse(4) present(ds,elem(:),s)
    do ie = 1 , nelemd
      do k = 1 , len
        do j = 1 , np
          do i = 1 , np
            dsdx00=0.0d0
            dsdy00=0.0d0
            do l = 1 , np
              dsdx00 = dsdx00 + deriv%Dvv(l,i)*s(l,j,k,ie)
              dsdy00 = dsdy00 + deriv%Dvv(l,j)*s(i,l,k,ie)
            enddo
            ds(i,j,1,k,ie) = ( elem(ie)%Dinv(1,1,i,j)*dsdx00 + elem(ie)%Dinv(2,1,i,j)*dsdy00 ) * rrearth
            ds(i,j,2,k,ie) = ( elem(ie)%Dinv(1,2,i,j)*dsdx00 + elem(ie)%Dinv(2,2,i,j)*dsdy00 ) * rrearth
          enddo
        enddo
      enddo
    enddo
  end function gradient_sphere

  subroutine neighbor_minmax(elem,hybrid,edgeMinMax,nets,nete,min_neigh,max_neigh)
    use hybrid_mod , only: hybrid_t
    use element_mod, only: element_t
    use perf_mod   , only: t_startf, t_stopf
    implicit none
    ! compute Q min&max over the element and all its neighbors
    integer :: nets,nete
    type (element_t)     , intent(in   ) :: elem(:)
    type (hybrid_t)      , intent(in   ) :: hybrid
    type (EdgeBuffer_t)  , intent(inout) :: edgeMinMax
    real (kind=real_kind), intent(inout) :: min_neigh(nlev,qsize,nelemd)
    real (kind=real_kind), intent(inout) :: max_neigh(nlev,qsize,nelemd)
    ! local
    integer :: ie,k,q,j,i
    ! compute Qmin, Qmax
    !$omp barrier
    !$omp master
    call minmax_pack(qmin_pack,qmax_pack,min_neigh,max_neigh)
    call edgeVpack(edgeMinMax,Qmin_pack,nlev*qsize,0         ,desc_putmapP,desc_reverse,1,1)
    call edgeVpack(edgeMinMax,Qmax_pack,nlev*qsize,nlev*qsize,desc_putmapP,desc_reverse,1,1)
    !$omp end master
    !$omp barrier

    call t_startf('nmm_bexchV')
    call bndry_exchangeV(hybrid,edgeMinMax)
    call t_stopf('nmm_bexchV')
       
    !$omp barrier
    !$omp master
    call edgeVunpackMin(edgeMinMax%buf,edgeMinMax%nlyr,Qmin_pack,nlev*qsize,0         ,desc_getmapP,1,1)
    call edgeVunpackMax(edgeMinMax%buf,edgeMinMax%nlyr,Qmax_pack,nlev*qsize,nlev*qsize,desc_getmapP,1,1)
    call minmax_reduce_corners(qmin_pack,qmax_pack,min_neigh,max_neigh)
    !$omp end master
    !$omp barrier
  end subroutine neighbor_minmax

  subroutine minmax_pack(qmin_pack,qmax_pack,min_neigh,max_neigh)
    implicit none
    real(kind=real_kind), intent(  out) :: qmin_pack(np,np,nlev,qsize,nelemd)
    real(kind=real_kind), intent(  out) :: qmax_pack(np,np,nlev,qsize,nelemd)
    real(kind=real_kind), intent(in   ) :: min_neigh(nlev,qsize,nelemd)
    real(kind=real_kind), intent(in   ) :: max_neigh(nlev,qsize,nelemd)
    integer :: ie,q,k,j,i
    !$acc parallel loop gang vector collapse(5) present(min_neigh,max_neigh,qmin_pack,qmax_pack)
    do ie=1,nelemd
      do q=1,qsize
        do k=1,nlev
          do j = 1 , np
            do i = 1 , np
              Qmin_pack(i,j,k,q,ie)=min_neigh(k,q,ie)
              Qmax_pack(i,j,k,q,ie)=max_neigh(k,q,ie)
            enddo
          enddo
        enddo
      enddo
    enddo
  end subroutine minmax_pack

  subroutine minmax_reduce_corners(qmin_pack,qmax_pack,min_neigh,max_neigh)
    implicit none
    real(kind=real_kind), intent(in   ) :: qmin_pack(np,np,nlev,qsize,nelemd)
    real(kind=real_kind), intent(in   ) :: qmax_pack(np,np,nlev,qsize,nelemd)
    real(kind=real_kind), intent(  out) :: min_neigh(nlev,qsize,nelemd)
    real(kind=real_kind), intent(  out) :: max_neigh(nlev,qsize,nelemd)
    integer :: ie,q,k
    !$acc parallel loop gang vector collapse(3) present(min_neigh,max_neigh,qmin_pack,qmax_pack)
    do ie=1,nelemd
      do q=1,qsize
        do k=1,nlev
          ! note: only need to consider the corners, since the data we packed was
          ! constant within each element
          min_neigh(k,q,ie)=min(qmin_pack(1,1,k,q,ie),qmin_pack(1,np,k,q,ie),qmin_pack(np,1,k,q,ie),qmin_pack(np,np,k,q,ie))
          min_neigh(k,q,ie)=max(min_neigh(k,q,ie),0d0)
          max_neigh(k,q,ie)=max(qmax_pack(1,1,k,q,ie),qmax_pack(1,np,k,q,ie),qmax_pack(np,1,k,q,ie),qmax_pack(np,np,k,q,ie))
        enddo
      enddo
    enddo
  end subroutine minmax_reduce_corners

  subroutine edgeVpack(edge,v,vlyr,kptr,putmapP,reverse,tdim,tl)
    use dimensions_mod, only : max_corner_elem
    use control_mod   , only : north, south, east, west, neast, nwest, seast, swest
    use perf_mod      , only : t_startf, t_stopf
    use parallel_mod  , only : haltmp
    use element_mod   , only : Element_t
    use edge_mod      , only : EdgeBuffer_t
    type (EdgeBuffer_t)    ,intent(inout) :: edge
    integer                ,intent(in   ) :: vlyr
    real (kind=real_kind)  ,intent(in   ) :: v(np,np,vlyr,tdim,nelemd)
    integer                ,intent(in   ) :: kptr
    integer(kind=int_kind) ,intent(in   ) :: putmapP(max_neigh_edges,nelemd)
    logical(kind=log_kind) ,intent(in   ) :: reverse(max_neigh_edges,nelemd)
    integer                ,intent(in   ) :: tdim,tl
    ! Local variables
    integer :: i,k,ir,ll,is,ie,in,iw,el
    call t_startf('edge_pack')
    if (edge%nlyr < (kptr+vlyr) ) call haltmp('edgeVpack: Buffer overflow: size of the vertical dimension must be increased!')
    !$acc parallel loop gang vector collapse(2) present(v,putmapP,reverse,edge%buf)
    do el = 1 , nelemd
      do k = 1 , vlyr
        do i = 1 , np
          edge%buf(kptr+k,putmapP(south,el)+i) = v(i ,1 ,k,tl,el)
          edge%buf(kptr+k,putmapP(east ,el)+i) = v(np,i ,k,tl,el)
          edge%buf(kptr+k,putmapP(north,el)+i) = v(i ,np,k,tl,el)
          edge%buf(kptr+k,putmapP(west ,el)+i) = v(1 ,i ,k,tl,el)
        enddo
        do i = 1 , np
          ir = np-i+1
          if(reverse(south,el)) edge%buf(kptr+k,putmapP(south,el)+ir) = v(i ,1 ,k,tl,el)
          if(reverse(east ,el)) edge%buf(kptr+k,putmapP(east ,el)+ir) = v(np,i ,k,tl,el)
          if(reverse(north,el)) edge%buf(kptr+k,putmapP(north,el)+ir) = v(i ,np,k,tl,el)
          if(reverse(west ,el)) edge%buf(kptr+k,putmapP(west ,el)+ir) = v(1 ,i ,k,tl,el)
        enddo
        do i = 1 , max_corner_elem
          ll = swest+0*max_corner_elem+i-1; if(putmapP(ll,el) /= -1) edge%buf(kptr+k,putmapP(ll,el)+1) = v(1 ,1 ,k,tl,el)
          ll = swest+1*max_corner_elem+i-1; if(putmapP(ll,el) /= -1) edge%buf(kptr+k,putmapP(ll,el)+1) = v(np,1 ,k,tl,el)
          ll = swest+2*max_corner_elem+i-1; if(putmapP(ll,el) /= -1) edge%buf(kptr+k,putmapP(ll,el)+1) = v(1 ,np,k,tl,el)
          ll = swest+3*max_corner_elem+i-1; if(putmapP(ll,el) /= -1) edge%buf(kptr+k,putmapP(ll,el)+1) = v(np,np,k,tl,el)
        enddo
      enddo
    enddo
    call t_stopf('edge_pack')
  end subroutine edgeVpack

  subroutine edgeVunpack(edgebuf,nlyr,v,vlyr,kptr,getmapP,tdim,tl)
    use dimensions_mod, only : np, max_corner_elem
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    use perf_mod, only: t_startf, t_stopf
    real(kind=real_kind)  , intent(in   ) :: edgebuf(nlyr,nbuf)
    integer               , intent(in   ) :: nlyr
    integer               , intent(in   ) :: vlyr
    real(kind=real_kind)  , intent(inout) :: v(np,np,vlyr,tdim,nelemd)
    integer               , intent(in   ) :: kptr
    integer(kind=int_kind), intent(in   ) :: getmapP(max_neigh_edges,nelemd)
    integer                ,intent(in   ) :: tdim,tl
    ! Local
    integer :: i,k,ll,is,ie,in,iw,el
    call t_startf('edge_unpack')
    !$acc parallel loop gang vector collapse(2) present(v,getmapP,edgebuf)
    do el = 1 , nelemd
      do k = 1 , vlyr
        do i = 1 , np
          v(i ,1 ,k,tl,el) = v(i ,1 ,k,tl,el) + edgebuf(kptr+k,getmapP(south,el)+i)
          v(np,i ,k,tl,el) = v(np,i ,k,tl,el) + edgebuf(kptr+k,getmapP(east ,el)+i)
          v(i ,np,k,tl,el) = v(i ,np,k,tl,el) + edgebuf(kptr+k,getmapP(north,el)+i)
          v(1 ,i ,k,tl,el) = v(1 ,i ,k,tl,el) + edgebuf(kptr+k,getmapP(west ,el)+i)
        enddo
        do i = 1 , max_corner_elem
          ll = swest+0*max_corner_elem+i-1; if(getmapP(ll,el) /= -1) v(1  ,1 ,k,tl,el) = v(1 ,1 ,k,tl,el) + edgebuf(kptr+k,getmapP(ll,el)+1)
          ll = swest+1*max_corner_elem+i-1; if(getmapP(ll,el) /= -1) v(np ,1 ,k,tl,el) = v(np,1 ,k,tl,el) + edgebuf(kptr+k,getmapP(ll,el)+1)
          ll = swest+2*max_corner_elem+i-1; if(getmapP(ll,el) /= -1) v(1  ,np,k,tl,el) = v(1 ,np,k,tl,el) + edgebuf(kptr+k,getmapP(ll,el)+1)
          ll = swest+3*max_corner_elem+i-1; if(getmapP(ll,el) /= -1) v(np ,np,k,tl,el) = v(np,np,k,tl,el) + edgebuf(kptr+k,getmapP(ll,el)+1)
        enddo
      enddo
    enddo
    call t_stopf('edge_unpack')
  end subroutine edgeVunpack

  subroutine edgeVunpackMin(edgebuf,nlyr,v,vlyr,kptr,getmapP,tdim,tl)
    use dimensions_mod, only : np, max_corner_elem
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    use perf_mod, only: t_startf, t_stopf
    real(kind=real_kind)  , intent(in   ) :: edgebuf(nlyr,nbuf)
    integer               , intent(in   ) :: nlyr
    integer               , intent(in   ) :: vlyr
    real(kind=real_kind)  , intent(inout) :: v(np,np,vlyr,tdim,nelemd)
    integer               , intent(in   ) :: kptr
    integer(kind=int_kind), intent(in   ) :: getmapP(max_neigh_edges,nelemd)
    integer                ,intent(in   ) :: tdim,tl
    ! Local
    integer :: i,k,ll,is,ie,in,iw,el
    call t_startf('edge_unpack_min')
    !$acc parallel loop gang vector collapse(2) present(v,getmapP,edgebuf)
    do el = 1 , nelemd
      do k = 1 , vlyr
        do i = 1 , np
          v(i ,1 ,k,tl,el) = min(v(i ,1 ,k,tl,el) , edgebuf(kptr+k,getmapP(south,el)+i))
          v(np,i ,k,tl,el) = min(v(np,i ,k,tl,el) , edgebuf(kptr+k,getmapP(east ,el)+i))
          v(i ,np,k,tl,el) = min(v(i ,np,k,tl,el) , edgebuf(kptr+k,getmapP(north,el)+i))
          v(1 ,i ,k,tl,el) = min(v(1 ,i ,k,tl,el) , edgebuf(kptr+k,getmapP(west ,el)+i))
        enddo
        do i = 1 , max_corner_elem
          ll = swest+0*max_corner_elem+i-1; if(getmapP(ll,el) /= -1) v(1  ,1 ,k,tl,el) = min(v(1 ,1 ,k,tl,el) , edgebuf(kptr+k,getmapP(ll,el)+1))
          ll = swest+1*max_corner_elem+i-1; if(getmapP(ll,el) /= -1) v(np ,1 ,k,tl,el) = min(v(np,1 ,k,tl,el) , edgebuf(kptr+k,getmapP(ll,el)+1))
          ll = swest+2*max_corner_elem+i-1; if(getmapP(ll,el) /= -1) v(1  ,np,k,tl,el) = min(v(1 ,np,k,tl,el) , edgebuf(kptr+k,getmapP(ll,el)+1))
          ll = swest+3*max_corner_elem+i-1; if(getmapP(ll,el) /= -1) v(np ,np,k,tl,el) = min(v(np,np,k,tl,el) , edgebuf(kptr+k,getmapP(ll,el)+1))
        enddo
      enddo
    enddo
    call t_stopf('edge_unpack_min')
  end subroutine edgeVunpackMin

  subroutine edgeVunpackMax(edgebuf,nlyr,v,vlyr,kptr,getmapP,tdim,tl)
    use dimensions_mod, only : np, max_corner_elem
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    use perf_mod, only: t_startf, t_stopf
    real(kind=real_kind)  , intent(in   ) :: edgebuf(nlyr,nbuf)
    integer               , intent(in   ) :: nlyr
    integer               , intent(in   ) :: vlyr
    real(kind=real_kind)  , intent(inout) :: v(np,np,vlyr,tdim,nelemd)
    integer               , intent(in   ) :: kptr
    integer(kind=int_kind), intent(in   ) :: getmapP(max_neigh_edges,nelemd)
    integer                ,intent(in   ) :: tdim,tl
    ! Local
    integer :: i,k,ll,is,ie,in,iw,el
    call t_startf('edge_unpack_max')
    !$acc parallel loop gang vector collapse(2) present(v,getmapP,edgebuf)
    do el = 1 , nelemd
      do k = 1 , vlyr
        do i = 1 , np
          v(i ,1 ,k,tl,el) = max(v(i ,1 ,k,tl,el) , edgebuf(kptr+k,getmapP(south,el)+i))
          v(np,i ,k,tl,el) = max(v(np,i ,k,tl,el) , edgebuf(kptr+k,getmapP(east ,el)+i))
          v(i ,np,k,tl,el) = max(v(i ,np,k,tl,el) , edgebuf(kptr+k,getmapP(north,el)+i))
          v(1 ,i ,k,tl,el) = max(v(1 ,i ,k,tl,el) , edgebuf(kptr+k,getmapP(west ,el)+i))
        enddo
        do i = 1 , max_corner_elem
          ll = swest+0*max_corner_elem+i-1; if(getmapP(ll,el) /= -1) v(1  ,1 ,k,tl,el) = max(v(1 ,1 ,k,tl,el) , edgebuf(kptr+k,getmapP(ll,el)+1))
          ll = swest+1*max_corner_elem+i-1; if(getmapP(ll,el) /= -1) v(np ,1 ,k,tl,el) = max(v(np,1 ,k,tl,el) , edgebuf(kptr+k,getmapP(ll,el)+1))
          ll = swest+2*max_corner_elem+i-1; if(getmapP(ll,el) /= -1) v(1  ,np,k,tl,el) = max(v(1 ,np,k,tl,el) , edgebuf(kptr+k,getmapP(ll,el)+1))
          ll = swest+3*max_corner_elem+i-1; if(getmapP(ll,el) /= -1) v(np ,np,k,tl,el) = max(v(np,np,k,tl,el) , edgebuf(kptr+k,getmapP(ll,el)+1))
        enddo
      enddo
    enddo
    call t_stopf('edge_unpack_max')
  end subroutine edgeVunpackMax

  subroutine bndry_exchangeV(hybrid,buffer)
    use hybrid_mod, only : hybrid_t
    use kinds, only : log_kind
    use edge_mod, only : Edgebuffer_t
    use schedule_mod, only : schedule_t, cycle_t, schedule
    use dimensions_mod, only: nelemd, np
    use parallel_mod, only : abortmp, status, srequest, rrequest, mpireal_t, mpiinteger_t, mpi_success
    implicit none
    type (hybrid_t)           :: hybrid
    type (EdgeBuffer_t)       :: buffer
    type (Schedule_t),pointer :: pSchedule
    type (Cycle_t),pointer    :: pCycle
    integer                   :: dest,length,tag
    integer                   :: icycle,ierr
    integer                   :: iptr,source,nlyr
    integer                   :: nSendCycles,nRecvCycles
    integer                   :: errorcode,errorlen
    character*(80) errorstring
    integer        :: i
    !$OMP BARRIER
    if(hybrid%ithr == 0) then 
#ifdef _PREDICT
      pSchedule => Schedule(iam)
#else
      pSchedule => Schedule(1)
#endif
      nlyr = buffer%nlyr
      nSendCycles = pSchedule%nSendCycles
      nRecvCycles = pSchedule%nRecvCycles
      !==================================================
      !  Post the Receives 
      !==================================================
      do icycle=1,nRecvCycles
        pCycle => pSchedule%RecvCycle(icycle)
        source =  pCycle%source - 1
        length =  nlyr * pCycle%lengthP
        tag    =  pCycle%tag
        iptr   =  pCycle%ptrP
        call MPI_Irecv(buffer%receive(1,iptr),length,MPIreal_t,source,tag,hybrid%par%comm,Rrequest(icycle),ierr)
        if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,'bndry_exchangeV: Error after call to MPI_Irecv: ',errorstring
        endif
      enddo    ! icycle
      !==================================================
      !  Fire off the sends
      !==================================================
      do icycle=1,nSendCycles
        pCycle => pSchedule%SendCycle(icycle)
        iptr   =  pCycle%ptrP
        if (pCycle%lengthP > 0) then
          call update_host(buffer%buf(1,iptr),nlyr*pCycle%lengthP)
        endif
      enddo
      !$acc wait
      do icycle=1,nSendCycles
        pCycle => pSchedule%SendCycle(icycle)
        dest   =  pCycle%dest - 1
        length =  nlyr * pCycle%lengthP
        tag    =  pCycle%tag
        iptr   =  pCycle%ptrP
        call MPI_Isend(buffer%buf(1,iptr),length,MPIreal_t,dest,tag,hybrid%par%comm,Srequest(icycle),ierr)
        if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,'bndry_exchangeV: Error after call to MPI_Isend: ',errorstring
        endif
      enddo    ! icycle
      !==================================================
      !  Wait for all the receives to complete
      !==================================================
      call MPI_Waitall(nSendCycles,Srequest,status,ierr)
      call MPI_Waitall(nRecvCycles,Rrequest,status,ierr)
      do icycle=1,nRecvCycles
        pCycle => pSchedule%RecvCycle(icycle)
        iptr   =  pCycle%ptrP
        do i=0,pCycle%lengthP-1
          buffer%buf(1:nlyr,iptr+i) = buffer%receive(1:nlyr,iptr+i)
        enddo
        call update_device(buffer%buf(1,iptr),nlyr*pCycle%lengthP)
      enddo   ! icycle
    endif  ! if (hybrid%ithr == 0)
    !$OMP BARRIER
  end subroutine bndry_exchangeV

  subroutine update_host(dat,len)
    implicit none
    real(kind=real_kind) :: dat(len)
    integer              :: len
    !$acc update host(dat(1:len))
  end subroutine update_host

  subroutine update_device(dat,len)
    implicit none
    real(kind=real_kind) :: dat(len)
    integer              :: len
    !$acc update device(dat(1:len))
  end subroutine update_device

  subroutine limiter2d_zero(Qdp,tdim,tl)
    ! mass conserving zero limiter (2D only).  to be called just before DSS
    !
    ! this routine is called inside a DSS loop, and so Q had already
    ! been multiplied by the mass matrix.  Thus dont include the mass
    ! matrix when computing the mass = integral of Q over the element
    !
    ! ps is only used when advecting Q instead of Qdp
    ! so ps should be at one timelevel behind Q
    implicit none
    real (kind=real_kind), intent(inout) :: Qdp(np,np,nlev,qsize,tdim,nelemd)
    integer              , intent(in   ) :: tdim
    integer              , intent(in   ) :: tl

    ! local
    real (kind=real_kind) :: mass,mass_new
    real (kind=real_kind) :: qtmp(np,np)
    integer i,j,k,q,ie

    !$acc parallel loop gang vector collapse(3) private(qtmp)
    do ie = 1 , nelemd
      do q = 1 , qsize
        do k = nlev , 1 , -1
          !$acc cache(qtmp)
          qtmp = Qdp(:,:,k,q,tl,ie)
          mass = 0
          do j = 1 , np
            do i = 1 , np
              mass = mass + qtmp(i,j)
            enddo
          enddo

          ! negative mass.  so reduce all postive values to zero 
          ! then increase negative values as much as possible
          if ( mass < 0 ) qtmp = -qtmp 
          mass_new = 0
          do j = 1 , np
            do i = 1 , np
              if ( qtmp(i,j) < 0 ) then
                qtmp(i,j) = 0
              else
                mass_new = mass_new + qtmp(i,j)
              endif
            enddo
          enddo

          ! now scale the all positive values to restore mass
          if ( mass_new > 0 ) qtmp = qtmp * abs(mass) / mass_new
          if ( mass     < 0 ) qtmp = -qtmp 
          Qdp(:,:,k,q,tl,ie) = qtmp
        enddo
      enddo
    enddo
  end subroutine limiter2d_zero

  subroutine limiter_optim_iter_full(ptens,elem,minp,maxp,dpmass)
    use element_mod, only: element_t
    use kinds         , only : real_kind
    use dimensions_mod, only : np, np, nlev
    implicit none
    !THIS IS A NEW VERSION OF LIM8, POTENTIALLY FASTER BECAUSE INCORPORATES KNOWLEDGE FROM
    !PREVIOUS ITERATIONS
    
    !The idea here is the following: We need to find a grid field which is closest
    !to the initial field (in terms of weighted sum), but satisfies the min/max constraints.
    !So, first we find values which do not satisfy constraints and bring these values
    !to a closest constraint. This way we introduce some mass change (addmass),
    !so, we redistribute addmass in the way that l2 error is smallest. 
    !This redistribution might violate constraints thus, we do a few iterations. 
    real (kind=real_kind), intent(inout) :: ptens (np*np,nlev,qsize,nelemd)
    type(element_t)      , intent(in   ) :: elem  (:)
    real (kind=real_kind), intent(inout) :: minp  (      nlev,qsize,nelemd)
    real (kind=real_kind), intent(inout) :: maxp  (      nlev,qsize,nelemd)
    real (kind=real_kind), intent(in   ) :: dpmass(np*np,nlev      ,nelemd)
 
    integer  k1, k, i, j, iter, i1, i2, q, ie
    integer :: whois_neg(np*np), whois_pos(np*np), neg_counter, pos_counter
    real (kind=real_kind) :: addmass, weightssum, mass
    real (kind=real_kind) :: x(np*np),c(np*np)
    real (kind=real_kind) :: al_neg(np*np), al_pos(np*np), howmuch
    real (kind=real_kind) :: tol_limiter = 1e-15
    integer, parameter :: maxiter = 5

    !$acc parallel loop gang vector collapse(3) present(ptens,elem(:),minp,maxp,dpmass) private(c,x,whois_neg,whois_pos,al_neg,al_pos)
    do ie = 1 , nelemd
      do q = 1 , qsize
        do k = 1 , nlev
          !$acc cache(c,x)
          do k1 = 1 , np*np
            i = modulo(k1-1,np)+1
            j = (k1-1)/np+1
            c(k1) = elem(ie)%spheremp(i,j) * dpmass(k1,k,ie)
            x(k1) = ptens(k1,k,q,ie) / dpmass(k1,k,ie)
          enddo

          mass = sum(c*x)

          ! relax constraints to ensure limiter has a solution:
          ! This is only needed if runnign with the SSP CFL>1 or 
          ! due to roundoff errors
          if( (mass / sum(c)) < minp(k,q,ie) ) then
            minp(k,q,ie) = mass / sum(c)
          endif
          if( (mass / sum(c)) > maxp(k,q,ie) ) then
            maxp(k,q,ie) = mass / sum(c)
          endif

          addmass = 0.0d0
          pos_counter = 0;
          neg_counter = 0;
          
          ! apply constraints, compute change in mass caused by constraints 
          do k1 = 1 , np*np
            if ( ( x(k1) >= maxp(k,q,ie) ) ) then
              addmass = addmass + ( x(k1) - maxp(k,q,ie) ) * c(k1)
              x(k1) = maxp(k,q,ie)
              whois_pos(k1) = -1
            else
              pos_counter = pos_counter+1;
              whois_pos(pos_counter) = k1;
            endif
            if ( ( x(k1) <= minp(k,q,ie) ) ) then
              addmass = addmass - ( minp(k,q,ie) - x(k1) ) * c(k1)
              x(k1) = minp(k,q,ie)
              whois_neg(k1) = -1
            else
              neg_counter = neg_counter+1;
              whois_neg(neg_counter) = k1;
            endif
          enddo
          
          ! iterate to find field that satifies constraints and is l2-norm closest to original 
          weightssum = 0.0d0
          if ( addmass > 0 ) then
            do i2 = 1 , maxIter
              weightssum = 0.0
              do k1 = 1 , pos_counter
                i1 = whois_pos(k1)
                weightssum = weightssum + c(i1)
                al_pos(i1) = maxp(k,q,ie) - x(i1)
              enddo
              
              if( ( pos_counter > 0 ) .and. ( addmass > tol_limiter * abs(mass) ) ) then
                do k1 = 1 , pos_counter
                  i1 = whois_pos(k1)
                  howmuch = addmass / weightssum
                  if ( howmuch > al_pos(i1) ) then
                    howmuch = al_pos(i1)
                    whois_pos(k1) = -1
                  endif
                  addmass = addmass - howmuch * c(i1)
                  weightssum = weightssum - c(i1)
                  x(i1) = x(i1) + howmuch
                enddo
                !now sort whois_pos and get a new number for pos_counter
                !here neg_counter and whois_neg serve as temp vars
                neg_counter = pos_counter
                whois_neg = whois_pos
                whois_pos = -1
                pos_counter = 0
                do k1 = 1 , neg_counter
                  if ( whois_neg(k1) .ne. -1 ) then
                    pos_counter = pos_counter+1
                    whois_pos(pos_counter) = whois_neg(k1)
                  endif
                enddo
              else
                exit
              endif
            enddo
          else
             do i2 = 1 , maxIter
               weightssum = 0.0
               do k1 = 1 , neg_counter
                 i1 = whois_neg(k1)
                 weightssum = weightssum + c(i1)
                 al_neg(i1) = x(i1) - minp(k,q,ie)
               enddo
               
               if ( ( neg_counter > 0 ) .and. ( (-addmass) > tol_limiter * abs(mass) ) ) then
                 do k1 = 1 , neg_counter
                   i1 = whois_neg(k1)
                   howmuch = -addmass / weightssum
                   if ( howmuch > al_neg(i1) ) then
                     howmuch = al_neg(i1)
                     whois_neg(k1) = -1
                   endif
                   addmass = addmass + howmuch * c(i1)
                   weightssum = weightssum - c(i1)
                   x(i1) = x(i1) - howmuch
                 enddo
                 !now sort whois_pos and get a new number for pos_counter
                 !here pos_counter and whois_pos serve as temp vars
                 pos_counter = neg_counter
                 whois_pos = whois_neg
                 whois_neg = -1
                 neg_counter = 0
                 do k1 = 1 , pos_counter
                   if ( whois_pos(k1) .ne. -1 ) then
                     neg_counter = neg_counter+1
                     whois_neg(neg_counter) = whois_pos(k1)
                   endif
                 enddo
               else
                 exit
               endif
             enddo
          endif
          ptens(:,k,q,ie) = x * dpmass(:,k,ie)
        enddo
      enddo
    enddo
  end subroutine limiter_optim_iter_full

  subroutine precompute_divdp_openacc( elem , hybrid , deriv , dt , nets , nete , n0_qdp )
    use element_mod   , only: element_t, derived_vn0, derived_divdp, derived_divdp_proj
    use hybrid_mod    , only: hybrid_t
    use derivative_mod, only: derivative_t
    use edge_mod      , only: edgeVpack, edgeVunpack
    use bndry_mod     , only: bndry_exchangeV
    use control_mod   , only: limiter_option
    implicit none
    type(element_t)      , intent(inout) :: elem(:)
    type (hybrid_t)      , intent(in   ) :: hybrid
    type (derivative_t)  , intent(in   ) :: deriv
    real(kind=real_kind) , intent(in   ) :: dt
    integer              , intent(in   ) :: nets , nete , n0_qdp
    integer :: ie , k
    real(kind=real_kind), pointer, dimension(:,:,:) :: DSSvar
    call copy_qdp1_h2d(elem,n0_qdp)
    !$omp barrier
    !$omp master
    !$acc update device(derived_vn0)
    call divergence_sphere(derived_vn0,deriv,elem,derived_divdp,nlev)
    call copy_arr(derived_divdp_proj,derived_divdp,product(shape(derived_divdp)))
    !$acc update host(derived_divdp,derived_divdp_proj)
    !$omp end master
    !$omp barrier
    do ie = nets , nete
      ! note: eta_dot_dpdn is actually dimension nlev+1, but nlev+1 data is
      ! all zero so we only have to DSS 1:nlev
      do k = 1 , nlev
        elem(ie)%derived%eta_dot_dpdn(:,:,k) = elem(ie)%spheremp(:,:) * elem(ie)%derived%eta_dot_dpdn(:,:,k) 
        elem(ie)%derived%omega_p(:,:,k)      = elem(ie)%spheremp(:,:) * elem(ie)%derived%omega_p(:,:,k)      
        elem(ie)%derived%divdp_proj(:,:,k)   = elem(ie)%spheremp(:,:) * elem(ie)%derived%divdp_proj(:,:,k)   
      enddo
      call edgeVpack( edgeAdv3 , elem(ie)%derived%eta_dot_dpdn(:,:,1:nlev) , nlev , 0      , elem(ie)%desc )
      call edgeVpack( edgeAdv3 , elem(ie)%derived%omega_p(:,:,1:nlev)      , nlev , nlev   , elem(ie)%desc )
      call edgeVpack( edgeAdv3 , elem(ie)%derived%divdp_proj(:,:,1:nlev)   , nlev , nlev*2 , elem(ie)%desc )
    enddo

    call bndry_exchangeV( hybrid , edgeAdv3   )

    do ie = nets , nete
      call edgeVunpack( edgeAdv3 , elem(ie)%derived%eta_dot_dpdn(:,:,1:nlev) , nlev , 0      , elem(ie)%desc )
      call edgeVunpack( edgeAdv3 , elem(ie)%derived%omega_p(:,:,1:nlev)      , nlev , nlev   , elem(ie)%desc )
      call edgeVunpack( edgeAdv3 , elem(ie)%derived%divdp_proj(:,:,1:nlev)   , nlev , nlev*2 , elem(ie)%desc )
      do k = 1 , nlev
        elem(ie)%derived%eta_dot_dpdn(:,:,k) = elem(ie)%rspheremp(:,:) * elem(ie)%derived%eta_dot_dpdn(:,:,k) 
        elem(ie)%derived%omega_p(:,:,k)      = elem(ie)%rspheremp(:,:) * elem(ie)%derived%omega_p(:,:,k)      
        elem(ie)%derived%divdp_proj(:,:,k)   = elem(ie)%rspheremp(:,:) * elem(ie)%derived%divdp_proj(:,:,k)   
      enddo
    enddo
  end subroutine precompute_divdp_openacc

  subroutine copy_arr(dest,src,len)
    implicit none
    real(kind=real_kind), intent(  out) :: dest(len)
    real(kind=real_kind), intent(in   ) :: src (len)
    integer             , intent(in   ) :: len
    integer :: i
    !$acc parallel loop gang vector present(src,dest)
    do i = 1 , len
      dest(i) = src(i)
    enddo
  end subroutine copy_arr

  subroutine divergence_sphere(v,deriv,elem,div,len)
!   input:  v = velocity in lat-lon coordinates
!   ouput:  div(v)  spherical divergence of v
    use element_mod   , only: element_t
    use derivative_mod, only: derivative_t
    use physical_constants, only: rrearth
    implicit none
    real(kind=real_kind), intent(in   ) :: v(np,np,2,len,nelemd)  ! in lat-lon coordinates
    type(derivative_t)  , intent(in   ) :: deriv
    type(element_t)     , intent(in   ) :: elem(:)
    real(kind=real_kind), intent(  out) :: div(np,np,len,nelemd)
    integer             , intent(in   ) :: len
    ! Local
    integer, parameter :: kchunk = 8
    integer :: i, j, l, k, ie, kc, kk
    real(kind=real_kind) ::  dudx00, dvdy00, gv(np,np,kchunk,2)
    ! convert to contra variant form and multiply by g
    !$acc parallel loop gang collapse(2) private(gv) present(v,deriv,div,elem(:))
    do ie = 1 , nelemd
      do kc = 1 , len/kchunk+1
        !$acc cache(gv)
        !$acc loop vector collapse(3) private(k)
        do kk = 1 , kchunk
          do j = 1 , np
            do i = 1 , np
              k = (kc-1)*kchunk+kk
              if (k <= len) then
                gv(i,j,kk,1)=elem(ie)%metdet(i,j)*(elem(ie)%Dinv(1,1,i,j)*v(i,j,1,k,ie) + elem(ie)%Dinv(1,2,i,j)*v(i,j,2,k,ie))
                gv(i,j,kk,2)=elem(ie)%metdet(i,j)*(elem(ie)%Dinv(2,1,i,j)*v(i,j,1,k,ie) + elem(ie)%Dinv(2,2,i,j)*v(i,j,2,k,ie))
              endif
            enddo
          enddo
        enddo
        ! compute d/dx and d/dy         
        !$acc loop vector collapse(3) private(dudx00,dvdy00,k)
        do kk = 1 , kchunk
          do j = 1 , np
            do i = 1 , np
              k = (kc-1)*kchunk+kk
              if (k <= len) then
                dudx00=0.0d0
                dvdy00=0.0d0
                do l = 1 , np
                  dudx00 = dudx00 + deriv%Dvv(l,i)*gv(l,j,kk,1)
                  dvdy00 = dvdy00 + deriv%Dvv(l,j)*gv(i,l,kk,2)
                enddo
                div(i,j,k,ie)=(dudx00+dvdy00)*(elem(ie)%rmetdet(i,j)*rrearth)
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
  end subroutine divergence_sphere

#endif
end module prim_advection_openacc_mod



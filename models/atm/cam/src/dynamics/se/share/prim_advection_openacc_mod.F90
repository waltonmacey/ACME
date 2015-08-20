!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! OpenACC implementations of some prim_advection_mod routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_advection_openacc_mod
#if USE_OPENACC
  use kinds          , only: real_kind
  use dimensions_mod , only: np,nlevp,nlev,qsize,qsize_d,max_corner_elem,max_neigh_edges,nelemd
  use element_mod    , only: timelevels
  use edge_mod       , only: EdgeBuffer_t
  implicit none
  private
  real(kind=real_kind), allocatable :: qmin(:,:,:), qmax(:,:,:)
  real(kind=real_kind), allocatable :: dp0(:)
  real(kind=real_kind), allocatable :: Qtens_biharmonic(:,:,:,:,:)
  real(kind=real_kind), allocatable :: Qmin_pack(:,:,:,:,:)
  real(kind=real_kind), allocatable :: Qmax_pack(:,:,:,:,:)
  type (EdgeBuffer_t) :: edgeAdv, edgeAdvQ3, edgeAdv_p1, edgeAdvQ2, edgeAdv1
  integer,parameter :: DSSeta = 1
  integer,parameter :: DSSomega = 2
  integer,parameter :: DSSdiv_vdp_ave = 3
  integer,parameter :: DSSno_var = -1

  public :: prim_advection_openacc_init
  public :: precompute_divdp_openacc
  public :: euler_step_openacc

contains

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
    real(kind=real_kind), pointer :: buf_ptr(:) => null()
    real(kind=real_kind), pointer :: receive_ptr(:) => null()
    integer :: k

    ! this might be called with qsize=0
    ! allocate largest one first
    ! Currently this is never freed. If it was, only this first one should
    ! be freed, as only it knows the true size of the buffer.
    call initEdgeBuffer(edgeAdvQ3,max(nlev,qsize*nlev*3), buf_ptr, receive_ptr)  ! Qtens,Qmin, Qmax

    ! remaining edge buffers can share %buf and %receive with edgeAdvQ3
    ! (This is done through the optional 1D pointer arguments.)
    call initEdgeBuffer(edgeAdv1,nlev,buf_ptr,receive_ptr)
    call initEdgeBuffer(edgeAdv,qsize*nlev,buf_ptr,receive_ptr)
    call initEdgeBuffer(edgeAdv_p1,qsize*nlev + nlev,buf_ptr,receive_ptr) 
    call initEdgeBuffer(edgeAdvQ2,qsize*nlev*2,buf_ptr,receive_ptr)  ! Qtens,Qmin, Qmax

    ! Don't actually want these saved, if this is ever called twice.
    nullify(buf_ptr)
    nullify(receive_ptr)

    ! this static array is shared by all threads, so dimension for all threads (nelemd), not nets:nete:
    allocate (qmin(nlev,qsize,nelemd))
    allocate (qmax(nlev,qsize,nelemd))
    allocate (Qtens_biharmonic(np,np,nlev,qsize,nelemd))
    allocate (qmin_pack(np,np,nlev,qsize,nelemd))
    allocate (qmax_pack(np,np,nlev,qsize,nelemd))
    allocate (dp0(nlev))

    do k = 1 , nlev
      dp0(k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
               ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0
    enddo

    !$OMP BARRIER
    if (hybrid%ithr == 0) then
      !$acc enter data pcreate(state_v,state_T,state_dp3d,state_lnps,state_ps_v,state_phis,state_Q,state_Qdp,derived_vn0,derived_vstar,derived_dpdiss_biharmonic,&
      !$acc&                   derived_dpdiss_ave,derived_phi,derived_omega_p,derived_eta_dot_dpdn,derived_grad_lnps,derived_zeta,derived_div,derived_dp,&
      !$acc&                   derived_divdp,derived_divdp_proj,qtens_biharmonic)
      !$acc enter data pcreate(qmin,qmax,qmin_pack,qmax_pack)
      !$acc enter data pcopyin(elem(1:nelemd),deriv,dp0)
      !$acc enter data pcopyin(edgeAdvQ3,edgeAdv1,edgeAdv,edgeAdv_p1,edgeAdvQ2)
      !$acc enter data pcopyin(edgeAdvQ3%buf,edgeAdv1%buf,edgeAdv%buf,edgeAdv_p1%buf,edgeAdvQ2%buf)
      !$acc enter data pcopyin(edgeAdvQ3%receive,edgeAdv1%receive,edgeAdv%receive,edgeAdv_p1%receive,edgeAdvQ2%receive)
    endif
    !$OMP BARRIER
  end subroutine prim_advection_openacc_init

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
  use derivative_mod , only: derivative_t, divergence_sphere, gradient_sphere, vorticity_sphere
  use edge_mod       , only: edgevpack, edgevunpack
  use bndry_mod      , only: bndry_exchangev
  use hybvcoord_mod  , only: hvcoord_t
  use control_mod    , only: limiter_option, nu_p, nu_q
  use perf_mod       , only: t_startf, t_stopf
  use viscosity_mod  , only: biharmonic_wk_scalar_minmax
  use element_mod    , only: derived_dp, derived_divdp_proj, state_qdp, derived_dpdiss_ave
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
  real(kind=real_kind), dimension(np,np       ) :: divdp, dpdiss
  real(kind=real_kind), dimension(np,np,2     ) :: gradQ
  real(kind=real_kind), dimension(np,np,2,nlev) :: Vstar
  real(kind=real_kind), dimension(np,np  ,nlev) :: Qtens
  real(kind=real_kind), dimension(np,np  ,nlev) :: dp,dp_star
  real(kind=real_kind), pointer, dimension(:,:,:) :: DSSvar
  real(kind=real_kind) :: tmp,mintmp,maxtmp
  integer :: ie,q,i,j,k
  integer :: rhs_viss

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
    !$acc update device(derived_dp,derived_divdp_proj,state_qdp)
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
    !$acc update device(qmin,qmax)
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
    !$acc update host(qmin,qmax,qtens_biharmonic)
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
        !$acc update device(qtens_biharmonic,derived_dpdiss_ave)
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
        !$acc update host(qtens_biharmonic)
        !$omp end master
        !$omp barrier
      endif
      call biharmonic_wk_scalar_minmax( elem , qtens_biharmonic(:,:,:,:,nets:nete) , deriv , edgeAdvQ3 , hybrid , nets , nete , qmin(:,:,nets:nete) , qmax(:,:,nets:nete) )
      !$omp barrier
      !$omp master
      !$acc update device(qtens_biharmonic)
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
      !$acc update host(qtens_biharmonic)
      !$omp end master
      !$omp barrier
    endif
  endif  ! compute biharmonic mixing term and qmin/qmax


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   2D Advection step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ie = nets , nete
    ! note: eta_dot_dpdn is actually dimension nlev+1, but nlev+1 data is
    ! all zero so we only have to DSS 1:nlev
    if ( DSSopt == DSSeta         ) DSSvar => elem(ie)%derived%eta_dot_dpdn(:,:,:)
    if ( DSSopt == DSSomega       ) DSSvar => elem(ie)%derived%omega_p(:,:,:)
    if ( DSSopt == DSSdiv_vdp_ave ) DSSvar => elem(ie)%derived%divdp_proj(:,:,:)

    ! Compute velocity used to advance Qdp 
    do k = 1 , nlev    !  Loop index added (AAM)
      ! derived variable divdp_proj() (DSS'd version of divdp) will only be correct on 2nd and 3rd stage
      ! but that's ok because rhs_multiplier=0 on the first stage:
      dp(:,:,k) = elem(ie)%derived%dp(:,:,k) - rhs_multiplier * dt * elem(ie)%derived%divdp_proj(:,:,k) 
      Vstar(:,:,1,k) = elem(ie)%derived%vn0(:,:,1,k) / dp(:,:,k)
      Vstar(:,:,2,k) = elem(ie)%derived%vn0(:,:,2,k) / dp(:,:,k)
    enddo

    ! advance Qdp
    do q = 1 , qsize
      do k = 1 , nlev  !  dp_star used as temporary instead of divdp (AAM)
        ! div( U dp Q), 
        gradQ(:,:,1) = Vstar(:,:,1,k) * elem(ie)%state%Qdp(:,:,k,q,n0_qdp)
        gradQ(:,:,2) = Vstar(:,:,2,k) * elem(ie)%state%Qdp(:,:,k,q,n0_qdp)
        dp_star(:,:,k) = divergence_sphere( gradQ , deriv , elem(ie) )
        Qtens(:,:,k) = elem(ie)%state%Qdp(:,:,k,q,n0_qdp) - dt * dp_star(:,:,k)
        ! optionally add in hyperviscosity computed above:
        if ( rhs_viss /= 0 ) Qtens(:,:,k) = Qtens(:,:,k) + Qtens_biharmonic(:,:,k,q,ie)
      enddo
         
      if ( limiter_option == 8 ) then
        do k = 1 , nlev  ! Loop index added (AAM)
          ! UN-DSS'ed dp at timelevel n0+1:  
          dp_star(:,:,k) = dp(:,:,k) - dt * elem(ie)%derived%divdp(:,:,k)  
          if ( nu_p > 0 .and. rhs_viss /= 0 ) then
            ! add contribution from UN-DSS'ed PS dissipation
!            dpdiss(:,:) = ( hvcoord%hybi(k+1) - hvcoord%hybi(k) ) * elem(ie)%derived%psdiss_biharmonic(:,:)
            dpdiss(:,:) = elem(ie)%derived%dpdiss_biharmonic(:,:,k)
            dp_star(:,:,k) = dp_star(:,:,k) - rhs_viss * dt * nu_q * dpdiss(:,:) / elem(ie)%spheremp(:,:)
          endif
        enddo
        ! apply limiter to Q = Qtens / dp_star 
        call limiter_optim_iter_full( Qtens(:,:,:) , elem(ie)%spheremp(:,:) , qmin(:,q,ie) , &
                                      qmax(:,q,ie) , dp_star(:,:,:) )
      endif

      ! apply mass matrix, overwrite np1 with solution:
      ! dont do this earlier, since we allow np1_qdp == n0_qdp 
      ! and we dont want to overwrite n0_qdp until we are done using it
      do k = 1 , nlev
        elem(ie)%state%Qdp(:,:,k,q,np1_qdp) = elem(ie)%spheremp(:,:) * Qtens(:,:,k) 
      enddo

      if ( limiter_option == 4 ) then
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
        ! sign-preserving limiter, applied after mass matrix
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
        call limiter2d_zero( elem(ie)%state%Qdp(:,:,:,q,np1_qdp) ) 
      endif
    enddo

    if ( DSSopt == DSSno_var ) then
      call edgeVpack(edgeAdv    , elem(ie)%state%Qdp(:,:,:,:,np1_qdp) , nlev*qsize , 0 , elem(ie)%desc )
    else
      call edgeVpack(edgeAdv_p1 , elem(ie)%state%Qdp(:,:,:,:,np1_qdp) , nlev*qsize , 0 , elem(ie)%desc )
      ! also DSS extra field
      do k = 1 , nlev
        DSSvar(:,:,k) = elem(ie)%spheremp(:,:) * DSSvar(:,:,k) 
      enddo
      call edgeVpack( edgeAdv_p1 , DSSvar(:,:,1:nlev) , nlev , nlev*qsize , elem(ie)%desc )
    endif
  enddo

  call t_startf('eus_bexchV')
  if ( DSSopt == DSSno_var ) then
    call bndry_exchangeV( hybrid , edgeAdv    )
  else
    call bndry_exchangeV( hybrid , edgeAdv_p1 )
  endif
  call t_stopf('eus_bexchV')

  do ie = nets , nete
    if ( DSSopt == DSSeta         ) DSSvar => elem(ie)%derived%eta_dot_dpdn(:,:,:)
    if ( DSSopt == DSSomega       ) DSSvar => elem(ie)%derived%omega_p(:,:,:)
    if ( DSSopt == DSSdiv_vdp_ave ) DSSvar => elem(ie)%derived%divdp_proj(:,:,:)

    if ( DSSopt == DSSno_var ) then
      call edgeVunpack( edgeAdv    , elem(ie)%state%Qdp(:,:,:,:,np1_qdp) , nlev*qsize , 0 , elem(ie)%desc )
      do q = 1 , qsize
        do k = 1 , nlev    !  Potential loop inversion (AAM)
          elem(ie)%state%Qdp(:,:,k,q,np1_qdp) = elem(ie)%rspheremp(:,:) * elem(ie)%state%Qdp(:,:,k,q,np1_qdp)
        enddo
      enddo
    else
      call edgeVunpack( edgeAdv_p1 , elem(ie)%state%Qdp(:,:,:,:,np1_qdp) , nlev*qsize , 0 , elem(ie)%desc )
      do q = 1 , qsize
        do k = 1 , nlev    !  Potential loop inversion (AAM)
          elem(ie)%state%Qdp(:,:,k,q,np1_qdp) = elem(ie)%rspheremp(:,:) * elem(ie)%state%Qdp(:,:,k,q,np1_qdp)
        enddo
      enddo
      call edgeVunpack( edgeAdv_p1 , DSSvar(:,:,1:nlev) , nlev , qsize*nlev , elem(ie)%desc )
       
      do k = 1 , nlev
        DSSvar(:,:,k) = DSSvar(:,:,k) * elem(ie)%rspheremp(:,:)
      enddo
    endif
  enddo
  end subroutine euler_step_openacc

  subroutine neighbor_minmax(elem,hybrid,edgeMinMax,nets,nete,min_neigh,max_neigh)
    use hybrid_mod , only: hybrid_t
    use element_mod, only: element_t
    use edge_mod   , only: edgeVunpackMin, edgeVunpackMax
    use perf_mod   , only: t_startf, t_stopf
    use bndry_mod  , only: bndry_exchangeV
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
    !$acc update device(min_neigh,max_neigh)
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
    !$acc update host(qmin_pack,qmax_pack)
    !$omp end master
    !$omp barrier
    call edgeVpack(edgeMinMax,Qmin_pack,nlev*qsize,0         ,elem(:),nets,nete)
    call edgeVpack(edgeMinMax,Qmax_pack,nlev*qsize,nlev*qsize,elem(:),nets,nete)

    call t_startf('nmm_bexchV')
    call bndry_exchangeV(hybrid,edgeMinMax)
    call t_stopf('nmm_bexchV')
       
    do ie=nets,nete
     ! WARNING - edgeVunpackMin/Max take second argument as input/ouput
      call edgeVunpackMin(edgeMinMax,Qmin_pack(:,:,:,:,ie),nlev*qsize,0,elem(ie)%desc)
      call edgeVunpackMax(edgeMinMax,Qmax_pack(:,:,:,:,ie),nlev*qsize,nlev*qsize,elem(ie)%desc)
    enddo
    !$omp barrier
    !$omp master
    !$acc update device(qmin_pack,qmax_pack)
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
    !$acc update host(min_neigh,max_neigh)
    !$omp end master
    !$omp barrier
  end subroutine neighbor_minmax

  subroutine edgeVpack(edge,v,vlyr,kptr,elem,nets,nete)
    use dimensions_mod, only : max_corner_elem
    use control_mod   , only : north, south, east, west, neast, nwest, seast, swest
    use perf_mod      , only : t_startf, t_stopf
    use parallel_mod  , only : haltmp
    use element_mod   , only : Element_t
    use edge_mod      , only : EdgeBuffer_t
    type (EdgeBuffer_t)    ,intent(inout) :: edge
    integer                ,intent(in   ) :: vlyr
    real (kind=real_kind)  ,intent(in   ) :: v(np,np,vlyr,nelemd)
    integer                ,intent(in   ) :: kptr
    type (Element_t)       ,intent(in   ) :: elem(:)
    integer                ,intent(in   ) :: nets,nete
    ! Local variables
    integer :: i,k,ir,ll,is,ie,in,iw,el
    call t_startf('edge_pack')
    if (edge%nlyr < (kptr+vlyr) ) call haltmp('edgeVpack: Buffer overflow: size of the vertical dimension must be increased!')
    !$acc update device(v)
    !$acc parallel loop gang vector collapse(3) present(v,elem(:),edge)
    do el = nets , nete
      do k = 1 , vlyr
        do i = 1 , np
          edge%buf(kptr+k,elem(el)%desc%putmapP(south)+i) = v(i ,1 ,k,el)
          edge%buf(kptr+k,elem(el)%desc%putmapP(east )+i) = v(np,i ,k,el)
          edge%buf(kptr+k,elem(el)%desc%putmapP(north)+i) = v(i ,np,k,el)
          edge%buf(kptr+k,elem(el)%desc%putmapP(west )+i) = v(1 ,i ,k,el)
        enddo
      enddo
    enddo
    !$acc update host(edge%buf)
    do el = nets , nete
      do k = 1 , vlyr
        do i = 1 , np
          ir = np-i+1
          if(elem(el)%desc%reverse(south)) edge%buf(kptr+k,elem(el)%desc%putmapP(south)+ir) = v(i ,1 ,k,el)
          if(elem(el)%desc%reverse(east )) edge%buf(kptr+k,elem(el)%desc%putmapP(east )+ir) = v(np,i ,k,el)
          if(elem(el)%desc%reverse(north)) edge%buf(kptr+k,elem(el)%desc%putmapP(north)+ir) = v(i ,np,k,el)
          if(elem(el)%desc%reverse(west )) edge%buf(kptr+k,elem(el)%desc%putmapP(west )+ir) = v(1 ,i ,k,el)
        enddo
      enddo
    enddo
    do el = nets , nete
      do k = 1 , vlyr
        do i = 1 , max_corner_elem
          ll = swest+0*max_corner_elem+i-1; if (elem(el)%desc%putmapP(ll) /= -1) edge%buf(kptr+k,elem(el)%desc%putmapP(ll)+1) = v(1 ,1 ,k,el)
          ll = swest+1*max_corner_elem+i-1; if (elem(el)%desc%putmapP(ll) /= -1) edge%buf(kptr+k,elem(el)%desc%putmapP(ll)+1) = v(np,1 ,k,el)
          ll = swest+2*max_corner_elem+i-1; if (elem(el)%desc%putmapP(ll) /= -1) edge%buf(kptr+k,elem(el)%desc%putmapP(ll)+1) = v(1 ,np,k,el)
          ll = swest+3*max_corner_elem+i-1; if (elem(el)%desc%putmapP(ll) /= -1) edge%buf(kptr+k,elem(el)%desc%putmapP(ll)+1) = v(np,np,k,el)
        enddo
      enddo
    enddo
    call t_stopf('edge_pack')
  end subroutine edgeVpack

  subroutine limiter2d_zero(Q)
  ! mass conserving zero limiter (2D only).  to be called just before DSS
  !
  ! this routine is called inside a DSS loop, and so Q had already
  ! been multiplied by the mass matrix.  Thus dont include the mass
  ! matrix when computing the mass = integral of Q over the element
  !
  ! ps is only used when advecting Q instead of Qdp
  ! so ps should be at one timelevel behind Q
  implicit none
  real (kind=real_kind), intent(inout) :: Q(np,np,nlev)

  ! local
  real (kind=real_kind) :: dp(np,np)
  real (kind=real_kind) :: mass,mass_new,ml
  integer i,j,k

  do k = nlev , 1 , -1
    mass = 0
    do j = 1 , np
      do i = 1 , np
        !ml = Q(i,j,k)*dp(i,j)*spheremp(i,j)  ! see above
        ml = Q(i,j,k)
        mass = mass + ml
      enddo
    enddo

    ! negative mass.  so reduce all postive values to zero 
    ! then increase negative values as much as possible
    if ( mass < 0 ) Q(:,:,k) = -Q(:,:,k) 
    mass_new = 0
    do j = 1 , np
      do i = 1 , np
        if ( Q(i,j,k) < 0 ) then
          Q(i,j,k) = 0
        else
          ml = Q(i,j,k)
          mass_new = mass_new + ml
        endif
      enddo
    enddo

    ! now scale the all positive values to restore mass
    if ( mass_new > 0 ) Q(:,:,k) = Q(:,:,k) * abs(mass) / mass_new
    if ( mass     < 0 ) Q(:,:,k) = -Q(:,:,k) 
  enddo
  end subroutine limiter2d_zero

  subroutine limiter_optim_iter_full(ptens,sphweights,minp,maxp,dpmass)
    !THIS IS A NEW VERSION OF LIM8, POTENTIALLY FASTER BECAUSE INCORPORATES KNOWLEDGE FROM
    !PREVIOUS ITERATIONS
    
    !The idea here is the following: We need to find a grid field which is closest
    !to the initial field (in terms of weighted sum), but satisfies the min/max constraints.
    !So, first we find values which do not satisfy constraints and bring these values
    !to a closest constraint. This way we introduce some mass change (addmass),
    !so, we redistribute addmass in the way that l2 error is smallest. 
    !This redistribution might violate constraints thus, we do a few iterations. 
    use kinds         , only : real_kind
    use dimensions_mod, only : np, np, nlev
    real (kind=real_kind), dimension(np*np,nlev), intent(inout)            :: ptens
    real (kind=real_kind), dimension(np*np     ), intent(in   )            :: sphweights
    real (kind=real_kind), dimension(      nlev), intent(inout)            :: minp
    real (kind=real_kind), dimension(      nlev), intent(inout)            :: maxp
    real (kind=real_kind), dimension(np*np,nlev), intent(in   ), optional  :: dpmass
 
    real (kind=real_kind), dimension(np*np,nlev) :: weights
    integer  k1, k, i, j, iter, i1, i2
    integer :: whois_neg(np*np), whois_pos(np*np), neg_counter, pos_counter
    real (kind=real_kind) :: addmass, weightssum, mass
    real (kind=real_kind) :: x(np*np),c(np*np)
    real (kind=real_kind) :: al_neg(np*np), al_pos(np*np), howmuch
    real (kind=real_kind) :: tol_limiter = 1e-15
    integer, parameter :: maxiter = 5

    do k = 1 , nlev
      weights(:,k) = sphweights(:) * dpmass(:,k)
      ptens(:,k) = ptens(:,k) / dpmass(:,k)
    enddo

    do k = 1 , nlev
      c = weights(:,k)
      x = ptens(:,k)

      mass = sum(c*x)

      ! relax constraints to ensure limiter has a solution:
      ! This is only needed if runnign with the SSP CFL>1 or 
      ! due to roundoff errors
      if( (mass / sum(c)) < minp(k) ) then
        minp(k) = mass / sum(c)
      endif
      if( (mass / sum(c)) > maxp(k) ) then
        maxp(k) = mass / sum(c)
      endif

      addmass = 0.0d0
      pos_counter = 0;
      neg_counter = 0;
      
      ! apply constraints, compute change in mass caused by constraints 
      do k1 = 1 , np*np
        if ( ( x(k1) >= maxp(k) ) ) then
          addmass = addmass + ( x(k1) - maxp(k) ) * c(k1)
          x(k1) = maxp(k)
          whois_pos(k1) = -1
        else
          pos_counter = pos_counter+1;
          whois_pos(pos_counter) = k1;
        endif
        if ( ( x(k1) <= minp(k) ) ) then
          addmass = addmass - ( minp(k) - x(k1) ) * c(k1)
          x(k1) = minp(k)
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
            al_pos(i1) = maxp(k) - x(i1)
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
             al_neg(i1) = x(i1) - minp(k)
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
      
      ptens(:,k) = x
    enddo
    
    do k = 1 , nlev
      ptens(:,k) = ptens(:,k) * dpmass(:,k)
    enddo
  end subroutine limiter_optim_iter_full

  subroutine precompute_divdp_openacc( elem , hybrid , deriv , dt , nets , nete , n0_qdp )
    use element_mod   , only: element_t, derived_vn0, derived_divdp, derived_divdp_proj
    use hybrid_mod    , only: hybrid_t
    use derivative_mod, only: derivative_t
    implicit none
    type(element_t)      , intent(inout) :: elem(:)
    type (hybrid_t)      , intent(in   ) :: hybrid
    type (derivative_t)  , intent(in   ) :: deriv
    real(kind=real_kind) , intent(in   ) :: dt
    integer              , intent(in   ) :: nets , nete , n0_qdp
    integer :: ie , k
    !$omp barrier
    !$omp master
    !$acc update device(derived_vn0)
    call divergence_sphere(derived_vn0,deriv,elem,derived_divdp,nlev)
    call copy_arr(derived_divdp_proj,derived_divdp,product(shape(derived_divdp)))
    !$acc update host(derived_divdp,derived_divdp_proj)
    !$omp end master
    !$omp barrier
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
              if (k <= nlev) then
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
              if (k <= nlev) then
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



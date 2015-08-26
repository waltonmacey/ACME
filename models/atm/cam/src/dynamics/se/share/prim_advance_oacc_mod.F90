#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
!#define _DBG_ print *,"File:",__FILE__," at ",__LINE__
!#define _DBG_ !DBG
!
#ifdef CUDA_DEBUG
#define _CHECK(line) ierr = cudaThreadSynchronize(); if (ierr .ne. 0) stop line
#else
#define _CHECK(line)
#endif

!
module prim_advance_oacc_mod
  use edge_mod, only : EdgeBuffer_t
  use kinds, only : real_kind, iulog
  use perf_mod, only: t_startf, t_stopf, t_barrierf ! _EXTERNAL
  use parallel_mod, only : abortmp
  implicit none
  private
  save
   public :: compute_and_apply_oacc_rhs
   public :: prim_advance_oacc_init, prim_advance_oacc_exit
!  public :: prim_advance_exp, prim_advance_si, prim_advance_init, preq_robert3,&
!       applyCAMforcing_dynamics, applyCAMforcing, smooth_phis, overwrite_SEdensity

  !device arrays:
 ! real (kind=real_kind),allocatable, dimension(:,:)      :: ps_oacc         ! surface pressure for current tiime level
  real (kind=real_kind),allocatable, dimension(:,:,:,:)  :: phi_oacc

contains

 subroutine prim_advance_oacc_init(elem, hvcoord, deriv,nets,nete)
  use kinds, only : real_kind
  use dimensions_mod, only : np, nc, nlev, ntrac, max_corner_elem
  use element_mod, only : element_t !uncomment after merge , PrintElem
  use derivative_mod, only : derivative_t, divergence_sphere, gradient_sphere, vorticity_sphere
  use hybvcoord_mod, only : hvcoord_t

 implicit none
 
    type(element_t)      , intent(in) :: elem(:)
    type (hvcoord_t)     , intent(in) :: hvcoord
    type (derivative_t)  , intent(in) :: deriv
    integer              , intent(in) :: nets,nete
    integer :: i, j, k, l, ie, nelem
    integer :: ierr 

   nelem = nete-nets

   allocate( phi_oacc  (np, np, nlev  , nete+1 ) , stat = ierr ); _CHECK(__LINE__)

   do ie=nets,nete
     do i=1,np
       do j=1,np
         do k=1,nlev
           phi_oacc(i,j,k,ie)=elem(ie)%derived%phi(i,j,k)
         enddo
       enddo
     enddo
   enddo


   !$acc enter data pcopyin( phi_oacc, hvcoord%hyai, hvcoord%hybi, hvcoord%hyam, hvcoord%hybm, &
   !$acc&                    elem(nets:nete))

 
 end subroutine prim_advance_oacc_init


 subroutine prim_advance_oacc_exit( )
 implicit none

   !$acc exit data
 end subroutine prim_advance_oacc_exit
  
  subroutine compute_and_apply_oacc_rhs(np1,nm1,n0,qn0,dt2,elem,hvcoord,hybrid,&
       deriv,nets,nete,compute_diagnostics,eta_ave_w, edge3p1)
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
  !
  ! ===================================
  use kinds, only : real_kind
  use dimensions_mod, only : np, nc, nlev, ntrac, max_corner_elem
  use hybrid_mod, only : hybrid_t
  use element_mod, only : element_t !uncomment after merge , PrintElem
  use derivative_mod, only : derivative_t, divergence_sphere, gradient_sphere, vorticity_sphere
!uncomment after merge  use derivative_mod, only : subcell_div_fluxes, subcell_dss_fluxes
  use edge_mod, only : edgevpack, edgevunpack, edgeDGVunpack
!uncomment after merge
!  use edgetype_mod, only : edgedescriptor_t
  use bndry_mod, only : bndry_exchangev
  use control_mod, only : moisture, qsplit, use_cpstar, rsplit, swest
  use hybvcoord_mod, only : hvcoord_t

  use physical_constants, only : cp, cpwater_vapor, Rgas, kappa, rrearth, Rwater_vapor
  use physics_mod, only : virtual_specific_heat, virtual_temperature
  use prim_si_mod, only : preq_vertadv, preq_omega_ps, preq_hydrostatic
#if ( defined CAM ) 
  use control_mod, only: se_met_nudge_u, se_met_nudge_p, se_met_nudge_t, se_met_tevolve
#endif

!uncomment after merge  use time_mod, only : tevolve 

  implicit none
  integer, intent(in) :: np1,nm1,n0,qn0,nets,nete
  real*8, intent(in) :: dt2
  logical, intent(in)  :: compute_diagnostics

  type (hvcoord_t)     , intent(in) :: hvcoord
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv
  type (EdgeBuffer_t)  , intent(in) :: edge3p1
  real (kind=real_kind) :: eta_ave_w  ! weighting for eta_dot_dpdn mean flux

  ! local
!  real (kind=real_kind), pointer, dimension(:,:)      :: ps         ! surface pressure for current tiime level
!  real (kind=real_kind), pointer, dimension(:,:,:)   :: phi

  real (kind=real_kind), dimension(np,np,nlev)   :: omega_p       
  real (kind=real_kind), dimension(np,np,nlev)   :: T_v         
  real (kind=real_kind), dimension(np,np,nlev)   :: divdp
  real (kind=real_kind), dimension(np,np,nlev+1)   :: eta_dot_dpdn  ! half level vertical velocity on p-grid
  real (kind=real_kind), dimension(np,np)      :: sdot_sum   ! temporary field
  real (kind=real_kind), dimension(np,np,2)    :: vtemp     ! generic gradient storage
  real (kind=real_kind), dimension(np,np,2,nlev):: vdp       !                                            
  real (kind=real_kind), dimension(np,np,2     ):: v         !  
  real (kind=real_kind), dimension(np,np)      :: vgrad_T    ! v.grad(T)
  real (kind=real_kind), dimension(np,np)      :: Ephi       ! kinetic energy + PHI term
  real (kind=real_kind), dimension(np,np,2)      :: grad_ps    ! lat-lon coord version
  real (kind=real_kind), dimension(np,np,2,nlev) :: grad_p
  real (kind=real_kind), dimension(np,np,2,nlev) :: grad_p_m_pmet  ! gradient(p - p_met)
  real (kind=real_kind), dimension(np,np,nlev)   :: vort       ! vorticity
  real (kind=real_kind), dimension(np,np,nlev)   :: p          ! pressure
  real (kind=real_kind), dimension(np,np,nlev)   :: dp         ! delta pressure
  real (kind=real_kind), dimension(np,np,nlev)   :: rdp        ! inverse of delta pressure
  real (kind=real_kind), dimension(np,np,nlev)   :: T_vadv     ! temperature vertical advection
  real (kind=real_kind), dimension(np,np,nlev)   :: vgrad_p    ! v.grad(p)
  real (kind=real_kind), dimension(np,np,nlev+1) :: ph               ! half level pressures on p-grid
  real (kind=real_kind), dimension(np,np,2,nlev) :: v_vadv   ! velocity vertical advection
  real (kind=real_kind), dimension(np+1,np+1,nlev) :: corners 
  real (kind=real_kind), dimension(2,2,2)          :: cflux 
  real (kind=real_kind) ::  kappa_star(np,np,nlev)
  real (kind=real_kind) ::  vtens1(np,np,nlev)
  real (kind=real_kind) ::  vtens2(np,np,nlev)
  real (kind=real_kind) ::  ttens(np,np,nlev)
  real (kind=real_kind) ::  stashdp3d (np,np,nlev)
  real (kind=real_kind) ::  tempdp3d  (np,np)
  real (kind=real_kind) ::  tempflux  (nc,nc,4)
  !uncomment after merge type (EdgeDescriptor_t) :: desc

  real (kind=real_kind) ::  cp2,cp_ratio,E,de,Qt,v1,v2
  real (kind=real_kind) ::  glnps1,glnps2,gpterm
  real(kind=real_kind) ::  dsdx00
  real(kind=real_kind) ::  dsdy00
  real(kind=real_kind) ::  v_tmp1(np,np),v_tmp2(np,np), gv(np,np,2)
  real(kind=real_kind) hkk,hkl, term          ! diagonal term of energy conversion matrix
  real(kind=real_kind), dimension(np,np,nlev) :: phii, phi      ! Geopotential at interfaces
  integer :: i,j,k,l,kptr,ie
  real (kind=real_kind) ::  u_m_umet, v_m_vmet, t_m_tmet


  call t_barrierf('sync_compute_and_apply_rhs', hybrid%par%comm) 
   !uncomment after merge type + comment line before   call t_adj_detailf(+1)

  call t_startf('compute_and_apply_oacc_rhs')
 !$OMP BARRIER
 if (hybrid%ithr == 0) then


!  !$acc update device( elem(nets:nete))
! !$acc parallel loop gang vector private(dsdx00, dsdy00) &
! !$acc&               present(phi_oacc) &
! !$acc&               pcopyin(np1, nm1, qn0,dt2, nets,nete, &
! !$acc&                       compute_diagnostics,eta_ave_w)

  do ie=nets,nete

     ! ==================================================
     ! compute pressure (p) on half levels from ps 
     ! using the hybrid coordinates relationship, i.e.
     ! e.g. equation (3.a.92) of the CCM-2 description, 
     ! (NCAR/TN-382+STR), June 1993, p. 24.
     ! ==================================================
     ! vertically eulerian only needs grad(ps)
     if (rsplit==0) then 
          !! inlined gradient_sphere below: grad_ps = gradient_sphere(elem(ie)%state%ps_v(:,:,n0),deriv,elem(ie)%Dinv)
          do j=1,np
            do l=1,np
             dsdx00=0.0d0
             dsdy00=0.0d0
              do i=1,np
               dsdx00 = dsdx00 + deriv%Dvv(i,l  )*elem(ie)%state%ps_v(i,j,n0)
               dsdy00 = dsdy00 + deriv%Dvv(i,l  )*elem(ie)%state%ps_v(j,i,n0)
              enddo
              v_tmp1(l  ,j  ) = dsdx00*rrearth
              v_tmp2(j  ,l  ) = dsdy00*rrearth
             enddo
           enddo
           do j=1,np
            do i=1,np
             grad_ps(i,j,1)=elem(ie)%Dinv(1,1,i,j)*v_tmp1(i,j) + elem(ie)%Dinv(2,1,i,j)*v_tmp2(i,j)
             grad_ps(i,j,2)=elem(ie)%Dinv(1,2,i,j)*v_tmp1(i,j) + elem(ie)%Dinv(2,2,i,j)*v_tmp2(i,j)
            enddo
           enddo
       endif

     ! ============================
     ! compute p and delta p
     ! ============================
!     do k=1,nlev+1
!       ph(:,:,k)   = hvcoord%hyai(k)*hvcoord%ps0 + hvcoord%hybi(k)*elem(ie)%state%ps_v(:,:,n0)
!     end do
    
     do k=1,nlev

      do j = 1 , np
        do i = 1 , np
          if (rsplit==0) then
            dp(i,j,k) = (hvcoord%hyai(k+1)*hvcoord%ps0 + hvcoord%hybi(k+1)*elem(ie)%state%ps_v(i,j,n0)) &
                        - (hvcoord%hyai(k)*hvcoord%ps0 + hvcoord%hybi(k)*elem(ie)%state%ps_v(i,j,n0))
            p(i,j,k)  = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem(ie)%state%ps_v(i,j,n0)
            grad_p(i,j,1,k) = hvcoord%hybm(k)*grad_ps(i,j,1)
            grad_p(i,j,2,k) = hvcoord%hybm(k)*grad_ps(i,j,2)
           else
           ! vertically lagrangian code: we advect dp3d instead of ps_v
           ! we also need grad(p) at all levels (not just grad(ps))
           !p(k)= hyam(k)*ps0 + hybm(k)*ps
           !    = .5*(hyai(k+1)+hyai(k))*ps0 + .5*(hybi(k+1)+hybi(k))*ps
           !    = .5*(ph(k+1) + ph(k) )  = ph(k) + dp(k)/2
           !
           ! p(k+1)-p(k) = ph(k+1)-ph(k) + (dp(k+1)-dp(k))/2
           !             = dp(k) + (dp(k+1)-dp(k))/2 = (dp(k+1)+dp(k))/2
             dp(i,j,k) = elem(ie)%state%dp3d(i,j,k,n0)
             if (k==1) then
                p(i,j,k)=hvcoord%hyai(k)*hvcoord%ps0 + dp(i,j,k)/2
             else
                p(i,j,k)=p(i,j,k-1) + dp(i,j,k-1)/2 + dp(i,j,k)/2
             endif
           endif
          enddo
        enddo

       if (rsplit/=0) then
          !! inline gradient_sphere here: grad_p(:,:,:,k) = gradient_sphere(p(:,:,k),deriv,elem(ie)%Dinv)
          ! ============================
          ! compute gradien_sphere 
          ! ============================
          do j=1,np
            do l=1,np
             dsdx00=0.0d0
             dsdy00=0.0d0
              do i=1,np
               dsdx00 = dsdx00 + deriv%Dvv(i,l  )*p(i,j,k)
               dsdy00 = dsdy00 + deriv%Dvv(i,l  )*p(j,i,k)
              enddo
              v_tmp1(l  ,j  ) = dsdx00*rrearth
              v_tmp2(j  ,l  ) = dsdy00*rrearth
             enddo
           enddo
           do j=1,np
            do i=1,np
             grad_p(i,j,1,k)=elem(ie)%Dinv(1,1,i,j)*v_tmp1(i,j) + elem(ie)%Dinv(2,1,i,j)*v_tmp2(i,j)
             grad_p(i,j,2,k)=elem(ie)%Dinv(1,2,i,j)*v_tmp1(i,j) + elem(ie)%Dinv(2,2,i,j)*v_tmp2(i,j)
            enddo
           enddo

        endif

        do j=1,np
         do i=1,np
           rdp(i,j,k) = 1.0D0/dp(i,j,k)
         enddo
        enddo

        ! ============================
        ! compute vgrad_lnps 
        ! ============================
        do j=1,np
           do i=1,np
              v1 = elem(ie)%state%v(i,j,1,k,n0)
              v2 = elem(ie)%state%v(i,j,2,k,n0)
!              vgrad_p(i,j,k) = &
!                   hvcoord%hybm(k)*(v1*grad_ps(i,j,1) + v2*grad_ps(i,j,2)) 
              vgrad_p(i,j,k) = (v1*grad_p(i,j,1,k) + v2*grad_p(i,j,2,k)) 
              vdp(i,j,1,k) = v1*dp(i,j,k)
              vdp(i,j,2,k) = v2*dp(i,j,k)
           end do
        end do
#if ( defined CAM ) 
      ! ============================
      ! compute grad(P-P_met) 
      ! ============================
      if (se_met_nudge_p.gt.0.D0) then 
         grad_p_m_pmet(:,:,:,k) = &
                grad_p(:,:,:,k) - &
                hvcoord%hybm(k) * &  
                gradient_sphere( elem(ie)%derived%ps_met(:,:)+tevolve*elem(ie)%derived%dpsdt_met(:,:), &
                                deriv,elem(ie)%Dinv)
      endif

#endif      
        ! ================================
        ! Accumulate mean Vel_rho flux in vn0
        ! ================================
        do j=1,np
           do i=1,np
              elem(ie)%derived%vn0(i,j,1,k)=elem(ie)%derived%vn0(i,j,1,k)+eta_ave_w*vdp(i,j,1,k)
              elem(ie)%derived%vn0(i,j,2,k)=elem(ie)%derived%vn0(i,j,2,k)+eta_ave_w*vdp(i,j,2,k)
           enddo
        enddo

        ! =========================================
        !
        ! Compute relative vorticity and divergence
        !
        ! =========================================
        ! inline divdp(:,:,k)=divergence_sphere(vtemp,deriv,elem(ie))

        ! convert to contra variant form and multiply by g
        do j=1,np
           do i=1,np
             gv(i,j,1)=elem(ie)%metdet(i,j)*(elem(ie)%Dinv(1,1,i,j)*vdp(i,j,1,k) + elem(ie)%Dinv(1,2,i,j)*vdp(i,j,2,k))
             gv(i,j,2)=elem(ie)%metdet(i,j)*(elem(ie)%Dinv(2,1,i,j)*vdp(i,j,1,k) + elem(ie)%Dinv(2,2,i,j)*vdp(i,j,2,k))
           enddo
         enddo
         ! compute d/dx and d/dy         
         do j=1,np
            do l=1,np
               dsdx00=0.0d0
               dsdy00=0.0d0
               do i=1,np
                  dsdx00 = dsdx00 + deriv%Dvv(i,l  )*gv(i,j  ,1)
                  dsdy00 = dsdy00 + deriv%Dvv(i,l  )*gv(j  ,i,2)
               end do
               divdp(l,j,k) = dsdx00
               v_tmp1(j,l)  = dsdy00
             end do
          end do
          do j=1,np
             do i=1,np
               divdp(i,j,k)=(divdp(i,j,k)+v_tmp1(i,j))*(elem(ie)%rmetdet(i,j)*rrearth)
             end do
          end do

        !inline  vort(:,:,k)=vorticity_sphere(elem(ie)%state%v(:,:,:,k,n0),deriv,elem(ie))
  
        ! convert to covariant form
        do j=1,np
          do i=1,np
             gv(i,j,1)=(elem(ie)%D(1,1,i,j)*elem(ie)%state%v(i,j,1,k,n0) + elem(ie)%D(2,1,i,j)*elem(ie)%state%v(i,j,2,k,n0))
             gv(i,j,2)=(elem(ie)%D(1,2,i,j)*elem(ie)%state%v(i,j,1,k,n0) + elem(ie)%D(2,2,i,j)*elem(ie)%state%v(i,j,2,k,n0))
          enddo
       enddo   

       do j=1,np
          do l=1,np
            dsdy00=0.0d0
            dsdx00=0.0d0
            do i=1,np
               dsdx00 = dsdx00 + deriv%Dvv(i,l  )*gv(i,j  ,2)
               dsdy00 = dsdy00 + deriv%Dvv(i,l  )*gv(j  ,i,1)
            enddo
            vort  (l  ,j  ,k  ) = dsdx00
            v_tmp1(j  ,l  ) = dsdy00
          enddo
        enddo

        do j=1,np
          do i=1,np
             vort(i,j,k)=(vort(i,j,k)-v_tmp1(i,j))*(elem(ie)%rmetdet(i,j)*rrearth)
          enddo
        enddo

     enddo
     
     ! compute T_v for timelevel n0
     !if ( moisture /= "dry") then
     if (qn0 == -1 ) then
        do k=1,nlev
           do j=1,np
              do i=1,np
                 T_v(i,j,k) = elem(ie)%state%T(i,j,k,n0)
                 kappa_star(i,j,k) = kappa
              end do
           end do
        end do
     else
        do k=1,nlev
           do j=1,np
              do i=1,np
                 ! Qt = elem(ie)%state%Q(i,j,k,1) 
                 Qt = elem(ie)%state%Qdp(i,j,k,1,qn0)/dp(i,j,k)
                 !T_v(i,j,k) = Virtual_Temperature(elem(ie)%state%T(i,j,k,n0),Qt)
                  T_v(i,j,k) = elem(ie)%state%T(i,j,k,n0)*(1_real_kind + (Rwater_vapor/Rgas - 1.0_real_kind)*Qt)
                 if (use_cpstar==1) then
                    !kappa_star(i,j,k) =  Rgas/Virtual_Specific_Heat(Qt)
                    kappa_star(i,j,k) =  Rgas/Cp*(1.0_real_kind + (Cpwater_vapor/Cp - 1.0_real_kind)*Qt)
                 else
                    kappa_star(i,j,k) = kappa
                 endif
              end do
           end do
        end do
     end if
     
     
     ! ====================================================
     ! Compute Hydrostatic equation, modeld after CCM-3
     ! ====================================================
     !call geopotential_t(p,dp,T_v,Rgas,phi)
    !inline here call preq_hydrostatic(phi,elem(ie)%state%phis,T_v,p,dp)

     do j=1,np   !   Loop inversion (AAM)

          do i=1,np
             hkk = dp(i,j,nlev)*0.5d0/p(i,j,nlev)
             hkl = 2*hkk
             phii(i,j,nlev)  = Rgas*T_v(i,j,nlev)*hkl
             phi_oacc(i,j,nlev,ie) = elem(ie)%state%phis(i,j) + Rgas*T_v(i,j,nlev)*hkk
          end do

          do k=nlev-1,2,-1
             do i=1,np
                ! hkk = dp*ckk
                hkk = dp(i,j,k)*0.5d0/p(i,j,k)
                hkl = 2*hkk
                phii(i,j,k) = phii(i,j,k+1) + Rgas*T_v(i,j,k)*hkl
                phi_oacc(i,j,k,ie) = elem(ie)%state%phis(i,j) + phii(i,j,k+1) + Rgas*T_v(i,j,k)*hkk
             end do
          end do

          do i=1,np
             ! hkk = dp*ckk
             hkk = 0.5d0*dp(i,j,1)/p(i,j,1)
             phi_oacc(i,j,1,ie) = elem(ie)%state%phis(i,j) + phii(i,j,2) + Rgas*T_v(i,j,1)*hkk
          end do

       end do
     


     ! ====================================================
     ! Compute omega_p according to CCM-3 
     ! ====================================================
    ! call preq_omega_ps(omega_p,hvcoord,p,vgrad_p,divdp)

    do j=1,np   !   Loop inversion (AAM)

          do i=1,np
             hkk = 0.5d0/p(i,j,1)
             term = divdp(i,j,1)
!             omega_p(i,j,1) = hvcoord%hybm(1)*vgrad_ps(i,j,1)/p(i,j,1)
             omega_p(i,j,1) = vgrad_p(i,j,1)/p(i,j,1)
             omega_p(i,j,1) = omega_p(i,j,1) - hkk*term
             v_tmp1(i,j) = term
          end do

          do k=2,nlev-1
             do i=1,np
                hkk = 0.5d0/p(i,j,k)
                hkl = 2*hkk
                term = divdp(i,j,k)
!                omega_p(i,j,k) = hvcoord%hybm(k)*vgrad_ps(i,j,k)/p(i,j,k)
                omega_p(i,j,k) = vgrad_p(i,j,k)/p(i,j,k)
                omega_p(i,j,k) = omega_p(i,j,k) - hkl*v_tmp1(i,j) - hkk*term
                v_tmp1(i,j) = v_tmp1(i,j) + term

             end do
          end do

          do i=1,np
             hkk = 0.5d0/p(i,j,nlev)
             hkl = 2*hkk
             term = divdp(i,j,nlev)
!             omega_p(i,j,nlev) = hvcoord%hybm(nlev)*vgrad_ps(i,j,nlev)/p(i,j,nlev)
             omega_p(i,j,nlev) = vgrad_p(i,j,nlev)/p(i,j,nlev)
             omega_p(i,j,nlev) = omega_p(i,j,nlev) - hkl*v_tmp1(i,j) - hkk*term
          end do

       end do


     
     ! ==================================================
     ! zero partial sum for accumulating sum
     !    (div(v_k) + v_k.grad(lnps))*dsigma_k = div( v dp )
     ! used by eta_dot_dpdn and lnps tendency
     ! ==================================================
     sdot_sum=0
     

     ! ==================================================
     ! Compute eta_dot_dpdn 
     ! save sdot_sum as this is the -RHS of ps_v equation
     ! ==================================================
     if (rsplit>0) then
        ! VERTICALLY LAGRANGIAN:   no vertical motion
        eta_dot_dpdn=0    
        T_vadv=0          
        v_vadv=0
     else
        do k=1,nlev
           ! ==================================================
           ! add this term to PS equation so we exactly conserve dry mass
           ! ==================================================
           do j=1,np
              do i=1,np
                 sdot_sum(i,j) = sdot_sum(i,j) + divdp(i,j,k) 
                 eta_dot_dpdn(i,j,k+1) = sdot_sum(i,j) 
              enddo
           enddo     
        enddo
     
     
        ! ===========================================================
        ! at this point, eta_dot_dpdn contains integral_etatop^eta[ divdp ]
        ! compute at interfaces:
        !    eta_dot_dpdn = -dp/dt - integral_etatop^eta[ divdp ]
        ! for reference: at mid layers we have:
        !    omega = v grad p  - integral_etatop^eta[ divdp ]
        ! ===========================================================
        do k=1,nlev-1
           do j=1,np
              do i=1,np
                 eta_dot_dpdn(i,j,k+1) = hvcoord%hybi(k+1)*sdot_sum(i,j) - eta_dot_dpdn(i,j,k+1)
              enddo
            enddo
        enddo
        
        do j=1,np
          do i=1,np
            eta_dot_dpdn(i,j,1     ) = 0.0D0
            eta_dot_dpdn(i,j,nlev+1) = 0.0D0
          enddo
        enddo
        
        ! ===========================================================
        ! Compute vertical advection of T and v from eq. CCM2 (3.b.1)
        ! ==============================================

        !inline   call preq_vertadv(elem(ie)%state%T(:,:,:,n0),elem(ie)%state%v(:,:,:,:,n0), eta_dot_dpdn,rdp,T_vadv,v_vadv)
                                    !                                                  eta_dot_dp_deta, rpdel, T_vadv, v_vadv)
         do j=1,np   !   Loop inversion (AAM)

          ! ===========================================================
          ! Compute vertical advection of T and v from eq. (3.b.1)
          !
          ! k = 1 case:
          ! ===========================================================

           k=1
            do i=1,np
               hkk             = (0.5_real_kind*rdp(i,j,k))*eta_dot_dpdn(i,j,k+1)
               T_vadv(i,j,k)   = hkk*(elem(ie)%state%T(i,j,k+1,n0)   - elem(ie)%state%T(i,j,k,n0))
               v_vadv(i,j,1,k) = hkk*(elem(ie)%state%v(i,j,1,k+1,n0) - elem(ie)%state%v(i,j,1,k,n0))
               v_vadv(i,j,2,k) = hkk*(elem(ie)%state%v(i,j,2,k+1,n0) - elem(ie)%state%v(i,j,2,k,n0))
            end do

          ! ===========================================================
          ! vertical advection
          !
          ! 1 < k < nlev case:
          ! ===========================================================

          do k=2,nlev-1
             do i=1,np
                hkk            = (0.5_real_kind*rdp(i,j,k))*eta_dot_dpdn(i,j,k+1)
                hkl            = (0.5_real_kind*rdp(i,j,k))*eta_dot_dpdn(i,j,k)
                T_vadv(i,j,k)   = hkk*(elem(ie)%state%T(i,j,k+1,n0)   - elem(ie)%state%T(i,j,k,n0)) + &
                                  hkl*(elem(ie)%state%T(i,j,k,n0)     - elem(ie)%state%T(i,j,k-1,n0))
                v_vadv(i,j,1,k) = hkk*(elem(ie)%state%v(i,j,1,k+1,n0) - elem(ie)%state%v(i,j,1,k,n0)) + &
                                  hkl*(elem(ie)%state%v(i,j,1,k,n0)   - elem(ie)%state%v(i,j,1,k-1,n0))
                v_vadv(i,j,2,k) = hkk*(elem(ie)%state%v(i,j,2,k+1,n0) - elem(ie)%state%v(i,j,2,k,n0)) + &
                                  hkl*(elem(ie)%state%v(i,j,2,k,n0)   - elem(ie)%state%v(i,j,2,k-1,n0))
             end do 
           end do

          ! ===========================================================
          ! vertical advection
          !
          ! k = nlev case:
          ! ===========================================================
 
          k=nlev
          do i=1,np
              hkk             = (0.5_real_kind*rdp(i,j,k))*eta_dot_dpdn(i,j,k)
              T_vadv(i,j,k)   = hkl*(elem(ie)%state%T(i,j,k,n0)- elem(ie)%state%T(i,j,k-1,n0))
              v_vadv(i,j,1,k) = hkl*(elem(ie)%state%v(i,j,1,k,n0)- elem(ie)%state%v(i,j,1,k-1,n0))
              v_vadv(i,j,2,k) = hkl*(elem(ie)%state%v(i,j,2,k,n0)- elem(ie)%state%v(i,j,2,k-1,n0))
          end do

        enddo

 

     endif


     ! ================================
     ! accumulate mean vertical flux:
     ! ================================
     do k=1,nlev  !  Loop index added (AAM)
        do j=1,np
           do i=1,np
               elem(ie)%derived%eta_dot_dpdn(i,j,k) = &
                       elem(ie)%derived%eta_dot_dpdn(i,j,k) + eta_ave_w*eta_dot_dpdn(i,j,k)
               elem(ie)%derived%omega_p(i,j,k) = &
                       elem(ie)%derived%omega_p(i,j,k) + eta_ave_w*omega_p(i,j,k)
           enddo
        enddo
     enddo

     do j=1,np
        do i=1,np
           elem(ie)%derived%eta_dot_dpdn(i,j,nlev+1) = &
                 elem(ie)%derived%eta_dot_dpdn(i,j,nlev+1) + eta_ave_w*eta_dot_dpdn(i,j,nlev+1)
        enddo
     enddo

     
     
     ! ==============================================
     ! Compute phi + kinetic energy term: 10*nv*nv Flops
     ! ==============================================
     vertloop: do k=1,nlev   
        do j=1,np
           do i=1,np
              v1     = elem(ie)%state%v(i,j,1,k,n0)
              v2     = elem(ie)%state%v(i,j,2,k,n0)
              E = 0.5D0*( v1*v1 + v2*v2 )
              Ephi(i,j)=E+phi_oacc(i,j,k,ie)+elem(ie)%derived%pecnd(i,j,k)
           end do
        end do
        ! ================================================
        ! compute gradp term (ps/p)*(dp/dps)*T
        ! ================================================
        !! inline gradient_sphere here vtemp(:,:,:)   = gradient_sphere(elem(ie)%state%T(:,:,k,n0),deriv,elem(ie)%Dinv)
          do j=1,np
            do l=1,np
             dsdx00=0.0d0
             dsdy00=0.0d0
              do i=1,np
               dsdx00 = dsdx00 + deriv%Dvv(i,l  )*elem(ie)%state%T(i,j,k,n0)
               dsdy00 = dsdy00 + deriv%Dvv(i,l  )*elem(ie)%state%T(j,i,k,n0)
              enddo
              v_tmp1(l  ,j  ) = dsdx00*rrearth
              v_tmp2(j  ,l  ) = dsdy00*rrearth
             enddo
           enddo
           do j=1,np
            do i=1,np
             vtemp(i,j,1)=elem(ie)%Dinv(1,1,i,j)*v_tmp1(i,j) + elem(ie)%Dinv(2,1,i,j)*v_tmp2(i,j)
             vtemp(i,j,2)=elem(ie)%Dinv(1,2,i,j)*v_tmp1(i,j) + elem(ie)%Dinv(2,2,i,j)*v_tmp2(i,j)
            enddo
           enddo 



        do j=1,np
           do i=1,np
              v1     = elem(ie)%state%v(i,j,1,k,n0)
              v2     = elem(ie)%state%v(i,j,2,k,n0)
              vgrad_T(i,j) =  v1*vtemp(i,j,1) + v2*vtemp(i,j,2) 
           end do
        end do
        
        
        ! vtemp = grad ( E + PHI )
        !inline vtemp = gradient_sphere(Ephi(:,:),deriv,elem(ie)%Dinv)
        do j=1,np
            do l=1,np
             dsdx00=0.0d0
             dsdy00=0.0d0
              do i=1,np
               dsdx00 = dsdx00 + deriv%Dvv(i,l  )*Ephi(i,j)
               dsdy00 = dsdy00 + deriv%Dvv(i,l  )*Ephi(j,i)
              enddo
              v_tmp1(l  ,j  ) = dsdx00*rrearth
              v_tmp2(j  ,l  ) = dsdy00*rrearth
             enddo
           enddo
           do j=1,np
            do i=1,np
             vtemp(i,j,1)=elem(ie)%Dinv(1,1,i,j)*v_tmp1(i,j) + elem(ie)%Dinv(2,1,i,j)*v_tmp2(i,j)
             vtemp(i,j,2)=elem(ie)%Dinv(1,2,i,j)*v_tmp1(i,j) + elem(ie)%Dinv(2,2,i,j)*v_tmp2(i,j)
            enddo
           enddo        

        do j=1,np
           do i=1,np
!              gpterm = hvcoord%hybm(k)*T_v(i,j,k)/p(i,j,k)
!              glnps1 = Rgas*gpterm*grad_ps(i,j,1)
!              glnps2 = Rgas*gpterm*grad_ps(i,j,2)
              gpterm = T_v(i,j,k)/p(i,j,k)
              glnps1 = Rgas*gpterm*grad_p(i,j,1,k)
              glnps2 = Rgas*gpterm*grad_p(i,j,2,k)
              
              v1     = elem(ie)%state%v(i,j,1,k,n0)
              v2     = elem(ie)%state%v(i,j,2,k,n0)
              
              vtens1(i,j,k) =   - v_vadv(i,j,1,k)                           &
                   + v2*(elem(ie)%fcor(i,j) + vort(i,j,k))        &
                   - vtemp(i,j,1) - glnps1   
              
              vtens2(i,j,k) =   - v_vadv(i,j,2,k)                            &
                   - v1*(elem(ie)%fcor(i,j) + vort(i,j,k))        &
                   - vtemp(i,j,2) - glnps2   
              
              ttens(i,j,k)  = - T_vadv(i,j,k) - vgrad_T(i,j) + kappa_star(i,j,k)*T_v(i,j,k)*omega_p(i,j,k)
              !
              ! phl: add forcing term to T
              !

#if ( defined CAM )

              if (se_prescribed_wind_2d) then
                 vtens1(i,j,k) = 0.D0
                 vtens2(i,j,k) = 0.D0
                 ttens(i,j,k) = 0.D0
              else
                 if(se_met_nudge_u.gt.0.D0)then
                    u_m_umet = v1 - &
                         elem(ie)%derived%u_met(i,j,k) - &
                         se_met_tevolve*tevolve*elem(ie)%derived%dudt_met(i,j,k)
                    v_m_vmet = v2 - &
                         elem(ie)%derived%v_met(i,j,k) - &
                         se_met_tevolve*tevolve*elem(ie)%derived%dvdt_met(i,j,k)

                    vtens1(i,j,k) =   vtens1(i,j,k) - se_met_nudge_u*u_m_umet * elem(ie)%derived%nudge_factor(i,j,k)

                    elem(ie)%derived%Utnd(i+(j-1)*np,k) = elem(ie)%derived%Utnd(i+(j-1)*np,k) &
                         + se_met_nudge_u*u_m_umet * elem(ie)%derived%nudge_factor(i,j,k)

                    vtens2(i,j,k) =   vtens2(i,j,k) - se_met_nudge_u*v_m_vmet * elem(ie)%derived%nudge_factor(i,j,k)

                    elem(ie)%derived%Vtnd(i+(j-1)*np,k) = elem(ie)%derived%Vtnd(i+(j-1)*np,k) &
                         + se_met_nudge_u*v_m_vmet * elem(ie)%derived%nudge_factor(i,j,k)

                 endif

                 if(se_met_nudge_p.gt.0.D0)then
                    vtens1(i,j,k) =   vtens1(i,j,k) - se_met_nudge_p*grad_p_m_pmet(i,j,1,k)  * elem(ie)%derived%nudge_factor(i,j,k)
                    vtens2(i,j,k) =   vtens2(i,j,k) - se_met_nudge_p*grad_p_m_pmet(i,j,2,k)  * elem(ie)%derived%nudge_factor(i,j,k)
                 endif

                 if(se_met_nudge_t.gt.0.D0)then
                    t_m_tmet = elem(ie)%state%T(i,j,k,n0) - &
                         elem(ie)%derived%T_met(i,j,k) - &
                         se_met_tevolve*tevolve*elem(ie)%derived%dTdt_met(i,j,k)
                    ttens(i,j,k)  = ttens(i,j,k) - se_met_nudge_t*t_m_tmet * elem(ie)%derived%nudge_factor(i,j,k)
                    elem(ie)%derived%Ttnd(i+(j-1)*np,k) = elem(ie)%derived%Ttnd(i+(j-1)*np,k) &
                         + se_met_nudge_t*t_m_tmet * elem(ie)%derived%nudge_factor(i,j,k)
                 endif
              endif
#endif


           end do
        end do

     end do vertloop 

#ifdef ENERGY_DIAGNOSTICS
     ! =========================================================
     !
     ! diagnostics
     ! recomputes some gradients that were not saved above
     ! uses:  sdot_sum(), eta_dot_dpdn(), grad_ps()
     ! grad_phi(), dp(), p(), T_vadv(), v_vadv(), divdp()
     ! =========================================================

     ! =========================================================
     ! (AAM) - This section has accumulations over vertical levels.
     !   Be careful if implementing OpenMP
     ! =========================================================

     if (compute_diagnostics) then
        elem(ie)%accum%KEhorz1=0
        elem(ie)%accum%KEhorz2=0
        elem(ie)%accum%IEhorz1=0
        elem(ie)%accum%IEhorz2=0
        elem(ie)%accum%IEhorz1_wet=0
        elem(ie)%accum%IEhorz2_wet=0
        elem(ie)%accum%KEvert1=0
        elem(ie)%accum%KEvert2=0
        elem(ie)%accum%IEvert1=0
        elem(ie)%accum%IEvert2=0
        elem(ie)%accum%IEvert1_wet=0
        elem(ie)%accum%IEvert2_wet=0
        elem(ie)%accum%T1=0
        elem(ie)%accum%T2=0
        elem(ie)%accum%T2_s=0
        elem(ie)%accum%S1=0
        elem(ie)%accum%S1_wet=0
        elem(ie)%accum%S2=0
        
        do j=1,np
           do i=1,np
              elem(ie)%accum%S2(i,j) = elem(ie)%accum%S2(i,j) - &
                   sdot_sum(i,j)*elem(ie)%state%phis(i,j)
           enddo
        enddo
        
        do k=1,nlev
           ! vtemp = grad_E(:,:,k)
           do j=1,np
              do i=1,np
                 v1     = elem(ie)%state%v(i,j,1,k,n0)
                 v2     = elem(ie)%state%v(i,j,2,k,n0)
                 Ephi(i,j)=0.5D0*( v1*v1 + v2*v2 )
              enddo
           enddo
           !inline vtemp = gradient_sphere(Ephi,deriv,elem(ie)%Dinv)
           do j=1,np
            do l=1,np
             dsdx00=0.0d0
             dsdy00=0.0d0
              do i=1,np
               dsdx00 = dsdx00 + deriv%Dvv(i,l  )*Ephi(i,j)
               dsdy00 = dsdy00 + deriv%Dvv(i,l  )*Ephi(j,i)
              enddo
              v_tmp1(l  ,j  ) = dsdx00*rrearth
              v_tmp2(j  ,l  ) = dsdy00*rrearth
             enddo
           enddo
           do j=1,np
            do i=1,np
             vtemp(i,j,1)=elem(ie)%Dinv(1,1,i,j)*v_tmp1(i,j) + elem(ie)%Dinv(2,1,i,j)*v_tmp2(i,j)
             vtemp(i,j,2)=elem(ie)%Dinv(1,2,i,j)*v_tmp1(i,j) + elem(ie)%Dinv(2,2,i,j)*v_tmp2(i,j)
            enddo
           enddo

           do j=1,np
              do i=1,np
                 ! dp/dn u dot grad(E)
                 v1     = elem(ie)%state%v(i,j,1,k,n0)
                 v2     = elem(ie)%state%v(i,j,2,k,n0)
                 elem(ie)%accum%KEhorz2(i,j) = elem(ie)%accum%KEhorz2(i,j) + &
                      (v1*vtemp(i,j,1)  + v2*vtemp(i,j,2))*dp(i,j,k)
                 ! E div( u dp/dn )
                 elem(ie)%accum%KEhorz1(i,j) = elem(ie)%accum%KEhorz1(i,j) + Ephi(i,j)*divdp(i,j,k)
                 
                 ! Cp T div( u dp/dn)   ! dry horizontal advection component
                 elem(ie)%accum%IEhorz1(i,j) = elem(ie)%accum%IEhorz1(i,j) + Cp*elem(ie)%state%T(i,j,k,n0)*divdp(i,j,k)
                 
                 
              enddo
           enddo
           
           
           ! vtemp = grad_phi(:,:,k)
           !inline vtemp = gradient_sphere(phi(:,:,k),deriv,elem(ie)%Dinv)
           do j=1,np
            do l=1,np
             dsdx00=0.0d0
             dsdy00=0.0d0
              do i=1,np
               dsdx00 = dsdx00 + deriv%Dvv(i,l  )*phi_oacc(i,j,k,ie)
               dsdy00 = dsdy00 + deriv%Dvv(i,l  )*phi_oacc(j,i,k,ie)
              enddo
              v_tmp1(l  ,j  ) = dsdx00*rrearth
              v_tmp2(j  ,l  ) = dsdy00*rrearth
             enddo
           enddo
           do j=1,np
            do i=1,np
             vtemp(i,j,1)=elem(ie)%Dinv(1,1,i,j)*v_tmp1(i,j) + elem(ie)%Dinv(2,1,i,j)*v_tmp2(i,j)
             vtemp(i,j,2)=elem(ie)%Dinv(1,2,i,j)*v_tmp1(i,j) + elem(ie)%Dinv(2,2,i,j)*v_tmp2(i,j)
            enddo
           enddo


           do j=1,np
              do i=1,np
                 v1     = elem(ie)%state%v(i,j,1,k,n0)
                 v2     = elem(ie)%state%v(i,j,2,k,n0)
                 E = 0.5D0*( v1*v1 + v2*v2 )
                 ! NOTE:  Cp_star = Cp + (Cpv-Cp)*q   
                 ! advection terms can thus be broken into two components: dry and wet
                 ! dry components cancel exactly
                 ! wet components should cancel exactly
                 ! 
                 ! some diagnostics
                 ! e = eta_dot_dpdn()    
                 de =  eta_dot_dpdn(i,j,k+1)-eta_dot_dpdn(i,j,k) 
                 ! Cp T de/dn, integral dn:
                 elem(ie)%accum%IEvert1(i,j)=elem(ie)%accum%IEvert1(i,j) + Cp*elem(ie)%state%T(i,j,k,n0)*de
                 ! E de/dn
                 elem(ie)%accum%KEvert1(i,j)=elem(ie)%accum%KEvert1(i,j) + E*de
                 ! Cp T_vadv dp/dn
                 elem(ie)%accum%IEvert2(i,j)=elem(ie)%accum%IEvert2(i,j) + Cp*T_vadv(i,j,k)*dp(i,j,k)
                 ! dp/dn V dot V_vadv
                 elem(ie)%accum%KEvert2(i,j)=elem(ie)%accum%KEvert2(i,j) + (v1*v_vadv(i,j,1,k) + v2*v_vadv(i,j,2,k)) *dp(i,j,k)
                 
                 ! IEvert1_wet():  (Cpv-Cp) T Qdp_vadv  (Q equation)
                 ! IEvert2_wet():  (Cpv-Cp) Qdp T_vadv   T equation
                 if (use_cpstar==1) then
                 elem(ie)%accum%IEvert2_wet(i,j)=elem(ie)%accum%IEvert2_wet(i,j) +&
                      (Cpwater_vapor-Cp)*elem(ie)%state%Q(i,j,k,1)*T_vadv(i,j,k)*dp(i,j,k)
                 endif

                 gpterm = T_v(i,j,k)/p(i,j,k)
                 elem(ie)%accum%T1(i,j) = elem(ie)%accum%T1(i,j) - &
                      Rgas*gpterm*(grad_p(i,j,1,k)*v1 + grad_p(i,j,2,k)*v2)*dp(i,j,k)
                 
                 elem(ie)%accum%T2(i,j) = elem(ie)%accum%T2(i,j) - &
                      (vtemp(i,j,1)*v1 + vtemp(i,j,2)*v2)*dp(i,j,k)
                 
                 ! S1 = < Cp_star dp/dn , RT omega_p/cp_star >  
                 elem(ie)%accum%S1(i,j) = elem(ie)%accum%S1(i,j) + &
                      Rgas*T_v(i,j,k)*omega_p(i,j,k)*dp(i,j,k)
                 
                 ! cp_star = cp + cp2
                 if (use_cpstar==1) then
                 cp2 = (Cpwater_vapor-Cp)*elem(ie)%state%Q(i,j,k,1)
                 cp_ratio = cp2/(cp+cp2)
                 elem(ie)%accum%S1_wet(i,j) = elem(ie)%accum%S1_wet(i,j) + &
                      cp_ratio*(Rgas*T_v(i,j,k)*omega_p(i,j,k)*dp(i,j,k))
                 endif
                 
                 elem(ie)%accum%CONV(i,j,:,k)=-Rgas*gpterm*grad_p(i,j,:,k)-vtemp(i,j,:)
              enddo
           enddo
           
           ! inline vtemp(:,:,:) = gradient_sphere(elem(ie)%state%phis(:,:),deriv,elem(ie)%Dinv)
           do j=1,np
            do l=1,np
             dsdx00=0.0d0
             dsdy00=0.0d0
              do i=1,np
               dsdx00 = dsdx00 + deriv%Dvv(i,l  )*elem(ie)%state%phis(i,j)
               dsdy00 = dsdy00 + deriv%Dvv(i,l  )*elem(ie)%state%phis(j,i)
              enddo
              v_tmp1(l  ,j  ) = dsdx00*rrearth
              v_tmp2(j  ,l  ) = dsdy00*rrearth
             enddo
           enddo
           do j=1,np
            do i=1,np
             vtemp(i,j,1)=elem(ie)%Dinv(1,1,i,j)*v_tmp1(i,j) + elem(ie)%Dinv(2,1,i,j)*v_tmp2(i,j)
             vtemp(i,j,2)=elem(ie)%Dinv(1,2,i,j)*v_tmp1(i,j) + elem(ie)%Dinv(2,2,i,j)*v_tmp2(i,j)
            enddo
           enddo

           do j=1,np
              do i=1,np
                 v1     = elem(ie)%state%v(i,j,1,k,n0)
                 v2     = elem(ie)%state%v(i,j,2,k,n0)
                 elem(ie)%accum%T2_s(i,j) = elem(ie)%accum%T2_s(i,j) - &
                      (vtemp(i,j,1)*v1 + vtemp(i,j,2)*v2)*dp(i,j,k)
              enddo
           enddo
           
           !inline vtemp(:,:,:)   = gradient_sphere(elem(ie)%state%T(:,:,k,n0),deriv,elem(ie)%Dinv)
           do j=1,np
            do l=1,np
             dsdx00=0.0d0
             dsdy00=0.0d0
              do i=1,np
               dsdx00 = dsdx00 + deriv%Dvv(i,l  )*elem(ie)%state%T(i,j,k,n0)
               dsdy00 = dsdy00 + deriv%Dvv(i,l  )*elem(ie)%state%T(j,i,k,n0)
              enddo
              v_tmp1(l  ,j  ) = dsdx00*rrearth
              v_tmp2(j  ,l  ) = dsdy00*rrearth
             enddo
           enddo
           do j=1,np
            do i=1,np
             vtemp(i,j,1)=elem(ie)%Dinv(1,1,i,j)*v_tmp1(i,j) + elem(ie)%Dinv(2,1,i,j)*v_tmp2(i,j)
             vtemp(i,j,2)=elem(ie)%Dinv(1,2,i,j)*v_tmp1(i,j) + elem(ie)%Dinv(2,2,i,j)*v_tmp2(i,j)
            enddo
           enddo

           do j=1,np
              do i=1,np
                 v1     = elem(ie)%state%v(i,j,1,k,n0)
                 v2     = elem(ie)%state%v(i,j,2,k,n0)
                 
                 ! Cp dp/dn u dot gradT
                 elem(ie)%accum%IEhorz2(i,j) = elem(ie)%accum%IEhorz2(i,j) + &
                      Cp*(v1*vtemp(i,j,1) + v2*vtemp(i,j,2))*dp(i,j,k)
                 
                 if (use_cpstar==1) then
                 elem(ie)%accum%IEhorz2_wet(i,j) = elem(ie)%accum%IEhorz2_wet(i,j) + &
                      (Cpwater_vapor-Cp)*elem(ie)%state%Q(i,j,k,1)*&
                      (v1*vtemp(i,j,1) + v2*vtemp(i,j,2))*dp(i,j,k)
                 endif
                 
              enddo
           enddo
           
        enddo
     endif
#endif
     ! =========================================================
     ! local element timestep, store in np1.  
     ! note that we allow np1=n0 or nm1
     ! apply mass matrix
     ! =========================================================
     if (dt2<0) then
        ! calling program just wanted DSS'd RHS, skip time advance
        do k=1,nlev
          do j=1,np
             do i=1,np
              elem(ie)%state%v(i,j,1,k,np1) = elem(ie)%spheremp(i,j)*vtens1(i,j,k) 
              elem(ie)%state%v(i,j,2,k,np1) = elem(ie)%spheremp(i,j)*vtens2(i,j,k) 
              elem(ie)%state%T(i,j,k,np1)   = elem(ie)%spheremp(i,j)*ttens(i,j,k)
              if (rsplit>0) &
                  elem(ie)%state%dp3d(i,j,k,np1) = -elem(ie)%spheremp(i,j)*&
                     (divdp(i,j,k) + eta_dot_dpdn(i,j,k+1)-eta_dot_dpdn(i,j,k)) 
              if (0<rsplit.and.0<ntrac.and.eta_ave_w.ne.0.) then
              v(i,j,1) =  elem(ie)%Dinv(i,j,1,1)*vdp(i,j,1,k) + elem(ie)%Dinv(i,j,1,2)*vdp(i,j,2,k)
              v(i,j,2) =  elem(ie)%Dinv(i,j,2,1)*vdp(i,j,1,k) + elem(ie)%Dinv(i,j,2,2)*vdp(i,j,2,k)
!uncomment after merge               tempflux =  eta_ave_w*subcell_div_fluxes(v, np, nc, elem(ie)%metdet)
!Irina TOFIX
!uncomment after merge              elem(ie)%sub_elem_mass_flux(i,j,:,k) = elem(ie)%sub_elem_mass_flux(i,j,:,k) - tempflux
           end if

             enddo
          enddo
        enddo
        do j=1,np
             do i=1,np
                 elem(ie)%state%ps_v(i,j,np1) = -elem(ie)%spheremp(i,j)*sdot_sum(i,j)
             enddo
        enddo
     else
        do k=1,nlev
          do j=1,np
             do i=1,np
               elem(ie)%state%v(i,j,1,k,np1) = elem(ie)%spheremp(i,j)*( elem(ie)%state%v(i,j,1,k,nm1) + dt2*vtens1(i,j,k) )
               elem(ie)%state%v(i,j,2,k,np1) = elem(ie)%spheremp(i,j)*( elem(ie)%state%v(i,j,2,k,nm1) + dt2*vtens2(i,j,k) )
               elem(ie)%state%T(i,j,k,np1)   = elem(ie)%spheremp(i,j)*( elem(ie)%state%T(i,j,k,nm1)   + dt2*ttens(i,j,k))
                if (rsplit>0) &
                    elem(ie)%state%dp3d(i,j,k,np1) = elem(ie)%spheremp(i,j)*&
                      (elem(ie)%state%dp3d(i,j,k,nm1)-dt2*&
                      (divdp(i,j,k) + eta_dot_dpdn(i,j,k+1)-eta_dot_dpdn(i,j,k)))
                if (0<rsplit.and.0<ntrac.and.eta_ave_w.ne.0.) then
                  v(i,j,1) =  elem(ie)%Dinv(i,j,1,1)*vdp(i,j,1,k) + elem(ie)%Dinv(i,j,1,2)*vdp(i,j,2,k)
                  v(i,j,2) =  elem(ie)%Dinv(i,j,2,1)*vdp(i,j,1,k) + elem(ie)%Dinv(i,j,2,2)*vdp(i,j,2,k)
!uncomment after merge                   tempflux =  eta_ave_w*subcell_div_fluxes(v, np, nc, elem(ie)%metdet)
!Irina TOFiX
!uncomment after merge                  elem(ie)%sub_elem_mass_flux(i,j,:,k) = elem(ie)%sub_elem_mass_flux(i,j,:,k) - tempflux
           end if
             enddo
          enddo
        enddo
        do j=1,np
           do i=1,np
                elem(ie)%state%ps_v(i,j,np1) = elem(ie)%spheremp(i,j)*( elem(ie)%state%ps_v(i,j,nm1) - dt2*sdot_sum(i,j) )
           enddo
        enddo

     endif

  enddo


  do ie=nets,nete     
     ! =========================================================
     !
     ! Pack ps(np1), T, and v tendencies into comm buffer
     !
     ! =========================================================
     kptr=0
!uncomment after merge: d foe all edgeVpack:
!call edgeVpack(edge3p1, elem(ie)%state%ps_v(:,:,np1),1,kptr,ie)
     call edgeVpack(edge3p1, elem(ie)%state%ps_v(:,:,np1),1,kptr,elem(ie)%desc)
     
     kptr=1
     call edgeVpack(edge3p1, elem(ie)%state%T(:,:,:,np1),nlev,kptr,elem(ie)%desc)
     
     kptr=nlev+1
     call edgeVpack(edge3p1, elem(ie)%state%v(:,:,:,:,np1),2*nlev,kptr,elem(ie)%desc)

     if (rsplit>0) then
        kptr=kptr+2*nlev
        call edgeVpack(edge3p1, elem(ie)%state%dp3d(:,:,:,np1),nlev,kptr,elem(ie)%desc)
     endif
  end do
  
  ! =============================================================
    ! Insert communications here: for shared memory, just a single
  ! sync is required
  ! =============================================================

  call t_startf('caar_bexchV_oacc')
  call bndry_exchangeV(hybrid,edge3p1)
  call t_stopf('caar_bexchV_oacc')

  do ie=nets,nete
     ! ===========================================================
     ! Unpack the edges for vgrad_T and v tendencies...
     ! ===========================================================
     kptr=0
!Irina TODO
     call edgeVunpack(edge3p1, elem(ie)%state%ps_v(:,:,np1), 1, kptr, elem(ie)%desc)
     
     kptr=1
     call edgeVunpack(edge3p1, elem(ie)%state%T(:,:,:,np1), nlev, kptr, elem(ie)%desc)
     
     kptr=nlev+1
     call edgeVunpack(edge3p1, elem(ie)%state%v(:,:,:,:,np1), 2*nlev, kptr, elem(ie)%desc)
     
     if (rsplit>0) then
         if (0<ntrac.and.eta_ave_w.ne.0.) then
          do k=1,nlev
             stashdp3d(:,:,k) = elem(ie)%state%dp3d(:,:,k,np1)/elem(ie)%spheremp(:,:)
          end do
        endif

        corners = 0.0d0
        corners(1:np,1:np,:) = elem(ie)%state%dp3d(:,:,:,np1)

        kptr=kptr+2*nlev
        call edgeVunpack(edge3p1, elem(ie)%state%dp3d(:,:,:,np1),nlev,kptr,elem(ie)%desc)
       
!uncomment after merge
!        if  (0<ntrac.and.eta_ave_w.ne.0.) then
!          desc = elem(ie)%desc
!          call edgeDGVunpack(edge3p1, corners, nlev, kptr, ie)
!          corners = corners/dt2

!          do k=1,nlev
!            tempdp3d = elem(ie)%rspheremp(:,:)*elem(ie)%state%dp3d(:,:,k,np1)
!            tempdp3d = tempdp3d - stashdp3d(:,:,k)
!            tempdp3d = tempdp3d/dt2

!            call distribute_flux_at_corners(cflux, corners(:,:,k), desc%getmapP)

!            cflux(1,1,:)   = elem(ie)%rspheremp(1,  1) * cflux(1,1,:)
!            cflux(2,1,:)   = elem(ie)%rspheremp(np, 1) * cflux(2,1,:)
!            cflux(1,2,:)   = elem(ie)%rspheremp(1, np) * cflux(1,2,:)
!            cflux(2,2,:)   = elem(ie)%rspheremp(np,np) * cflux(2,2,:)

!            tempflux =  eta_ave_w*subcell_dss_fluxes(tempdp3d, np, nc, elem(ie)%metdet, cflux)
!            elem(ie)%sub_elem_mass_flux(:,:,:,k) = elem(ie)%sub_elem_mass_flux(:,:,:,k) + tempflux
!          end do
!        end if
        

     endif
     
     ! ====================================================
     ! Scale tendencies by inverse mass matrix
     ! ====================================================

     do k=1,nlev
       do j=1,np
          do i=1,np
             elem(ie)%state%T(i,j,k,np1)   = elem(ie)%rspheremp(i,j)*elem(ie)%state%T(i,j,k,np1)
             elem(ie)%state%v(i,j,1,k,np1) = elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,1,k,np1)
             elem(ie)%state%v(i,j,2,k,np1) = elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,2,k,np1)
          enddo
        enddo
     enddo

     if (rsplit>0) then
        ! vertically lagrangian: complete dp3d timestep:
        do k=1,nlev
          do j=1,np
             do i=1,np
                 elem(ie)%state%dp3d(i,j,k,np1)= elem(ie)%rspheremp(i,j)*elem(ie)%state%dp3d(i,j,k,np1)
             enddo
          enddo
        enddo
        ! when debugging: also update ps_v
        !elem(ie)%state%ps_v(:,:,np1) = elem(ie)%rspheremp(:,:)*elem(ie)%state%ps_v(:,:,np1)
     else
        ! vertically eulerian: complete ps_v timestep:
        do j=1,np
           do i=1,np
              elem(ie)%state%ps_v(i,j,np1) = elem(ie)%rspheremp(i,j)*elem(ie)%state%ps_v(i,j,np1)
           enddo
        enddo
     endif

  end do

 !$OMP BARRIER
 endif !hybrid%ithr
  call t_stopf('compute_and_apply_oacc_rhs')
  end subroutine compute_and_apply_oacc_rhs


end module prim_advance_oacc_mod


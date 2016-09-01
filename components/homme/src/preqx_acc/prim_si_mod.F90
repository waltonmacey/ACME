
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


module prim_si_mod
  use prim_si_mod_base, only: preq_omegap, preq_omega_lnps, geopotential_t, preq_pressure, preq_hydrostatic, preq_omega_ps, preq_vertadv
#if USE_OPENACC
  implicit none
  private

  public :: preq_omegap, preq_omega_lnps, geopotential_t, preq_pressure, preq_hydrostatic, preq_omega_ps, preq_vertadv
  public :: preq_hydrostatic_openacc
  public :: preq_omega_ps_openacc
  public :: preq_vertadv_openacc


contains



  subroutine preq_vertadv_openacc(elem, eta_dot_dp_deta, rpdel, T_vadv, v_vadv, nets , nete, tl)
    use kinds         , only : real_kind
    use dimensions_mod, only : nlev, np, nlevp, nelemd
    use element_mod   , only : element_t, state_t, state_v
    implicit none
    type(element_t)      , intent(in)  :: elem(:)
    real (kind=real_kind), intent(in)  :: eta_dot_dp_deta(np,np  ,nlevp,nelemd)
    real (kind=real_kind), intent(in)  :: rpdel          (np,np  ,nlev ,nelemd)
    real (kind=real_kind), intent(out) :: T_vadv         (np,np  ,nlev ,nelemd)
    real (kind=real_kind), intent(out) :: v_vadv         (np,np,2,nlev ,nelemd)
    integer              , intent(in)  :: nets, nete, tl
    ! Local Variables
    integer :: i,j,k,ie
    real(kind=real_kind) :: facp(np,np,nlev), facm(np,np,nlev)
    !$acc parallel loop gang private(facm,facp) present(elem,eta_dot_dp_deta,rpdel,t_vadv,v_vadv,state_t,state_v)
    do ie = nets , nete
      !$acc cache(facm,facp)
      !$acc loop vector
      do k = 1 , nlev
        !$acc loop vector
        do j = 1 , np
          !$acc loop vector
          do i = 1 , np
            facm(i,j,k) = 0.5 * rpdel(i,j,k,ie)*eta_dot_dp_deta(i,j,k  ,ie)
            facp(i,j,k) = 0.5 * rpdel(i,j,k,ie)*eta_dot_dp_deta(i,j,k+1,ie)
          enddo
        enddo
      enddo
      !$acc loop vector
      do j = 1 , np
        !$acc loop vector
        do i = 1 , np 
          !T_vadv(i,j  ,1,ie) = facp(i,j,1)*( elem(ie)%state%T(i,j  ,2,tl) - elem(ie)%state%T(i,j  ,1,tl) )
          !v_vadv(i,j,1,1,ie) = facp(i,j,1)*( elem(ie)%state%v(i,j,1,2,tl) - elem(ie)%state%v(i,j,1,1,tl) )
          !v_vadv(i,j,2,1,ie) = facp(i,j,1)*( elem(ie)%state%v(i,j,2,2,tl) - elem(ie)%state%v(i,j,2,1,tl) )
          T_vadv(i,j  ,1,ie) = facp(i,j,1)*( state_T(i,j  ,2,tl,ie) - state_T(i,j  ,1,tl,ie) )
          v_vadv(i,j,1,1,ie) = facp(i,j,1)*( state_v(i,j,1,2,tl,ie) - state_v(i,j,1,1,tl,ie) )
          v_vadv(i,j,2,1,ie) = facp(i,j,1)*( state_v(i,j,2,2,tl,ie) - state_v(i,j,2,1,tl,ie) )
          do k = 2 , nlev-1
            !T_vadv(i,j  ,k,ie) = facp(i,j,k)*( elem(ie)%state%T(i,j  ,k+1,tl) - elem(ie)%state%T(i,j  ,k,tl) ) + facm(i,j,k)*( elem(ie)%state%T(i,j  ,k,tl) - elem(ie)%state%T(i,j  ,k-1,tl) )
            !v_vadv(i,j,1,k,ie) = facp(i,j,k)*( elem(ie)%state%v(i,j,1,k+1,tl) - elem(ie)%state%v(i,j,1,k,tl) ) + facm(i,j,k)*( elem(ie)%state%v(i,j,1,k,tl) - elem(ie)%state%v(i,j,1,k-1,tl) )
            !v_vadv(i,j,2,k,ie) = facp(i,j,k)*( elem(ie)%state%v(i,j,2,k+1,tl) - elem(ie)%state%v(i,j,2,k,tl) ) + facm(i,j,k)*( elem(ie)%state%v(i,j,2,k,tl) - elem(ie)%state%v(i,j,2,k-1,tl) )
            T_vadv(i,j  ,k,ie) = facp(i,j,k)*( state_T(i,j  ,k+1,tl,ie) - state_T(i,j  ,k,tl,ie) ) + facm(i,j,k)*( state_T(i,j  ,k,tl,ie) - state_T(i,j  ,k-1,tl,ie) )
            v_vadv(i,j,1,k,ie) = facp(i,j,k)*( state_v(i,j,1,k+1,tl,ie) - state_v(i,j,1,k,tl,ie) ) + facm(i,j,k)*( state_v(i,j,1,k,tl,ie) - state_v(i,j,1,k-1,tl,ie) )
            v_vadv(i,j,2,k,ie) = facp(i,j,k)*( state_v(i,j,2,k+1,tl,ie) - state_v(i,j,2,k,tl,ie) ) + facm(i,j,k)*( state_v(i,j,2,k,tl,ie) - state_v(i,j,2,k-1,tl,ie) )
          enddo
          !T_vadv(i,j  ,nlev,ie) = facm(i,j,nlev)*( elem(ie)%state%T(i,j  ,nlev,tl) - elem(ie)%state%T(i,j  ,nlev-1,tl) )
          !v_vadv(i,j,1,nlev,ie) = facm(i,j,nlev)*( elem(ie)%state%v(i,j,1,nlev,tl) - elem(ie)%state%v(i,j,1,nlev-1,tl) )
          !v_vadv(i,j,2,nlev,ie) = facm(i,j,nlev)*( elem(ie)%state%v(i,j,2,nlev,tl) - elem(ie)%state%v(i,j,2,nlev-1,tl) )
          T_vadv(i,j  ,nlev,ie) = facm(i,j,nlev)*( state_T(i,j  ,nlev,tl,ie) - state_T(i,j  ,nlev-1,tl,ie) )
          v_vadv(i,j,1,nlev,ie) = facm(i,j,nlev)*( state_v(i,j,1,nlev,tl,ie) - state_v(i,j,1,nlev-1,tl,ie) )
          v_vadv(i,j,2,nlev,ie) = facm(i,j,nlev)*( state_v(i,j,2,nlev,tl,ie) - state_v(i,j,2,nlev-1,tl,ie) )
        enddo
      enddo
    enddo
  end subroutine preq_vertadv_openacc



  subroutine preq_hydrostatic_openacc(elem,T_v,p,dp,nets,nete,asyncid_in)
    use kinds             , only : real_kind
    use dimensions_mod    , only : np, nlev, nelemd
    use element_mod       , only : element_t
    use physical_constants, only : rgas
    use openacc_utils_mod, only: acc_async_sync
    implicit none
    type(element_t)     , intent(inout) :: elem(:)
    real(kind=real_kind), intent(in)    :: T_v(np,np,nlev,nelemd)
    real(kind=real_kind), intent(in)    :: p  (np,np,nlev,nelemd)
    real(kind=real_kind), intent(in)    :: dp (np,np,nlev,nelemd)
    integer             , intent(in)    :: nets,nete
    integer, optional   , intent(in)    :: asyncid_in
    !---------------------------Local workspace-----------------------------
    integer i,j,k,ie                                       ! longitude, level indices
    real(kind=real_kind) :: phii    ! Geopotential at interfaces
    real(kind=real_kind) :: tmp
    integer :: asyncid
    if (present(asyncid_in)) then
      asyncid = asyncid_in
    else
      asyncid = acc_async_sync
    endif
    !$acc parallel loop gang vector collapse(3) private(tmp,phii) present(elem,t_v,p,dp) async(asyncid)
    do ie = nets , nete
      do j = 1 , np
        do i = 1 , np
          tmp = Rgas*T_v(i,j,nlev,ie)*dp(i,j,nlev,ie)*0.5d0/p(i,j,nlev,ie)
          elem(ie)%derived%phi(i,j,nlev) = elem(ie)%state%phis(i,j) + tmp
          phii = tmp*2
          do k = nlev-1 , 2 , -1
            tmp = Rgas*T_v(i,j,k,ie)*dp(i,j,k,ie)*0.5d0/p(i,j,k,ie)
            elem(ie)%derived%phi(i,j,k) = elem(ie)%state%phis(i,j) + phii + tmp
            phii = phii + tmp*2
          enddo
          tmp = Rgas*T_v(i,j,1,ie)*dp(i,j,1,ie)*0.5d0/p(i,j,1,ie)
          elem(ie)%derived%phi(i,j,1) = elem(ie)%state%phis(i,j) + phii + tmp
        enddo
      enddo
    enddo
  end subroutine preq_hydrostatic_openacc



  subroutine preq_omega_ps_openacc(omega_p,p,vgrad_p,divdp,nets,nete,asyncid_in)
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nelemd
    use openacc_utils_mod, only: acc_async_sync
    implicit none
    real(kind=real_kind), intent(in) :: divdp  (np,np,nlev,nelemd)      ! divergence
    real(kind=real_kind), intent(in) :: vgrad_p(np,np,nlev,nelemd) ! v.grad(p)
    real(kind=real_kind), intent(in) :: p      (np,np,nlev,nelemd)     ! layer thicknesses (pressure)
    real(kind=real_kind), intent(out):: omega_p(np,np,nlev,nelemd)   ! vertical pressure velocity
    integer             , intent(in) :: nets , nete
    integer, optional   , intent(in) :: asyncid_in
    !---------------------------Local workspace-----------------------------
    integer :: i,j,k,ie                   ! longitude, level indices
    real(kind=real_kind) :: suml      ! partial sum over l = (1, k-1)
    integer :: asyncid
    if (present(asyncid_in)) then
      asyncid = asyncid_in
    else
      asyncid = acc_async_sync
    endif
    !$acc parallel loop gang vector collapse(3) private(suml) present(divdp,vgrad_p,p,omega_p) async(asyncid)
    do ie = nets , nete
      do j = 1 , np
        do i = 1 , np
          omega_p(i,j,1,ie) = ( vgrad_p(i,j,1,ie) - 0.5d0*divdp(i,j,1,ie) ) / p(i,j,1,ie)
          suml = divdp(i,j,1,ie)
          do k = 2 , nlev-1
            omega_p(i,j,k,ie) = ( ( vgrad_p(i,j,k,ie) - 0.5d0*divdp(i,j,k,ie) ) - suml ) / p(i,j,k,ie)
            suml = suml + divdp(i,j,k,ie)
          enddo
          omega_p(i,j,nlev,ie) = ( ( vgrad_p(i,j,nlev,ie) - 0.5d0*divdp(i,j,nlev,ie) ) - suml ) / p(i,j,nlev,ie)
        enddo
      enddo
    enddo
  end subroutine preq_omega_ps_openacc



#endif
end module prim_si_mod





#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_driver_mod

  use prim_driver_base, only:&
    apply_forcing, &
    flt_advection, &
    get_diagnostic_state, &
    get_subcycle_timesteps, &
    leapfrog_bootstrap, &
    prim_energy_fixer, &
    prim_init1,&
    prim_init2 ,&
    prim_finalize,&
    prim_run,&
    prim_run_subcycle_diags

  use control_mod,        only: statefreq, energy_fixer, ftype, qsplit, rsplit, test_cfldep, disable_diagnostics
  use dimensions_mod,     only: np, nlev, nlevp, nelem, nelemd, qsize, nc, ntrac
  use element_mod,        only: element_t, timelevels,  allocate_element_desc
  use fvm_control_volume_mod, only: n0_fvm, fvm_struct
  use hybrid_mod,         only: hybrid_t
  use hybvcoord_mod,      only: hvcoord_t
  use kinds,              only: real_kind, iulog
  use parallel_mod,       only: abortmp
  use perf_mod,           only: t_startf, t_stopf
  use prim_advance_mod,   only: applycamforcing, applycamforcing_dynamics, compute_and_apply_rhs  
  use prim_state_mod,     only: prim_printstate, prim_diag_scalars, prim_energy_halftimes
  use reduction_mod,      only: parallelmax
  use test_mod,           only: apply_test_forcing
  use time_mod,           only: timeLevel_t, timelevel_update, timelevel_qdp, nsplit

  implicit none
  contains

! customized: prim_run_subcycle , prim_step

!_______________________________________________________________________
subroutine prim_run_subcycle(elem,fvm,hybrid,nets,nete,dt,tl,hvcoord,nsubstep)

    type (element_t) ,    intent(inout) :: elem(:)
    type (fvm_struct),    intent(inout) :: fvm(:)
    type (hybrid_t),      intent(in)    :: hybrid                       ! distributed parallel structure (shared)
    type (hvcoord_t),     intent(inout) :: hvcoord                      ! hybrid vertical coordinate struct
    integer,              intent(in)    :: nets                         ! starting thread element number (private)
    integer,              intent(in)    :: nete                         ! ending thread element number   (private)
    real(kind=real_kind), intent(in)    :: dt                           ! "timestep dependent" timestep
    type (TimeLevel_t),   intent(inout) :: tl
    integer,              intent(in)    :: nsubstep                     ! nsubstep = 1 .. nsplit

    logical :: compute_diagnostics, compute_energy
    real(kind=real_kind) :: dt_q, dt_remap
    integer :: n0_qdp,np1_qdp,r,nstep_end

    ! get dt_q, dt_remp, and nstep_end from dt
    call get_subcycle_timesteps(dt_q,dt_remap,nstep_end,dt,tl)

    ! enable/disable diagnositics for this time-step
    call get_diagnostic_state(compute_diagnostics,compute_energy,tl,nstep_end)

    ! compute scalar diagnostics if currently active
    if (compute_diagnostics) call prim_diag_scalars(elem,hvcoord,tl,4,.true.,nets,nete)

    ! apply cam forcing or stand-alone test forcing to state variables
    call apply_forcing(elem,fvm,hybrid,hvcoord,tl,dt_remap,nets,nete)

    ! get E(1): energy after CAM forcing
    if (compute_energy) call prim_energy_halftimes(elem,hvcoord,tl,1,.true.,nets,nete)

    ! get qmass and variance, using Q(n0),Qdp(n0)
    if (compute_diagnostics) call prim_diag_scalars(elem,hvcoord,tl,1,.true.,nets,nete)

    ! loop over rsplit vertically lagrangian timesteps
    call prim_step(elem, fvm, hybrid,nets,nete, dt, tl, hvcoord,compute_diagnostics,1)
    do r=2,rsplit
       call TimeLevel_update(tl,"leapfrog")
       call prim_step(elem, fvm, hybrid,nets,nete, dt, tl, hvcoord,.false.,r)
    enddo

    ! update diagnostic variables
    call prim_run_subcycle_diags(elem,hvcoord,tl,nets,nete)

    ! apply energy fixer, if active
    if (compute_diagnostics)  call prim_diag_scalars    (elem,hvcoord,tl,2,.false.,nets,nete)
    if (compute_energy)       call prim_energy_halftimes(elem,hvcoord,tl,2,.false.,nets,nete)
    if (energy_fixer > 0)     call prim_energy_fixer    (elem,hvcoord,hybrid,tl,nets,nete,nsubstep)
    if (compute_diagnostics)  call prim_diag_scalars    (elem,hvcoord,tl,3,.false.,nets,nete)
    if (compute_energy)       call prim_energy_halftimes(elem,hvcoord,tl,3,.false.,nets,nete)

    ! update dynamics time level pointers
    call TimeLevel_update(tl,"leapfrog")

    ! print some diagnostic information
    if (compute_diagnostics) call prim_printstate(elem, tl, hybrid,hvcoord,nets,nete, fvm)

  end subroutine prim_run_subcycle

!_______________________________________________________________________
subroutine prim_step(elem, fvm, hybrid,nets,nete,dt,tl,hvcoord,compute_diagnostics,rstep)

    use control_mod,        only: statefreq, integration, ftype, qsplit, nu_p, test_cfldep, rsplit
    use control_mod,        only: use_semi_lagrange_transport, tracer_transport_type
    use control_mod,        only: tracer_grid_type, TRACER_GRIDTYPE_GLL
    use derivative_mod,     only: subcell_integration
    use fvm_bsp_mod,        only: get_boomerang_velocities_gll, get_solidbody_velocities_gll
    use fvm_control_volume_mod, only : fvm_supercycling
    use fvm_mod,            only: fvm_ideal_test, IDEAL_TEST_OFF, IDEAL_TEST_ANALYTICAL_WINDS
    use fvm_mod,            only: fvm_test_type, IDEAL_TEST_BOOMERANG, IDEAL_TEST_SOLIDBODY
    use hybvcoord_mod,      only : hvcoord_t
    use parallel_mod,       only: abortmp
    use prim_advance_mod,   only: prim_advance_exp, overwrite_SEdensity
    use prim_advection_mod, only: prim_advec_tracers_fvm, prim_advec_tracers_remap, deriv
    use reduction_mod,      only: parallelmax
    use time_mod,           only: time_at,TimeLevel_t, timelevel_update, nsplit

    type(element_t),      intent(inout) :: elem(:)
    type(fvm_struct),     intent(inout) :: fvm(:)
    type(hybrid_t),       intent(in)    :: hybrid   ! distributed parallel structure (shared)
    type(hvcoord_t),      intent(inout) :: hvcoord  ! hybrid vertical coordinate struct
    integer,              intent(in)    :: nets     ! starting thread element number (private)
    integer,              intent(in)    :: nete     ! ending thread element number   (private)
    real(kind=real_kind), intent(in)    :: dt       ! "timestep dependent" timestep
    type(TimeLevel_t),    intent(inout) :: tl
    integer,              intent(in)    :: rstep    ! vertical remap subcycling step
    logical,              intent(in)    :: compute_diagnostics

    real(kind=real_kind) :: st, st1, dp, dt_q
    integer :: ie, t, q,k,i,j,n, n_Q, qn0
integer :: nm1,n0,np1,nstep
real (kind=real_kind) ::  eta_ave_w


    real (kind=real_kind) :: maxcflx, maxcfly
    real (kind=real_kind) ::  tempdp3d(np,np), x
    real (kind=real_kind) ::  tempmass(nc,nc)
    real (kind=real_kind) ::  tempflux(nc,nc,4)
    real (kind=real_kind) :: dp_np1(np,np)

    ! Clear derived quantities

    call t_startf("prim_step_init")
    dt_q = dt*qsplit
    do ie=nets,nete
      elem(ie)%derived%eta_dot_dpdn = 0         ! mean vertical mass flux
      elem(ie)%derived%vn0                  = 0 ! mean horizontal mass flux
      elem(ie)%derived%omega_p              = 0
      if (nu_p>0) then
         elem(ie)%derived%dpdiss_ave        = 0
         elem(ie)%derived%dpdiss_biharmonic = 0
      endif
    enddo
    call t_stopf("prim_step_init")

    ! Take qsplit dynamics steps

    call t_startf("prim_step_dyn")
    n_Q = tl%n0  ! n_Q = timelevel of FV tracers at time t. 

    qn0=n_Q
    nm1   = tl%nm1
    n0    = tl%n0
    np1   = tl%np1
    nstep = tl%nstep
    eta_ave_w = 1.0d0

    call t_startf("RK3_timestep")
       call compute_and_apply_rhs(np1,n0,n0,qn0,dt/3,elem,hvcoord,hybrid,deriv(hybrid%ithr),nets,nete,compute_diagnostics,0d0)
       ! u2 = u0 + dt/2 RHS(u1)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt/2,elem,hvcoord,hybrid,deriv(hybrid%ithr),nets,nete,.false.,0d0)
       ! u3 = u0 + dt RHS(u2)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,deriv(hybrid%ithr),nets,nete,.false.,eta_ave_w)
       call t_stopf("RK3_timestep")
    !call compute_and_apply_rhs(np1,nm1,n0,n_Q,dt,elem,hvcoord,hybrid, deriv(hybrid%ithr),nets,nete,compute_diagnostics,eta_ave_w)
    !call prim_advance_exp(elem, deriv(hybrid%ithr), hvcoord, hybrid, dt, tl, nets, nete, compute_diagnostics)
    call t_stopf("prim_step_dyn")

    ! Advect tracers

  end subroutine prim_step

end module


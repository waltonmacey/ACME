!This is where all of the PGI CUDA FORTRAN code will go, and these routines will be called from prim_advection_mod.
!This is compiled regardless, but PGI-specific calls are always wrapped in the USE_CUDA_FORTRAN ifdefs that are automagically
!activated when -Mcuda is specified during compilation with a PGI compiler. Thus, it will be ignored unless explicitly
!activated by the user
!
!As a general rule, all of the routines in here will be called within a threaded context (assuming ELEMENT_OPENMP is not
!deifned), and therefore, we enforce BARRIERS, MASTERS, and SINGLES from within these routines rather than outside them.
!This is to minimize the visible code impacts on the existing CPU code.

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!#define CUDA_DEBUG

#ifdef CUDA_DEBUG
#define _CHECK(line) ierr = cudaThreadSynchronize(); if (ierr .ne. 0) stop line
#else
#define _CHECK(line)
#endif


module cuda_mod
#if USE_CUDA_FORTRAN

! NP > 4 is not supported due to shared memory constraints
#if NP > 4
#error CUDA Fortran build only supported with NP <= 4
#endif

#define PAD 1

!Put everything CUDA-specific in here so it doesn't get compiled without -Mcuda enabled on a PGI compiler
  use cudafor
  use kinds          , only: real_kind
  use dimensions_mod , only: np,nlevp,nlev,qsize,qsize_d,max_corner_elem,max_neigh_edges,nelemd
  use element_mod    , only: timelevels
  use edge_mod       , only: EdgeBuffer_t
  implicit none
  private
  save

  !First listed are all externally accescible routines
  public :: cuda_mod_init
  public :: euler_step_cuda 
  public :: qdp_time_avg_cuda
  public :: advance_hypervis_scalar_cuda
  public :: copy_qdp_d2h
  public :: copy_qdp_h2d
  public :: cudaThreadSynchronize_wrap
  public :: precompute_divdp_cuda

  integer,parameter :: numk_eul = 6
  integer,parameter :: numk_hyp = 6
  integer,parameter :: numk_lim2d = 15

  !This is from prim_advection_mod.F90
  type(EdgeBuffer_t) :: edgeAdv, edgeAdv1, edgeAdvQ3, edgeAdvQ2, edgeAdv_p1, edgeAdvDSS
  integer,parameter :: DSSeta = 1
  integer,parameter :: DSSomega = 2
  integer,parameter :: DSSdiv_vdp_ave = 3
  integer,parameter :: DSSno_var = -1

  !Device arrays
  real (kind=real_kind),device,allocatable,dimension(:,:,:,:,:,:) :: qdp_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:,:)     :: qdp1_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:,:,:)   :: qtens_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:,:,:)   :: qtens_biharmonic_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:)       :: spheremp_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:,:)     :: dpdiss_ave_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:)       :: rspheremp_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:,:,:)   :: dinv_d
  real (kind=real_kind),device,allocatable,dimension(:,:)         :: deriv_dvv_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:)       :: variable_hyperviscosity_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:)       :: metdet_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:,:)     :: dpdiss_biharmonic_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:)       :: rmetdet_d
  real (kind=real_kind),device,allocatable,dimension(:,:)         :: edgebuf_d
  real (kind=real_kind),device,allocatable,dimension(:,:)         :: edgebuf_Q2_d
  real (kind=real_kind),device,allocatable,dimension(:,:)         :: edgebuf_Q3_d
  real (kind=real_kind),device,allocatable,dimension(:,:)         :: recvbuf_Q2_d
  real (kind=real_kind),device,allocatable,dimension(:,:)         :: recvbuf_Q3_d
  real (kind=real_kind),device,allocatable,dimension(:)           :: hyai_d
  real (kind=real_kind),device,allocatable,dimension(:)           :: hybi_d
  logical              ,device,allocatable,dimension(:,:)         :: reverse_d
  integer              ,device,allocatable,dimension(:,:)         :: putmapP_d
  integer              ,device,allocatable,dimension(:,:)         :: getmapP_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:,:,:)   :: vstar_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:,:)     :: divdp_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:,:)     :: divdp_proj_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:,:)     :: dp_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:,:)     :: dp_star_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:)       :: qmin_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:)       :: qmax_d
  integer              ,device,allocatable,dimension(:)           :: recv_internal_indices_d
  integer              ,device,allocatable,dimension(:)           :: recv_external_indices_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:)       :: recvbuf_d
  real(kind=real_kind) ,device,allocatable,dimension(:,:,:,:)     :: dpo_d
  real(kind=real_kind) ,device,allocatable,dimension(:,:,:,:,:)   :: ppmdx_d
  real(kind=real_kind) ,device,allocatable,dimension(:,:,:,:)     :: z1_d
  real(kind=real_kind) ,device,allocatable,dimension(:,:,:,:)     :: z2_d
  integer              ,device,allocatable,dimension(:,:,:,:)     :: kid_d
  real(kind=real_kind) ,device,allocatable,dimension(:)           :: dp0_d
  real(kind=real_kind) ,device,allocatable,dimension(:,:,:,:,:)   :: qmin_exch
  real(kind=real_kind) ,device,allocatable,dimension(:,:,:,:,:)   :: qmax_exch

  !PINNED Host arrays
  real(kind=real_kind) ,pinned,allocatable,dimension(:,:,:)       :: spheremp_h
  real(kind=real_kind) ,pinned,allocatable,dimension(:,:,:,:,:,:) :: qdp_h
  real(kind=real_kind) ,pinned,allocatable,dimension(:,:,:,:)     :: qdp1_h
  real(kind=real_kind) ,pinned,allocatable,dimension(:,:,:,:,:)   :: vstar_h
  real(kind=real_kind) ,pinned,allocatable,dimension(:,:,:,:,:)   :: qtens_h
  real(kind=real_kind) ,pinned,allocatable,dimension(:,:,:,:)     :: dp_h
  real(kind=real_kind) ,pinned,allocatable,dimension(:,:,:,:)     :: divdp_h
  real(kind=real_kind) ,pinned,allocatable,dimension(:,:,:,:)     :: divdp_proj_h
  real(kind=real_kind) ,pinned,allocatable,dimension(:,:,:,:)     :: dp_star_h
  real(kind=real_kind) ,pinned,allocatable,dimension(:,:,:,:)     :: dp_np1_h
  real(kind=real_kind) ,pinned,allocatable,dimension(:,:,:,:)     :: dpdiss_ave_h
  real(kind=real_kind) ,pinned,allocatable,dimension(:,:,:)       :: sendbuf_h
  real(kind=real_kind) ,pinned,allocatable,dimension(:,:,:)       :: recvbuf_h
  real(kind=real_kind) ,pinned,allocatable,dimension(:,:,:,:,:)   :: qtens_biharmonic_h
  real(kind=real_kind) ,pinned,allocatable,dimension(:,:,:)       :: variable_hyperviscosity_h
  real(kind=real_kind) ,pinned,allocatable,dimension(:,:,:)       :: qmin
  real(kind=real_kind) ,pinned,allocatable,dimension(:,:,:)       :: qmax
  real(kind=real_kind) ,pinned,allocatable,dimension(:,:,:,:)     :: dpdiss_biharmonic_h
  real(kind=real_kind) ,pinned,allocatable,dimension(:,:,:,:)     :: dpo_h
  real(kind=real_kind) ,pinned,allocatable,dimension(:,:,:,:,:)   :: ppmdx_h
  real(kind=real_kind) ,pinned,allocatable,dimension(:,:,:,:)     :: z1_h
  real(kind=real_kind) ,pinned,allocatable,dimension(:,:,:,:)     :: z2_h
  real(kind=real_kind) ,pinned,allocatable,dimension(:,:)         :: edgebuf_Q2
  integer              ,pinned,allocatable,dimension(:,:,:,:)     :: kid_h
  real(kind=real_kind) ,pinned,allocatable,dimension(:)           :: dp0_h
  logical              ,pinned,allocatable,dimension(:,:)         :: reverse_h
  integer              ,pinned,allocatable,dimension(:,:)         :: putmapP_h
  integer              ,pinned,allocatable,dimension(:,:)         :: getmapP_h
  real (kind=real_kind),pinned,allocatable,dimension(:,:)         :: edgebuf_Q2_h
  real (kind=real_kind),pinned,allocatable,dimension(:,:)         :: edgebuf_Q3_h
  real (kind=real_kind),pinned,allocatable,dimension(:,:)         :: recvbuf_Q2_h
  real (kind=real_kind),pinned,allocatable,dimension(:,:)         :: recvbuf_Q3_h

  !OpenACC data structures
  real (kind=real_kind),allocatable,dimension(:,:)         :: deriv_dvv_acc
  real (kind=real_kind),allocatable,dimension(:,:,:,:,:)   :: dinv_acc
  real (kind=real_kind),allocatable,dimension(:,:,:,:,:)   :: tensorvisc_acc

  !Normal Host arrays
  integer,allocatable,dimension(:)   :: send_nelem
  integer,allocatable,dimension(:)   :: recv_nelem
  integer,allocatable,dimension(:,:) :: send_indices
  integer,allocatable,dimension(:,:) :: recv_indices
  integer,allocatable,dimension(:)   :: recv_internal_indices
  integer,allocatable,dimension(:)   :: recv_external_indices
  integer :: recv_external_nelem
  integer :: recv_internal_nelem
  logical :: old_peu
  logical,allocatable,dimension(:)   :: d2h_done
  logical,allocatable,dimension(:)   :: msg_sent
  logical,allocatable,dimension(:)   :: msg_rcvd
  logical,allocatable,dimension(:)   :: h2d_done
  integer, parameter :: south_px = 1
  integer, parameter :: east_px  = 3
  integer, parameter :: north_px = 4
  integer, parameter :: west_px  = 2
  integer :: cuda_streams
  integer(kind=cuda_stream_kind), allocatable :: streams(:)
  integer(kind=cuda_stream_kind), allocatable :: streams2(:)
  integer            :: nbuf
  integer            :: nmsg_rcvd
  integer            :: nmsg_sent
  integer, device :: max_neigh_edges_d, max_corner_elem_d

  type(cudaEvent) :: timer1, timer2
  integer, parameter :: gs = 2
  real(kind=real_kind), device :: rrearth_d


contains



  !The point of this is to initialize any data required in other routines of this module as well
  !as to run one initial CUDA kernel just to get those overheads out of the way so that subsequent
  !timing routines are accurage.
  subroutine cuda_mod_init(elem,deriv,hvcoord)
    use edge_mod      , only: initEdgeBuffer
    use schedule_mod  , only: schedule_t, cycle_t, schedule
    use edge_mod      , only: Edgebuffer_t
    use element_mod   , only: element_t
    use derivative_mod, only: derivative_t
    use hybvcoord_mod, only: hvcoord_t
    use physical_constants, only: rrearth
    implicit none
    type(element_t)   , intent(in) :: elem(:)
    type(derivative_t), intent(in) :: deriv
    type(hvcoord_t)   , intent(in) :: hvcoord

    type (Cycle_t),pointer    :: pCycle
    type (Schedule_t),pointer :: pSchedule
    integer                   :: ie , ierr , icycle , iPtr , rank , nSendCycles , nRecvCycles , nlyr , mx_send_len , mx_recv_len , n, k
    real(kind=real_kind)      :: dinv_t(np , np , 2 , 2)
    type (dim3)               :: griddim , blockdim
    logical,allocatable,dimension(:,:) :: send_elem_mask
    logical,allocatable,dimension(:,:) :: recv_elem_mask
    logical,allocatable,dimension(:)   :: elem_computed
    integer :: total_work
    real(kind=real_kind), pointer :: buf_ptr(:) => null()
    real(kind=real_kind), pointer :: receive_ptr(:) => null()

    rrearth_d = rrearth

    max_neigh_edges_d = max_neigh_edges
    max_corner_elem_d = max_corner_elem

    cuda_streams = nelemd*8
    allocate(streams (0:cuda_streams))
    allocate(streams2(0:cuda_streams))

#if (defined ELEMENT_OPENMP)
    write(*,*) 'ERROR: Do not use ELEMENT_OPENMP and CUDA FORTRAN'
    stop
#endif
!$OMP BARRIER
!$OMP MASTER
    write(*,*) "cuda_mod_init"

    write(*,*) "allocate arrays on device & host"
#ifdef _PREDICT
    pSchedule => Schedule(iam)
#else
    pSchedule => Schedule(1)
#endif
    nlyr = edgeAdv%nlyr
    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles
    mx_send_len = 0
    mx_recv_len = 0
    do icycle=1,nSendCycles
      if (pSchedule%SendCycle(icycle)%lengthP > mx_send_len) mx_send_len = pSchedule%SendCycle(icycle)%lengthP
    enddo
    do icycle=1,nRecvCycles
      if (pSchedule%RecvCycle(icycle)%lengthP > mx_recv_len) mx_recv_len = pSchedule%RecvCycle(icycle)%lengthP 
    enddo
    nbuf=4*(np+max_corner_elem)*nelemd  !inlined from edge_mod.F90, initEdgeBuffer()

    !Allocate the host and device arrays
    allocate( qmin_d                   (nlev,qsize_d                 ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( qmax_d                   (nlev,qsize_d                 ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( qdp_d                    (np,np,nlev,qsize_d,timelevels,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( qdp1_d                   (np,np,nlev                   ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( dpdiss_ave_d             (np,np,nlev                   ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( qtens_d                  (np,np,nlev,qsize_d           ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( qtens_biharmonic_d       (np,np,nlev,qsize_d           ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( spheremp_d               (np,np                        ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( rspheremp_d              (np,np                        ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( dinv_d                   (np,np,2,2                    ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( deriv_dvv_d              (np,np                               ) , stat = ierr ); _CHECK(__LINE__)
    allocate( variable_hyperviscosity_d(np,np                        ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( metdet_d                 (np,np                        ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( rmetdet_d                (np,np                        ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( dpdiss_biharmonic_d      (np,np,nlev                   ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( vstar_d                  (np,np,nlev,2                 ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( divdp_d                  (np,np,nlev                   ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( divdp_proj_d             (np,np,nlev                   ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( dp_d                     (np,np,nlev                   ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( hyai_d                   (      nlev+1                        ) , stat = ierr ); _CHECK(__LINE__)
    allocate( hybi_d                   (      nlev+1                        ) , stat = ierr ); _CHECK(__LINE__)
    allocate( dp_star_d                (np,np,nlev                   ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( reverse_d                (max_neigh_edges              ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( putmapP_d                (max_neigh_edges              ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( getmapP_d                (max_neigh_edges              ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( edgebuf_d                (nlev*qsize_d,nbuf                   ) , stat = ierr ); _CHECK(__LINE__)
    allocate( edgebuf_Q2_d             (nlev*qsize_d*2,nbuf                 ) , stat = ierr ); _CHECK(__LINE__)
    allocate( edgebuf_Q3_d             (nlev*qsize_d*3,nbuf                 ) , stat = ierr ); _CHECK(__LINE__)
    allocate( recvbuf_Q2_d             (nlev*qsize_d*2,nbuf                 ) , stat = ierr ); _CHECK(__LINE__)
    allocate( recvbuf_Q3_d             (nlev*qsize_d*3,nbuf                 ) , stat = ierr ); _CHECK(__LINE__)
    allocate( recv_internal_indices_d  (nelemd                              ) , stat = ierr ); _CHECK(__LINE__)
    allocate( recv_external_indices_d  (nelemd                              ) , stat = ierr ); _CHECK(__LINE__)
    allocate( dpo_d                    (np,np,   1-gs:nlev+gs        ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( ppmdx_d                  (np,np,10,   0:nlev+1         ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( z1_d                     (np,np,        nlev           ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( z2_d                     (np,np,        nlev           ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( kid_d                    (np,np,        nlev           ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( dp0_d                    (              nlev                  ) , stat = ierr ); _CHECK(__LINE__)

    allocate( spheremp_h               (np,np                        ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( qdp_h                    (np,np,nlev,qsize_d,timelevels,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( qdp1_h                   (np,np,nlev                   ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( dpdiss_ave_h             (np,np,nlev                   ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( qmin                     (nlev,qsize_d                 ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( qmax                     (nlev,qsize_d                 ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( qmin_exch                (np,np,nlev,qsize_d           ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( qmax_exch                (np,np,nlev,qsize_d           ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( vstar_h                  (np,np,nlev,2                 ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( qtens_h                  (np,np,nlev,qsize_d           ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( qtens_biharmonic_h       (np,np,nlev,qsize_d           ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( variable_hyperviscosity_h(np,np                        ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( dp_h                     (np,np,nlev                   ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( divdp_h                  (np,np,nlev                   ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( divdp_proj_h             (np,np,nlev                   ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( dp_star_h                (np,np,nlev                   ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( dpdiss_biharmonic_h      (np,np,nlev                   ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( dp_np1_h                 (np,np,nlev                   ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( dpo_h                    (np,np,   1-gs:nlev+gs        ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( ppmdx_h                  (np,np,10,   0:nlev+1         ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( z1_h                     (np,np,        nlev           ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( z2_h                     (np,np,        nlev           ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( kid_h                    (np,np,        nlev           ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( dp0_h                    (              nlev                  ) , stat = ierr ); _CHECK(__LINE__)
    allocate( edgebuf_Q2               (nlev*qsize_d*2,nbuf                 ) , stat = ierr ); _CHECK(__LINE__)
    allocate( reverse_h                (max_neigh_edges              ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( putmapP_h                (max_neigh_edges              ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( getmapP_h                (max_neigh_edges              ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( edgebuf_Q2_h             (nlev*qsize_d*2,nbuf                 ) , stat = ierr ); _CHECK(__LINE__)
    allocate( edgebuf_Q3_h             (nlev*qsize_d*3,nbuf                 ) , stat = ierr ); _CHECK(__LINE__)
    allocate( recvbuf_Q2_h             (nlev*qsize_d*2,nbuf                 ) , stat = ierr ); _CHECK(__LINE__)
    allocate( recvbuf_Q3_h             (nlev*qsize_d*3,nbuf                 ) , stat = ierr ); _CHECK(__LINE__)

    !OpenACC arrays
    allocate( deriv_dvv_acc            (np,np                               ) , stat = ierr ); _CHECK(__LINE__)
    allocate( dinv_acc                 (np,np,2,2                    ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( tensorvisc_acc           (2,2,np,np                    ,nelemd) , stat = ierr ); _CHECK(__LINE__)

    ! The PGI compiler with cuda enabled errors when allocating arrays of zero
    !   size - here when using only one MPI task
    if (nSendCycles > 0) then
      allocate( sendbuf_h                (nlev*qsize_d,mx_send_len,nSendCycles) , stat = ierr ); _CHECK(__LINE__)
      allocate( send_elem_mask           (nelemd,nSendCycles                  ) , stat = ierr ); _CHECK(__LINE__)
      allocate( send_nelem               (       nSendCycles                  ) , stat = ierr ); _CHECK(__LINE__)
      allocate( send_indices             (nelemd,nSendCycles                  ) , stat = ierr ); _CHECK(__LINE__)
      allocate( d2h_done                 (nSendCycles                         ) , stat = ierr ); _CHECK(__LINE__)
      allocate( msg_sent                 (nSendCycles                         ) , stat = ierr ); _CHECK(__LINE__)
    endif

    if (nRecvCycles > 0) then
      allocate( recvbuf_h                (nlev*qsize_d,mx_recv_len,nRecvCycles) , stat = ierr ); _CHECK(__LINE__)
      allocate( recvbuf_d                (nlev*qsize_d,mx_recv_len,nRecvCycles) , stat = ierr ); _CHECK(__LINE__)
      allocate( recv_elem_mask           (nelemd,nRecvCycles                  ) , stat = ierr ); _CHECK(__LINE__)
      allocate( recv_nelem               (       nRecvCycles                  ) , stat = ierr ); _CHECK(__LINE__)
      allocate( recv_indices             (nelemd,nRecvCycles                  ) , stat = ierr ); _CHECK(__LINE__)
      allocate( h2d_done                 (nRecvCycles                         ) , stat = ierr ); _CHECK(__LINE__)
      allocate( msg_rcvd                 (nRecvCycles                         ) , stat = ierr ); _CHECK(__LINE__)
    endif 

    allocate( recv_internal_indices    (nelemd                                ) , stat = ierr ); _CHECK(__LINE__)
    allocate( recv_external_indices    (nelemd                                ) , stat = ierr ); _CHECK(__LINE__)
    allocate( elem_computed            (nelemd                                ) , stat = ierr ); _CHECK(__LINE__)
    write(*,*) "send data from host to device"
    !Copy over data to the device
    ierr = cudaMemcpy( deriv_dvv_d   , deriv%dvv    , size( deriv%dvv    ) , cudaMemcpyHostToDevice ); _CHECK(__LINE__)
    deriv_dvv_acc = deriv%dvv
    ierr = cudaMemcpy( hyai_d        , hvcoord%hyai , size( hvcoord%hyai ) , cudaMemcpyHostToDevice ); _CHECK(__LINE__)
    ierr = cudaMemcpy( hybi_d        , hvcoord%hybi , size( hvcoord%hybi ) , cudaMemcpyHostToDevice ); _CHECK(__LINE__)
    do ie = 1,nelemd
      dinv_t(:,:,1,1) = elem(ie)%dinv(1,1,:,:)
      dinv_t(:,:,1,2) = elem(ie)%dinv(1,2,:,:)
      dinv_t(:,:,2,1) = elem(ie)%dinv(2,1,:,:)
      dinv_t(:,:,2,2) = elem(ie)%dinv(2,2,:,:)
      tensorvisc_acc(:,:,:,:,ie) = elem(ie)%tensorvisc
      ierr = cudaMemcpy( qdp_d                    (1,1,1,1,1,ie) , elem(ie)%state%Qdp               , size(elem(ie)%state%Qdp              ) , cudaMemcpyHostToDevice ); _CHECK(__LINE__)
      ierr = cudaMemcpy( spheremp_d               (1,1      ,ie) , elem(ie)%spheremp                , size(elem(ie)%spheremp               ) , cudaMemcpyHostToDevice ); _CHECK(__LINE__)
      ierr = cudaMemcpy( rspheremp_d              (1,1      ,ie) , elem(ie)%rspheremp               , size(elem(ie)%rspheremp              ) , cudaMemcpyHostToDevice ); _CHECK(__LINE__)
      ierr = cudaMemcpy( dinv_d                   (1,1,1,1  ,ie) , dinv_t                           , size(dinv_t                          ) , cudaMemcpyHostToDevice ); _CHECK(__LINE__)
      ierr = cudaMemcpy( metdet_d                 (1,1      ,ie) , elem(ie)%metdet                  , size(elem(ie)%metdet                 ) , cudaMemcpyHostToDevice ); _CHECK(__LINE__)
      ierr = cudaMemcpy( rmetdet_d                (1,1      ,ie) , elem(ie)%rmetdet                 , size(elem(ie)%rmetdet                ) , cudaMemcpyHostToDevice ); _CHECK(__LINE__)
      ierr = cudaMemcpy( putmapP_d                (1        ,ie) , elem(ie)%desc%putmapP            , size(elem(ie)%desc%putmapP           ) , cudaMemcpyHostToDevice ); _CHECK(__LINE__)
      ierr = cudaMemcpy( getmapP_d                (1        ,ie) , elem(ie)%desc%getmapP            , size(elem(ie)%desc%getmapP           ) , cudaMemcpyHostToDevice ); _CHECK(__LINE__)
      ierr = cudaMemcpy( reverse_d                (1        ,ie) , elem(ie)%desc%reverse            , size(elem(ie)%desc%reverse           ) , cudaMemcpyHostToDevice ); _CHECK(__LINE__)
      ierr = cudaMemcpy( variable_hyperviscosity_d(1,1      ,ie) , elem(ie)%variable_hyperviscosity , size(elem(ie)%variable_hyperviscosity) , cudaMemcpyHostToDevice ); _CHECK(__LINE__)
    enddo
    do k = 1 , nlev    
      dp0_h(k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                 ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0
    enddo
    ierr = cudaMemcpy( dp0_d , dp0_h , size(dp0_h) , cudaMemcpyHostToDevice ); _CHECK(__LINE__)

    ierr = cudaMemcpy( putmapP_h , putmapP_d , size(putmapP_h) , cudaMemcpyDeviceToHost ); _CHECK(__LINE__)
    ierr = cudaMemcpy( reverse_h , reverse_d , size(reverse_h) , cudaMemcpyDeviceToHost ); _CHECK(__LINE__)
    ierr = cudaMemcpy( getmapP_h , getmapP_d , size(getmapP_h) , cudaMemcpyDeviceToHost ); _CHECK(__LINE__)
    ierr = cudaMemcpy( variable_hyperviscosity_h , variable_hyperviscosity_d , size(variable_hyperviscosity_h) , cudaMemcpyDeviceToHost ); _CHECK(__LINE__)
    ierr = cudaMemcpy( spheremp_h , spheremp_d , size(spheremp_h) , cudaMemcpyDeviceToHost ); _CHECK(__LINE__)
    ierr = cudaMemcpy( dinv_acc , dinv_d , size(dinv_d) , cudaMemcpyDeviceToHost ); _CHECK(__LINE__)

    write(*,*) "edgebuffers"
    !These have to be in a threaded region or they complain and die
!$OMP END MASTER
    call initEdgeBuffer(edgeAdvDSS,      3*nlev                              )
    call initEdgeBuffer(edgeAdvQ3 ,max(nlev,qsize*nlev*3),buf_ptr,receive_ptr)  ! Qtens,Qmin, Qmax
    call initEdgeBuffer(edgeAdv1  ,nlev                  ,buf_ptr,receive_ptr)
    call initEdgeBuffer(edgeAdv   ,qsize*nlev            ,buf_ptr,receive_ptr)
    call initEdgeBuffer(edgeAdv_p1,qsize*nlev + nlev     ,buf_ptr,receive_ptr) 
    call initEdgeBuffer(edgeAdvQ2 ,qsize*nlev*2                              )  ! Qtens,Qmin, Qmax
    nullify(buf_ptr)
    nullify(receive_ptr)
!$OMP MASTER

    write(*,*) "initial kernel"
    !This needs to run because we need accurate timing during out cuda profiler runs. Initial kernel runs incur overheads, so we do this here
    blockdim = dim3(1,1,1)
    griddim  = dim3(1,1,1)
    call warmup <<< griddim , blockdim >>> ( ie ); _CHECK(__LINE__)
    ierr = cudaThreadSynchronize()

    do n = 0 , cuda_streams
      ierr = cudaStreamCreate(streams(n)); _CHECK(__LINE__)
      ierr = cudaStreamCreate(streams2(n)); _CHECK(__LINE__)
    enddo
    ierr = cudaDeviceSetCacheConfig(cudaFuncCachePreferShared); _CHECK(__LINE__)
    ierr = cudaEventCreate(timer1); _CHECK(__LINE__)
    ierr = cudaEventCreate(timer2); _CHECK(__LINE__)

    write(*,*) "Dividing elements among cycles in which they participate"
    !For efficient MPI, PCI-e, packing, and unpacking, we need to separate out the cycles by dependence. Once on cycle has packed, then stage the PCI-e D2H, MPI, PCI-e H2D, & internal unpack
    !We begin by testing what elements contribute to packing in what cycle's MPI data.
    do ie = 1,nelemd
      if (nRecvCycles > 0) recv_elem_mask(ie,:) = .false.
      do icycle = 1 , nRecvCycles
        do n = 1 , max_neigh_edges
          if ( elem(ie)%desc%getmapP(n) >= pSchedule%RecvCycle(icycle)%ptrP .and. &
               elem(ie)%desc%getmapP(n) <= pSchedule%RecvCycle(icycle)%ptrP + pSchedule%RecvCycle(icycle)%lengthP-1 ) then
            recv_elem_mask(ie,icycle) = .true.
          endif
        enddo
      enddo
    enddo
    edgebuf_d = 0.

    elem_computed = .false.
    !This pass accumulates for each cycle incides participating in the MPI_Irecv
    do icycle = 1 , nRecvCycles
      recv_nelem(icycle) = 0
      do ie = 1 , nelemd
        if ( recv_elem_mask(ie,icycle) .and. ( .not. elem_computed(ie) ) ) then
          recv_nelem(icycle) = recv_nelem(icycle) + 1
          recv_indices(recv_nelem(icycle),icycle) = ie
          elem_computed(ie) = .true.
        endif
      enddo
    enddo
    !This pass accumulates all elements from all cycles participating in MPI_Irecv into the recv_external_indices array
    recv_external_nelem = 0
    do icycle = 1 , nRecvCycles
      do ie = 1 , recv_nelem(icycle)
        recv_external_nelem = recv_external_nelem + 1
        recv_external_indices(recv_external_nelem) = recv_indices(ie,icycle)
      enddo
    enddo
    !This pass goes through all elements, and distributes evenly the elements not participating in MPI_Irecv 
    recv_internal_nelem = 0
    do ie = 1 , nelemd
      if ( .not. elem_computed(ie) ) then
        recv_internal_nelem = recv_internal_nelem + 1
        recv_internal_indices(recv_internal_nelem) = ie
      endif
    enddo

    write(*,*) "Sending element & cycle informationt to device"
    ierr = cudaMemcpy(recv_internal_indices_d,recv_internal_indices,size(recv_internal_indices),cudaMemcpyHostToDevice); _CHECK(__LINE__)
    ierr = cudaMemcpy(recv_external_indices_d,recv_external_indices,size(recv_external_indices),cudaMemcpyHostToDevice); _CHECK(__LINE__)

    write(*,*)"done cuda_mod_init"
    !$acc enter data pcreate(deriv_dvv_acc,variable_hyperviscosity_h,dinv_acc,spheremp_h,tensorvisc_acc)
    !$acc update device(deriv_dvv_acc,variable_hyperviscosity_h,dinv_acc,spheremp_h,tensorvisc_acc)
!$OMP END MASTER
!$OMP BARRIER
  end subroutine cuda_mod_init



  !Meaningless kernel just to get initial kernel overheads out of the way.
  attributes(global) subroutine warmup(a)
    integer,value :: a
    a = 2.0 * a
  end subroutine warmup



  function cudaThreadSynchronize_wrap() result(ierr)
    integer :: ierr
    ierr = cudaThreadSynchronize()
  end function cudaThreadSynchronize_wrap



  subroutine copy_qdp_h2d( elem , nt )
    use element_mod, only: element_t
    use perf_mod, only: t_startf, t_stopf
    implicit none
    type(element_t), intent(in   ) :: elem(:)
    integer        , intent(in   ) :: nt
    integer :: ierr , ie
    !$OMP BARRIER
    !$OMP MASTER
    ierr = cudaThreadSynchronize()
    call t_startf('CUDA QDP H2D')
    do ie = 1,nelemd
      ierr = cudaMemcpy( qdp_d(1,1,1,1,nt,ie) , elem(ie)%state%qdp(1,1,1,1,nt) , size(elem(ie)%state%qdp(:,:,:,:,nt)) , cudaMemcpyHostToDevice ); _CHECK(__LINE__)
    enddo
    ierr = cudaThreadSynchronize()
    call t_stopf('CUDA QDP H2D')
    !$OMP END MASTER
    !$OMP BARRIER
  end subroutine copy_qdp_h2d



  subroutine copy_qdp_d2h( elem , nt )
    use element_mod, only: element_t
    use perf_mod, only: t_startf, t_stopf
    implicit none
    type(element_t), intent(in   ) :: elem(:)
    integer        , intent(in   ) :: nt
    integer :: ierr , ie
    !$OMP BARRIER
    !$OMP MASTER
    ierr = cudaThreadSynchronize()
    call t_startf('CUDA QDP D2H')
    do ie = 1,nelemd
      ierr = cudaMemcpy( elem(ie)%state%qdp(1,1,1,1,nt) , qdp_d(1,1,1,1,nt,ie) , size(elem(ie)%state%qdp(:,:,:,:,nt)) , cudaMemcpyDeviceToHost ); _CHECK(__LINE__)
    enddo
    ierr = cudaThreadSynchronize()
    call t_stopf('CUDA QDP D2H')
    !$OMP END MASTER
    !$OMP BARRIER
  end subroutine copy_qdp_d2h



  subroutine precompute_divdp_cuda( elem , hybrid , deriv , dt , nets , nete , n0_qdp )
    use element_mod, only: element_t
    use derivative_mod, only: derivative_t, divergence_sphere
    use hybrid_mod, only: hybrid_t
    use edge_mod, only: edgeVpack,edgeVunpack
    use bndry_mod, only: bndry_exchangeV
    use perf_mod, only: t_startf, t_stopf
    implicit none
    type(element_t)      , intent(inout) :: elem(:)
    type (hybrid_t)      , intent(in   ) :: hybrid
    type (derivative_t)  , intent(in   ) :: deriv
    real(kind=real_kind) , intent(in   ) :: dt
    integer              , intent(in   ) :: nets , nete , n0_qdp
    integer :: ie , k
    integer :: ierr
    type(dim3) :: blockdim , griddim

    call t_startf('precompute_divdp_cuda')

    do ie = nets , nete
      qdp1_h(:,:,:,ie) = elem(ie)%state%Qdp(:,:,:,1,n0_qdp)
      dpdiss_ave_h(:,:,:,ie) = elem(ie)%derived%dpdiss_ave(:,:,:)
    enddo
    !$OMP BARRIER
    !$OMP MASTER
    ierr = cudaMemcpyAsync( dpdiss_ave_d , dpdiss_ave_h , size( dpdiss_ave_h ) , cudaMemcpyHostToDevice , streams(1) ); _CHECK(__LINE__)
    ierr = cudaMemcpyAsync( qdp1_d , qdp1_h , size( qdp1_h ) , cudaMemcpyHostToDevice , streams(1) ); _CHECK(__LINE__)
    blockdim = dim3( np     , np , nlev )
    griddim  = dim3( nelemd , 1  , 1    )
    call unpack_qdp1<<<griddim,blockdim,0,streams(1)>>>( qdp1_d , qdp_d , n0_qdp , nets , nete ); _CHECK(__LINE__)
    !$OMP END MASTER

    do ie = nets , nete
      do k = 1 , nlev   ! div( U dp Q), 
        elem(ie)%derived%divdp(:,:,k) = divergence_sphere(elem(ie)%derived%vn0(:,:,:,k),deriv,elem(ie))
        elem(ie)%derived%eta_dot_dpdn(:,:,k) = elem(ie)%spheremp(:,:) * elem(ie)%derived%eta_dot_dpdn(:,:,k)
        elem(ie)%derived%omega_p     (:,:,k) = elem(ie)%spheremp(:,:) * elem(ie)%derived%omega_p     (:,:,k)
        elem(ie)%derived%divdp_proj  (:,:,k) = elem(ie)%spheremp(:,:) * elem(ie)%derived%divdp       (:,:,k)
      enddo
      call edgeVpack(edgeAdvDSS,elem(ie)%derived%eta_dot_dpdn(:,:,1:nlev),nlev,0*nlev,elem(ie)%desc)
      call edgeVpack(edgeAdvDSS,elem(ie)%derived%omega_p     (:,:,1:nlev),nlev,1*nlev,elem(ie)%desc)
      call edgeVpack(edgeAdvDSS,elem(ie)%derived%divdp_proj  (:,:,1:nlev),nlev,2*nlev,elem(ie)%desc)
    enddo
    call t_startf('eus_cuda_bexchV')
    call bndry_exchangeV(hybrid,edgeAdvDSS)
    call t_stopf('eus_cuda_bexchV')
    do ie = nets , nete
      call edgeVunpack(edgeAdvDSS,elem(ie)%derived%eta_dot_dpdn(:,:,1:nlev),nlev,0*nlev,elem(ie)%desc)
      call edgeVunpack(edgeAdvDSS,elem(ie)%derived%omega_p     (:,:,1:nlev),nlev,1*nlev,elem(ie)%desc)
      call edgeVunpack(edgeAdvDSS,elem(ie)%derived%divdp_proj  (:,:,1:nlev),nlev,2*nlev,elem(ie)%desc)
      do k = 1 , nlev
        elem(ie)%derived%eta_dot_dpdn(:,:,k) = elem(ie)%derived%eta_dot_dpdn(:,:,k) * elem(ie)%rspheremp(:,:)
        elem(ie)%derived%omega_p     (:,:,k) = elem(ie)%derived%omega_p     (:,:,k) * elem(ie)%rspheremp(:,:)
        elem(ie)%derived%divdp_proj  (:,:,k) = elem(ie)%derived%divdp_proj  (:,:,k) * elem(ie)%rspheremp(:,:)
      enddo
    enddo

    do ie = nets , nete
      divdp_h(:,:,:,ie) = elem(ie)%derived%divdp(:,:,:)
      dpdiss_biharmonic_h(:,:,:,ie) = elem(ie)%derived%dpdiss_biharmonic(:,:,:)
    enddo
    ierr = cudaMemcpyAsync( divdp_d , divdp_h , size( divdp_h ) , cudaMemcpyHostToDevice , streams(1) ); _CHECK(__LINE__)
    ierr = cudaMemcpyAsync( dpdiss_biharmonic_d , dpdiss_biharmonic_h , size( dpdiss_biharmonic_h ) , cudaMemcpyHostToDevice , streams(1) ); _CHECK(__LINE__)

    call t_stopf('precompute_divdp_cuda')


  end subroutine precompute_divdp_cuda



  subroutine euler_step_cuda( np1_qdp , n0_qdp , dt , elem , hvcoord , hybrid , deriv , nets , nete , DSSopt , rhs_multiplier )
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
  use kinds             , only: real_kind
  use dimensions_mod    , only: np, npdg, nlev, qsize
  use hybrid_mod        , only: hybrid_t
  use element_mod       , only: element_t
  use derivative_mod    , only: derivative_t, divergence_sphere, gradient_sphere, vorticity_sphere
  use edge_mod          , only: edgevpack, edgevunpack
  use bndry_mod         , only: bndry_exchangev
  use hybvcoord_mod     , only: hvcoord_t
  use control_mod       , only: nu_q, nu_p, limiter_option
  use perf_mod          , only: t_startf, t_stopf  ! _EXTERNAL
  use viscosity_mod     , only: biharmonic_wk_scalar, biharmonic_wk_scalar_minmax, neighbor_minmax
  use parallel_mod      , only: iam
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
  real(kind=real_kind) :: dp0,tmp,tmp_min,tmp_max,qtens_tmp(np,np)
  real(kind=real_kind) :: large,small
  integer :: ie,q,i,j,k
  integer :: rhs_viss = 0

  integer :: ierr
  type(dim3) :: blockdim , griddim

  large =  huge(large)
  small = (huge(small)-1)*-1

  !call t_startf('euler_step')
  call t_startf('euler_step_cuda')

  do ie = nets , nete
    ! Compute velocity used to advance Qdp 
    ! derived variable divdp_proj() (DSS'd version of divdp) will only be correct on 2nd and 3rd stage
    ! but that's ok because rhs_multiplier=0 on the first stage:
    dp_h(:,:,:,ie) = elem(ie)%derived%dp(:,:,:) - rhs_multiplier * dt * elem(ie)%derived%divdp_proj(:,:,:) 
    vstar_h(:,:,:,1,ie) = elem(ie)%derived%vn0(:,:,1,:) / dp_h(:,:,:,ie)
    vstar_h(:,:,:,2,ie) = elem(ie)%derived%vn0(:,:,2,:) / dp_h(:,:,:,ie)
  enddo
  !$OMP BARRIER
  if (hybrid%ithr == 0) then
    ierr = cudaMemcpyAsync( dp_d , dp_h , size( dp_h ) , cudaMemcpyHostToDevice , streams(1) ); _CHECK(__LINE__)
    ierr = cudaStreamSynchronize(streams(1))
  endif

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
    !$OMP BARRIER
    if (hybrid%ithr == 0) then
      !$acc parallel loop gang vector collapse(5) deviceptr(qdp_d,dp_d,qtens_biharmonic_d) async(1)
      do ie = 1 , nelemd
        do q = 1 , qsize
          do k = 1 , nlev    
            do j = 1 , np
              do i = 1 , np
                qtens_biharmonic_d(i,j,k,q,ie) = qdp_d(i,j,k,q,n0_qdp,ie) / dp_d(i,j,k,ie)
              enddo
            enddo
          enddo
        enddo
      enddo
      !$acc parallel loop gang vector collapse(3) deviceptr(qmin_d,qmax_d,qtens_biharmonic_d) private(tmp_min,tmp_max) async(1)
      do ie = 1 , nelemd
        do q = 1 , qsize
          do k = 1 , nlev    
            if ( rhs_multiplier == 0 .or. rhs_multiplier == 2 ) then
              tmp_min = large
              tmp_max = small
            else
              tmp_min = qmin_d(k,q,ie)
              tmp_max = qmax_d(k,q,ie)
            endif
            do j = 1 , np
              do i = 1 , np
                tmp_min = min(tmp_min,qtens_biharmonic_d(i,j,k,q,ie))
                tmp_max = max(tmp_max,qtens_biharmonic_d(i,j,k,q,ie))
              enddo
            enddo
            qmin_d(k,q,ie) = max(tmp_min,0d0)
            qmax_d(k,q,ie) = tmp_max
          enddo
        enddo
      enddo
      if ( nu_p > 0 .and. rhs_multiplier == 2 ) then
        !$acc parallel loop gang vector collapse(5) deviceptr(qtens_biharmonic_d,dpdiss_ave_d,dp0_d) private(tmp_min,tmp_max) async(1)
        do ie = 1 , nelemd
          do q = 1 , qsize
            do k = 1 , nlev
              do j = 1 , np
                do i = 1 , np
                  ! two scalings depending on nu_p:
                  ! nu_p=0:    qtens_biharmonic *= dp0                   (apply viscsoity only to q)
                  ! nu_p>0):   qtens_biharmonic *= elem()%psdiss_ave      (for consistency, if nu_p=nu_q)
                  ! NOTE: divide by dp0 since we multiply by dp0 below
                  qtens_biharmonic_d(i,j,k,q,ie) = qtens_biharmonic_d(i,j,k,q,ie)*dpdiss_ave_d(i,j,k,ie)/dp0_d(k)
                enddo
              enddo
            enddo
          enddo
        enddo
      endif
    endif
    if ( rhs_multiplier == 0 ) call neighbor_minmax_loc(elem,hybrid,nets,nete,qmin_d,qmax_d)   ! update qmin/qmax based on neighbor data for lim8
    if ( rhs_multiplier == 2 ) call biharmonic_wk_scalar_minmax_loc( elem , qtens_biharmonic_d , deriv , hybrid , nets , nete , qmin_d , qmax_d )
    if (hybrid%ithr == 0) then
      if ( rhs_multiplier == 2 ) then
        !compute biharmonic mixing term
        rhs_viss = 3
        !$acc parallel loop gang vector collapse(5) deviceptr(spheremp_d,dp0_d) async(1)
        do ie = 1 , nelemd
          do q = 1 , qsize
            do k = 1 , nlev
              do j = 1 , np
                do i = 1 , np
                  qtens_biharmonic_d(i,j,k,q,ie) = -rhs_viss*dt*nu_q*dp0_d(k)*qtens_biharmonic_d(i,j,k,q,ie) / spheremp_d(i,j,ie)
                enddo
              enddo
            enddo
          enddo
        enddo
      endif
      !$acc wait(1)
    endif
    !$OMP BARRIER
  endif  ! compute biharmonic mixing term and qmin/qmax

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   2D Advection step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !$OMP BARRIER
  !$OMP MASTER
  ierr = cudaMemcpyAsync( vstar_d , vstar_h , size( vstar_h ) , cudaMemcpyHostToDevice , streams(1) ); _CHECK(__LINE__)
  blockdim = dim3( np*np*numk_eul , 1 , 1 )
  griddim  = dim3( int(ceiling(dble(nlev)/numk_eul))*qsize_d*nelemd , 1 , 1 )
  call euler_step_kernel1<<<griddim,blockdim,0,streams(1)>>>( qdp_d , qtens_d , spheremp_d , vstar_d , metdet_d , rmetdet_d , dinv_d , deriv_dvv_d , &
                                                              n0_qdp , np1_qdp , rhs_viss , dt , 1 , nelemd ); _CHECK(__LINE__)

  ierr = cudathreadsynchronize()
  if ( limiter_option == 8 ) then
    !$acc parallel loop gang vector collapse(5) deviceptr(qtens_d,dp_star_d,dp_d,divdp_d,spheremp_d,qtens_biharmonic_d) private(tmp) async(1)
    do ie = 1 , nelemd
      do q = 1 , qsize
        do k = 1 , nlev
          do j = 1 , np
            do i = 1 , np
              if ( rhs_viss /= 0 ) Qtens_d(i,j,k,q,ie) = Qtens_d(i,j,k,q,ie) + Qtens_biharmonic_d(i,j,k,q,ie)
              tmp = dp_d(i,j,k,ie) - dt * divdp_d(i,j,k,ie)    ! UN-DSS'ed dp at timelevel n0+1:   
              if ( nu_p > 0 .and. rhs_viss /= 0 ) tmp = tmp - rhs_viss * dt * nu_q * dpdiss_biharmonic_d(i,j,k,ie) / spheremp_d(i,j,ie)
              dp_star_d(i,j,k,ie) = tmp
            enddo
          enddo
        enddo
      enddo
    enddo
    call limiter_optim_iter_full_loc( Qtens_d , spheremp_d , qmin_d , qmax_d , dp_star_d )  ! apply limiter to Q = Qtens / dp_star 
    !$acc wait(1)
  endif
  blockdim = dim3( np*np*numk_eul , 1 , 1 )
  griddim  = dim3( int(ceiling(dble(nlev)/numk_eul))*qsize_d*nelemd , 1 , 1 )
  call euler_step_kernel2<<<griddim,blockdim,0,streams(1)>>>( qdp_d , qtens_d , spheremp_d , np1_qdp , 1 , nelemd ); _CHECK(__LINE__)
  if ( limiter_option == 4 ) then
    blockdim = dim3( np*np*numk_lim2d , 1 , 1 )
    griddim  = dim3( int(ceiling(dble(nlev)/numk_lim2d))*qsize_d*nelemd , 1 , 1 )
    call limiter2d_zero_kernel<<<griddim,blockdim,0,streams(1)>>>( qdp_d , 1 , nelemd , np1_qdp ); _CHECK(__LINE__)
  endif
  call t_startf('eus_cuda_peus')
  call pack_exchange_unpack_stage(np1_qdp,hybrid,qdp_d,timelevels)
  call t_stopf('eus_cuda_peus')
  blockdim = dim3( np*np   * nlev  , 1 , 1 )
  griddim  = dim3( qsize_d , nelemd , 1 )
  call euler_hypervis_kernel_last<<<griddim,blockdim,0,streams(1)>>>( qdp_d , rspheremp_d , 1 , nelemd , np1_qdp ); _CHECK(__LINE__)
  !$OMP END MASTER
  !$OMP BARRIER

  call t_stopf('euler_step_cuda')
  !call t_stopf('euler_step')
end subroutine euler_step_cuda



subroutine qdp_time_avg_cuda( elem , rkstage , n0_qdp , np1_qdp , limiter_option , nu_p , nets , nete )
  use element_mod, only: element_t
  implicit none
  type(element_t)     , intent(inout) :: elem(:)
  integer             , intent(in   ) :: rkstage , n0_qdp , np1_qdp , nets , nete, limiter_option
  real(kind=real_kind), intent(in   ) :: nu_p
  integer :: ie
  type(dim3) :: griddim , blockdim
  integer :: ierr
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$OMP BARRIER
!$OMP MASTER
  blockdim = dim3( np      , np     , nlev )
  griddim  = dim3( qsize_d , nelemd , 1    )
  call qdp_time_avg_kernel<<<griddim,blockdim,0,streams(1)>>>( qdp_d , rkstage , n0_qdp , np1_qdp , 1 , nelemd ); _CHECK(__LINE__)
  if (limiter_option == 8) then
    blockdim = dim3( np , np , nlev )
    griddim  = dim3( nelemd , 1 , 1 )
    call pack_qdp1<<<griddim,blockdim,0,streams(1)>>>( qdp1_d , qdp_d , np1_qdp , 1 , nelemd ); _CHECK(__LINE__)
    ierr = cudaMemcpyAsync( qdp1_h , qdp1_d , size( qdp1_h ) , cudaMemcpyDeviceToHost , streams(1) ); _CHECK(__LINE__)
    ierr = cudaStreamSynchronize(streams(1))
  endif
!$OMP END MASTER
!$OMP BARRIER
  if (limiter_option == 8) then
    do ie = nets , nete
      elem(ie)%state%Qdp(:,:,:,1,np1_qdp) = qdp1_h(:,:,:,ie)
    enddo
  endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine qdp_time_avg_cuda



subroutine advance_hypervis_scalar_cuda( edgeAdv , elem , hvcoord , hybrid , deriv , nt , nt_qdp , nets , nete , dt2 )
  !  hyperviscsoity operator for foward-in-time scheme
  !  take one timestep of:  
  !          Q(:,:,:,np) = Q(:,:,:,np) +  dt2*nu*laplacian**order ( Q )
  !
  !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
  use kinds          , only : real_kind
  use dimensions_mod , only : np, nlev
  use hybrid_mod     , only : hybrid_t
  use hybvcoord_mod  , only : hvcoord_t
  use element_mod    , only : element_t
  use derivative_mod , only : derivative_t
  use edge_mod       , only : EdgeBuffer_t, edgevpack, edgevunpack
  use bndry_mod      , only : bndry_exchangev
  use perf_mod       , only : t_startf, t_stopf, t_barrierf
  use control_mod    , only : rsplit, nu_q, hypervis_order, hypervis_subcycle_q , nu_p
  use parallel_mod   , only : iam
  implicit none
  type (EdgeBuffer_t)  , intent(inout)         :: edgeAdv
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
  integer :: k , kptr , i , j , ie , ic , q , ierr
  real (kind=real_kind) :: dt
  type(dim3) :: griddim , blockdim
  if ( nu_q           == 0 ) return
  if ( hypervis_order /= 2 ) return
  call t_barrierf('sync_advance_hypervis_scalar', hybrid%par%comm)
  call t_startf('advance_hypervis_scalar_cuda')

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  hyper viscosity  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dt = dt2 / hypervis_subcycle_q
  do ic = 1 , hypervis_subcycle_q
    ! Qtens = Q/dp   (apply hyperviscsoity to dp0 * Q, not Qdp)
    do ie = nets , nete
      do k = 1 , nlev
        dp_h(:,:,k,ie) = elem(ie)%derived%dp(:,:,k) - dt2*elem(ie)%derived%divdp_proj(:,:,k)
      enddo
    enddo
!$OMP BARRIER
!$OMP MASTER
    ierr = cudaMemcpyAsync( dp_d         , dp_h         , size( dp_h )         , cudaMemcpyHostToDevice , streams(1) ); _CHECK(__LINE__)

    !KERNEL 1
    blockdim = dim3( np*np*numk_hyp , 1 , 1 )
    griddim  = dim3( int(ceiling(dble(nlev)/numk_hyp))*qsize_d*nelemd , 1 , 1 )
    call hypervis_kernel1<<<griddim,blockdim,0,streams(1)>>>( qdp_d , qtens_d , dp_d , dinv_d , variable_hyperviscosity_d , dpdiss_ave_d , spheremp_d , &
                                                              deriv_dvv_d , hyai_d , hybi_d , hvcoord%ps0 , 1 , nelemd , dt , nt_qdp , nu_p ); _CHECK(__LINE__)
    call t_startf('ahs_cuda_peus1')
    call pack_exchange_unpack_stage(1,hybrid,qtens_d,1)
    call t_stopf('ahs_cuda_peus1')

    !KERNEL 2
    blockdim = dim3( np*np*numk_hyp , 1 , 1 )
    griddim  = dim3( int(ceiling(dble(nlev)/numk_hyp))*qsize_d*nelemd , 1 , 1 )
    call hypervis_kernel2<<<griddim,blockdim,0,streams(1)>>>( qdp_d , qtens_d , dp_d , dinv_d , variable_hyperviscosity_d , spheremp_d , rspheremp_d , deriv_dvv_d , &
                                                              nu_q , 1 , nelemd , dt , nt_qdp ); _CHECK(__LINE__)
    blockdim = dim3( np*np*numk_lim2d , 1 , 1 )
    griddim  = dim3( int(ceiling(dble(nlev)/numk_lim2d))*qsize_d*nelemd , 1 , 1 )
    call limiter2d_zero_kernel<<<griddim,blockdim,0,streams(1)>>>(qdp_d,nets,nete,nt_qdp); _CHECK(__LINE__)

    call t_startf('ahs_cuda_peus2')
    call pack_exchange_unpack_stage(nt_qdp,hybrid,qdp_d,timelevels)
    call t_stopf('ahs_cuda_peus2')

    !KERNEL 3
    blockdim = dim3( np * np * nlev , 1 , 1 )
    griddim  = dim3( qsize_d , nelemd , 1 )
    call euler_hypervis_kernel_last<<<griddim,blockdim,0,streams(1)>>>( qdp_d , rspheremp_d , nets , nete , nt_qdp ); _CHECK(__LINE__)
    blockdim = dim3( np , np , nlev )
    griddim  = dim3( nelemd , 1 , 1 )
    call pack_qdp1<<<griddim,blockdim,0,streams(1)>>>( qdp1_d , qdp_d , nt_qdp , 1 , nelemd ); _CHECK(__LINE__)
    ierr = cudaMemcpyAsync( qdp1_h , qdp1_d , size( qdp1_h ) , cudaMemcpyDeviceToHost , streams(1) ); _CHECK(__LINE__)
    ierr = cudaStreamSynchronize(streams(1))
!$OMP END MASTER
!$OMP BARRIER
    do ie = nets , nete
      elem(ie)%state%Qdp(:,:,:,1,nt_qdp) = qdp1_h(:,:,:,ie)
    enddo
  enddo

  call t_stopf('advance_hypervis_scalar_cuda')

end subroutine advance_hypervis_scalar_cuda



attributes(global) subroutine pack_qdp1( qdp1 , qdp , nt , nets , nete )
  implicit none
  real(kind=real_kind), intent(  out) :: qdp1( np , np , nlev                        , nets:nete )
  real(kind=real_kind), intent(in   ) :: qdp ( np , np , nlev , qsize_d , timelevels , nets:nete )
  integer, value      , intent(in   ) :: nt , nets , nete
  integer :: i,j,k,ie
  i  = threadidx%x
  j  = threadidx%y
  k  = threadidx%z
  ie = blockidx%x
  qdp1(i,j,k,ie) = qdp(i,j,k,1,nt,ie)
end subroutine pack_qdp1



attributes(global) subroutine unpack_qdp1( qdp1 , qdp , nt , nets , nete )
  implicit none
  real(kind=real_kind), intent(in   ) :: qdp1( np , np , nlev                        , nets:nete )
  real(kind=real_kind), intent(  out) :: qdp ( np , np , nlev , qsize_d , timelevels , nets:nete )
  integer, value      , intent(in   ) :: nt , nets , nete
  integer :: i,j,k,ie
  i  = threadidx%x
  j  = threadidx%y
  k  = threadidx%z
  ie = blockidx%x
  qdp(i,j,k,1,nt,ie) = qdp1(i,j,k,ie)
end subroutine unpack_qdp1



attributes(global) subroutine qdp_time_avg_kernel( qdp , rkstage , n0_qdp , np1_qdp , nets , nete )
  implicit none
  real(kind=real_kind), intent(inout) :: qdp( np , np , nlev , qsize_d , timelevels , nets:nete )
  integer, value      , intent(in   ) :: rkstage , n0_qdp , np1_qdp , nets , nete
  integer :: i, j, k, q, ie
  i  = threadidx%x
  j  = threadidx%y
  k  = threadidx%z
  q  = blockidx%x
  ie = blockidx%y
  qdp(i,j,k,q,np1_qdp,ie) = ( qdp(i,j,k,q,n0_qdp,ie) + (rkstage-1) * qdp(i,j,k,q,np1_qdp,ie) ) / rkstage
end subroutine qdp_time_avg_kernel



attributes(global) subroutine euler_step_kernel1( Qdp , qtens , spheremp , vstar , metdet , rmetdet , dinv , deriv_dvv , n0_qdp , np1_qdp , rhs_viss , dt , nets , nete )
  implicit none
  real(kind=real_kind), dimension(np,np,nlev,qsize_d,timelevels,nets:nete), intent(in   ) :: Qdp
  real(kind=real_kind), dimension(np,np,nlev,qsize_d           ,nets:nete), intent(  out) :: qtens
  real(kind=real_kind), dimension(np,np                        ,nets:nete), intent(in   ) :: spheremp
  real(kind=real_kind), dimension(np,np,nlev,2                 ,nets:nete), intent(in   ) :: vstar
  real(kind=real_kind), dimension(np,np                        ,nets:nete), intent(in   ) :: metdet
  real(kind=real_kind), dimension(np,np                        ,nets:nete), intent(in   ) :: rmetdet
  real(kind=real_kind), dimension(np,np,2,2                    ,nets:nete), intent(in   ) :: dinv
  real(kind=real_kind), dimension(np,np                                  ), intent(in   ) :: deriv_dvv
  real(kind=real_kind), value                                             , intent(in   ) :: dt
  integer, value                                                          , intent(in   ) :: n0_qdp, np1_qdp, rhs_viss, nets, nete
  integer :: ks
  integer :: i , j , k , kk , q , ie , ij
  real(kind=real_kind), shared :: deriv_dvv_s(np*np+1)
  real(kind=real_kind) :: qtmp

  ks = int(ceiling(dble(nlev)/numk_eul))

  !Define the indices
  i  = modulo( threadidx%x-1    ,np)+1
  j  = modulo((threadidx%x-1)/np,np)+1
  kk = (threadidx%x-1)/(np*np)+1
  k  = modulo(blockidx%x-1,ks)*numk_eul + kk
  q  = modulo((blockidx%x-1)/ks,qsize_d)+1
  ie =       ((blockidx%x-1)/ks)/qsize_d+1
  ij   =              (j-1)*np+i
  if (k  <= nlev .and. q  <= qsize_d .and. ie <= nete ) then
    if (kk == 1) deriv_dvv_s(ij) = deriv_dvv(i,j)
    qtmp = Qdp(i,j,k,q,n0_qdp,ie)
    qtens(i,j,k,q,ie) = qtmp - dt * divergence_sphere( i , j , k , q , ie , kk , ij , vstar , qtmp , metdet , rmetdet , dinv , deriv_dvv_s , nets , nete )
  endif
end subroutine euler_step_kernel1



attributes(global) subroutine euler_step_kernel2( Qdp , qtens , spheremp , np1_qdp , nets , nete )
  implicit none
  real(kind=real_kind), dimension(np,np,nlev,qsize_d,timelevels,nets:nete), intent(  out) :: Qdp
  real(kind=real_kind), dimension(np,np,nlev,qsize_d           ,nets:nete), intent(in   ) :: qtens
  real(kind=real_kind), dimension(np,np                        ,nets:nete), intent(in   ) :: spheremp
  integer, value                                                          , intent(in   ) :: np1_qdp, nets, nete
  integer :: i,j,k,q,ie,kk,ks
  ks = int(ceiling(dble(nlev)/numk_eul))
  !Define the indices
  i  = modulo( threadidx%x-1    ,np)+1
  j  = modulo((threadidx%x-1)/np,np)+1
  kk = (threadidx%x-1)/(np*np)+1
  k  = modulo(blockidx%x-1,ks)*numk_eul + kk
  q  = modulo((blockidx%x-1)/ks,qsize_d)+1
  ie =       ((blockidx%x-1)/ks)/qsize_d+1
  if (k  > nlev   ) return
  if (q  > qsize_d) return
  if (ie > nete   ) return
  Qdp(i,j,k,q,np1_qdp,ie) = spheremp(i,j,ie) * Qtens(i,j,k,q,ie)
end subroutine euler_step_kernel2



attributes(device) function divergence_sphere(i,j,k,q,ie,kk,ij,Vstar,qtmp,metdet,rmetdet,dinv,deriv_dvv,nets,nete) result(dp_star)
  implicit none
  integer,              intent(in) :: i, j, k, q, ie, kk, ij, nets, nete
  real(kind=real_kind), intent(in) :: Dinv     (np,np,2,2,nets:nete)
  real(kind=real_kind), intent(in) :: metdet   (np,np,nets:nete)
  real(kind=real_kind), intent(in) :: rmetdet  (np,np,nets:nete)
  real(kind=real_kind), intent(in) :: deriv_dvv(np*np+1)
  real(kind=real_kind), intent(in) :: Vstar    (np,np,nlev,2,nets:nete)
  real(kind=real_kind), intent(in), value :: qtmp
  real(kind=real_kind)             :: dp_star
  real(kind=real_kind), shared :: gv(np*np+1,numk_eul,2)
  integer :: s
  real(kind=real_kind) :: vvtemp, divtemp, vs1tmp, vs2tmp
  vs1tmp = vstar(i,j,k,1,ie) * metdet(i,j,ie) * qtmp
  vs2tmp = vstar(i,j,k,2,ie) * metdet(i,j,ie) * qtmp
  gv(ij,kk,1) = dinv(i,j,1,1,ie) * vs1tmp + dinv(i,j,1,2,ie) * vs2tmp
  gv(ij,kk,2) = dinv(i,j,2,1,ie) * vs1tmp + dinv(i,j,2,2,ie) * vs2tmp
  call syncthreads()
  divtemp = 0.0d0
  vvtemp   = 0.0d0
  do s = 1 , np
    divtemp = divtemp + deriv_dvv((i-1)*np+s) * gv((j-1)*np+s,kk,1)
    vvtemp  = vvtemp  + deriv_dvv((j-1)*np+s) * gv((s-1)*np+i,kk,2)
  end do
  dp_star = ( divtemp + vvtemp ) * ( rmetdet(i,j,ie) * rrearth_d )
end function divergence_sphere



attributes(global) subroutine limiter2d_zero_kernel(Qdp,nets,nete,np1)
  use kinds, only : real_kind
  use dimensions_mod, only : np, nlev
  implicit none
  real (kind=real_kind), intent(inout) :: Qdp(np,np,nlev,qsize_d,timelevels,nets:nete)
  integer, value       , intent(in   ) :: nets,nete,np1
  integer :: i, j, k, kk, q, ie, jj, ij, ijk, ks
  real (kind=real_kind) :: mass,mass_new
  real (kind=real_kind), shared :: Qdp_shared(np*np+PAD,nlev)
  real (kind=real_kind), shared :: mass_shared(nlev)
  real (kind=real_kind), shared :: mass_new_shared(nlev)

  ks = int(ceiling(dble(nlev)/numk_lim2d))

  i  = modulo( threadidx%x-1    ,np)+1
  j  = modulo((threadidx%x-1)/np,np)+1
  kk = (threadidx%x-1)/(np*np)+1
  k  = modulo(blockidx%x-1,ks)*numk_lim2d + kk
  q  = modulo((blockidx%x-1)/ks,qsize_d)+1
  ie =       ((blockidx%x-1)/ks)/qsize_d+1
  ij  = (j-1)*np+i
  ijk = (kk-1)*np*np+ij

  Qdp_shared(ij,kk) = Qdp(i,j,k,q,np1,ie)
  call syncthreads()

  if ( ijk <= numk_lim2d ) then
    mass = 0.
    do jj = 1 , np*np
      mass = mass + Qdp_shared(jj,ijk)
    enddo
    mass_shared(ijk) = mass
  endif
  call syncthreads()

  if ( mass_shared(kk)   < 0 ) Qdp_shared(ij,kk) = -Qdp_shared(ij,kk)
  if ( Qdp_shared(ij,kk) < 0 ) Qdp_shared(ij,kk) = 0
  call syncthreads()

  if ( ijk <= numk_lim2d ) then
    mass = 0.
    do jj = 1 , np*np
      mass = mass + Qdp_shared(jj,ijk)
    enddo
    mass_new_shared(ijk) = mass
  endif
  call syncthreads()

  ! now scale the all positive values to restore mass
  if ( mass_new_shared(kk) > 0 ) Qdp_shared(ij,kk) =  Qdp_shared(ij,kk) * abs(mass_shared(kk)) / mass_new_shared(kk)
  if ( mass_shared    (kk) < 0 ) Qdp_shared(ij,kk) = -Qdp_shared(ij,kk)
  Qdp(i,j,k,q,np1,ie) = Qdp_shared(ij,kk)
end subroutine limiter2d_zero_kernel



attributes(global) subroutine euler_hypervis_kernel_last( Qdp , rspheremp , nets , nete , np1 )
  implicit none
  real(kind=real_kind), dimension(np*np*nlev,qsize_d,timelevels,nets:nete), intent(inout) :: Qdp
  real(kind=real_kind), dimension(np*np                        ,nets:nete), intent(in   ) :: rspheremp
  integer, value                                                          , intent(in   ) :: nets , nete , np1
  integer :: ijk, ij, q, ie
  ijk  = threadidx%x
  ij   = modulo(threadidx%x-1,np*np)+1
  q  = blockidx%x
  ie = blockidx%y
  Qdp(ijk,q,np1,ie) = rspheremp(ij,ie) * Qdp(ijk,q,np1,ie)
end subroutine euler_hypervis_kernel_last



subroutine pack_exchange_unpack_stage(np1,hybrid,array_in,tl_in)
  use hybrid_mod, only : hybrid_t
  use schedule_mod, only : schedule_t, schedule, cycle_t
  use parallel_mod, only : abortmp, status, srequest, rrequest, mpireal_t, mpiinteger_t, iam
  use perf_mod, only: t_startf, t_stopf
  implicit none
  include 'mpif.h'
  type(hybrid_t)              , intent(in   ) :: hybrid
  real(kind=real_kind), device, intent(inout) :: array_in(np,np,nlev,qsize_d,tl_in,nelemd)
  integer, value              , intent(in   ) :: np1 , tl_in
  ! local
  type(dim3)                :: griddim6,blockdim6
  integer                   :: icycle,nSendCycles,nRecvCycles,n, ierr
  type (Schedule_t),pointer :: pSchedule
  type (Cycle_t),pointer    :: pCycle
  integer                   :: dest,length,tag,iptr,source,nlyr,query_sum, npacked
  logical :: recvflag, internal_unpacked
  real :: time_milli
#ifdef _PREDICT
  pSchedule => Schedule(iam)
#else
  pSchedule => Schedule(1)
#endif
  nlyr = edgeAdv%nlyr
  nSendCycles = pSchedule%nSendCycles
  nRecvCycles = pSchedule%nRecvCycles

  ierr = cudaStreamSynchronize(streams(1))
  call t_startf('CUDA PEU NEWER')
  do icycle=1,nRecvCycles
    pCycle => pSchedule%RecvCycle(icycle)
    source  = pCycle%source - 1
    length  = nlyr * pCycle%lengthP
    tag     = pCycle%tag
    call MPI_Irecv(recvbuf_h(1,1,icycle),length,MPIreal_t,source,tag,hybrid%par%comm,Rrequest(icycle),ierr)
  enddo
  if ( recv_external_nelem > 0 ) then
    blockdim6 = dim3( max(np,max_corner_elem_d) , nlev , 1 )
    griddim6  = dim3( qsize_d , recv_external_nelem , 1    )
    call edgeVpack_kernel_stage<<<griddim6,blockdim6,0,streams(1)>>>(edgebuf_d,array_in,putmapP_d,reverse_d,nbuf,0,1,nelemd,np1,recv_external_indices_d,tl_in); _CHECK(__LINE__)
  endif
  do icycle = 1 , nSendCycles
    pCycle => pSchedule%SendCycle(icycle)
    iptr   =  pCycle%ptrP
    ierr = cudaMemcpyAsync(sendbuf_h(1,1,icycle),edgebuf_d(1,iptr),size(sendbuf_h(1:nlyr,1:pCycle%lengthP,icycle)),cudaMemcpyDeviceToHost,streams(1)); _CHECK(__LINE__)
  enddo
  if ( recv_internal_nelem > 0 ) then
    blockdim6 = dim3( max(np,max_corner_elem_d) , nlev , 1 )
    griddim6  = dim3( qsize_d , recv_internal_nelem , 1    )
    call edgeVpack_kernel_stage<<<griddim6,blockdim6,0,streams(2)>>>(edgebuf_d,array_in,putmapP_d,reverse_d,nbuf,0,1,nelemd,np1,recv_internal_indices_d,tl_in); _CHECK(__LINE__)
    blockdim6 = dim3( max(np,max_corner_elem_d) , nlev , 1 )
    griddim6  = dim3( qsize_d , recv_internal_nelem , 1    )
    call edgeVunpack_kernel_stage<<<griddim6,blockdim6,0,streams(2)>>>(edgebuf_d,array_in,getmapP_d,nbuf,0,1,nelemd,np1,recv_internal_indices_d,tl_in); _CHECK(__LINE__)
  endif
  do icycle = 1 , nSendCycles
    pCycle => pSchedule%SendCycle(icycle)
    iptr   =  pCycle%ptrP
    ierr = cudaMemcpyAsync(sendbuf_h(1,1,icycle),edgebuf_d(1,iptr),size(sendbuf_h(1:nlyr,1:pCycle%lengthP,icycle)),cudaMemcpyDeviceToHost,streams(1)); _CHECK(__LINE__)
  enddo
  ierr = cudaStreamSynchronize(streams(1))
  do icycle = 1 , nSendCycles
    pCycle => pSchedule%SendCycle(icycle)
    dest   =  pCycle%dest - 1
    iptr   =  pCycle%ptrP
    length =  nlyr * pCycle%lengthP
    tag    =  pCycle%tag
    call MPI_Isend(sendbuf_h(1,1,icycle),length,MPIreal_t,dest,tag,hybrid%par%comm,Srequest(icycle),ierr)
  enddo
  call MPI_WaitAll(nRecvCycles,Rrequest,status,ierr)
  !When this cycle's MPI transfer is compliete, then call the D2H memcopy asynchronously
  do icycle = 1 , nRecvCycles
    pCycle => pSchedule%RecvCycle(icycle)
    iptr   =  pCycle%ptrP
    ierr = cudaMemcpyAsync(edgebuf_d(1,iptr),recvbuf_h(1,1,icycle),size(recvbuf_h(1:nlyr,1:pCycle%lengthP,icycle)),cudaMemcpyHostToDevice,streams(1)); _CHECK(__LINE__)
  enddo
  call MPI_WaitAll(nSendCycles,Srequest,status,ierr)
  if ( recv_external_nelem > 0 ) then
    blockdim6 = dim3( max(np,max_corner_elem_d) , nlev , 1 )
    griddim6  = dim3( qsize_d , recv_external_nelem , 1    )
    call edgeVunpack_kernel_stage<<<griddim6,blockdim6,0,streams(1)>>>(edgebuf_d,array_in,getmapP_d,nbuf,0,1,nelemd,np1,recv_external_indices_d,tl_in); _CHECK(__LINE__)
  endif
  ierr = cudaStreamSynchronize(streams(2))
  call t_stopf('CUDA PEU NEWER')

end subroutine pack_exchange_unpack_stage



attributes(global) subroutine edgeVpack_kernel_stage(edgebuf,v,putmapP,reverse,nbuf,kptr,nets,nete,nt,indices,tl_in)
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest
  implicit none
  real(kind=real_kind)   ,intent(in   ) :: v(np*np*nlev,qsize_d,tl_in,nets:nete)
  real(kind=real_kind)   ,intent(  out) :: edgebuf(nlev*qsize_d,nbuf)
  integer, value         ,intent(in   ) :: nbuf
  integer, value         ,intent(in   ) :: kptr
  integer, value         ,intent(in   ) :: nt
  integer                ,intent(in   ) :: putmapP(max_neigh_edges_d,nets:nete)
  logical                ,intent(in   ) :: reverse(max_neigh_edges_d,nets:nete)
  integer                ,intent(in   ) :: indices(nets:nete)
  integer, value         ,intent(in   ) :: tl_in
  integer, value         ,intent(in   ) :: nets
  integer, value         ,intent(in   ) :: nete
  integer :: i,k,ir,ll,kq,ie,q,ii,ik,koff,ind,ind_k,ind_ij
  real(kind=real_kind), shared :: v_s((np*np+1)*nlev)
  integer             , shared :: put_s(20)
  logical             , shared :: rev_s(4)

  i  = threadidx%x
  k  = threadidx%y
  q  = blockidx%x
  ie = indices(blockidx%y)
  kq = (q-1)*nlev+k + kptr
  ik = (k-1)*np+i
  koff = (k-1)*(np*np+1)

  !Efficiently load v, reverse, and putmapP into shared memory to re-use and to efficiently retreive DRAM memory.
  if (i <= np) then
    do ii = 1 , np  !Though this looping structure is strange, it allows the kernel to access contiguous np*nlev chunks of v.
      ind = (ii-1)*np*nlev+ik
      ind_k = (ind-1)/(np*np) + 1
      ind_ij = modulo(ind-1,np*np) + 1
      v_s((ind_k-1)*(np*np+1) + ind_ij) = v(ind,q,nt,ie)
    enddo
  endif
  if (ik <= 4                ) then
    rev_s(ik) = reverse(ik,ie)
  elseif (ik <= max_neigh_edges_d+4) then
    put_s(ik-4) = putmapP(ik-4,ie)
  endif
  call syncthreads()

  !Begin the packing process
  if (i <= np) then
    edgebuf(kq,put_s(west )+i) = v_s(koff+(i -1)*np+1 )
    edgebuf(kq,put_s(east )+i) = v_s(koff+(i -1)*np+np)
    edgebuf(kq,put_s(south)+i) = v_s(koff          +i )
    edgebuf(kq,put_s(north)+i) = v_s(koff+(np-1)*np+i )
  endif
  call syncthreads()
  if (i <= np) then
    ir = np-i+1
    if (rev_s(south)) edgebuf(kq,put_s(south)+ir) = v_s(koff          +i )
    if (rev_s(east )) edgebuf(kq,put_s(east )+ir) = v_s(koff+(i -1)*np+np)
    if (rev_s(north)) edgebuf(kq,put_s(north)+ir) = v_s(koff+(np-1)*np+i )
    if (rev_s(west )) edgebuf(kq,put_s(west )+ir) = v_s(koff+(i -1)*np+1 )
  endif
  if (i <= max_corner_elem_d) then
    ll = put_s(swest+0*max_corner_elem_d+i-1); if (ll /= -1) edgebuf(kq,ll+1) = v_s(koff          +1 )
    ll = put_s(swest+1*max_corner_elem_d+i-1); if (ll /= -1) edgebuf(kq,ll+1) = v_s(koff          +np)
    ll = put_s(swest+2*max_corner_elem_d+i-1); if (ll /= -1) edgebuf(kq,ll+1) = v_s(koff+(np-1)*np+1 )
    ll = put_s(swest+3*max_corner_elem_d+i-1); if (ll /= -1) edgebuf(kq,ll+1) = v_s(koff+(np-1)*np+np)
  endif
end subroutine edgeVpack_kernel_stage



attributes(global) subroutine edgeVunpack_kernel_stage(edgebuf,v,getmapP,nbuf,kptr,nets,nete,nt,recv_indices,tl_in)
  use edge_mod, only: EdgeDescriptor_t, EdgeBuffer_t
  use dimensions_mod, only : np, max_corner_elem
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest
  implicit none
  real(kind=real_kind)   ,intent(in   ) :: edgebuf(nlev*qsize_d,nbuf)
  real(kind=real_kind)   ,intent(inout) :: v(np*np*nlev,qsize_d,tl_in,nets:nete)
  integer                ,intent(in   ) :: getmapP(max_neigh_edges_d,nets:nete)
  integer, value         ,intent(in   ) :: nbuf
  integer, value         ,intent(in   ) :: kptr
  integer, value         ,intent(in   ) :: nets
  integer, value         ,intent(in   ) :: nete
  integer, value         ,intent(in   ) :: nt
  integer                ,intent(in   ) :: recv_indices(nets:nete)
  integer, value         ,intent(in   ) :: tl_in
  real(kind=real_kind), shared :: v_s((np*np+1)*nlev)
  integer             , shared :: get_s(20)
  integer :: i,k,ll,q,ie,kq,ik,koff,ind,ind_k,ind_ij,ii
  i  = threadidx%x
  k  = threadidx%y
  q  = blockidx%x
  ie = recv_indices(blockidx%y)
  kq = (q-1)*nlev+k + kptr
  ik = (k-1)*np+i
  koff = (k-1)*(np*np+1)

  !Efficiently load v, and getmapP into shared memory to re-use and to efficiently retreive DRAM memory.
  if (i <= np) then
    do ii = 1 , np  !Though this looping structure is strange, it allows the kernel to access contiguous np*nlev chunks of v.
      ind = (ii-1)*np*nlev+ik
      ind_k = (ind-1)/(np*np) + 1
      ind_ij = modulo(ind-1,np*np) + 1
      v_s((ind_k-1)*(np*np+1) + ind_ij) = v(ind,q,nt,ie)
    enddo
  endif
  if (ik <= max_neigh_edges_d) get_s(ik) = getmapP(ik,ie)
  call syncthreads()

  if (i <= np) then
    v_s(koff+(1 -1)*np+i ) = v_s(koff+(1 -1)*np+i ) + edgebuf(kptr+kq,get_s(south)+i)
    v_s(koff+(np-1)*np+i ) = v_s(koff+(np-1)*np+i ) + edgebuf(kptr+kq,get_s(north)+i)
  endif
  call syncthreads()
  if (i <= np) then
    v_s(koff+(i -1)*np+1 ) = v_s(koff+(i -1)*np+1 ) + edgebuf(kptr+kq,get_s(west )+i)
    v_s(koff+(i -1)*np+np) = v_s(koff+(i -1)*np+np) + edgebuf(kptr+kq,get_s(east )+i)
  endif
  call syncthreads()
  if (i <= max_corner_elem_d) then
    ll = get_s(swest+0*max_corner_elem_d+i-1); if(ll /= -1) v_s(koff+(1 -1)*np+1 ) = v_s(koff+(1 -1)*np+1 ) + edgebuf(kptr+kq,ll+1)
    ll = get_s(swest+1*max_corner_elem_d+i-1); if(ll /= -1) v_s(koff+(1 -1)*np+np) = v_s(koff+(1 -1)*np+np) + edgebuf(kptr+kq,ll+1)
    ll = get_s(swest+2*max_corner_elem_d+i-1); if(ll /= -1) v_s(koff+(np-1)*np+1 ) = v_s(koff+(np-1)*np+1 ) + edgebuf(kptr+kq,ll+1)
    ll = get_s(swest+3*max_corner_elem_d+i-1); if(ll /= -1) v_s(koff+(np-1)*np+np) = v_s(koff+(np-1)*np+np) + edgebuf(kptr+kq,ll+1)
  endif
  call syncthreads()

  !Efficiently store v_s into DRAM memory.
  if (i <= np) then
    do ii = 1 , np  !Though this looping structure is strange, it allows the kernel to access contiguous np*nlev chunks of v.
      ind = (ii-1)*np*nlev+ik
      ind_k = (ind-1)/(np*np) + 1
      ind_ij = modulo(ind-1,np*np) + 1
      v(ind,q,nt,ie) = v_s((ind_k-1)*(np*np+1) + ind_ij)
    enddo
  endif
end subroutine edgeVunpack_kernel_stage



attributes(global) subroutine hypervis_kernel1( Qdp , qtens , dp , dinv , variable_hyperviscosity , dpdiss_ave , spheremp , deriv_dvv , hyai , hybi , ps0 , nets , nete , dt , nt , nu_p )
  implicit none
  real(kind=real_kind), dimension(np,np,nlev,qsize_d,timelevels,nets:nete), intent(in   ) :: Qdp
  real(kind=real_kind), dimension(np,np,nlev,qsize_d           ,nets:nete), intent(  out) :: Qtens
  real(kind=real_kind), dimension(np,np,nlev                   ,nets:nete), intent(in   ) :: dp
  real(kind=real_kind), dimension(np,np,4                      ,nets:nete), intent(in   ) :: dinv
  real(kind=real_kind), dimension(np,np                                  ), intent(in   ) :: deriv_dvv
  real(kind=real_kind), dimension(      nlev+1                           ), intent(in   ) :: hyai
  real(kind=real_kind), dimension(      nlev+1                           ), intent(in   ) :: hybi
  real(kind=real_kind), dimension(np,np                        ,nets:nete), intent(in   ) :: variable_hyperviscosity
  real(kind=real_kind), dimension(np,np,nlev                   ,nets:nete), intent(in   ) :: dpdiss_ave
  real(kind=real_kind), dimension(np,np                        ,nets:nete), intent(in   ) :: spheremp
  real(kind=real_kind), value                                             , intent(in   ) :: dt , ps0 , nu_p
  integer, value                                                          , intent(in   ) :: nets , nete , nt
  real(kind=real_kind), dimension(np*np+1,numk_hyp,2), shared :: s
  real(kind=real_kind), dimension(np*np+1,4  ), shared :: dinv_s
  real(kind=real_kind), dimension(np,np), shared :: spheremp_s
  real(kind=real_kind), dimension(np,np), shared :: deriv_dvv_s
  real(kind=real_kind) :: dp0
  integer :: i, j, k, q, ie, kk, ks, ij

  ks = int(ceiling(dble(nlev)/numk_hyp))

  i  = modulo( threadidx%x-1    ,np)+1
  j  = modulo((threadidx%x-1)/np,np)+1
  kk = (threadidx%x-1)/(np*np)+1
  k  = modulo(blockidx%x-1,ks)*numk_hyp + kk
  q  = modulo((blockidx%x-1)/ks,qsize_d)+1
  ie =       ((blockidx%x-1)/ks)/qsize_d+1
  ij = (j-1)*np+i

  if (k  > nlev   ) return
  if (q  > qsize_d) return
  if (ie > nete   ) return

  if (kk == 1) then
    dinv_s(ij,:) = dinv(i,j,:,ie)
    spheremp_s(i,j) = spheremp(i,j,ie)
    deriv_dvv_s(i,j) = deriv_dvv(i,j)
  endif
  if (nu_p>0) then
    s(ij,kk,1) = dpdiss_ave(i,j,k,ie) * Qdp(i,j,k,q,nt,ie) / dp(i,j,k,ie)
  else
    dp0 = ( ( hyai(k+1) - hyai(k) )*ps0 + ( hybi(k+1) - hybi(k) )*ps0 )
    s(ij,kk,1) = dp0 * Qdp(i,j,k,q,nt,ie) / dp(i,j,k,ie)
  endif
  call syncthreads()
  qtens(i,j,k,q,ie) = laplace_sphere_wk(i,j,ie,kk,s,dinv_s,spheremp_s,variable_hyperviscosity,deriv_dvv_s,nets,nete)
end subroutine hypervis_kernel1



attributes(global) subroutine hypervis_kernel2( Qdp , qtens , dp , dinv , variable_hyperviscosity , spheremp , rspheremp , deriv_dvv , nu_q , nets , nete , dt , nt )
  implicit none
  real(kind=real_kind), dimension(np,np,nlev,qsize_d,timelevels,nets:nete), intent(inout) :: Qdp
  real(kind=real_kind), dimension(np,np,nlev,qsize_d           ,nets:nete), intent(in   ) :: Qtens
  real(kind=real_kind), dimension(np,np,nlev                   ,nets:nete), intent(in   ) :: dp
  real(kind=real_kind), dimension(np,np                                  ), intent(in   ) :: deriv_dvv
  real(kind=real_kind), dimension(np,np,4                      ,nets:nete), intent(in   ) :: dinv
  real(kind=real_kind), dimension(np,np                        ,nets:nete), intent(in   ) :: variable_hyperviscosity
  real(kind=real_kind), dimension(np,np                        ,nets:nete), intent(in   ) :: spheremp
  real(kind=real_kind), dimension(np,np                        ,nets:nete), intent(in   ) :: rspheremp
  real(kind=real_kind), value                                             , intent(in   ) :: dt, nu_q
  integer, value                                                          , intent(in   ) :: nets , nete , nt
  real(kind=real_kind), dimension(np*np+1,numk_hyp,2), shared :: s
  real(kind=real_kind), dimension(np*np+1,4  ), shared :: dinv_s
  real(kind=real_kind), dimension(np,np), shared :: spheremp_s
  real(kind=real_kind), dimension(np,np), shared :: deriv_dvv_s
  integer :: i, j, k, q, ie, kk, ks, ij

  ks = int(ceiling(dble(nlev)/numk_hyp))

  i  = modulo( threadidx%x-1    ,np)+1
  j  = modulo((threadidx%x-1)/np,np)+1
  kk = (threadidx%x-1)/(np*np)+1
  k  = modulo(blockidx%x-1,ks)*numk_hyp + kk
  q  = modulo((blockidx%x-1)/ks,qsize_d)+1
  ie =       ((blockidx%x-1)/ks)/qsize_d+1
  ij = (j-1)*np+i

  if (k  > nlev   ) return
  if (q  > qsize_d) return
  if (ie > nete   ) return

  if (kk == 1) then
    dinv_s(ij,:) = dinv(i,j,:,ie)
    spheremp_s(i,j) = spheremp(i,j,ie)
    deriv_dvv_s(i,j) = deriv_dvv(i,j)
  endif
  s(ij,kk,1) = rspheremp(i,j,ie)*qtens(i,j,k,q,ie)
  call syncthreads()
  Qdp(i,j,k,q,nt,ie) = Qdp(i,j,k,q,nt,ie)*spheremp_s(i,j)-dt*nu_q*laplace_sphere_wk(i,j,ie,kk,s,dinv_s,spheremp_s,variable_hyperviscosity,deriv_dvv_s,nets,nete)
end subroutine hypervis_kernel2



attributes(device) function laplace_sphere_wk(i,j,ie,k,s,dinv,spheremp,variable_hyperviscosity,deriv_dvv,nets,nete) result(lapl)
  implicit none
  integer,                                              intent(in) :: nets, nete, i, j, ie, k
  real(kind=real_kind), dimension(np*np+1,numk_hyp,2) , intent(inout) :: s
  real(kind=real_kind), dimension(np*np+1,2,2        ), intent(in) :: dinv
  real(kind=real_kind), dimension(np,np              ), intent(in) :: deriv_dvv
  real(kind=real_kind), dimension(np,np,nets:nete    ), intent(in) :: variable_hyperviscosity
  real(kind=real_kind), dimension(np,np              ), intent(in) :: spheremp
  real(kind=real_kind)                                             :: lapl
  integer :: l
  real(kind=real_kind) :: dsdx00 , dsdy00, tmp1, tmp2
  real(kind=real_kind), dimension(2) :: ds
  ds = gradient_sphere(i,j,ie,k,s,dinv,variable_hyperviscosity,deriv_dvv,nets,nete)
  lapl = divergence_sphere_wk(i,j,ie,k,ds,s,dinv,spheremp,deriv_dvv,nets,nete)
end function laplace_sphere_wk



attributes(device) function gradient_sphere(i,j,ie,k,s,dinv,variable_hyperviscosity,deriv_dvv,nets,nete) result(ds)
  implicit none
  integer,                                              intent(in) :: nets, nete, i, j, ie, k
  real(kind=real_kind), dimension(np*np+1,numk_hyp,2) , intent(in) :: s
  real(kind=real_kind), dimension(np*np+1,2,2        ), intent(in) :: dinv
  real(kind=real_kind), dimension(np,np              ), intent(in) :: deriv_dvv
  real(kind=real_kind), dimension(np,np,nets:nete    ), intent(in) :: variable_hyperviscosity
  real(kind=real_kind), dimension(2)                               :: ds
  integer :: l
  real(kind=real_kind) :: dsdx00 , dsdy00, tmp1, tmp2
  dsdx00 = 0.0d0
  dsdy00 = 0.0d0
  do l = 1 , np
    dsdx00 = dsdx00 + deriv_dvv(l,i)*s((j-1)*np+l,k,1)
    dsdy00 = dsdy00 + deriv_dvv(l,j)*s((l-1)*np+i,k,1)
  enddo
  ds(1) = ( dinv((j-1)*np+i,1,1)*dsdx00 + dinv((j-1)*np+i,2,1)*dsdy00 ) * rrearth_d * variable_hyperviscosity(i,j,ie)
  ds(2) = ( dinv((j-1)*np+i,1,2)*dsdx00 + dinv((j-1)*np+i,2,2)*dsdy00 ) * rrearth_d * variable_hyperviscosity(i,j,ie)
end function gradient_sphere



attributes(device) function divergence_sphere_wk(i,j,ie,k,tmp,s,dinv,spheremp,deriv_dvv,nets,nete) result(lapl)
  implicit none
  integer,                                              intent(in   ) :: nets, nete, i, j, ie, k
  real(kind=real_kind), dimension(2),                   intent(in   ) :: tmp
  real(kind=real_kind), dimension(np*np+1,numk_hyp,2) , intent(inout) :: s
  real(kind=real_kind), dimension(np*np+1,2,2        ), intent(in   ) :: dinv
  real(kind=real_kind), dimension(np,np              ), intent(in   ) :: deriv_dvv
  real(kind=real_kind), dimension(np,np              ), intent(in   ) :: spheremp
  real(kind=real_kind)                                                :: lapl
  integer :: l, ij
  ij = (j-1)*np+i
  s(ij,k,1) = ( dinv((j-1)*np+i,1,1)*tmp(1) + dinv((j-1)*np+i,1,2)*tmp(2) ) * spheremp(i,j)
  s(ij,k,2) = ( dinv((j-1)*np+i,2,1)*tmp(1) + dinv((j-1)*np+i,2,2)*tmp(2) ) * spheremp(i,j)
  call syncthreads()
  lapl = 0.0d0
  do l = 1 , np
    lapl = lapl - (s((j-1)*np+l,k,1)*deriv_dvv(i,l)+&
                   s((l-1)*np+i,k,2)*deriv_dvv(j,l)) * rrearth_d
  enddo
end function divergence_sphere_wk



  subroutine limiter_optim_iter_full_loc(ptens,sphweights,minp,maxp,dpmass)
    !THIS IS A NEW VERSION OF LIM8, POTENTIALLY FASTER BECAUSE INCORPORATES KNOWLEDGE FROM
    !PREVIOUS ITERATIONS
    
    !The idea here is the following: We need to find a grid field which is closest
    !to the initial field (in terms of weighted sum), but satisfies the min/max constraints.
    !So, first we find values which do not satisfy constraints and bring these values
    !to a closest constraint. This way we introduce some mass change (addmass),
    !so, we redistribute addmass in the way that l2 error is smallest. 
    !This redistribution might violate constraints thus, we do a few iterations. 
    implicit none
    real (kind=real_kind), dimension(np*np,nlev,qsize,nelemd), intent(inout) , device :: ptens
    real (kind=real_kind), dimension(np*np           ,nelemd), intent(in   ) , device :: sphweights
    real (kind=real_kind), dimension(      nlev,qsize,nelemd), intent(inout) , device :: minp
    real (kind=real_kind), dimension(      nlev,qsize,nelemd), intent(inout) , device :: maxp
    real (kind=real_kind), dimension(np*np,nlev      ,nelemd), intent(in   ) , device :: dpmass
 
    integer :: k1, i, j, k, iter, i1, i2, q, ie
    integer :: whois_neg(np*np), whois_pos(np*np)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!! OPENACC WORKAROUND !!!!!!!!!!!!!!!!!!!!!!!
    ! These should be integers, but for some strange reason, the
    ! OpenACC kernel gives wrong answers unless you make them reals
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real (kind=real_kind) :: neg_counter, pos_counter
    real (kind=real_kind) :: addmass, weightssum, mass,sumc,tmp,min_tmp,max_tmp
    real (kind=real_kind) :: x(np*np),c(np*np)
    real (kind=real_kind) :: al_neg(np*np), al_pos(np*np), howmuch
    real (kind=real_kind) :: tol_limiter = 1e-15
    integer :: maxiter = 5

    !$acc  parallel loop gang vector collapse(3) &
    !$acc&   private(k1,i,j,k,iter,i1,i2,q,ie,whois_neg,whois_pos,neg_counter, &
    !$acc&           pos_counter,addmass,weightssum,mass,x,c,al_neg,al_pos,howmuch,sumc,tmp) &
    !$acc&   deviceptr(ptens,sphweights,minp,maxp,dpmass) vector_length(512) async(1)
    do ie = 1 , nelemd
      do q = 1 , qsize
        do k = 1 , nlev
          !$acc cache(c,x,al_neg,al_pos,whois_neg,whois_pos)
          min_tmp = minp(k,q,ie)
          max_tmp = maxp(k,q,ie)
          c = sphweights(:,ie) * dpmass(:,k,ie)
          x = ptens(:,k,q,ie) / dpmass(:,k,ie)

          sumc = 0d0
          mass = 0d0
          !$acc loop seq
          do i1 = 1 , np*np
            mass = mass + c(i1)*x(i1)
            sumc = sumc + c(i1)
          enddo

          ! relax constraints to ensure limiter has a solution:
          ! This is only needed if runnign with the SSP CFL>1 or 
          ! due to roundoff errors
          if( (mass / sumc) < min_tmp ) then
            min_tmp = mass / sumc
          endif
          if( (mass / sumc) > max_tmp ) then
            max_tmp = mass / sumc
          endif

          addmass = 0.0d0
          pos_counter = 0
          neg_counter = 0
          
          ! apply constraints, compute change in mass caused by constraints 
          !$acc loop seq
          do k1 = 1 , np*np
            if ( ( x(k1) >= max_tmp ) ) then
              addmass = addmass + ( x(k1) - max_tmp ) * c(k1)
              x(k1) = max_tmp
              whois_pos(k1) = -1
            else
              pos_counter = pos_counter+1
              whois_pos(pos_counter) = k1
            endif
            if ( ( x(k1) <= min_tmp ) ) then
              addmass = addmass - ( min_tmp - x(k1) ) * c(k1)
              x(k1) = min_tmp
              whois_neg(k1) = -1
            else
              neg_counter = neg_counter+1
              whois_neg(neg_counter) = k1
            endif
          enddo
          
          ! iterate to find field that satifies constraints and is l2-norm closest to original 
          weightssum = 0.0d0
          if ( addmass > 0 ) then
            !$acc loop seq
            do i2 = 1 , maxIter
              weightssum = 0.0
              !$acc loop seq
              do k1 = 1 , pos_counter
                i1 = whois_pos(k1)
                weightssum = weightssum + c(i1)
                al_pos(i1) = max_tmp - x(i1)
              enddo
              
              if( ( pos_counter > 0 ) .and. ( addmass > tol_limiter * abs(mass) ) ) then
                !$acc loop seq
                do k1 = 1 , pos_counter
                  i1 = whois_pos(k1)
                  howmuch = addmass / weightssum
                  if ( howmuch > al_pos(i1) ) then
                    howmuch = al_pos(i1)
                    whois_pos(k1) = -1
                  endif
                  tmp = addmass - howmuch * c(i1)
                  addmass = tmp
                  weightssum = weightssum - c(i1)
                  x(i1) = x(i1) + howmuch
                enddo
                !now sort whois_pos and get a new number for pos_counter
                !here neg_counter and whois_neg serve as temp vars
                neg_counter = pos_counter
                whois_neg = whois_pos
                whois_pos = -1
                pos_counter = 0
                !$acc loop seq
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
            !$acc loop seq
            do i2 = 1 , maxIter
              weightssum = 0.0
              !$acc loop seq
              do k1 = 1 , neg_counter
                i1 = whois_neg(k1)
                weightssum = weightssum + c(i1)
                al_neg(i1) = x(i1) - min_tmp
              enddo
              
              if ( ( neg_counter > 0 ) .and. ( (-addmass) > tol_limiter * abs(mass) ) ) then
                !$acc loop seq
                do k1 = 1 , neg_counter
                  i1 = whois_neg(k1)
                  howmuch = -addmass / weightssum
                  if ( howmuch > al_neg(i1) ) then
                    howmuch = al_neg(i1)
                    whois_neg(k1) = -1
                  endif
                  tmp = addmass + howmuch * c(i1)
                  addmass = tmp
                  weightssum = weightssum - c(i1)
                  x(i1) = x(i1) - howmuch
                enddo
                !now sort whois_pos and get a new number for pos_counter
                !here pos_counter and whois_pos serve as temp vars
                pos_counter = neg_counter
                whois_pos = whois_neg
                whois_neg = -1
                neg_counter = 0
                !$acc loop seq
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
          minp(k,q,ie) = min_tmp
          maxp(k,q,ie) = max_tmp
        enddo
      enddo
    enddo
  end subroutine limiter_optim_iter_full_loc



  subroutine neighbor_minmax_loc(elem,hybrid,nets,nete,min_neigh,max_neigh)
    ! compute Q min&max over the element and all its neighbors
    use element_mod, only: element_t
    use hybrid_mod, only: hybrid_t
    use edge_mod, only: edgeVpack, edgeVunpackMin, edgeVunpackMax
    use bndry_mod, only: bndry_exchangeV
    use perf_mod, only: t_startf, t_stopf
    implicit none
    integer                      , intent(in   ) :: nets,nete
    type (element_t)             , intent(in   ) :: elem(:)
    type (hybrid_t)              , intent(in   ) :: hybrid
    real (kind=real_kind), device, intent(inout) :: min_neigh(nlev,qsize,nelemd)
    real (kind=real_kind), device, intent(inout) :: max_neigh(nlev,qsize,nelemd)
    ! local
    integer :: ie,kk,q,ks,i,j,k
    integer, parameter :: kchunk = 8
    real (kind=real_kind) :: min_tmp(kchunk), max_tmp(kchunk)
    if (hybrid%ithr == 0) then
      !$acc parallel loop gang collapse(3) deviceptr(min_neigh,max_neigh,qmin_exch,qmax_exch) private(min_tmp,max_tmp) async(1) vector_length(np*np*kchunk)
      do ie=1,nelemd
        do q=1,qsize
          do ks=1,nlev/kchunk+1
            !$acc cache(min_tmp,max_tmp)
            !$acc loop vector
            do kk=1,kchunk
              k = (ks-1)*kchunk+kk
              if (k <= nlev) then
                min_tmp(kk) = min_neigh(k,q,ie)
                max_tmp(kk) = max_neigh(k,q,ie)
              endif
            enddo
            !$acc loop collapse(3) vector
            do kk=1,kchunk
              do j=1,np
                do i=1,np
                  k = (ks-1)*kchunk+kk
                  if (k <= nlev) then
                    qmin_exch(i,j,k,q,ie) = min_tmp(kk)
                    qmax_exch(i,j,k,q,ie) = max_tmp(kk)
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      call edgeVpack_gpu(edgebuf_Q2_d,qmin_exch,nlev*qsize*2,nlev*qsize,0         ,putmapP_d,reverse_d)
      call edgeVpack_gpu(edgebuf_Q2_d,qmax_exch,nlev*qsize*2,nlev*qsize,nlev*qsize,putmapP_d,reverse_d)
      !$acc wait(1)
    endif
  
    call t_startf('nmm_bexchV')
    call bndry_exchangeV_loc(hybrid,edgebuf_Q2_d,recvbuf_Q2_d,edgebuf_Q2_h,recvbuf_Q2_h,nlev*qsize*2)
    call t_stopf('nmm_bexchV')
       
    if (hybrid%ithr == 0) then
      call edgeVunpackMin_gpu(edgebuf_Q2_d,qmin_exch,nlev*qsize*2,nlev*qsize,0         ,getmapP_d)
      call edgeVunpackMax_gpu(edgebuf_Q2_d,qmax_exch,nlev*qsize*2,nlev*qsize,nlev*qsize,getmapP_d)
      !$acc parallel loop gang vector collapse(3) deviceptr(min_neigh,max_neigh,qmin_exch,qmax_exch) async(1)
      do ie=1,nelemd
        do q=1,qsize
          do k=1,nlev
            ! note: only need to consider the corners, since the data we packed was constant within each element
            min_neigh(k,q,ie)=max(min(qmin_exch(1,1,k,q,ie),qmin_exch(1,np,k,q,ie),qmin_exch(np,1,k,q,ie),qmin_exch(np,np,k,q,ie)),0d0)
            max_neigh(k,q,ie)=    max(qmax_exch(1,1,k,q,ie),qmax_exch(1,np,k,q,ie),qmax_exch(np,1,k,q,ie),qmax_exch(np,np,k,q,ie))
          enddo
        enddo
      enddo
    endif
  end subroutine neighbor_minmax_loc



  subroutine biharmonic_wk_scalar_minmax_loc(elem,qtens,deriv,hybrid,nets,nete,emin,emax)
    use element_mod, only: element_t
    use derivative_mod, only: derivative_t
    use hybrid_mod, only: hybrid_t
    use derivative_mod, only: laplace_sphere_wk
    use control_mod, only : hypervis_scaling, hypervis_power
    use edge_mod, only: edgeVpack, edgeVunpack, edgeVunpackMin, edgeVunpackMax
    use perf_mod, only: t_startf, t_stopf
    use bndry_mod, only: bndry_exchangeV
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! compute weak biharmonic operator
    !    input:  qtens = Q
    !    output: qtens = weak biharmonic of Q and Q element min/max
    !    note: emin/emax must be initialized with Q element min/max.  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    type (element_t)             , intent(inout), target :: elem(:)
    real (kind=real_kind), device, intent(inout)         :: qtens(np,np,nlev,qsize,nelemd)
    type (derivative_t)          , intent(in   )         :: deriv
    type (hybrid_t)              , intent(in   )         :: hybrid
    integer                      , intent(in   )         :: nets,nete
    real (kind=real_kind), device, intent(  out)         :: emin(nlev,qsize,nelemd)
    real (kind=real_kind), device, intent(  out)         :: emax(nlev,qsize,nelemd)
    integer :: k,kptr,i,j,ie,ic,q,ks,kk
    real (kind=real_kind) :: lap_p(np,np)
    logical :: var_coef1
    integer, parameter :: kchunk = 8
    real (kind=real_kind) :: min_tmp(kchunk), max_tmp(kchunk)
    !$acc routine(laplace_sphere_wk_loc) seq
    !if tensor hyperviscosity with tensor V is used, then biharmonic operator is (\grad\cdot V\grad) (\grad \cdot \grad) 
    !so tensor is only used on second call to laplace_sphere_wk
    var_coef1 = .true.
    if(hypervis_scaling > 0) var_coef1 = .false.
    if (hybrid%ithr == 0) then
      !$acc enter data pcopyin(hypervis_power,hypervis_scaling,var_coef1)
      !$acc parallel loop gang collapse(3) deviceptr(emin,emax,qmin_exch,qmax_exch) private(min_tmp,max_tmp) async(1) vector_length(np*np*kchunk)
      do ie=1,nelemd
        do q=1,qsize
          do ks=1,nlev/kchunk+1
            !$acc cache(min_tmp,max_tmp)
            !$acc loop vector
            do kk=1,kchunk
              k = (ks-1)*kchunk+kk
              if (k <= nlev) then
                min_tmp(kk) = emin(k,q,ie)
                max_tmp(kk) = emax(k,q,ie)
              endif
            enddo
            !$acc loop collapse(3) vector
            do kk=1,kchunk
              do j=1,np
                do i=1,np
                  k = (ks-1)*kchunk+kk
                  if (k <= nlev) then
                    qmin_exch(i,j,k,q,ie) = min_tmp(kk)
                    qmax_exch(i,j,k,q,ie) = max_tmp(kk)
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      !$acc  parallel loop gang vector collapse(3) deviceptr(emin,emax,qtens,qmin_exch,qmax_exch) private(lap_p) &
      !$acc& present(deriv_dvv_acc,variable_hyperviscosity_h,dinv_acc,tensorvisc_acc,spheremp_h,hypervis_power,hypervis_scaling,var_coef1) async(1)
      do ie=1,nelemd
        do q=1,qsize      
          do k=1,nlev    !  Potential loop inversion (AAM)
            lap_p = qtens(:,:,k,q,ie)
            call laplace_sphere_wk_loc(lap_p,deriv_dvv_acc,dinv_acc(:,:,:,:,ie),spheremp_h(:,:,ie),tensorvisc_acc(:,:,:,:,ie),hypervis_power,hypervis_scaling,variable_hyperviscosity_h,var_coef1,lap_p)
            qtens(:,:,k,q,ie) = lap_p
          enddo
        enddo
      enddo
      call edgeVpack_gpu(edgebuf_Q3_d,    qtens,nlev*qsize*3,nlev*qsize,           0,putmapP_d,reverse_d)
      call edgeVpack_gpu(edgebuf_Q3_d,qmin_exch,nlev*qsize*3,nlev*qsize,  nlev*qsize,putmapP_d,reverse_d)
      call edgeVpack_gpu(edgebuf_Q3_d,qmax_exch,nlev*qsize*3,nlev*qsize,2*nlev*qsize,putmapP_d,reverse_d)
      !$acc wait(1)
    endif

    call t_startf('biwkscmm_bexchV')
    call bndry_exchangeV_loc(hybrid,edgebuf_Q3_d,recvbuf_Q3_d,edgebuf_Q3_h,recvbuf_Q3_h,nlev*qsize*3)
    call t_stopf('biwkscmm_bexchV')

    if (hybrid%ithr == 0) then
      call edgeVunpack_gpu   (edgebuf_Q3_d,    qtens,nlev*qsize*3,qsize*nlev,           0,getmapP_d)
      call edgeVunpackMin_gpu(edgebuf_Q3_d,qmin_exch,nlev*qsize*3,qsize*nlev,  qsize*nlev,getmapP_d)
      call edgeVunpackMax_gpu(edgebuf_Q3_d,qmax_exch,nlev*qsize*3,qsize*nlev,2*qsize*nlev,getmapP_d)
      !$acc  parallel loop gang vector collapse(3) deviceptr(emin,emax,qtens,qmin_exch,qmax_exch,rspheremp_d) private(lap_p) &
      !$acc& present(deriv_dvv_acc,variable_hyperviscosity_h,dinv_acc,tensorvisc_acc,spheremp_h,hypervis_power,hypervis_scaling,var_coef1) async(1)
      do ie=1,nelemd
        do q=1,qsize      
          do k=1,nlev
            lap_p(:,:)=rspheremp_d(:,:,ie)*qtens(:,:,k,q,ie)
            call laplace_sphere_wk_loc(lap_p,deriv_dvv_acc,dinv_acc(:,:,:,:,ie),spheremp_h(:,:,ie),tensorvisc_acc(:,:,:,:,ie),hypervis_power,hypervis_scaling,variable_hyperviscosity_h,.true.,lap_p)
            qtens(:,:,k,q,ie) = lap_p
            emin(k,q,ie)=max(min(qmin_exch(1,1,k,q,ie),qmin_exch(1,np,k,q,ie),qmin_exch(np,1,k,q,ie),qmin_exch(np,np,k,q,ie)),0d0)
            emax(k,q,ie)=    max(qmax_exch(1,1,k,q,ie),qmax_exch(1,np,k,q,ie),qmax_exch(np,1,k,q,ie),qmax_exch(np,np,k,q,ie))
          enddo
        enddo
      enddo
    endif
  end subroutine biharmonic_wk_scalar_minmax_loc



  subroutine edgeVpack_gpu(edge_buf,v,nlyr,vlyr,kptr,putmapP,reverse)
    use dimensions_mod, only : np, max_corner_elem, max_neigh_edges
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    use edge_mod, only: edgeDescriptor_t
    use perf_mod, only: t_startf, t_stopf
    implicit none
    real (kind=real_kind), device , intent(inout) :: edge_buf(nlyr,nbuf)
    real (kind=real_kind), device , intent(in   ) :: v(np,np,vlyr,nelemd)
    integer                       , intent(in   ) :: nlyr
    integer                       , intent(in   ) :: vlyr
    integer                       , intent(in   ) :: kptr
    integer              , device , intent(in   ) :: putmapP(max_neigh_edges,nelemd)
    logical              , device , intent(in   ) :: reverse(max_neigh_edges,nelemd)
    integer :: i,k,ir,ll,ie
    call t_startf('edge_pack')
    !$acc parallel loop gang vector collapse(2) deviceptr(putmapP,reverse,v,edge_buf) private(i,ll,ir) async(1)
    do ie=1,nelemd
      do k=1,vlyr
        do i=1,np
          edge_buf(kptr+k,putmapP(south,ie)+i) = v(i ,1 ,k,ie)
          edge_buf(kptr+k,putmapP(east ,ie)+i) = v(np,i ,k,ie)
          edge_buf(kptr+k,putmapP(north,ie)+i) = v(i ,np,k,ie)
          edge_buf(kptr+k,putmapP(west ,ie)+i) = v(1 ,i ,k,ie)
        enddo
        do i=1,np
          ir = np-i+1
          if(reverse(south,ie)) edge_buf(kptr+k,putmapP(south,ie)+ir) = v(i ,1 ,k,ie)
          if(reverse(east ,ie)) edge_buf(kptr+k,putmapP(east ,ie)+ir) = v(np,i ,k,ie)
          if(reverse(north,ie)) edge_buf(kptr+k,putmapP(north,ie)+ir) = v(i ,np,k,ie)
          if(reverse(west ,ie)) edge_buf(kptr+k,putmapP(west ,ie)+ir) = v(1 ,i ,k,ie)
        enddo
        do i = 0 , max_corner_elem-1
          ll = swest+0*max_corner_elem+i ; if (putmapP(ll,ie) /= -1) edge_buf(kptr+k,putmapP(ll,ie)+1) = v(1  ,1 ,k,ie)
          ll = swest+1*max_corner_elem+i ; if (putmapP(ll,ie) /= -1) edge_buf(kptr+k,putmapP(ll,ie)+1) = v(np ,1 ,k,ie)
          ll = swest+2*max_corner_elem+i ; if (putmapP(ll,ie) /= -1) edge_buf(kptr+k,putmapP(ll,ie)+1) = v(1  ,np,k,ie)
          ll = swest+3*max_corner_elem+i ; if (putmapP(ll,ie) /= -1) edge_buf(kptr+k,putmapP(ll,ie)+1) = v(np ,np,k,ie)
        enddo
      enddo
    enddo
    call t_stopf('edge_pack')
  end subroutine edgeVpack_gpu



  subroutine edgeVunpack_gpu(edge_buf,v,nlyr,vlyr,kptr,getmapP)
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    use perf_mod, only: t_startf, t_stopf
    implicit none
    real (kind=real_kind), device , intent(in   ) :: edge_buf(nlyr,nbuf)
    real (kind=real_kind), device , intent(inout) :: v(np,np,vlyr,nelemd)
    integer                       , intent(in   ) :: nlyr
    integer                       , intent(in   ) :: vlyr
    integer                       , intent(in   ) :: kptr
    integer              , device , intent(in   ) :: getmapP(max_neigh_edges,nelemd)
    integer :: i,k,ll,ie,kc,kk,j,loc_ind,glob_k,ij
    integer, parameter :: kchunk = 64
    real(kind=real_kind) :: vtmp(np*np+1,kchunk)
    call t_startf('edge_unpack')
    !$acc parallel loop gang collapse(2) deviceptr(getmapP,v,edge_buf) private(ll,i,vtmp) async(1) vector_length(kchunk)
    do ie = 1 , nelemd
      do kc = 1 , vlyr/kchunk+1
        !$acc cache(vtmp)
        !$acc loop vector
        do kk = 1 , kchunk
          do ll = 1 , np*np
            loc_ind = (ll-1)*kchunk+kk-1
            i = modulo(loc_ind,np)+1
            j = modulo(loc_ind/np,np)+1
            ij = (j-1)*np+i
            k = loc_ind/np/np+1
            glob_k = (kc-1)*kchunk+k
            if (glob_k <= vlyr) vtmp(ij,k) = v(i,j,glob_k,ie)
          enddo
        enddo
        !$acc loop vector
        do kk = 1 , kchunk
          k = (kc-1)*kchunk+kk
          if (k <= vlyr) then
            do i = 1 , np
              vtmp(i +(1 -1)*np,kk) = vtmp(i +(1 -1)*np,kk) + edge_buf(kptr+k,getmapP(south,ie)+i)
              vtmp(np+(i -1)*np,kk) = vtmp(np+(i -1)*np,kk) + edge_buf(kptr+k,getmapP(east ,ie)+i)
              vtmp(i +(np-1)*np,kk) = vtmp(i +(np-1)*np,kk) + edge_buf(kptr+k,getmapP(north,ie)+i)
              vtmp(1 +(i -1)*np,kk) = vtmp(1 +(i -1)*np,kk) + edge_buf(kptr+k,getmapP(west ,ie)+i)
            enddo
            do i = 0 , max_corner_elem-1
              ll = swest+0*max_corner_elem+i; if(getmapP(ll,ie) /= -1) vtmp(1 +(1 -1)*np,kk) = vtmp(1 +(1 -1)*np,kk) + edge_buf(kptr+k,getmapP(ll,ie)+1)
              ll = swest+1*max_corner_elem+i; if(getmapP(ll,ie) /= -1) vtmp(np+(1 -1)*np,kk) = vtmp(np+(1 -1)*np,kk) + edge_buf(kptr+k,getmapP(ll,ie)+1)
              ll = swest+2*max_corner_elem+i; if(getmapP(ll,ie) /= -1) vtmp(1 +(np-1)*np,kk) = vtmp(1 +(np-1)*np,kk) + edge_buf(kptr+k,getmapP(ll,ie)+1)
              ll = swest+3*max_corner_elem+i; if(getmapP(ll,ie) /= -1) vtmp(np+(np-1)*np,kk) = vtmp(np+(np-1)*np,kk) + edge_buf(kptr+k,getmapP(ll,ie)+1)
            enddo
          endif
        enddo
        !$acc loop vector
        do kk = 1 , kchunk
          do ll = 1 , np*np
            loc_ind = (ll-1)*kchunk+kk-1
            i = modulo(loc_ind,np)+1
            j = modulo(loc_ind/np,np)+1
            ij = (j-1)*np+i
            k = loc_ind/np/np+1
            glob_k = (kc-1)*kchunk+k
            if (glob_k <= vlyr) v(i,j,glob_k,ie) = vtmp(ij,k)
          enddo
        enddo
      enddo
    enddo
    call t_stopf('edge_unpack')
  end subroutine edgeVunpack_gpu



  subroutine edgeVunpackMIN_gpu(edge_buf,v,nlyr,vlyr,kptr,getmapP)
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    implicit none
    real (kind=real_kind), device , intent(in   ) :: edge_buf(nlyr,nbuf)
    real (kind=real_kind), device , intent(inout) :: v(np,np,vlyr,nelemd)
    integer                       , intent(in   ) :: nlyr
    integer                       , intent(in   ) :: vlyr
    integer                       , intent(in   ) :: kptr
    integer              , device , intent(in   ) :: getmapP(max_neigh_edges,nelemd)
    integer :: i,k,ll,ie
    !$acc parallel loop gang vector collapse(2) deviceptr(getmapP,v,edge_buf) private(ll,i) async(1)
    do ie = 1 , nelemd
      do k = 1 , vlyr
        do i = 1 , np
          v(i ,1 ,k,ie) = min( v(i ,1 ,k,ie) , edge_buf(kptr+k,getmapP(south,ie)+i) )
          v(np,i ,k,ie) = min( v(np,i ,k,ie) , edge_buf(kptr+k,getmapP(east ,ie)+i) )
          v(i ,np,k,ie) = min( v(i ,np,k,ie) , edge_buf(kptr+k,getmapP(north,ie)+i) )
          v(1 ,i ,k,ie) = min( v(1 ,i ,k,ie) , edge_buf(kptr+k,getmapP(west ,ie)+i) )
        enddo
        do i = 0 , max_corner_elem-1
          ll = swest+0*max_corner_elem+i; if(getmapP(ll,ie) /= -1) v(1 ,1 ,k,ie) = min( v(1 ,1 ,k,ie) , edge_buf(kptr+k,getmapP(ll,ie)+1) )
          ll = swest+1*max_corner_elem+i; if(getmapP(ll,ie) /= -1) v(np,1 ,k,ie) = min( v(np,1 ,k,ie) , edge_buf(kptr+k,getmapP(ll,ie)+1) )
          ll = swest+2*max_corner_elem+i; if(getmapP(ll,ie) /= -1) v(1 ,np,k,ie) = min( v(1 ,np,k,ie) , edge_buf(kptr+k,getmapP(ll,ie)+1) )
          ll = swest+3*max_corner_elem+i; if(getmapP(ll,ie) /= -1) v(np,np,k,ie) = min( v(np,np,k,ie) , edge_buf(kptr+k,getmapP(ll,ie)+1) )
        enddo
      enddo
    enddo
  end subroutine edgeVunpackMIN_gpu



  subroutine edgeVunpackMAX_gpu(edge_buf,v,nlyr,vlyr,kptr,getmapP)
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    implicit none
    real (kind=real_kind), device , intent(in   ) :: edge_buf(nlyr,nbuf)
    real (kind=real_kind), device , intent(inout) :: v(np,np,vlyr,nelemd)
    integer                       , intent(in   ) :: nlyr
    integer                       , intent(in   ) :: vlyr
    integer                       , intent(in   ) :: kptr
    integer              , device , intent(in   ) :: getmapP(max_neigh_edges,nelemd)
    integer :: i,k,ll,ie
    !$acc parallel loop gang vector collapse(2) deviceptr(getmapP,v,edge_buf) private(ll,i) async(1)
    do ie = 1 , nelemd
      do k = 1 , vlyr
        do i = 1 , np
          v(i ,1 ,k,ie) = max( v(i ,1 ,k,ie) , edge_buf(kptr+k,getmapP(south,ie)+i) )
          v(np,i ,k,ie) = max( v(np,i ,k,ie) , edge_buf(kptr+k,getmapP(east ,ie)+i) )
          v(i ,np,k,ie) = max( v(i ,np,k,ie) , edge_buf(kptr+k,getmapP(north,ie)+i) )
          v(1 ,i ,k,ie) = max( v(1 ,i ,k,ie) , edge_buf(kptr+k,getmapP(west ,ie)+i) )
        enddo
        do i = 0 , max_corner_elem-1
          ll = swest+0*max_corner_elem+i; if(getmapP(ll,ie) /= -1) v(1 ,1 ,k,ie) = max( v(1 ,1 ,k,ie) , edge_buf(kptr+k,getmapP(ll,ie)+1) )
          ll = swest+1*max_corner_elem+i; if(getmapP(ll,ie) /= -1) v(np,1 ,k,ie) = max( v(np,1 ,k,ie) , edge_buf(kptr+k,getmapP(ll,ie)+1) )
          ll = swest+2*max_corner_elem+i; if(getmapP(ll,ie) /= -1) v(1 ,np,k,ie) = max( v(1 ,np,k,ie) , edge_buf(kptr+k,getmapP(ll,ie)+1) )
          ll = swest+3*max_corner_elem+i; if(getmapP(ll,ie) /= -1) v(np,np,k,ie) = max( v(np,np,k,ie) , edge_buf(kptr+k,getmapP(ll,ie)+1) )
        enddo
      enddo
    enddo
  end subroutine edgeVunpackMAX_gpu



  subroutine laplace_sphere_wk_loc(s,deriv_dvv,dinv,spheremp,tensorvisc,hypervis_power,hypervis_scaling,variable_hyperviscosity,var_coef,laplace)
    use element_mod, only: element_t
    use derivative_mod, only: derivative_t
    implicit none
    !$acc routine seq
!   input:  s = scalar
!   ouput:  -< grad(PHI), grad(s) >   = weak divergence of grad(s)
!     note: for this form of the operarad(s) does not need to be made C0
    real(kind=real_kind), intent(in)   :: s(np,np) 
    real(kind=real_kind), intent(in)   :: dinv(np,np,2,2)
    real(kind=real_kind), intent(in)   :: tensorvisc(2,2,np,np)
    real(kind=real_kind), intent(in)   :: spheremp(np,np)
    real(kind=real_kind), intent(in)   :: hypervis_power,hypervis_scaling
    real(kind=real_kind), intent(in)   :: variable_hyperviscosity(np,np)
    logical             , intent(in)   :: var_coef
    real(kind=real_kind), intent(in)   :: deriv_dvv(np,np)
    real(kind=real_kind), intent(out)  :: laplace(np,np)
    real(kind=real_kind) :: laplace2(np,np)
    integer i,j
    real(kind=real_kind) :: grads(np,np,2), oldgrads(np,np,2)
    !$acc routine(gradient_sphere_loc) seq
    !$acc routine(divergence_sphere_wk_loc) seq
    call gradient_sphere_loc(s,deriv_dvv,dinv,grads)
    if (var_coef) then
       if (hypervis_power/=0 ) then
          ! scalar viscosity with variable coefficient
          do j = 1 , np
            do i = 1 , np
              grads(i,j,1) = grads(i,j,1)*variable_hyperviscosity(i,j)
              grads(i,j,2) = grads(i,j,2)*variable_hyperviscosity(i,j)
            enddo
          enddo
       else if (hypervis_scaling /=0 ) then
          ! tensor hv, (3)
          do j=1,np
             do i=1,np
                grads(i,j,1) = sum(grads(i,j,:)*tensorVisc(1,:,i,j))
                grads(i,j,2) = sum(grads(i,j,:)*tensorVisc(2,:,i,j))
             end do
          end do
       else
          ! do nothing: constant coefficient viscsoity
       endif
    endif
    ! note: divergnece_sphere and divergence_sphere_wk are identical *after* bndry_exchange
    ! if input is C_0.  Here input is not C_0, so we should use divergence_sphere_wk().  
    call divergence_sphere_wk_loc(grads,deriv_dvv,dinv,spheremp,laplace)
  end subroutine laplace_sphere_wk_loc



  subroutine divergence_sphere_wk_loc(v,deriv_dvv,dinv,spheremp,div)
    use physical_constants, only: rrearth
    use element_mod, only: element_t
    use derivative_mod, only: derivative_t
    implicit none
    !$acc routine seq
!   input:  v = velocity in lat-lon coordinates
!   ouput:  div(v)  spherical divergence of v, integrated by parts
!   Computes  -< grad(psi) dot v > 
!   (the integrated by parts version of < psi div(v) > )
!   note: after DSS, divergence_sphere() and divergence_sphere_wk() 
!   are identical to roundoff, as theory predicts.
    real(kind=real_kind), intent(in   ) :: v(np,np,2)  ! in lat-lon coordinates
    real(kind=real_kind), intent(in   ) :: deriv_dvv(np,np)
    real(kind=real_kind), intent(in   ) :: Dinv(np,np,2,2)
    real(kind=real_kind), intent(in   ) :: spheremp(np,np)
    real(kind=real_kind), intent(  out) :: div(np,np)
    integer i,j,m,n
    real(kind=real_kind) :: vtemp(np,np,2)
    real(kind=real_kind) :: ggtemp(np,np,2)
    real(kind=real_kind) :: gtemp(np,np,2)
    real(kind=real_kind) :: psi(np,np)
    real(kind=real_kind) :: xtmp
    do j=1,np
       do i=1,np
          vtemp(i,j,1)=(Dinv(i,j,1,1)*v(i,j,1) + Dinv(i,j,1,2)*v(i,j,2))
          vtemp(i,j,2)=(Dinv(i,j,2,1)*v(i,j,1) + Dinv(i,j,2,2)*v(i,j,2))
       enddo
    enddo
    do n=1,np
       do m=1,np
          div(m,n)=0
          do j=1,np
             div(m,n)=div(m,n)-(  spheremp(j,n)*vtemp(j,n,1)*deriv_dvv(m,j) &
                                + spheremp(m,j)*vtemp(m,j,2)*deriv_dvv(n,j) ) * rrearth
          enddo
       end do
    end do
  end subroutine divergence_sphere_wk_loc



  subroutine gradient_sphere_loc(s,deriv_dvv,dinv,ds)
    use physical_constants, only: rrearth
    use derivative_mod, only: derivative_t
    implicit none
    !$acc routine seq
!   input s:  scalar
!   output  ds: spherical gradient of s, lat-lon coordinates
    real(kind=real_kind), intent(in   ) :: deriv_dvv(np,np)
    real(kind=real_kind), intent(in   ) :: Dinv(np,np,2,2)
    real(kind=real_kind), intent(in   ) :: s(np,np)
    real(kind=real_kind), intent(  out) :: ds(np,np,2)
    integer :: i, j, l
    real(kind=real_kind) ::  dsdx00
    real(kind=real_kind) ::  dsdy00
    real(kind=real_kind) ::  v1(np,np),v2(np,np)
    do j=1,np
       do l=1,np
          dsdx00=0.0d0
          dsdy00=0.0d0
          do i=1,np
             dsdx00 = dsdx00 + deriv_dvv(i,l)*s(i,j)
             dsdy00 = dsdy00 + deriv_dvv(i,l)*s(j,i)
          end do
          v1(l,j) = dsdx00*rrearth
          v2(j,l) = dsdy00*rrearth
       end do
    end do
    ! convert covarient to latlon
    do j=1,np
       do i=1,np
          ds(i,j,1)=Dinv(i,j,1,1)*v1(i,j) + Dinv(i,j,2,1)*v2(i,j)
          ds(i,j,2)=Dinv(i,j,1,2)*v1(i,j) + Dinv(i,j,2,2)*v2(i,j)
       enddo
    enddo
  end subroutine gradient_sphere_loc



  subroutine bndry_exchangeV_loc(hybrid,sendbuf_dev,recvbuf_dev,sendbuf_hst,recvbuf_hst,nlyr)
    use hybrid_mod, only : hybrid_t
    use kinds, only : log_kind
    use schedule_mod, only : schedule_t, cycle_t, schedule
    use dimensions_mod, only: nelemd, np
#ifdef _MPI
    use parallel_mod, only : abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success
#else
    use parallel_mod, only : abortmp
#endif
    implicit none
    type (hybrid_t)              , intent(in   ) :: hybrid
    real (kind=real_kind), device, intent(inout) :: sendbuf_dev(nlyr,nbuf)
    real (kind=real_kind), device, intent(inout) :: recvbuf_dev(nlyr,nbuf)
    real (kind=real_kind)        , intent(inout) :: sendbuf_hst(nlyr,nbuf)
    real (kind=real_kind)        , intent(inout) :: recvbuf_hst(nlyr,nbuf)
    integer                      , intent(in   ) :: nlyr
    type (Schedule_t),pointer     :: pSchedule
    type (Cycle_t),pointer        :: pCycle
    integer                       :: dest,length,tag
    integer                       :: icycle,ierr
    integer                       :: iptr,source
    integer                       :: nSendCycles,nRecvCycles
    integer                       :: errorcode,errorlen
    character*(80) :: errorstring
    integer        :: i
    if(hybrid%ithr == 0) then 
#ifdef _MPI
      ! Setup the pointer to proper Schedule
#ifdef _PREDICT
      pSchedule => Schedule(iam)
#else
      pSchedule => Schedule(1)
#endif
      nSendCycles = pSchedule%nSendCycles
      nRecvCycles = pSchedule%nRecvCycles

      !==================================================
      !  Post the Receives 
      !==================================================
      do icycle=1,nRecvCycles
        pCycle      => pSchedule%RecvCycle(icycle)
        source      =  pCycle%source - 1
        length      =  nlyr * pCycle%lengthP
        tag         =  pCycle%tag
        iptr        =  pCycle%ptrP
        call MPI_Irecv(recvbuf_hst(1,iptr),length,MPIreal_t, source,tag,hybrid%par%comm,Rrequest(icycle),ierr)
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
        pCycle      => pSchedule%SendCycle(icycle)
        iptr        =  pCycle%ptrP
        if (pCycle%lengthP > 0) ierr = cudaMemCpyAsync(sendbuf_hst(1,iptr),sendbuf_dev(1,iptr),size(sendbuf_hst(1:nlyr,iptr:iptr+pCycle%lengthP-1)),cudaMemcpyDeviceToHost,streams(1))
      enddo
      ierr = cudaStreamSynchronize(streams(1))
      do icycle=1,nSendCycles
        pCycle      => pSchedule%SendCycle(icycle)
        dest        =  pCycle%dest - 1
        length      =  nlyr * pCycle%lengthP
        tag         =  pCycle%tag
        iptr        =  pCycle%ptrP
        call MPI_Isend(sendbuf_hst(1,iptr),length,MPIreal_t,dest,tag,hybrid%par%comm,Srequest(icycle),ierr)
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
        pCycle         => pSchedule%RecvCycle(icycle)
        length         =  pCycle%lengthP
        iptr           =  pCycle%ptrP
        if (length > 0) ierr = cudaMemCpyAsync(sendbuf_dev(1,iptr),recvbuf_hst(1,iptr),size(recvbuf_hst(1:nlyr,iptr:iptr+length-1)),cudaMemcpyHostToDevice,streams(1))
      enddo   ! icycle
      ierr = cudaStreamSynchronize(streams(1))
#endif
    endif  ! if (hybrid%ithr == 0)
    !$OMP BARRIER
  end subroutine bndry_exchangeV_loc



#endif
end module cuda_mod



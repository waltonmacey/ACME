
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module openacc_utils_mod
  use kinds, only: real_kind
  use dimensions_mod, only: nelemd
  implicit none
  private

  public :: copy_qdp_h2d
  public :: copy_qdp_d2h
  public :: copy_arr
  public :: arch_init2

contains

  subroutine arch_init2( elem , deriv )
    use element_mod, only: element_t, state_qdp, derived_vn0, derived_divdp, derived_divdp_proj
    use derivative_mod, only: derivative_t
    implicit none
    type(element_t)   , intent(in) :: elem(:)
    type(derivative_t), intent(in) :: deriv
    !$omp barrier
    !$omp master

    !$acc enter data pcreate(state_Qdp,derived_vn0,derived_divdp,derived_divdp_proj)
    !$acc enter data pcopyin(elem(1:nelemd),deriv)

    !$omp end master
    !$omp barrier
  end subroutine arch_init2

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

end module openacc_utils_mod


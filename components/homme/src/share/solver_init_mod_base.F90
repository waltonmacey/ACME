#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module solver_init_mod_base
  implicit none
  private

  public :: solver_init2


contains


  subroutine solver_init2( elem , hvcoord , deriv )
    use element_mod, only: element_t
    use derivative_mod, only: derivative_t
    use hybvcoord_mod, only : hvcoord_t
    implicit none
    type(element_t)   , intent(in) :: elem(:)
    type(derivative_t), intent(in) :: deriv
    type(hvcoord_t)   , intent(in) :: hvcoord
    !do nothing
  end subroutine solver_init2


end module solver_init_mod_base

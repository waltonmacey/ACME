
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module arch_switch_mod
#if USE_OPENACC
  use prim_advection_openacc_mod, only: prim_advec_init1, prim_advec_init2, Prim_Advec_Tracers_remap_rk2
#else
  use prim_advection_mod, only: prim_advec_init1, prim_advec_init2, Prim_Advec_Tracers_remap_rk2
#endif
  implicit none
end module arch_switch_mod


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_advance_mod
  use prim_advance_mod_base, only: prim_advance_exp, prim_advance_si, prim_advance_init, preq_robert3, &
                                   applyCAMforcing_dynamics, applyCAMforcing, smooth_phis, overwrite_SEdensity, prim_step_init, init_dp3d
  implicit none
end module prim_advance_mod


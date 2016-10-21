#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_advance_mod

  use prim_advance_base, only:&
    applyCAMforcing_dynamics, &
    applyCAMforcing,          &
    compute_and_apply_rhs,    &
    init_dynamics,            &
    overwrite_SEdensity,      &
    preq_robert3,             &
    prim_advance_exp,         &
    prim_advance_init,        &
    prim_advance_si,          &
    smooth_phis

end module prim_advance_mod

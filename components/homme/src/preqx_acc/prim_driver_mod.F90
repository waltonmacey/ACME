#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_driver_mod

  use prim_driver_base, only:&
    leapfrog_bootstrap, &
    prim_init1,&
    prim_init2 ,&
    prim_finalize,&
    prim_run, &
    prim_run_subcycle

  implicit none
  contains
end module

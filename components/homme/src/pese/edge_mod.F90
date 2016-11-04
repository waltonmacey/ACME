#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module edge_mod
  use edge_mod_base, only: &
    FreeEdgeBuffer, &
    FreeGhostBuffer3D, &
    FreeGhostBufferTR,&
    FreeLongEdgeBuffer, &
    LongEdgeVpack, &
    LongEdgeVunpackMIN, &
    buffermap,&
    edgeDGVpack,&
    edgeDGVunpack,&
    edgeDefaultVal,&
    edgeSpack, &
    edgeSunpackMax, &
    edgeSunpackMin, &
    edgeVpack, &
    edgeVunpack,&
    edgeVunpackMAX,&
    edgeVunpackMIN,&
    edgeVunpackVert,&
    edgerotate,&
    ghostVpack, &
    ghostVpack2d, &
    ghostVpack2d_level,&
    ghostVpack2d_single, &
    ghostVpack3d, &
    ghostVpackR, &
    ghostVpack_unoriented, &
    ghostVpackfull,&
    ghostVunpack, &
    ghostVunpack2d,&
    ghostVunpack2d_level,&
    ghostVunpack2d_single, &
    ghostVunpack3d, &
    ghostVunpackR, &
    ghostVunpack_unoriented,&
    ghostVunpackfull,&
    initEdgeBuffer, &
    initEdgeSBuffer, &
    initGhostBuffer3D,&
    initGhostBufferTR,&
    initLongEdgeBuffer
  implicit none
end module edge_mod

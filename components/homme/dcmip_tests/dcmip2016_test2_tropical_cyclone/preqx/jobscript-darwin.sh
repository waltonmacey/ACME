#!/bin/bash
#
#   Jobscript for launching dcmip2016 test2 on a mac running Darwin
#

EXEC=../../../test_execs/preqx-nlev30-interp/preqx-nlev30-interp        # set name of executable
openmpiexec -n 6 $EXEC < ./namelist-lowres.nl                          # launch simulation

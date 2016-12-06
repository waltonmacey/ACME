#!/bin/bash
#
#   Jobscript for launching dcmip2016 test2 on a mac running Darwin
#

EXEC=../../../test_execs/pese-nlev30/pese-nlev30       # set name of executable
openmpiexec -n 6 $EXEC < ./namelist-lowres.nl                          # launch simulation

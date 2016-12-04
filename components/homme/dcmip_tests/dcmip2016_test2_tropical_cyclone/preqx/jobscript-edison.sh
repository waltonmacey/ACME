#!/bin/bash
#
#   Jobscript for launching dcmip2016 test2 on the NERSC Edison machine
#
#SBATCH -J dcmip16-2          # job name
#SBATCH -o out_dcmip16-1.o%j  # output and error file name (%j expands to jobID)
#SBATCH -n 600                # total number of mpi tasks requested
#SBATCH -p debug              # queue (partition) -- normal, development, etc.
#SBATCH -t 00:10:00           # run time (hh:mm:ss)
#SBATCH -A acme               # charge hours to account 1

EXEC=../../../test_execs/preqx-nlev30-interp/preqx-nlev30-interp        # set name of executable
srun -n 600 $EXEC < ./namelist-default.nl                              # launch simulation


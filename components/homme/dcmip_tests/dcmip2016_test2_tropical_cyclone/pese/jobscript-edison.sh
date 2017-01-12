#!/bin/bash
#
#   Jobscript for launching dcmip2016 test2 on the NERSC Edison machine
#
#SBATCH -J dcmip16-2          # job name
#SBATCH -o out_dcmip16-1.o%j  # output and error file name (%j expands to jobID)
#SBATCH -n 900                # total number of mpi tasks requested
#SBATCH -p regular            # queue (partition) -- normal, development, etc.
#SBATCH --qos=premium
#SBATCH -t 01:30:00           # run time (hh:mm:ss)
#SBATCH -A m2618 #acme        # charge hours to account 1

EXEC=../../../test_execs/pese-nlev30/pese-nlev30                        # set name of executable
srun -n 900 $EXEC < ./namelist-default.nl                              # launch simulation


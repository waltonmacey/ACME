#!/bin/bash
cwd=`pwd`
cd ../../..
  echo "make -j 4 preqx-nlev30-interp"
  make -j 4 preqx-nlev30-interp
cd $cwd
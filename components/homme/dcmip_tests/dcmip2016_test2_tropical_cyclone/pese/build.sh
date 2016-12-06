#!/bin/bash
cwd=`pwd`
cd ../../..
  echo "make -j 4 pese-nlev30"
  make -j 4 pese-nlev30
cd $cwd
#!/bin/bash
# script to run toulbar2 in parallel within a debugger for each process
# usage: mpirun.sh [nproc] ./bin/Linux/toulbar2 file.wcsp [other options]"

if [ "$1" -eq "$1" ]
then
  NPROC=$1
  shift
  mpirun -n $NPROC xterm -e gdb -ex run --args $*
else
  mpirun xterm -e gdb -ex run --args $*
fi

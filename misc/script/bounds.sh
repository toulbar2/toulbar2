#!/bin/tcsh

# usage: toulbar2 problem | bounds.sh 1

${0:h}/ts "%.s" - | awk -v N=$1 -f ${0:h}/bounds.awk

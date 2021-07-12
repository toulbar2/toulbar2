#!/bin/sh

# usage: toulbar2 problem | bounds.sh 1

./misc/script/ts "%.s" - | awk -v N=$1 -f ./misc/script/bounds.awk


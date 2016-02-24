#!/bin/sh

# usage: toulbar2 problem | bounds.sh 1

ts "%.s" - | awk -v N=$1 -f ./misc/script/bounds.awk


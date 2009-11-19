#!/bin/sh

# Provides a variable elimination ordering to ToulBar2

# Usage: td.sh problem.wcsp toulbar2_options

# by default, uses Maximum Cardinality Search heuristic
peo $1 1 > order
toulbar2 $1 $2 Oorder

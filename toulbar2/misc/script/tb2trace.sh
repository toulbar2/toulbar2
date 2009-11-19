#!/bin/sh

# usage: tb2trace.sh trace.txt

# warning! still buggy...

# for debugging purposes,
# extract from toulbar2 verbose trace the messages 
# corresponding to the current search path only.

# for instance, the following input trace:

#[1,26,35,44] Try 1 == 1
#[2,29,35,40] Try 21 == 0
#[3,32,35,36] Try 9 == 1
#[3,32,35,36] Try 9 == 0
#[2,29,35,40] Try 21 == 1
#[3,31,35,36] Try 9 == 1
#[4,31,35,32] Try 7 == 0
#[5,33,35,26] Try 4 == 0
#New solution: 34 (1 backtracks, 8 nodes, depth 5)

# will result in the following output:

#[1][1,26,35,44] Try 1 == 1
#[2][2,29,35,40] Try 21 == 1
#[3][3,31,35,36] Try 9 == 1
#[4][4,31,35,32] Try 7 == 0
#[5][5,33,35,26] Try 4 == 0
#[5]New solution: 34 (1 backtracks, 8 nodes, depth 5)

awk 'BEGIN{depth=1} /^[[][0-9]*/{d=$1;sub("[[]","",d);depth=0+d} {print "[" depth "]" $0} /contradiction/{depth--}' $1 | tac | awk -f tb2trace.awk | tac

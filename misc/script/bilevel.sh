#!/bin/tcsh

# compute all solutions of problem1.cfn and problem2.cfn and return their combination (tuples joined on the first N fields) sorted by decreasing cost of sol1 minus cost of sol2

#usage:
# ./misc/scripts/bilevel.sh problem1.cfn problem2.cfn N

./bin/Linux/toulbar2 $1 -s=3 -a -svo -nopre -hbfs: -k=0 | awk -v N=$3 '/solution.*[)]:/{sub("solution[(]","",$2); for (i=1;i<=N;i++) printf(":%s",$(2+i)); for (j=3+N;j<=NF;j++) printf(" %s",$j); printf(" %d",0+$2); print ""}' | sort >! ${1:r}.all

set T = `head -1 ${1:r}.all | wc -w`

./bin/Linux/toulbar2 $2 -s=3 -a -svo -nopre -hbfs: -k=0 | awk -v N=$3 '/solution.*[)]:/{sub("solution[(]","",$2); for (i=1;i<=N;i++) printf(":%s",$(2+i)); for (j=3+N;j<=NF;j++) printf(" %s",$j); printf(" %d",0+$2); print ""}' | sort >! ${2:r}.all

join ${1:r}.all ${2:r}.all | awk '{print $'$T'-$NF,$0}' | sed 's/:/ /g' | sort -g -r -k 1,1


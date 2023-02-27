#!/bin/tcsh

# compute all solutions of problem1.cfn and problem2.cfn and return their combination (tuples joined on the first N fields) sorted by decreasing cost of problem1 minus the minimum cost of problem2

#usage:
# ./misc/scripts/bilevel.sh problem1.cfn problem2.cfn N

#note: if N is negative then output all combinations without selecting the best one of problem2

if ( $# != 3 ) then
   echo "usage: $0 problem1.cfn problem2.cfn NumberOfCommonVariables"
   exit(0)
endif

set N = 0
if ( $3 < 0 ) then
   @ N -= $3
else
   @ N += $3
endif

toulbar2 $1 -s=3 -a -svo -nopre -hbfs: -k=0 | awk -v N=$N '/solution.*[)]:/{sub("solution[(]","",$2); for (i=1;i<=N;i++) printf(":%s",$(2+i)); for (j=3+N;j<=NF;j++) printf(" %s",$j); printf(" %d",0+$2); print ""}' | sort >! ${1:r}.all

set T = `head -1 ${1:r}.all | wc -w`

toulbar2 $2 -s=3 -a -svo -nopre -hbfs: -k=0 | awk -v N=$N '/solution.*[)]:/{sub("solution[(]","",$2); for (i=1;i<=N;i++) printf(":%s",$(2+i)); for (j=3+N;j<=NF;j++) printf(" %s",$j); printf(" %d",0+$2); print ""}' | sort >! ${2:r}.all

set F = `head -1 ${2:r}.all | wc -w`

set P1 = $T
@ P1 += $N
@ P1 -= 1

set P2 = $P1
@ P2 += $F

if ( $3 < 0 ) then
  join ${1:r}.all ${2:r}.all | awk '{print $'$T' - $NF, $0}' | sed 's/:/ /g' | sort -g -r -k 1,1
else
  join ${1:r}.all ${2:r}.all | awk '{print $'$T' - $NF, $0}' | sed 's/:/ /g' | sort -k 2,$P1 -k $P2,${P2}gr | awk '{o=$1;$1="";s=$0;sub(" [-0-9]+.*","",s);res[s]=o "" $0} END{for (e in res) {print res[e]}}' | sort -g -r -k 1,1
endif


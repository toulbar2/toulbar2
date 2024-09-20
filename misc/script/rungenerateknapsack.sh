#!/bin/bash

# Usage:
# rungenerateknapsack.sh "options"

nerr=0
ninstances=100
nend=10
bctr=3
tctr=0
nary=1
tight=80
n=5
d=3
K=1

while (( $n < $nend )) ; do
  seed=0
  echo n:$n d:$d tight:$tight%  binary:$bctr  ternary:$tctr  quintary:$nary  
  while (( $seed < $ninstances )) ; do  
    # test toulbar2
    rm -f toulbar2_opt
    rm -f toulbar2_verif
    rm -f sol
    randomfile="knapsack-$n-$d-$tight-$bctr-$tctr-0-$nary"
#    echo $randomfile
    ./toulbar2 -random=$randomfile -seed=$seed -C=$K -nopre -k=0 -ub=1000000 -z=1 > /dev/null
    python2 ./misc/script/wcsp2lp-support.py problem.wcsp problem.lp > /dev/null
    ./toulbar2 problem.wcsp "${@:1}" -w | awk 'BEGIN{opt=-1;} /No solution/{opt=-2} /^Optimum: /{opt=$2;} END{printf("%d",opt);}' > toulbar2_opt
    ub1=`awk 'BEGIN{opt=-1;} {opt=$1} END{printf("%d", opt)}' toulbar2_opt`
    if [[ $ub1 -gt -2 ]] ; then
    if [[ $ub1 -lt 0 ]] ; then
      echo "error $nerr found negative toulbar2 optimum! $ub1"
      mv problem.wcsp error$nerr.wcsp
      nerr=`expr $nerr + 1`
      seed=`expr $seed + 1`
	  continue
    fi
    ./toulbar2 problem.wcsp -x | awk 'BEGIN{opt=-1;} /nb. of unassigned variables: 0/{ sub("[[]","",$0); opt=$4;} END{printf("%d",opt);}' > toulbar2_verif
    ub1b=`awk 'BEGIN{opt=-1;} {opt=$1} END{printf("%d", opt)}' toulbar2_verif`	
    if [[ $ub1b -lt 0 ]] ; then
      echo "error $nerr found negative certificate upperbound! $ub1b"
      mv problem.wcsp error$nerr.wcsp
      nerr=`expr $nerr + 1`
      seed=`expr $seed + 1`
	  continue
    fi
    if [[ $ub1 != $ub1b ]] ; then
      echo "error $nerr found diff between optimum and certificate! $ub1 $ub1b"
      mv problem.wcsp error$nerr.wcsp
      nerr=`expr $nerr + 1`
      seed=`expr $seed + 1`
	  continue
    fi
    fi

    # compare with cplex
    rm -f toulbar2_cplex
    rm -f cplex.log
    cplex -c read cplex.prm read problem.lp opt | awk 'BEGIN{opt=-1;} /Integer infeasible/{opt=-2} /Integer optimal solution:/{opt=0+$NF} END{printf("%d", opt)}' > toulbar2_cplex
    ub2=`awk 'BEGIN{opt=-1;} {opt=$1} END{printf("%d", opt)}' toulbar2_cplex`
    if [[ $ub2 -gt -2 ]] ; then
    if [[ $ub2 -lt 0 ]] ; then
      echo "error $nerr found negative cplex optimum! $ub2"
      mv problem.lp error$nerr.lp
      nerr=`expr $nerr + 1`
      seed=`expr $seed + 1`
	  continue
    fi
    fi
    if [[ $ub1 != $ub2 ]] ; then
      echo "error $nerr found diff between optimum of toulbar2 and cplex! $ub1 $ub2"
      mv problem.wcsp error$nerr.wcsp
      mv problem.lp error$nerr.lp
      nerr=`expr $nerr + 1`
      seed=`expr $seed + 1`
	  continue
    fi
    seed=`expr $seed + 1`
  done	
  nary=`expr $nary + 1` 
  tctr=`expr $tctr + 1`  
  bctr=`expr $bctr + 3`  
  n=`expr $n + 1`
done

rm -f problem.wcsp
rm -f problem.lp
rm -f toulbar2_opt
rm -f sol
rm -f problem.dot
rm -f problem.degree
rm -f toulbar2_cplex
rm -f cplex.log
rm -f toulbar2_opt
rm -f toulbar2_verif

#!/bin/bash

# Usage:
# rungeneratealldiff.sh "options"

nerr=0
ninstances=100
nend=10
bctr=10
tctr=0
nary="1"
tight=80
n=4
d=4
K=1

while (( $n < $nend )) ; do
  seed=0
  echo n:$n d:$d tight:$tight%  binary:$bctr  ternary:$tctr  quatary:$nary  
  while (( $seed < $ninstances )) ; do
  
    rm -f toulbar2_opt
    rm -f toulbar2_verif
    rm -f sol
    # tests hungarian-based propagation
    randomfile="alldiff-$n-$d-$tight-$bctr-$tctr-$nary-$seed"
    ./Debug/bin/Linux/toulbar2 -random=$randomfile -C=$K -nopre -k=0 -z -v=-1 > /dev/null
    cp problem.wcsp problemHUN.wcsp
    ./Debug/bin/Linux/toulbar2 problem.wcsp "${@:2}" -w | awk 'BEGIN{opt=-1;} /^Optimum: /{opt=$2;} END{printf("%d",opt);}' > toulbar2_opt
    ub1=`awk 'BEGIN{opt=-1;} {opt=$1} END{printf("%d", opt)}' toulbar2_opt`
    ./Debug/bin/Linux/toulbar2 problem.wcsp -x | awk 'BEGIN{opt=-1;} /nb. of unassigned variables: 0/{ sub("[[]","",$0); opt=$4;} END{printf("%d",opt);}' > toulbar2_verif
    ub1b=`awk 'BEGIN{opt=-1;} {opt=$1} END{printf("%d", opt)}' toulbar2_verif`
    if [[ $ub1 -lt 0 ]] ; then
      echo "error found $ub1 < 0"
      mv problem.wcsp error$nerr.wcsp
      nerr=`expr $nerr + 1`
    fi
    if [[ $ub1b -lt 0 ]] ; then
      echo "error found $ub1b < 0"
      mv problem.wcsp error$nerr.wcsp
      nerr=`expr $nerr + 1`
    fi
    if [[ $ub1 != $ub1b ]] ; then
      echo "error found $ub1 != $ub1b"
      mv problem.wcsp error$nerr.wcsp
      nerr=`expr $nerr + 1`
    fi

    # tests knapsack-based propagation
    rm -f toulbar2_opt
    rm -f toulbar2_verif
    rm -f sol
    randomfile="salldiffkp-$n-$d-$tight-$bctr-$tctr-$nary-$seed"
    ./Debug/bin/Linux/toulbar2 -random=$randomfile -C=$K -nopre -k=0 -z > /dev/null
    cp problem.wcsp problemKP.wcsp
    ./Debug/bin/Linux/toulbar2 problem.wcsp $2 -w | awk 'BEGIN{opt=-1;} /^Optimum: /{opt=$2;} END{printf("%d",opt);}' > toulbar2_opt
    ub2=`awk 'BEGIN{opt=-1;} {opt=$1} END{printf("%d", opt)}' toulbar2_opt`
    ./Debug/bin/Linux/toulbar2 problem.wcsp -x | awk 'BEGIN{opt=-1;} /nb. of unassigned variables: 0/{ sub("[[]","",$0); opt=$4;} END{printf("%d",opt);}' > toulbar2_verif
    ub2b=`awk 'BEGIN{opt=-1;} {opt=$1} END{printf("%d", opt)}' toulbar2_verif`
    if [[ $ub2 -lt 0 ]] ; then
      echo "error found $ub2 < 0"
      mv problem.wcsp errorKP$nerr.wcsp
      nerr=`expr $nerr + 1`
    fi
    if [[ $ub2b -lt 0 ]] ; then
      echo "error found $ub2b < 0"
      mv problem.wcsp errorKP$nerr.wcsp
      nerr=`expr $nerr + 1`
    fi
    if [[ $ub2 != $ub2b ]] ; then
      echo "error found $ub2 != $ub2b"
      mv problem.wcsp errorKP$nerr.wcsp
      nerr=`expr $nerr + 1`
    fi

    if [[ $ub1 != $ub2 ]] ; then
      echo "error found between hungarian and knapsack at seed $seed"
    fi

    if [[ $ub1 != $ub2 ]] ; then
      mv problemHUN.wcsp error$nerr.wcsp
      mv problemKP.wcsp errorKP$nerr.wcsp
      nerr=`expr $nerr + 1`
    fi
    
    seed=`expr $seed + 1`
  done	
  nary="0-$nary" 
  tctr=`expr $tctr + 5`  
  bctr=`expr $bctr + 10`  
  n=`expr $n + 1`
  d=$n
done

rm -f problem.wcsp
rm -f toulbar2_opt
rm -f toolbar_opt
rm -f toolbar_sol
rm -f sol
rm -f problem.dot
rm -f problem.degree

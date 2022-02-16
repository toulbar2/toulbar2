#!/bin/bash

nerr=0
ninstances=100
nend=20
bctr=20
tctr=0
nary=1
nary5=1
tight=90  # replace by 4% in order to generate clauses, i.e., with at most one non-zero tuple per 4-ary and 5-ary constraints
n=5
d=3
K=1

while (( $n < $nend )) ; do
  seed=0
  echo n:$n d:$d tight:$tight%  binary:$bctr  ternary:$tctr  quatary:$nary  quintary: $nary5
  while (( $seed < $ninstances )) ; do
    rm -f toulbar2_opt
    rm -f toolbar_opt
    rm -f toolbar_sol
    rm -f sol
    randomfile=nary-$n-$d-$tight-$bctr-$tctr-$nary-$nary5-$seed 
    ./toulbar2 -random=$randomfile -C=$K -z > /dev/null
    toolbar problem.wcsp  | awk 'BEGIN{opt="-";} /^Optimum: /{opt=$2;}  END{printf("%d \n",opt); }' > toolbar_opt
    ub0=`cat toolbar_opt`
    ub=`expr $ub0 + 2`
    ./toulbar2 problem.wcsp "$@" -w -ub=$ub | awk 'BEGIN{opt="-";} /Read [0-9]* variables/{n=$2;} / unassigned variables/{f=$1} /^Optimum: /{opt=$2;}  END{printf("%d %d %d",opt,n,f);}'  > toulbar2_opt
    toolbar problem.wcsp  -csol  | awk 'BEGIN{opt="-";} /^Total cost /{opt=$4;}  END{printf("%d \n",opt); }' > toolbar_sol
    ub1=`awk '{printf("%d", $1)}' toulbar2_opt`
    ub2=`awk '{printf("%d", $1)}' toolbar_opt`
    ub3=`awk '{printf("%d", $1)}' toolbar_sol`
    sz1=`awk '{printf("%d", $2)}' toulbar2_opt`
    sz2=`awk '{printf("%d", $3)}' toulbar2_opt`
#    if  [[ $sz1 != $sz2 ]] ; then
#      cat toulbar2_opt
#      echo " " $randomfile
#    fi
    if [[ $ub1 != $ub2 ]] ; then
      echo "error found"
      mv problem.wcsp error$nerr.wcsp
      nerr=`expr $nerr + 1`
    fi
    if [[ $ub1 != $ub3 ]] ; then
      echo "solution error found"
    fi
    seed=`expr $seed + 1`
  done	
  nary=`expr $nary + 2` 
  nary5=`expr $nary5 + 1` 
  tctr=`expr $tctr + 3`  
  bctr=`expr $bctr + 5`  
  n=`expr $n + 1`
done

rm -f problem.wcsp
rm -f toulbar2_opt
rm -f toolbar_opt
rm -f toolbar_sol
rm -f sol
rm -f problem.dot
rm -f problem.degree

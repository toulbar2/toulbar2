#!/bin/bash

nerr=0
ninstances=100
nend=10
bctr=10
tctr=1
nary=1
tight=40
n=6
d=5
K=100

while (( $n < $nend )) ; do
  seed=0
  echo n:$n m:$d tight:$tight%  binary:$bctr  ternary:$tctr  quatary:$nary  
  while (( $seed < $ninstances )) ; do
    rm -f toulbar2_opt
    rm -f toolbar_opt
    rm -f toolbar_sol
    rm -f order
    rm -f sol
    randomfile=nary-$n-$d-$tight-$bctr-$tctr-$nary-$seed 
    toulbar2 $randomfile C$K z > /dev/null
    ./misc/bin/Linux/peo problem.wcsp 1 > order
    ./bin/Linux/toulbar2 problem.wcsp $1 w Oorder | awk 'BEGIN{opt="-";} /^Optimum: /{opt=$2;}  END{printf("%d ",opt); }'   > toulbar2_opt
    toolbar problem.wcsp  | awk 'BEGIN{opt="-";} /^Optimum: /{opt=$2;}  END{printf("%d \n",opt); }' > toolbar_opt
    toolbar problem.wcsp  -csol  | awk 'BEGIN{opt="-";} /^Total cost /{opt=$4;}  END{printf("%d \n",opt); }' > toolbar_sol
    ub1=`cat toulbar2_opt`
    ub2=`cat toolbar_opt`
    ub3=`cat toolbar_sol`
    if [[ $ub1 != $ub2 ]] ; then
      echo "error found"
      mv problem.wcsp error$nerr.wcsp
      mv order        o$nerr
      nerr=`expr $nerr + 1`
    fi
    if [[ $ub1 != $ub3 ]] ; then
      echo "solution error found"
    fi
    seed=`expr $seed + 1`
  done	
  nary=`expr $nary + 1` 
  tctr=`expr $tctr + 1`  
  bctr=`expr $bctr + 10`  
  n=`expr $n + 1`
done

rm -f problem.wcsp
rm -f toulbar2_opt
rm -f toolbar_opt
rm -f toolbar_sol
rm -f order
rm -f sol

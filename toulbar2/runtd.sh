#!/bin/sh


nerr=0
ninstances=50
nend=99
bctr=4
tctr=0
nary=0
tight=80
n=5
m=4
h=2


while (( $n < $nend )) ; do
  seed=0
  echo n:$n m:$m tight:$tight%  binary:$bctr  ternary:$tctr  quatary:$nary  quintary:$nary 
  #echo  n:$n m:$m tight:$tight%  h:$h 
  while (( $seed < $ninstances )) ; do
    rm -f toulbar2_opt
    rm -f toolbar_opt
    randomfile=nary-$n-$m-$tight-$bctr-$tctr-$nary-$nary-$seed 
    toulbar2 $randomfile z > /dev/null
    peo problem.wcsp 1 > order
    toulbar2 problem.wcsp $1 Oorder | awk 'BEGIN{opt="-";} /^Optimum: /{opt=$2;}  END{printf("%d ",opt); }'   > toulbar2_opt
    toolbar problem.wcsp | awk 'BEGIN{opt="-";} /^Optimum: /{opt=$2;}  END{printf("%d \n",opt); }' > toolbar_opt
    ub1=`cat toulbar2_opt`
    ub2=`cat toolbar_opt`
    if [[ $ub1 != $ub2 ]] ; then
      echo "error found"
      mv problem.wcsp error$nerr.wcsp
      mv order        o$nerr
      nerr=`expr $nerr + 1`
    fi
    seed=`expr $seed + 1`
  done	
  nary=`expr $nary + 2` 
  bctr=`expr $bctr + 3` 
  n=`expr $n + 1`
  #h=`expr $h + 1`
done

rm -f problem.wcsp
rm -f toulbar2_opt
rm -f toolbar_opt
  
 

#!/bin/sh


nerr=0
ninstances=50
nend=100
bctr=50
tctr=50
tight=40
n=3
m=5


while (( $n < $nend )) ; do
  seed=0
  while (( $seed < $ninstances )) ; do
    rm -f toulbar2_opt
    rm -f toolbar_opt
    randomfile=ternary$n-$m-$tight-$bctr-$tctr-$seed 
    /tmp/toulbar2 $randomfile z > /dev/null
    /tmp/toulbar2 problem.wcsp $1 | awk 'BEGIN{opt="-";} /^Optimum: /{opt=$2;}  END{printf("%d ",opt); }'   > toulbar2_opt
    /tmp/toolbar  problem.wcsp    | awk 'BEGIN{opt="-";} /^Optimum: /{opt=$2;}  END{printf("%d \n",opt); }' > toolbar_opt
    ub1=`cat toulbar2_opt`
    ub2=`cat toolbar_opt`
    if [[ $ub1 != $ub2 ]] ; then
      echo "error found"
      mv problem.wcsp error$nerr.wcsp
      nerr=`expr $nerr + 1`
    fi
    seed=`expr $seed + 1`
  done	
  tctr=`expr $tctr + 20` 
  n=`expr $n + 1`
done

rm -f problem.wcsp
rm -f toulbar2_opt
rm -f toolbar_opt
  
 

#!/bin/sh


nerr=0
ninstances=50
nend=100
bctr=10
tctr=10
nary=1
tight=30
n=4
m=4


while (( $n < $nend )) ; do
  seed=0
  echo n:$n m:$m tight:$tight%  binary:$bctr  ternary:$tctr  quatary:$nary  quintary:$nary 
  while (( $seed < $ninstances )) ; do
    rm -f toulbar2_opt
    rm -f toolbar_opt
    randomfile=nary-$n-$m-$tight-$bctr-$tctr-$nary-$nary-$seed 
<<<<<<< rungenerate.sh
    ./toulbar2 $randomfile z > /dev/null
    ./toulbar2 problem.wcsp $1 | awk 'BEGIN{opt="-";} /^Optimum: /{opt=$2;}  END{printf("%d ",opt); }'   > toulbar2_opt
    ./toolbar  problem.wcsp    | awk 'BEGIN{opt="-";} /^Optimum: /{opt=$2;}  END{printf("%d \n",opt); }' > toolbar_opt
=======
    toulbar2 $randomfile z > /dev/null
    toulbar2 problem.wcsp $1 | awk 'BEGIN{opt="-";} /^Optimum: /{opt=$2;}  END{printf("%d ",opt); }'   > toulbar2_opt
    toolbar  problem.wcsp    | awk 'BEGIN{opt="-";} /^Optimum: /{opt=$2;}  END{printf("%d \n",opt); }' > toolbar_opt
>>>>>>> 1.4
    ub1=`cat toulbar2_opt`
    ub2=`cat toolbar_opt`
    if [[ $ub1 != $ub2 ]] ; then
      echo "error found"
      mv problem.wcsp error$nerr.wcsp
      nerr=`expr $nerr + 1`
    fi
    seed=`expr $seed + 1`
  done	
  nary=`expr $nary + 15` 
  n=`expr $n + 1`
done

rm -f problem.wcsp
rm -f toulbar2_opt
rm -f toolbar_opt
  
 

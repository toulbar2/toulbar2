#!/bin/bash

nerr=0
ninstances=100
nend=6
bctr=20
tctr=0
nary=1
nary5=1
tight=90  # replace by 4% in order to generate clauses, i.e., with at most one non-zero tuple per 4-ary and 5-ary constraints
n=3
d=3
K=1
m=$n

while (( $n <= $nend )) ; do
  seed=0
  echo n:$n d:$d tight:$tight%  binary:$bctr  ternary:$tctr  quatary:$nary  quintary: $nary5
  while (( $seed < $ninstances )) ; do
    rm -f toulbar2_opt
    rm -f problem1_opt
    rm -f problem2_opt
    rm -f bilevel_opt
    rm -f sol
#    randomfile=nary-$n-$d-$tight-$bctr-$tctr-$nary-$nary5-$seed 
#    randomfile=tern-$n-$d-$tight-$bctr-$tctr-$seed 
    randomfile=bin-$n-$d-$tight-$bctr-$seed 
    ./toulbar2 -random=$randomfile -C=$K -z=3 -z=problem1.cfn > /dev/null
    m=`expr $n + $n`
    seed2=`expr $seed + 1000000`
#    randomfile=nary-$m-$d-$tight-$bctr-$tctr-$nary-$nary5-$seed2
#    randomfile=tern-$m-$d-$tight-$bctr-$tctr-$seed2
    randomfile=bin-$m-$d-$tight-$bctr-$seed2
    ./toulbar2 -random=$randomfile -C=$K -z=3 -z=problem2.cfn > /dev/null
    ./toulbar2 -bilevel problem1.cfn problem2.cfn "$@" -w=3 | awk 'BEGIN{opt=1e9;} /^Optimum: /{opt=$2;} END{printf("%s",opt);}' > toulbar2_opt
    ./toulbar2 problem1.cfn -x | awk 'BEGIN{opt=-1e9} /Input solution cost:/ && /of unassigned variables: 0/{opt=$4;} END{printf("%s",opt);}' > problem1_opt
    ./toulbar2 problem2.cfn -x | awk 'BEGIN{opt=-1e9} /Input solution cost:/ && /of unassigned variables: 0/{opt=$4;} END{printf("%s",opt);}' > problem2_opt
    ./misc/script/bilevel.sh problem1.cfn problem2.cfn $n | tail -1 | awk '{printf("%s",$1);}' > bilevel_opt
    ub1=`awk '{printf("%d", $1)}' toulbar2_opt`
    problem1=`awk '{printf("%d", $1)}' problem1_opt`
    problem2=`awk '{printf("%d", $1)}' problem2_opt`
    ub2=`expr $problem1 - $problem2`
    ub3=`awk '{printf("%d", $1)}' bilevel_opt`
    #echo "$ub1 $problem1 $problem2 $ub2 $ub3"
    if [[ $ub1 != $ub2 ]] ; then
      echo "solution error found"
    fi
    if [[ $ub1 != $ub3 ]] ; then
      echo "error found"
      mv problem1.cfn error${nerr}_leader.cfn
      mv problem2.cfn error${nerr}_follower.cfn
      nerr=`expr $nerr + 1`
    fi
    seed=`expr $seed + 1`
  done	
  nary=`expr $nary + 2` 
  nary5=`expr $nary5 + 1` 
  tctr=`expr $tctr + 3`  
  bctr=`expr $bctr + 5`  
  n=`expr $n + 1`
done

rm -f problem1.cfn
rm -f problem2.cfn
rm -f problem1.all
rm -f problem2.all
rm -f toulbar2_opt
rm -f problem1_opt
rm -f problem2_opt
rm -f bilevel_opt
rm -f sol


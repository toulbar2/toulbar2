#!/bin/bash

nerr=0
ninstances=100
nend=10
bctr=20
tctr=0
nary=10
tight=80
n=4
d=5
K=1

while (( $n < $nend )) ; do
  seed=0
  echo n:$n d:$d tight:$tight%  binary:$bctr  ternary:$tctr  quatary:$nary  
  while (( $seed < $ninstances )) ; do
  
    # tests salldiff
    rm -f toulbar2_opt
    rm -f toulbar2_verif
    rm -f sol
    randomfile="salldiff-$n-$d-$tight-$bctr-$tctr-$nary-$seed"
    ./toulbar2 -random=$randomfile -C=$K -nopre -k=0 -z > /dev/null
    ./toulbar2 problem.wcsp $1 -w | awk 'BEGIN{opt=-1;} /^Optimum: /{opt=$2;} END{printf("%d",opt);}' > toulbar2_opt
    ub1=`awk 'BEGIN{opt=-1;} {opt=$1} END{printf("%d", opt)}' toulbar2_opt`
    ./toulbar2 problem.wcsp -x | awk 'BEGIN{opt=-1;} /nb. of unassigned variables: 0/{ sub("[[]","",$0); opt=$3;} END{printf("%d",opt);}' > toulbar2_verif
    ub1b=`awk 'BEGIN{opt=-1;} {opt=$1} END{printf("%d", opt)}' toulbar2_verif`
   	if [[ $ub1 -lt 0 ]] ; then
      echo "error found"
      mv problem.wcsp errorflow$nerr.wcsp
      nerr=`expr $nerr + 1`
    fi
   	if [[ $ub1b -lt 0 ]] ; then
      echo "error found"
      mv problem.wcsp errorflow$nerr.wcsp
      nerr=`expr $nerr + 1`
    fi
   	if [[ $ub1 != $ub1b ]] ; then
      echo "error found"
      mv problem.wcsp errorflow$nerr.wcsp
      nerr=`expr $nerr + 1`
    fi
  
    # tests walldiff
    rm -f toulbar2_opt
    rm -f toulbar2_verif
    rm -f sol
    randomfile="walldiff-$n-$d-$tight-$bctr-$tctr-$nary-$seed"
    ./toulbar2 -random=$randomfile -C=$K -nopre -k=0 -z > /dev/null
    ./toulbar2 problem.wcsp $1 -w | awk 'BEGIN{opt=-1;} /^Optimum: /{opt=$2;} END{printf("%d",opt);}' > toulbar2_opt
    ub2=`awk 'BEGIN{opt=-1;} {opt=$1} END{printf("%d", opt)}' toulbar2_opt`
    ./toulbar2 problem.wcsp -x | awk 'BEGIN{opt=-1;} /nb. of unassigned variables: 0/{ sub("[[]","",$0); opt=$3;} END{printf("%d",opt);}' > toulbar2_verif
    ub2b=`awk 'BEGIN{opt=-1;} {opt=$1} END{printf("%d", opt)}' toulbar2_verif`
   	if [[ $ub2 -lt 0 ]] ; then
      echo "error found"
      mv problem.wcsp errornetwork$nerr.wcsp
      nerr=`expr $nerr + 1`
    fi
   	if [[ $ub2b -lt 0 ]] ; then
      echo "error found"
      mv problem.wcsp errornetwork$nerr.wcsp
      nerr=`expr $nerr + 1`
    fi
   	if [[ $ub2 != $ub2b ]] ; then
      echo "error found"
      mv problem.wcsp errornetwork$nerr.wcsp
      nerr=`expr $nerr + 1`
    fi

   	if [[ $ub1 != $ub2 ]] ; then
      echo "error found between flow and network"
      mv problem.wcsp error$nerr.wcsp
      nerr=`expr $nerr + 1`
    fi

    seed=`expr $seed + 1`
  done	
  nary=`expr $nary + 5` 
  tctr=`expr $tctr + 1`  
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

#!/bin/bash

# Usage:
# ./runall.sh ../validation

solver=./toulbar2
timelimit=3600
vmemlimit=16000000
K=1

rm -f outall

for e in `find $1 -regex ".*[.]wcsp" -print | sort` ; do
    dir=`dirname $e`
    base=`basename $e .wcsp`
    file=$dir/$base

    ubfile=${dir}/${base}.ub

    rm -f out
    rm -f outsolver
    rm -f usedtime

    if [[ -e $ubfile ]] ; then
	    ub=`cat $ubfile`
	    ub=`expr $ub \* $K`
        echo -n $ub > out
    else
        ub=
        echo -n "-" > out
    fi
   
    echo -n $file " "
    echo -n $file " " >> outall

    ulimit -t $timelimit > /dev/null
    ulimit -v $vmemlimit > /dev/null

#    (/usr/bin/time -f "%U user %S sys" $solver $file.wcsp -ub=$ub "${@:2}" -C=$K >> outsolver) 2> usedtime
#    cat outsolver | awk -v UB=$ub -f ./misc/script/runall.awk >> out ; cat out
# UNCOMMENT PREVIOUS *OR* NEXT TWO LINES
    (/usr/bin/time -f "%U user %S sys" $solver $file.wcsp "${@:2}" -C=$K >> outsolver) 2> usedtime
    cat outsolver | awk -v UB=-1 -f ./misc/script/runall.awk >> out ; cat out

    cat usedtime | awk '/ user /{ printf("%.2f",0.0+$1+$3); }'
   
    if [[ -e $ubfile ]] ; then
        cat out  | awk '/ /{ if($1 != $2) printf("           *******ERROR optimal cost"); }'
    fi
    
    cat out >> outall  
    cat usedtime | awk '/ user /{ printf(" %.2f",0.0+$1+$3); }' >> outall
    echo " " >> outall
    echo " " ;
done

cat outall | awk -f ./misc/script/evalresults.awk | sort

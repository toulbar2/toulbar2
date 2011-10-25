#!/bin/bash

# Usage:
# ./runall.sh ../validation

solver=./toulbar2
timelimit=600
K=1

rm -f outall

for e in `find $1 -regex ".*[.]wcsp" -print | sort` ; do
    dir=`dirname $e`
    base=`basename $e .wcsp`
    file=$dir/$base

# COMMENT OUT IF USING INITIAL UPPER BOUND
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
    (/usr/bin/time -f "%U user %S sys" $solver $file.wcsp -ub=$ub $2 -C=$K >> outsolver) 2> usedtime

    cat outsolver | awk -v UB=$ub -f runall.awk >> out ; cat out

    cat usedtime | awk '/ user /{ printf("%.2f",0.0+$1+$3); }'
   
    cat out  | awk '/ /{ if($1 != $2) printf("           *******ERROR optimal cost"); }'
    
    cat out >> outall  
    cat usedtime | awk '/ user /{ printf(" %.2f",0.0+$1+$3); }' >> outall
    echo " " >> outall
    echo " " ;
done

cat outall | awk -f ./misc/script/evalresults.awk | sort

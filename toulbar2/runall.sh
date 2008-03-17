#!/bin/sh

timelimit=1800
K=10


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
    fi
    
    echo -n $file " "
    echo -n $file " " >> outall
   
    ulimit -t $timelimit > /dev/null
    (/usr/bin/time -f "%U user %S sys" ./toulbar2 $file.wcsp $2C$K >> outsolver) 2> usedtime

    cat outsolver | awk 'BEGIN{opt="-";nodes=0;} /^Initial upperbound: /{if(ubini<0) ubini=$3;} /^Optimum: /{opt=$2; nodes=$7;}  /^No solution /{opt=ubini; nodes=$7;}  END{printf(" %s %d ",opt,nodes); }' >> out ; cat out

    cat usedtime | awk '/ user /{ printf("%.2f",0.0+$1+$3); }'
   

    cat out  | awk '/ /{ if($1 != $2) printf("           *******ERROR optimal cost"); }'
    
    cat out >> outall  
    cat usedtime | awk '/ user /{ printf(" %.2f",0.0+$1+$3); }' >> outall
    echo " " >> outall
    echo " "
done


cat outall | awk -f evalresults.awk | sort
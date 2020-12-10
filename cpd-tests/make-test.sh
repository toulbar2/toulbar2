#!/bin/bash

i=1
n=3
d=6
t=2

python3 ../bilevel-generator.py -n $n -d $d -t $2

zcat bilevel.cfn.gz | head -1 | cut -d ' ' -f8 > bilevel.opt

python3 ../make-bilevel-seqvars.py --pos p1.cfn.gz --neg p2.cfn.gz
shift=`zcat bilevel.cfn.gz | head -1`
shift2=`echo $shift | cut -d' ' -f7 `
shiftneg2=`echo $shift | cut -d' ' -f11 `

../../toulbar2 bilevel.cfn.gz bilevel.cov -B=1 -bilevel -hbfs: -nopre -shift2=$shift2 -shiftneg2=$shiftneg2 > bilevel.tb2

cat bilevel.opt
grep "Optimum" bilevel.tb2

#!/bin/bash

solver=./toulbar2
instance=$1
outbase=${1%.*}
shift
params="$@"
options="${params// /_}"

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo Instance is ${instance}
echo Options are ${options}
echo -n load:
uptime
taskset -c -p $$

/usr/bin/time --verbose --output=${outbase}_${options}.time ${solver} ${instance} "$@" | ./misc/script/bounds.sh ${outbase}_${options} > ${outbase}_${options}.out


#! /bin/tcsh

#usage: ./misc/script/wcsp2dzn.sh problem.wcsp MaximumForbiddenCost CostDivide
awk -f ./misc/script/wcsp2dzn.awk $1 1000000000 1 > ${1:r}.dzn
cp ./misc/script/wcsp.mzn ${1:h}/


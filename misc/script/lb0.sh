#!/bin/tcsh

foreach f ( $* )
echo $f:r
#awk 'BEGIN{lb=0} NF==2{lb=$2; exit} END{print lb}' ${f:r}.lb > ${f:r}.lb0
awk 'BEGIN{best=10**20} NF==2{if ($2<best) best=$2} END{print best}' ${f:r}_*.lb >! ${f:r}.lb0
end


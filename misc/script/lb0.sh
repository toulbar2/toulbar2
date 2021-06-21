#!/bin/tcsh

foreach f ( $* )
echo $f:r
awk 'BEGIN{lb=0} NF==2{lb=$2; exit} END{print lb}' ${f:r}.lb > ${f:r}.lb0
end


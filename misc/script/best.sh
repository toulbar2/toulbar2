#!/bin/tcsh

foreach f ( $* )
awk 'BEGIN{best=10**20} NF==2{if ($2<best) best=$2} END{print best}' ${f:r}_*.ub >! ${f:r}.best
echo $f:r
end

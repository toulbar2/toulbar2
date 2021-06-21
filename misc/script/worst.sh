#!/bin/tcsh

foreach f ( $* )
awk 'BEGIN{worst=0} NF==2{if ($2>worst) worst=$2} END{print worst}' ${f:r}_*.ub >! ${f:r}.worst
echo $f:r

end

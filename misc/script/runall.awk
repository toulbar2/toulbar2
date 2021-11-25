BEGIN{opt="-";initlb=0;bt=0;nodes=0;ub=UB;maxdepth=1;totaldepth=1} 

/^Optimum: /{opt=$2; bt=$4; nodes=$7;} 

/^No solution /{opt=ub; bt=$4; nodes=$7;} 

/Initial lower and upper bounds: /{gsub("[[]","",$0); gsub("[]]","",$0); gsub("[,]"," ",$0); initlb=$6; ub=$7;} 

/Depth:/{maxdepth=0+$2;totaldepth=0+$4}

END{printf(" %s %d %d %d %d %d ",opt,initlb,bt,nodes,maxdepth,totaldepth);}


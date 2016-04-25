BEGIN{opt="-";initlb=0;bt=0;nodes=0;} 

/^Optimum: /{opt=$2; bt=$4; nodes=$7;} 

/^No solution /{opt=UB; bt=$4; nodes=$7;} 

/Initial lower and upper bounds: /{gsub("[[]","",$0); gsub("[,]"," ",$0); initlb=$6;} 

END{printf(" %s %d %d %d ",opt,initlb,bt,nodes);}

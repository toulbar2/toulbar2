BEGIN{opt="-";nodes=0;initlb=0;tern=0;einit=0;elast=0;} 

/^Optimum: /{opt=$2; nodes=$7;} 

/^No solution /{opt=UB; nodes=$7;} 

/ ternary ctrs/{tern=$2} 

/ with maximum arity /{rmax=$16;einit=$10} 

/ current domains /{elast = $11} 

/Initial lower and upper bounds/{gsub("[[]","",$0); gsub("[,]"," ",$0); initlb=$6} 

END{printf(" %s %d %d %d %d %d ",opt,nodes,initlb,tern,einit,elast); }
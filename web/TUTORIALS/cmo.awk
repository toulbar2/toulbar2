FNR==1{n=0+$1;print "nodes=" n+1;printf "edges=["} NF==2{if (e) printf ","; e++;printf "(" 1+$1 "," 1+$2 ")"} END{print "]"; print "nbedges=" e}

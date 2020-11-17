
# Translate toulbar2 solution in CFN format (-s=3 or -w=3) to block modeling with M matrix and k clusters

# python3 blockmodel.py simple.mat 3 | toulbar2 --stdin=cfn -s=3 | awk -f ./sol2block.awk

BEGIN {
 k = 0;
}

{print $0}

/([A-Z][_0-9]+=[v]?[0-9]+)+/ {
 gsub("=v","=",$0);
 sub("^ +","",$0);
 s = $0;
 gsub("[~M_]","",s);
 k2 = length($0) - length(s);
 if (index($0, "M=")>0) k2--;
 k2 /=3;
 gsub("_"," ",$0);
 gsub("="," ",$0);
 for (i=2;i<=k2*4;i+=4) {
   M[0+$i, 0+$(i+1)] = 0+$(i+2);
   if (1+$i > k) k = 1+$i;
 }
 print "M" k "*" k " = ";
 for (i=0;i<k;i++) {
   for (j=0;j<k;j++) {
     if ((i,j) in M) {
       printf " " M[i,j];
     } else {
       printf " " M[j,i];
     }
   }
   print ""
 }
 for (c=0;c<k;c++) {
   printf " {"
   for (i=k2*4+1;i<=NF;i+=2) {
     if ($(i+1) == c) printf " " $i;
   }
   printf " }"
 }
 print "";
}


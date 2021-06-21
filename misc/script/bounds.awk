BEGIN {
  print "" > N ".ub";
  print "" > N ".lb";
}

FNR==1{start=$1}
{$1=$1-start}

/Initial lower and upper bounds:/ {
  gsub("[[]","",$0);
  gsub(","," ",$0);
#  print $1,0+$(NF-1) >> N ".ub";
  print $1,0+$(NF-2) >> N ".lb";
  fflush(N ".ub");
  fflush(N ".lb");
}

/New solution:/ && /backtracks/ {
  if (match($0,"log10like") || match($0,"energy")) {
    print $1,0+$(NF-12) >> N ".ub";
  } else {
    print $1,0+$(NF-8) >> N ".ub";
  }
  fflush(N ".ub");
  fflush(N ".lb");
} 

/New solution:/ && / in .* seconds[.]$/ {
  if (match($0,"log10like") || match($0,"energy")) {
    print $1,0+$(NF-7) >> N ".ub";
  } else {
    print $1,0+$(NF-3) >> N ".ub";
  }
  fflush(N ".ub");
  fflush(N ".lb");
} 

/Optimality gap:/{
  gsub("[[]","",$0);
  gsub(","," ",$0);
  print $1,0+$(NF-9) >> N ".lb";
  fflush(N ".ub");
  fflush(N ".lb");
} 

/Optimum:/{
 if (match($0,"log10like") || match($0,"energy")) {
  if (match($0,"DEE")) {
    print $1,0+$(NF-18) >> N ".ub";
    print $1,0+$(NF-18) >> N ".lb";
  } else {
    print $1,0+$(NF-13) >> N ".ub";
    print $1,0+$(NF-13) >> N ".lb";
  }
 } else {
  if (match($0,"DEE")) {
    print $1,0+$(NF-14) >> N ".ub";
    print $1,0+$(NF-14) >> N ".lb";
  } else {
    print $1,0+$(NF-9) >> N ".ub";
    print $1,0+$(NF-9) >> N ".lb";
  }
 }
  fflush(N ".ub");
  fflush(N ".lb");
}

/[0-9]+ .[0-9]+,[0-9]+.[/][0-9]+[/][0-9]+[/][0-9]+ [0-9.]+% [0-9.]+/ {
  gsub("[[]","",$0);
  gsub(","," ",$0);
  gsub("/"," ",$0);
  print $1,0+$3 >> N ".lb";
  fflush(N ".ub");
  fflush(N ".lb");
}

{print $0;}


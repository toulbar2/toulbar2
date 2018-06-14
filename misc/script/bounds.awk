BEGIN {
  print "" > "ub" N;
  print "" > "lb" N;
  print "" > "open" N;
}

FNR==1{start=$1}
{$1=$1-start}

/Initial lower and upper bounds:/ {
  gsub("[[]","",$0);
  gsub(","," ",$0);
#  print $1,$(NF-1) >> "ub" N;
  print $1,$(NF-2) >> "lb" N;
  fflush("ub" N);
  fflush("lb" N);
}

/New solution:/ {
  if (match($0,"log10like") || match($0,"energy")) {
    print $1,$(NF-10) >> "ub" N;
  } else {
    print $1,$(NF-6) >> "ub" N;
  }
  fflush("ub" N);
  fflush("lb" N);
} 

/Optimality gap:/{
  print $1,$(NF-9) >> "lb" N;
  fflush("ub" N);
  fflush("lb" N);
} 

/Optimum:/{
 if (match($0,"log10like") || match($0,"energy")) {
  if (match($0,"DEE")) {
    print $1,$(NF-18) >> "ub" N;
    print $1,$(NF-18) >> "lb" N;
  } else {
    print $1,$(NF-13) >> "ub" N;
    print $1,$(NF-13) >> "lb" N;
  }
 } else {
  if (match($0,"DEE")) {
    print $1,$(NF-14) >> "ub" N;
    print $1,$(NF-14) >> "lb" N;
  } else {
    print $1,$(NF-9) >> "ub" N;
    print $1,$(NF-9) >> "lb" N;
  }
 }
  fflush("ub" N);
  fflush("lb" N);
}

/[0-9]+ .[0-9]+,[0-9]+.[/][0-9]+[/][0-9]+[/][0-9]+ [0-9.]+% [0-9.]+/ {
  gsub("[[]","",$0);
  gsub(","," ",$0);
  gsub("/"," ",$0);
  print $1,$3 >> "lb" N;
  print $1,$5 >> "open" N;
  fflush("ub" N);
  fflush("lb" N);
  fflush("open" N);
}

{print $0;}

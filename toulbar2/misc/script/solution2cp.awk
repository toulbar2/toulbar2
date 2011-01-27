
# Translate the solution found by toolbar in cp format

# Usage: toolbar -v problem.wcsp | awk -f solution2cp.awk problem.cp -

# Note: it can save the solution in a file whose name is problem.sol (see awk variable SAVESOLUTION)
# SAVESOLUTION = 1 (default value) :  save the solution in a file
# SAVESOLUTION = 0                 :  do not save the solution

# Verification of the solution cost: 
# cat problem.sol problem.cp | awk -f cp2wcsp.awk | ../../toolbar/toolbar
# or a safer (solver independent) verification:
# cat problem.sol problem.cp | awk -f cp2wcsp.awk | awk '/^0 / && NF == 3{cost+=$2} END{print "--- solution cost = " cost}'

BEGIN {
  SAVESOLUTION = 1;
  first = 1;
  tuplemode = 0;
  idx = 0;
}

# outside a constraint defined by a list of tuples
FNR == NR && !first && tuplemode && !/^\#/ && NF >= 1 && !/^ *((-)?[0-9]+ +)*((-)?[0-9]+) *$/ {
    tuplemode = 0;
}

# a new constraint defined by a list of tuples
FNR == NR && !first && !tuplemode && !/^\#/ && NF >= 2 && /^ *([a-zA-Z_][a-zA-Z0-9_]* +)+((-)?[0-9]+) *$/ && (NF >= 3 || ($1 in domainsize)) {
    tuplemode = 1;
}

# a new variable
FNR == NR && !first && !tuplemode && !/^\#/ && NF >= 1 && /^ *([a-zA-Z_][a-zA-Z0-9_]*)( +(-)?[0-9]+)+ *$/ {
  if (NF - 1 >= 2) { # constants are ignored
    var = $1;
    domainsize[var] = NF - 1;
    varname[idx] = var;
    for (i = 2; i <= NF; i++) {
      domains[idx, i - 2] = $i;
    }
    idx++;
  }
}
  
# problem name and possibly global upper bound
FNR == NR && first && !/^\#/ && NF >= 1 {
  first = 0;
}

FNR != NR && SAVESOLUTION && /^Filename = / {
  if (match($3,".*[.]")) filename = substr($3,1,RLENGTH) "sol";
  else filename = $3 ".sol";
  print "--- output solution in file ",filename;
  print $3 > filename;
}

FNR != NR && SAVESOLUTION && /^Optimum: / {
  print "# optimum = ",$2 >> filename;
}

FNR != NR && SAVESOLUTION && /^Best bound: / {
  print "# solution cost = " $3 >> filename;
}

FNR != NR && $2 == "solution:" {
  for (i=3; i<=NF; i++) {
    $i = domains[i-3, $i];
    if (SAVESOLUTION) {
#      print "hard(",varname[i-3],"==",$i,")" >> filename;
      print varname[i-3],$i >> filename;
    }
  }
}

FNR != NR {
  print $0;
}

END {
  if (SAVESOLUTION) {
    printf("#") >> filename;
    close(filename);
  }
}

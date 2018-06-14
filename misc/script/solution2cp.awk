
# Translate the solution found by toulbar2 in cp format

# Usage: toulbar2 problem.wcsp -s | awk -f solution2cp.awk problem.cp -

# Note: it can save the solution in a file whose name is problem.sol (see awk variable SAVESOLUTION)
# SAVESOLUTION = 1 (default value) :  save the solution in a file
# SAVESOLUTION = 0                 :  do not save the solution

# Note bis: it can read a problem solution instead of the solver output
# SOLVEROUTPUT = 1 (default value) :  read solutions from solver output
# SOLVEROUTPUT = 0                 :  read from problem solution file

# Verification of the solution cost: 
#  gawk -f cp2wcsp.awk problem.sol problem.cp > problem.wcsp
#  toulbar2 problem.wcsp
# or using a safer (solver independent) verification:
#  gawk -f cp2wcsp.awk problem.sol problem.cp | awk '/^0 / && NF == 3{cost+=$2} END{print "# solution cost = " cost}'
# (warning! these verifications work only if all the cost functions are defined by a formula, see e.g. golomb4.c
#  and an initial upper bound is given in problem.cp)

# Example:

# awk -f cp2wcsp.awk celar6sub0.cp > ! celar6sub0.wcsp
# toulbar2 celar6sub0.wcsp -s | awk -f ./solution2cp.awk celar6sub0.cp -
# gawk -f cp2wcsp.awk celar6sub0.sol celar6sub0.cp | awk '/^0 / && NF == 3{cost+=$2} END{print "# solution cost = " cost}'

BEGIN {
  SAVESOLUTION = 1;
  SOLVEROUTPUT = 1;
  first = 1;
  tuplemode = 0;
  idx = 0;
  oksol = 0;
}

# outside a constraint defined by a list of tuples
FNR == NR && !first && tuplemode && !/^ *\#/ && NF >= 1 && !/^ *((-)?[0-9]+ +)*((-)?[0-9]+) *$/ {
    tuplemode = 0;
}

# a new constraint defined by a list of tuples
FNR == NR && !first && !tuplemode && !/^ *\#/ && NF >= 2 && /^ *([a-zA-Z_][a-zA-Z0-9_]* +)+((-)?[0-9]+) *$/ && (NF >= 3 || ($1 in domainsize)) {
    tuplemode = 1;
}

# a new variable
FNR == NR && !first && !tuplemode && !/^ *\#/ && NF >= 1 && /^ *([a-zA-Z_][a-zA-Z0-9_]*)( +(-)?[0-9]+)+ *$/ {
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
FNR == NR && first && !/^ *\#/ && NF >= 1 {
  first = 0;
}

!SOLVEROUTPUT && FNR != NR && FNR==1 {
	filename = FILENAME ".cp";
	print "" > filename;
}

SOLVEROUTPUT && FNR != NR && SAVESOLUTION && /^loading wcsp file:/ {
  if (match($4,".*[.]")) filename = substr($4,1,RLENGTH) "sol";
  else filename = $4 ".sol";
  print "--- output solution in file ",filename;
  print $4 > filename;
}

SOLVEROUTPUT && FNR != NR && SAVESOLUTION && /^Optimum: / {
  print "# optimum = ",$2 >> filename;
}

FNR != NR && (oksol || !SOLVEROUTPUT || /^[0-9]+ solution[:] /) {
  for (i=($2=="solution:")?3:1; i<=NF; i++) {
	$i = domains[i-(($2=="solution:")?3:1), $i];
    if (SAVESOLUTION) {
#      print "hard(",varname[i-(($2=="solution:")?3:1)],"==",$i,")" >> filename;
      print varname[i-(($2=="solution:")?3:1)],$i >> filename;
    }
  }
  oksol = 0;
}

SOLVEROUTPUT && FNR != NR && /New solution:/ {
	if (SAVESOLUTION) print "#",$0 > filename;
	oksol = 1;
}

SOLVEROUTPUT && FNR != NR {
	print $0;
}

END {
  if (SAVESOLUTION) {
#    printf("#") >> filename;
    close(filename);
  }
}

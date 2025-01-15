
# Translate a Weighted Max-SAT problem in wcnf format into Minizinc

# convention: 1: true , 0: false

# See the file cnf_format.ps for a description of the cnf format

# Warning! Each clause must be defined on a single line with a zero at the end

# Preprocessing:
# - Tautologies ("x and not(x)") are removed
# - Duplicate literals are simplified ("x and x" replaced by "x")

# Remark: some variables may be not involved in any clauses.

# Usage: 
# awk -f wcnf2dzn.awk file.wcnf > file.dzn
# minizinc wcnf.mzn file.dzn

BEGIN {
    ok = 0;
    costs = "";
    top = 100000000;
    first = 1;
}

ok && !/^c/ && NF >= 2 && $NF == 0 {
    delete vars;
    removed = 0;
    scope = "";
    arity = 0;
    weight=$1;
    for (i=2; i<NF; i++) {
		if ($i >= 1) {
			if ($i in vars) {
				if (vars[$i] == -1) {
					removed = 1; # remove a tautology
				} # else do not duplicate the variable
			} else {
				vars[$i] = 1;
				arity++;
				scope = scope "," $i;
			}
		} else {
			if ((- $i) in vars) {
				if (vars[(- $i)] == 1) {
					removed = 1; # remove a tautology
				} # else do not duplicate the variable
			} else {
				vars[(- $i)] = -1;
				arity++;
				scope = scope "," $i;
			}
		}
    }
    if (!removed && (arity > 0 || weight < top)) {
		nbclause++;
		if (first) {
			printf("{");
			first = 0;
		} else {
			printf(",{");
		}
		printf("%s",substr(scope,2));
		if (weight < top) {
			nbcost++;
			costs = costs "," weight;
			printf(",%d",nbvar+nbcost);
		}
		printf("}");
    }
}

/^p/ {
	nbcost = 0;
    nbvar = $3;
    nbclause = 0;
	if ($5 < top) top = $5;
    ok = 1;
	printf("formula = [");
}

END {
	print "];"
	print "costs = [" substr(costs,2) "];"
	print "nb_variables = " nbvar ";"
	print "nb_clauses = " nbclause ";"
	print "nb_costs = " nbcost ";"
	print "max_cost = " top ";"
}

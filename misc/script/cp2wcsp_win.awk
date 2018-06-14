
# Translator from cp format to wcsp format

# gawk program to transform WCSP problems with cost functions defined by formulae 
# into cost functions defined by tuples

# Usage: gawk -f cp2wcsp_win.awk problem.cp > problem.wcsp

# Note that negative costs are replaced by the sum of the maximum costs per cost function plus one
# or by a user-defined global upper bound (second field after the problem name)

# For faster evaluations of the cost function formulae, 
# Use shared cost functions to save memory and reduce formulae evaluation time

# Warning! all variables in wcsp format will have an enumerated domain (no interval) 

# Example:

# gawk -f ./cp2wcsp_win.awk golomb4.cp > golomb4.wcsp
# toulbar2 golomb4.wcsp -s | gawk -f ./solution2cp.awk golomb4.cp -

function perror(i,message) {
	print message;
	print "line " NR ": " $0;
	error = i;
	exit(i);
}

function min(x,y) {
    if (x < y) return(x); 
    else return(y);
}

function max(x,y) {
    if (x > y) return(x); 
    else return(y);
}

# generate all the possible tuples with their associated cost value
# i is a local variable
function enumerate(assignment,tuple,depth,i) {
    if (depth>arity[nbc]) {
		print assignment "print (" formula ");" | command " > " tempfile;
		nbt++;
		tuples[nbc,nbt] = tuple;
    } else if (domainsize[orderedvars[nbc,depth]] == 1) {
		enumerate(assignment "" orderedvars[nbc,depth] "=" domains[orderedvars[nbc,depth],1] ";", tuple, depth + 1);
    } else {
		for (i=1;i<=domainsize[orderedvars[nbc,depth]];i++) {
			enumerate(assignment "" orderedvars[nbc, depth] "=" domains[orderedvars[nbc, depth],i] ";", tuple "" i-1 " ", depth + 1);
		}
    }
}

# intialization
BEGIN {
	MAXCOST = 10000000;
    error = 0;
    first = 1;
    nbc = 0;
    nbv = 0;
    nbt = 0;
    maxdom = 0;
    name = FILENAME;
    ub = -1;
    tuplemode = 0;
	tempfile = ("tmpawk." PROCINFO["pid"]);
    command="gawk --assign ub=" ub " -f libcp.awk -f -";
    infinity = 0;
    maxcost = 0;
    nbshared = 0;
}

# tuple outside a cost function!
!first && !tuplemode && !/^ *\#/ && NF >= 1 && /^ *((-)?[0-9]+ +)+((-)?[0-9]+) *$/ {
	perror(1,"Error: tuple defined outside a cost function!");
}

# bad arity for the current tuple!
!first && tuplemode && !/^ *\#/ && NF >= 1 && /^ *((-)?[0-9]+ +)*((-)?[0-9]+) *$/ && NF - 1 != arity[nbc] {
	perror(2,"Error: tuple with a bad arity! Should be equal to " arity[nbc] " + 1 (cost)");
}

# inside a cost function defined by a list of tuples
!first && tuplemode && !/^ *\#/ && NF >= 1 && /^ *((-)?[0-9]+ +)*((-)?[0-9]+) *$/ && NF - 1 == arity[nbc] && $NF != defaultvalue[nbc] {
    nbt++;
    tuplesize[nbc]++;
    realtuplesize[nbc]++;
    tuple = "";
    for (i=1; i <= NF-1 ;i++) {
        if (!((orderedvars[nbc,i] SUBSEP $i) in domainsindex)) {
			nbt--;
			tuplesize[nbc]--;
			realtuplesize[nbc]--;
			next; # to be compatible with solution2cp.awk
			perror(3,"Error: value " $i " not defined in the domain of " orderedvars[nbc,i] "!");
        }
		if (domainsize[orderedvars[nbc,i]] > 1) {
			tuple = tuple "" domainsindex[orderedvars[nbc,i], $i] - 1 " ";
		}
    }
    tuples[nbc,nbt] = tuple;
    costs[nbc,nbt] = $NF;
    maxcost = max(maxcost,$NF);
}

# outside a cost function defined by a list of tuples
!first && tuplemode && !/^ *\#/ && NF >= 1 && !/^ *((-)?[0-9]+ +)*((-)?[0-9]+) *$/ {
    tuplemode = 0;
    if (!(nbc in definedby) && realarity[nbc] == 0) {
		if (tuplesize[nbc] > 1) {
			perror(9, "Error: zero arity cost function with more than one tuple!");
		}
		if (tuplesize[nbc] == 1) {
			defaultvalue[nbc] = costs[nbc,1];
			realtuplesize[nbc] = 0;
		}
    }
    if (nbc in definedby) {
		if (!(definedby[nbc] in sharedmaxcost)) {
			perror(10,"Error: shared n-ary cost function undefined!");
		}
		infinity += sharedmaxcost[definedby[nbc]];
    } else {
		infinity += maxcost;
    }
}

# a new cost function defined by a list of tuples
!first && !tuplemode && !/^ *\#/ && NF >= 2 && /^ *([a-zA-Z_][a-zA-Z0-9_]* +)+((-)?[0-9]+) *$/ && (NF >= 3 || ($1 in domainsize)) {
    nbc++;
    tuplemode = 1;

    if (match($0, " +defined +by +")) {
		if ($NF > nbshared) perror(11,"shared n-ary cost function " $NF " undefined!");
		definedby[nbc] = $NF;
		sub(" +defined +by +"," ");
    }

    arity[nbc] = NF - 1;
    realarity[nbc] = 0;
    
    for (i=1; i <= NF-1 ;i++) {
		var = $i;
		if (!(var in domainsize)) {
			perror(4,"Error: variable " var " not defined before used in the following cost function!");
		}
		if (domainsize[var] > 1) realarity[nbc]++;
		orderedvars[nbc,i] = var;
    }
    nbt = 0;
    defaultvalue[nbc] = $NF;
    tuplesize[nbc] = 0;
    realtuplesize[nbc] = 0;
    maxcost = max(0,defaultvalue[nbc]);
}

# a new variable
!first && !tuplemode && !/^ *\#/ && NF >= 1 && /^ *([a-zA-Z_][a-zA-Z0-9_]*)( +(-)?[0-9]+)+ *$/ {
    var = $1;
    if (var in domainsize) {
		next; # to be compatible with solution2cp.awk
#      perror(5, "Error: variable " var " already defined!");
    }
    domainsize[var] = NF - 1;
    if (domainsize[var] > maxdom) maxdom = domainsize[var];
    if (domainsize[var] > 1) {
		nbv++;
		variables[nbv] = var;
		varnum[var] = nbv;
    }
    for (i = 2; i <= NF; i++) {
		domains[var,"" i - 1] = $i;
		domainsindex[var, $i] = i - 1;
    }
}

# a new cost function (at least on two variables) defined by a specific wcsp cost function keyword
!first && !tuplemode && !/^ *\#/ && NF >= 1 && !/^ *([a-zA-Z_][a-zA-Z0-9_]*)( +(-)?[0-9]+)+ *$/ && /^ *([a-zA-Z_][a-zA-Z0-9_]*)( +[a-zA-Z_][a-zA-Z0-9_]*)+ +(-1) +([a-zA-Z<>=][a-zA-Z0-9_<>=]*)/ {
	nbc++;
	copy = $0;
	_arity = 0;
	for (i=1;i<NF;i++) {
		if ($i == "-1") break; 
		else _arity++;
	}
    arity[nbc] = _arity;
    realarity[nbc] = _arity;
    
    for (i=1; i<=_arity; i++) {
		var = $i;
		if (!(var in domainsize)) {
			perror(4,"Error: variable " var " not defined before used in the following cost function!");
		}
		if (domainsize[var] == 1) {
			perror(12,"Error: variable " var " has an initial domain size equal to 1! It is not allowed by wcsp cost functions in intention (" copy "). You should declare an additional unary constraint to express this unit domain.");
		}
		if (i>=2) {
			for (j=1; j<=domainsize[var]; j++) {
				if (domains[var,j] != domains[orderedvars[nbc,i-1],j]) {
					perror(13,"Error: variable " var " has an initial domain different from another variable in the scope of a wcsp cost function in intention (" copy "). You should declare all the variables in the scope by using the same domain and add additional unary constraints to take into account specific restrictions on some domains.");
				}
				if ((match($(_arity+2),"[<>=]") || match($(_arity+2),"disj")) && j>=2 && domains[var,j] != domains[var,j-1]+1) {
					perror(14,"Error: variable " var " has an initial domain with non-consecutive or non-increasing integer values. It is not allowed by arithmetic and disjunctive wcsp cost functions in intention (" copy "). You should declare all the variables in the scope by using the same domain of increasing consecutive integer values and add additional unary constraints to take into account specific restrictions on some domains.");
				}
			}
		}			
		orderedvars[nbc,i] = var;
		$i = varnum[var] -1;
    }
	if (match($(_arity+2),"[<>=]")) infinity += $(_arity+4);
	if ($(_arity+2) == "disj") {
		if ($(_arity+5) < 0) $(_arity+5) = ((ub>0)?ub:MAXCOST);
		else infinity += $(_arity+5);
	}
	if ($(_arity+2) == "sdisj") {
		$(_arity+5) -= domains[orderedvars[nbc,1],"" 1];
		$(_arity+6) -= domains[orderedvars[nbc,2],"" 1];
		if ($(_arity+7) < 0) $(_arity+7) = ((ub>0)?ub:MAXCOST);
		else infinity += $(_arity+7);
		if ($(_arity+8) < 0) $(_arity+8) = ((ub>0)?ub:MAXCOST);
		else infinity += $(_arity+8);
	}
	if ($(_arity+2) == "salldiff") {
		if ($(_arity+4) < 0) $(_arity+4) = ((ub>0)?ub:MAXCOST);
		else {
			if ($(_arity+3) == "var") infinity += $(_arity+4) * _arity;
			else infinity += $(_arity+4) * _arity * (_arity - 1) / 2;
		}
	}
	if ($(_arity+2) == "sgcc") {
		for (j=_arity+6; j<=NF; j +=3) {
			$j -= domains[orderedvars[nbc,1],"" 1];
		}
		if ($(_arity+4) < 0) $(_arity+4) = ((ub>0)?ub:MAXCOST);
		else {
			if ($(_arity+3) == "var") infinity += $(_arity+4) * _arity;
			else infinity += $(_arity+4) * _arity * domainsize[orderedvars[nbc,1]];
		}
	}
	if ($(_arity+2) == "ssame") {
		for (j=_arity+6; j<=NF; j++) {
			if (!($j in domainsize)) {
				perror(4,"Error: variable " $j " not defined before used in the following cost function!");
			}

			$j = varnum[$j] - 1;
		}
		if ($(_arity+3) < 0) $(_arity+3) = ((ub>0)?ub:MAXCOST);
		else infinity += $(_arity+3) * $(_arity+4);
	}
	if ($(_arity+2) == "sregular") {
		for (j=_arity+10+($(_arity+6))+($(_arity+7+$(_arity+6))); j<=NF; j+=3) {
			$j -= domains[orderedvars[nbc,1],"" 1];
		}
		if ($(_arity+4) < 0) $(_arity+4) = ((ub>0)?ub:MAXCOST);
		else {
			infinity += $(_arity+4) * _arity;
		}
	}

	wcspfunctions[nbc] = _arity " " $0;
	next;
}	

# a new cost function defined by a formula
!first && !tuplemode && !/^ *\#/ && NF >= 1 && !/^ *([a-zA-Z_][a-zA-Z0-9_]*)( +(-)?[0-9]+)+ *$/ {
    nbc++;
    arity[nbc] = 0;
    realarity[nbc] = 0;

    if (match($0, "^ *shared[(].*[)] *$")) {
		nbshared++;
		shared[nbc] = 1;
		sub("^ *shared[(]","");
		sub("[)] *$","");	
    }

    delete vars;

    string = $0;
    while (start = match(string, "[a-zA-Z_][a-zA-Z0-9_]*")) {
		size = match(substr(string, start + 1), "[^a-zA-Z0-9_]");
		if (size == 0) size = length(string) + 1;
		var = substr(string, start, size)
		if (!(var in vars) && substr(string, start + size, 1) != "(") {
			arity[nbc]++;
			if (!(var in domainsize)) {
				perror(6, "Error: variable " var " not defined before used in the following cost function!");
			}
			if (domainsize[var] > 1) realarity[nbc]++;
			vars[var] = arity[nbc];
			orderedvars[nbc,arity[nbc]] = var;
		}
		string = substr(string, start + size);
    }
    nbt = 0;
    formula = $0;
    print "BEGIN {" | command " > " tempfile;
    enumerate("","",1);
    print "}" | command " > " tempfile;
    close(command " > " tempfile);
    delete freq;
    maxfreq = 0;
    defaultvalue[nbc] = 0;
    maxcost = 0;
    for (i=1; i<=nbt; i++) {
		getline costs[nbc,i] < tempfile;
		maxcost = max(maxcost,costs[nbc,i]);
		freq[costs[nbc,i]]++;
		if (freq[costs[nbc,i]] > maxfreq) {
			maxfreq = freq[costs[nbc,i]];
			defaultvalue[nbc] = costs[nbc,i];
		}
    }
    close(tempfile);
    tuplesize[nbc] = nbt;
    realtuplesize[nbc] = nbt - maxfreq;
    if (realarity[nbc] == 0 && tuplesize[nbc] != 1) {
		perror(8, "Error: zero arity cost function without just one tuple!");
    }
    if (nbc in shared) {
		sharedmaxcost[nbshared] = maxcost;
    }
    infinity += maxcost;
}

# problem name and possibly global upper bound
first && !/^ *\#/ && NF >= 1 {
    if (NF >= 3) {
		perror(7, "Error: the first line does not define the problem name!");
    }
    first = 0;
    name = $1;
    if (NF >= 2) {
        # global upper bound defined by the user
		ub = $2;
		command="gawk --assign ub=" ub " -f libcp.awk -f -";
    }
}

# output file result
END {
    if (!error) {
        if (!first && tuplemode) { # do not forget to close a cost function in extension
			if (!(nbc in definedby) && realarity[nbc] == 0) {
				if (tuplesize[nbc] > 1) {
					perror(10, "Error: zero arity cost function with more than one tuple!");
				}
				if (tuplesize[nbc] == 1) {
					defaultvalue[nbc] = costs[nbc,1];
					realtuplesize[nbc] = 0;
				}
			}
			if (nbc in definedby) {
				if (!(definedby[nbc] in sharedmaxcost)) {
					perror(10,"Error: shared n-ary cost function undefined!");
				}
				infinity += sharedmaxcost[definedby[nbc]];
			} else {
				infinity += maxcost;
			}
		}
		# compute the upper bound
		infinity++;
		if (ub < 0) ub = infinity;
		# print name and problem sizes
		print name,nbv,maxdom,nbc,ub;
		# print every domain size
		if (nbv > 0) {
			printf("%d", domainsize[variables[1]]);
			for (i=2; i <= nbv; i++) {
				printf(" %d", domainsize[variables[i]]);
			}
			print "";
		}
		# print all the cost functions with a list of tuples with their costs
		for (i=1; i <= nbc; i++) {
			if (i in wcspfunctions) {
				print wcspfunctions[i];
				continue;
			}
			if (i in shared) printf("-%d", realarity[i]);
			else printf("%d", realarity[i]);
			for (j=1; j <= arity[i]; j++) {
				if (domainsize[orderedvars[i,j]] > 1) printf(" %d", varnum[orderedvars[i,j]] - 1);
			}
			print " " ((defaultvalue[i] < 0)?ub:defaultvalue[i]),((i in definedby)?-definedby[i]:realtuplesize[i]);
			for (k=1; k <= tuplesize[i]; k++) {
				if (costs[i,k] != defaultvalue[i]) {
					print tuples[i,k] "" ((costs[i,k] < 0)?ub:costs[i,k]);
				}
			}
		}
    }
}


# Translate in dzn format from an original binary WCSP instance including linear constraints (knapsackp)

# Usage: awk -f wcsp2dzn.awk smallrandom.wcsp TOP DIV > smallrandom.dzn

# Maximum cost is given by TOP and all costs (including initial upper bound) are divided by DIV

# Warning! wcsp problem file must have only unary and binary cost functions in extension with all tuples defined

function cost2int(x)
{
	if (x<TOP) return(int(x/DIV));
	else return(int(TOP/DIV));
}

BEGIN {
	TOP = 0+ARGV[2];
	DIV = 0+ARGV[3];
	ARGV[2]="";
	ARGV[3]="";

	tuple = 0;
	binary = 0;
	n = 0;
	c1 = 0;
	cum1 = 0;
	c2 = 0;
	cum2 = 0;
	maxcost1 = 0;
	maxcost2 = 0;
	i = 0;
	j = 0;
	doms = "";
	costs1 = "";
	costs2 = "";
	nbtuples1 = "";
	nbtuples2 = "";
	cumtuples1 = "";
	cumtuples2 = "";
	func1x = "";
	func2x = "";
	func2y = "";
	e = 0;
	numkp = 0;
	maxkp = 0;
	totalw = 0;
	capacities = "";
	numw = "";
	cumw = "";
	kpvar = "";
	kpval = "";
	kpw = "";
}

# read wcsp header
FNR == 1 {
	n =$2;
	maxd = $3;
	c = $4;
	ub = $5;
	lb = 0;
}

# read domain sizes
FNR == 2 {
	for (i=1; i<=NF; i++) {
		domainsize[i] = $i;
		doms = doms "," $i;
	}
}

# read zero-artiy, unary and binary cost functions ONLY
FNR >= 3 {
  if (tuple==0) {
	  # read header of current unary or binary cost function
	  constraint = "";
	  if ($(3+$1) == "knapsackp") {
		  numkp++;
		  capacities = capacities "," 0+$(4+$1);
		  w = 0;
		  pos = 5+$1;
		  for (i=1; i<=0+$1; i++) {
		    k = 0+$pos;
		    pos++;
		    for (j=1; j<=k; j++) {
		      #print $(1+i), $(pos), $(pos+1);
		      kpvar = kpvar "," 1+$(1+i);
		      kpval = kpval "," $(pos);
		      kpw = kpw "," $(pos+1);
		      pos += 2;
		      w++;
		    }
		  }
		  numw = numw "," w;
		  cumw = cumw "," totalw;
		  totalw += w;
		  maxkp = totalw;
		  tuple = 0;
	  } else if ($1 == 2) {
		  c2++;
		  nbtuples2 = nbtuples2 "," $NF;
		  cumtuples2 = cumtuples2 "," cum2;
		  cum2 += $NF;
		  func2x = func2x "," 1+$2;
		  func2y = func2y "," 1+$3;
		  binary = 1;
		  tuple = $NF;
	  } else if ($1 == 1) {
		  c1++;
		  nbtuples1 = nbtuples1 "," $NF;
		  cumtuples1 = cumtuples1 "," cum1;
		  cum1 += $NF;
		  func1x = func1x "," 1+$2;
		  binary = 0;
		  tuple = $NF;
	  } else if ($1 == 0) {
		  lb += $(NF-1);
		  tuple = $NF;
	  } else {
		  print "ERROR bad arity" > "/dev/stderr";
		  err = 1;
		  exit(err);
	  }
          e++;
  } else {
	  if (binary) {
		  if ($NF > maxcost2) maxcost2=$NF;
		  constraint = constraint "," cost2int($NF) "," $1 "," $2;
		  tuple--;
		  if (tuple==0) {
			  costs2 = costs2 "," substr(constraint,2) "";
		  }
	  } else  {
		  if ($NF > maxcost1) maxcost1=$NF;
		  constraint = constraint "," cost2int($NF) "," $1;
		  tuple--;
		  if (tuple==0) {
			  costs1 = costs1 "," substr(constraint,2) "";
		  }
	  }
  }
}

END {
	print "num_variables= " n ";";
	print "max_domain = " maxd ";";
	print "domains = [" substr(doms,2) "];";
	print "top = " cost2int(ub) ";";
	print "cost0 = " cost2int(lb) ";";
	print "num_constraints1= " c1 ";";
	print "max_constraints1 = " cum1*2 ";";
	print "max_costs1 = " cost2int(maxcost1) ";";
	print "num_constraints2= " c2 ";";
	print "max_constraints2 = " cum2*3 ";";
	print "max_costs2 = " cost2int(maxcost2) ";";
	print "num_knapsack= " numkp ";";
	print "max_knapsack = " maxkp ";";
	print "func1x = [" substr(func1x,2) "];";
	print "num_tuples1 = [" substr(nbtuples1,2) "];";
	print "cum_tuples1 = [" substr(cumtuples1,2) "];";
	print "costs1 = [" substr(costs1,2) "];";
	print "func2x = [" substr(func2x,2) "];";
	print "func2y = [" substr(func2y,2) "];";
	print "num_tuples2 = [" substr(nbtuples2,2) "];";
	print "cum_tuples2 = [" substr(cumtuples2,2) "];";
	print "costs2 = [" substr(costs2,2) "];";
	print "knapsack_capacity = [" substr(capacities,2) "];";
	print "num_weights = [" substr(numw,2) "];";
	print "cum_weights = [" substr(cumw,2) "];";
	print "knapsack_variable = [" substr(kpvar,2) "];";
	print "knapsack_value = [" substr(kpval,2) "];";
	print "knapsack_weight = [" substr(kpw,2) "];";
}

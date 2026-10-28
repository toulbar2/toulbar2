
# translate flatzinc format to python/Numberjack 

BEGIN {
	MAX_COST_FZN = 2147483647 * 1024;
	MAXINT = 10000; # maximum integer domain value
	#VERBOSITY = -1;
	#SHOWSOLUTIONS = 0;
	#DUMPMODEL = 0;	
	# ALLSOLUTIONS = 10000; # needed for running unit tests (automatically assigned by ./fzn2py with extra option -a)
	#FREESEARCH = 1;
	TIMELIMIT = 3600;
	print "#! /usr/bin/env python"
	print "import time"
	print "import pytoulbar2 as tb2"
	system("cat " MZNNJ_DIR "/fzn2py.py");
	print "";
	print "def get_model():"
	print "  global model";
	print "  global DelayedLinEq";
	print "  model = tb2.CFN(ubinit=" int((MAX_COST_FZN / 2) / 3 / 3) ",verbose=" VERBOSITY ")";
	print "  try:";
	parameter = 1;
	error = 0;
}

# { print $0 } # for debugging purposes only

/::_output_var/ {
	gsub("::_output_var","");
	$0 = $0 " ::_output_";
}

/::_output_array/ {
	regexp = "::_output_array[(][[](range[(][0-9]+,1[+][0-9]+[)][ ,]*)+[]][)]";
	if (match($0,regexp)) {
	  array = substr($0, RSTART, RLENGTH);
	  sub(regexp,"");
	  gsub("[[]","",array);
	  gsub("[]]","",array);
	  gsub(" ","",array);
	  gsub("range[(]","",array);
	  gsub(",1[+]","..",array);
	  gsub("[)]","",array);
	  nbdim = array;
	  gsub("[^,]","",nbdim);
	  dim = length(nbdim)+1;
	  sub("::_output_array", "array" dim "d",array);
	  $0 = $0 " ::_output_" array;
	} else {
		print "ERROR WRONG OUTPUT ARRAY DEFINITION",$0;
		error = 3;
		exit(2);
	}
}

/^var / {
	parameter = 0;
}

/ of var / {
	parameter = 0;
}

parameter {
	sub(".*: ","");
	if ($2 == "=") print "    " $0; # it might be a predicate instead if no =
}

/^var bool:/{
	print "    " $3 " = Variable(0, 1, '" $3 "')";
	if (match($0,"::_output_")) output[$3] = 1;
	if ($4 == "=") print "    int_eq(" $3 " , " $5 ")";
}

/^var int:/{
	print "    " $3 " = Variable(-" MAXINT "," MAXINT ",'" $3 "')";
	if (match($0,"::_output_")) output[$3] = 1;
	if ($4 == "=") print "    int_eq(" $3 " , " $5 ")";
}

/^var range[(][-]*[0-9]+,1[+][-]*[0-9]+[)]:/{
	sub("range[(]","",$2);
	sub("[)]:","",$2);
	sub(",1[+]",",",$2);
	print "    " $3 " = Variable(" $2 ",'" $3 "')";
	if (match($0,"::_output_")) output[$3] = 1;
	if ($4 == "=") print "    int_eq(" $3 " , " $5 ")";
}

/^var {[-]*[0-9]+(,[-]*[0-9]+)*}:/{
	sub("{","[",$2);
	sub("}:","]",$2);
	print "    " $3 " = VariableInDomain(" $2 ",'" $3 "')";
	if (match($0,"::_output_")) output[$3] = 1;
	if ($4 == "=") print "    int_eq(" $3 " , " $5 ")";
}

/^array [[]range[(][-]*[0-9]+,1[+][-]*[0-9]+[)][]] of var bool:/{
	sub("[[]range[(]","",$2);
	sub("[)][]]","",$2);
	sub(",1[+]",",",$2);
	iinf = $2;
	sub(",.*","",iinf);
	iinf=int(iinf);
	if (iinf != 1) {
		print "ERROR ASSUMING MINIMUM INDEX VALUE EQUAL TO ONE IN ARRAY DEFINITION",$0;
		error = 1;
		exit(1);
	}
	isup = $2;
	sub(".*,","",isup);
	isup=int(isup);
	if (isup < iinf) {
		print "ERROR WRONG RANGE INDEX VALUE IN ARRAY DEFINITION",$0;
		error = 2;
		exit(2);
	}
	name = $6;
	if (match($0,"::_output_")) {
		output[ name ] = isup;
		outputstring[ name ] = substr($0, RSTART+RLENGTH);
	}
	if ($7 == "=") {
		printf "    " name " = ";
		sub(".*= ","");
		sub("::_output_.*","");
		gsub(" ","");
		print $0;
	} else {
		print "    " name " = VarArrayBoolean(" isup ",'" name "')";
	}
}

/^array [[]range[(][-]*[0-9]+,1[+][-]*[0-9]+[)][]] of var int:/{
	sub("[[]range[(]","",$2);
	sub("[)][]]","",$2);
	sub(",1[+]",",",$2);
	iinf = $2;
	sub(",.*","",iinf);
	iinf=int(iinf);
	if (iinf != 1) {
		print "ERROR ASSUMING MINIMUM INDEX VALUE EQUAL TO ONE IN ARRAY DEFINITION",$0;
		error = 1;
		exit(1);
	}
	isup = $2;
	sub(".*,","",isup);
	isup=int(isup);
	if (isup < iinf) {
		print "ERROR WRONG RANGE INDEX VALUE IN ARRAY DEFINITION",$0;
		error = 2;
		exit(2);
	}
	name = $6;
	
	if (match($0,"::_output_")) {
		output[ name ] = isup;
		outputstring[ name ] = substr($0, RSTART+RLENGTH);
	}
	if ($7 == "=") {
		printf "    " name " = ";
		sub(".*= ","");
		sub("::_output_.*","");
		gsub(" ","");
		print $0;
	} else {
		print "    " name " = VarArray(" isup ",-" MAXINT "," MAXINT ",'" name "')";
	}
}

/^array [[]range[(][-]*[0-9]+,1[+][-]*[0-9]+[)][]] of var range[(][-]*[0-9]+,1[+][-]*[0-9]+[)]:/{
	sub("[[]range[(]","",$2);
	sub("[)][]]","",$2);
	sub(",1[+]",",",$2);
	iinf = $2;
	sub(",.*","",iinf);
	iinf=int(iinf);
	if (iinf != 1) {
		print "ERROR ASSUMING MINIMUM INDEX VALUE EQUAL TO ONE IN ARRAY DEFINITION",$0;
		error = 1;
		exit(1);
	}
	isup = $2;
	sub(".*,","",isup);
	isup=int(isup);
	if (isup < iinf) {
		print "ERROR WRONG RANGE INDEX VALUE IN ARRAY DEFINITION",$0;
		error = 2;
		exit(2);
	}
	sub("range[(]","",$5);
	sub("[)]:","",$5);
	sub(",1[+]",",",$5);
	name = $6;
	if (match($0,"::_output_")) {
		output[ name ] = isup;
		outputstring[ name ] = substr($0, RSTART+RLENGTH);
	}
	if ($7 == "=") {
		printf "    " name " = ";
		sub(".*= ","");
		sub("::_output_.*","");
		gsub(" ","");
		print $0;
	} else {
		print "    " name " = VarArray(" isup "," $5 ",'" name "')";
	}
}

/^array [[]range[(][-]*[0-9]+,1[+][-]*[0-9]+[)][]] of var {[-]*[0-9]+(,[-]*[0-9]+)*}:/{
	sub("[[]range[(]","",$2);
	sub("[)][]]","",$2);
	sub(",1[+]",",",$2);
	iinf = $2;
	sub(",.*","",iinf);
	iinf=int(iinf);
	if (iinf != 1) {
		print "ERROR ASSUMING MINIMUM INDEX VALUE EQUAL TO ONE IN ARRAY DEFINITION",$0;
		error = 1;
		exit(1);
	}
	isup = $2;
	sub(".*,","",isup);
	isup=int(isup);
	if (isup < iinf) {
		print "ERROR WRONG RANGE INDEX VALUE IN ARRAY DEFINITION",$0;
		error = 2;
		exit(2);
	}
	sub("{","[",$5);
	sub("}:","]",$5);
	name = $6;
	if (match($0,"::_output_")) {
		output[ name ] = isup;
		outputstring[ name ] = substr($0, RSTART+RLENGTH);
	}
	if ($7 == "=") {
		printf "    " name " = ";
		sub(".*= ","");
		sub("::_output_.*","");
		gsub(" ","");
		print $0;
	} else {
		print "    " name " = [VariableInDomain(" $5 ",'" name "_" i "') for i in range(" isup ")]";
	}
}

/^constraint /{
	$1 = "   ";
    	$NF = $NF "";
	print $0;
}

/^solve / && / minimize /{
	print "    Minimize(" $NF ")";
	objective = $NF;
}

/^solve / && / maximize /{
	print "    Maximize(" $NF ")";
	objective = $NF;
}

END {
	if (!error) {
	output_vars = "";

	n = asorti(output,varnames);
	for (i=1; i<=n; i++) {
		e = varnames[i];
		if(e != objective) {
			if(length(output_vars) > 0) output_vars = output_vars ", ";
			output_vars = output_vars e;
		}
	}
	if(objective) {
		 if(length(output_vars) > 0) output_vars = output_vars ", ";
		 output_vars = output_vars objective;
	}
	print "    finalize_int_lin_eq(" output_vars ")";
	print "  except model.Contradiction:";
	print "    return None, None";
	print "  output_vars = (" output_vars ")";
	print "  return model, output_vars";

	print "";
	print "model, output_vars = get_model()";
	print "if model is None:";
	print "    print('=====UNSATISFIABLE=====')";
	print "    exit(0)";
	if(output_vars) print output_vars " = output_vars";
	if(DUMPMODEL) print "model.Dump('fzn.cfn')";
	if(!FREESEARCH) {
	    print "model.Option.Static_variable_ordering = True";
	    print "model.Option.QueueComplexity = True";
	}
	print "result = model.Solve(timeLimit = " TIMELIMIT ", showSolutions= " SHOWSOLUTIONS ", allSolutions= " ALLSOLUTIONS ")";
	print "if result is not None:";
	print "    sols = model.GetSolutions()";
	print "    for sol in sols:";
	#print "        print('')";
	for (i=1; i<=n; i++) {
		e = varnames[i];
		if (output[e]==1) {
			if (e != objective && e != "objective" && e != "obj") print "        print('" e " = ',get_value(sol[1]," e "),';')";
			else print "        print('" e " = ', sol[0],';')";
		} else {
			print "        print('" e " = " outputstring[e] ",',get_values(sol[1]," e "),');')";
		}
	}
	print "        print('----------')";
	if (objective){
		print "    if model.GetDDualBound() == result[1]:";
		print "        print('==========')";
	} else {
		print "    print('==========')";
	}
	print "    print('% Solutions', result[2])";
	#print "elif solver.is_unsat():"
	#print "    print('=====UNSATISFIABLE=====')"
	print "else:";
	print "    print('=====UNKNOWN=====')";
	print "    print('% Solutions', 0)";
		# print "print '% SolveTime', model.getTime()";
	print "print('% Nodes', model.GetNbNodes())";
	print "print('% Failures', model.GetNbBacktracks())";
	if (objective) print "if result is not None: print('% Objective', int(result[1]))";
	}
}

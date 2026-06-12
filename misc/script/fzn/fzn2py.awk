
# translate flatzinc format to python/Numberjack 

BEGIN {
	MAXINT = 10000;
	VERBOSITY = 0;
	SHOWSOLUTIONS = 3;
	TIMELIMIT = 1200;
	print "#! /usr/bin/env python"
	print "import time"
	print "import pytoulbar2 as tb2"
	system("cat " MZNNJ_DIR "/fzn2py.py");
	print "";
	print "def get_model():"
	print "    global model";
	print "    model = tb2.CFN(verbose=" VERBOSITY ")";
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
	print "    " $3 " = Variable(" $2 ",'" $3 "')";
	if (match($0,"::_output_")) output[$3] = 1;
	if ($4 == "=") print "    int_eq(" $3 " , " $5 ")";
}

/^array [[]range[(][-]*[0-9]+,1[+][-]*[0-9]+[)][]] of var bool:/{
	sub("[[]range[(]","",$2);
	sub("[)][]]","",$2);
	sub(",1[+]",",",$2);
	iinf = $2;
	sub(",.*","",iinf);
	if (iinf != 1) {
		print "ERROR ASSUMING MINIMUM INDEX VALUE EQUAL TO ONE IN ARRAY DEFINITION",$0;
		error = 1;
		exit(1);
	}
	isup = $2;
	sub(".*,","",isup);
	if (isup < iinf) {
		print "ERROR WRONG RANGE INDEX VALUE IN ARRAY DEFINITION",$0;
		error = 2;
		exit(2);
	}
	name = $6;
	print "    " name " = VarArray(" isup ",'" name "')";
	if (match($0,"::_output_")) {
		output[ name ] = isup;
		outputstring[ name ] = substr($0, RSTART+RLENGTH);
	}
	if ($7 == "=") {
		sub(".*= ","");
		sub("::_output_.*","");
		gsub(" ","");
		for (i=0; i<isup; i++) {
			print "    int_eq(" name "[" i "] , " $0 "[" i "])";
		}
	}
}

/^array [[]range[(][-]*[0-9]+,1[+][-]*[0-9]+[)][]] of var int:/{
	sub("[[]range[(]","",$2);
	sub("[)][]]","",$2);
	sub(",1[+]",",",$2);
	iinf = $2;
	sub(",.*","",iinf);
	if (iinf != 1) {
		print "ERROR ASSUMING MINIMUM INDEX VALUE EQUAL TO ONE IN ARRAY DEFINITION",$0;
		error = 1;
		exit(1);
	}
	isup = $2;
	sub(".*,","",isup);
	if (isup < iinf) {
		print "ERROR WRONG RANGE INDEX VALUE IN ARRAY DEFINITION",$0;
		error = 2;
		exit(2);
	}
	name = $6;
	print "    " name " = VarArray(" isup ",-" MAXINT "," MAXINT ",'" name "')";
	if (match($0,"::_output_")) {
		output[ name ] = isup;
		outputstring[ name ] = substr($0, RSTART+RLENGTH);
	}
	if ($7 == "=") {
		sub(".*= ","");
		sub("::_output_.*","");
		gsub(" ","");
		for (i=0; i<isup; i++) {
			print "    int_eq(" name "[" i "] , " $0 "[" i "])";
		}
	}
}

/^array [[]range[(][-]*[0-9]+,1[+][-]*[0-9]+[)][]] of var range[(][-]*[0-9]+,1[+][-]*[0-9]+[)]:/{
	sub("[[]range[(]","",$2);
	sub("[)][]]","",$2);
	sub(",1[+]",",",$2);
	iinf = $2;
	sub(",.*","",iinf);
	if (iinf != 1) {
		print "ERROR ASSUMING MINIMUM INDEX VALUE EQUAL TO ONE IN ARRAY DEFINITION",$0;
		error = 1;
		exit(1);
	}
	isup = $2;
	sub(".*,","",isup);
	if (isup < iinf) {
		print "ERROR WRONG RANGE INDEX VALUE IN ARRAY DEFINITION",$0;
		error = 2;
		exit(2);
	}
	sub("range[(]","",$5);
	sub("[)]:","",$5);
	sub(",1[+]",",",$5);
	name = $6;
	print "    " name " = VarArray(" isup "," $5 ",'" name "')";
	if (match($0,"::_output_")) {
		output[ name ] = isup;
		outputstring[ name ] = substr($0, RSTART+RLENGTH);
	}
	if ($7 == "=") {
		sub(".*= ","");
		sub("::_output_.*","");
		gsub(" ","");
		for (i=0; i<isup; i++) {
			print "    int_eq(" name "[" i "] , " $0 "[" i "])";
		}
	}
}

/^array [[]range[(][-]*[0-9]+,1[+][-]*[0-9]+[)][]] of var {[-]*[0-9]+(,[-]*[0-9]+)*}:/{
	sub("[[]range[(]","",$2);
	sub("[)][]]","",$2);
	sub(",1[+]",",",$2);
	iinf = $2;
	sub(",.*","",iinf);
	if (iinf != 1) {
		print "ERROR ASSUMING MINIMUM INDEX VALUE EQUAL TO ONE IN ARRAY DEFINITION",$0;
		error = 1;
		exit(1);
	}
	isup = $2;
	sub(".*,","",isup);
	if (isup < iinf) {
		print "ERROR WRONG RANGE INDEX VALUE IN ARRAY DEFINITION",$0;
		error = 2;
		exit(2);
	}
	sub("{","[",$5);
	sub("}:","]",$5);
	name = $6;
	print "    " name " = [Variable(" $5 ",'" name "_" i "') for i in range(" isup ")]";
	if (match($0,"::_output_")) {
		output[ name ] = isup;
		outputstring[ name ] = substr($0, RSTART+RLENGTH);
	}
	if ($7 == "=") {
		sub(".*= ","");
		sub("::_output_.*","");
		gsub(" ","");
		for (i=0; i<isup; i++) {
			print "    int_eq(" name "[" i "] , " $0 "[" i "])";
		}
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
	print "    output_vars = (" output_vars ")";
	print "    return model, output_vars";

	print "";
	print "model, output_vars = get_model()";
	if(output_vars) print output_vars " = output_vars";
	print "#model.Dump('fzn.cfn')";
	print "result = model.Solve(timeLimit = " TIMELIMIT ", showSolutions= "SHOWSOLUTIONS ")";
	print "if result is not None:";
	print "    print('')";
	for (i=1; i<=n; i++) {
		e = varnames[i];
		if (output[e]==1) {
			if (e != objective && e != "objective" && e != "obj") print "    print('" e " = ', result[0][" i "])";
			else print "    print('" e " = ', result[1])";
		} else {
			print "    print('" e " = " outputstring[e] ",',get_values(result[0]," e "),')')";
		}
	}
	print "    print('----------')";
	if (objective){
		print "    if model.GetDDualBound() == result[1]:"
		print "        print('==========')"
	} else {
		print "    print('==========')"
	}
	#print "elif solver.is_unsat():"
	#print "    print('=====UNSATISFIABLE=====')"
        print "else:"
	print "    print('=====UNKNOWN=====')"
	# print "print '% SolveTime', model.getTime()";
	print "print('% Nodes', model.GetNbNodes())";
	print "print('% Failures', model.GetNbBacktracks())";
	if (objective) print "if result is not None: print('% Objective', int(result[1]))";
	}
}

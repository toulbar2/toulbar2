/*
 * **************** Main function ***********************************
 */

#include "tb2solver.hpp"
#include "tb2pedigree.hpp"
#include "tb2bep.hpp"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>

//* definition of path separtor depending of OS '/'  => Unix ;'\' ==> windows
#ifdef WINDOWS
#define PATH_SEP_CHR '\\'
#define PATH_DELIM ";"
#else
#define PATH_SEP_CHR '/'
#define PATH_DELIM ":"
#endif

// USED FOR DEBUGGING PURPOSE
// under gdb: p ((BinaryConstraint *) constrs[13])->dump
// under gdb: p $2(constrs[13], myCout)
ostream myCout(cout.rdbuf());

#ifdef PARETOPAIR_COST
void initCosts(Cost c)
{
    if (ToulBar2::LcLevel > LC_FDAC) {cerr << "EDAC not implemented on Paretopair => force to FDAC." << endl; ToulBar2::LcLevel = LC_FDAC;}
    if (ToulBar2::vac) {cerr << "VAC not implemented on Paretopair." << endl; ToulBar2::vac = false; ToulBar2::minsumDiffusion = false;}
	if (ToulBar2::elimDegree >= 0 || ToulBar2::elimDegree_preprocessing >= 0) {cerr << "Variable elimination not implemented on Paretopair." << endl; ToulBar2::elimDegree = -1; ToulBar2::elimDegree_preprocessing = -1;}
    if (ToulBar2::btdMode>=1) {cerr << "BTD-like methods not implemented on Paretopair." << endl; ToulBar2::btdMode = 0;}
}
#endif


/* return current binary path extract from or argv[0] or from the env var $path the env var $path */

char* find_bindir(const char* bin_name, char* buffer, size_t buflen) {
        struct stat st;
        char* path, * tok;
        if(!stat(bin_name, &st)) {
                char* end = (char*) strrchr(bin_name, PATH_SEP_CHR);
		 static char bin_path[256];
                if(end) {
                        *end=0;
                        strncpy(buffer, bin_name, buflen);
                        sprintf(bin_path,"%s%c", buffer,PATH_SEP_CHR);
                } else {
                        strcpy(buffer, ".");
                    //path separator added to the path value
                        sprintf(bin_path,"%s%c", buffer,PATH_SEP_CHR);
                }
                return(bin_path);
        }
        path = strdup(getenv("PATH"));
        tok = strtok(path, PATH_DELIM);
        while(tok) {
                snprintf(buffer, buflen, "%s%c%s", tok, PATH_SEP_CHR, bin_name);
                if(!stat(buffer, &st)) {
                    static char bin_path[256];
                    strncpy(buffer, tok, buflen);
                    free(path);
                    sprintf(bin_path,"%s%c", buffer,PATH_SEP_CHR);
                    return bin_path;


                }
                tok = strtok(NULL, PATH_DELIM);
        }
        buffer[0]=0;
        return NULL;
}

bool localSearch(char *filename, Cost *upperbound, char *CurrentBinaryPath)
{
    string keyword;
// narcycsp file name
    static char fich[512];
    char * tname ;
    tname=tmpnam (NULL);
    strcpy(fich,tname);
    cout << "incop narycsp output file :" << fich << endl;

    char line[1024];
    //int pid = get_pid();
    Cost tmpUB=MIN_COST;
    int varIndex;
    int value;
    map<int, int> solution;
    map<int, int> bestsol;
    Cost bestcost = MAX_COST;
    Cost solcost = MAX_COST;
    // * get narycsp path into system variable call Narycsp_path
    char* Narycsp_Path = getenv ("NARYCSP");

  // run local search solver INCOP from Bertrand Neveu on wcsp format
    if (Narycsp_Path!=NULL)
    {
	printf ("env variable NARYCPS founded: narycsp path = %s \n",Narycsp_Path);
   sprintf(line,"%s%cnarycsp %s %s 0 1 5 idwa 100000 cv v 0 200 1 0 0", Narycsp_Path,PATH_SEP_CHR, fich, filename);
   // sprintf(line,"%s%cautonarycsp %s %s 0 1 3 idwa 200 cv v 0 10 10 0 0 100", Narycsp_Path,PATH_SEP_CHR, fich, filename);

    } else {

    printf ("narycsp default path : %s \n",CurrentBinaryPath);

 sprintf(line,"%snarycsp %s %s 0 1 5 idwa 100000 cv v 0 200 1 0 0", CurrentBinaryPath, fich, filename);

 // sprintf(line,"%s%cautonarycsp %s %s 0 1 3 idwa 200 cv v 0 10 10 0 0 100", CurrentBinaryPath,PATH_SEP_CHR, fich, filename);

    }

    if(system(line))
    {
        printf ("exec error from narcycsp call: %s\n",line);//erreur
    }


    // open the file
    ifstream file(fich);
    if (!file) {
        return false;
    }

    while (file) {
        file >> keyword;
        if (keyword == "Meilleur") {
            file >> keyword;
            if (keyword == "essai") {
                file >> keyword; // skip ":"
                file >> tmpUB;
                if (tmpUB>=MIN_COST && tmpUB < *upperbound) *upperbound = tmpUB;
                if (ToulBar2::writeSolution && *upperbound == bestcost) {
                    ofstream sol("sol");
                    if (file) {
                        for (map<int,int>::iterator iter = bestsol.begin(); iter != bestsol.end(); ++iter) {
                            sol << " " << (*iter).second;
                        }
                        sol << endl;
                    }
                }
                break;
            }
        }
        if (keyword == "variable") {
            file >> varIndex;
            file >> keyword; // skip ":"
            file >> value;
            solution[varIndex] = value;
        }
        if (keyword == "verification") {
            file >> solcost;
            if (solcost < bestcost) {
                bestcost = solcost;
                bestsol = solution;
            }
        }
    }
    file.close();

/* deletion of tmp file created by narycsp system call */

 if( remove(fich ) != 0 )
    perror( "narycsp File: Error deleting file" );
  else
    puts( "narycsp tmp File successfully deleted" );

      return true;
}

void help_msg(char *toulbar2filename)
{
        cerr << "*************************" << endl;
        cerr << "* ToulBar2 Help Message *" << endl;
        cerr << "*************************" << endl;
        cerr << toulbar2filename << " version " << ToulBar2::version << endl;
        cerr << endl;
        cerr << "Command line is:" << endl;
        cerr << toulbar2filename << " filename [upperbound] [options]" << endl;
        cerr << endl;
#ifndef MENDELSOFT
        cerr << "Available problem formats (specified by the filename extension) are:" << endl;
        cerr << "   *.wcsp : Weighted CSP format (see SoftCSP web site)" << endl;
#ifdef XMLFLAG
        cerr << "   *.xml : CSP and weighted CSP in XML format XCSP 2.1";
#ifdef MAXCSP
        cerr << " (Max-CSP only)";
#endif
		cerr << endl;
#endif
        cerr << "   *.uai : Bayesian network and Markov Random Field format (see UAI'08 evaluation workshop) followed by an optional evidence filename (perform MPE task)" << endl;
        cerr << "   *.pre  : pedigree format (see doc/MendelSoft.txt for Mendelian error correction)" << endl;
        cerr << "   *.bep  : satellite scheduling format (CHOCO benchmark)" << endl << endl;

        cerr << "Alternatively one can call the random problem generator: " << endl;
		cerr << "     bin-{n}-{m}-{p1}-{p2}-{seed}        p1 is the tightness in percentage %" << endl;
		cerr << "                                         p2 is the num of binary cost functions to include" << endl;
		cerr << "                                         the seed parameter is optional" << endl;

		cerr << "or:                                                                               " << endl;
		cerr << "     binsub-{n}-{m}-{p1}-{p2}-{p3}-{seed} binary random & submodular cost functions" << endl;
		cerr << "                                         p1 is the tightness in percentage % of random cost functions" << endl;
		cerr << "                                         p2 is the num of binary cost functions to include" << endl;
		cerr << "                                         p3 is the percentage % of submodular cost functions among p2 cost functions" << endl;
		cerr << "                                         (plus 10 permutations of two randomly-chosen values for each domain)" << endl;
		cerr << "or:                                                                               " << endl;
		cerr << "     tern-{n}-{m}-{p1}-{p2}-{p3}-{seed}  p3 is the num of ternary cost functions" << endl;
        cerr << "or:                                                                               " << endl;
		cerr << "     nary-{n}-{m}-{p1}-{p2}-{p3}...{pn}-{seed}  pn is the num of n-ary cost functions" << endl;
        cerr << endl;
#endif
        cerr << "An initial upperbound is optional (default value is " << MAX_COST << ")" << endl;
        cerr << endl;
        cerr << "Available options are (WARNING! each letter is meaningful, do not use a \"-\" before each option, use symbol \":\" before an option to remove a default option):" << endl;
        cerr << "   v : verbosity (repeat this option letter to increase the verbosity level)" << endl;
        cerr << "   s : show each solution found" << endl;
#ifndef MENDELSOFT
		cerr << "   w : write last solution found in filename \"sol\"" << endl;
#else
		cerr << "   w[mode] : write last solution found" << endl;
		cerr << "               mode=0: save pedigree with erroneous genotypings removed" << endl;
		cerr << "               mode=1: save pedigree with erroneous genotypings corrected" << endl;
		cerr << "               mode=2: save pedigree with erroneous genotypings corrected and missing genotypes of informative individuals inferred" << endl;
        cerr << "   g : sort pedigree by increasing generation number and if equal by increasing individual number" << endl;
        cerr << "   u[integer] : add a penalty weight (must use option y also) on genotyped individuals depending on the number of their genotyped children in order to penalize genotyping removals if the number of genotyped children is strictly greater than a given threshold" << endl;
        cerr << "   y [genotypinpErrorRate probabilityPrecision genotypePriorMode]  : pedigree solved by Bayesian MPE" << endl;
		cerr << "        y must be the last option in the command line possibly followed by three arguments:" << endl;
        cerr << "               genotypingErrorRate is a prior uniform probability of genotyping errors (default value is " << ToulBar2::errorg << ")" << endl;
        cerr << "               probabilityPrecision is a conversion factor (a power of ten) for representing fixed point numbers (default value is " << ToulBar2::resolution << ")" << endl;
        cerr << "               genotypePriorMode selects the prior mode for allele probability distribution (default value is " << ToulBar2::foundersprob_class << ")" << endl;
        cerr << "                   = 0 : uniform allele probability distribution" << endl;
        cerr << "                   = 1 : allele probability distribution read from pedigree data" << endl;
        cerr << "                   = 2 p1 p2 p3...: allele probability distribution given explicitely in the command line" << endl << endl;
#endif
#ifndef MENDELSOFT
        cerr << "   b : perform binary branching always instead of binary branching for interval domains and n-ary branching for enumerated domains";
		if (ToulBar2::binaryBranching) cerr << " (default option)";
		cerr << endl;
        cerr << "   c : perform binary branching with last conflict backjumping variable ordering heuristic";
		if (ToulBar2::lastConflict) cerr << " (default option)";
		cerr << endl;
        cerr << "   q : weighted degree variable ordering heuristic";
		if (ToulBar2::weightedDegree) cerr << " (default option)";
		cerr << endl;
        cerr << "   d : dichotomic branching instead of binary branching when current domain size is strictly greater than " << ToulBar2::dichotomicBranchingSize;
		if (ToulBar2::dichotomicBranching) cerr << " (default option)";
		cerr << endl;
        cerr << "   e[integer] : boosting search with variable elimination of small degree (less than or equal to 3) (default value is " << ToulBar2::elimDegree << ")" << endl;
        cerr << "   p[integer] : preprocessing only: general variable elimination of degree less than or equal to the given value (default value is " << ToulBar2::elimDegree_preprocessing << ")" << endl;
        cerr << "   t : preprocessing only: project ternary cost functions on binary cost functions and apply 3-consistency (can be very slow)";
		if (ToulBar2::preprocessTernary) cerr << " (default option)";
		cerr << endl;
        cerr << "   h : preprocessing only: project ternary cost functions on binary cost functions following a heuristic (to be used in conjunction with option \"t\")";
		if (ToulBar2::preprocessTernaryHeuristic) cerr << " (default option)";
		cerr << endl;
        cerr << "   m : preprocessing only: minimum degree (if compiled with BOOST)/user-specified re-ordering of variables (in conjunction with options \"p\" and \"O\")";
		if (ToulBar2::elimOrderType != ELIM_NONE) cerr << " (default option)";
		cerr << endl;
        cerr << "   o : ensure optimal worst-case time complexity of DAC and EAC (can be costly in practice)";
		if (ToulBar2::QueueComplexity) cerr << " (default option)";
		cerr << endl;
        cerr << "   k[integer] : soft local consistency level (NC=0, AC=1, DAC=2, FDAC=3, EDAC=4) (default value is " << ToulBar2::LcLevel << ")" << endl;
        cerr << "   l : limited discrepancy search";
		if (ToulBar2::lds) cerr << " (default option)";
		cerr << endl;
#ifndef WINDOWS
        cerr << "   i : initial upperbound found by INCOP local search solver (filename \"./misc/bin/linux/narycsp\")" << endl;
#endif
        cerr << "   z : save problem in wcsp format in filename \"problem.wcsp\" (repeat this option letter to save the current problem after preprocessing)" << endl;
        cerr << "   Z : debug mode (save problem at each node if verbosity option set!)" << endl;
        cerr << "   x : load a solution from filename \"sol\"" << endl << endl;
        cerr << "   M[integer] : preprocessing only: Min Sum Diffusion algorithm (default number of iterations is " << ToulBar2::minsumDiffusion << ")" << endl;
        cerr << "   A[integer] : enforce VAC at each search node with a search depth less than a given value (default value is " << ToulBar2::vac << ")" << endl;
        cerr << "   T[integer] : threshold cost value for VAC (default value is " << ToulBar2::costThreshold << ")" << endl;
        cerr << "   P[integer] : threshold cost value for VAC during the preprocessing phase (default value is " << ToulBar2::costThresholdPre << ")" << endl;
        cerr << "   C[integer] : multiply all costs by this number (default value is " << ToulBar2::costMultiplier << ")" << endl;
        cerr << "   S : preprocessing only: perform singleton consistency (only in conjunction with option \"A\")";
		if (ToulBar2::singletonConsistency) cerr << " (default option)";
		cerr << endl;
        cerr << "   V : VAC-based value ordering heuristic";
		if (ToulBar2::vacValueHeuristic) cerr << " (default option)";
		cerr << endl << endl;

        cerr << "   B[integer] : (0) DFBB, (1) BTD, (2) RDS-BTD, (3) RDS-BTD with path decomposition instead of tree decomposition (default value is " << ToulBar2::btdMode << ")" << endl;
		cerr << "   O[filename ] : read a variable elimination order from a file in order to build a tree decomposition (WARNING! DO NOT PUT A SPACE BETWEEN \"O\" AND FILENAME)" << endl;
		cerr << "                  (if not specified, then use the variable order in which variables appear in the problem file)" << endl;
        cerr << "   j[integer] : split large clusters into a chain of smaller embedded clusters with a number of proper variables less than this number" << endl;
        cerr << "                (use options \"B3j1\" for pure RDS, use value 0 for no splitting) (default value is " << ToulBar2::splitClusterMaxSize << ")" << endl;
        cerr << "   r[integer] : limit on maximum cluster separator size (merge cluster with its father otherwise, use a negative value for no limit) (default value is " << ToulBar2::maxSeparatorSize << ")" << endl;
        cerr << "   X[integer] : limit on minimum number of proper variables in a cluster (merge cluster with its father otherwise, use a zero for no limit) (default value is " << ToulBar2::minProperVarSize << ")" << endl;
        cerr << "   E : merge leaf clusters with their fathers if small local treewidth (in conjunction with option \"e\")";
		if (ToulBar2::boostingBTD) cerr << " (default option)";
		cerr << endl;
        cerr << "   R[integer] : choice for a specific root cluster number" << endl;
        cerr << "   I[integer] : choice for solving only a particular rooted cluster subtree" << endl << endl;
        cerr << "   a : find all solutions (or count the number of zero-cost satisfiable solutions in conjunction with BTD)";
		if (ToulBar2::allSolutions) cerr << " (default option)";
		cerr << endl;
        cerr << "   D : approximate satisfiable solution count with BTD";
		if (ToulBar2::approximateCountingBTD) cerr << " (default option)";
		cerr << endl;
#endif
}

int main(int argc, char **argv)
{
    bool localsearch = false;
    bool certificate = false;
    char buf [512];
    char* CurrentBinaryPath = find_bindir(argv[0], buf, 512); // current binary path search
    cout << CurrentBinaryPath<<"toulbar2"<<"  version : " << ToulBar2::version << ", copyright (c) INRA 2010"<<endl;
//     cout << "Toulbar2 Path="<<CurrentBinaryPath<<'\n';

    if (argc <= 1) {
        cerr << "Missing a problem filename as first argument!" << endl;
        cerr << endl;
	    help_msg(argv[0]);
        exit(0);
    }

    string strext;
	string strfile(argv[1]);
    unsigned int pos = strfile.find_last_of(".");
    if(pos < strfile.size()) strext = strfile.substr(pos,strfile.size());

	int iniarg = 2;
    if (strstr(strext.c_str(),".xml")) { ToulBar2::xmlflag = true; ToulBar2::writeSolution = true; }
    if (strstr(strext.c_str(),".uai")) {
    	ToulBar2::uai = true; ToulBar2::bayesian = true; ToulBar2::writeSolution = true;
		if(argc > 2 && strstr(argv[2],".evid")) {
			ToulBar2::evidence_file = string(argv[2]);
			iniarg = 3;
		}
    }


    char* ch;
    ToulBar2::verbose = 0;

    for (int i=iniarg; i<argc; i++) {
        if ( (ch = strchr(argv[i],'O')) ) {
        	char buf[80];
        	sprintf(buf,"%s",&ch[1]);
        	if(ToulBar2::varOrder) delete [] ToulBar2::varOrder;
        	ToulBar2::varOrder = new char [ strlen(buf) + 1 ];
	       	sprintf(ToulBar2::varOrder, "%s",buf);
			continue; // skip current argument in order to not search for other options inside filename
    	}
        if ( (ch = strchr(argv[i],'B')) ) {
        	int mode = atoi(&ch[1]);
        	if(mode >= 0) {
			  ToulBar2::btdMode = mode;
			}
        }
        if ( (ch = strchr(argv[i],'R')) ) {
        	int root = atoi(&ch[1]);
        	if(root > 0) ToulBar2::btdRootCluster = root;
        }
        if ((ch = strchr(argv[i],'I'))) {
			int subcluster = atoi(&ch[1]);
		    if(subcluster >= 1) ToulBar2::btdSubTree = subcluster;
		}
        if ( (ch = strchr(argv[i],'j')) ) {
        	int cmaxsize = atoi(&ch[1]);
        	if(cmaxsize >= 1) ToulBar2::splitClusterMaxSize = cmaxsize;
        }
        if ( (ch = strchr(argv[i],'E')) ) {
		  ToulBar2::boostingBTD = true;
		}
        if ( (ch = strchr(argv[i],'r')) ) {
        	int sepmaxsize = atoi(&ch[1]);
        	if(sepmaxsize >= -1) ToulBar2::maxSeparatorSize = sepmaxsize;
        }
        if ( (ch = strchr(argv[i],'X')) ) {
        	int minpvarsize = atoi(&ch[1]);
        	if(minpvarsize >= 0) ToulBar2::minProperVarSize = minpvarsize;
        }
        for (int j=0; argv[i][j] != 0; j++) if (argv[i][j]=='v') ToulBar2::verbose++;
        for (int j=0; argv[i][j] != 0; j++) if (argv[i][j]=='Z') ToulBar2::debug++;
        if (strchr(argv[i],'s')) ToulBar2::showSolutions = true;
        if ( (ch = strchr(argv[i],'w')) ) {
		  ToulBar2::writeSolution = true;
		  int correct = atoi(&ch[1]);
		  ToulBar2::pedigreeCorrectionMode = 0;
		  if((correct > 0) && (correct <= 2)) ToulBar2::pedigreeCorrectionMode = correct;
        }
        if ( (ch = strchr(argv[i],'u')) ) {
		  ToulBar2::pedigreePenalty = 1;
		  int penaltyThreshold = atoi(&ch[1]);
		  if(penaltyThreshold >= 1) ToulBar2::pedigreePenalty = penaltyThreshold;
		}
		if (strchr(argv[i],'a')) ToulBar2::allSolutions = true;
		if (strchr(argv[i],'D')) {ToulBar2::approximateCountingBTD = true; ToulBar2::allSolutions = true; ToulBar2::btdMode = 1;}
        if ((ch = strchr(argv[i],'b'))) {if (ch[-1]==':') { ToulBar2::binaryBranching = false; } else { ToulBar2::binaryBranching = true; }}
        if ((ch = strchr(argv[i],'c'))) {if (ch[-1]==':') { ToulBar2::lastConflict = false; } else { ToulBar2::binaryBranching = true; ToulBar2::lastConflict = true; }}
        if ((ch = strchr(argv[i],'d'))) {if (ch[-1]==':') { ToulBar2::dichotomicBranching = false; } else  { ToulBar2::binaryBranching = true; ToulBar2::dichotomicBranching = true; }}
        if ((ch = strchr(argv[i],'q'))) {if (ch[-1]==':') ToulBar2::weightedDegree = false; else ToulBar2::weightedDegree = true;}
        if ( (ch = strchr(argv[i],'e')) ) {
        	if (ch[-1]==':') ToulBar2::elimDegree = -1; else ToulBar2::elimDegree = 3;
			char *tmpendchr = NULL;
        	int ndegree = strtol(&ch[1], &tmpendchr, 10);
        	if((tmpendchr != &ch[1]) && (ndegree >= 0) && (ndegree <= 3)) ToulBar2::elimDegree = ndegree;
        }
        if ( (ch = strchr(argv[i],'p')) ) {
        	if (ch[-1]==':') ToulBar2::elimDegree_preprocessing = -1; else ToulBar2::elimDegree_preprocessing = 3;
        	int ndegree = atoi(&ch[1]);
        	if(ndegree > 0) ToulBar2::elimDegree_preprocessing = ndegree;
        }
        if ( (ch = strchr(argv[i],'M')) ) {
        	if (!ToulBar2::vac) ToulBar2::vac = 1;
         	ToulBar2::minsumDiffusion = 1000;
        	int nit = atoi(&ch[1]);
        	if(nit > 0) ToulBar2::minsumDiffusion = nit;
        }
        if ( (ch = strchr(argv[i],'A')) ) {
		  if (ch[-1]==':') ToulBar2::vac = 0; else ToulBar2::vac = 1;
		  int depth = atoi(&ch[1]);
		  if(depth >= 1) ToulBar2::vac = depth;
		}
        if ( (ch = strchr(argv[i],'T')) ) {
        	Cost ct = string2Cost(&ch[1]);
        	if(ct > UNIT_COST) ToulBar2::costThreshold = ct;
        }
        if ( (ch = strchr(argv[i],'P')) ) {
        	Cost ct = string2Cost(&ch[1]);
        	if(ct > UNIT_COST) ToulBar2::costThresholdPre = ct;
        }
        /*if ( (ch = strchr(argv[i],'R')) ) {
        	Cost ct = string2Cost(&ch[1]);
        	if(ct > UNIT_COST) ToulBar2::relaxThreshold = ct;
        }*/
        if ( (ch = strchr(argv[i],'C')) ) {
        	Cost co = string2Cost(&ch[1]);
        	if(co > MIN_COST) ToulBar2::costMultiplier = co;
        }
        if (strchr(argv[i],'S')) ToulBar2::singletonConsistency = true;
        if (strchr(argv[i],'V')) ToulBar2::vacValueHeuristic = true;
        if (strchr(argv[i],'t')) ToulBar2::preprocessTernary = true;
        if (strchr(argv[i],'h')) { ToulBar2::preprocessTernary = true; ToulBar2::preprocessTernaryHeuristic = true; }
        if (strchr(argv[i],'m')) ToulBar2::elimOrderType = MIN_DEGREE;
        if (strchr(argv[i],'o')) ToulBar2::QueueComplexity = true;
        if (strchr(argv[i],'l')) ToulBar2::lds = true;
        if (strchr(argv[i],'i')) localsearch = true;
        if ( (ch = strchr(argv[i],'k')) ) {
        	ToulBar2::LcLevel = LC_EDAC;
        	LcLevelType lclevel = (LcLevelType) atoi(&ch[1]);
        	if((lclevel >= LC_NC) && (lclevel < LC_THEMAX)) ToulBar2::LcLevel = lclevel;
        }
        for (int j=0; argv[i][j] != 0; j++) if (argv[i][j]=='z') ToulBar2::dumpWCSP++;
        if (strchr(argv[i],'x')) certificate = true;
        if (strchr(argv[i],'g')) ToulBar2::generation = true;
        if (strchr(argv[i],'y')) {
        	ToulBar2::bayesian = true;
        	float f;
            int pos = i + 1;
        	if(argc > pos) { sscanf(argv[pos++],"%f",&f); ToulBar2::errorg = f; }
        	if(argc > pos) sscanf(argv[pos++],"%d",&ToulBar2::resolution);
        	if(argc > pos) sscanf(argv[pos++],"%d",&ToulBar2::foundersprob_class); // 0: 			equal frequencies
								   										           // 1: 			probs depending on the frequencies found in the current pedigree problem
										 									       // otherwise:    read probability distribution from command line
            while (argc > pos) { sscanf(argv[pos++],"%f",&f); ToulBar2::allelefreqdistrib.push_back(f); }
        }
    }

    if (ToulBar2::elimDegree_preprocessing > 0 && (ToulBar2::showSolutions || ToulBar2::writeSolution)) {
	  cout << "Warning! Cannot show/save solutions if general variable elimination used in preprocessing." << endl;
	  ToulBar2::showSolutions = false;
	  ToulBar2::writeSolution = false;
	}
    if (ToulBar2::elimDegree_preprocessing > 2 && ToulBar2::varOrder && ToulBar2::elimOrderType == ELIM_NONE) {
	  cout << "Warning! Re-order variables from file \"" << ToulBar2::varOrder << "\"" << endl;
	  ToulBar2::elimOrderType = MIN_DEGREE;
	}
    if (ToulBar2::approximateCountingBTD && ToulBar2::btdMode != 1) {
		  cout << "Warning! Cannot find an approximation of solution count without BTD." << endl;
		  ToulBar2::approximateCountingBTD = false;
		  ToulBar2::allSolutions = false;
		}
	if (ToulBar2::allSolutions && ToulBar2::btdMode > 1) {
	  cout << "Warning! Cannot find all solutions with RDS-like search methods." << endl;
	  ToulBar2::allSolutions = false;
	}
	if (ToulBar2::allSolutions && ToulBar2::elimDegree >= 0) {
	  //	  cout << "Warning! Cannot count all solutions with variable elimination during search (except with degree 0 for #BTD)" << endl;
	  if (!ToulBar2::btdMode) {
		ToulBar2::elimDegree = -1;
	  } else {
		ToulBar2::elimDegree = 0;
	  }
	}
	if (ToulBar2::lds && ToulBar2::btdMode >= 1) {
	  cout << "Warning! Limited Discrepancy Search not compatible with BTD-like search methods." << endl;
	  ToulBar2::lds = false;
	}
	if (!ToulBar2::binaryBranching && ToulBar2::btdMode >= 1) {
	  cout << "Warning! n-ary branching not implemented with BTD-like search methods => force binary branching." << endl;
	  ToulBar2::binaryBranching = true;
	}
	if (ToulBar2::vac > 1 && ToulBar2::btdMode >= 1) {
	  cout << "Warning! VAC not implemented with BTD-like search methods during search => VAC in preprocessing only." << endl;
	  ToulBar2::vac = 1;
	}

	Cost c = (argc >= 3)?string2Cost(argv[2]):MAX_COST;
    if (c <= MIN_COST) c = MAX_COST;

    if (localsearch && !strstr(strext.c_str(),".pre")) {
#ifndef WINDOWS
    if (localSearch(argv[1],&c,CurrentBinaryPath)) {
            cout << "Initial upperbound: " << c << endl;

        } else cerr << "INCOP solver narycsp not found in:"<<CurrentBinaryPath << endl;

#else
          cerr << "Initial upperbound: INCOP not compliant with Windows OS" << endl;

#endif
    }
    if (c==MIN_COST) {
        cout << "Initial upperbound equal to zero!" << endl;
        cout << "No solution found by initial propagation!" << endl;
        cout << "end." << endl;
        return 0;
    }

    initCosts(c);
    Solver solver(STORE_SIZE, c);

    bool randomproblem = false;
	bool forceSubModular = false;

    int n = 10;
    int m = 2;
    int seed = 3;
    vector<int> p;
#ifdef MENDELSOFT
    ToulBar2::pedigree = new Pedigree;
#else
    if(!strchr(argv[1],'.') && !strstr(argv[1],"bEpInstance")) {
    	int pn[10];
    	int narities = 0;
    	if(strstr(argv[1],"bin"))  { randomproblem = true; sscanf(argv[1], "bin-%d-%d-%d-%d-%d", &n, &m, &pn[0], &pn[1],&seed); narities = 2; }
    	if(strstr(argv[1],"binsub"))  { forceSubModular = true; randomproblem = true; sscanf(argv[1], "binsub-%d-%d-%d-%d-%d-%d", &n, &m, &pn[0], &pn[1], &pn[2], &seed); narities = 2; }
    	if(strstr(argv[1],"tern")) { randomproblem = true; sscanf(argv[1], "tern-%d-%d-%d-%d-%d-%d", &n, &m, &pn[0], &pn[1], &pn[2],&seed); narities = 3; }
    	if(strstr(argv[1],"nary")) {
    		randomproblem = true;
       		char* pch = strtok (argv[1],"-");
       		pch = strtok (NULL, "-"); n = atoi(pch);
       		pch = strtok (NULL, "-"); m = atoi(pch);

			while (pch != NULL) {
			   pch = strtok (NULL, "-");
			   if(pch != NULL) {
				   pn[narities] = atoi(pch);
				   narities++;
			   }
			}
			narities--;
			seed = pn[narities];
    	}
		if(pn[0] > 100) { cout << pn[0] << " tightness is a percentage" << endl; pn[0] = 100; }
		for(int i=0;i<narities;i++) p.push_back( pn[i] );
		if (forceSubModular) p.push_back( pn[narities] );
		if(narities == 0) {
		  help_msg(argv[0]);
		  exit(0);
		}
    }
    if (strstr(strext.c_str(),".pre")) ToulBar2::pedigree = new Pedigree;
    if (strstr(strext.c_str(),".bep") || strstr(argv[1],"bEpInstance")) ToulBar2::bep = new BEP;
#endif
    try {
        if(randomproblem)    solver.read_random(n,m,p,seed,forceSubModular);
        else 		         solver.read_wcsp(argv[1]);

		ToulBar2::startCpuTime = cpuTime();

        if (certificate) solver.read_solution("sol");
        if (ToulBar2::dumpWCSP==1) solver.dump_wcsp("problem.wcsp");
		else if (!certificate || ToulBar2::btdMode>=2) {
		  if (CSP(solver.getWCSP()->getLb(), solver.getWCSP()->getUb())) {
			ToulBar2::LcLevel = LC_AC;
		  }
		  solver.solve();
		}
    } catch (Contradiction) {
        cout << "No solution found by initial propagation!" << endl;
    }
    if(!ToulBar2::xmlflag && !ToulBar2::uai) cout << "end." << endl;


    // for the competition it was necessary to write a file with the optimal sol
	/*char line[80];
    string strfile(argv[1]);
    int pos = strfile.find_last_of(".");
    string strfilewcsp = strfile.substr(0,pos) + ".ub";
    sprintf(line,"echo %d > %s",(int)solver.getWCSP()->getUb(),strfilewcsp.c_str());
    system(line); */

    return 0;
}

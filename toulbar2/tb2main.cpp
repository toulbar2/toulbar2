/*
 * **************** Main function ***********************************
 */

#include "tb2solver.hpp"
#include "tb2pedigree.hpp"
#include "tb2bep.hpp"
#include <string.h>

// used for debugging purpose.
// under gdb: p ((BinaryConstraint *) constrs[13])->dump
// under gdb: p $2(constrs[13], myCout)
ostream myCout(cout.rdbuf());

#ifdef PARETOPAIR_COST
void initCosts(Cost c)
{
    if (ToulBar2::LcLevel > LC_FDAC) {cerr << "EDAC not implemented on Paretopair => force to FDAC." << endl; ToulBar2::LcLevel = LC_FDAC;}
    if (ToulBar2::vac) {cerr << "VAC not implemented on Paretopair." << endl; ToulBar2::vac = false; ToulBar2::minsumDiffusion = false;}
	if (ToulBar2::elimDegree >= 0 || ToulBar2::elimDegree_preprocessing >= 0) {cerr << "Variable elimination not implemented on Paretopair." << endl; ToulBar2::elimDegree = -1; ToulBar2::elimDegree_preprocessing = -1;}
    if (ToulBar2::btdMode>=1) {cerr << "BTD-like methods not implemented on Paretopair." << endl; ToulBar2::btdMode = 0; ToulBar2::varOrder = NULL;}
}
#endif

bool localSearch(char *filename, Cost *upperbound)
{
    string keyword;
    const char *fich = "resincop";
    char line[1024];
    Cost tmpUB=MIN_COST;
    int varIndex;
    int value;
    map<int, int> solution;
    map<int, int> bestsol;
    Cost bestcost = MAX_COST;
    Cost solcost = MAX_COST;
    
    // run local search solver INCOP from Bertrand Neveu on wcsp format
    sprintf(line,"./narycsp %s %s 0 1 5 idwa 100000 cv v 0 200 1 0 0", fich, filename);
    system(line);

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
    sprintf(line,"rm -f %s",fich);
    system(line);
    
    return true;
}

int main(int argc, char **argv)
{
    bool localsearch = false;
    bool saveproblem = false;
    bool certificate = false;
    if (argc <= 1) {
        cerr << argv[0] << " version " << ToulBar2::version << endl;
        cerr << "Missing a problem filename as first argument!" << endl;
        cerr << endl;
        cerr << argv[0] << " filename upperbound options" << endl;
        cerr << endl;
#ifndef MENDELSOFT
        cerr << "Available problem formats (specified by the filename extension):" << endl;
        cerr << "   *.wcsp : wcsp format" << endl;
#ifdef XMLFLAG
        cerr << "   *.xml : XML format XCSP 2.1";
#ifdef MAXCSP
        cerr << " (Max-CSP only)";
#endif
		cerr << endl;
#endif
        cerr << "   *.uai : Bayesian network (UAI-08) format followed by an optional evidence filename (MPE task)" << endl;
        cerr << "   *.pre  : pedigree format" << endl << endl;
        cerr << "   *.bep  : BEP (satellite scheduling) format" << endl << endl;

        cerr << "Alternatively one call the random problem generator: " << endl; 
		cerr << "     bin-{n}-{m}-{p1}-{p2}-{seed}        p1 is the tightness in percentage %" << endl; 
		cerr << "                                         p2 is the num of binary constraints to include" << endl;
		cerr << "                                         the seed parameter is optional" << endl;
		 
		cerr << "or:                                                                               " << endl;            
		cerr << "     binsub-{n}-{m}-{p1}-{p2}-{seed}     binary submodular cost functions (p1 unused)" << endl; 
		cerr << "                                         (plus 10 permutations of two randomly-chosen values for each domain)" << endl; 
		cerr << "or:                                                                               " << endl;            
		cerr << "     tern-{n}-{m}-{p1}-{p2}-{p3}-{seed}  p3 is the num of ternary constraints" << endl; 
        cerr << "or:                                                                               " << endl;            
		cerr << "     nary-{n}-{m}-{p1}-{p2}-{p3}...{pn}-{seed}   " << endl; 
        cerr << endl;
#endif
        cerr << "Initial upperbound is optional (default value: " << MAX_COST << ")" << endl;
        cerr << endl;
        cerr << "Available options (use symbol \":\" before an option letter to remove a default option):" << endl;
        cerr << "   v : verbosity (repeat this option to increase the verbosity level)" << endl;
        cerr << "   s : show each solution found" << endl;
        cerr << "   g : sort pedigree by increasing generation number and if equal by increasing individual number" << endl;
        cerr << "   u[integer] : add a penalty weight (must use option y also) on genotyped individuals depending on the number of their genotyped children in order to penalize genotyping removals if the number of genotyped children is strictly greater than a given threshold" << endl;
		cerr << "   w[mode] : write last solution found" << endl;
		cerr << "     and for pedigree problems:" << endl;
		cerr << "               mode=0: save pedigree with erroneous genotypings removed" << endl;
		cerr << "               mode=1: save pedigree with erroneous genotypings corrected" << endl;
		cerr << "               mode=2: save pedigree with erroneous genotypings corrected and missing genotypes of informative individuals inferred" << endl;
        cerr << "   y [genotypinpErrorRate probabilityPrecision genotypePriorMode]  : pedigree solved by Bayesian MPE" << endl;
		cerr << "        y must be the last option in the command line followed by three arguments:" << endl;
        cerr << "               genotypingErrorRate is a prior uniform probability of genotyping errors (default value: " << ToulBar2::errorg << ")" << endl;
        cerr << "               probabilityPrecision is a conversion factor (a power of ten) for representing fixed point numbers (default value: " << ToulBar2::resolution << ")" << endl;
        cerr << "               genotypePriorMode selects the prior mode for allele probability distribution (default value: " << ToulBar2::foundersprob_class << ")" << endl;
        cerr << "                   = 0 : uniform allele probability distribution" << endl;
        cerr << "                   = 1 : allele probability distribution read from pedigree data" << endl;
        cerr << "                   = 2 p1 p2 p3...: allele probability distribution given explicitely in the command line" << endl;
#ifndef MENDELSOFT
        cerr << "   a : find all solutions" << endl;
        cerr << "   b : binary branching always instead of binary branching for interval domains and n-ary branching for enumerated domains (default option)" << endl;
        cerr << "   c : binary branching with conflict-directed variable ordering heuristic (default option)" << endl;
        cerr << "   q : weighted degree variable ordering heuristic" << endl;
        cerr << "   d : dichotomic branching instead of binary branching when current domain size is strictly greater than " << ToulBar2::dichotomicBranchingSize << " (default option)" << endl;
        cerr << "   e[integer] : boosting search with variable elimination of small degree (less than or equal to 3)" << " (default option)" << endl;
        cerr << "   p[integer] : preprocessing only: variable elimination of degree less than or equal to the given value" << endl;
        cerr << "   t : preprocessing only: project ternary constraints on binary constraints and apply 3-consistency" << endl;
        cerr << "   h : preprocessing only: project ternary constraints on binary constraints following a heuristic" << endl;
#ifdef BOOST
        cerr << "   m : preprocessing only: minimum degree re-ordering of variables" << endl;
#endif
        cerr << "   o : ensure optimal worst-case time complexity of DAC and EAC (can be costly in practice)" << endl;
        cerr << "   k[integer] : soft local consistency level (NC=0, AC=1, DAC=2, FDAC=3, EDAC=4)" << endl;
        cerr << "   l : limited discrepancy search" << endl;
        cerr << "   i : initial upperbound found by INCOP local search solver" << endl;
        cerr << "   z : save current problem in wcsp format" << endl;
        cerr << "   Z : debug mode (save problem at each node if verbosity option set!)" << endl;
        cerr << "   x : load a solution from a file" << endl;
        cerr << "   M[integer] : Min Sum Diffusion on preprocessing" << endl;
        cerr << "   A[integer] : enforce VAC at search nodes with depth less than a given threshold" << endl;
        cerr << "   T[integer] : threshold cost value for VAC" << endl;
        cerr << "   P[integer] : threshold cost value for VAC during the preprocessing phase" << endl;
        cerr << "   C[integer] : multiply all costs by this number" << endl;
        cerr << "   S : singleton consistency on preprocessing" << endl;
        cerr << "   V : VAC-based value ordering heuristic" << endl << endl;

        cerr << "   B[integer] : (0) DFBB, (1) BTD, (2) RDS-BTD, (3) RDS-BTD with path decomposition (default value: 0)" << endl;
		cerr << "   O[filename] : read variable elimination order from a file (if not specified, then use variable order in which variables appear in the problem file)" << endl;
        cerr << "   j[integer] : split larger clusters into a chain of smaller embedded clusters with a number of proper variables less than this number (B3j1 for pure RDS)" << endl;
        cerr << "   r[integer] : limit on maximum cluster separator size" << endl;
        cerr << "   E : merge leaf clusters with their fathers if small local treewidth (in conjunction with option e)" << endl;
        cerr << "   R[integer] : choice for root cluster" << endl;
        cerr << "   I[integer] : choice solving only a particular subtree" << endl;
#endif
        cerr << endl;
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
        	char buf[80];
        	int mode = atoi(&ch[1]);
        	if(mode > 0) {
			  ToulBar2::btdMode = mode;
			  sprintf(buf,"%s","default");
			  ToulBar2::varOrder = new char [ strlen(buf) + 1 ];
			  sprintf(ToulBar2::varOrder, "%s",buf);
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
        if ((ch = strchr(argv[i],'b'))) {if (ch[-1]==':') { ToulBar2::binaryBranching = false; } else { ToulBar2::binaryBranching = true; }}
        if ((ch = strchr(argv[i],'c'))) {if (ch[-1]==':') { ToulBar2::lastConflict = false; } else { ToulBar2::binaryBranching = true; ToulBar2::lastConflict = true; }}
        if ((ch = strchr(argv[i],'d'))) {if (ch[-1]==':') { ToulBar2::dichotomicBranching = false; } else  { ToulBar2::binaryBranching = true; ToulBar2::dichotomicBranching = true; }}
        if (strchr(argv[i],'q')) { ToulBar2::weightedDegree = true; }
        if ( (ch = strchr(argv[i],'e')) ) {
        	if (ch[-1]==':') ToulBar2::elimDegree = -1; else ToulBar2::elimDegree = 3;
        	int ndegree = atoi(&ch[1]);
        	if((ndegree > 0) && (ndegree <= 3)) ToulBar2::elimDegree = ndegree;
        }
        if ( (ch = strchr(argv[i],'p')) ) { 
        	ToulBar2::elimDegree_preprocessing = 3; 
        	int ndegree = atoi(&ch[1]); 
        	if(ndegree > 0) ToulBar2::elimDegree_preprocessing = ndegree;
        }
        if ( (ch = strchr(argv[i],'M')) ) {
        	if (!ToulBar2::vac) ToulBar2::vac = 1; 
         	ToulBar2::minsumDiffusion = 1000;
        	int nit = atoi(&ch[1]);
        	if(nit > 0) ToulBar2::minsumDiffusion = nit;
        }
        if ((ch = strchr(argv[i],'A'))) { 
		  ToulBar2::vac = 1;
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
        if (strchr(argv[i],'z')) saveproblem = true;
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
	if (ToulBar2::allSolutions && ToulBar2::btdMode >= 1) {
	  cout << "Warning! Cannot find all solutions with BTD-like search methods." << endl;
	  ToulBar2::allSolutions = false;
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
        if (localSearch(argv[1],&c)) {
            cout << "Initial upperbound: " << c << endl;
            
        } else cerr << "INCOP solver ./narycsp not found!" << endl;
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
    	if(strstr(argv[1],"binsub"))  { forceSubModular = true; randomproblem = true; sscanf(argv[1], "binsub-%d-%d-%d-%d-%d", &n, &m, &pn[0], &pn[1],&seed); narities = 2; }  
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
		if(narities == 0) cout << "Random problem incorrect, use:   bin{n}-{m}-{%}-{n. of bin ctrs}  or  tern{n}-{m}-{%}-{num bin}-{num tern}" << endl;  
    } 
    if (strstr(strext.c_str(),".pre")) ToulBar2::pedigree = new Pedigree;
    if (strstr(strext.c_str(),".bep") || strstr(argv[1],"bEpInstance")) ToulBar2::bep = new BEP;
#endif
    try {
        if(randomproblem)    solver.read_random(n,m,p,seed,forceSubModular);
        else 		         solver.read_wcsp(argv[1]);

		ToulBar2::startCpuTime = cpuTime();
        
        if (certificate) solver.read_solution("sol");
        if (saveproblem) solver.dump_wcsp("problem.wcsp");
        else if (!certificate || ToulBar2::btdMode>=2) solver.solve();
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

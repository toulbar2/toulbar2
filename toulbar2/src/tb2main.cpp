/*
 * **************** Main function ***********************************
 */

#include "toulbar2lib.hpp"
#include "tb2pedigree.hpp"
#include "tb2haplotype.hpp"
#include "tb2bep.hpp"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>

const int maxdiscrepancy = 4;
const Long maxrestarts =  10000;
const Long hbfsgloballimit = 10000;

// INCOP default command line option
const string Incop_cmd = "0 1 3 idwa 100000 cv v 0 200 1 0 0";

//* definition of path separtor depending of OS '/'  => Unix ;'\' ==> windows
#ifdef WINDOWS
#define PATH_SEP_CHR '\\'
#define PATH_DELIM ";"
#else
#define PATH_SEP_CHR '/'
#define PATH_DELIM ":"
#endif
//*definition of  windows include for command line.
#if defined(_MSC_VER)||_WIN32
# include <windows.h>
# include <tchar.h>
#else
# define TCHAR          char
# define _T(x)          x
# define _tprintf       printf
# define _tmain         main
# define _ttoi      atoi
#endif
#include "SimpleOpt.h"
#include "SimpleGlob.h"

// used for debugging purpose.
// under gdb: p ((BinaryConstraint *) constrs[13])->dump
// under gdb: p $2(constrs[13], myCout)
ostream myCout(cout.rdbuf());

#ifdef PARETOPAIR_COST
void initCosts(Cost c)
{
    if (ToulBar2::LcLevel > LC_FDAC)
    {
        cerr << "EDAC not implemented on Paretopair => force to FDAC." << endl;
        ToulBar2::LcLevel = LC_FDAC;
    }
    if (ToulBar2::vac)
    {
        cerr << "VAC not implemented on Paretopair." << endl;
        ToulBar2::vac = 0;
        ToulBar2::minsumDiffusion = 0;
    }
    if (ToulBar2::elimDegree >= 0 || ToulBar2::elimDegree_preprocessing >= 0 || ToulBar2::elimDegree_preprocessing<-1)
    {
        cerr << "Variable elimination not implemented on Paretopair." << endl;
        ToulBar2::elimDegree = -1;
        ToulBar2::elimDegree_preprocessing = -1;
    }
    if (ToulBar2::btdMode>=1)
    {
        cerr << "BTD-like methods not implemented on Paretopair." << endl;
        ToulBar2::btdMode = 0;
    }
}
#endif

/* read upper bound from given file */
string read_UB( char* ubfilename)
{
    string ubs;
    ifstream ubfile(ubfilename);
    if (ubfile.is_open())
    {
        ubfile >> ubs;
        cout << "UB read from "<< ubfilename<<"=" << ubs << endl;
        ubfile.close();

        return ubs;

    } else {
        cout << "Unable to open file:" << ubfilename << endl;
        return NULL;
    }

}

/* get ending sbsrting */
//bool getExt(char* fileName, const char* ext)
bool check_file_ext( const string fileName, const string  ext)
{
    // Finds the last persiod character of the string
    int period = fileName.find_last_of(ext[0]);
    if (period < 0) return false;
    string found_ext =  fileName.substr(period);
    if(found_ext == ext ) {
        return true;
    } else {
        //cout <<"ext " << ext << "not found in "<< fileName << endl;
        return false;
    }
}



enum { 
    OPT_verbose = 0,
    OPT_debug,
    OPT_dumpWCSP,
    OPT_HELP ,

    // file extension option
    OPT_wcsp_ext,
    OPT_wcspXML_ext,
    OPT_order_ext,
    OPT_uai_ext,
    OPT_uai_log_ext,
    OPT_evid_ext,
    OPT_map_ext,
    OPT_sol_ext,
    OPT_bep_ext,
    OPT_ub_ext,
    OPT_pre_ext,
    OPT_wcnf_ext,
    OPT_cnf_ext,
    OPT_qpbo_ext,

    // search option
    OPT_SEARCH_METHOD,
    OPT_btdRootCluster,
    OPT_btdSubTree,
    OPT_splitClusterMaxSize,
    OPT_maxSeparatorSize,
    OPT_boostingBTD,
    OPT_minProperVarSize,
    OPT_varOrder,
    OPT_problemsaved_filename,
    OPT_PARTIAL_ASSIGNMENT,
    OpT_showSolutions,
    OPT_writeSolution,
    OPT_pedigreePenalty,
    OPT_allSolutions,
    OPT_binaryBranching,
    NO_OPT_binaryBranching,
    OPT_approximateCountingBTD,
    OPT_Static_variable_ordering,
    NO_OPT_Static_variable_ordering,
    OPT_lastConflict,
    NO_OPT_lastConflict,
    OPT_dichotomicBranching,
    NO_OPT_dichotomicBranching,
    OPT_sortDomains,
    NO_OPT_sortDomains,
    OPT_weightedDegree,
    NO_OPT_weightedDegree,
    OPT_weightedTightness,
    NO_OPT_weightedTightness,
    OPT_nbDecisionVars,
    OPT_elimDegree,
    NO_OPT_elimDegree,
    OPT_elimDegree_preprocessing,
    NO_OPT_elimDegree_preprocessing,
    OPT_showSolutions,
    OPT_costfuncSeparate,
    NO_OPT_costfuncSeparate,
    OPT_nopre,

    // VAC OPTION
    OPT_minsumDiffusion,
    OPT_vac,
    NO_OPT_vac,
    OPT_costThreshold,
    OPT_costThresholdPre,
    OPT_costMultiplier,
    OPT_singletonConsistency,
    NO_OPT_singletonConsistency,
    OPT_vacValueHeuristic,
    NO_OPT_vacValueHeuristic,
    OPT_preprocessTernary,
    NO_OPT_preprocessTernary,
    OPT_preprocessFunctional,
    NO_OPT_preprocessFunctional,
    OPT_preprocessNary,
    NO_OPT_preprocessNary,
    OPT_QueueComplexity,
    NO_OPT_QueueComplexity,
    OPT_MSTDAC,
    NO_OPT_MSTDAC,
    OPT_DEE,
    NO_OPT_DEE,
    OPT_lds,
    NO_OPT_lds,
    OPT_restart,
    NO_OPT_restart,
    OPT_hbfs,
    NO_OPT_hbfs,
    OPT_open,
    OPT_localsearch,
    NO_OPT_localsearch,
    OPT_EDAC,
    OPT_ub,
    OPT_Z,
    OPT_epsilon,
    OPT_learning,
    OPT_timer,
#ifndef NDEBUG
    OPT_verifyopt,
#endif
    // MEDELESOFT OPTION
    OPT_generation=99,
    MENDEL_OPT_genotypingErrorRate=100,
    MENDEL_OPT_resolution=101,
    OPT_pedigree_by_MPE=102,
    MENDEL_OPT_EQUAL_FREQ=103,
    MENDEL_OPT_ESTIMAT_FREQ=104,
    MENDEL_OPT_ALLOCATE_FREQ=105,

    // random generator
    OPT_seed,
    OPT_random=1000
};

string getExt(string FileName)
{
    // Finds the last persiod character of the string
    int period = FileName.find_last_of(".");
    // I use  + 1 because I don't really need the to include the period
    string ext = FileName.substr(period + 1);
    return ext;
}

CSimpleOpt::SOption g_rgOptions[] =
{
    { OPT_HELP,                     (char*) "-h",                    SO_NONE         }, // boolean help
    { OPT_HELP,                     (char*) "-?",                    SO_NONE         }, // boolean help
    { OPT_HELP,                		(char*) "-help",          			SO_NONE     	}, // boolean help
    { OPT_HELP,                		(char*) "--help",         			SO_NONE     	}, // boolean help
    { OPT_verbose,  			 (char*) "-v", 				SO_OPT			}, // verbose level
    { OPT_debug,  				 (char*) "-Z", 				SO_OPT			},  // debug level
    { OPT_dumpWCSP,  			 (char*) "-z", 				SO_OPT			}, // dump wcsp

    // file extension
    { OPT_wcsp_ext,      	          (char*) "--wcsp_ext",                                 SO_REQ_SEP      },
    { OPT_wcspXML_ext,                 (char*) "--wcspXML_ext",                              SO_REQ_SEP      },
    { OPT_order_ext, 	                 (char*) "--order_ext",                                SO_REQ_SEP      },
    { OPT_uai_ext,          	         (char*) "--uai_ext",                                  SO_REQ_SEP      },
    { OPT_uai_log_ext,                       (char*) "--uai_log_ext",                                  SO_REQ_SEP      },
    { OPT_evid_ext,                 	 (char*) "--evid_ext",                                 SO_REQ_SEP      },
    { OPT_map_ext,                     (char*) "--map_ext",                                  SO_REQ_SEP      },
    { OPT_sol_ext,                     (char*) "--sol_ext",                                  SO_REQ_SEP      },
    { OPT_bep_ext,                     (char*) "--bep_ext",                                  SO_REQ_SEP      },
    { OPT_ub_ext,                      (char*) "--ub_ext",                                   SO_REQ_SEP      },
    { OPT_pre_ext,                     (char*) "--pre_ext",                                  SO_REQ_SEP      },
    { OPT_wcnf_ext,                 (char*) "--wcnf_ext",                              SO_REQ_SEP      },
    { OPT_cnf_ext,                 (char*) "--cnf_ext",                              SO_REQ_SEP      },
    { OPT_qpbo_ext,                 (char*) "--qpbo_ext",                              SO_REQ_SEP      },


    { OPT_SEARCH_METHOD,       		(char*) "-B",             			SO_REQ_SEP  	}, // -B [0,1,2] search method
    { OPT_SEARCH_METHOD,       		(char*) "--search",       			SO_REQ_SEP 	},
    { OPT_btdRootCluster,      		(char*) "-R",             			SO_REQ_SEP  	}, // root cluster used in BTD
    { OPT_btdRootCluster,       		(char*) "--RootCluster",  			SO_REQ_CMB  	},
    { OPT_btdSubTree,          		(char*) "-I",             			SO_REQ_SEP  	}, // btd sub tree
    { OPT_splitClusterMaxSize, 		(char*) "-j",             			SO_REQ_SEP  	},
    { OPT_maxSeparatorSize,     		(char*) "-r",             			SO_REQ_SEP  	},
    { OPT_maxSeparatorSize,     		(char*) "--maxSepSize",   			SO_REQ_CMB  	},

    { OPT_minProperVarSize,     		(char*) "-X",             			SO_REQ_SEP 	},
    { OPT_PARTIAL_ASSIGNMENT,  		(char*) "-x",             			SO_OPT 	 	},
    { OPT_boostingBTD,          		(char*) "-E",             			SO_NONE    	},
    { OPT_varOrder,		    		(char*) "-O",             			SO_REQ_SEP 	}, // filename of variable order
    { OPT_problemsaved_filename,        (char*) "--save",                       SO_REQ_SEP  }, // filename of saved problem
    { OPT_showSolutions,         		(char*) "-s",             			SO_NONE    	},//print solution found
    { OPT_showSolutions,         		(char*) "--show",          			SO_NONE    	},//print solution found
    { OPT_writeSolution,        		(char*) "-w",       			  	SO_OPT  	}, //  depending of value write last solution found in filename "sol" in different format

    { OPT_pedigreePenalty,       		(char*) "-u",       			   	SO_REQ_SEP  	}, // int ..
    { OPT_allSolutions,	  		(char*) "-a",       			  	SO_NONE		}, // counting option ...print solution found
    { OPT_approximateCountingBTD,	 	(char*) "-D",           			SO_NONE    	}, //approximate counting
    { OPT_binaryBranching,			(char*) "-b",       			  	SO_NONE    	},
    { OPT_binaryBranching,			(char*) "-binaryBranching",         		SO_NONE    	},
    { NO_OPT_binaryBranching,		(char*) "-b:",             		   	SO_NONE    	},
    { NO_OPT_binaryBranching,		(char*) "-no--b",             		SO_NONE    	},
    { NO_OPT_binaryBranching,		(char*) "-no--binaryBranching",             	SO_NONE    	},
    { OPT_Static_variable_ordering,		(char*) "-svo", 			SO_NONE 		},
    { NO_OPT_Static_variable_ordering, 		(char*) "-svo:", 			SO_NONE 	},
    { OPT_lastConflict, 			(char*) "-c", 				SO_NONE 	},
    { NO_OPT_lastConflict, 			(char*) "-c:", 				SO_NONE 	},
    { NO_OPT_lastConflict, 			(char*) "-no--c", 				SO_NONE 	},
    { NO_OPT_lastConflict, 			(char*) "--lastConflict--off", 		SO_NONE 	},
    { OPT_dichotomicBranching,		(char*) "-d", 				SO_OPT 	},
    { NO_OPT_dichotomicBranching,		(char*) "-d:", 				SO_NONE 	},
    { OPT_sortDomains,		(char*) "-sortd", 				SO_NONE 	},
    { NO_OPT_sortDomains,		(char*) "-sortd:", 				SO_NONE 	},
    { OPT_weightedDegree,			(char*) "-q", 				SO_OPT 	},
    { NO_OPT_weightedDegree, 		(char*) "-q:", 				SO_NONE 	},
    { OPT_weightedTightness,	    (char*) "-m",       			  	SO_OPT    	},
    { NO_OPT_weightedTightness,		(char*) "-m:",             		   	SO_NONE    	},
    { OPT_nbDecisionVars,			(char*) "-var", 				SO_REQ_SEP		},

    { OPT_elimDegree,			(char*) "-e", 				SO_OPT			},
    { NO_OPT_elimDegree,		 	(char*) "-e:", 				SO_NONE 		},
    { OPT_elimDegree_preprocessing, 	(char*) "-p", 				SO_OPT	 		},
    { NO_OPT_elimDegree_preprocessing,	(char*) "-p:", 				SO_NONE 		},
    { OPT_costfuncSeparate,			(char*) "-dec", 			SO_NONE 		},
    { NO_OPT_costfuncSeparate,		(char*) "-dec:",			SO_NONE			},
    { OPT_nopre,			(char*) "-nopre", 			SO_NONE 		},
    // vac option
    { OPT_vac,				(char*) "-A", 				SO_OPT			},
    { NO_OPT_vac,				(char*) "-A:", 				SO_NONE			},
    { OPT_vacValueHeuristic,		(char*) "-V", 				SO_NONE 		},
    { OPT_costThreshold,			(char*) "-T", 				SO_REQ_SEP		},
    { OPT_costThresholdPre,	 		(char*) "-P", 				SO_REQ_SEP		},
    { OPT_costMultiplier,	 		(char*) "-C", 				SO_REQ_SEP		},

    //preprocessing
    { OPT_minsumDiffusion,	 		(char*) "-M", 				SO_REQ_SEP		},
    { OPT_singletonConsistency,		(char*) "-S", 				SO_NONE			},
    { OPT_preprocessTernary,		(char*) "-t", 				SO_OPT		},
    { NO_OPT_preprocessTernary,		(char*) "-t:", 				SO_NONE		},
    { OPT_preprocessFunctional,		(char*) "-f", 				SO_OPT		},
    { NO_OPT_preprocessFunctional,	(char*) "-f:", 		 	    SO_NONE		},
    { OPT_preprocessNary,			(char*) "-n", 				SO_OPT		},
    { NO_OPT_preprocessNary,		(char*) "-n:", 				SO_NONE			},

    { OPT_QueueComplexity,			(char*) "-o", 				SO_NONE		},
    { OPT_MSTDAC,			        (char*) "-mst", 			SO_NONE		},
    { NO_OPT_MSTDAC,			(char*) "-mst:", 			SO_NONE		},
    { OPT_DEE,	     			(char*) "-dee", 			SO_OPT 		},
    { NO_OPT_DEE,				(char*) "-dee:",			SO_OPT			},
    { OPT_lds,				(char*) "-l", 				SO_OPT			},
    { NO_OPT_lds,	 			(char*) "-l:", 				SO_NONE			},
    { OPT_restart,	 			(char*) "-L", 				SO_OPT			},
    { NO_OPT_restart,			(char*) "-L:", 				SO_NONE			},
    { OPT_hbfs,              (char*) "-hbfs",               SO_OPT          },
    { OPT_hbfs,              (char*) "-bfs",               SO_OPT          },
    { NO_OPT_hbfs,           (char*) "-hbfs:",              SO_NONE         },
    { NO_OPT_hbfs,           (char*) "-bfs:",              SO_NONE         },
    { OPT_open,              (char*) "-open",               SO_REQ_SEP          },
    { OPT_localsearch,			(char*) "-i", 				SO_OPT			},// incop option default or string for narycsp argument
    { OPT_EDAC,				(char*) "-k", 				SO_REQ_SEP		},
    { OPT_ub, 	 			(char*) "-ub", 				SO_REQ_SEP		}, // init upper bound in cli
    // MEDELSOFT
    { OPT_generation,			(char*) "-g", 				SO_NONE			}, //sort pedigree by increasing generation number and if equal by increasing individual number
    //	{ OPT_pedigree_by_MPE,  		(char*) "-y", 				SO_OPT			}, // bayesian flag
    { MENDEL_OPT_genotypingErrorRate,	(char*) "-genoError", 			SO_REQ_SEP		},
    { MENDEL_OPT_resolution,		(char*) "-precision", 			SO_REQ_SEP		},
    { OPT_pedigree_by_MPE,  		(char*) "-bayes", 				SO_NONE 		}, // bayesian flag
    { MENDEL_OPT_EQUAL_FREQ,  		(char*) "-pequal", 				SO_NONE			}, // allocate equal frequencies to all allele
    { MENDEL_OPT_ESTIMAT_FREQ,		(char*) "-probdata", 			SO_NONE			}, //  probs depending on the frequencies found in the current pedigree problem
    { MENDEL_OPT_ALLOCATE_FREQ,		(char*) "-problist", 			SO_MULTI		}, // read probability distribution from command line

    { OPT_Z,  				(char*) "-logz", 				SO_NONE			},  // compute log partition function (log Z)
    { OPT_epsilon,				(char*) "-epsilon", 			SO_REQ_SEP		}, // approximation parameter for computing Z

    { OPT_learning,                         	(char*) "-learning",                    SO_NONE }, // pseudoboolean learning during search
    #ifndef NDEBUG
    { OPT_verifyopt,                (char*) "-opt",             SO_NONE }, // for debugging purposes, checks the given optimal solution (problem.sol) is not pruned during search
    #endif
    { OPT_timer,              (char*) "-timer",             SO_REQ_SEP      }, // CPU timer

    // random generator
    { OPT_seed,			         (char*) "-seed", 				SO_REQ_SEP},
    { OPT_random, 	 			 (char*) "-random", 			SO_REQ_SEP }, // init upper bound in cli


    SO_END_OF_OPTIONS
};

void ShowFiles(int argc, TCHAR ** argv)
{
    // glob files to catch expand wildcards
    CSimpleGlob glob(SG_GLOB_NODOT|SG_GLOB_NOCHECK);
    if (SG_SUCCESS != glob.Add(argc, argv))
    {
        _tprintf(_T("Error while globbing files\n"));
        return;
    }

    for (int n = 0; n < glob.FileCount(); ++n)
    {
        _tprintf(_T("file %2d: '%s'\n"), n, glob.File(n));
    }
}

static const TCHAR * GetLastErrorText( int a_nError)
{
    switch (a_nError)
    {
    case SO_SUCCESS:
        return _T("Success");
    case SO_OPT_INVALID:
        return _T("Unrecognized option");
    case SO_OPT_MULTIPLE:
        return _T("Option matched multiple strings");
    case SO_ARG_INVALID:
        return _T("Option does not accept argument");
    case SO_ARG_INVALID_TYPE:
        return _T("Invalid argument format");
    case SO_ARG_MISSING:
        return _T("\n Required argument is missing");
    case SO_ARG_INVALID_DATA:
        return _T("Invalid argument data");
    default:
        return _T("Unknown error");
    }
}

// processing of  multi argument -problist option 
// -problist #arguments arg1 arg2 arg2 ...#argument

static void Pedi_Args( CSimpleOpt &    args, int nMultiArgs)
{
    TCHAR ** rgpszArg = NULL;
    int nMultiBackup = nMultiArgs ;

    // get the number of arguments if necessary
    if (nMultiArgs == -1) {
        // first arg is a count of how many we have
        rgpszArg = args.MultiArg(1);

        if (!rgpszArg) {
            _tprintf(
                    _T("%s: '%s' (missing argument or please use --help to get command line help)\n \n"),
                    GetLastErrorText(args.LastError()), args.OptionText());
            return;

        }


        nMultiArgs = _ttoi(rgpszArg[0]);

        cout << endl;
        cout << nMultiArgs << "Alleles probability will be read from command line " << endl;
    }

    // get the arguments to follow
    rgpszArg = args.MultiArg(nMultiArgs);

    if(nMultiBackup < 0 ) {
        for (int n = 0; n < nMultiArgs; ++n) {
            if(ToulBar2::verbose >= 0) _tprintf(_T("MultiArg %d: %s\n"), n, rgpszArg[n]);
            float f2;
            sscanf(rgpszArg[n],"%f",&f2) ;
            ToulBar2::allelefreqdistrib.push_back(f2);
        }
        if(ToulBar2::verbose >= 0)
        {
            for(unsigned int x=0; x<=sizeof(ToulBar2::allelefreqdistrib); x++) cout <<"Allele prob used "<< x << "=" <<  ToulBar2::allelefreqdistrib[x] << endl;
        }
        _tprintf(_T("%s: expecting %d args\n"), args.OptionText(), nMultiArgs);

    }

    if (!rgpszArg) {
        _tprintf(
                _T("%s: '%s' (use --help to get command line help)\n"),
                GetLastErrorText(args.LastError()), args.OptionText());
        exit(-1);
    }




}


/* return current binary path extract from or argv[0] or from the env var $path the env var $path */

char* find_bindir(const char* bin_name, char* buffer, size_t buflen)
{
    struct stat st;
    char* path, * tok;
    if (!stat(bin_name, &st))
    {
        char* end = (char*) strrchr(bin_name, PATH_SEP_CHR);
        static char bin_path[512];
        if (end)
        {
            *end=0;
            strncpy(buffer, bin_name, buflen);
            sprintf(bin_path,"%s%c", buffer,PATH_SEP_CHR);
        }
        else
        {
            strcpy(buffer, ".");
            //path separator added to the path value
            sprintf(bin_path,"%s%c", buffer,PATH_SEP_CHR);
        }
        return(bin_path);
    }
    path = strdup(getenv("PATH"));
    tok = strtok(path, PATH_DELIM);
    while (tok)
    {
        snprintf(buffer, buflen, "%s%c%s", tok, PATH_SEP_CHR, bin_name);
        if (!stat(buffer, &st))
        {
            static char bin_path[512];
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

//  current unused option letters: 	f F G H J K n N Q U W Y
void help_msg(char *toulbar2filename)
{
    cout << "*************************" << endl;
    cout << "* ToulBar2 Help Message *" << endl;
    cout << "*************************" << endl;
    cout << endl;
    cout << "Command line is:" << endl;
    cout << toulbar2filename << " problem_filename [options]" << endl;
    cout << endl;
#ifndef MENDELSOFT
    cout << "Available problem formats (specified by the filename extension) are:" << endl;
    cout << "   *.wcsp : Weighted CSP format (see SoftCSP web site)" << endl;
    cout << "   *.wcnf : Weighted Partial Max-SAT format (see Max-SAT Evaluation)" << endl;
    cout << "   *.cnf : (Max-)SAT format" << endl;
    cout << "   *.qpbo : quadratic pseudo-Boolean optimization (unconstrained quadratic programming) format" << endl;
#ifdef XMLFLAG
    cout << "   *.xml : CSP and weighted CSP in XML format XCSP 2.1";
#ifdef MAXCSP
    cout << " (Max-CSP only)";
#endif
    cout << endl;
#endif
    cout << "   *.uai : Bayesian network and Markov Random Field format (see UAI'08 Evaluation) followed by an optional evidence filename (performs MPE task, see -logz for PR task)" << endl;
    cout << "   *.LG : Bayesian network and Markov Random Field format using logarithms instead of probabilities" << endl;
    cout << "   *.pre : pedigree format (see doc/MendelSoft.txt for Mendelian error correction)" << endl;
    cout << "   *.pre *.map : pedigree and genetic map formats (see doc/HaplotypeHalfSib.txt for haplotype reconstruction in half-sib families)" << endl;
    cout << "   *.bep  : satellite scheduling format (CHOCO benchmark)" << endl << endl;
    cout << "   *.order  : variable elimination order" << endl;
    cout << "   *.sol  : solution/certificate for the problem" << endl << endl;
    cout << "Warning! a New file extension can be enforced using --foo_ext=\".myext\" ex: --wcsp_ext='.test' --sol_ext='.sol2'  " << endl << endl;
#endif
    cout << "Available options are (use symbol \":\" after an option to remove a default option):" << endl;
    cout << "   -help : shows this help message" << endl;
    cout << "   -ub=[integer] : initial problem upperbound (default value is " << MAX_COST << ")" << endl;
    cout << "   -v=[integer] : verbosity level" << endl;
    cout << "   -s : shows each solution found" << endl;
#ifndef MENDELSOFT
    cout << "   -w=[filename] : writes last solution found in filename (or \"sol\" if no parameter is given)" << endl;
    cout << "   -precision=[integer] : probability/real precision is a conversion factor (a power of ten) for representing fixed point numbers (default value is " << ToulBar2::resolution << ")" << endl;
#else
    cout << "   -w=[mode] : writes last solution found" << endl;
    cout << "               mode=0: saves pedigree with erroneous genotypings removed" << endl;
    cout << "               mode=1: saves pedigree with erroneous genotypings corrected" << endl;
    cout << "               mode=2: saves pedigree with erroneous genotypings corrected and missing genotypes of informative individuals inferred" << endl;
    cout << "   --save=[filename] : saves pedigree in filename (or \"pedigree_corrected.pre\" if no parameter is given)" << endl;
    cout << "   -g : sorts pedigree by increasing generation number and if equal by increasing individual number" << endl;
    cout << "   -u=[integer] : adds a penalty weight (must use option y also) on genotyped individuals depending on the number of their genotyped children in order to penalize genotyping removals if the number of genotyped children is strictly greater than a given threshold" << endl;


    cout << "   -bayes : pedigree solved by Bayesian MPE . the following option can be tune" << endl;
    cout << "               -genoError [real]<=> genotyping Error Rate is a prior uniform probability of genotyping errors (default value is " << ToulBar2::errorg << ")" << endl;
    cout << "               -precision [int]<=> probability Precision is a conversion factor (a power of ten) for representing fixed point numbers (default value is " << ToulBar2::resolution << ")" << endl;

    cout << "         the command line possibly followed by three exclusive options:" << endl;
    cout << "           -pequal	  : uniform allele probability distribution (default mode) " << endl;
    cout << "           -probd    : allele probability distribution read from pedigree data" << endl;
    cout << "           -problist [nbre of prob] p1 p2 p3... : allele probability distribution given explicitely in the command line" << endl << endl;
#endif
#ifndef MENDELSOFT
#ifdef LINUX
    cout << "   -timer=[integer] : CPU time limit in seconds" << endl;
#endif
    cout << "   -var=[integer] : searches by branching only on the first -the given value- decision variables, assuming the remaining variables are intermediate variables completely assigned by the decision variables (use a zero if all variables are decision variables) (default value is " << ToulBar2::nbDecisionVars << ")" << endl;
    cout << "   -b : searches using binary branching always instead of binary branching for interval domains and n-ary branching for enumerated domains";
    if (ToulBar2::binaryBranching) cout << " (default option)";
    cout << endl;
    cout << "   -svo : searches using a static variable ordering heuristic (same order as DAC)";
    if (ToulBar2::Static_variable_ordering) cout << " (default option)";
    cout << endl;
    cout << "   -c : searches using binary branching with last conflict backjumping variable ordering heuristic";
    if (ToulBar2::lastConflict) cout << " (default option)";
    cout << endl;
    cout << "   -q=[integer] : weighted degree variable ordering heuristic if the number of cost functions is less than the given value (default value is " << ToulBar2::weightedDegree << ")" << endl;
    cout << "   -m=[integer] : variable ordering heuristic based on mean (m=1) or median (m=2) costs (in conjunction with weighted degree heuristic -q) (default value is " << ToulBar2::weightedTightness << ")" << endl;
    cout << "   -d=[integer] : searches using dichotomic branching (d=1 splitting in the middle of domain range, d=2 splitting in the middle of sorted unary costs) instead of binary branching when current domain size is strictly greater than " << ToulBar2::dichotomicBranchingSize << " (default value is " << ToulBar2::dichotomicBranching << ")" << endl;
    cout << "   -sortd : sorts domains based on increasing unary costs (warning! works only for binary WCSPs)";
    if (ToulBar2::sortDomains) cout << " (default option)";
    cout << endl;
    cout << "   -e=[integer] : boosting search with variable elimination of small degree (less than or equal to 3) (default value is " << ToulBar2::elimDegree << ")" << endl;
    cout << "   -p=[integer] : preprocessing only: general variable elimination of degree less than or equal to the given value (default value is " << ToulBar2::elimDegree_preprocessing << ")" << endl;
    cout << "   -t=[integer] : preprocessing only: simulates restricted path consistency by adding ternary cost functions on triangles of binary cost functions within a given maximum space limit (in MB)";
    if (ToulBar2::preprocessTernaryRPC) cout << " (" << ToulBar2::preprocessTernaryRPC << " MB)";
    cout << endl;
    cout << "   -f=[integer] : preprocessing only: variable elimination of functional (f=1) (resp. bijective (f=2)) variables (default value is " << ToulBar2::preprocessFunctional << ")" << endl;
    cout << "   -dec : preprocessing only: pairwise decomposition of cost functions with arity >=3 into smaller arity cost functions";
    if (ToulBar2::costfuncSeparate) cout << " (default option)";
    cout << endl;
    cout << "   -n=[integer] : preprocessing only: projects n-ary cost functions on all binary cost functions if n is lower than the given value (default value is " << ToulBar2::preprocessNary << ")" << endl;
#ifdef BOOST
    cout << "   -mst : maximum spanning tree DAC ordering";
    if (ToulBar2::MSTDAC) cout << " (default option)";
    cout << endl;
#endif
    cout << "   -nopre : removes all preprocessing options (equivalent to -e: -p: -t: -f: -dec: -n: -mst: -dee:)" << endl;
    cout << "   -o : ensures optimal worst-case time complexity of DAC and EAC (can be slower in practice)";
    if (ToulBar2::QueueComplexity) cout << " (default option)";
    cout << endl;
    cout << "   -k=[integer] : soft local consistency level (NC with Strong NIC for global cost functions=0, (G)AC=1, D(G)AC=2, FD(G)AC=3, (weak) ED(G)AC=4) (default value is " << ToulBar2::LcLevel << ")" << endl;
    cout << "   -dee=[integer] : restricted dead-end elimination (value pruning by dominance rule from EAC value (dee>=1 and dee<=3)) and soft neighborhood substitutability (in preprocessing (dee=2 or dee=4) or during search (dee=3)) (default value is " << ToulBar2::DEE << ")" << endl;
    cout << "   -l=[integer] : limited discrepancy search, use a negative value to stop the search after the given absolute number of discrepancies has been explored (discrepancy bound = " << maxdiscrepancy << " by default)";
    if (ToulBar2::lds) cout << " (default option)";
    cout << endl;
    cout << "   -L=[integer] : randomized (quasi-random variable ordering) search with restart (maximum number of nodes = " << maxrestarts << " by default)";
    if (ToulBar2::restart>=0) cout << " (default option)";
    cout << endl;
    cout << "   -i=[\"string\"] : initial upperbound found by INCOP local search solver." << endl;
    cout << "       string parameter is optional, using \"" << Incop_cmd << "\" by default with the following meaning:" << endl;
    cout << "       stoppinglowerbound randomseed nbiterations method nbmoves neighborhoodchoice neighborhoodchoice2 minnbneighbors maxnbneighbors neighborhoodchoice3 autotuning tracemode"<< endl;

    cout << "   -z=[filename] : saves problem in wcsp format in filename (or \"problem.wcsp\"  if no parameter is given)" << endl;
    cout << "                   writes also the  graphviz dot file  and the degree distribution of the input problem" << endl;
    cout << "   -z=[integer] : 1: saves original instance (by default), 2: saves after preprocessing" << endl;
    cout << "   -Z=[integer] : debug mode (save problem at each node if verbosity option -v=num >= 1 and -Z=num >=3)" << endl;
#ifndef NDEBUG
    cout << "   -opt filename.sol : checks a given optimal solution (given as input filename with \".sol\" extension) is never pruned by propagation (works only if compiled with debug)" << endl;
#endif
    cout << "   -x=[(,i=a)*] : assigns variable of index i to value a (multiple assignments are separated by a comma and no space) (without any argument, a complete assignment -- used as initial upper bound and as value heuristic -- read from default file \"sol\" or given as input filename with \".sol\" extension)" << endl << endl;
    cout << "   -M=[integer] : preprocessing only: Min Sum Diffusion algorithm (default number of iterations is " << ToulBar2::minsumDiffusion << ")" << endl;
    cout << "   -A=[integer] : enforces VAC at each search node with a search depth less than a given value (default value is " << ToulBar2::vac << ")" << endl;
    cout << "   -T=[integer] : threshold cost value for VAC (default value is " << ToulBar2::costThreshold << ")" << endl;
    cout << "   -P=[integer] : threshold cost value for VAC during the preprocessing phase (default value is " << ToulBar2::costThresholdPre << ")" << endl;
    cout << "   -C=[integer] : multiplies all costs by this number (default value is " << ToulBar2::costMultiplier << ")" << endl;
    cout << "   -S : preprocessing only: performs singleton consistency (only in conjunction with option \"-A\")";
    if (ToulBar2::singletonConsistency) cout << " (default option)";
    cout << endl;
    cout << "   -V : VAC-based value ordering heuristic";
    if (ToulBar2::vacValueHeuristic) cout << " (default option)";
    cout << endl << endl;

    cout << "   -B=[integer] : (0) DFBB, (1) BTD, (2) RDS-BTD, (3) RDS-BTD with path decomposition instead of tree decomposition (default value is " << ToulBar2::btdMode << ")" << endl;
    cout << "   -O=[filename] : reads a variable elimination order from a file in order to build a tree decomposition (if BTD-like and/or variable elimination methods are used) and also a compatible DAC ordering" << endl;
#ifdef BOOST
    cout << "   -O=[negative integer] : build a tree decomposition (if BTD-like and/or variable elimination methods are used) and also a compatible DAC ordering using" << endl;
    cout << "                           (-1) maximum cardinality search ordering, (-2) minimum degree ordering, (-3) minimum fill-in ordering," << endl;
    cout << "                           (-4) maximum spanning tree ordering (see -mst), (-5) reverse Cuthill-Mckee ordering, (-6) approximate minimum degree ordering" << endl;
#endif
    cout << "                  (if not specified, then use the variable order in which variables appear in the problem file)" << endl;
    cout << "   -j=[integer] : splits large clusters into a chain of smaller embedded clusters with a number of proper variables less than this number" << endl;
    cout << "                (use options \"-B=3 -j=1 -svo -k=1\" for pure RDS, use value 0 for no splitting) (default value is " << ToulBar2::splitClusterMaxSize << ")" << endl;
    cout << "   -r=[integer] : limit on maximum cluster separator size (merge cluster with its father otherwise, use a negative value for no limit) (default value is " << ToulBar2::maxSeparatorSize << ")" << endl;
    cout << "   -X=[integer] : limit on minimum number of proper variables in a cluster (merge cluster with its father otherwise, use a zero for no limit) (default value is " << ToulBar2::minProperVarSize << ")" << endl;
    cout << "   -E : merges leaf clusters with their fathers if small local treewidth (in conjunction with option \"-e\")";
    if (ToulBar2::boostingBTD) cout << " (default option)";
    cout << endl;
    cout << "   -R=[integer] : choice for a specific root cluster number" << endl;
    cout << "   -I=[integer] : choice for solving only a particular rooted cluster subtree (with RDS-BTD only)" << endl << endl;
    cout << "   -a : finds all solutions (or count the number of zero-cost satisfiable solutions in conjunction with BTD)";
    if (ToulBar2::allSolutions) cout << " (default option)";
    cout << endl;
    cout << "   -D : approximate satisfiable solution count with BTD";
    if (ToulBar2::approximateCountingBTD) cout << " (default option)";
    cout << endl;
    cout << "   -logz : computes log of probability of evidence (i.e. log partition function or log(Z) or PR task) for graphical models only (problem file extension .uai)" << endl;
    cout << "   -epsilon=[float] : approximation factor for computing the partition function (default value is " << Exp(-ToulBar2::logepsilon) << ")" << endl;
    cout << endl;
    cout << "   -hbfs=[integer] : hybrid best-first search, restarting from the root after a given number of backtracks (default value is " << hbfsgloballimit << ")" << endl;
    cout << "   -open=[integer] : hybrid best-first search limit on the number of open nodes (default value is " << ToulBar2::hbfsOpenNodeLimit << ")" << endl;
    cout << "---------------------------" << endl;
    cout << "Alternatively one can call the random problem generator with the following options: " << endl;
    cout << endl;
    cout << "   -random=[bench profile]  : bench profile must be specified as follow :" << endl;
    cout << "                         n and m are respectively the number of variable and the maximum domain size  of the random problem." << endl;
    cout << "			"<< endl;
    cout << "       bin-{n}-{m}-{p1}-{p2}-{seed}       :p1 is the tightness in percentage %" << endl;
    cout << "                                          :p2 is the num of binary cost functions to include" << endl;
    cout << "                                          :the seed parameter is optional" << endl;

    cout << "   or:                                                                               " << endl;
    cout << "       binsub-{n}-{m}-{p1}-{p2}-{p3}-{seed} binary random & submodular cost functions" << endl;
    cout << "                                          p1 is the tightness in percentage % of random cost functions" << endl;
    cout << "                                          p2 is the num of binary cost functions to include" << endl;
    cout << "                                          p3 is the percentage % of submodular cost functions among p2 cost functions" << endl;
    cout << "                                           (plus 10 permutations of two randomly-chosen values for each domain)" << endl;
    cout << " or:                                                                               " << endl;
    cout << "      tern-{n}-{m}-{p1}-{p2}-{p3}-{seed}  p3 is the num of ternary cost functions" << endl;
    cout << " or:                                                                               " << endl;
    cout << "      nary-{n}-{m}-{p1}-{p2}-{p3}...{pn}-{seed}  pn is the num of n-ary cost functions" << endl;
    cout << "---------------------------" << endl;
    cout << "			"<< endl;

    cout << endl;
#endif
}

int _tmain(int argc, TCHAR * argv[])
{
    tb2init();

    setlocale( LC_ALL, "C" );
    bool certificate = false;
    char *certificateFilename = NULL;
    char *certificateString = NULL;
    char buf [512];
    char* CurrentBinaryPath = find_bindir(argv[0], buf, 512); // current binary path search
    Cost ub = MAX_COST;
    int timeout = 0;

    // Configuration for MaxSAT Evaluation
    //	ToulBar2::maxsateval = true;
    //	ToulBar2::verbose = -1;
    //	ToulBar2::binaryBranching = false;
    //	ToulBar2::lds = 1;

    // Configuration for UAI Evaluation
    //	ToulBar2::uaieval = true;
    //  ToulBar2::verbose = 0;
    //  ToulBar2::lds = 1;
    //  ToulBar2::incop_cmd = "0 1 3 idwa 100000 cv v 0 200 1 0 0";

    char *random_desc = NULL ; // benchmark description set from command line;

    //default file extension : can be enforced using --foo_ext option in command line

    std::map<std::string, string> file_extension_map;
    file_extension_map["wcsp_ext"]=".wcsp";
    file_extension_map["wcspXML_ext"]=".xml";
    file_extension_map["order_ext"]=".order";
    file_extension_map["ub_ext"]=".ub";
    file_extension_map["sol_ext"]=".sol";
    file_extension_map["uai_ext"]=".uai";
    file_extension_map["uai_log_ext"]=".LG";
    file_extension_map["evid_ext"]=".evid";
    file_extension_map["bep_ext"]=".bep";
    file_extension_map["pre_ext"]=".pre";
    file_extension_map["map_ext"]=".map";
    file_extension_map["wcnf_ext"]=".wcnf";
    file_extension_map["cnf_ext"]=".cnf";
    file_extension_map["qpbo_ext"]=".qpbo";


    assert(cout << "Warning! toulbar2 was compiled in debug mode and it can be very slow..." << endl);
    if (ToulBar2::verbose >= 0) cout << "c " << CurrentBinaryPath<<"toulbar2"<<"  version : " << ToulBar2::version << ", copyright (c) INRA 2015"<<endl;

    // --------------------------simple opt ----------------------

    // declare our options parser, pass in the arguments from main
    // as well as our array of valid options.
    CSimpleOpt args(argc, argv, g_rgOptions);


    while (args.Next())
    {

        if (args.LastError() == SO_SUCCESS)
        {

            //if (check_file_ext(to_string(args.OptionText()),"_ext") )
            if (strstr(args.OptionText(),"_ext") )
            {
                string force_extension = args.OptionText();
                force_extension.erase(force_extension.begin(),force_extension.begin()+2);
                if (ToulBar2::debug) cout << "Old extension " << file_extension_map[force_extension] << " --> ";
                //			cout << " extension " << force_extension << " forced with " << args.OptionArg() << endl;
                file_extension_map[force_extension]=args.OptionArg();

                if (ToulBar2::debug) cout <<" New extension " << file_extension_map[force_extension] << endl;

            }


            //  search algo
            if (args.OptionId() == OPT_SEARCH_METHOD)
            {
                int mode = 0;
                if(args.OptionArg() != NULL ) mode = atoi(args.OptionArg());
                if (mode >= 0)
                {
                    ToulBar2::btdMode = mode;
                } else ToulBar2::btdMode = 0;
                if (ToulBar2::debug) cout << "Search Method used =  " << mode << endl;
            }
            // BTD root cluster
            if (args.OptionId() == OPT_btdRootCluster)
            {
                int root = atoi(args.OptionArg());
                if (root > 0) ToulBar2::btdRootCluster = root;
            }
            // btd SubTree initialisation sub cluster

            if (args.OptionId() == OPT_btdSubTree)
            {
                int subcluster =atoi(args.OptionArg());
                if (subcluster >= 1) ToulBar2::btdSubTree = subcluster;
            }
            //cluster Max size
            if (args.OptionId() == OPT_splitClusterMaxSize)
            {
                int cmaxsize = atoi(args.OptionArg());
                if (cmaxsize >= 1) ToulBar2::splitClusterMaxSize = cmaxsize;
            }

            // E : merge leaf clusters with their fathers if small local treewidth
            if (args.OptionId() == OPT_boostingBTD)
            {
                ToulBar2::boostingBTD = true;
            }

            //   -r=[integer] : limit on maximum cluster separator size (merge cluster with its father otherwise, use a negative value for no limit)
            if (args.OptionId() == OPT_maxSeparatorSize)
            {
                int sepmaxsize = atoi(args.OptionArg());
                if (sepmaxsize >= -1) ToulBar2::maxSeparatorSize = sepmaxsize;
            }
            // minimal number of Proper variable included into each cluster of BTD
            if (args.OptionId() == OPT_minProperVarSize)
            {
                int minpvarsize = atoi(args.OptionArg());
                if (minpvarsize >= 0) ToulBar2::minProperVarSize = minpvarsize;
            }
            // -help print command line HELP
            if (args.OptionId() == OPT_HELP) {
                //	ShowUsage();
                help_msg(argv[0]);
                return 0;
            }


            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            // filename containing variable order

            if (args.OptionId() == OPT_varOrder)
            {
                int varElimOrder = atoi(args.OptionArg());
                if (varElimOrder>=0) {
                    char buf[512];
                    sprintf(buf,"%s",args.OptionArg());
                    //				if (ToulBar2::varOrder) delete [] ToulBar2::varOrder;
                    ToulBar2::varOrder = new char [ strlen(buf) + 1 ];
                    sprintf(ToulBar2::varOrder, "%s",buf);
                    if (ToulBar2::debug) cout << "variable order read from file " << args.OptionArg() << endl;
                } else {
                    ToulBar2::varOrder = reinterpret_cast<char *>(-varElimOrder);
                }
            }

            if (args.OptionId() == OPT_problemsaved_filename)
            {
                char buf[512];
                sprintf(buf,"%s",args.OptionArg());
                ToulBar2::problemsaved_filename = to_string(buf);
                //                if (!ToulBar2::dumpWCSP) ToulBar2::dumpWCSP = 1;
                if (ToulBar2::debug) cout << "saved problem into file " << ToulBar2::problemsaved_filename << endl;
            }

            // filename of solution
            if (args.OptionId() == OPT_PARTIAL_ASSIGNMENT) {
                if (args.OptionArg()!= NULL)	{
                    certificate = true;
                    certificateString = args.OptionArg() ;
                    if (ToulBar2::debug) cout << "partial assignment to be checked ..." << certificateString << endl;
                } else {
                    certificate = true;
                    certificateFilename = (char *) "sol";
                    if (ToulBar2::debug) cout << "certificate of solution read in file: ./" << certificateFilename << endl;

                }
            }

            // show Solutions
            if (args.OptionId() == OPT_showSolutions ) {
                ToulBar2::showSolutions = true;
            }



            //#############################################
            if (args.OptionId() == OPT_writeSolution )

            {
                ToulBar2::writeSolution = (char *) "sol";
                if(args.OptionArg() != NULL) {
                    char *tmpFile = new char[strlen(args.OptionArg())+1];
                    strcpy(tmpFile,args.OptionArg());
                    if (strlen(tmpFile) == 1 && (tmpFile[0] == '0' || tmpFile[0] == '1' || tmpFile[0] == '2')) ToulBar2::pedigreeCorrectionMode = atoi(tmpFile);
                    else ToulBar2::writeSolution = tmpFile;
                }
            }


            if ( args.OptionId() == OPT_pedigreePenalty)
            {
                ToulBar2::pedigreePenalty = 1;
                int penaltyThreshold = atoi(args.OptionArg());
                if (penaltyThreshold >= 1) ToulBar2::pedigreePenalty = penaltyThreshold;
            }
            // counting
            if (args.OptionId() == OPT_allSolutions)
            {
                ToulBar2::allSolutions = true;
            }


            // approximate conting

            if (args.OptionId() == OPT_approximateCountingBTD)
            {
                ToulBar2::approximateCountingBTD = true;
                ToulBar2::allSolutions = true;
                ToulBar2::btdMode = 1;
            }


            if (args.OptionId() == OPT_binaryBranching	) { ToulBar2::binaryBranching = true; }
            if (args.OptionId() == NO_OPT_binaryBranching   ) { ToulBar2::binaryBranching = false; }

            // static variable ordering

            if (args.OptionId() == OPT_Static_variable_ordering) {
                ToulBar2::Static_variable_ordering = true;
            } else if (args.OptionId() == NO_OPT_Static_variable_ordering) { ToulBar2::Static_variable_ordering= false; }

            // last conflict

            if (args.OptionId() == OPT_lastConflict ) {
                ToulBar2::binaryBranching = true;
                ToulBar2::lastConflict = true;
            } else if (args.OptionId() == NO_OPT_lastConflict ) { ToulBar2::lastConflict = false; }

            if (args.OptionId() == OPT_dichotomicBranching ) {
                if(args.OptionArg() == NULL ) {
                    ToulBar2::dichotomicBranching = 1;
                } else {
                    int dico = atoi(args.OptionArg());
                    if (dico > 0) ToulBar2::dichotomicBranching = dico;
                    cout << "dichotomicBranching=" << ToulBar2::dichotomicBranching << endl;
                }
            } else if (args.OptionId() == NO_OPT_dichotomicBranching ) {
                ToulBar2::dichotomicBranching = 0;
            }

            if (args.OptionId() == OPT_sortDomains ) {
                ToulBar2::sortDomains = true;
            } else if (args.OptionId() == NO_OPT_sortDomains ) {
                ToulBar2::sortDomains = false;
            }

            if (args.OptionId() == OPT_weightedTightness) {
                if (args.OptionArg() != NULL) {
                    int weightedtight = atol(args.OptionArg());
                    if (weightedtight >= 0) ToulBar2::weightedTightness = weightedtight;
                } else {
                    ToulBar2::weightedTightness = 2;
                }
                if (ToulBar2::weightedTightness && !ToulBar2::weightedDegree) ToulBar2::weightedDegree = 10000;
            } else if (args.OptionId() == NO_OPT_weightedTightness ) { ToulBar2::weightedTightness = 0; }

            // weitghted Degree (var ordering )
            if (args.OptionId() == OPT_weightedDegree and args.OptionArg() != NULL ) {
                int weighteddegree = atol(args.OptionArg());
                if (weighteddegree > 0) ToulBar2::weightedDegree = weighteddegree;
            }   else if ( args.OptionId() == NO_OPT_weightedDegree )
            {
                ToulBar2::weightedDegree = 0;
                ToulBar2::weightedTightness = 0;
                if (ToulBar2::debug) cout << "ToulBar2::weightedDegree = false" << endl;
            }

            // LIMIT BRANCHING ON FIRST nbDecisionVars VARIABLES OPTION
            if ( args.OptionId() == OPT_nbDecisionVars )
            {
                ToulBar2::nbDecisionVars = atoi(args.OptionArg());
            }

            //////////////////////////////////////////////////////

            if ( args.OptionId() == NO_OPT_elimDegree )
            {
                ToulBar2::elimDegree = -1;
                if (ToulBar2::debug) cout << "elimDegree OFF " << ToulBar2::elimDegree << endl;

            } else if ( args.OptionId() == OPT_elimDegree  and args.OptionArg() != NULL ) {
                int ndegree =-1;
                ndegree = atol(args.OptionArg());
                if ((ndegree >= 0) && (ndegree <= 3)) {ToulBar2::elimDegree = ndegree;}
                if (ToulBar2::debug) cout << "elimDegree ON " << ToulBar2::elimDegree << endl;
            }
            ////////////////////////////////////////////

            // p[integer]: preprocessing only: general variable elimination of degree less than or equal to the given value
            // elimination degree in preprocessing
            if ( args.OptionId() == OPT_elimDegree_preprocessing)
            {
                int ndegree;

                if(args.OptionArg() !=NULL) {
                    ndegree = atoi(args.OptionArg());
                    if (ndegree != 0) ToulBar2::elimDegree_preprocessing = ndegree;
                    if (ndegree<0) ToulBar2::elimSpaceMaxMB = 128;
                } else ToulBar2::elimDegree_preprocessing = 3;
                if (ToulBar2::debug) cout << "elimDegree_preprocessing ON: "  << ToulBar2::elimDegree_preprocessing << endl;

            }   else if ( args.OptionId() == NO_OPT_elimDegree_preprocessing)
            {
                ToulBar2::elimDegree_preprocessing = -1;
                if (ToulBar2::debug) cout << "elimDegree_preprocessing OFF: " << ToulBar2::elimDegree_preprocessing << endl;
            }



            // VAC PARAMETER
            if ( args.OptionId() == OPT_minsumDiffusion)

            {
                if (!ToulBar2::vac) ToulBar2::vac = 1;
                ToulBar2::minsumDiffusion = 1000;
                int nit = atoi(args.OptionArg());
                if (nit > 0) ToulBar2::minsumDiffusion = nit;
            }

            if ( args.OptionId() == OPT_vac)
            {
                ToulBar2::vac = 1;
                if(args.OptionArg() == NULL)
                { ToulBar2::vac = 1 ;}
                else {
                    int depth =  atoi(args.OptionArg());
                    if (depth >= 1) ToulBar2::vac = depth;
                }
                if (ToulBar2::debug) cout <<"VAC propagation ON"<< endl;

            } else if( args.OptionId() == NO_OPT_vac)
            {
                if (ToulBar2::debug) cout << "VAC propagation OFF"<< endl;

                ToulBar2::vac = 0;
            }

            if ( args.OptionId() == OPT_costThreshold)
            {
                Cost ct = string2Cost(args.OptionArg());
                if (ct > UNIT_COST) ToulBar2::costThreshold = ct;
            }

            if ( args.OptionId() == OPT_costThresholdPre)
            {
                Cost ct = string2Cost(args.OptionArg());
                if (ct > UNIT_COST) ToulBar2::costThresholdPre = ct;
            }
            /*if ( (ch = strchr(argv[i],'R')) ) {
			  Cost ct = string2Cost(&ch[1]);
			  if(ct > UNIT_COST) ToulBar2::relaxThreshold = ct;
			  }*/


            if ( args.OptionId() == OPT_costMultiplier)
            {
                double co = atof(args.OptionArg());
                if (co > MIN_COST) ToulBar2::costMultiplier = co;
            }

            if ( args.OptionId() == OPT_singletonConsistency)  ToulBar2::singletonConsistency = true;
            if ( args.OptionId() == OPT_vacValueHeuristic) ToulBar2::vacValueHeuristic = true;
            if ( args.OptionId() == OPT_preprocessTernary) {
                if (args.OptionArg() != NULL) {
                    int size = atol(args.OptionArg());
                    if (size>=0) ToulBar2::preprocessTernaryRPC = size;
                } else ToulBar2::preprocessTernaryRPC = 128;
            } else if ( args.OptionId() == NO_OPT_preprocessTernary) {
                if (ToulBar2::debug) cout <<"preprocess triangles of binary cost functions into ternary cost functions OFF" << endl;
                ToulBar2::preprocessTernaryRPC = 0;
            }

            // elimination of functional variables
            if ( args.OptionId() == OPT_preprocessFunctional)
            {
                if(args.OptionArg() == NULL ) {
                    ToulBar2::preprocessFunctional = 1;
                } else {
                    int func = atoi(args.OptionArg());
                    if (func > 0) ToulBar2::preprocessFunctional = func;
                }
                if (ToulBar2::debug) cout << "elimination of functional variables ON "<< ToulBar2::preprocessFunctional << endl;

            } else if( args.OptionId() == NO_OPT_preprocessFunctional)
            {
                if (ToulBar2::debug) cout <<"elimination of functional variables OFF" << endl;
                ToulBar2::preprocessFunctional  = 0;

            }

            if (args.OptionId() == OPT_costfuncSeparate)
            {
                if (ToulBar2::debug) cout << "decomposition of cost functions" << endl;
                ToulBar2::costfuncSeparate = true;
            } else if (args.OptionId() == NO_OPT_costfuncSeparate)
            {
                if (ToulBar2::debug) cout << "decomposition of cost functions OFF" << endl;
                ToulBar2::costfuncSeparate = false;
            }

            if (args.OptionId() == NO_OPT_DEE)
            {
                if (ToulBar2::debug) cout << "dead-end elimination OFF" << endl;
                ToulBar2::DEE = 0;
            } else if (args.OptionId() == OPT_DEE) {
                ToulBar2::DEE = 1;
                if (args.OptionArg() != NULL) {
                    int dee = atoi(args.OptionArg());
                    if (dee >= 0) ToulBar2::DEE = dee;
                }
                if (ToulBar2::debug) cout << "dead-end elimination: " << ToulBar2::DEE << endl;
            }

            // pre projection of nary cost functions
            if ( args.OptionId() == OPT_preprocessNary)
            {
                if(args.OptionArg() == NULL ) {
                    ToulBar2::preprocessNary = 10;
                } else {
                    int maxnary = atoi(args.OptionArg());
                    if (maxnary > 0) ToulBar2::preprocessNary = maxnary;
                }
                if (ToulBar2::debug) cout <<"preproject cost functions with arity lower than "<< ToulBar2::preprocessNary << " ON" << endl;

            } else if( args.OptionId() == NO_OPT_preprocessNary)
            {
                if (ToulBar2::debug) cout <<"preproject of n-ary cost functions OFF" << endl;
                ToulBar2::preprocessNary  = 0;
            }

            if ( args.OptionId() == OPT_QueueComplexity) ToulBar2::QueueComplexity = true;
            if ( args.OptionId() == OPT_MSTDAC) ToulBar2::MSTDAC = true;
            else if ( args.OptionId() == NO_OPT_MSTDAC) ToulBar2::MSTDAC = false;

            // LDS
            if ( args.OptionId() == OPT_lds)
            {
                string comment ;
                if(args.OptionArg() == NULL ) {
                    ToulBar2::lds = maxdiscrepancy;
                    comment = " (default value) " ;
                } else {
                    int maxlds =atoi(args.OptionArg());
                    if (maxlds > 0 || maxlds < 0) ToulBar2::lds = maxlds;
                }
                if (ToulBar2::debug) cout << "LDS ON #iter = " << ToulBar2::lds << comment << endl ;

            } else if (args.OptionId() ==NO_OPT_lds) {

                ToulBar2::lds = 0;
                if (ToulBar2::debug) cout << "LDS OFF iter = " << ToulBar2::lds << endl ;
            }


            // restart option
            if ( args.OptionId() == OPT_restart)
            {
                string comment ;
                if( args.OptionArg() == NULL ) {
                    comment = " (default value) " ;
                    ToulBar2::restart = maxrestarts;
                } else {
                    Long maxbt = atoll(args.OptionArg());
                    if (maxbt >= 0) ToulBar2::restart = maxbt;
                }
                if (ToulBar2::debug) cout << "restart ON #iter = " << ToulBar2::restart << comment << endl;
            } else if (args.OptionId() ==NO_OPT_restart)
            {
                if (ToulBar2::debug) cout <<"restart OFF" << endl;
                ToulBar2::restart = -1;
            }

            // hybrid BFS option
            if ( args.OptionId() == OPT_hbfs)
            {
                string comment ;
                if( args.OptionArg() == NULL ) {
                    comment = " (default value) " ;
                    ToulBar2::hbfsGlobalLimit = hbfsgloballimit;
                } else {
                    Long maxbt = atoll(args.OptionArg());
                    if (maxbt > 0) ToulBar2::hbfsGlobalLimit = maxbt;
                    else ToulBar2::hbfsGlobalLimit = hbfsgloballimit;
                }
                ToulBar2::hbfs = 1; // initial value to perform a greedy search exploration before visiting a new open search node
                if (ToulBar2::debug) cout << "hybrid BFS ON with global backtrack limit = " << ToulBar2::hbfsGlobalLimit << comment << endl;
            } else if (args.OptionId() ==NO_OPT_hbfs)
            {
                if (ToulBar2::debug) cout <<"hybrid BFS OFF" << endl;
                ToulBar2::hbfs = 0;
                ToulBar2::hbfsGlobalLimit = 0;
            }
            if ( args.OptionId() == OPT_open)
            {
                ToulBar2::hbfsOpenNodeLimit = OPEN_NODE_LIMIT;
                Long openlimit = atoll(args.OptionArg());
                if (openlimit>0) ToulBar2::hbfsOpenNodeLimit = openlimit;
                if (ToulBar2::hbfsGlobalLimit==0) ToulBar2::hbfsGlobalLimit = hbfsgloballimit;
                ToulBar2::hbfs = 1; // initial value to perform a greedy search exploration before visiting a new open search node
                if (ToulBar2::debug) cout << "hybrid BFS ON with open node limit = " << ToulBar2::hbfsOpenNodeLimit << endl;
            }

            // local search INCOP
            if ( args.OptionId() == OPT_localsearch)  {
                if (args.OptionArg() != NULL) {
                    ToulBar2::incop_cmd = args.OptionArg();
                } else {
                    ToulBar2::incop_cmd = Incop_cmd;
                }
            }

            // EDAC OPTION
            if ( args.OptionId() == OPT_EDAC )
            {
                ToulBar2::LcLevel = LC_EDAC;
                LcLevelType lclevel = (LcLevelType) atoi(args.OptionArg());
                if ((lclevel >= LC_NC) && (lclevel < LC_THEMAX)) ToulBar2::LcLevel = lclevel;
            }

            ////////////////////////////////////////////////////////////////////////////////////////
            ////////		MEDELSOFT command line processing 	////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////

            // g  sort pedigree by increasing generation number and if equal by increasing individual number

            if ( args.OptionId() == OPT_generation)
            {
                ToulBar2::generation = true;
                if (ToulBar2::debug) cout << "-g flag ==> sort pedigree option ON" <<  endl;

            }

            // bayesien -y [genotypinpErrorRate probabilityPrecision genotypePriorMode]  : pedigree solved by Bayesian MPE

            if ( args.OptionId() == MENDEL_OPT_genotypingErrorRate)
            {
                if(args.OptionArg() != NULL) {
                    float f;
                    sscanf(args.OptionArg(),"%f",&f) ;
                    ToulBar2::errorg=f;
                    if (ToulBar2::debug) cout << "New assignment for genotyping Error Rate = " << ToulBar2::errorg  <<  endl;
                }

            }

            if ( args.OptionId() == MENDEL_OPT_resolution)
            {
                if(args.OptionArg() != NULL) {
                    ToulBar2::resolution=atoi(args.OptionArg());
                    if (ToulBar2::debug) cout << "New assignment for precision = " << ToulBar2::resolution  <<  endl;
                }

            }

            if ( args.OptionId() == OPT_pedigree_by_MPE)
            {
                ToulBar2::bayesian = true;
                if (ToulBar2::debug) cout << endl ;
                if (ToulBar2::debug) cout << "<<< Bayesian mode ON >>>" << endl;

                if(ToulBar2::verbose >= 1) {
                    cout << endl ;
                    cout << "genotypingErrorRate = " << ToulBar2::errorg << " (default value)" <<  endl;
                    cout << "resolution = " << ToulBar2::resolution << " (default value)" <<  endl;
                    cout << endl ;
                }



            }

            if ( args.OptionId() == MENDEL_OPT_EQUAL_FREQ && ToulBar2::bayesian)
            {

                ToulBar2::foundersprob_class=0;
                if (ToulBar2::debug) cout << "equal frequencies used ( default mode)" << ToulBar2::foundersprob_class << endl;

            } else if ( args.OptionId() == MENDEL_OPT_ESTIMAT_FREQ && ToulBar2::bayesian)
            {
                ToulBar2::foundersprob_class=1;
                if (ToulBar2::debug) cout << " => prob depending on the frequencies found in the current pedigree probleme used" << endl;


            } else if ( args.OptionId() == MENDEL_OPT_ALLOCATE_FREQ && ToulBar2::bayesian)
            {
                if (ToulBar2::debug) cout << " => prob frequencies read from command line" << endl;
                ToulBar2::foundersprob_class=2;
                Pedi_Args(args,-1);
            }

            ////////////////////////////////////////////////////////////////////////////////////////
            /////				DEBUG and verbose LEVEL  management + DUMP
            ////////////////////////////////////////////////////////////////////////////////////////
            //verbose mode
            //   v: verbosity level
            if ( args.OptionId() == OPT_verbose) {
                if (args.OptionArg() != NULL) { ToulBar2::verbose= atoi(args.OptionArg()); }
                else { ToulBar2::verbose=0; }
                if (ToulBar2::debug) cout <<"verbose level = " << ToulBar2::verbose << endl;


            }

            //  z: save problem in wcsp format in filename \"problem.wcsp\" (1:before, 2:current problem after preprocessing)

            if ( args.OptionId() == OPT_dumpWCSP) {
                if (!ToulBar2::dumpWCSP) ToulBar2::dumpWCSP=1;
                if (args.OptionArg() != NULL) {
                    char *tmpFile = new char[strlen(args.OptionArg())+1];
                    strcpy(tmpFile,args.OptionArg());
                    if (strlen(tmpFile) == 1 && (tmpFile[0] == '1' || tmpFile[0] == '2')) ToulBar2::dumpWCSP = atoi(tmpFile);
                    else ToulBar2::problemsaved_filename = to_string(tmpFile);
                }

                if(ToulBar2::dumpWCSP <= 1 ) {
                    if (ToulBar2::debug) cout <<"original problem dump in problem.wcsp (see also Graphviz and degree distribution)"  << endl;
                } else {
                    if (ToulBar2::debug) cout <<"dump after preprocessing in problem.wcsp (see also Graphviz and degree distribution)"  << endl;

                }


            }
            //   Z: debug mode (save problem at each node if verbosity option set!)

            //		for (int j=0; argv[i][j] != 0; j++) if (argv[i][j]=='Z') ToulBar2::debug++;

            if ( args.OptionId() == OPT_debug) {
                if (args.OptionArg() != NULL) { ToulBar2::debug= atoi(args.OptionArg());} else ToulBar2::debug=1;
                if (ToulBar2::debug) cout <<"debug level = " << ToulBar2::debug << endl;


            }

            // discrete integration for computing the partition function Z
            if ( args.OptionId() == OPT_Z) ToulBar2::isZ = true;

            if ( args.OptionId() == OPT_epsilon)
            {
                if(args.OptionArg() != NULL) {
                    ToulBar2::logepsilon=-Log(atol(args.OptionArg()));
                    if (ToulBar2::debug) cout << "New assignment for epsilon = " << Exp(-ToulBar2::logepsilon)  <<  endl;
                }

            }

            if ( args.OptionId() == OPT_learning )
            {
                ToulBar2::learning = true;
            }

#ifndef NDEBUG
            if ( args.OptionId() == OPT_verifyopt )
                ToulBar2::verifyOpt = true;
#endif

            // upper bound initialisation from command line
            if ( args.OptionId() == OPT_ub) {

                if (args.OptionArg() != NULL) {
                    ub = (args.OptionArg())?string2Cost(args.OptionArg()):MAX_COST;
                }
                if (ToulBar2::debug) cout <<"UB =" << ub << " passed in  command line" << endl;
            }

            // CPU timer
            if ( args.OptionId() == OPT_timer)
            {
                if(args.OptionArg() != NULL) {
                    timeout=atoi(args.OptionArg());
                    if (ToulBar2::debug) cout << "TimeOut = " << timeout <<  endl;
                }

            }


            //////////RANDOM GENERATOR///////
            if ( args.OptionId() == OPT_seed ) {
                mysrand(atol(args.OptionArg()));
            }

            if ( args.OptionId() == OPT_random) {

                random_desc = args.OptionArg();


            }

            if (args.OptionId() == OPT_nopre)
            {
                ToulBar2::elimDegree = -1;
                if (ToulBar2::debug) cout << "elimDegree OFF " << ToulBar2::elimDegree << endl;
                ToulBar2::elimDegree_preprocessing = -1;
                if (ToulBar2::debug) cout << "elimDegree_preprocessing OFF: " << ToulBar2::elimDegree_preprocessing << endl;
                if (ToulBar2::debug) cout <<"preprocess triangles of binary cost functions into ternary cost functions OFF" << endl;
                ToulBar2::preprocessTernaryRPC = 0;
                if (ToulBar2::debug) cout <<"elimination of functional variables OFF" << endl;
                ToulBar2::preprocessFunctional  = 0;
                if (ToulBar2::debug) cout <<"preproject of n-ary cost functions OFF" << endl;
                ToulBar2::preprocessNary  = 0;
                if (ToulBar2::debug) cout << "decomposition of cost functions OFF" << endl;
                ToulBar2::costfuncSeparate = false;
                if (ToulBar2::debug) cout << "maximum spanning tree DAC ordering OFF" << endl;
                ToulBar2::MSTDAC = false;
                if (ToulBar2::debug) cout << "dead-end elimination OFF" << endl;
                ToulBar2::DEE = 0;
            }

        }


        else
        {
            cout <<  "<<<< ERROR >>>>  : "<< endl ;
            _tprintf(
                    _T("%s: '%s' (use --help to get command line help)\n"),
                    GetLastErrorText(args.LastError()), args.OptionText());

            _tprintf(_T("invalid Option :%s or Invalid argument: %s \n please check help using --help option "), args.OptionText(),  args.OptionArg() ? args.OptionArg() : "");

            return 1 ;
        }
    }
    // process any files that were passed to us on the command line.
    // send them to the globber so that all wildcards are expanded
    // into valid filenames (e.g. *.cpp -> a.cpp, b.cpp, c.cpp, etc)
    // See the SimpleGlob.h header file for details of the flags.
    CSimpleGlob glob(SG_GLOB_NODOT|SG_GLOB_NOCHECK);
    if (SG_SUCCESS != glob.Add(args.FileCount(), args.Files())) {
        _tprintf(_T("Error while globbing files\n"));
        return 1;
    }

    // dump all of the details, the script that was passed on the
    // command line and the expanded file names
    if (ToulBar2::verbose > 0 )  {
        cout << "cmd line parsing ---> "<< glob.FileCount() << " Filename(s) found in command line"<< endl;
    }
    string strfile;
    string strext;



    if(random_desc ==NULL)
    {
        for (int n = 0; n < glob.FileCount(); ++n) {
            if (ToulBar2::verbose > 0 )   _tprintf(_T("file %d: '%s'\n"), n, glob.File(n));

            // wcsp input file  == first file
            //	if (strstr(glob.File(n),".wcsp"))

            if(check_file_ext(glob.File(n),file_extension_map["wcsp_ext"]) )
            {
                if (ToulBar2::verbose >= 0) cout <<  "loading wcsp file: "<< glob.File(n)  << endl;
                strext = ".wcsp";
                strfile = glob.File(n);
            }
            // uai  file
            if(check_file_ext(glob.File(n),file_extension_map["uai_ext"]) ) {
                strfile = glob.File(n);
                strext = ".uai";
                cout <<  "loading uai file:  "<< glob.File(n) << endl;
                ToulBar2::uai = 1;
                ToulBar2::bayesian = true;
            }
            // uai log file
            if(check_file_ext(glob.File(n),file_extension_map["uai_log_ext"]) ) {
                strfile = glob.File(n);
                strext = ".LG";
                cout <<  "loading uai log file:  "<< glob.File(n) << endl;
                ToulBar2::uai = 2;
                ToulBar2::bayesian = true;
            }
            // UAI evidence file
            if(check_file_ext(glob.File(n),file_extension_map["evid_ext"]) ) {
                cout <<  "loading evidence file:  "<< glob.File(n) << endl;
                ToulBar2::evidence_file = string(glob.File(n));
            }

            // xml file
            if(check_file_ext(glob.File(n),file_extension_map["wcspXML_ext"]) ) {
                cout <<  "loading xml file:" << glob.File(n) << endl;

                ToulBar2::xmlflag = true;
                if (!ToulBar2::writeSolution) ToulBar2::writeSolution = (char *) "sol";
                strext = ".xml";
                strfile = glob.File(n);
            }

            // wcnf or cnf file

            if(check_file_ext(glob.File(n),file_extension_map["wcnf_ext"]) ) {
                if (ToulBar2::verbose >= 0) cout <<  "loading wcnf file:" << glob.File(n) << endl;
                ToulBar2::wcnf = true;
                strext = ".wcnf";
                strfile = glob.File(n);
            } else if(check_file_ext(glob.File(n),file_extension_map["cnf_ext"]) ) {
                if (ToulBar2::verbose >= 0) cout <<  "loading cnf file:" << glob.File(n) << endl;
                ToulBar2::wcnf = true;
                strext = ".cnf";
                strfile = glob.File(n);
            }

            // unconstrained quadratic programming file

            if(check_file_ext(glob.File(n),file_extension_map["qpbo_ext"]) ) {
                cout <<  "loading quadratic pseudo-Boolean optimization file:" << glob.File(n) << endl;
                ToulBar2::qpbo = true;
                strext = ".qpbo";
                strfile = glob.File(n);
            }

            // upperbound file

            if(check_file_ext(glob.File(n),file_extension_map["ub_ext"]) ) {
                cout <<  "loading upper bound from file: "<< glob.File(n)  << endl;
                string ubs;
                ubs=read_UB(glob.File(n));
                if(ubs.c_str() != NULL ) {
                    ub  = string2Cost((char*) ubs.c_str() );
                } else {
                    cerr << "error reading UB in " << glob.File(n) << endl;
                    exit (-1);
                }

            }



            // bep input file
            if(check_file_ext(glob.File(n),file_extension_map["bep_ext"]) ) {
                cout <<  "loading BEP file: "<< glob.File(n)  << endl;
                strext = ".bep";
                strfile = glob.File(n);
            }

            //////////////////////Mendelian-error analysis and haplotype reconstruction ////////////////////////////////////
            // map file

            if(check_file_ext(glob.File(n),file_extension_map["map_ext"]) ) {
                ToulBar2::map_file = string(glob.File(n));
                ToulBar2::haplotype = new Haplotype;
                cout << "loading map file: " << ToulBar2::map_file << endl;

                if( glob.FileCount() <2) {
                    cerr << "pedigree file is missing (.pre): " << endl;
                    exit(-1);


                }
            }

            // pre file
            if(check_file_ext(glob.File(n),file_extension_map["pre_ext"]) ) {
                strfile = glob.File(n);
                strext = ".pre";
                cout << "loading pre file: " << glob.File(n) << endl;
                if( glob.FileCount() <2) ToulBar2::pedigree = new Pedigree;
            }

            //////////////////////VARIABLE ORDERING ////////////////////////////////////
            // filename containing variable order
            //if (strstr(glob.File(n),".order"))
            if(check_file_ext(glob.File(n),file_extension_map["order_ext"]) )
            {
                cout << "loading variable order in file: " << glob.File(n) << endl;
                char buf[80];
                sprintf(buf,"%s",glob.File(n));
                //				if (ToulBar2::varOrder) delete [] ToulBar2::varOrder;
                ToulBar2::varOrder = new char [ strlen(buf) + 1 ];
                sprintf(ToulBar2::varOrder, "%s",buf);
            }

            // read assignment in file or filename of solution
            if(check_file_ext(glob.File(n),file_extension_map["sol_ext"]) )
            {
                if (certificateString && strcmp(certificateString, "")!=0)
                {
                    cerr << "\n COMMAND LINE ERROR cannot read a solution if a partial assignment is given in the command line using -x= argument " << endl ;
                    exit(-1) ;
                }
                cout << "loading solution in file: " << glob.File(n) << endl;

                certificate = true;
                certificateFilename = new char[256];
                certificateFilename= glob.File(n);
                certificateString = (char *) "";  // ensure the search will continue starting from this solution
            }




        }

    }

    // ----------simple opt end ----------------------

    //  command line => check existencise of problem filname argument

    // PB FILENAME

    if (argc <= 1)
    {
        cerr << "Problem filename is missing as command line argument!" << endl;
        cerr << endl;
        help_msg(argv[0]);
        exit(0);
    }


    //------------------------------tb2 option --------------
    /////////////////////////////////////////
    // test on initial ub cost value;
    /////////////////////////////////////////

    ub = tb2checkOptions(ub);
    if (ToulBar2::verifyOpt && (!certificate || certificateFilename == NULL)) {
        cout << "Warning! An optimal solution file is missing in the command line. Verifying the optimal solution is disabled." << endl;
        ToulBar2::verifyOpt = false;
    }

    if (ToulBar2::uai || ToulBar2::uaieval) {
        char *tmpPath = new char[strlen(argv[0])+1];
        strcpy(tmpPath,argv[0]);
        char *tmpFile = new char[strlen(strfile.c_str())+1];
        strcpy(tmpFile,strfile.c_str());
        string filename(tmpPath);
        filename += "/";
        filename += basename(tmpFile);
        size_t wcsppos = string::npos;
        if (ToulBar2::uaieval && (wcsppos = filename.rfind( ".wcsp" )) != string::npos) filename.replace( wcsppos, 5, ".uai" );
        filename += ".";
        if (ToulBar2::isZ) filename += "PR";
        else filename += "MPE";
        ToulBar2::solution_uai_filename = filename;
        if (!ToulBar2::uaieval) {
            ToulBar2::solution_file.open(ToulBar2::solution_uai_filename.c_str());
            ToulBar2::solution_file << ((ToulBar2::isZ)?"PR":"MPE") << endl;
            ToulBar2::solution_file.flush();
        }
        delete [] tmpPath;
        delete [] tmpFile;
    }

    //TODO: If --show_options then dump ToulBar2 object here

    ToulBar2::startCpuTime = cpuTime();

    initCosts(ub);
    WeightedCSPSolver *solver = WeightedCSPSolver::makeWeightedCSPSolver(STORE_SIZE, ub);

    bool randomproblem = false;
    bool forceSubModular = false;

    int n = 10;
    int m = 2;
    int seed = 3;
    vector<int> p;
#ifndef MENDELSOFT
    if (random_desc !=NULL)
    {
        cout <<"random test ON"<<endl;
        int pn[10];
        int narities = 0;
        if (strstr(random_desc,"bin"))
        {
            randomproblem = true;
            sscanf(random_desc, "bin-%d-%d-%d-%d-%d", &n, &m, &pn[0], &pn[1],&seed);
            narities = 2;
        }
        if (strstr(random_desc,"binsub"))
        {
            forceSubModular = true;
            randomproblem = true;
            sscanf(random_desc, "binsub-%d-%d-%d-%d-%d-%d", &n, &m, &pn[0], &pn[1], &pn[2], &seed);
            narities = 2;
        }
        if (strstr(random_desc,"tern"))
        {
            randomproblem = true;
            sscanf(random_desc, "tern-%d-%d-%d-%d-%d-%d", &n, &m, &pn[0], &pn[1], &pn[2],&seed);
            narities = 3;
        }
        if (strstr(random_desc,"nary"))
        {
            randomproblem = true;
            char* pch = strtok (random_desc,"-");
            pch = strtok (NULL, "-");
            n = atoi(pch);
            pch = strtok (NULL, "-");
            m = atoi(pch);

            while (pch != NULL)
            {
                pch = strtok (NULL, "-");
                if (pch != NULL)
                {
                    pn[narities] = atoi(pch);
                    narities++;
                }
            }
            narities--;
            seed = pn[narities];
        }
        if (pn[0] > 100)
        {
            cout << pn[0] << " tightness is a percentage" << endl;
            pn[0] = 100;
        }
        for (int i=0;i<narities;i++) p.push_back( pn[i] );
        if (forceSubModular) p.push_back( pn[narities] );
        if (narities == 0)
        {
            //help_msg(argv[0]);
            //exit(0);
        }
    }
    if (strstr(strext.c_str(),".bep") || strstr((char*)strfile.c_str(),"bEpInstance")) ToulBar2::bep = new BEP;
#endif
    try
    {
        if (randomproblem)    solver->read_random(n,m,p,seed,forceSubModular);
        else 		         solver->read_wcsp((char*)strfile.c_str());

        if (certificate)
        {
            if (certificateFilename!= NULL) solver->read_solution(certificateFilename);
            else solver->parse_solution(certificateString);
        }
        if (ToulBar2::dumpWCSP==1) {
            string problemname = ToulBar2::problemsaved_filename;
            if (ToulBar2::uaieval) {
                problemname = ToulBar2::solution_uai_filename;
                problemname.replace( problemname.rfind( ".uai.MPE" ), 8, ".wcsp" );
            }
            solver->dump_wcsp(problemname.c_str());
        }
        else if (!certificate || certificateString!=NULL || ToulBar2::btdMode>=2)
        {
#ifdef LINUX
            signal(SIGINT,  timeOut);
            if (timeout>0) timer(timeout);
#endif
            solver->solve();
        }
    }
    catch (Contradiction)
    {
        if (ToulBar2::verbose >= 0) cout << "No solution found by initial propagation!" << endl;
        if (ToulBar2::isZ) {
            if (ToulBar2::uai) {
                if (ToulBar2::uai_firstoutput) ToulBar2::uai_firstoutput = false;
                else ToulBar2::solution_file << "-BEGIN-" << endl;
                ToulBar2::solution_file << "1" << endl;
                ToulBar2::solution_file << -numeric_limits<TProb>::infinity() << endl;
                ToulBar2::solution_file.flush();
            }
            cout << "Log(Z)= ";
            cout << -numeric_limits<TProb>::infinity() << endl;
        }
        if (ToulBar2::maxsateval) {
            cout << "o " << solver->getWCSP()->getUb() << endl;
            cout << "s UNSATISFIABLE" << endl;
        }
    }
    if (ToulBar2::verbose >= 0) cout << "end." << endl;
    if (ToulBar2::uai || ToulBar2::uaieval) ToulBar2::solution_file.close();

    // for the competition it was necessary to write a file with the optimal sol
    /*char line[80];
	  string strfile(argv[1]);
	  int pos = strfile.find_last_of(".");
	  string strfilewcsp = strfile.substr(0,pos) + ".ub";
	  sprintf(line,"echo %d > %s",(int)solver->getWCSP()->getUb(),strfilewcsp.c_str());
	  system(line); */

    return 0;
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */


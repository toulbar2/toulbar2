/*
 * **************** Main function ***********************************
 */

#include "toulbar2lib.hpp"
#include "applis/tb2pedigree.hpp"
#include "applis/tb2haplotype.hpp"
#include "applis/tb2bep.hpp"
#include "vns/tb2vnsutils.hpp"
#include "vns/tb2dgvns.hpp"
#ifdef OPENMPI
#include "vns/tb2cpdgvns.hpp"
#include "vns/tb2rpdgvns.hpp"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include <cfenv>

const int maxdiscrepancy = 4;
const Long maxrestarts = 10000;
const Long hbfsgloballimit = 10000;
const int raspsangle = 10;
const Long raspsbacktracks = 1000;
const double relativegap = 0.0001;
const int maxdivnbsol = 1000;

// INCOP default command line option
const string Incop_cmd = "0 1 3 idwa 100000 cv v 0 200 1 0 0";

//* definition of path separtor depending of OS '/'  => Unix ;'\' ==> windows
#ifdef __WIN32__
#define PATH_SEP_CHR '\\'
#define PATH_DELIM ";"
#else
#define PATH_SEP_CHR '/'
#define PATH_DELIM ":"
#endif
//*definition of  windows include for command line.
#ifdef __WIN32__
#include <windows.h>
#include <tchar.h>
#else
#define TCHAR char
#define _T(x) x
#define _tprintf printf
#define _tmain main
#define _ttoi atoi
#endif
#include "utils/SimpleOpt.h"
#include "utils/SimpleGlob.h"

// used for debugging purpose.
// under gdb: p ((BinaryConstraint *) constrs[13])->dump
// under gdb: p $2(constrs[13], myCout)
ostream myCout(cout.rdbuf());

void conflict() {}

#ifdef PARETOPAIR_COST
void initCosts()
{
    if (ToulBar2::LcLevel > LC_FDAC) {
        cerr << "EDAC not implemented on Paretopair => force to FDAC." << endl;
        ToulBar2::LcLevel = LC_FDAC;
    }
    if (ToulBar2::vac) {
        cerr << "VAC not implemented on Paretopair." << endl;
        ToulBar2::vac = 0;
        ToulBar2::minsumDiffusion = 0;
    }
    if (ToulBar2::elimDegree >= 0 || ToulBar2::elimDegree_preprocessing >= 0 || ToulBar2::elimDegree_preprocessing < -1) {
        cerr << "Variable elimination not implemented on Paretopair." << endl;
        ToulBar2::elimDegree = -1;
        ToulBar2::elimDegree_preprocessing = -1;
    }
    if (ToulBar2::btdMode >= 1) {
        cerr << "BTD-like methods not implemented on Paretopair." << endl;
        ToulBar2::btdMode = 0;
    }
}
#endif

/* read upper bound from given file */
string read_UB(const char* ubfilename)
{
    string ubs;
    ifstream ubfile(ubfilename);
    if (ubfile.is_open()) {
        ubfile >> ubs;
        cout << "UB read from " << ubfilename << "=" << ubs << endl;
        ubfile.close();

        return ubs;

    } else {
        cout << "Unable to open file:" << ubfilename << endl;
        return NULL;
    }
}

// chek if a filename end as ext
bool check_file_ext(const string fileName, const string ext)
{
    size_t extLen = ext.length();
    return ((extLen <= fileName.length()) && (fileName.substr(fileName.length() - extLen) == ext));
}

enum {
    OPT_verbose = 0,
    OPT_debug,
    OPT_dumpWCSP,
    OPT_HELP,
    OPT_stdin,

    // file extension option
    OPT_wcsp_ext,
    OPT_wcspgz_ext,
    OPT_wcspxz_ext,
    OPT_wcspXML_ext,
    OPT_cfn_ext,
    OPT_cfngz_ext,
    OPT_cfnxz_ext,
    OPT_order_ext,
    OPT_uai_ext,
    OPT_uaigz_ext,
    OPT_uaixz_ext,
    OPT_uai_log_ext,
    OPT_uaigz_log_ext,
    OPT_uaixz_log_ext,
    OPT_evid_ext,
    OPT_map_ext,
    OPT_sol_ext,
    OPT_bep_ext,
    OPT_ub_ext,
    OPT_pre_ext,
    OPT_wcnf_ext,
    OPT_wcnfgz_ext,
    OPT_wcnfxz_ext,
    OPT_cnf_ext,
    OPT_cnfgz_ext,
    OPT_cnfxz_ext,
    OPT_qpbo_ext,
    OPT_qpbogz_ext,
    OPT_qpboxz_ext,
    OPT_opb_ext,
    OPT_opbgz_ext,
    OPT_opbxz_ext,
    OPT_treedec_ext,
    OPT_clusterdec_ext,

    OPT_qpbo_mult,
    // search option
    OPT_SEARCH_METHOD,
    OPT_btdRootCluster,
    OPT_btdSubTree,
    OPT_splitClusterMaxSize,
    OPT_maxSeparatorSize,
    OPT_boostingBTD,
    NO_OPT_boostingBTD,
    OPT_minProperVarSize,
    OPT_varOrder,
    OPT_problemsaved_filename,
    OPT_PARTIAL_ASSIGNMENT,
    NO_OPT_PARTIAL_ASSIGNMENT,
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
    OPT_solutionBasedPhaseSaving,
    NO_OPT_solutionBasedPhaseSaving,
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
    OPT_trwsAccuracy,
    OPT_trwsOrder,
    NO_OPT_trwsOrder,
    OPT_trwsNIter,
    OPT_trwsNIterNoChange,
    OPT_trwsNIterComputeUb,
    NO_OPT_trws,
    OPT_costMultiplier,
    OPT_deltaUbAbsolute,
    OPT_deltaUbRelative,

    OPT_VACINT,
    NO_OPT_VACINT,
    OPT_VACthreshold,
    NO_OPT_VACthreshold,
    OPT_RASPS,
    NO_OPT_RASPS,
    OPT_RASPSangle,
    NO_OPT_RASPSangle,
    OPT_RASPSreset,
    NO_OPT_RASPSreset,
    OPT_RASPSlds,

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
    OPT_btlimit,
    NO_OPT_btlimit,
    OPT_hbfs,
    NO_OPT_hbfs,
    OPT_open,
    OPT_localsearch,
    NO_OPT_localsearch,
    OPT_EDAC,
    OPT_ub,
    OPT_divDist,
    OPT_divWidth,
    OPT_divMethod,
    OPT_divRelax,
    OPT_Z,
    OPT_epsilon,
    OPT_learning,
    OPT_timer,
#ifndef NDEBUG
    OPT_verifyopt,
#endif
    // MENDELESOFT OPTION
    OPT_generation,
    MENDEL_OPT_genotypingErrorRate,
    MENDEL_OPT_resolution,
    OPT_pedigree_by_MPE,
    MENDEL_OPT_EQUAL_FREQ,
    MENDEL_OPT_ESTIMAT_FREQ,
    MENDEL_OPT_ALLOCATE_FREQ,

    // random generator
    OPT_seed,
    OPT_random,

// VNS Methods
#ifdef BOOST
    OPT_VNS_search,
#ifdef OPENMPI
    OPT_CPDGVNS_search,
    OPT_RADGVNS_search,
    OPT_RSDGVNS_search,
    OPT_plimit,
#endif
    OPT_TREEDEC_search,
    OPT_decfile,
    OPT_vns_output,
    OPT_vns_sol_init,
    OPT_lds_min,
    OPT_lds_max,
    OPT_lds_inc,
    OPT_k_min,
    OPT_k_max,
    OPT_k_inc,
    //    OPT_vns_restart_strategy,
    //    OPT_vns_var_heur,
    OPT_neighbor_change,
    OPT_neighbor_synch,
    OPT_optimum
#endif
};

string getExt(string FileName)
{
    // Finds the last persiod character of the string
    int period = FileName.find_last_of(".");
    // I use  + 1 because I don't really need to include the period
    string ext = FileName.substr(period + 1);
    return ext;
}

CSimpleOpt::SOption g_rgOptions[] = {
    { OPT_HELP, (char*)"-h", SO_NONE }, // boolean help
    { OPT_HELP, (char*)"-?", SO_NONE }, // boolean help
    { OPT_HELP, (char*)"-help", SO_NONE }, // boolean help
    { OPT_HELP, (char*)"--help", SO_NONE }, // boolean help
    { OPT_verbose, (char*)"-v", SO_OPT }, // verbose level
    { OPT_debug, (char*)"-Z", SO_OPT }, // debug level
    { OPT_dumpWCSP, (char*)"-z", SO_OPT }, // dump wcsp
    //stdin format
    { OPT_stdin, (char*)"--stdin", SO_REQ_SEP },

    // file extension
    { OPT_wcsp_ext, (char*)"--wcsp_ext", SO_REQ_SEP },
    { OPT_wcspXML_ext, (char*)"--wcspXML_ext", SO_REQ_SEP },
    { OPT_cfn_ext, (char*)"--cfn_ext", SO_REQ_SEP },
    { OPT_cfngz_ext, (char*)"--cfngz_ext", SO_REQ_SEP },
    { OPT_cfnxz_ext, (char*)"--cfnxz_ext", SO_REQ_SEP },
    { OPT_order_ext, (char*)"--order_ext", SO_REQ_SEP },
    { OPT_uai_ext, (char*)"--uai_ext", SO_REQ_SEP },
    { OPT_uaigz_ext, (char*)"--uaigz_ext", SO_REQ_SEP },
    { OPT_uaixz_ext, (char*)"--uaixz_ext", SO_REQ_SEP },
    { OPT_uai_log_ext, (char*)"--uai_log_ext", SO_REQ_SEP },
    { OPT_uaigz_log_ext, (char*)"--uaigz_log_ext", SO_REQ_SEP },
    { OPT_uaixz_log_ext, (char*)"--uaixz_log_ext", SO_REQ_SEP },
    { OPT_evid_ext, (char*)"--evid_ext", SO_REQ_SEP },
    { OPT_map_ext, (char*)"--map_ext", SO_REQ_SEP },
    { OPT_sol_ext, (char*)"--sol_ext", SO_REQ_SEP },
    { OPT_bep_ext, (char*)"--bep_ext", SO_REQ_SEP },
    { OPT_ub_ext, (char*)"--ub_ext", SO_REQ_SEP },
    { OPT_pre_ext, (char*)"--pre_ext", SO_REQ_SEP },
    { OPT_wcnf_ext, (char*)"--wcnf_ext", SO_REQ_SEP },
    { OPT_wcnfgz_ext, (char*)"--wcnfgz_ext", SO_REQ_SEP },
    { OPT_wcnfxz_ext, (char*)"--wcnfxz_ext", SO_REQ_SEP },
    { OPT_cnf_ext, (char*)"--cnf_ext", SO_REQ_SEP },
    { OPT_cnfgz_ext, (char*)"--cnfgz_ext", SO_REQ_SEP },
    { OPT_cnfxz_ext, (char*)"--cnfxz_ext", SO_REQ_SEP },
    { OPT_qpbo_ext, (char*)"--qpbo_ext", SO_REQ_SEP },
    { OPT_qpbogz_ext, (char*)"--qpbogz_ext", SO_REQ_SEP },
    { OPT_qpboxz_ext, (char*)"--qpboxz_ext", SO_REQ_SEP },
    { OPT_opb_ext, (char*)"--opb_ext", SO_REQ_SEP },
    { OPT_opbgz_ext, (char*)"--opbgz_ext", SO_REQ_SEP },
    { OPT_opbxz_ext, (char*)"--opbxz_ext", SO_REQ_SEP },
    { OPT_treedec_ext, (char*)"--treedec_ext", SO_REQ_SEP },
    { OPT_clusterdec_ext, (char*)"--clusterdec_ext", SO_REQ_SEP },

    { OPT_qpbo_mult, (char*)"-qpmult", SO_REQ_SEP },
    { OPT_SEARCH_METHOD, (char*)"-B", SO_REQ_SEP }, // -B [0,1,2] search method
    { OPT_SEARCH_METHOD, (char*)"--search", SO_REQ_SEP },
    { OPT_btdRootCluster, (char*)"-R", SO_REQ_SEP }, // root cluster used in BTD
    { OPT_btdRootCluster, (char*)"--RootCluster", SO_REQ_CMB },
    { OPT_btdSubTree, (char*)"-I", SO_REQ_SEP }, // btd sub tree
    { OPT_splitClusterMaxSize, (char*)"-j", SO_REQ_SEP },
    { OPT_maxSeparatorSize, (char*)"-r", SO_REQ_SEP },
    { OPT_maxSeparatorSize, (char*)"--maxSepSize", SO_REQ_CMB },

    { OPT_minProperVarSize, (char*)"-X", SO_REQ_SEP },
    { OPT_PARTIAL_ASSIGNMENT, (char*)"-x", SO_OPT },
    { NO_OPT_PARTIAL_ASSIGNMENT, (char*)"-x:", SO_NONE },
    { OPT_boostingBTD, (char*)"-E", SO_OPT },
    { NO_OPT_boostingBTD, (char*)"-E:", SO_NONE },
    { OPT_varOrder, (char*)"-O", SO_REQ_SEP }, // filename of variable order
    { OPT_problemsaved_filename, (char*)"--save", SO_REQ_SEP }, // filename of saved problem
    { OPT_showSolutions, (char*)"-s", SO_OPT }, //print solution found
    { OPT_showSolutions, (char*)"--show", SO_OPT }, //print solution found
    { OPT_writeSolution, (char*)"-w", SO_OPT }, //  write last/all solutions found in file (default filename "sol")

    { OPT_pedigreePenalty, (char*)"-u", SO_REQ_SEP }, // int ..
    { OPT_allSolutions, (char*)"-a", SO_OPT }, // counting option ...print solution found
    { OPT_approximateCountingBTD, (char*)"-D", SO_NONE }, //approximate counting
    { OPT_binaryBranching, (char*)"-b", SO_NONE },
    { OPT_binaryBranching, (char*)"-binaryBranching", SO_NONE },
    { NO_OPT_binaryBranching, (char*)"-b:", SO_NONE },
    { NO_OPT_binaryBranching, (char*)"-no--b", SO_NONE },
    { NO_OPT_binaryBranching, (char*)"-no--binaryBranching", SO_NONE },
    { OPT_Static_variable_ordering, (char*)"-svo", SO_NONE },
    { NO_OPT_Static_variable_ordering, (char*)"-svo:", SO_NONE },
    { OPT_lastConflict, (char*)"-c", SO_NONE },
    { NO_OPT_lastConflict, (char*)"-c:", SO_NONE },
    { NO_OPT_lastConflict, (char*)"-no--c", SO_NONE },
    { NO_OPT_lastConflict, (char*)"--lastConflict--off", SO_NONE },
    { OPT_dichotomicBranching, (char*)"-d", SO_OPT },
    { NO_OPT_dichotomicBranching, (char*)"-d:", SO_NONE },
    { OPT_sortDomains, (char*)"-sortd", SO_NONE },
    { NO_OPT_sortDomains, (char*)"-sortd:", SO_NONE },
    { OPT_solutionBasedPhaseSaving, (char*)"-solr", SO_NONE },
    { NO_OPT_solutionBasedPhaseSaving, (char*)"-solr:", SO_NONE },
    { OPT_weightedDegree, (char*)"-q", SO_OPT },
    { NO_OPT_weightedDegree, (char*)"-q:", SO_NONE },
    { OPT_weightedTightness, (char*)"-m", SO_OPT },
    { NO_OPT_weightedTightness, (char*)"-m:", SO_NONE },
    { OPT_nbDecisionVars, (char*)"-var", SO_REQ_SEP },

    { OPT_elimDegree, (char*)"-e", SO_OPT },
    { NO_OPT_elimDegree, (char*)"-e:", SO_NONE },
    { OPT_elimDegree_preprocessing, (char*)"-p", SO_OPT },
    { NO_OPT_elimDegree_preprocessing, (char*)"-p:", SO_NONE },
    { OPT_costfuncSeparate, (char*)"-dec", SO_NONE },
    { NO_OPT_costfuncSeparate, (char*)"-dec:", SO_NONE },
    { OPT_nopre, (char*)"-nopre", SO_NONE },
    // vac option
    { OPT_vac, (char*)"-A", SO_OPT },
    { NO_OPT_vac, (char*)"-A:", SO_NONE },
    { OPT_vacValueHeuristic, (char*)"-V", SO_NONE },
    { NO_OPT_vacValueHeuristic, (char*)"-V:", SO_NONE },
    { OPT_costThreshold, (char*)"-T", SO_REQ_SEP },
    { OPT_costThresholdPre, (char*)"-P", SO_REQ_SEP },
    { OPT_costMultiplier, (char*)"-C", SO_REQ_SEP },

    { OPT_VACINT, (char*)"-vacint", SO_OPT },
    { OPT_VACINT, (char*)"-strictAC", SO_OPT }, //deprecated
    { OPT_VACINT, (char*)"-sac", SO_OPT }, //deprecated
    { NO_OPT_VACINT, (char*)"-vacint:", SO_NONE },
    { NO_OPT_VACINT, (char*)"-sac:", SO_NONE }, //deprecated
    { OPT_VACthreshold, (char*)"-vacthr", SO_NONE },
    { OPT_VACthreshold, (char*)"-VACthreshold", SO_NONE }, //deprecated
    { NO_OPT_VACthreshold, (char*)"-vacthr:", SO_NONE },
    { OPT_RASPS, (char*)"-rasps", SO_OPT },
    { OPT_RASPS, (char*)"-RINS", SO_OPT }, //deprecated
    { OPT_RASPS, (char*)"-rins", SO_OPT }, //deprecated
    { NO_OPT_RASPS, (char*)"-rasps:", SO_NONE },
    { NO_OPT_RASPS, (char*)"-RINS:", SO_NONE }, //deprecated
    { NO_OPT_RASPS, (char*)"-rins:", SO_NONE }, //deprecated
    { OPT_RASPSangle, (char*)"-raspsdeg", SO_REQ_SEP },
    { OPT_RASPSangle, (char*)"-RINSangle", SO_OPT }, //deprecated
    { OPT_RASPSangle, (char*)"-auto", SO_OPT }, //deprecated
    { NO_OPT_RASPSangle, (char*)"-raspsdeg:", SO_NONE },
    { NO_OPT_RASPSangle, (char*)"-auto:", SO_NONE }, //deprecated
    { OPT_RASPSreset, (char*)"-raspsini", SO_NONE },
    { OPT_RASPSreset, (char*)"-RINSreset:", SO_NONE }, //deprecated
    { NO_OPT_RASPSreset, (char*)"-raspsini:", SO_NONE },
    { OPT_RASPSlds, (char*)"-raspslds", SO_OPT },

    { OPT_deltaUbAbsolute, (char*)"-agap", SO_REQ_SEP },
    { OPT_deltaUbRelative, (char*)"-rgap", SO_OPT },
    { NO_OPT_trws, (char*)"-trws:", SO_NONE },
    { OPT_trwsAccuracy, (char*)"-trws", SO_OPT },
    { OPT_trwsAccuracy, (char*)"--trws-accuracy", SO_REQ_SEP },
    { OPT_trwsOrder, (char*)"--trws-order", SO_NONE },
    { NO_OPT_trwsOrder, (char*)"--trws-order:", SO_NONE },
    { OPT_trwsNIter, (char*)"--trws-n-iters", SO_REQ_SEP },
    { OPT_trwsNIterNoChange, (char*)"--trws-n-iters-no-change", SO_REQ_SEP },
    { OPT_trwsNIterComputeUb, (char*)"--trws-n-iters-compute-ub", SO_REQ_SEP },

    //preprocessing
    { OPT_minsumDiffusion, (char*)"-M", SO_REQ_SEP },
    { OPT_singletonConsistency, (char*)"-S", SO_NONE },
    { OPT_preprocessTernary, (char*)"-t", SO_OPT },
    { NO_OPT_preprocessTernary, (char*)"-t:", SO_NONE },
    { OPT_preprocessFunctional, (char*)"-f", SO_OPT },
    { NO_OPT_preprocessFunctional, (char*)"-f:", SO_NONE },
    { OPT_preprocessNary, (char*)"-n", SO_OPT },
    { NO_OPT_preprocessNary, (char*)"-n:", SO_NONE },

    { OPT_QueueComplexity, (char*)"-o", SO_NONE },
    { OPT_MSTDAC, (char*)"-mst", SO_NONE },
    { NO_OPT_MSTDAC, (char*)"-mst:", SO_NONE },
    { OPT_DEE, (char*)"-dee", SO_OPT },
    { NO_OPT_DEE, (char*)"-dee:", SO_OPT },
    { OPT_lds, (char*)"-l", SO_OPT },
    { NO_OPT_lds, (char*)"-l:", SO_NONE },
    { OPT_restart, (char*)"-L", SO_OPT },
    { NO_OPT_restart, (char*)"-L:", SO_NONE },
    { OPT_btlimit, (char*)"-bt", SO_OPT },
    { NO_OPT_btlimit, (char*)"-bt:", SO_NONE },
    { OPT_hbfs, (char*)"-hbfs", SO_OPT },
    { OPT_hbfs, (char*)"-bfs", SO_OPT },
    { NO_OPT_hbfs, (char*)"-hbfs:", SO_NONE },
    { NO_OPT_hbfs, (char*)"-bfs:", SO_NONE },
    { OPT_open, (char*)"-open", SO_REQ_SEP },
    { OPT_localsearch, (char*)"-i", SO_OPT }, // incop option default or string for narycsp argument
    { OPT_EDAC, (char*)"-k", SO_REQ_SEP },
    { OPT_ub, (char*)"-ub", SO_REQ_SEP }, // init upper bound in cli
    { OPT_divDist, (char*)"-div", SO_REQ_SEP }, // distance between solutions
    { OPT_divWidth, (char*)"-divwidth", SO_REQ_SEP }, // max relaxed MDD width
    { OPT_divWidth, (char*)"-mdd", SO_REQ_SEP }, // max relaxed MDD width
    { OPT_divRelax, (char*)"-divrelax", SO_REQ_SEP }, // relaxation method
    { OPT_divRelax, (char*)"-mddh", SO_REQ_SEP }, // relaxation method
    { OPT_divMethod, (char*)"-divmethod", SO_REQ_SEP }, // encoding of diversity constraint
    { OPT_divMethod, (char*)"-divm", SO_REQ_SEP }, // encoding of diversity constraint
    // MENDELSOFT
    { OPT_generation, (char*)"-g", SO_NONE }, //sort pedigree by increasing generation number and if equal by increasing individual number
    //	{ OPT_pedigree_by_MPE,  		(char*) "-y", 				SO_OPT			}, // bayesian flag
    { MENDEL_OPT_genotypingErrorRate, (char*)"-genoError", SO_REQ_SEP },
    { MENDEL_OPT_resolution, (char*)"-precision", SO_REQ_SEP },
    { OPT_pedigree_by_MPE, (char*)"-bayes", SO_NONE }, // bayesian flag
    { MENDEL_OPT_EQUAL_FREQ, (char*)"-pequal", SO_NONE }, // allocate equal frequencies to all allele
    { MENDEL_OPT_ESTIMAT_FREQ, (char*)"-probdata", SO_NONE }, //  probs depending on the frequencies found in the current pedigree problem
    { MENDEL_OPT_ALLOCATE_FREQ, (char*)"-problist", SO_MULTI }, // read probability distribution from command line

    { OPT_Z, (char*)"-logz", SO_NONE }, // compute log partition function (log Z)
    { OPT_epsilon, (char*)"-epsilon", SO_REQ_SEP }, // approximation parameter for computing Z

    { OPT_learning, (char*)"-learning", SO_NONE }, // pseudoboolean learning during search
#ifndef NDEBUG
    { OPT_verifyopt, (char*)"-opt", SO_NONE }, // for debugging purposes, checks the given optimal solution (problem.sol) is not pruned during search
#endif
    { OPT_timer, (char*)"-timer", SO_REQ_SEP }, // CPU timer

    // random generator
    { OPT_seed, (char*)"-seed", SO_REQ_SEP },
    { OPT_random, (char*)"-random", SO_REQ_SEP }, // init upper bound in cli

// VNS Methods
#ifdef BOOST
    { OPT_VNS_search, (char*)"-vns", SO_NONE },
    { OPT_VNS_search, (char*)"--vns", SO_NONE },
    { OPT_VNS_search, (char*)"-dgvns", SO_NONE },
    { OPT_VNS_search, (char*)"--dgvns", SO_NONE },
#ifdef OPENMPI
    { OPT_CPDGVNS_search, (char*)"--cpdgvns", SO_NONE },
    { OPT_RADGVNS_search, (char*)"-radgvns", SO_NONE },
    { OPT_RADGVNS_search, (char*)"--radgvns", SO_NONE },
    { OPT_RSDGVNS_search, (char*)"--rsdgvns", SO_NONE },
    { OPT_plimit, (char*)"--plimit", SO_NONE },
#endif
    { OPT_TREEDEC_search, (char*)"--treedec", SO_NONE },
    { OPT_decfile, (char*)"--decfile", SO_REQ_SEP },
    { OPT_vns_output, (char*)"--foutput", SO_REQ_SEP },
    { OPT_vns_sol_init, (char*)"-vnsini", SO_REQ_SEP },
    { OPT_vns_sol_init, (char*)"-ldsini", SO_REQ_SEP },
    { OPT_vns_sol_init, (char*)"--solution-init", SO_REQ_SEP },
    { OPT_lds_min, (char*)"-ldsmin", SO_REQ_SEP },
    { OPT_lds_max, (char*)"-ldsmax", SO_REQ_SEP },
    { OPT_lds_inc, (char*)"-ldsinc", SO_REQ_SEP },
    { OPT_k_min, (char*)"-kmin", SO_REQ_SEP },
    { OPT_k_min, (char*)"--kinit", SO_REQ_SEP },
    { OPT_k_max, (char*)"-kmax", SO_REQ_SEP },
    { OPT_k_max, (char*)"--kmax", SO_REQ_SEP },
    { OPT_k_inc, (char*)"-kinc", SO_REQ_SEP },
    //    { OPT_vns_var_heur, (char*) "--variable-heuristic", SO_REQ_SEP },
    { OPT_neighbor_change, (char*)"--strategy", SO_NONE },
    { OPT_neighbor_synch, (char*)"--synch", SO_NONE },
    { OPT_optimum, (char*)"-best", SO_REQ_SEP },
    { OPT_optimum, (char*)"--best", SO_REQ_SEP },
#endif
    SO_END_OF_OPTIONS
};

void ShowFiles(int argc, TCHAR** argv)
{
    // glob files to catch expand wildcards
    CSimpleGlob glob(SG_GLOB_NODOT | SG_GLOB_NOCHECK);
    if (SG_SUCCESS != glob.Add(argc, argv)) {
        _tprintf(_T("Error while globbing files\n"));
        return;
    }

    for (int n = 0; n < glob.FileCount(); ++n) {
        _tprintf(_T("file %2d: '%s'\n"), n, glob.File(n));
    }
}

static const TCHAR* GetLastErrorText(int a_nError)
{
    switch (a_nError) {
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

static void Pedi_Args(CSimpleOpt& args, int nMultiArgs)
{
    TCHAR** rgpszArg = NULL;
    int nMultiBackup = nMultiArgs;

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

    if (nMultiBackup < 0) {
        for (int n = 0; n < nMultiArgs; ++n) {
            if (ToulBar2::verbose >= 0)
                _tprintf(_T("MultiArg %d: %s\n"), n, rgpszArg[n]);
            float f2;
            sscanf(rgpszArg[n], "%f", &f2);
            ToulBar2::allelefreqdistrib.push_back(f2);
        }
        if (ToulBar2::verbose >= 0) {
            for (unsigned int x = 0; x <= sizeof(ToulBar2::allelefreqdistrib); x++)
                cout << "Allele prob used " << x << "=" << ToulBar2::allelefreqdistrib[x] << endl;
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

/* return current binary path extracted from argv[0] or from the env var $path */

char* find_bindir(const char* bin_name, char* buffer, size_t buflen)
{
    struct stat st;
    char *path, *tok;
    if (!stat(bin_name, &st)) {
        char* end = (char*)strrchr(bin_name, PATH_SEP_CHR);
        static char bin_path[512];
        if (end) {
            *end = 0;
            strncpy(buffer, bin_name, buflen);
            sprintf(bin_path, "%s%c", buffer, PATH_SEP_CHR);
        } else {
            strcpy(buffer, ".");
            //path separator added to the path value
            sprintf(bin_path, "%s%c", buffer, PATH_SEP_CHR);
        }
        return (bin_path);
    }
    path = strdup(getenv("PATH"));
    tok = strtok(path, PATH_DELIM);
    while (tok) {
        snprintf(buffer, buflen, "%s%c%s", tok, PATH_SEP_CHR, bin_name);
        if (!stat(buffer, &st)) {
            static char bin_path[512];
            strncpy(buffer, tok, buflen);
            free(path);
            sprintf(bin_path, "%s%c", buffer, PATH_SEP_CHR);
            return bin_path;
        }
        tok = strtok(NULL, PATH_DELIM);
    }
    free(path);
    buffer[0] = 0;
    return NULL;
}

//  current unused option letters: 	f F G H J K n N Q U W Y
void help_msg(char* toulbar2filename)
{
    cout << "*************************" << endl;
#ifdef MENDELSOFT
    cout << "* MendelSoft Help Message *" << endl;
#else
    cout << "* ToulBar2 Help Message *" << endl;
#endif
    cout << "*************************" << endl;
    cout << endl;
    cout << "Command line is:" << endl;
    cout << toulbar2filename << " problem_filename [options]" << endl;
    cout << endl;
#ifndef MENDELSOFT
    cout << "Available problem formats (specified by the filename extension) are:" << endl;
#ifdef BOOST
    cout << "   *.cfn : Cost Function Network format (see toulbar2 web site)" << endl;
#endif
    cout << "   *.wcsp : Weighted CSP format (see toulbar2 web site)" << endl;
    cout << "   *.wcnf : Weighted Partial Max-SAT format (see Max-SAT Evaluation)" << endl;
    cout << "   *.cnf : (Max-)SAT format" << endl;
    cout << "   *.qpbo : quadratic pseudo-Boolean optimization (unconstrained quadratic programming) format (see also option -qpmult)" << endl;
    cout << "   *.opb : pseudo-Boolean optimization format" << endl;
    cout << "   *.uai : Bayesian network and Markov Random Field format (see UAI'08 Evaluation) followed by an optional evidence filename (performs MPE task, see -logz for PR task, and write its solution in file .MPE or .PR using the same directory as toulbar2)" << endl;
    cout << "   *.LG : Bayesian network and Markov Random Field format using logarithms instead of probabilities" << endl;
#ifdef XMLFLAG
    cout << "   *.xml : CSP and weighted CSP in XML format XCSP 2.1 (constraints in extension only)";
#ifdef MAXCSP
    cout << " (Max-CSP only)";
#endif
    cout << endl;
#endif
    cout << "   *.pre : pedigree format (see doc/MendelSoft.txt for Mendelian error correction)" << endl;
    cout << "   *.pre *.map : pedigree and genetic map formats (see doc/HaplotypeHalfSib.txt for haplotype reconstruction in half-sib families)" << endl;
    cout << "   *.bep  : satellite scheduling format (CHOCO benchmark)" << endl
         << endl;
    cout << "   *.order  : variable elimination order" << endl;
    cout << "   *.cov  : tree decomposition given by a list of clusters in topological order of a rooted forest," << endl;
    cout << "      each line contains a cluster number, then a cluster parent number with -1 for the root(s) cluster(s), followed by a list of variable indexes" << endl;
    cout << "   *.dec  : a list of overlapping clusters without the running intersection property used by VNS-like methods," << endl;
    cout << "      each line contains a list of variable indexes" << endl;
    cout << "   *.sol  : initial solution for the problem (given as initial upperbound plus one and as default value heuristic, or only as initial upperbound if option -x: is added)" << endl
         << endl;
#ifdef BOOST
    cout << "Note: cfn, cnf, LG, qpbo, opb, uai, wcnf, wcsp formats can be read in gzip'd or xz compressed format, e.g., toulbar2 problem.cfn.xz" << endl;
#endif
    cout << "Warning! File formats are recognized by filename extensions. To change the default file format extension, use option --old_ext=\".new\" Examples: --cfn_ext='.json' --wcspgz_ext='.wgz' --sol_ext='.sol2'  " << endl;
    cout << endl;
#endif
    cout << "Available options are (use symbol \":\" after an option to remove a default option):" << endl;
    cout << "   -help : shows this help message" << endl;
    cout << "   -ub=[decimal] : initial problem upperbound (default value is " << MAX_COST << ")" << endl;
    cout << "   -agap=[decimal] : stop search if the absolute optimality gap reduces below the given value (provides guaranteed approximation) (default value is " << ToulBar2::deltaUbS << ")" << endl;
    cout << "   -rgap=[double] : stop search if the relative optimality gap reduces below the given value (provides guaranteed approximation) (default value is " << ToulBar2::deltaUbRelativeGap << ")" << endl;
    cout << "   -v=[integer] : verbosity level" << endl;
    cout << "   -s=[integer] : shows each solution found. 1 prints value numbers, 2 prints value names, 3 prints also variable names (default 1)" << endl;
#ifndef MENDELSOFT
    cout << "   -w=[filename] : writes last/all solutions in filename (or \"sol\" if no parameter is given)" << endl;
    cout << "   -w=[integer] : 1 writes value numbers, 2 writes value names, 3 writes also variable names (default 1)" << endl;
    cout << "   -precision=[integer] defines the number of digits that should be representable on probabilities in uai/pre files (default value is " << ToulBar2::resolution << ")" << endl;
    cout << "   -qpmult=[double] defines coefficient multiplier for quadratic terms (default value is " << ToulBar2::qpboQuadraticCoefMultiplier << ")" << endl;
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
    cout << "           -problist [nbre of prob] p1 p2 p3... : allele probability distribution given explicitely in the command line" << endl
         << endl;
#endif
#ifndef MENDELSOFT
#ifndef __WIN32__
    cout << "   -timer=[integer] : CPU time limit in seconds" << endl;
#endif
    cout << "   -bt=[integer] : limit on the number of backtracks (" << ToulBar2::backtrackLimit << " by default)" << endl;
    cout << "   -seed=[integer] : random seed non-negative value or use current time if a negative value is given (default value is " << ToulBar2::seed << ")" << endl;
    cout << "   --stdin=[format] : read file from pipe ; e.g., cat example.wcsp | toulbar2 --stdin=wcsp" << endl;
    cout << "   -var=[integer] : searches by branching only on the first -the given value- decision variables, assuming the remaining variables are intermediate variables completely assigned by the decision variables (use a zero if all variables are decision variables) (default value is " << ToulBar2::nbDecisionVars << ")" << endl;
    cout << "   -b : searches using binary branching always instead of binary branching for interval domains and n-ary branching for enumerated domains";
    if (ToulBar2::binaryBranching)
        cout << " (default option)";
    cout << endl;
    cout << "   -svo : searches using a static variable ordering heuristic (same order as DAC)";
    if (ToulBar2::Static_variable_ordering)
        cout << " (default option)";
    cout << endl;
    cout << "   -c : searches using binary branching with last conflict backjumping variable ordering heuristic";
    if (ToulBar2::lastConflict)
        cout << " (default option)";
    cout << endl;
    cout << "   -q=[integer] : weighted degree variable ordering heuristic if the number of cost functions is less than the given value (default value is " << ToulBar2::weightedDegree << ")" << endl;
    cout << "   -m=[integer] : variable ordering heuristic based on mean (m=1) or median (m=2) costs (in conjunction with weighted degree heuristic -q) (default value is " << ToulBar2::weightedTightness << ")" << endl;
    cout << "   -d=[integer] : searches using dichotomic branching (d=1 splitting in the middle of domain range, d=2 splitting in the middle of sorted unary costs) instead of binary branching when current domain size is strictly greater than " << ToulBar2::dichotomicBranchingSize << " (default value is " << ToulBar2::dichotomicBranching << ")" << endl;
    cout << "   -sortd : sorts domains based on increasing unary costs (warning! works only for binary WCSPs)";
    if (ToulBar2::sortDomains)
        cout << " (default option)";
    cout << endl;
    cout << "   -solr : solution-based phase saving";
    if (ToulBar2::solutionBasedPhaseSaving)
        cout << " (default option)";
    cout << endl;
    cout << "   -e=[integer] : boosting search with variable elimination of small degree (less than or equal to 3) (default value is " << ToulBar2::elimDegree << ")" << endl;
    cout << "   -p=[integer] : preprocessing only: general variable elimination of degree less than or equal to the given value (default value is " << ToulBar2::elimDegree_preprocessing << ")" << endl;
    cout << "   -t=[integer] : preprocessing only: simulates restricted path consistency by adding ternary cost functions on triangles of binary cost functions within a given maximum space limit (in MB)";
    if (ToulBar2::preprocessTernaryRPC)
        cout << " (" << ToulBar2::preprocessTernaryRPC << " MB)";
    cout << endl;
    cout << "   -f=[integer] : preprocessing only: variable elimination of functional (f=1) (resp. bijective (f=2)) variables (default value is " << ToulBar2::preprocessFunctional << ")" << endl;
    cout << "   -dec : preprocessing only: pairwise decomposition of cost functions with arity >=3 into smaller arity cost functions";
    if (ToulBar2::costfuncSeparate)
        cout << " (default option)";
    cout << endl;
    cout << "   -n=[integer] : preprocessing only: projects n-ary cost functions on all binary cost functions if n is lower than the given value (default value is " << ToulBar2::preprocessNary << ")" << endl;
#ifdef BOOST
    cout << "   -mst : maximum spanning tree DAC ordering";
    if (ToulBar2::MSTDAC)
        cout << " (default option)";
    cout << endl;
#endif
    cout << "   -nopre : removes all preprocessing options (equivalent to -e: -p: -t: -f: -dec: -n: -mst: -dee: -trws:)" << endl;
    cout << "   -o : ensures optimal worst-case time complexity of DAC and EAC (can be slower in practice)";
    if (ToulBar2::QueueComplexity)
        cout << " (default option)";
    cout << endl;
    cout << "   -k=[integer] : soft local consistency level (NC with Strong NIC for global cost functions=0, (G)AC=1, D(G)AC=2, FD(G)AC=3, (weak) ED(G)AC=4) (default value is " << ToulBar2::LcLevel << ")" << endl;
    cout << "   -dee=[integer] : restricted dead-end elimination (value pruning by dominance rule from EAC value (dee>=1 and dee<=3)) and soft neighborhood substitutability (in preprocessing (dee=2 or dee=4) or during search (dee=3)) (default value is " << ToulBar2::DEE << ")" << endl;
    cout << "   -l=[integer] : limited discrepancy search, use a negative value to stop the search after the given absolute number of discrepancies has been explored (discrepancy bound = " << maxdiscrepancy << " by default)";
    if (ToulBar2::lds)
        cout << " (default option)";
    cout << endl;
    cout << "   -L=[integer] : randomized (quasi-random variable ordering) search with restart (maximum number of nodes/VNS restarts = " << maxrestarts << " by default)";
    if (ToulBar2::restart >= 0)
        cout << " (default option)";
    cout << endl;
    cout << "   -i=[\"string\"] : initial upperbound found by INCOP local search solver." << endl;
    cout << "       string parameter is optional, using \"" << Incop_cmd << "\" by default with the following meaning:" << endl;
    cout << "       stoppinglowerbound randomseed nbiterations method nbmoves neighborhoodchoice neighborhoodchoice2 minnbneighbors maxnbneighbors neighborhoodchoice3 autotuning tracemode" << endl;
#ifdef BOOST
    cout << "   -vns : unified decomposition guided variable neighborhood search (a problem decomposition can be given as *.dec, *.cov, or *.order input files or using tree decomposition options such as -O)";
#ifdef OPENMPI
    //    cout << "   -cpdgvns : initial upperbound found by cooperative parallel DGVNS (usage: \"mpirun -n [NbOfProcess] toulbar2 -cpdgvns problem.wcsp\")" << endl;
    //    cout << "   -rsdgvns : initial upperbound found by replicated synchronous DGVNS (usage: \"mpirun -n [NbOfProcess] toulbar2 -rsdgvns problem.wcsp\")" << endl;
    cout << " (usage for parallel version: \"mpirun -n [NbOfProcess] toulbar2 -vns problem.wcsp\")";
#endif
    cout << endl;
    cout << "   -vnsini=[integer] : initial solution for VNS-like methods found (-1) at random, (-2) min domain values, (-3) max domain values, (-4) first solution found by a complete method, (k=0 or more) tree search with k discrepancy max (" << ToulBar2::vnsInitSol << " by default)" << endl;
    cout << "   -ldsmin=[integer] : minimum discrepancy for VNS-like methods (" << ToulBar2::vnsLDSmin << " by default)" << endl;
    cout << "   -ldsmax=[integer] : maximum discrepancy for VNS-like methods (number of problem variables multiplied by maximum domain size -1 by default)" << endl;
    cout << "   -ldsinc=[integer] : discrepancy increment strategy for VNS-like methods using (1) Add1, (2) Mult2, (3) Luby operator (" << ToulBar2::vnsLDSinc << " by default)" << endl;
    cout << "   -kmin=[integer] : minimum neighborhood size for VNS-like methods (" << ToulBar2::vnsKmin << " by default)" << endl;
    cout << "   -kmax=[integer] : maximum neighborhood size for VNS-like methods (number of problem variables by default)" << endl;
    cout << "   -kinc=[integer] : neighborhood size increment strategy for VNS-like methods using (1) Add1, (2) Mult2, (3) Luby operator (4) Add1/Jump (" << ToulBar2::vnsKinc << " by default)" << endl;
    cout << "   -best=[integer] : stop VNS-like methods if a better solution is found (default value is " << ToulBar2::vnsOptimum << ")" << endl;
    cout << endl;
#endif
    cout << "   -z=[filename] : saves problem in wcsp (by default) or cfn format (see below) in filename (or \"problem.wcsp/.cfn\"  if no parameter is given)" << endl;
    cout << "                   writes also the  graphviz dot file  and the degree distribution of the input problem (wcsp format only)" << endl;
    cout << "   -z=[integer] : 1 or 3: saves original instance in 1-wcsp or 3-cfn format (1 by default), 2 or 4: saves after preprocessing in 2-wcsp or 4-cfn format (this option can be combined with the previous one)" << endl;
    cout << "   -Z=[integer] : debug mode (save problem at each node if verbosity option -v=num >= 1 and -Z=num >=3)" << endl;
#ifndef NDEBUG
    cout << "   -opt filename.sol : checks a given optimal solution (given as input filename with \".sol\" extension) is never pruned by propagation (works only if compiled with debug)" << endl;
#endif
    cout << "   -x=[(,i[=#<>]a)*] : performs an elementary operation ('=':assign, '#':remove, '<':decrease, '>':increase) with value a on variable of index i (multiple operations are separated by a comma and no space) (without any argument, a complete assignment -- used as initial upper bound and as value heuristic -- read from default file \"sol\" taken as a certificate or given as input filename with \".sol\" extension)" << endl;
    cout << endl;
    cout << "   -M=[integer] : preprocessing only: Min Sum Diffusion algorithm (default number of iterations is " << ToulBar2::minsumDiffusion << ")" << endl;
    cout << "   -A=[integer] : enforces VAC at each search node with a search depth less than the absolute value of a given value, if negative value then VAC is not performed inside depth-first search of hybrid best-first search (default value is " << ToulBar2::vac << ")" << endl;
    cout << "   -T=[decimal] : threshold cost value for VAC (default value is " << ToulBar2::costThreshold << ")" << endl;
    cout << "   -P=[decimal] : threshold cost value for VAC during the preprocessing phase (default value is " << ToulBar2::costThresholdPre << ")" << endl;
    cout << "   -C=[float] : multiplies all costs internally by this number when loading the problem (default value is " << ToulBar2::costMultiplier << ")" << endl;
    cout << "   -S : preprocessing only: performs singleton consistency (only in conjunction with option \"-A\")";
    if (ToulBar2::singletonConsistency)
        cout << " (default option)";
    cout << endl;
    cout << "   -V : VAC-based value ordering heuristic";
    if (ToulBar2::vacValueHeuristic)
        cout << " (default option)";
    cout << endl;
    cout << "   -vacint : VAC-integrality/Full-EAC variable ordering heuristic";
    if (ToulBar2::FullEAC)
        cout << " (default option)";
    cout << endl;
    cout << "   -vacthr : automatic threshold cost value selection for VAC during search";
    if (ToulBar2::VACthreshold)
        cout << " (default option)";
    cout << endl;
    cout << "   -rasps=[integer] : VAC-based upper bound probing heuristic (0: disable, >0: max. nb. of backtracks) (default value is " << ((ToulBar2::useRASPS) ? ToulBar2::RASPSnbBacktracks : 0) << ")" << endl;
    cout << "   -raspslds=[integer] : VAC-based upper bound probing heuristic using LDS instead of DFS (0: DFS, >0: max. discrepancy) (default value is " << ((ToulBar2::useRASPS > 1) ? (ToulBar2::useRASPS - 1) : 0) << ")" << endl;
    cout << "   -raspsdeg=[integer] : automatic threshold cost value selection for probing heuristic (default value is " << ToulBar2::RASPSangle << "°)" << endl;
    cout << "   -raspsini : reset weighted degree variable ordering heuristic after doing upper bound probing";
    if (ToulBar2::RASPSreset)
        cout << " (default option)";
    cout << endl;
    cout << "   -trws=[float] : enforces TRW-S in preprocessing until a given precision is reached (default value is " << ToulBar2::trwsAccuracy << ")" << endl;
    cout << "   --trws-order : replaces DAC order by Kolmogorov's TRW-S order";
    if (ToulBar2::trwsOrder)
        cout << " (default option)";
    cout << endl;
    cout << "   --trws-n-iters=[integer] : enforce at most N iterations of TRW-S (default value is " << ToulBar2::trwsNIter << ")" << endl;
    cout << "   --trws-n-iters-no-change=[integer] : stop TRW-S when N iterations did not change the lower bound up the given precision (default value is " << ToulBar2::trwsNIterNoChange << ", -1=never)" << endl;
    cout << "   --trws-n-iters-compute-ub=[integer] : compute UB every N steps in TRW-S (default value is " << ToulBar2::trwsNIterComputeUb << ")" << endl;
    cout << endl;

    cout << "   -B=[integer] : (0) DFBB, (1) BTD, (2) RDS-BTD, (3) RDS-BTD with path decomposition instead of tree decomposition (default value is " << ToulBar2::btdMode << ")" << endl;
    cout << "   -O=[filename] : reads a variable elimination order or directly a valid tree decomposition (given by a list of clusters in topological order of a rooted forest, each line contains a cluster number, " << endl;
    cout << "      followed by a cluster parent number with -1 for the root(s) cluster(s), followed by a list of variable indexes) from a file used for BTD-like and variable elimination methods, and also DAC ordering" << endl;
#ifdef BOOST
    cout << "   -O=[negative integer] : build a tree decomposition (if BTD-like and/or variable elimination methods are used) and also a compatible DAC ordering using" << endl;
    cout << "                           (-" << MAX_CARD << ") maximum cardinality search ordering, (-" << MIN_DEGREE << ") minimum degree ordering, (-" << MIN_FILL << ") minimum fill-in ordering," << endl;
    cout << "                           (-" << ELIM_MST << ") maximum spanning tree ordering (see -mst), (-" << CUTHILL_MCKEE << ") reverse Cuthill-Mckee ordering, (-" << APPROX_MIN_DEGREE << ") approximate minimum degree ordering," << endl;
    cout << "                           (-" << ELIM_FILE_ORDER << ") default file ordering (the same if this option is missing, i.e. use the variable order in which variables appear in the problem file)" << endl;
#endif
    cout << "   -j=[integer] : splits large clusters into a chain of smaller embedded clusters with a number of proper variables less than this number" << endl;
    cout << "                (use options \"-B=3 -j=1 -svo -k=1\" for pure RDS, use value 0 for no splitting) (default value is " << ToulBar2::splitClusterMaxSize << ")" << endl;
    cout << "   -r=[integer] : limit on maximum cluster separator size (merge cluster with its father otherwise, use a negative value for no limit) (default value is " << ToulBar2::maxSeparatorSize << ")" << endl;
    cout << "   -X=[integer] : limit on minimum number of proper variables in a cluster (merge cluster with its father otherwise, use a zero for no limit) (default value is " << ToulBar2::minProperVarSize << ")" << endl;
    cout << "   -E=[float] : merges leaf clusters with their fathers if small local treewidth (in conjunction with option \"-e\" and positive threshold value) or ratio of number of separator variables by number of cluster variables above a given threshold (in conjunction with option \"-vns\") (default value is " << ToulBar2::boostingBTD << ")" << endl;
    cout << "   -R=[integer] : choice for a specific root cluster number" << endl;
    cout << "   -I=[integer] : choice for solving only a particular rooted cluster subtree (with RDS-BTD only)" << endl
         << endl;
    cout << "   -a=[integer] : finds at most a given number of solutions with a cost strictly lower than the initial upper bound and stops, or if no integer is given, finds all solutions (or counts the number of zero-cost satisfiable solutions in conjunction with BTD)";
    if (ToulBar2::allSolutions)
        cout << " (default value is " << ToulBar2::allSolutions << ")";
    cout << endl;
    cout << "   -div=[integer] : minimum Hamming distance between diverse solutions (use in conjunction with -a=integer with a limit of " << maxdivnbsol << " solutions) (default value is " << ToulBar2::divBound << ")" << endl;
    cout << "   -divm=[integer] : diversity encoding method: 0:Dual 1:Hidden 2:Ternary (default value is " << ToulBar2::divMethod << ")" << endl;
    cout << "   -mdd=[integer] : maximum relaxed MDD width for diverse solution global constraint (default value is " << ToulBar2::divWidth << ")" << endl;
    cout << "   -mddh=[integer] : MDD relaxation heuristic: 0: random, 1: high div, 2: small div, 3: high unary costs (default value is " << ToulBar2::divRelax << ")" << endl;
    cout << "   -D : approximate satisfiable solution count with BTD";
    if (ToulBar2::approximateCountingBTD)
        cout << " (default option)";
    cout << endl;
    cout << "   -logz : computes log of probability of evidence (i.e. log partition function or log(Z) or PR task) for graphical models only (problem file extension .uai)" << endl;
    cout << "   -epsilon=[float] : approximation factor for computing the partition function (greater than 1, default value is " << Exp(-ToulBar2::logepsilon) << ")" << endl;
    cout << endl;
    cout << "   -hbfs=[integer] : hybrid best-first search, restarting from the root after a given number of backtracks (default value is " << hbfsgloballimit << ")" << endl;
    cout << "   -open=[integer] : hybrid best-first search limit on the number of open nodes (default value is " << ToulBar2::hbfsOpenNodeLimit << ")" << endl;
    cout << "---------------------------" << endl;
    cout << "Alternatively one can call the random problem generator with the following options: " << endl;
    cout << endl;
    cout << "   -random=[bench profile]  : bench profile must be specified as follow :" << endl;
    cout << "                         n and d are respectively the number of variable and the maximum domain size  of the random problem." << endl;
    cout << "			" << endl;
    cout << "       bin-{n}-{d}-{t1}-{p2}-{seed}       :t1 is the tightness in percentage %of random binary cost functions" << endl;
    cout << "                                          :p2 is the num of binary cost functions to include" << endl;
    cout << "                                          :the seed parameter is optional (and will overwrite -seed)" << endl;

    cout << "   or:                                                                               " << endl;
    cout << "       binsub-{n}-{d}-{t1}-{p2}-{p3}-{seed} binary random & submodular cost functions" << endl;
    cout << "                                          t1 is the tightness in percentage % of random cost functions" << endl;
    cout << "                                          p2 is the num of binary cost functions to include" << endl;
    cout << "                                          p3 is the percentage % of submodular cost functions among p2 cost functions" << endl;
    cout << "                                           (plus 10 permutations of two randomly-chosen values for each domain)" << endl;
    cout << " or:                                                                               " << endl;
    cout << "      tern-{n}-{d}-{t1}-{p2}-{p3}-{seed}  p3 is the num of ternary cost functions" << endl;
    cout << " or:                                                                               " << endl;
    cout << "      nary-{n}-{d}-{t1}-{p2}-{p3}...-{pn}-{seed}  pn is the num of n-ary cost functions" << endl;
    cout << " or:                                                                               " << endl;
    cout << "      salldiff-{n}-{d}-{t1}-{p2}-{p3}...-{pn}-{seed}  pn is the num of salldiff global cost functions (p2 and p3 still being used for the number of random binary and ternary cost functions)" << endl;
    cout << "---------------------------" << endl;
    cout << "			" << endl;

    cout << endl;
#endif
}

int _tmain(int argc, TCHAR* argv[])
{
    //#pragma STDC FENV_ACCESS ON
    std::fesetround(FE_TONEAREST);

#ifdef OPENMPI
    MPIEnv env0;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &env0.ntasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &env0.myrank);
#endif
    tb2init();
#ifdef OPENMPI
    if (env0.myrank != 0)
        ToulBar2::verbose = -1;
#endif

    setlocale(LC_ALL, "C");
    bool certificate = false;
    char* certificateFilename = NULL;
    char* certificateString = NULL;
    char* solutionFileName = NULL;
    char buf[512];
    char* CurrentBinaryPath = find_bindir(argv[0], buf, 512); // current binary path search
    int timeout = 0;
    bool updateValueHeuristic = true;

    // Configuration for MaxSAT Evaluation
    //	ToulBar2::maxsateval = true;
    //	ToulBar2::verbose = -1;
    //	ToulBar2::binaryBranching = false;
    //	ToulBar2::lds = 1;

    // Configuration for UAI Evaluation
    // ToulBar2::uaieval = (env0.myrank == 0);
    //  ToulBar2::verbose = 0;
    //  ToulBar2::lds = 1;
    //  ToulBar2::incop_cmd = "0 1 3 idwa 100000 cv v 0 200 1 0 0";

    char* random_desc = NULL; // benchmark description set from command line;

    //default file extension : can be enforced using --foo_ext option in command line

    std::map<std::string, string> file_extension_map;
    file_extension_map["wcsp_ext"] = ".wcsp";
    file_extension_map["wcspgz_ext"] = ".wcsp.gz";
    file_extension_map["wcspxz_ext"] = ".wcsp.xz";
    file_extension_map["cfn_ext"] = ".cfn";
    file_extension_map["cfngz_ext"] = ".cfn.gz";
    file_extension_map["cfnxz_ext"] = ".cfn.xz";
    file_extension_map["wcspXML_ext"] = ".xml";
    file_extension_map["order_ext"] = ".order";
    file_extension_map["ub_ext"] = ".ub";
    file_extension_map["sol_ext"] = ".sol";
    file_extension_map["uai_ext"] = ".uai";
    file_extension_map["uaigz_ext"] = ".uai.gz";
    file_extension_map["uaixz_ext"] = ".uai.xz";
    file_extension_map["uai_log_ext"] = ".LG";
    file_extension_map["uaigz_log_ext"] = ".LG.gz";
    file_extension_map["uaixz_log_ext"] = ".LG.xz";
    file_extension_map["evid_ext"] = ".evid";
    file_extension_map["bep_ext"] = ".bep";
    file_extension_map["pre_ext"] = ".pre";
    file_extension_map["map_ext"] = ".map";
    file_extension_map["wcnf_ext"] = ".wcnf";
    file_extension_map["wcnfgz_ext"] = ".wcnf.gz";
    file_extension_map["wcnfxz_ext"] = ".wcnf.xz";
    file_extension_map["cnf_ext"] = ".cnf";
    file_extension_map["cnfgz_ext"] = ".cnf.gz";
    file_extension_map["cnfxz_ext"] = ".cnf.xz";
    file_extension_map["qpbo_ext"] = ".qpbo";
    file_extension_map["qpbogz_ext"] = ".qpbo.gz";
    file_extension_map["qpboxz_ext"] = ".qpbo.xz";
    file_extension_map["opb_ext"] = ".opb";
    file_extension_map["opbgz_ext"] = ".opb.gz";
    file_extension_map["opbxz_ext"] = ".opb.xz";
    file_extension_map["treedec_ext"] = ".cov";
    file_extension_map["clusterdec_ext"] = ".dec";

    assert(cout << "Warning! toulbar2 was compiled in debug mode and it can be very slow..." << endl);
    if (ToulBar2::verbose >= 0) {
        cout << "c " << CurrentBinaryPath;
#ifdef MENDELSOFT
        cout << "mendelsoft";
#else
        cout << "toulbar2";
#endif
        cout << "  version : " << ToulBar2::version << ", copyright (c) 2006-2020, toulbar2 team" << endl;
    }

    // --------------------------simple opt ----------------------

    // declare our options parser, pass in the arguments from main
    // as well as our array of valid options.
    CSimpleOpt args(argc, argv, g_rgOptions);

    while (args.Next()) {

        if (args.LastError() == SO_SUCCESS) {

            //if (check_file_ext(to_string(args.OptionText()),"_ext") )
            if (strstr(args.OptionText(), "_ext")) {
                string force_extension = args.OptionText();
                force_extension.erase(force_extension.begin(), force_extension.begin() + 2);
                if (ToulBar2::debug)
                    cout << "Old extension " << file_extension_map[force_extension] << " --> ";
                //			cout << " extension " << force_extension << " forced with " << args.OptionArg() << endl;
                file_extension_map[force_extension] = args.OptionArg();

                if (ToulBar2::debug)
                    cout << " New extension " << file_extension_map[force_extension] << endl;
            }

            //  search algorithm
            if (args.OptionId() == OPT_SEARCH_METHOD) {
                int mode = 0;
                if (args.OptionArg() != NULL)
                    mode = atoi(args.OptionArg());
                if (mode >= 0)
                    ToulBar2::btdMode = mode;
                else
                    ToulBar2::btdMode = 0;
                if (ToulBar2::debug)
                    cout << "Search Method used =  " << mode << endl;
            }

// VNS
#ifdef BOOST
            if (args.OptionId() == OPT_VNS_search) {
                //                ToulBar2::searchMethod = VNS;
                //                ToulBar2::vnsNeighborVarHeur = RANDOMVAR;
                ToulBar2::lds = maxdiscrepancy;
                ToulBar2::restart = maxrestarts;
#ifdef OPENMPI
                if (env0.ntasks > 1) {
                    ToulBar2::searchMethod = RPDGVNS;
                    ToulBar2::vnsParallel = true;
                    ToulBar2::vnsNeighborVarHeur = MASTERCLUSTERRAND;
                    ToulBar2::vnsParallelSync = false;
                } else {
                    ToulBar2::searchMethod = DGVNS;
                    ToulBar2::vnsNeighborVarHeur = CLUSTERRAND;
                }
#else
                ToulBar2::searchMethod = DGVNS;
                ToulBar2::vnsNeighborVarHeur = CLUSTERRAND;
#endif
            }
#ifdef OPENMPI
            if (args.OptionId() == OPT_CPDGVNS_search) {
                ToulBar2::searchMethod = CPDGVNS;
                ToulBar2::vnsParallel = true;
                ToulBar2::vnsNeighborVarHeur = MASTERCLUSTERRAND;
            }
            if (args.OptionId() == OPT_RADGVNS_search) {
                ToulBar2::searchMethod = RPDGVNS;
                ToulBar2::vnsParallel = true;
                ToulBar2::vnsNeighborVarHeur = MASTERCLUSTERRAND;
                ToulBar2::vnsParallelSync = false;
            }
            if (args.OptionId() == OPT_RSDGVNS_search) {
                ToulBar2::searchMethod = RPDGVNS;
                ToulBar2::vnsParallel = true;
                ToulBar2::vnsNeighborVarHeur = MASTERCLUSTERRAND;
                ToulBar2::vnsParallelSync = true;
            }
#endif
            if (args.OptionId() == OPT_TREEDEC_search) {
                ToulBar2::searchMethod = TREEDEC;
            }
            if (args.OptionId() == OPT_decfile) {
                ToulBar2::clusterFile = args.OptionArg();
                ifstream decfile(ToulBar2::clusterFile.c_str());
                if (!decfile) {
                    cerr << "File " << ToulBar2::clusterFile << " not found!" << endl;
                    exit(EXIT_FAILURE);
                }
            }

            if (args.OptionId() == OPT_vns_output) {
#ifdef OPENMPI
                if (env0.myrank == 0) {
#endif
                    ToulBar2::vnsOutput.clear();
                    ToulBar2::vnsOutput.open(args.OptionArg(),
                        ios::out | ios::trunc);
                    if (!ToulBar2::vnsOutput) {
                        cerr << "File " << args.OptionArg() << " cannot be open!" << endl;
#ifdef OPENMPI
                        for (int rank = 1; rank < env0.ntasks; ++rank) {
                            MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
                        }
#endif
                        exit(EXIT_FAILURE);
                    }
#ifdef OPENMPI
                }
#endif
            }
            if (args.OptionId() == OPT_vns_sol_init) {
                if (args.OptionArg() != NULL)
                    ToulBar2::vnsInitSol = static_cast<VNSSolutionInitMethod>(atoi(args.OptionArg()));
            }
            if (args.OptionId() == OPT_lds_min) {
                if (args.OptionArg() != NULL)
                    ToulBar2::vnsLDSmin = atoi(args.OptionArg());
            }
            if (args.OptionId() == OPT_lds_max) {
                if (args.OptionArg() != NULL)
                    ToulBar2::vnsLDSmax = atoi(args.OptionArg());
            }
            if (args.OptionId() == OPT_lds_inc) {
                if (args.OptionArg() != NULL)
                    ToulBar2::vnsLDSinc = static_cast<VNSInc>(atoi(args.OptionArg()));
            }
            if (args.OptionId() == OPT_k_min) {
                if (args.OptionArg() != NULL)
                    ToulBar2::vnsKmin = atoi(args.OptionArg());
            }
            if (args.OptionId() == OPT_k_max) {
                if (args.OptionArg() != NULL)
                    ToulBar2::vnsKmax = atoi(args.OptionArg());
            }
            if (args.OptionId() == OPT_k_inc) {
                if (args.OptionArg() != NULL)
                    ToulBar2::vnsKinc = static_cast<VNSInc>(atoi(args.OptionArg()));
            }
            //            if (args.OptionId() == OPT_vns_restart_strategy) {
            //                if (args.OptionArg() != NULL) {
            //                    string type = args.OptionArg();
            //                    if (type.compare("0") == 0) {
            //                        ToulBar2::vnsRestart = VNS_NORESTART;
            //                    } else if (type.compare("1") == 0) {
            //                        ToulBar2::vnsRestart = VNS_RESET;
            //                    } else if (type.compare("2") == 0) {
            //                        ToulBar2::vnsRestart = VNS_RESET_INC;
            //                    } else if (type.compare("3") == 0) {
            //                        ToulBar2::vnsRestart = VNS_FULL_RESTART;
            //                    } else {
            //                        if (env0.myrank == 0) {
            //                            cout << "Error : No implementation found for the given strategy"
            //                                 << endl;
            //                            cout << "Program will exit" << endl;
            //                        }
            //                        exit(EXIT_FAILURE);
            //                    }
            //                } else {
            //                    cout << "Warning : The strategy for local search method is NoRestart"
            //                         << endl;
            //                    cout << "specify for example --reset=2, to use ResetIncr strategy"
            //                         << endl;
            //                }
            //            }
            //            if (args.OptionId() == OPT_vns_var_heur) {
            //                if (args.OptionArg() != NULL)
            //                    ToulBar2::vnsNeighborVarHeur = static_cast<VNSVariableHeuristic>(atoi(args.OptionArg()));
            //            }
            if (args.OptionId() == OPT_neighbor_change) {
                ToulBar2::vnsNeighborChange = true;
            }
            if (args.OptionId() == OPT_neighbor_synch) {
                ToulBar2::vnsNeighborSizeSync = true;
            }
#ifdef OPENMPI
            if (args.OptionId() == OPT_plimit) {
                ToulBar2::vnsParallelLimit = true;
            }
#endif
            if (args.OptionId() == OPT_optimum) {
                if (args.OptionArg() != NULL)
                    //                    ToulBar2::vnsOptimum = atoll(args.OptionArg());
                    ToulBar2::vnsOptimumS = args.OptionArg();
            }
#endif

            if (args.OptionId() == OPT_stdin) {
                // stdin format reading by default stdin type is cfn format
                ToulBar2::stdin_format = args.OptionArg();
                if (ToulBar2::stdin_format.length() == 0) {
                    ToulBar2::stdin_format = "cfn";
                } else {
                    if (ToulBar2::stdin_format.compare("bep") == 0 || ToulBar2::stdin_format.compare("map") == 0 || ToulBar2::stdin_format.compare("pre") == 0) {
                        cerr << "Error: cannot read this " << ToulBar2::stdin_format << " format using stdin option!" << endl;
                        exit(EXIT_FAILURE);
                    }
                }
                //      cout << "pipe STDIN on waited FORMAT : " << ToulBar2::stdin_format<<endl;
            }

            // BTD root cluster
            if (args.OptionId() == OPT_btdRootCluster) {
                int root = atoi(args.OptionArg());
                if (root > 0)
                    ToulBar2::btdRootCluster = root;
            }
            // btd SubTree initialisation sub cluster

            if (args.OptionId() == OPT_btdSubTree) {
                int subcluster = atoi(args.OptionArg());
                if (subcluster >= 1)
                    ToulBar2::btdSubTree = subcluster;
            }
            //cluster Max size
            if (args.OptionId() == OPT_splitClusterMaxSize) {
                int cmaxsize = atoi(args.OptionArg());
                if (cmaxsize >= 1)
                    ToulBar2::splitClusterMaxSize = cmaxsize;
            }

            // E : merge leaf clusters with their fathers if small local treewidth
            if (args.OptionId() == OPT_boostingBTD) {
                ToulBar2::boostingBTD = 0.7;
                if (args.OptionArg() != NULL) {
                    double ratio = atof(args.OptionArg());
                    if (fabs(ratio) > 0. && fabs(ratio) <= 1.)
                        ToulBar2::boostingBTD = ratio;
                }
                if (ToulBar2::debug)
                    cout << "boostingBTD ON: " << ToulBar2::boostingBTD << endl;
            } else if (args.OptionId() == NO_OPT_boostingBTD) {
                ToulBar2::boostingBTD = 0.;
                if (ToulBar2::debug)
                    cout << "boostingBTD OFF " << ToulBar2::boostingBTD << endl;
            }

            //   -r=[integer] : limit on maximum cluster separator size (merge cluster with its father otherwise, use a negative value for no limit)
            if (args.OptionId() == OPT_maxSeparatorSize) {
                int sepmaxsize = atoi(args.OptionArg());
                if (sepmaxsize >= -1)
                    ToulBar2::maxSeparatorSize = sepmaxsize;
            }
            // minimal number of Proper variable included into each cluster of BTD
            if (args.OptionId() == OPT_minProperVarSize) {
                int minpvarsize = atoi(args.OptionArg());
                if (minpvarsize >= 0)
                    ToulBar2::minProperVarSize = minpvarsize;
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

            if (args.OptionId() == OPT_varOrder) {
                int varElimOrder = atoi(args.OptionArg());
                if (varElimOrder >= 0) {
                    char buf[512];
                    sprintf(buf, "%s", args.OptionArg());
                    //				if (ToulBar2::varOrder) delete [] ToulBar2::varOrder;
                    ToulBar2::varOrder = new char[strlen(buf) + 1];
                    sprintf(ToulBar2::varOrder, "%s", buf);
                    if (ToulBar2::debug)
                        cout << "variable order read from file " << args.OptionArg() << endl;
                } else {
                    ToulBar2::varOrder = reinterpret_cast<char*>(-varElimOrder);
                }
            }

            if (args.OptionId() == OPT_problemsaved_filename) {
                char buf[512];
                sprintf(buf, "%s", args.OptionArg());
                ToulBar2::problemsaved_filename = to_string(buf);
                //                if (!ToulBar2::dumpWCSP) ToulBar2::dumpWCSP = 1;
                if (ToulBar2::debug)
                    cout << "saved problem into file " << ToulBar2::problemsaved_filename << endl;
            }

            // filename of solution
            if (args.OptionId() == OPT_PARTIAL_ASSIGNMENT) {
                if (args.OptionArg() != NULL) {
                    certificate = true;
                    certificateString = args.OptionArg();
                    if (ToulBar2::debug)
                        cout << "partial assignment to be checked ..." << certificateString << endl;
                } else {
                    certificate = true;
                    certificateFilename = (char*)"sol";
                    if (ToulBar2::debug)
                        cout << "certificate of solution read in file: ./" << certificateFilename << endl;
                }
            } else if (args.OptionId() == NO_OPT_PARTIAL_ASSIGNMENT) {
                updateValueHeuristic = false;
            }

            // show Solutions
            if (args.OptionId() == OPT_showSolutions) {
                if (args.OptionArg() != NULL) {
                    int showType = atoi(args.OptionArg());
                    if (showType > 0 && showType < 4)
                        ToulBar2::showSolutions = showType;
                } else
                    ToulBar2::showSolutions = 1;
            }

            //#############################################
            if (args.OptionId() == OPT_writeSolution) {
                if (!ToulBar2::writeSolution)
                    ToulBar2::writeSolution = 1;
                if (solutionFileName == NULL)
                    solutionFileName = (char*)"sol";

                if (args.OptionArg() != NULL) {
                    char* tmpFile = new char[strlen(args.OptionArg()) + 1];
                    strcpy(tmpFile, args.OptionArg());
                    if (strlen(tmpFile) == 1 && (tmpFile[0] == '0' || tmpFile[0] == '1' || tmpFile[0] == '2' || tmpFile[0] == '3')) {
                        if (atoi(tmpFile) <= 2)
                            ToulBar2::pedigreeCorrectionMode = atoi(tmpFile);
                        if (atoi(tmpFile) >= 1)
                            ToulBar2::writeSolution = atoi(tmpFile);
                    } else {
                        solutionFileName = tmpFile;
                    }
                }
            }

            if (args.OptionId() == OPT_pedigreePenalty) {
                ToulBar2::pedigreePenalty = 1;
                int penaltyThreshold = atoi(args.OptionArg());
                if (penaltyThreshold >= 1)
                    ToulBar2::pedigreePenalty = penaltyThreshold;
            }
            // counting
            if (args.OptionId() == OPT_allSolutions) {
                ToulBar2::allSolutions = LONGLONG_MAX;
                if (args.OptionArg() != NULL) {
                    Long nbsol = atoll(args.OptionArg());
                    if (nbsol > 0)
                        ToulBar2::allSolutions = nbsol;
                }
            }

            // approximate counting

            if (args.OptionId() == OPT_approximateCountingBTD) {
                ToulBar2::approximateCountingBTD = true;
                ToulBar2::allSolutions = LONGLONG_MAX;
                ToulBar2::btdMode = 1;
            }

            if (args.OptionId() == OPT_binaryBranching) {
                ToulBar2::binaryBranching = true;
            }
            if (args.OptionId() == NO_OPT_binaryBranching) {
                ToulBar2::binaryBranching = false;
            }

            // static variable ordering

            if (args.OptionId() == OPT_Static_variable_ordering) {
                ToulBar2::Static_variable_ordering = true;
            } else if (args.OptionId() == NO_OPT_Static_variable_ordering) {
                ToulBar2::Static_variable_ordering = false;
            }

            // last conflict

            if (args.OptionId() == OPT_lastConflict) {
                ToulBar2::binaryBranching = true;
                ToulBar2::lastConflict = true;
            } else if (args.OptionId() == NO_OPT_lastConflict) {
                ToulBar2::lastConflict = false;
            }

            if (args.OptionId() == OPT_dichotomicBranching) {
                if (args.OptionArg() == NULL) {
                    ToulBar2::dichotomicBranching = 1;
                } else {
                    int dico = atoi(args.OptionArg());
                    if (dico > 0)
                        ToulBar2::dichotomicBranching = dico;
                    cout << "dichotomicBranching=" << ToulBar2::dichotomicBranching << endl;
                }
            } else if (args.OptionId() == NO_OPT_dichotomicBranching) {
                ToulBar2::dichotomicBranching = 0;
            }

            if (args.OptionId() == OPT_sortDomains) {
                ToulBar2::sortDomains = true;
            } else if (args.OptionId() == NO_OPT_sortDomains) {
                ToulBar2::sortDomains = false;
            }

            if (args.OptionId() == OPT_solutionBasedPhaseSaving) {
                ToulBar2::solutionBasedPhaseSaving = true;
            } else if (args.OptionId() == NO_OPT_solutionBasedPhaseSaving) {
                ToulBar2::solutionBasedPhaseSaving = false;
            }

            if (args.OptionId() == OPT_weightedTightness) {
                if (args.OptionArg() != NULL) {
                    int weightedtight = atol(args.OptionArg());
                    if (weightedtight >= 0)
                        ToulBar2::weightedTightness = weightedtight;
                } else {
                    ToulBar2::weightedTightness = 2;
                }
                if (ToulBar2::weightedTightness && !ToulBar2::weightedDegree)
                    ToulBar2::weightedDegree = 1000000;
            } else if (args.OptionId() == NO_OPT_weightedTightness) {
                ToulBar2::weightedTightness = 0;
            }

            // weitghted Degree (var ordering )
            if (args.OptionId() == OPT_weightedDegree and args.OptionArg() != NULL) {
                int weighteddegree = atol(args.OptionArg());
                if (weighteddegree > 0)
                    ToulBar2::weightedDegree = weighteddegree;
            } else if (args.OptionId() == NO_OPT_weightedDegree) {
                ToulBar2::weightedDegree = 0;
                ToulBar2::weightedTightness = 0;
                if (ToulBar2::debug)
                    cout << "ToulBar2::weightedDegree = false" << endl;
            }

            // LIMIT BRANCHING ON FIRST nbDecisionVars VARIABLES OPTION
            if (args.OptionId() == OPT_nbDecisionVars) {
                ToulBar2::nbDecisionVars = atoi(args.OptionArg());
            }

            //////////////////////////////////////////////////////

            if (args.OptionId() == NO_OPT_elimDegree) {
                ToulBar2::elimDegree = -1;
                if (ToulBar2::debug)
                    cout << "elimDegree OFF " << ToulBar2::elimDegree << endl;

            } else if (args.OptionId() == OPT_elimDegree and args.OptionArg() != NULL) {
                int ndegree = -1;
                ndegree = atol(args.OptionArg());
                if ((ndegree >= 0) && (ndegree <= 3)) {
                    ToulBar2::elimDegree = ndegree;
                }
                if (ToulBar2::debug)
                    cout << "elimDegree ON " << ToulBar2::elimDegree << endl;
            }
            ////////////////////////////////////////////

            // p[integer]: preprocessing only: general variable elimination of degree less than or equal to the given value
            // elimination degree in preprocessing
            if (args.OptionId() == OPT_elimDegree_preprocessing) {
                int ndegree;

                if (args.OptionArg() != NULL) {
                    ndegree = atoi(args.OptionArg());
                    if (ndegree != 0)
                        ToulBar2::elimDegree_preprocessing = ndegree;
                    if (ndegree < 0)
                        ToulBar2::elimSpaceMaxMB = 128;
                } else
                    ToulBar2::elimDegree_preprocessing = 3;
                if (ToulBar2::debug)
                    cout << "elimDegree_preprocessing ON: " << ToulBar2::elimDegree_preprocessing << endl;

            } else if (args.OptionId() == NO_OPT_elimDegree_preprocessing) {
                ToulBar2::elimDegree_preprocessing = -1;
                if (ToulBar2::debug)
                    cout << "elimDegree_preprocessing OFF: " << ToulBar2::elimDegree_preprocessing << endl;
            }

            // VAC PARAMETER
            if (args.OptionId() == OPT_minsumDiffusion)

            {
                if (!ToulBar2::vac)
                    ToulBar2::vac = 1;
                ToulBar2::minsumDiffusion = 1000;
                int nit = atoi(args.OptionArg());
                if (nit > 0)
                    ToulBar2::minsumDiffusion = nit;
            }

            if (args.OptionId() == OPT_vac) {
                ToulBar2::vac = 1;
                if (args.OptionArg() == NULL) {
                    ToulBar2::vac = 1;
                } else {
                    int depth = atoi(args.OptionArg());
                    if (depth != 0)
                        ToulBar2::vac = depth;
                }
                if (ToulBar2::debug)
                    cout << "VAC propagation ON" << endl;

            } else if (args.OptionId() == NO_OPT_vac) {
                if (ToulBar2::debug)
                    cout << "VAC propagation OFF" << endl;

                ToulBar2::vac = 0;
            }

            if (args.OptionId() == OPT_costThreshold) {
                //Cost ct = string2Cost(args.OptionArg());
                //if (ct > UNIT_COST)
                ToulBar2::costThresholdS = args.OptionArg();
            }

            if (args.OptionId() == OPT_costThresholdPre) {
                //Cost ct = string2Cost(args.OptionArg());
                //if (ct > UNIT_COST)
                ToulBar2::costThresholdPreS = args.OptionArg();
            }

            if (args.OptionId() == OPT_costMultiplier) {
                double co = atof(args.OptionArg());
                if (co > MIN_COST)
                    ToulBar2::costMultiplier = co;
            }

            if (args.OptionId() == OPT_qpbo_mult) {
                double co = atof(args.OptionArg());
                if (co != 0.)
                    ToulBar2::qpboQuadraticCoefMultiplier = co;
            }

            if (args.OptionId() == OPT_singletonConsistency)
                ToulBar2::singletonConsistency = true;
            if (args.OptionId() == OPT_vacValueHeuristic)
                ToulBar2::vacValueHeuristic = true;
            else if (args.OptionId() == NO_OPT_vacValueHeuristic)
                ToulBar2::vacValueHeuristic = false;
            if (args.OptionId() == OPT_preprocessTernary) {
                if (args.OptionArg() != NULL) {
                    int size = atol(args.OptionArg());
                    if (size >= 0)
                        ToulBar2::preprocessTernaryRPC = size;
                } else
                    ToulBar2::preprocessTernaryRPC = 128;
            } else if (args.OptionId() == NO_OPT_preprocessTernary) {
                if (ToulBar2::debug)
                    cout << "preprocess triangles of binary cost functions into ternary cost functions OFF" << endl;
                ToulBar2::preprocessTernaryRPC = 0;
            }

            if (args.OptionId() == OPT_trwsAccuracy) {
                if (args.OptionArg() == NULL) {
                    ToulBar2::trwsAccuracy = 0.00001;
                } else {
                    double co = atof(args.OptionArg());
                    if (co >= 0.)
                        ToulBar2::trwsAccuracy = co;
                    else
                        ToulBar2::trwsAccuracy = -1.;
                }
            } else if (args.OptionId() == NO_OPT_trws) {
                ToulBar2::trwsAccuracy = -1.;
            }
            if (args.OptionId() == OPT_trwsOrder) {
                ToulBar2::trwsOrder = true;
            } else if (args.OptionId() == NO_OPT_trwsOrder) {
                ToulBar2::trwsOrder = false;
            }
            if (args.OptionId() == OPT_trwsNIter) {
                ToulBar2::trwsNIter = atol(args.OptionArg());
            }
            if (args.OptionId() == OPT_trwsNIterNoChange) {
                ToulBar2::trwsNIterNoChange = atol(args.OptionArg());
            }
            if (args.OptionId() == OPT_trwsNIterComputeUb) {
                ToulBar2::trwsNIterComputeUb = atol(args.OptionArg());
            }

            // elimination of functional variables
            if (args.OptionId() == OPT_preprocessFunctional) {
                if (args.OptionArg() == NULL) {
                    ToulBar2::preprocessFunctional = 1;
                } else {
                    int func = atoi(args.OptionArg());
                    if (func > 0)
                        ToulBar2::preprocessFunctional = func;
                }
                if (ToulBar2::debug)
                    cout << "elimination of functional variables ON " << ToulBar2::preprocessFunctional << endl;

            } else if (args.OptionId() == NO_OPT_preprocessFunctional) {
                if (ToulBar2::debug)
                    cout << "elimination of functional variables OFF" << endl;
                ToulBar2::preprocessFunctional = 0;
            }

            if (args.OptionId() == OPT_costfuncSeparate) {
                if (ToulBar2::debug)
                    cout << "decomposition of cost functions" << endl;
                ToulBar2::costfuncSeparate = true;
            } else if (args.OptionId() == NO_OPT_costfuncSeparate) {
                if (ToulBar2::debug)
                    cout << "decomposition of cost functions OFF" << endl;
                ToulBar2::costfuncSeparate = false;
            }

            if (args.OptionId() == NO_OPT_DEE) {
                if (ToulBar2::debug)
                    cout << "dead-end elimination OFF" << endl;
                ToulBar2::DEE = 0;
            } else if (args.OptionId() == OPT_DEE) {
                ToulBar2::DEE = 1;
                if (args.OptionArg() != NULL) {
                    int dee = atoi(args.OptionArg());
                    if (dee >= 0)
                        ToulBar2::DEE = dee;
                }
                if (ToulBar2::debug)
                    cout << "dead-end elimination: " << ToulBar2::DEE << endl;
            }

            // pre projection of nary cost functions
            if (args.OptionId() == OPT_preprocessNary) {
                if (args.OptionArg() == NULL) {
                    ToulBar2::preprocessNary = 10;
                } else {
                    int maxnary = atoi(args.OptionArg());
                    if (maxnary > 0)
                        ToulBar2::preprocessNary = maxnary;
                }
                if (ToulBar2::debug)
                    cout << "preproject cost functions with arity lower than " << ToulBar2::preprocessNary << " ON" << endl;

            } else if (args.OptionId() == NO_OPT_preprocessNary) {
                if (ToulBar2::debug)
                    cout << "preproject of n-ary cost functions OFF" << endl;
                ToulBar2::preprocessNary = 0;
            }

            if (args.OptionId() == OPT_QueueComplexity)
                ToulBar2::QueueComplexity = true;
            if (args.OptionId() == OPT_MSTDAC)
                ToulBar2::MSTDAC = true;
            else if (args.OptionId() == NO_OPT_MSTDAC)
                ToulBar2::MSTDAC = false;

            // LDS
            if (args.OptionId() == OPT_lds) {
                string comment;
                if (args.OptionArg() == NULL) {
                    ToulBar2::lds = maxdiscrepancy;
                    comment = " (default value) ";
                } else {
                    int maxlds = atoi(args.OptionArg());
                    if (maxlds > 0 || maxlds < 0) {
                        ToulBar2::lds = maxlds;
                        ToulBar2::vnsLDSmin = abs(maxlds);
                        ToulBar2::vnsLDSmax = abs(maxlds);
                    }
                }
                if (ToulBar2::debug)
                    cout << "LDS ON #iter = " << ToulBar2::lds << comment << endl;

            } else if (args.OptionId() == NO_OPT_lds) {
                ToulBar2::lds = 0;
                if (ToulBar2::debug)
                    cout << "LDS OFF iter = " << ToulBar2::lds << endl;
            }

            // restart option
            if (args.OptionId() == OPT_restart) {
                string comment;
                if (args.OptionArg() == NULL) {
                    comment = " (default value) ";
                    ToulBar2::restart = maxrestarts;
                } else {
                    Long maxbt = atoll(args.OptionArg());
                    if (maxbt >= 0)
                        ToulBar2::restart = maxbt;
                }
                if (ToulBar2::debug)
                    cout << "restart ON #iter = " << ToulBar2::restart << comment << endl;
            } else if (args.OptionId() == NO_OPT_restart) {
                if (ToulBar2::debug)
                    cout << "restart OFF" << endl;
                ToulBar2::restart = -1;
            }

            // backtrack limit option
            if (args.OptionId() == OPT_btlimit) {
                string comment;
                if (args.OptionArg() == NULL) {
                    comment = " (default value) ";
                    ToulBar2::backtrackLimit = maxrestarts;
                } else {
                    Long maxbt = atoll(args.OptionArg());
                    if (maxbt >= 0)
                        ToulBar2::backtrackLimit = maxbt;
                }
                if (ToulBar2::debug)
                    cout << "backtrack limit ON #bt = " << ToulBar2::backtrackLimit << comment << endl;
            } else if (args.OptionId() == NO_OPT_btlimit) {
                if (ToulBar2::debug)
                    cout << "backtrack limit OFF" << endl;
                ToulBar2::backtrackLimit = LONGLONG_MAX;
            }

            // hybrid BFS option
            if (args.OptionId() == OPT_hbfs) {
                string comment;
                if (args.OptionArg() == NULL) {
                    comment = " (default value) ";
                    ToulBar2::hbfsGlobalLimit = hbfsgloballimit;
                } else {
                    Long maxbt = atoll(args.OptionArg());
                    if (maxbt > 0)
                        ToulBar2::hbfsGlobalLimit = maxbt;
                    else
                        ToulBar2::hbfsGlobalLimit = hbfsgloballimit;
                }
                ToulBar2::hbfs = 1; // initial value to perform a greedy search exploration before visiting a new open search node
                if (ToulBar2::debug)
                    cout << "hybrid BFS ON with global backtrack limit = " << ToulBar2::hbfsGlobalLimit << comment << endl;
            } else if (args.OptionId() == NO_OPT_hbfs) {
                if (ToulBar2::debug)
                    cout << "hybrid BFS OFF" << endl;
                ToulBar2::hbfs = 0;
                ToulBar2::hbfsGlobalLimit = 0;
            }
            if (args.OptionId() == OPT_open) {
                ToulBar2::hbfsOpenNodeLimit = OPEN_NODE_LIMIT;
                Long openlimit = atoll(args.OptionArg());
                if (openlimit > 0)
                    ToulBar2::hbfsOpenNodeLimit = openlimit;
                if (ToulBar2::hbfsGlobalLimit == 0)
                    ToulBar2::hbfsGlobalLimit = hbfsgloballimit;
                ToulBar2::hbfs = 1; // initial value to perform a greedy search exploration before visiting a new open search node
                if (ToulBar2::debug)
                    cout << "hybrid BFS ON with open node limit = " << ToulBar2::hbfsOpenNodeLimit << endl;
            }

            // local search INCOP
            if (args.OptionId() == OPT_localsearch) {
                if (args.OptionArg() != NULL) {
                    ToulBar2::incop_cmd = args.OptionArg();
                } else {
                    ToulBar2::incop_cmd = Incop_cmd;
                }
            }

            // EDAC OPTION
            if (args.OptionId() == OPT_EDAC) {
                ToulBar2::LcLevel = LC_EDAC;
                LcLevelType lclevel = (LcLevelType)atoi(args.OptionArg());
                if ((lclevel >= LC_NC) && (lclevel < LC_THEMAX))
                    ToulBar2::LcLevel = lclevel;
            }

            ////////////////////////////////////////////////////////////////////////////////////////
            ////////		MEDELSOFT command line processing 	////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////

            // g  sort pedigree by increasing generation number and if equal by increasing individual number

            if (args.OptionId() == OPT_generation) {
                ToulBar2::generation = true;
                if (ToulBar2::debug)
                    cout << "-g flag ==> sort pedigree option ON" << endl;
            }

            // bayesien -y [genotypinpErrorRate probabilityPrecision genotypePriorMode]  : pedigree solved by Bayesian MPE

            if (args.OptionId() == MENDEL_OPT_genotypingErrorRate) {
                if (args.OptionArg() != NULL) {
                    float f;
                    sscanf(args.OptionArg(), "%f", &f);
                    ToulBar2::errorg = f;
                    if (ToulBar2::debug)
                        cout << "New assignment for genotyping Error Rate = " << ToulBar2::errorg << endl;
                }
            }

            if (args.OptionId() == MENDEL_OPT_resolution) {
                if (args.OptionArg() != NULL) {
                    ToulBar2::resolution = atoi(args.OptionArg());
                    if (ToulBar2::debug)
                        cout << "New assignment for precision = " << ToulBar2::resolution << endl;
                }
            }

            if (args.OptionId() == OPT_pedigree_by_MPE) {
                ToulBar2::bayesian = true;
                if (ToulBar2::debug)
                    cout << endl;
                if (ToulBar2::debug)
                    cout << "<<< Bayesian mode ON >>>" << endl;

                if (ToulBar2::verbose >= 1) {
                    cout << endl;
                    cout << "genotypingErrorRate = " << ToulBar2::errorg << " (default value)" << endl;
                    cout << "resolution = " << ToulBar2::resolution << " (default value)" << endl;
                    cout << endl;
                }
            }

            if (args.OptionId() == MENDEL_OPT_EQUAL_FREQ && ToulBar2::bayesian) {

                ToulBar2::foundersprob_class = 0;
                if (ToulBar2::debug)
                    cout << "equal frequencies used ( default mode)" << ToulBar2::foundersprob_class << endl;

            } else if (args.OptionId() == MENDEL_OPT_ESTIMAT_FREQ && ToulBar2::bayesian) {
                ToulBar2::foundersprob_class = 1;
                if (ToulBar2::debug)
                    cout << " => prob depending on the frequencies found in the current pedigree probleme used" << endl;

            } else if (args.OptionId() == MENDEL_OPT_ALLOCATE_FREQ && ToulBar2::bayesian) {
                if (ToulBar2::debug)
                    cout << " => prob frequencies read from command line" << endl;
                ToulBar2::foundersprob_class = 2;
                Pedi_Args(args, -1);
            }

            ////////////////////////////////////////////////////////////////////////////////////////
            /////				DEBUG and verbose LEVEL  management + DUMP
            ////////////////////////////////////////////////////////////////////////////////////////
            //verbose mode
            //   v: verbosity level
            if (args.OptionId() == OPT_verbose) {
                if (args.OptionArg() != NULL) {
                    ToulBar2::verbose = atoi(args.OptionArg());
                } else {
                    ToulBar2::verbose = 0;
                }
                if (ToulBar2::debug)
                    cout << "verbose level = " << ToulBar2::verbose << endl;
            }

            //  z: save problem in wcsp/cfn format in filename \"problem.wcsp/cfn\" (1/3:before, 2/4:current problem after preprocessing)

            if (args.OptionId() == OPT_dumpWCSP) {
                if (args.OptionArg() != NULL) {
                    if (strlen(args.OptionArg()) == 1 && args.OptionArg()[0] >= '1' && args.OptionArg()[0] <= '4') {
                        ToulBar2::dumpWCSP = atoi(args.OptionArg());
                        if (ToulBar2::problemsaved_filename.rfind(".cfn") != string::npos && static_cast<ProblemFormat>((ToulBar2::dumpWCSP >> 1)+(ToulBar2::dumpWCSP & 1)) != CFN_FORMAT) {
                            cerr << "Error: filename extension .cfn not compatible with option -z=" << ToulBar2::dumpWCSP << endl;
                            exit(EXIT_FAILURE);
                        }
                        if (ToulBar2::problemsaved_filename.rfind(".wcsp") != string::npos && static_cast<ProblemFormat>((ToulBar2::dumpWCSP >> 1)+(ToulBar2::dumpWCSP & 1)) != WCSP_FORMAT) {
                            cerr << "Error: filename extension .wcsp not compatible with option -z=" << ToulBar2::dumpWCSP << endl;
                            exit(EXIT_FAILURE);
                        }
                    } else {
                        ToulBar2::problemsaved_filename = to_string(args.OptionArg());
                        if (ToulBar2::problemsaved_filename.rfind(".cfn") != string::npos && !ToulBar2::dumpWCSP)
                            ToulBar2::dumpWCSP = 3;
                        else
                            ToulBar2::dumpWCSP = 1;
                    }
                } else if (!ToulBar2::dumpWCSP)
                    ToulBar2::dumpWCSP = 1;

                if (ToulBar2::debug) {
                    cout << "Problem will be saved in " << ToulBar2::problemsaved_filename;
                    cout << "after " << ((ToulBar2::dumpWCSP % 2) ? "loading." : "preprocessing.") << endl;
                }
            }
            //   Z: debug mode (save problem at each node if verbosity option set!)
            //		for (int j=0; argv[i][j] != 0; j++) if (argv[i][j]=='Z') ToulBar2::debug++;

            if (args.OptionId() == OPT_debug) {
                if (args.OptionArg() != NULL) {
                    ToulBar2::debug = atoi(args.OptionArg());
                } else
                    ToulBar2::debug = 1;
                if (ToulBar2::debug)
                    cout << "debug level = " << ToulBar2::debug << endl;
            }

            // discrete integration for computing the partition function Z
            if (args.OptionId() == OPT_Z)
                ToulBar2::isZ = true;

            if (args.OptionId() == OPT_epsilon) {
                if (args.OptionArg() != NULL) {
                    ToulBar2::logepsilon = -Log(atof(args.OptionArg()));
                    if (ToulBar2::debug)
                        cout << "New assignment for epsilon = " << Exp(-ToulBar2::logepsilon) << endl;
                }
            }

            if (args.OptionId() == OPT_learning) {
                ToulBar2::learning = true;
            }

#ifndef NDEBUG
            if (args.OptionId() == OPT_verifyopt)
                ToulBar2::verifyOpt = true;
#endif

            // diversity
            if (args.OptionId() == OPT_divDist) {
                if (args.OptionArg() != NULL) {
                    ToulBar2::divBound = atoi(args.OptionArg());
                    ToulBar2::divNbSol = maxdivnbsol;
                    if (ToulBar2::debug)
                        cout << "Diversity distance = " << ToulBar2::divBound << endl;
                }
            }

            if (args.OptionId() == OPT_divWidth) {
                if (args.OptionArg() != NULL) {
                    ToulBar2::divWidth = atoi(args.OptionArg());
                    if (ToulBar2::debug)
                        cout << "Diversity MDD maximum width = " << ToulBar2::divWidth << endl;
                }
            }

            if (args.OptionId() == OPT_divMethod) {
                if (args.OptionArg() != NULL) {
                    ToulBar2::divMethod = atoi(args.OptionArg());
                    if (ToulBar2::debug)
                        cout << "Diversity method = " << ToulBar2::divMethod << endl;
                }
            }

            if (args.OptionId() == OPT_divRelax) {
                if (args.OptionArg() != NULL) {
                    ToulBar2::divRelax = atoi(args.OptionArg());
                    if (ToulBar2::debug)
                        cout << "Diversity MDD relaxation method = " << ToulBar2::divRelax << endl;
                }
            }

            // upper bound initialisation from command line
            if (args.OptionId() == OPT_ub) {
                ToulBar2::externalUB = args.OptionArg();
                rtrim(ToulBar2::externalUB);
            }

            if (args.OptionId() == OPT_deltaUbAbsolute) {
                ToulBar2::deltaUbS = args.OptionArg();
                rtrim(ToulBar2::deltaUbS);
            }

            if (args.OptionId() == OPT_deltaUbRelative) {
                ToulBar2::deltaUbRelativeGap = relativegap;
                if (args.OptionArg() != NULL) {
                    ToulBar2::deltaUbRelativeGap = atof(args.OptionArg());
                }
            }

            // CPU timer
            if (args.OptionId() == OPT_timer) {
                if (args.OptionArg() != NULL) {
                    timeout = atoi(args.OptionArg());
                    if (ToulBar2::debug)
                        cout << "TimeOut = " << timeout << endl;
                }
            }

            //////////RANDOM GENERATOR///////
            if (args.OptionId() == OPT_seed) {
                int seed = atoi(args.OptionArg());
                ToulBar2::seed = seed;
            }

            if (args.OptionId() == OPT_random) {

                random_desc = args.OptionArg();
            }

            if (args.OptionId() == OPT_nopre) {
                ToulBar2::elimDegree = -1;
                if (ToulBar2::debug)
                    cout << "elimDegree OFF " << ToulBar2::elimDegree << endl;
                ToulBar2::elimDegree_preprocessing = -1;
                if (ToulBar2::debug)
                    cout << "elimDegree_preprocessing OFF: " << ToulBar2::elimDegree_preprocessing << endl;
                if (ToulBar2::debug)
                    cout << "preprocess triangles of binary cost functions into ternary cost functions OFF" << endl;
                ToulBar2::preprocessTernaryRPC = 0;
                if (ToulBar2::debug)
                    cout << "elimination of functional variables OFF" << endl;
                ToulBar2::preprocessFunctional = 0;
                if (ToulBar2::debug)
                    cout << "preproject of n-ary cost functions OFF" << endl;
                ToulBar2::preprocessNary = 0;
                if (ToulBar2::debug)
                    cout << "decomposition of cost functions OFF" << endl;
                ToulBar2::costfuncSeparate = false;
                if (ToulBar2::debug)
                    cout << "maximum spanning tree DAC ordering OFF" << endl;
                ToulBar2::MSTDAC = false;
                if (ToulBar2::debug)
                    cout << "dead-end elimination OFF" << endl;
                ToulBar2::DEE = 0;
                if (ToulBar2::debug)
                    cout << "TRW-S OFF" << endl;
                ToulBar2::trwsAccuracy = -1.;
            }

            if (args.OptionId() == OPT_VACINT) {
                ToulBar2::FullEAC = true;
                if (ToulBar2::debug)
                    cout << "Strict AC ON" << endl;
            } else if (args.OptionId() == NO_OPT_VACINT) {
                ToulBar2::FullEAC = false;
                if (ToulBar2::debug)
                    cout << "Strict AC OFF" << endl;
            }

            if (args.OptionId() == OPT_VACthreshold) {
                ToulBar2::VACthreshold = true;
                if (ToulBar2::debug)
                    cout << "VAC iterations during search will go until the threshold calculated by RASPS approach." << endl;
            } else if (args.OptionId() == NO_OPT_VACthreshold) {
                ToulBar2::VACthreshold = false;
                if (ToulBar2::debug)
                    cout << "VAC automatic threshold OFF" << endl;
            }

            if (args.OptionId() == OPT_RASPS) {
                if (args.OptionArg() == NULL) {
                    if (ToulBar2::useRASPS == 0)
                        ToulBar2::useRASPS = 1;
                    ToulBar2::RASPSnbBacktracks = raspsbacktracks;
                } else {
                    Long limit = atoll(args.OptionArg());
                    if (ToulBar2::useRASPS == 0)
                        ToulBar2::useRASPS = 1;
                    ToulBar2::RASPSnbBacktracks = limit;
                }
                if (ToulBar2::debug)
                    cout << "RAPS ON " << ToulBar2::RASPSnbBacktracks << endl;
            } else if (args.OptionId() == NO_OPT_RASPS) {
                if (ToulBar2::debug)
                    cout << "RASPS OFF" << endl;
                ToulBar2::useRASPS = 0;
            }

            if (args.OptionId() == OPT_RASPSangle) {
                if (args.OptionArg() == NULL)
                    ToulBar2::RASPSangle = raspsangle;
                else {
                    int n = atoi(args.OptionArg());
                    ToulBar2::RASPSangle = n;
                }
                if (ToulBar2::debug)
                    cout << "RASPS angle set to " << ToulBar2::RASPSangle << "°" << endl;
            } else if (args.OptionId() == NO_OPT_RASPSangle) {
                if (ToulBar2::debug)
                    cout << "RASPS angle set to 90°" << endl;
                ToulBar2::RASPSangle = 90;
            }

            if (args.OptionId() == OPT_RASPSreset) {
                ToulBar2::RASPSreset = true;
                if (ToulBar2::debug)
                    cout << "RASPS reset ON" << endl;
            } else if (args.OptionId() == NO_OPT_RASPSreset) {
                ToulBar2::RASPSreset = false;
                if (ToulBar2::debug)
                    cout << "RASPS reset OFF" << endl;
            }

            if (args.OptionId() == OPT_RASPSlds) {
                if (args.OptionArg() == NULL)
                    ToulBar2::useRASPS = maxdiscrepancy + 1;
                else {
                    int lds = atoi(args.OptionArg());
                    if (lds > 0)
                        ToulBar2::useRASPS = lds + 1;
                }
                if (ToulBar2::debug)
                    cout << "RASPSlds" << ToulBar2::useRASPS << endl;
            }
        }

        else {
            cout << "<<<< ERROR >>>>  : " << endl;
            _tprintf(
                _T("%s: '%s' (use --help to get command line help)\n"),
                GetLastErrorText(args.LastError()), args.OptionText());

            _tprintf(_T("invalid Option :%s or Invalid argument: %s \n please check help using --help option "), args.OptionText(), args.OptionArg() ? args.OptionArg() : "");

            return 1;
        }
    }
    // process any files that were passed to us on the command line.
    // send them to the globber so that all wildcards are expanded
    // into valid filenames (e.g. *.cpp -> a.cpp, b.cpp, c.cpp, etc)
    // See the SimpleGlob.h header file for details of the flags.
    CSimpleGlob glob(SG_GLOB_NODOT | SG_GLOB_NOCHECK);
    if (SG_SUCCESS != glob.Add(args.FileCount(), args.Files())) {
        _tprintf(_T("Error while globbing files\n"));
        return 1;
    }

    // dump all of the details, the script that was passed on the
    // command line and the expanded file names
    if (ToulBar2::verbose > 0) {
        cout << "cmd line parsing ---> " << glob.FileCount() << " Filename(s) found in command line" << endl;
    }
    vector<string> strfile;
    set<string> strext;

    if (random_desc == NULL) {
        for (int n = 0; n < glob.FileCount() + ((ToulBar2::stdin_format.size() > 0) ? 1 : 0); ++n) {
            string problem = "";
            if (n < glob.FileCount())
                problem = to_string(glob.File(n));

            if (check_file_ext(problem, file_extension_map["wcsp_ext"]) || ToulBar2::stdin_format.compare("wcsp") == 0) {
                if (ToulBar2::verbose >= 0)
                    cout << "loading wcsp file: " << problem << endl;
                strext.insert(".wcsp");
                strfile.push_back(problem);
            }
            if (check_file_ext(problem, file_extension_map["wcspgz_ext"])) {
                if (ToulBar2::verbose >= 0)
                    cout << "loading gzip'd wcsp file: " << problem << endl;
                strext.insert(".wcsp.gz");
                strfile.push_back(problem);
                ToulBar2::gz = true;
            }
            if (check_file_ext(problem, file_extension_map["wcspxz_ext"])) {
                if (ToulBar2::verbose >= 0)
                    cout << "loading xz compressed wcsp file: " << problem << endl;
                strext.insert(".wcsp.xz");
                strfile.push_back(problem);
                ToulBar2::xz = true;
            }
            // CFN file
            if (check_file_ext(problem, file_extension_map["cfn_ext"]) || ToulBar2::stdin_format.compare("cfn") == 0) {
                if (ToulBar2::verbose >= 0)
                    cout << "loading cfn file: " << problem << endl;
                strext.insert(".cfn");
                strfile.push_back(problem);
                ToulBar2::cfn = true;
            }
            // CFN gzip'd file
            if (check_file_ext(problem, file_extension_map["cfngz_ext"])) {
                if (ToulBar2::verbose >= 0)
                    cout << "loading gzip'd cfn file: " << problem << endl;
                strext.insert(".cfn.gz");
                strfile.push_back(problem);
                ToulBar2::cfn = true;
                ToulBar2::gz = true;
            }
            // CFN xz file
            if (check_file_ext(problem, file_extension_map["cfnxz_ext"])) {
                if (ToulBar2::verbose >= 0)
                    cout << "loading xz compressed cfn file: " << problem << endl;
                strext.insert(".cfn.xz");
                strfile.push_back(problem);
                ToulBar2::cfn = true;
                ToulBar2::xz = true;
            }
            // uai  file
            if (check_file_ext(problem, file_extension_map["uai_ext"]) || ToulBar2::stdin_format.compare("uai") == 0) {
                strfile.push_back(problem);
                strext.insert(".uai");
                if (ToulBar2::verbose >= 0)
                    cout << "loading uai file:  " << problem << endl;
                ToulBar2::uai = 1;
                ToulBar2::bayesian = true;
            }
            // uai  gzip'd file
            if (check_file_ext(problem, file_extension_map["uaigz_ext"])) {
                strfile.push_back(problem);
                strext.insert(".uai.gz");
                if (ToulBar2::verbose >= 0)
                    cout << "loading gzip'd uai file:  " << problem << endl;
                ToulBar2::uai = 1;
                ToulBar2::bayesian = true;
                ToulBar2::gz = true;
            }
            // uai xz compressed file
            if (check_file_ext(problem, file_extension_map["uaixz_ext"])) {
                strfile.push_back(problem);
                strext.insert(".uai.xz");
                if (ToulBar2::verbose >= 0)
                    cout << "loading xz compressed uai file:  " << problem << endl;
                ToulBar2::uai = 1;
                ToulBar2::bayesian = true;
                ToulBar2::xz = true;
            }
            // uai log file
            if (check_file_ext(problem, file_extension_map["uai_log_ext"]) || ToulBar2::stdin_format.compare("LG") == 0) {
                strfile.push_back(problem);
                strext.insert(".LG");
                if (ToulBar2::verbose >= 0)
                    cout << "loading uai log file:  " << problem << endl;
                ToulBar2::uai = 2;
                ToulBar2::bayesian = true;
            }
            // uai log file
            if (check_file_ext(problem, file_extension_map["uaigz_log_ext"])) {
                strfile.push_back(problem);
                strext.insert(".LG.gz");
                if (ToulBar2::verbose >= 0)
                    cout << "loading gzip'd uai log file:  " << problem << endl;
                ToulBar2::uai = 2;
                ToulBar2::bayesian = true;
                ToulBar2::gz = true;
            }
            // uai log file
            if (check_file_ext(problem, file_extension_map["uaixz_log_ext"])) {
                strfile.push_back(problem);
                strext.insert(".LG.xz");
                if (ToulBar2::verbose >= 0)
                    cout << "loading xz compressed uai log file:  " << problem << endl;
                ToulBar2::uai = 2;
                ToulBar2::bayesian = true;
                ToulBar2::xz = true;
            }
            // UAI evidence file
            if (check_file_ext(problem, file_extension_map["evid_ext"])) {
                if (ToulBar2::verbose >= 0)
                    cout << "loading evidence file:  " << problem << endl;
                ToulBar2::evidence_file = string(problem);
            }

            // xml file
            if (check_file_ext(problem, file_extension_map["wcspXML_ext"])) {
                if (ToulBar2::verbose >= 0)
                    cout << "loading xml file:" << problem << endl;

                ToulBar2::xmlflag = true;
                if (!ToulBar2::writeSolution) {
                    ToulBar2::writeSolution = 1;
                    solutionFileName = (char*)"sol";
                }
                strext.insert(".xml");
                strfile.push_back(problem);
            }

            // wcnf or cnf file
            if (check_file_ext(problem, file_extension_map["wcnf_ext"]) || ToulBar2::stdin_format.compare("wcnf") == 0 || ToulBar2::stdin_format.compare("cnf") == 0) {
                if (ToulBar2::verbose >= 0)
                    cout << "loading wcnf file:" << problem << endl;
                ToulBar2::wcnf = true;
                strext.insert(".wcnf");
                strfile.push_back(problem);
            } else if (check_file_ext(problem, file_extension_map["cnf_ext"])) {
                if (ToulBar2::verbose >= 0)
                    cout << "loading cnf file:" << problem << endl;
                ToulBar2::wcnf = true;
                strext.insert(".cnf");
                strfile.push_back(problem);
            }
            if (check_file_ext(problem, file_extension_map["wcnfgz_ext"])) {
                if (ToulBar2::verbose >= 0)
                    cout << "loading gzip'd wcnf file:" << problem << endl;
                ToulBar2::wcnf = true;
                strext.insert(".wcnf.gz");
                strfile.push_back(problem);
                ToulBar2::gz = true;
            } else if (check_file_ext(problem, file_extension_map["cnfgz_ext"])) {
                if (ToulBar2::verbose >= 0)
                    cout << "loading gzip'd cnf file:" << problem << endl;
                ToulBar2::wcnf = true;
                strext.insert(".cnf.gz");
                strfile.push_back(problem);
                ToulBar2::gz = true;
            }
            if (check_file_ext(problem, file_extension_map["wcnfxz_ext"])) {
                if (ToulBar2::verbose >= 0)
                    cout << "loading xz compressed wcnf file:" << problem << endl;
                ToulBar2::wcnf = true;
                strext.insert(".wcnf.xz");
                strfile.push_back(problem);
                ToulBar2::xz = true;
            } else if (check_file_ext(problem, file_extension_map["cnfxz_ext"])) {
                if (ToulBar2::verbose >= 0)
                    cout << "loading xz compressed cnf file:" << problem << endl;
                ToulBar2::wcnf = true;
                strext.insert(".cnf.xz");
                strfile.push_back(problem);
                ToulBar2::xz = true;
            }

            // unconstrained quadratic programming file
            if (check_file_ext(problem, file_extension_map["qpbo_ext"]) || ToulBar2::stdin_format.compare("qpbo") == 0) {
                if (ToulBar2::verbose >= 0)
                    cout << "loading quadratic pseudo-Boolean optimization file:" << problem << endl;
                ToulBar2::qpbo = true;
                strext.insert(".qpbo");
                strfile.push_back(problem);
            }
            if (check_file_ext(problem, file_extension_map["qpbogz_ext"])) {
                if (ToulBar2::verbose >= 0)
                    cout << "loading gzip'd quadratic pseudo-Boolean optimization file:" << problem << endl;
                ToulBar2::qpbo = true;
                strext.insert(".qpbo.gz");
                strfile.push_back(problem);
                ToulBar2::gz = true;
            }
            if (check_file_ext(problem, file_extension_map["qpboxz_ext"])) {
                if (ToulBar2::verbose >= 0)
                    cout << "loading xz compressed quadratic pseudo-Boolean optimization file:" << problem << endl;
                ToulBar2::qpbo = true;
                strext.insert(".qpbo.xz");
                strfile.push_back(problem);
                ToulBar2::xz = true;
            }

            // pseudo-Boolean optimization file
            if (check_file_ext(problem, file_extension_map["opb_ext"]) || ToulBar2::stdin_format.compare("opb") == 0) {
                if (ToulBar2::verbose >= 0)
                    cout << "loading pseudo-Boolean optimization file:" << problem << endl;
                ToulBar2::opb = true;
                strext.insert(".opb");
                strfile.push_back(problem);
            }
            if (check_file_ext(problem, file_extension_map["opbgz_ext"])) {
                if (ToulBar2::verbose >= 0)
                    cout << "loading gzip'd pseudo-Boolean optimization file:" << problem << endl;
                ToulBar2::opb = true;
                strext.insert(".opb.gz");
                strfile.push_back(problem);
                ToulBar2::gz = true;
            }
            if (check_file_ext(problem, file_extension_map["opbxz_ext"])) {
                if (ToulBar2::verbose >= 0)
                    cout << "loading xz compressed pseudo-Boolean optimization file:" << problem << endl;
                ToulBar2::opb = true;
                strext.insert(".opb.xz");
                strfile.push_back(problem);
                ToulBar2::xz = true;
            }

            // upperbound file

            if (check_file_ext(problem, file_extension_map["ub_ext"])) {
                if (ToulBar2::verbose >= 0)
                    cout << "loading upper bound from file: " << problem << endl;
                string ubstring = read_UB(problem.c_str());
                if (ubstring.c_str() != NULL) {
                    if (ToulBar2::externalUB.length() != 0) {
                        cerr << "Error: cannot set upper bound from command line and file simultaneously." << endl;
                        exit(EXIT_FAILURE);
                    } else {
                        ToulBar2::externalUB = ubstring;
                    }
                } else {
                    cerr << "error reading UB in " << problem << endl;
                    exit(-1);
                }
            }

            // bep input file
            if (check_file_ext(problem, file_extension_map["bep_ext"])) {
                if (ToulBar2::verbose >= 0)
                    cout << "loading BEP file: " << problem << endl;
                strext.insert(".bep");
                strfile.push_back(problem);
            }

            //////////////////////Mendelian-error analysis and haplotype reconstruction ////////////////////////////////////
            // map file

            if (check_file_ext(problem, file_extension_map["map_ext"])) {
                ToulBar2::map_file = string(problem);
                ToulBar2::haplotype = new Haplotype;
                if (ToulBar2::verbose >= 0)
                    cout << "loading map file: " << ToulBar2::map_file << endl;

                if (glob.FileCount() < 2) {
                    cerr << "pedigree file is missing (.pre): " << endl;
                    exit(-1);
                }
            }

            // pre file
            if (check_file_ext(problem, file_extension_map["pre_ext"])) {
                strfile.push_back(problem);
                strext.insert(".pre");
                if (ToulBar2::verbose >= 0)
                    cout << "loading pre file: " << problem << endl;
                if (glob.FileCount() < 2)
                    ToulBar2::pedigree = new Pedigree;
            }

            //////////////////////VARIABLE ORDERING ////////////////////////////////////
            // filename containing variable order
            //if (strstr(problem,".order"))
            if (check_file_ext(problem, file_extension_map["order_ext"])) {
                if (ToulBar2::verbose >= 0)
                    cout << "loading variable order in file: " << problem << endl;
                //				if (ToulBar2::varOrder) delete [] ToulBar2::varOrder;
                ToulBar2::varOrder = new char[problem.length() + 1];
                sprintf(ToulBar2::varOrder, "%s", problem.c_str());
            }

            //////////////////////TREE DECOMPOSITION AND VARIABLE ORDERING ////////////////////////////////////
            // filename containing list of clusters in topological ordering
            //if (strstr(problem,".cov"))
            if (check_file_ext(problem, file_extension_map["treedec_ext"])) {
                if (ToulBar2::verbose >= 0)
                    cout << "loading tree decomposition in file: " << problem << endl;
                //              if (ToulBar2::varOrder) delete [] ToulBar2::varOrder;
                ToulBar2::varOrder = new char[problem.length() + 1];
                sprintf(ToulBar2::varOrder, "%s", problem.c_str());
                if (!WCSP::isAlreadyTreeDec(ToulBar2::varOrder)) {
                    cerr << "Input tree decomposition file is not valid! (first cluster must be a root, i.e., parentID=-1)" << endl;
                    exit(EXIT_FAILURE);
                }
                if (ToulBar2::btdMode == 0 && ToulBar2::searchMethod == DFBB)
                    ToulBar2::btdMode = 1;
                ToulBar2::clusterFile = string(ToulBar2::varOrder);
            }

            // filename containing list of clusters without the running intersection property
            //if (strstr(problem,".dec"))
            if (check_file_ext(problem, file_extension_map["clusterdec_ext"])) {
                if (ToulBar2::verbose >= 0)
                    cout << "loading cluster decomposition in file: " << problem << endl;
                ToulBar2::clusterFile = problem;
                ifstream decfile(ToulBar2::clusterFile.c_str());
                if (!decfile) {
                    cerr << "File " << ToulBar2::clusterFile << " not found!" << endl;
                    exit(EXIT_FAILURE);
                }
            }

            // read assignment in file or filename of solution
            if (check_file_ext(problem, file_extension_map["sol_ext"])) {
                if (certificateString && strcmp(certificateString, "") != 0) {
                    cerr << "\n COMMAND LINE ERROR cannot read a solution if a partial assignment is given in the command line using -x= argument " << endl;
                    exit(-1);
                }
                if (ToulBar2::verbose >= 0)
                    cout << "loading solution in file: " << problem << endl;

                certificate = true;
                certificateFilename = new char[256];
                sprintf(certificateFilename, "%s", problem.c_str());
                certificateString = (char*)""; // ensure the search will continue starting from this solution
            }
        }
    }

    // ----------simple opt end ----------------------

    //  command line => check existencise of problem filname argument

    // PB FILENAME

    if (argc <= 1) {
        cerr << "Problem filename is missing as command line argument!" << endl;
        cerr << endl;
        help_msg(argv[0]);
        exit(EXIT_FAILURE);
    }

    //------------------------------tb2 option --------------
    /////////////////////////////////////////
    // test on initial ub cost value;
    /////////////////////////////////////////

    if (ToulBar2::verifyOpt && (!certificate || certificateFilename == NULL)) {
        cerr << "Error: no optimal solution file given. Cannot verify the optimal solution." << endl;
        exit(EXIT_FAILURE);
    }

    //TODO: If --show_options then dump ToulBar2 object here

    ToulBar2::startCpuTime = cpuTime();

    initCosts();
    Cost globalUb = MAX_COST;
    WeightedCSPSolver* solver = WeightedCSPSolver::makeWeightedCSPSolver(MAX_COST);

    bool randomproblem = false;
    bool forceSubModular = false;
    string randomglobal;

    int n = 10;
    int m = 2;
    vector<int> p;
#ifndef MENDELSOFT
    if (random_desc != NULL) {
        int pn[10];
        int narities = 0;
        if (strstr(random_desc, "binsub")) {
            forceSubModular = true;
            randomproblem = true;
            sscanf(random_desc, "binsub-%d-%d-%d-%d-%d-%d", &n, &m, &pn[0], &pn[1], &pn[2], &ToulBar2::seed);
            narities = 2;
        } else if (strstr(random_desc, "bin")) {
            randomproblem = true;
            sscanf(random_desc, "bin-%d-%d-%d-%d-%d", &n, &m, &pn[0], &pn[1], &ToulBar2::seed);
            narities = 2;
        } else if (strstr(random_desc, "tern")) {
            randomproblem = true;
            sscanf(random_desc, "tern-%d-%d-%d-%d-%d-%d", &n, &m, &pn[0], &pn[1], &pn[2], &ToulBar2::seed);
            narities = 3;
        } else // if (strstr(random_desc,"salldiff"))
        {
            randomproblem = true;
            char* pch = strtok(random_desc, "-");
            randomglobal = string(pch);
            pch = strtok(NULL, "-");
            n = atoi(pch);
            pch = strtok(NULL, "-");
            m = atoi(pch);

            while (pch != NULL) {
                pch = strtok(NULL, "-");
                if (pch != NULL) {
                    pn[narities] = atoi(pch);
                    narities++;
                }
            }
            narities--;
            ToulBar2::seed = pn[narities];
        }
        if (pn[0] > 100) {
            cout << pn[0] << " tightness is a percentage" << endl;
            pn[0] = 100;
        }
        for (int i = 0; i < narities; i++)
            p.push_back(pn[i]);
        if (forceSubModular)
            p.push_back(pn[narities]);
        if (narities == 0) {
            cerr << "Random problem generator with no cost functions! (no arity)" << endl;
            exit(-1);
        }
    }
    if (strext.count(".bep") || (strfile.size() > 0 && strstr((char*)strfile.back().c_str(), "bEpInstance")))
        ToulBar2::bep = new BEP;
#endif

    if (ToulBar2::seed < 0) { // initialize seed using current time
        ToulBar2::seed = abs((int)time(NULL) * getpid() * ToulBar2::seed);
        if (ToulBar2::verbose >= 0)
            cout << "Initial random seed is " << ToulBar2::seed << endl;
    }
    mysrand(ToulBar2::seed);
    if (ToulBar2::incop_cmd.size() > 0 && ToulBar2::seed != 1 && ToulBar2::incop_cmd.find("0 1 ") == 0) {
        string sseed = to_string(ToulBar2::seed);
        ToulBar2::incop_cmd.replace(2, 1, sseed);
    }

    tb2checkOptions();
    try {
        if (randomproblem)
            solver->read_random(n, m, p, ToulBar2::seed, forceSubModular, randomglobal);
        else {
            if (strfile.size() == 0) {
                cerr << "No problem file given as input!" << endl;
                exit(EXIT_FAILURE);
            } else if (strfile.size() == 1) {
                globalUb = solver->read_wcsp((char*)strfile.back().c_str());
                if (globalUb <= MIN_COST) {
                    THROWCONTRADICTION;
                }
            } else {
                if (strext.size() > 1) {
                    cerr << "Sorry, multiple problem files must have the same file extension!" << endl;
                    exit(EXIT_FAILURE);
                }
                if (strext.begin()->find(".wcsp") == string::npos && strext.begin()->find(".cfn") == string::npos && strext.begin()->find(".xml") == string::npos) {
                    cerr << "Sorry, multiple problem files must have a file extension which contains either '.wcsp' or '.cfn' or '.xml'!" << endl;
                    exit(EXIT_FAILURE);
                }
                for (auto f : strfile) {
                    globalUb = solver->read_wcsp((char*)f.c_str());
                    if (globalUb <= MIN_COST) {
                        THROWCONTRADICTION;
                    }
                }
            }
        }

        //TODO: If --show_options then dump ToulBar2 object here

        if (certificate) {
            if (certificateFilename != NULL)
                solver->read_solution(certificateFilename, updateValueHeuristic);
            else
                solver->parse_solution(certificateString);
        }

#ifdef OPENMPI
        if (env0.myrank == 0) {
#endif
            if (ToulBar2::writeSolution) {
                ToulBar2::solutionFile = fopen(solutionFileName, "w");
                if (!ToulBar2::solutionFile) {
                    cerr << "Could not open file " << solutionFileName << endl;
                    exit(EXIT_FAILURE);
                }
            }
            if (ToulBar2::uaieval) {
                char* tmpPath = new char[strlen(argv[0]) + 1];
                strcpy(tmpPath, argv[0]);
                if (strcmp(tmpPath, "toulbar2") == 0)
                    strcpy(tmpPath, ".");
                char* tmpFile = new char[strlen(strfile.back().c_str()) + 1];
                strcpy(tmpFile, strfile.back().c_str());
                string filename(tmpPath);
                filename += "/";
                filename += basename(tmpFile);
                size_t wcsppos = string::npos;
                if (ToulBar2::uaieval && (wcsppos = filename.rfind(".wcsp")) != string::npos)
                    filename.replace(wcsppos, 5, ".uai");
                filename += ".";
                if (ToulBar2::isZ)
                    filename += "PR";
                else
                    filename += "MPE";
                ToulBar2::solution_uai_filename = filename;
                ToulBar2::solution_uai_file = fopen(ToulBar2::solution_uai_filename.c_str(), "w");
                if (!ToulBar2::solution_uai_file) {
                    cerr << "Could not open file " << ToulBar2::solution_uai_filename << endl;
                    exit(EXIT_FAILURE);
                }
                delete[] tmpPath;
                delete[] tmpFile;
            }
#ifdef OPENMPI
        }
#endif

        if (ToulBar2::problemsaved_filename.empty())
            ToulBar2::problemsaved_filename = ((static_cast<ProblemFormat>((ToulBar2::dumpWCSP >> 1)+(ToulBar2::dumpWCSP & 1)) == CFN_FORMAT) ? "problem.cfn" : "problem.wcsp");

        if (ToulBar2::dumpWCSP % 2) {
            string problemname = ToulBar2::problemsaved_filename;
            if (ToulBar2::uaieval) {
                problemname = ToulBar2::solution_uai_filename;
                problemname.replace(problemname.rfind(".uai.MPE"), 8, (static_cast<ProblemFormat>((ToulBar2::dumpWCSP >> 1)+(ToulBar2::dumpWCSP & 1)) == CFN_FORMAT) ? ".cfn" : ".wcsp");
            }
            solver->dump_wcsp(problemname.c_str(), true, static_cast<ProblemFormat>((ToulBar2::dumpWCSP >> 1)+(ToulBar2::dumpWCSP & 1)));
        } else if (!certificate || certificateString != NULL || ToulBar2::btdMode >= 2) {
#ifndef __WIN32__
            signal(SIGINT, timeOut);
            signal(SIGTERM, timeOut);
            if (timeout > 0)
                timer(timeout);
#endif
            solver->solve();
        }
    } catch (const Contradiction&) {
        if (ToulBar2::verbose >= 0)
            cout << "No solution found by initial propagation!" << endl;
        if (ToulBar2::isZ) {
            if (ToulBar2::uaieval) {
                rewind(ToulBar2::solution_uai_file);
                fprintf(ToulBar2::solution_uai_file, "PR\n");
                fprintf(ToulBar2::solution_uai_file, PrintFormatProb, -numeric_limits<TProb>::infinity());
                fprintf(ToulBar2::solution_uai_file, "\n");
            }
            cout << "Log(Z)= ";
            cout << -numeric_limits<TProb>::infinity() << endl;
        }
        if (ToulBar2::maxsateval) {
            cout << "o " << solver->getWCSP()->getUb() << endl;
            cout << "s UNSATISFIABLE" << endl;
        }
    }
    if (ToulBar2::verbose >= 0)
        cout << "end." << endl;

    if (ToulBar2::solutionFile != NULL) {
        if (ftruncate(fileno(ToulBar2::solutionFile), ftell(ToulBar2::solutionFile)))
            exit(EXIT_FAILURE);
        fclose(ToulBar2::solutionFile);
    }
    if (ToulBar2::solution_uai_file != NULL) {
        if (ftruncate(fileno(ToulBar2::solution_uai_file), ftell(ToulBar2::solution_uai_file)))
            exit(EXIT_FAILURE);
        fclose(ToulBar2::solution_uai_file);
    }

    // for the competition it was necessary to write a file with the optimal sol
    /*char line[1024];
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

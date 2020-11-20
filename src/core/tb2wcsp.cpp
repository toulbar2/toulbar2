/*
 * ****** Weighted constraint satisfaction problem modeling and local reasoning ********
 *
 * Contains also ToulBar2 options expressed by global variable definitions
 */

#include "tb2wcsp.hpp"
#include "tb2enumvar.hpp"
#include "tb2intervar.hpp"
#include "tb2binconstr.hpp"
#include "tb2ternaryconstr.hpp"
#include "tb2naryconstr.hpp"
#include "tb2arithmetic.hpp"
#include "applis/tb2pedigree.hpp"
#include "applis/tb2haplotype.hpp"
#include "tb2vac.hpp"
#include "search/tb2clusters.hpp"

#include "tb2globaldecomposable.hpp"
#include "globals/tb2globalconstr.hpp"
#ifdef ILOGCPLEX
#include "global/tb2lpsconstr.hpp"
#endif
#include "globals/tb2flowbasedconstr.hpp"
#include "globals/tb2alldiffconstr.hpp"
#include "globals/tb2globalcardinalityconstr.hpp"
#include "globals/tb2sameconstr.hpp"
#include "globals/tb2regularflowconstr.hpp"
#include "globals/tb2amongconstr.hpp"
#include "globals/tb2regulardpconstr.hpp"
#include "globals/tb2grammarconstr.hpp"
#include "globals/tb2treeconstr.hpp"
#include "globals/tb2maxconstr.hpp"
#include "tb2clause.hpp"
#include "tb2clqcover.hpp"
#include "tb2knapsack.hpp"

/*
 * Global variables with their default value
 *
 */

int Store::depth = 0;
StoreStack<BTList<Value>, DLink<Value>*> Store::storeDomain(STORE_SIZE);
StoreStack<BTList<ConstraintLink>, DLink<ConstraintLink>*> Store::storeConstraint(STORE_SIZE);
StoreStack<BTList<Variable*>, DLink<Variable*>*> Store::storeVariable(STORE_SIZE);
StoreStack<BTList<Separator*>, DLink<Separator*>*> Store::storeSeparator(STORE_SIZE);

int WCSP::wcspCounter = 0;

int ToulBar2::verbose;

bool ToulBar2::FullEAC;
bool ToulBar2::VACthreshold;
int ToulBar2::nbTimesIsVAC;
int ToulBar2::nbTimesIsVACitThresholdMoreThanOne;
bool ToulBar2::RASPS;
int ToulBar2::useRASPS;
bool ToulBar2::RASPSreset;
int ToulBar2::RASPSnbStrictACVariables;
Cost ToulBar2::RASPSlastitThreshold;
bool ToulBar2::RASPSsaveitThresholds;
vector<pair<Cost, double>> ToulBar2::RASPSitThresholds;
int ToulBar2::RASPSangle;
Long ToulBar2::RASPSnbBacktracks;
int ToulBar2::debug;
string ToulBar2::externalUB;
int ToulBar2::showSolutions;
int ToulBar2::writeSolution;
FILE* ToulBar2::solutionFile;
long ToulBar2::solutionFileRewindPos;
Long ToulBar2::allSolutions;
int ToulBar2::dumpWCSP;
bool ToulBar2::approximateCountingBTD;
int ToulBar2::elimDegree;
int ToulBar2::elimDegree_preprocessing;
int ToulBar2::elimDegree_;
int ToulBar2::elimDegree_preprocessing_;
int ToulBar2::elimSpaceMaxMB;
int ToulBar2::preprocessTernaryRPC;
int ToulBar2::preprocessFunctional;
bool ToulBar2::costfuncSeparate;
int ToulBar2::preprocessNary;
LcLevelType ToulBar2::LcLevel;
bool ToulBar2::QueueComplexity;
bool ToulBar2::binaryBranching;
bool ToulBar2::lastConflict;
int ToulBar2::dichotomicBranching;
unsigned int ToulBar2::dichotomicBranchingSize;
bool ToulBar2::sortDomains;
map<int, ValueCost*> ToulBar2::sortedDomains;
bool ToulBar2::solutionBasedPhaseSaving;
int ToulBar2::lds;
bool ToulBar2::limited;
Long ToulBar2::restart;
Long ToulBar2::backtrackLimit;
bool ToulBar2::generation;
int ToulBar2::minsumDiffusion;
bool ToulBar2::Static_variable_ordering;
int ToulBar2::weightedDegree;
int ToulBar2::weightedTightness;
bool ToulBar2::MSTDAC;
int ToulBar2::DEE;
int ToulBar2::DEE_;
int ToulBar2::nbDecisionVars;
bool ToulBar2::singletonConsistency;
bool ToulBar2::vacValueHeuristic;

externalevent ToulBar2::setvalue;
externalevent ToulBar2::setmin;
externalevent ToulBar2::setmax;
externalevent ToulBar2::removevalue;
externalcostevent ToulBar2::setminobj;
externalsolution ToulBar2::newsolution;
Pedigree* ToulBar2::pedigree;
Haplotype* ToulBar2::haplotype;

bool ToulBar2::cfn;
bool ToulBar2::gz;
bool ToulBar2::xz;
bool ToulBar2::bayesian;
int ToulBar2::uai;
string ToulBar2::evidence_file;
string ToulBar2::stdin_format;
FILE* ToulBar2::solution_uai_file;
string ToulBar2::solution_uai_filename;
string ToulBar2::problemsaved_filename;
TLogProb ToulBar2::markov_log;
bool ToulBar2::xmlflag;
string ToulBar2::map_file;
bool ToulBar2::maxsateval;
bool ToulBar2::uaieval;

int ToulBar2::resolution;
TProb ToulBar2::errorg;
TLogProb ToulBar2::NormFactor;
/// Allele frequencies of founders
/// - 0: 			equal frequencies
/// - 1: 			probs depending on the frequencies found in the problem
/// - otherwise:  read probability distribution from command line
int ToulBar2::foundersprob_class;
vector<TProb> ToulBar2::allelefreqdistrib;
bool ToulBar2::consecutiveAllele;
int ToulBar2::pedigreeCorrectionMode;
int ToulBar2::pedigreePenalty;

int ToulBar2::vac;
Cost ToulBar2::costThreshold;
Cost ToulBar2::costThresholdPre;
string ToulBar2::costThresholdS;
string ToulBar2::costThresholdPreS;
double ToulBar2::trwsAccuracy;
bool ToulBar2::trwsOrder;
unsigned int ToulBar2::trwsNIter;
unsigned int ToulBar2::trwsNIterNoChange;
unsigned int ToulBar2::trwsNIterComputeUb;
double ToulBar2::costMultiplier;
unsigned int ToulBar2::decimalPoint;
string ToulBar2::deltaUbS;
Cost ToulBar2::deltaUb;
Cost ToulBar2::deltaUbAbsolute;
Double ToulBar2::deltaUbRelativeGap;

unsigned int ToulBar2::divNbSol;
unsigned int ToulBar2::divBound;
unsigned int ToulBar2::divWidth;
unsigned int ToulBar2::divMethod;
unsigned int ToulBar2::divRelax;

BEP* ToulBar2::bep;
bool ToulBar2::wcnf;
bool ToulBar2::qpbo;
double ToulBar2::qpboQuadraticCoefMultiplier;
bool ToulBar2::opb;

char* ToulBar2::varOrder;
int ToulBar2::btdMode;
int ToulBar2::btdSubTree;
int ToulBar2::btdRootCluster;

double ToulBar2::startCpuTime;

int ToulBar2::splitClusterMaxSize;
double ToulBar2::boostingBTD;
int ToulBar2::maxSeparatorSize;
int ToulBar2::minProperVarSize;

int ToulBar2::smallSeparatorSize;

bool ToulBar2::isZ;
TLogProb ToulBar2::logZ;
TLogProb ToulBar2::logU;
TLogProb ToulBar2::logepsilon;
int ToulBar2::Berge_Dec = 0; // berge decomposition flag  > 0 if wregular found in the problem

externalfunc ToulBar2::timeOut;
bool ToulBar2::interrupted;

bool ToulBar2::learning;

int ToulBar2::seed;

string ToulBar2::incop_cmd;

string ToulBar2::clusterFile;
ofstream ToulBar2::vnsOutput;

SearchMethod ToulBar2::searchMethod;

VNSSolutionInitMethod ToulBar2::vnsInitSol;
int ToulBar2::vnsLDSmin;
int ToulBar2::vnsLDSmax;
VNSInc ToulBar2::vnsLDSinc;
int ToulBar2::vnsKmin;
int ToulBar2::vnsKmax;
VNSInc ToulBar2::vnsKinc;

int ToulBar2::vnsLDScur;
int ToulBar2::vnsKcur;
VNSVariableHeuristic ToulBar2::vnsNeighborVarHeur;
bool ToulBar2::vnsNeighborChange;
bool ToulBar2::vnsNeighborSizeSync;
bool ToulBar2::vnsParallelLimit;
bool ToulBar2::vnsParallelSync;
string ToulBar2::vnsOptimumS;
Cost ToulBar2::vnsOptimum;
bool ToulBar2::vnsParallel;

Long ToulBar2::hbfs;
Long ToulBar2::hbfsGlobalLimit;
Long ToulBar2::hbfsAlpha; // inverse of minimum node redundancy goal limit
Long ToulBar2::hbfsBeta; // inverse of maximum node redundancy goal limit
ptrdiff_t ToulBar2::hbfsCPLimit; // limit on the number of choice points stored inside open node list
ptrdiff_t ToulBar2::hbfsOpenNodeLimit; // limit on the number of open nodes

bool ToulBar2::verifyOpt;
Cost ToulBar2::verifiedOptimum;

/// \brief initialization of ToulBar2 global variables needed by numberjack/toulbar2
void tb2init()
{
    Store::depth = 0;

    ToulBar2::stdin_format = "";
    ToulBar2::externalUB = "";
    ToulBar2::verbose = 0;

    ToulBar2::FullEAC = false;
    ToulBar2::VACthreshold = false;
    ToulBar2::nbTimesIsVAC = 0;
    ToulBar2::nbTimesIsVACitThresholdMoreThanOne = 0;
    ToulBar2::RASPS = false;
    ToulBar2::useRASPS = 0;
    ToulBar2::RASPSreset = false;
    ToulBar2::RASPSnbStrictACVariables = 0;
    ToulBar2::RASPSlastitThreshold = 1;
    ToulBar2::RASPSsaveitThresholds = false;
    ToulBar2::RASPSangle = 10;
    ToulBar2::RASPSnbBacktracks = 1000;
    ToulBar2::RASPSitThresholds.clear();

    ToulBar2::debug = 0;
    ToulBar2::showSolutions = 0;
    ToulBar2::writeSolution = 0;
    ToulBar2::solutionFile = NULL;
    ToulBar2::solutionFileRewindPos = 0L;
    ToulBar2::allSolutions = 0;
    ToulBar2::dumpWCSP = 0;
    ToulBar2::approximateCountingBTD = false;
    ToulBar2::elimDegree = 3;
    ToulBar2::elimDegree_preprocessing = -1;
    ToulBar2::elimDegree_ = -1;
    ToulBar2::elimDegree_preprocessing_ = -1;
    ToulBar2::elimSpaceMaxMB = 0;
    ToulBar2::preprocessTernaryRPC = 0;
    ToulBar2::preprocessFunctional = 1;
    ToulBar2::costfuncSeparate = true;
    ToulBar2::preprocessNary = 10;
    ToulBar2::LcLevel = LC_EDAC;
    ToulBar2::QueueComplexity = false;
    ToulBar2::binaryBranching = true;
    ToulBar2::lastConflict = true;
    ToulBar2::dichotomicBranching = 1;
    ToulBar2::dichotomicBranchingSize = 10;
    ToulBar2::sortDomains = false;
    ToulBar2::solutionBasedPhaseSaving = true;
    ToulBar2::lds = 0;
    ToulBar2::limited = false;
    ToulBar2::restart = -1;
    ToulBar2::backtrackLimit = LONGLONG_MAX;
    ToulBar2::generation = false;
    ToulBar2::minsumDiffusion = 0;
    ToulBar2::Static_variable_ordering = false;
    ToulBar2::weightedDegree = 1000000;
    ToulBar2::weightedTightness = 0;
    ToulBar2::MSTDAC = false;
    ToulBar2::DEE = 1;
    ToulBar2::DEE_ = 0;
    ToulBar2::nbDecisionVars = 0;
    ToulBar2::singletonConsistency = false;
    ToulBar2::vacValueHeuristic = true;

    ToulBar2::setvalue = NULL;
    ToulBar2::setmin = NULL;
    ToulBar2::setmax = NULL;
    ToulBar2::removevalue = NULL;
    ToulBar2::setminobj = NULL;
    ToulBar2::newsolution = NULL;
    ToulBar2::pedigree = NULL;
    ToulBar2::haplotype = NULL;

    ToulBar2::cfn = false;
    ToulBar2::gz = false;
    ToulBar2::xz = false;
    ToulBar2::bayesian = false;
    ToulBar2::uai = 0;
    ToulBar2::solution_uai_file = NULL;
    ToulBar2::solution_uai_filename = "sol";
    ToulBar2::problemsaved_filename = "";
    ToulBar2::markov_log = 0;
    ToulBar2::xmlflag = false;
    ToulBar2::maxsateval = false;
    ToulBar2::uaieval = false;

    ToulBar2::resolution = 7;
    ToulBar2::errorg = 0.05;
    ToulBar2::NormFactor = 1;
    ToulBar2::foundersprob_class = 0;
    ToulBar2::consecutiveAllele = false;
    ToulBar2::pedigreeCorrectionMode = 0;
    ToulBar2::pedigreePenalty = 0;
    ToulBar2::allelefreqdistrib.clear();

    ToulBar2::vac = 0;
    ToulBar2::costThresholdS = "";
    ToulBar2::costThresholdPreS = "";
    ToulBar2::costThreshold = UNIT_COST;
    ToulBar2::costThresholdPre = UNIT_COST;
    ToulBar2::trwsAccuracy = -1; // 0.001;
    ToulBar2::trwsOrder = false;
    ToulBar2::trwsNIter = 1000;
    ToulBar2::trwsNIterNoChange = 5;
    ToulBar2::trwsNIterComputeUb = 100;
    ToulBar2::costMultiplier = UNIT_COST;
    ToulBar2::decimalPoint = 0;
    ToulBar2::deltaUbS = "0";
    ToulBar2::deltaUb = MIN_COST;
    ToulBar2::deltaUbAbsolute = MIN_COST;
    ToulBar2::deltaUbRelativeGap = 0.;

    ToulBar2::divNbSol = 0;
    ToulBar2::divBound = 0;
    ToulBar2::divWidth = 0;
    ToulBar2::divMethod = 0;
    ToulBar2::divRelax = 0;

    ToulBar2::bep = NULL;
    ToulBar2::wcnf = false;
    ToulBar2::qpbo = false;
    ToulBar2::qpboQuadraticCoefMultiplier = 2.;
    ToulBar2::opb = false;

    ToulBar2::varOrder = NULL;
    ToulBar2::btdMode = 0;
    ToulBar2::btdSubTree = -1;
    ToulBar2::btdRootCluster = -1;

    ToulBar2::startCpuTime = 0;

    ToulBar2::splitClusterMaxSize = 0;
    ToulBar2::boostingBTD = 0.;
    ToulBar2::maxSeparatorSize = -1;
    ToulBar2::minProperVarSize = 0;

    ToulBar2::smallSeparatorSize = 4;

    ToulBar2::isZ = false;
    ToulBar2::logZ = -numeric_limits<TLogProb>::infinity();
    ToulBar2::logU = -numeric_limits<TLogProb>::infinity();
    ToulBar2::logepsilon = -numeric_limits<TLogProb>::infinity();
    ToulBar2::Berge_Dec = 0;

    ToulBar2::timeOut = NULL;
    ToulBar2::interrupted = false;

    ToulBar2::learning = false;

    ToulBar2::seed = 1;

    ToulBar2::incop_cmd = "";

    ToulBar2::searchMethod = DFBB;

    ToulBar2::clusterFile = "";
    ToulBar2::vnsOutput.setstate(std::ios::failbit);

    ToulBar2::vnsInitSol = LS_INIT_DFBB;
    ToulBar2::vnsLDSmin = 1;
    ToulBar2::vnsLDSmax = -1;
    ToulBar2::vnsLDSinc = VNS_MULT2;
    ToulBar2::vnsKmin = 4;
    ToulBar2::vnsKmax = 0;
    ToulBar2::vnsKinc = VNS_ADD1JUMP;

    ToulBar2::vnsLDScur = -1;
    ToulBar2::vnsKcur = 0;
    ToulBar2::vnsNeighborVarHeur = RANDOMVAR;
    ToulBar2::vnsNeighborChange = false;
    ToulBar2::vnsNeighborSizeSync = false;
    ToulBar2::vnsParallelLimit = false;
    ToulBar2::vnsParallelSync = false;
    ToulBar2::vnsOptimumS = "";
    ToulBar2::vnsOptimum = MIN_COST;
    ToulBar2::vnsParallel = false;

    ToulBar2::hbfs = 1;
    ToulBar2::hbfsGlobalLimit = 10000;
    ToulBar2::hbfsAlpha = 20LL; // i.e., alpha = 1/20 = 0.05
    ToulBar2::hbfsBeta = 10LL; // i.e., beta = 1/10 = 0.1
    ToulBar2::hbfsCPLimit = CHOICE_POINT_LIMIT;
    ToulBar2::hbfsOpenNodeLimit = OPEN_NODE_LIMIT;

    ToulBar2::verifyOpt = false;
    ToulBar2::verifiedOptimum = MAX_COST;
}

/// \brief checks compatibility between selected options of ToulBar2 needed by numberjack/toulbar2
void tb2checkOptions()
{
    if (ToulBar2::divBound >= 1 && ToulBar2::divNbSol == 0) {
        cerr << "Error: ask for zero diverse solutions!" << endl;
        exit(1);
    }
    if (ToulBar2::divNbSol >= 1 && ToulBar2::allSolutions > 0) {
        if (ToulBar2::verbose >= 0 && ToulBar2::allSolutions > ToulBar2::divNbSol)
            cout << "Warning! A limit of " << ToulBar2::divNbSol << " diverse solutions has been applied! (try a smaller value with option -a)" << endl;
        ToulBar2::divNbSol = min((Long)ToulBar2::divNbSol, ToulBar2::allSolutions);
        ToulBar2::allSolutions = 0;
    }
    if (ToulBar2::costMultiplier != UNIT_COST && (ToulBar2::uai || ToulBar2::qpbo || ToulBar2::opb)) {
        cerr << "Error: cost multiplier cannot be used with UAI, PBO, and QPBO formats. Use option -precision instead." << endl;
        exit(1);
    }
    if (ToulBar2::costMultiplier != UNIT_COST && (ToulBar2::haplotype || ToulBar2::pedigree || ToulBar2::bep || ToulBar2::xmlflag)) {
        cerr << "Error: cost multiplier not implemented for this file format." << endl;
        exit(1);
    }
    if (ToulBar2::searchMethod != DFBB && ToulBar2::btdMode >= 1) {
        cerr << "Error: BTD-like search methods are compatible with VNS. Deactivate either '-B' or '-vns'" << endl;
        exit(1);
    }
    if (ToulBar2::searchMethod != DFBB && ToulBar2::restart < 1) {
        ToulBar2::restart = 1; // Force random variable selection during (LDS) search within variable neighborhood search methods
    }
    if ((ToulBar2::allSolutions || ToulBar2::isZ) && ToulBar2::searchMethod != DFBB) {
        cerr << "Error: cannot find all solutions or compute a partition function with VNS. Deactivate either option." << endl;
        exit(1);
    }
    if (ToulBar2::divNbSol > 1 && ToulBar2::searchMethod != DFBB) {
        cerr << "Error: cannot find diverse solutions with VNS. Deactivate either option." << endl;
        exit(1);
    }
    if (ToulBar2::approximateCountingBTD && ToulBar2::searchMethod != DFBB) {
        cerr << "Error: cannot compute an approximate solution count with VNS. Deactivate '-vns' for counting." << endl;
        exit(1);
    }
    if (ToulBar2::searchMethod == RPDGVNS && !ToulBar2::vnsParallelSync && ToulBar2::vnsKinc == VNS_LUBY) {
        cerr << "Error: Luby operator not implemented for neighborhood growth strategy in asynchronous parallel VNS-like methods, use Add1 instead." << endl;
        exit(1);
    }
    if (ToulBar2::searchMethod == RPDGVNS && !ToulBar2::vnsParallelSync && ToulBar2::vnsLDSinc == VNS_LUBY) {
        cerr << "Error: Luby operator not implemented for  discrepancy growth strategy in asynchronous parallel VNS-like methods, use Add1 instead." << endl;
        exit(1);
    }
    if (ToulBar2::approximateCountingBTD && ToulBar2::btdMode != 1) {
        cerr << "Error: BTD search mode required for approximate solution counting (use '-B=1')." << endl;
        exit(1);
    }
    if (ToulBar2::allSolutions && ToulBar2::btdMode > 1) {
        cerr << "Error: RDS-like method cannot currently enumerate solutions. Use DFS/HBFS search or BTD (feasibility only)." << endl;
        exit(1);
    }
    if (ToulBar2::divNbSol > 1 && ToulBar2::btdMode >= 1) {
        cerr << "Error: BTD-like methods cannot currently find diverse solutions. Use DFS/HBFS search." << endl;
        exit(1);
    }
    if (ToulBar2::allSolutions && ToulBar2::btdMode == 1 && ToulBar2::elimDegree > 0) {
        //    if (!ToulBar2::uai || ToulBar2::debug) cout << "Warning! Cannot count all solutions with variable elimination during search (except with degree 0 for #BTD)" << endl;
        ToulBar2::elimDegree = 0;
    }
    if (ToulBar2::allSolutions && ToulBar2::btdMode != 1 && ToulBar2::elimDegree >= 0) {
        //    if (!ToulBar2::uai || ToulBar2::debug) cout << "Warning! Cannot count all solutions with variable elimination during search (except with degree 0 for #BTD)" << endl;
        ToulBar2::elimDegree = -1;
    }
    if (ToulBar2::allSolutions && ToulBar2::elimDegree_preprocessing >= 0) {
        //    if (!ToulBar2::uai || ToulBar2::debug) cout << "Warning! Cannot count all solutions with generic variable elimination" << endl;
        ToulBar2::elimDegree_preprocessing = -1;
    }
    if (ToulBar2::allSolutions || ToulBar2::isZ) {
        ToulBar2::DEE = 0;
        ToulBar2::FullEAC = false;
    }
    if (ToulBar2::lds && ToulBar2::btdMode >= 1) {
        cerr << "Error: Limited Discrepancy Search not compatible with BTD-like search methods." << endl;
        exit(1);
    }
    if (ToulBar2::lds && ToulBar2::hbfs) {
        // cout << "Warning! Hybrid best-first search not compatible with Limited Discrepancy Search." << endl;
        ToulBar2::hbfs = 0;
    }
    if (ToulBar2::lds && ToulBar2::solutionBasedPhaseSaving) {
        // cout << "Warning! Solution based phase saving is not recommended with Limited Discrepancy Search." << endl;
        ToulBar2::solutionBasedPhaseSaving = false;
    }
    if (ToulBar2::hbfs && ToulBar2::btdMode >= 2) {
        cout << "Warning! Hybrid best-first search not compatible with RDS-like search methods." << endl;
        ToulBar2::hbfs = 0;
    }
    if (ToulBar2::restart >= 0 && ToulBar2::btdMode >= 1) {
        cerr << "Error: Randomized search with restart not compatible with BTD-like search methods." << endl;
        exit(1);
    }
    if (!ToulBar2::binaryBranching && ToulBar2::btdMode >= 1) {
        cout << "Warning! N-ary branching not implemented with BTD-like search methods (remove -b: or -B option)." << endl;
        exit(1);
    }
    if (ToulBar2::btdSubTree >= 0 && ToulBar2::btdMode <= 1) {
        cerr << "Error: cannot restrict solving to a problem rooted at a subtree, use RDS (-B=2)." << endl;
        exit(1);
    }
    if (abs(ToulBar2::vac) > 1 && ToulBar2::btdMode >= 1) { /// \warning VAC supports can break EAC supports (e.g. SPOT5 404.wcsp)
        cerr << "Error: VAC during search not implemented with BTD-like search methods (use -A only or unset -B)." << endl;
        exit(1);
    }
    if (ToulBar2::FullEAC && ToulBar2::btdMode >= 1) {
        cerr << "Error: VAC-based variable ordering heuristic not implemented with BTD-like search methods (remove -vacint option)." << endl;
        exit(1);
    }
    if (ToulBar2::FullEAC && ToulBar2::LcLevel != LC_EDAC) { /// \warning VAC-integral assumes EAC supports
        cerr << "Error: VAC-based variable ordering heuristic requires EDAC local consistency (select EDAC using -k option)." << endl;
        exit(1);
    }
    if (ToulBar2::useRASPS && ToulBar2::btdMode >= 1) {
        cerr << "Error: VAC-based upper bound probing heuristic not implemented with BTD-like search methods (remove -rasps option)." << endl;
        exit(1);
    }
    if (ToulBar2::useRASPS && !ToulBar2::vac) {
        cerr << "Error: VAC-based upper bound probing heuristic requires VAC at least in preprocessing (add -A option)." << endl;
        exit(1);
    }
    if (ToulBar2::VACthreshold && !ToulBar2::vac) {
        cerr << "Error: VAC threshold heuristic requires VAC during search (add -A option)." << endl;
        exit(1);
    }
    if (ToulBar2::vac && (ToulBar2::LcLevel == LC_NC || ToulBar2::LcLevel == LC_DAC)) { /// \warning VAC assumes AC supports
        cerr << "Error: VAC requires at least AC local consistency (select AC, FDAC, or EDAC using -k option)." << endl;
        exit(1);
    }
    if (ToulBar2::vac && ToulBar2::FullEAC && !ToulBar2::vacValueHeuristic) { /// \warning VAC must update EAC supports in order to make new FullEAC supports based on VAC-integrality
        ToulBar2::vacValueHeuristic = true;
    }
    if (ToulBar2::preprocessFunctional > 0 && (ToulBar2::LcLevel == LC_NC || ToulBar2::LcLevel == LC_DAC)) {
        cerr << "Error: functional elimination requires at least AC local consistency (select AC, FDAC, or EDAC using -k option)." << endl;
        exit(1);
    }
    if (ToulBar2::learning && ToulBar2::elimDegree >= 0) {
        cout << "Warning! Cannot perform variable elimination during search with pseudo-boolean learning." << endl;
        ToulBar2::elimDegree = -1;
    }
    if (ToulBar2::incop_cmd.size() > 0 && (ToulBar2::allSolutions || ToulBar2::isZ)) {
        cout << "Error: Cannot use INCOP local search for (weighted) counting (remove -i option)." << endl;
        exit(1);
    }
    if (!ToulBar2::binaryBranching && ToulBar2::hbfs) {
        cout << "Error: hybrid best-first search restricted to binary branching (remove -b: or add -hbfs: options)." << endl;
        exit(1);
    }
    if (ToulBar2::dichotomicBranching >= 2 && ToulBar2::hbfs) {
        cout << "Error: general dichotomic branching not implemented with hybrid best-first search (use simple dichotomic branching or add -hbfs: parameter)." << endl;
        exit(1);
    }
    if (ToulBar2::verifyOpt && (ToulBar2::elimDegree >= 0 || ToulBar2::elimDegree_preprocessing >= 0)) {
        cout << "Warning! Cannot perform variable elimination while verifying that the optimal solution is preserved." << endl;
        ToulBar2::elimDegree = -1;
        ToulBar2::elimDegree_preprocessing = -1;
    }
    if (ToulBar2::verifyOpt && ToulBar2::preprocessFunctional > 0) {
        cout << "Warning! Cannot perform functional elimination while verifying that the optimal solution is preserved." << endl;
        ToulBar2::preprocessFunctional = 0;
    }
    if (ToulBar2::verifyOpt && ToulBar2::DEE >= 1) {
        cout << "Warning! Cannot perform dead-end elimination while verifying that the optimal solution is preserved." << endl;
        ToulBar2::DEE = 0;
    }
}

/*
 * WCSP constructors
 *
 */

/// \note isDelayedNaryCtr should be false if toulbar2 is used within numberjack
WCSP::WCSP(Cost upperBound, void* _solver_)
    : solver(_solver_)
    , lb(MIN_COST)
    , ub(upperBound)
    , negCost(MIN_COST)
    , solutionCost(MAX_COST)
    , NCBucketSize(cost2log2gub(upperBound) + 1)
    , NCBuckets(NCBucketSize, VariableList(&Store::storeVariable))
    , PendingSeparator(&Store::storeSeparator)
    , objectiveChanged(false)
    , nbNodes(0)
    , nbDEE(0)
    , lastConflictConstr(NULL)
    , maxdomainsize(0)
    ,
#ifdef NUMBERJACK
    isDelayedNaryCtr(false)
    ,
#else
    isDelayedNaryCtr(true)
    ,
#endif
    isPartOfOptimalSolution(0)
    , elimOrder(0)
    , elimBinOrder(0)
    , elimTernOrder(0)
    , maxDegree(-1)
    , elimSpace(0)
{
    instance = wcspCounter++;
    if (ToulBar2::vac)
        vac = new VACExtension(this);
    else
        vac = NULL;

    td = NULL;
}

WCSP::~WCSP()
{
    if (vars.size())
        for (unsigned int i = 0; i < vars.size(); i++)
            delete vars[i];
    if (constrs.size())
        for (unsigned int i = 0; i < constrs.size() - 1; i++)
            delete constrs[i]; // Warning! The last constraint may be badly allocated due to an exception occuring in its constructor (because of propagate) // If there is no constraint then (constrs.size()-1) overflow!
    if (elimBinConstrs.size())
        for (unsigned int i = 0; i < elimBinConstrs.size(); i++)
            delete elimBinConstrs[i];
    if (elimTernConstrs.size())
        for (unsigned int i = 0; i < elimTernConstrs.size(); i++)
            delete elimTernConstrs[i];
}

WeightedCSP* WeightedCSP::makeWeightedCSP(Cost upperBound, void* solver)
{
    WeightedCSP* W = new WCSP(upperBound, solver);
    return W;
}

/// \brief create an enumerated variable with its domain bounds
int WCSP::makeEnumeratedVariable(string n, Value iinf, Value isup)
{
    EnumeratedVariable* x;
    if (!ToulBar2::vac) {
        x = new EnumeratedVariable(this, n, iinf, isup);
    } else {
        x = new VACVariable(this, n, iinf, isup);
    }
    if (maxdomainsize < isup - iinf + 1)
        maxdomainsize = isup - iinf + 1;
    listofsuccessors.push_back(vector<int>()); // add new variable in the topological order list;
    return x->wcspIndex;
}

/// \brief create an enumerated variable with its domain values
int WCSP::makeEnumeratedVariable(string n, Value* d, int dsize)
{
    EnumeratedVariable* x;
    if (!ToulBar2::vac) {
        x = new EnumeratedVariable(this, n, d, dsize);
    } else {
        x = new VACVariable(this, n, d, dsize);
    }
    if (maxdomainsize < dsize)
        maxdomainsize = dsize;
    listofsuccessors.push_back(vector<int>()); // add new variable in the topological order list;
    return x->wcspIndex;
}

void WCSP::addValueName(int xIndex, const string& name)
{
    Variable* x = getVar(xIndex);
    if (x->enumerated()) {
        ((EnumeratedVariable*)x)->addValueName(name);
    }
}

/// \brief create an interval variable with its domain bounds
int WCSP::makeIntervalVariable(string n, Value iinf, Value isup)
{
    if (ToulBar2::vac) {
        cerr << "VAC not implemented on interval variables!" << endl;
        ToulBar2::vac = 0;
        ToulBar2::minsumDiffusion = 0;
    }
    IntervalVariable* x = new IntervalVariable(this, n, iinf, isup);
    if (maxdomainsize < isup - iinf + 1)
        maxdomainsize = isup - iinf + 1;
    listofsuccessors.push_back(vector<int>()); // add new variable in the topological order list;
    return x->wcspIndex;
}

/// \brief create a binary cost function from a vector of costs
/// \param xIndex index of enumerated variable x as returned by makeEnumeratedVariable
/// \param yIndex index of enumerated variable y
/// \param costs a flat vector of costs (y indexes moving first)
///
/// \note It looks for an existing constraint
/// (even not connected constraints). It also allocates
/// memory for a new constraint. For these two reasons it should
/// ONLY be called before search.
///
/// \warning Vector costs must have the same size as Cartesian product of original domains.
int WCSP::postBinaryConstraint(int xIndex, int yIndex, vector<Cost>& costs)
{
    assert(xIndex != yIndex);
    EnumeratedVariable* x = (EnumeratedVariable*)vars[xIndex];
    EnumeratedVariable* y = (EnumeratedVariable*)vars[yIndex];

    assert(costs.size() == (x->getDomainInitSize() * y->getDomainInitSize()));
    if (ToulBar2::vac) {
        for (unsigned int a = 0; a < x->getDomainInitSize(); a++) {
            for (unsigned int b = 0; b < y->getDomainInitSize(); b++) {
                Cost c = costs[a * y->getDomainInitSize() + b];
                histogram(c);
            }
        }
    }

    BinaryConstraint* ctr = x->getConstr(y);
    if (ctr) {
        ctr->reconnect();
        ctr->addCosts(x, y, costs);
        ctr->propagate();
    } else {
        if (!ToulBar2::vac) {
            ctr
                = new BinaryConstraint(this, (EnumeratedVariable*)vars[xIndex], (EnumeratedVariable*)vars[yIndex], costs);
        } else {
            ctr
                = new VACBinaryConstraint(this, (EnumeratedVariable*)vars[xIndex], (EnumeratedVariable*)vars[yIndex], costs);
        }
    }

    return ctr->wcspIndex;
}

/// \brief create a binary cost function from a vector of floating point values that will be approximated to the ToulBar2:decimalPoint precision
/// \param xIndex index of enumerated variable x as returned by makeEnumeratedVariable
/// \param yIndex index of enumerated variable y
/// \param costs a flat vector of floating point costs (y indexes moving first)
/// \param incremental true if this cost function is automatically removed when backtracking
///
/// \note This must ONLY be called before search.
///
/// \warning The cost Vector must have the same size as Cartesian product of original domains.
int WCSP::postBinaryConstraint(int xIndex, int yIndex, vector<Double>& dcosts, bool incremental)
{
    assert(xIndex != yIndex);
    EnumeratedVariable* x = (EnumeratedVariable*)vars[xIndex];
    EnumeratedVariable* y = (EnumeratedVariable*)vars[yIndex];

    assert(dcosts.size() == (x->getDomainInitSize() * y->getDomainInitSize()));

    long double minCost = std::numeric_limits<long double>::infinity();
    for (long double cost : dcosts) {
        minCost = min(minCost, cost);
    }

    vector<Cost> icosts;
    icosts.resize(dcosts.size());
    for (unsigned int a = 0; a < x->getDomainInitSize(); a++) {
        for (unsigned int b = 0; b < y->getDomainInitSize(); b++) {
            icosts[a * y->getDomainInitSize() + b] = (Cost)(round((dcosts[a * y->getDomainInitSize() + b] - minCost) * pow(10, ToulBar2::decimalPoint)));
        }
    }
    negCost -= (Cost)(round(minCost * pow(10, ToulBar2::decimalPoint)));
    if (incremental) {
        return postIncrementalBinaryConstraint(xIndex, yIndex, icosts);
    } else {
        return postBinaryConstraint(xIndex, yIndex, icosts);
    }
}

/// \brief create a binary cost function from a vector of floating point values that will be approximated to the ToulBar2:decimalPoint precision
/// \param xIndex index of enumerated variable x as returned by makeEnumeratedVariable
/// \param yIndex index of enumerated variable y
/// \param costs a flat vector of floating point costs (y indexes moving first)
/// \param incremental true if this cost function is automatically removed when backtracking
///
/// \note This must ONLY be called before search.
///
/// \warning The cost Vector must have the same size as Cartesian product of original domains.
int WCSP::postTernaryConstraint(int xIndex, int yIndex, int zIndex, vector<Double>& dcosts, bool incremental)
{
    assert(xIndex != yIndex);
    assert(xIndex != zIndex);
    assert(yIndex != zIndex);
    EnumeratedVariable* x = (EnumeratedVariable*)vars[xIndex];
    EnumeratedVariable* y = (EnumeratedVariable*)vars[yIndex];
    EnumeratedVariable* z = (EnumeratedVariable*)vars[zIndex];

    assert(dcosts.size() == (x->getDomainInitSize() * y->getDomainInitSize() * z->getDomainInitSize()));

    long double minCost = std::numeric_limits<long double>::infinity();
    for (long double cost : dcosts) {
        minCost = min(minCost, cost);
    }

    vector<Cost> icosts;
    icosts.resize(dcosts.size());
    for (unsigned int a = 0; a < x->getDomainInitSize(); a++) {
        for (unsigned int b = 0; b < y->getDomainInitSize(); b++) {
            for (unsigned int c = 0; c < z->getDomainInitSize(); c++) {
                icosts[a * y->getDomainInitSize() * z->getDomainInitSize() + b * z->getDomainInitSize() + c] = (Cost)(round((dcosts[a * y->getDomainInitSize() * z->getDomainInitSize() + b * z->getDomainInitSize() + c] - minCost) * pow(10, ToulBar2::decimalPoint)));
            }
        }
    }
    negCost -= (Cost)(round(minCost * pow(10, ToulBar2::decimalPoint)));
    if (incremental) {
        return postIncrementalTernaryConstraint(xIndex, yIndex, zIndex, icosts);
    } else {
        return postTernaryConstraint(xIndex, yIndex, zIndex, icosts);
    }
}

/// \brief create a ternary cost function from a flat vector of costs (z indexes moving first)
int WCSP::postTernaryConstraint(int xIndex, int yIndex, int zIndex, vector<Cost>& costs)
{
    assert(xIndex != yIndex && xIndex != zIndex && yIndex != zIndex);
    EnumeratedVariable* x = (EnumeratedVariable*)vars[xIndex];
    EnumeratedVariable* y = (EnumeratedVariable*)vars[yIndex];
    EnumeratedVariable* z = (EnumeratedVariable*)vars[zIndex];

    if (ToulBar2::vac) {
        for (unsigned int a = 0; a < x->getDomainInitSize(); a++) {
            for (unsigned int b = 0; b < y->getDomainInitSize(); b++) {
                for (unsigned int c = 0; c < z->getDomainInitSize(); c++) {
                    Cost co = costs[a * y->getDomainInitSize() * z->getDomainInitSize() + b * z->getDomainInitSize() + c];
                    histogram(co);
                }
            }
        }
    }

    TernaryConstraint* ctr = x->getConstr(y, z);

    if (!ctr) {
        unsigned int a, b;
        vector<Cost> zerocostsxy;
        vector<Cost> zerocostsxz;
        vector<Cost> zerocostsyz;

        for (a = 0; a < x->getDomainInitSize(); a++) {
            for (b = 0; b < y->getDomainInitSize(); b++) {
                zerocostsxy.push_back(MIN_COST);
            }
        }
        for (a = 0; a < x->getDomainInitSize(); a++) {
            for (b = 0; b < z->getDomainInitSize(); b++) {
                zerocostsxz.push_back(MIN_COST);
            }
        }
        for (a = 0; a < y->getDomainInitSize(); a++) {
            for (b = 0; b < z->getDomainInitSize(); b++) {
                zerocostsyz.push_back(MIN_COST);
            }
        }

        BinaryConstraint* xy = x->getConstr(y);
        BinaryConstraint* xz = x->getConstr(z);
        BinaryConstraint* yz = y->getConstr(z);

        if (!ToulBar2::vac) {
            if (!xy) {
                xy = new BinaryConstraint(this, x, y, zerocostsxy);
                xy->deconnect(true);
            }
            if (!xz) {
                xz = new BinaryConstraint(this, x, z, zerocostsxz);
                xz->deconnect(true);
            }
            if (!yz) {
                yz = new BinaryConstraint(this, y, z, zerocostsyz);
                yz->deconnect(true);
            }
        } else {
            if (!xy) {
                xy = new VACBinaryConstraint(this, x, y, zerocostsxy);
                xy->deconnect(true);
            }
            if (!xz) {
                xz = new VACBinaryConstraint(this, x, z, zerocostsxz);
                xz->deconnect(true);
            }
            if (!yz) {
                yz = new VACBinaryConstraint(this, y, z, zerocostsyz);
                yz->deconnect(true);
            }
        }

        ctr = new TernaryConstraint(this, x, y, z, xy, xz, yz, costs);
    } else {
        ctr->addCosts(x, y, z, costs);
        ctr->propagate();
    }

    return ctr->wcspIndex;
}

/// \brief create a global cost function using a default cost (tuples with a different cost will be enter later using WCSP::postNaryConstraintTuple)
/// \param scopeIndex array of enumerated variable indexes (as returned by makeEnumeratedVariable)
/// \param arity size of scopeIndex
/// \param defval default cost for any tuple
/// \param nbtuples number of tuples to be inserted after creating this global cost function (optional parameter, zero if unknown)
/// \param forcenary if false then it may create a specialized global cost function (e.g., a clause) instead of a basic NaryConstraint object
/// \note should not be used for unary or binary or ternary cost functions
/// \warning do not forget to initially propagate the global cost function using WCSP::postNaryConstraintEnd
int WCSP::postNaryConstraintBegin(int* scopeIndex, int arity, Cost defval, Long nbtuples, bool forcenary)
{
#ifndef NDEBUG
    for (int i = 0; i < arity; i++)
        for (int j = i + 1; j < arity; j++)
            assert(scopeIndex[i] != scopeIndex[j]);
#endif
    EnumeratedVariable** scopeVars = new EnumeratedVariable*[arity];
    bool binary = true;
    for (int i = 0; i < arity; i++) {
        scopeVars[i] = (EnumeratedVariable*)vars[scopeIndex[i]];
        if (scopeVars[i]->getDomainInitSize() != 2)
            binary = false;
    }
    AbstractNaryConstraint* ctr = NULL;
    if (!forcenary && binary && nbtuples == 1 && defval == MIN_COST && arity > NARYPROJECTIONSIZE) {
        ctr = new WeightedClause(this, scopeVars, arity);
    } else {
        ctr = new NaryConstraint(this, scopeVars, arity, defval, nbtuples);
    }
    if (arity > NARYPROJECTIONSIZE) {
        if (isDelayedNaryCtr)
            delayedNaryCtr.push_back(ctr->wcspIndex);
        else {
            BinaryConstraint* bctr;
            TernaryConstraint* tctr = new TernaryConstraint(this);
            elimTernConstrs.push_back(tctr);
            for (int j = 0; j < 3; j++) {
                if (!ToulBar2::vac)
                    bctr = new BinaryConstraint(this);
                else
                    bctr = new VACBinaryConstraint(this);
                elimBinConstrs.push_back(bctr);
            }
        }
    }
    delete[] scopeVars;
    return ctr->wcspIndex;
}

/// \brief set one tuple with a specific cost for global cost function in extension
/// \param ctrindex index of cost function as returned by WCSP::postNaryConstraintBegin
/// \param tuple array of values assigned to variables ordered by following the original scope order
/// \param arity size of the array
/// \param cost new cost for this tuple
/// \warning valid only for global cost function in extension
void WCSP::postNaryConstraintTuple(int ctrindex, Value* tuple, int arity, Cost cost)
{
    static Tuple s;
    if (ToulBar2::vac)
        histogram(cost);
    Constraint* ctr = getCtr(ctrindex);
    //    assert(ctr->extension()); // must be an NaryConstraint or WeightedClause
    assert(arity == ctr->arity());
    s.resize(arity);
    for (int i = 0; i < arity; i++)
        s[i] = ((EnumeratedVariable*)ctr->getVar(i))->toIndex(tuple[i]);
    ctr->setTuple(s, cost);
}

/// \brief set one tuple with a specific cost for global cost function in extension
/// \param ctrindex index of cost function as returned by WCSP::postNaryConstraintBegin
/// \param tuple Tuple encoding of values assigned to variables ordered by following the original scope order
/// \param cost new cost for this tuple
/// \warning valid only for global cost function in extension
/// \warning string encoding of tuples is for advanced users only!
void WCSP::postNaryConstraintTuple(int ctrindex, const Tuple& tuple, Cost cost)
{
    if (ToulBar2::vac)
        histogram(cost);
    Constraint* ctr = getCtr(ctrindex);
    //    assert(ctr->extension()); // must be an NaryConstraint or WeightedClause
    ctr->setTuple(tuple, cost);
}

void WCSP::postNaryConstraintEnd(int ctrindex)
{
    AbstractNaryConstraint* ctr = (AbstractNaryConstraint*)getCtr(ctrindex);
    if (ctr->arity() <= NARYPROJECTIONSIZE)
        ctr->projectNaryBeforeSearch();
    else if (!isDelayedNaryCtr)
        ctr->propagate();
}

// Add a temporary (backtrackable) binary constraint for incremental search (like "on the fly ElimVar")
int WCSP::postIncrementalBinaryConstraint(int xIndex, int yIndex, vector<Cost>& costs)
{
    assert(getTreeDec() == NULL);
    EnumeratedVariable* x = (EnumeratedVariable*)getVar(xIndex);
    EnumeratedVariable* y = (EnumeratedVariable*)getVar(yIndex);
    BinaryConstraint* xy = x->getConstr(y);

    initElimConstr();
    BinaryConstraint* xynew = newBinaryConstr(x, y, NULL, NULL);
    elimBinOrderInc();
    for (auto iterx = x->begin(); iterx != x->end(); ++iterx) {
        for (auto itery = y->begin(); itery != y->end(); ++itery) {
            xynew->setcost(*iterx, *itery, costs[x->toIndex(*iterx) * y->getDomainInitSize() + y->toIndex(*itery)]);
        }
    }
    if (xy) {
        xy->addCosts(xynew);
        assert(xynew->deconnected());
        if (x->unassigned() && y->unassigned())
            xy->reconnect();
    } else {
        xy = xynew;
        xy->reconnect();
    }
    xy->propagate();
    return xy->wcspIndex;
}

// Add a temporary (backtrackable) ternary constraint for incremental search (like "on the fly ElimVar")
int WCSP::postIncrementalTernaryConstraint(int xIndex, int yIndex, int zIndex, vector<Cost>& costs)
{
    assert(getTreeDec() == NULL);
    BinaryConstraint* bctr;
    TernaryConstraint* xyz = new TernaryConstraint(this);
    elimTernConstrs.push_back(xyz);
    for (int j = 0; j < 3; j++) {
        if (!ToulBar2::vac)
            bctr = new BinaryConstraint(this);
        else
            bctr = new VACBinaryConstraint(this);
        elimBinConstrs.push_back(bctr);
    }

    EnumeratedVariable* z = (EnumeratedVariable*)getVar(zIndex);
    EnumeratedVariable* y = (EnumeratedVariable*)getVar(yIndex);
    EnumeratedVariable* x = (EnumeratedVariable*)getVar(xIndex);

    xyz = newTernaryConstr(x, y, z);
    elimTernOrderInc();
    for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
        for (EnumeratedVariable::iterator itery = y->begin(); itery != y->end(); ++itery) {
            for (EnumeratedVariable::iterator iterz = z->begin(); iterz != z->end(); ++iterz) {
                xyz->setcost(*iterx, *itery, *iterz, costs[x->toIndex(*iterx) * y->getDomainInitSize() * z->getDomainInitSize() + y->toIndex(*itery) * z->getDomainInitSize() + z->toIndex(*iterz)]);
            }
        }
    }

    TernaryConstraint* ctr = x->getConstr(y, z);
    if (!ctr) {
        xyz->fillElimConstrBinaries();
        xyz->reconnect();
    } else {
        ctr->addCosts(xyz);
        assert(ctr->connected());
        assert(xyz->deconnected());
        xyz = ctr;
    }
    xyz->propagate();
    return xyz->wcspIndex;
}

// Dual representation, Hamming, the dual variables have 2 (divBound + 1) values
// each value represent a (d,d) or (d,d+1) transition except for the last two
// that means we are beyond divBound.
void WCSP::addDivConstraint(const vector<Value> solution, int sol_j, Cost cost)
{
    if (ToulBar2::verbose >= 1)
        cout << "adding diversity constraint (dual)" << endl;

    // add diversity constraint from solution sol_id
    vector<Cost> vc;
    EnumeratedVariable* ex;

    EnumeratedVariable* c;
    EnumeratedVariable* cp = NULL;
    int cId;
    int cpId;

    bool first_pos = true;

    for (Variable* x : getDivVariables()) {

        ex = (EnumeratedVariable*)x;
        int xId = x->wcspIndex; //index of variable x

        // Add constraint between x and c_j_x
        cId = divVarsId[sol_j][xId]; //index of variable c
        c = (EnumeratedVariable*)getVar(cId);
        vc.clear();
        for (unsigned val_x = 0; val_x < ex->getDomainInitSize(); val_x++) { //val_x = value index ?
            for (unsigned val_c = 0; val_c < c->getDomainInitSize(); val_c++) { // Si les domaines sont réduits, ça ne marche pas!
                // val_c = delta(divBound+1) + qp
                unsigned delta = val_c / (ToulBar2::divBound + 1); // the first divBound+1 values of c (delta = 0) are (d,d) transitions, the rest is (d,d+1)
                vc.push_back(((val_x != (unsigned)solution[xId]) == delta) ? MIN_COST : getUb());
            }
        }
        postIncrementalBinaryConstraint(xId, cId, vc);

        // Add constraint between c_j_x and c_j_x-1
        vc.clear();
        if (first_pos) {
            for (unsigned val_c = 0; val_c < c->getDomainInitSize(); val_c++) {
                vc.push_back(((val_c % (ToulBar2::divBound + 1)) == 0) ? MIN_COST : getUb()); // allow only (0,0) and (0,1)
            }
            postIncrementalUnaryConstraint(cId, vc);
            first_pos = false;
        } else {
            // add binary constraint between cp and c
            for (unsigned val_cp = 0; val_cp < cp->getDomainInitSize(); val_cp++) {
                for (unsigned val_c = 0; val_c < c->getDomainInitSize(); val_c++) {
                    unsigned deltap = val_cp / (ToulBar2::divBound + 1); // 1 iff (d,d+1) on previous
                    unsigned qp = val_cp % (ToulBar2::divBound + 1); // (d,.) on previous, so deltap+qp is (.,d/d+1)
                    unsigned q = val_c % (ToulBar2::divBound + 1); // (d,.) on current
                    vc.push_back((q == min(ToulBar2::divBound, qp + deltap)) ? MIN_COST : getUb()); // why min? (Thomas)
                }
            }
            postIncrementalBinaryConstraint(cpId, cId, vc);
        }
        cp = c;
        cpId = cId;
    }

    //Unary constraint on last var_c to ensure diversity
    if (!getDivVariables().empty()) {
        vc.clear();
        for (unsigned val_c = 0; val_c < c->getDomainInitSize(); val_c++) {
            unsigned q = val_c % (ToulBar2::divBound + 1);
            unsigned delta = val_c / (ToulBar2::divBound + 1);
            vc.push_back(q + delta >= ToulBar2::divBound ? MIN_COST : getUb());
        }
        postIncrementalUnaryConstraint(cId, vc);
    }
}

// Hidden representation: dual + state variables
void WCSP::addHDivConstraint(const vector<Value> solution, int sol_j, Cost cost)
{
    if (ToulBar2::verbose >= 1)
        cout << "adding diversity constraint (hidden)" << endl;

    // add diversity constraint from solution sol_id
    vector<Cost> vc;
    EnumeratedVariable* ex;

    EnumeratedVariable* c = NULL;
    EnumeratedVariable* h = NULL;
    EnumeratedVariable* hp = NULL;

    int cId = -1;
    int hId = -1;
    int hpId = -1;

    for (unsigned divVarPos = 0; divVarPos < divVariables.size(); divVarPos++) {
        ex = (EnumeratedVariable*)divVariables[divVarPos];
        int xId = ex->wcspIndex; //wcsp index of current divVariable

        // Add constraint between ex and c_j_x
        cId = divVarsId[sol_j][xId]; //index of variable c
        c = (EnumeratedVariable*)getVar(cId);
        vc.clear();
        for (unsigned val_ex = 0; val_ex < ex->getDomainInitSize(); val_ex++) { // val_x = value index ?
            for (unsigned val_c = 0; val_c < c->getDomainInitSize(); val_c++) {
                unsigned delta = val_c / (ToulBar2::divBound + 1);
                vc.push_back(((val_ex != (unsigned)solution[xId]) == delta) ? MIN_COST : getUb());
            }
        }
        postIncrementalBinaryConstraint(xId, cId, vc);

        // Add constraint between hp_j_x and c_j_x
        vc.clear();
        if (divVarPos == 0) {
            for (unsigned val_c = 0; val_c < c->getDomainInitSize(); val_c++) {
                vc.push_back(((val_c % (ToulBar2::divBound + 1)) == 0) ? MIN_COST : getUb());
            }
            postIncrementalUnaryConstraint(cId, vc);
        } else {
            // add binary constraint between previous hp_j_x and c_j_x
            for (unsigned val_hp = 0; val_hp < hp->getDomainInitSize(); val_hp++) {
                for (unsigned val_c = 0; val_c < c->getDomainInitSize(); val_c++) {
                    unsigned q = val_c % (ToulBar2::divBound + 1);
                    vc.push_back((val_hp == q) ? MIN_COST : getUb());
                }
            }
            postIncrementalBinaryConstraint(hpId, cId, vc);
        }

        // Add constraint between c_j_x and h_j_x a (except last var, done after loop)
        vc.clear();
        if (divVarPos < divVariables.size() - 1) {
            hId = divHVarsId[sol_j][xId]; //index of variable h
            h = (EnumeratedVariable*)getVar(hId);
            for (unsigned val_c = 0; val_c < c->getDomainInitSize(); val_c++) {
                for (unsigned val_h = 0; val_h < h->getDomainInitSize(); val_h++) {
                    unsigned q = val_c % (ToulBar2::divBound + 1);
                    unsigned delta = val_c / (ToulBar2::divBound + 1);
                    vc.push_back((val_h == min(ToulBar2::divBound, q + delta)) ? MIN_COST : getUb());
                }
            }
            postIncrementalBinaryConstraint(cId, hId, vc);
        }
        hp = h;
        hpId = hId;
    }
    //Unary constraint on last var_c to ensure diversity
    if (!divVariables.empty()) {
        vc.clear();
        for (unsigned val_c = 0; val_c < c->getDomainInitSize(); val_c++) {
            unsigned q = val_c % (ToulBar2::divBound + 1);
            unsigned delta = val_c / (ToulBar2::divBound + 1);
            vc.push_back(q + delta >= ToulBar2::divBound ? MIN_COST : getUb());
        }
        postIncrementalUnaryConstraint(cId, vc);
    }
}

// Ternary representation
void WCSP::addTDivConstraint(const vector<Value> solution, int sol_j, Cost cost)
{
    if (ToulBar2::verbose >= 1)
        cout << "adding diversity constraint (ternary decomposition)" << endl;

    // add diversity constraint from solution sol_id
    vector<Cost> vc;
    EnumeratedVariable* ex;

    EnumeratedVariable* h = NULL;
    EnumeratedVariable* hp = NULL;

    int hId = -1;
    int hpId = -1;

    for (unsigned divVarPos = 0; divVarPos < divVariables.size(); divVarPos++) {
        ex = (EnumeratedVariable*)divVariables[divVarPos];
        int xId = ex->wcspIndex; //wcsp index of current divVariable

        // Add constraint between hp ex and h (except on extremities)
        vc.clear();
        if (divVarPos == 0) { // first position, no hp, only distance 0
            hId = divHVarsId[sol_j][xId]; //index of variable h
            h = (EnumeratedVariable*)getVar(hId);
            for (unsigned val_ex = 0; val_ex < ex->getDomainInitSize(); val_ex++) {
                for (unsigned val_h = 0; val_h < h->getDomainInitSize(); val_h++) {
                    vc.push_back(((val_ex != (unsigned)solution[xId]) == val_h) ? MIN_COST : getUb());
                }
            }
            postIncrementalBinaryConstraint(xId, hId, vc);
        } else if (divVarPos + 1 == divVariables.size()) { // last position, no h
            for (unsigned val_hp = 0; val_hp < hp->getDomainInitSize(); val_hp++) {
                for (unsigned val_ex = 0; val_ex < ex->getDomainInitSize(); val_ex++) {
                    vc.push_back((val_hp + (val_ex != (unsigned)solution[xId]) >= ToulBar2::divBound) ? MIN_COST : getUb());
                }
            }
            postIncrementalBinaryConstraint(hpId, xId, vc);
        } else {
            hId = divHVarsId[sol_j][xId]; //index of variable h
            h = (EnumeratedVariable*)getVar(hId);
            for (unsigned val_hp = 0; val_hp < hp->getDomainInitSize(); val_hp++) {
                for (unsigned val_ex = 0; val_ex < ex->getDomainInitSize(); val_ex++) {
                    for (unsigned val_h = 0; val_h < h->getDomainInitSize(); val_h++) {
                        //cout << "hp " << val_hp << " x " << val_ex << " h " << val_h << "(" << ((val_hp + (val_ex != (unsigned)solution[xId]) == val_h) ? MIN_COST : getUb()) << ")" << endl;
                        vc.push_back((min(ToulBar2::divBound, val_hp + (val_ex != (unsigned)solution[xId])) == val_h) ? MIN_COST : getUb());
                    }
                }
            }
            postIncrementalTernaryConstraint(hpId, xId, hId, vc);
        }
        hp = h;
        hpId = hId;
    }
}

void WCSP::addMDDConstraint(Mdd mdd, int relaxed)
{ //sol_j: to recognize the set of variables to use
    if (ToulBar2::verbose >= 1)
        cout << "adding relaxed mdd constraint" << endl;

    vector<Variable*> varReverse;
    for (unsigned v = 0; v < getDivVariables().size(); v++) {
        varReverse.push_back(getDivVariables()[getDivVariables().size() - v - 1]);
    }
    int nLayers = varReverse.size();
    vector<Cost> vc;
    bool first_pos = true;
    EnumeratedVariable* x;

    EnumeratedVariable* c;
    EnumeratedVariable* cp = NULL;
    int cId;
    int cpId = -1;

    unsigned source;
    unsigned target;
    for (int layer = 0; layer < nLayers; layer++) {
        x = (EnumeratedVariable*)varReverse[layer];
        int xId = x->wcspIndex; //index of variable x

        // Add constraint between x and c_j_x
        cId = divVarsId[relaxed][xId]; //index of variable c
        c = (EnumeratedVariable*)getVar(cId);
        vc.clear();
        for (unsigned val_x = 0; val_x < x->getDomainInitSize(); val_x++) {
            for (unsigned val_c = 0; val_c < c->getDomainInitSize(); val_c++) {
                source = val_c / (ToulBar2::divWidth);
                target = val_c % ToulBar2::divWidth;
                vc.push_back((source < mdd[layer].size() && target < mdd[layer][source].size()) ? mdd[layer][source][target][val_x] : getUb());
            }
        }
        postIncrementalBinaryConstraint(xId, cId, vc);

        vc.clear();
        for (unsigned val_c = 0; val_c < c->getDomainInitSize(); val_c++) {
            source = val_c / (ToulBar2::divWidth);
            target = val_c % ToulBar2::divWidth;
            vc.push_back((source < mdd[layer].size() && target < mdd[layer][source].size()) ? MIN_COST : getUb());
        }
        postIncrementalUnaryConstraint(cId, vc);
        vc.clear();
        if (first_pos) {
            first_pos = false;
        } else {
            // add binary constraint between cp and c
            for (unsigned val_cp = 0; val_cp < cp->getDomainInitSize(); val_cp++) {
                for (unsigned val_c = 0; val_c < c->getDomainInitSize(); val_c++) {
                    target = val_cp % (ToulBar2::divWidth);
                    source = val_c / (ToulBar2::divWidth);
                    vc.push_back((target == source) ? MIN_COST : getUb());
                }
            }
            postIncrementalBinaryConstraint(cpId, cId, vc);
            vc.clear();
        }
        cp = c;
        cpId = cId;
    }
}

void WCSP::addHMDDConstraint(Mdd mdd, int relaxed)
{
    if (ToulBar2::verbose >= 1)
        cout << "adding relaxed mdd constraint - hidden decomposition" << endl;

    vector<Variable*> varReverse;
    for (unsigned v = 0; v < getDivVariables().size(); v++) {
        varReverse.push_back(getDivVariables()[getDivVariables().size() - v - 1]);
    }

    unsigned nLayers = varReverse.size();
    vector<Cost> vc;
    EnumeratedVariable* ex;

    EnumeratedVariable* c;
    EnumeratedVariable* h = NULL;
    EnumeratedVariable* hp = NULL;

    int cId;
    int hId = -1;
    int hpId = -1;
    unsigned source;
    unsigned target;

    for (unsigned layer = 0; layer < nLayers; layer++) {
        ex = (EnumeratedVariable*)varReverse[layer];
        int xId = ex->wcspIndex; //wcsp index of current divVariable

        hId = divHVarsId[relaxed][xId]; //index of variable h
        h = (EnumeratedVariable*)getVar(hId);
        cId = divVarsId[relaxed][xId];
        c = (EnumeratedVariable*)getVar(cId);
        vc.clear();
        //add constraint between hp and c, except on layer 0
        if (layer == 0) {
            for (unsigned val_c = 0; val_c < c->getDomainInitSize(); val_c++) {
                source = val_c / (ToulBar2::divWidth);
                vc.push_back((source == 0) ? MIN_COST : getUb());
            }
            postIncrementalUnaryConstraint(cId, vc);
        } else {
            for (unsigned val_c = 0; val_c < c->getDomainInitSize(); val_c++) {
                source = val_c / (ToulBar2::divWidth);
                for (unsigned val_hp = 0; val_hp < hp->getDomainInitSize(); val_hp++) {
                    vc.push_back((val_hp == source) ? MIN_COST : getUb());
                }
            }
            postIncrementalBinaryConstraint(cId, hpId, vc);
        }
        vc.clear();
        //add constraint between c and h
        if (layer != nLayers - 1) {
            for (unsigned val_c = 0; val_c < c->getDomainInitSize(); val_c++) {
                source = val_c / (ToulBar2::divWidth);
                target = val_c % (ToulBar2::divWidth);
                for (unsigned val_h = 0; val_h < h->getDomainInitSize(); val_h++) {
                    vc.push_back((val_h == target) ? MIN_COST : getUb());
                }
            }
            postIncrementalBinaryConstraint(cId, hId, vc);
        }
        vc.clear();
        //Add constraint between c and ex
        for (unsigned val_c = 0; val_c < c->getDomainInitSize(); val_c++) {
            source = val_c / (ToulBar2::divWidth);
            target = val_c % (ToulBar2::divWidth);
            for (unsigned val_ex = 0; val_ex < ex->getDomainInitSize(); val_ex++) {
                vc.push_back((source < mdd[layer].size() && target < mdd[layer][source].size()) ? mdd[layer][source][target][val_ex] : getUb());
            }
        }
        postIncrementalBinaryConstraint(cId, xId, vc);
        hp = h;
        hpId = hId;
    }
}

void WCSP::addTMDDConstraint(Mdd mdd, int relaxed)
{
    if (ToulBar2::verbose >= 1)
        cout << "adding relaxed mdd constraint - ternary decomposition" << endl;

    vector<Variable*> varReverse;
    for (unsigned v = 0; v < getDivVariables().size(); v++) {
        varReverse.push_back(getDivVariables()[getDivVariables().size() - v - 1]);
    }

    unsigned nLayers = varReverse.size();
    vector<Cost> vc;
    EnumeratedVariable* ex;

    EnumeratedVariable* h = NULL;
    EnumeratedVariable* hp = NULL;

    int hId = -1;
    int hpId = -1;

    for (unsigned layer = 0; layer < nLayers; layer++) {
        ex = (EnumeratedVariable*)varReverse[layer];
        int xId = ex->wcspIndex; //wcsp index of current divVariable

        vc.clear();
        // Add constraint between hp ex and h (except on extremities)
        if (layer == 0) { // first position, no hp, only distance 0
            hId = divHVarsId[relaxed][xId]; //index of variable h
            h = (EnumeratedVariable*)getVar(hId);
            for (unsigned val_ex = 0; val_ex < ex->getDomainInitSize(); val_ex++) {
                for (unsigned val_h = 0; val_h < h->getDomainInitSize(); val_h++) {
                    vc.push_back((val_h < mdd[0][0].size()) ? mdd[0][0][val_h][val_ex] : getUb());
                }
            }
            postIncrementalBinaryConstraint(xId, hId, vc);
        } else if (layer + 1 == nLayers) { // last position, no h
            for (unsigned val_hp = 0; val_hp < hp->getDomainInitSize(); val_hp++) {
                for (unsigned val_ex = 0; val_ex < ex->getDomainInitSize(); val_ex++) {
                    vc.push_back((val_hp < mdd[layer].size()) ? mdd[layer][val_hp][0][val_ex] : getUb());
                }
            }
            postIncrementalBinaryConstraint(hpId, xId, vc);
        } else {
            hId = divHVarsId[relaxed][xId]; //index of variable h
            h = (EnumeratedVariable*)getVar(hId);
            for (unsigned val_hp = 0; val_hp < hp->getDomainInitSize(); val_hp++) {
                for (unsigned val_ex = 0; val_ex < ex->getDomainInitSize(); val_ex++) {
                    for (unsigned val_h = 0; val_h < h->getDomainInitSize(); val_h++) {
                        vc.push_back((val_hp < mdd[layer].size() && val_h < mdd[layer][val_hp].size()) ? mdd[layer][val_hp][val_h][val_ex] : getUb());
                    }
                }
            }
            postIncrementalTernaryConstraint(hpId, xId, hId, vc);
        }
        hp = h;
        hpId = hId;
    }
}

void WCSP::postWSum(int* scopeIndex, int arity, string semantics, Cost baseCost, string comparator, int rightRes)
{
#ifndef NDEBUG
    for (int i = 0; i < arity; i++)
        for (int j = i + 1; j < arity; j++)
            assert(scopeIndex[i] != scopeIndex[j]);
#endif
    string gcname = "wsum";
    WeightedSum* decomposableGCF = new WeightedSum(arity, scopeIndex);
    decomposableGCF->setSemantics(semantics);
    decomposableGCF->setBaseCost(baseCost);
    decomposableGCF->setComparator(comparator);
    decomposableGCF->setRightRes(rightRes);
    decomposableGCF->addToCostFunctionNetwork(this);
}

void WCSP::postWVarSum(int* scopeIndex, int arity, string semantics, Cost baseCost, string comparator, int varIndex)
{
#ifndef NDEBUG
    for (int i = 0; i < arity; i++)
        for (int j = i + 1; j < arity; j++)
            assert(scopeIndex[i] != scopeIndex[j]);
#endif
    string gcname = "wvarsum";
    WeightedVarSum* decomposableGCF = new WeightedVarSum(arity, scopeIndex);
    decomposableGCF->setSemantics(semantics);
    decomposableGCF->setBaseCost(baseCost);
    decomposableGCF->setComparator(comparator);
    decomposableGCF->setIndex(varIndex);
    decomposableGCF->addToCostFunctionNetwork(this);
}

void WCSP::postWAmong(int* scopeIndex, int arity, string semantics, Cost baseCost, Value* values, int nbValues, int lb, int ub)
{
#ifndef NDEBUG
    for (int i = 0; i < arity; i++)
        for (int j = i + 1; j < arity; j++)
            assert(scopeIndex[i] != scopeIndex[j]);
#endif
    WeightedAmong* decomposableGCF = new WeightedAmong(arity, scopeIndex);
    decomposableGCF->setSemantics(semantics);
    decomposableGCF->setBaseCost(baseCost);
    for (int i = 0; i < nbValues; i++) {
        decomposableGCF->addValue(values[i]);
    }
    decomposableGCF->setBounds(lb, ub);
    decomposableGCF->addToCostFunctionNetwork(this);
    delete[] values;
    //delete [] decomposableGCF;
}

void WCSP::postWVarAmong(int* scopeIndex, int arity, const string& semantics, Cost baseCost, Value* values, int nbValues, int varIndex)
{
#ifndef NDEBUG
    for (int i = 0; i < arity; i++)
        for (int j = i + 1; j < arity; j++)
            assert(scopeIndex[i] != scopeIndex[j]);
#endif
    WeightedVarAmong* decomposableGCF = new WeightedVarAmong(arity, scopeIndex);
    decomposableGCF->setSemantics(semantics);
    decomposableGCF->setBaseCost(baseCost);
    for (int i = 0; i < nbValues; i++) {
        decomposableGCF->addValue(values[i]);
    }
    decomposableGCF->setIndex(varIndex);
    decomposableGCF->addToCostFunctionNetwork(this);
    delete[] values;
}

void WCSP::postWRegular(int* scopeIndex, int arity, int nbStates, vector<pair<int, Cost>> initial_States, vector<pair<int, Cost>> accepting_States, int** Wtransitions,
    vector<Cost> transitionsCosts)
{
#ifndef NDEBUG
    for (int i = 0; i < arity; i++)
        for (int j = i + 1; j < arity; j++)
            assert(scopeIndex[i] != scopeIndex[j]);
#endif
    WFA* automaton = new WFA(nbStates);
    for (unsigned int i = 0; i < initial_States.size(); i++) {
        automaton->getInitialStates().push_back(initial_States[i]);
    }
    for (unsigned int i = 0; i < accepting_States.size(); i++) {
        automaton->getAcceptingStates().push_back(accepting_States[i]);
    }
    for (unsigned int i = 0; i < transitionsCosts.size(); i++) {
        automaton->getTransitions().push_back(new WTransition(Wtransitions[i][0], Wtransitions[i][1], Wtransitions[i][2], transitionsCosts[i]));
    }
    WeightedRegular regular(arity, scopeIndex);
    regular.setWFA(automaton);
    regular.addToCostFunctionNetwork(this);
}

void WCSP::postWGcc(int* scopeIndex, int arity, string semantics, Cost baseCost, Value* values, int nbValues, int* lb, int* ub)
{
#ifndef NDEBUG
    for (int i = 0; i < arity; i++)
        for (int j = i + 1; j < arity; j++)
            assert(scopeIndex[i] != scopeIndex[j]);
#endif
    WeightedGcc* decomposableGCF = new WeightedGcc(arity, scopeIndex);
    decomposableGCF->setSemantics(semantics);
    decomposableGCF->setBaseCost(baseCost);
    decomposableGCF->setNbValue(nbValues);
    for (int i = 0; i < nbValues; i++) {
        decomposableGCF->setBounds(values[i], lb[i], ub[i]);
    }
    decomposableGCF->addToCostFunctionNetwork(this);
    delete[] values;
    delete[] lb;
    delete[] ub;
}

void WCSP::postWSame(int* scopeIndex, int arity, string semantics, Cost baseCost)
{
#ifndef NDEBUG
    for (int i = 0; i < arity / 2; i++)
        for (int j = i + 1; j < arity / 2; j++)
            assert(scopeIndex[i] != scopeIndex[j]);
    for (int i = arity / 2; i < arity; i++)
        for (int j = i + 1; j < arity; j++)
            assert(scopeIndex[i] != scopeIndex[j]);
    for (int i = 0; i < arity / 2; i++)
        for (int j = arity / 2; j < arity; j++)
            assert(scopeIndex[i] != scopeIndex[j]);
#endif
    WeightedSame* decomposableGCF = new WeightedSame(arity, scopeIndex);
    decomposableGCF->setSemantics(semantics);
    decomposableGCF->setBaseCost(baseCost);
    decomposableGCF->addToCostFunctionNetwork(this);
}

void WCSP::postWAllDiff(int* scopeIndex, int arity, string semantics, Cost baseCost)
{
#ifndef NDEBUG
    for (int i = 0; i < arity; i++)
        for (int j = i + 1; j < arity; j++)
            assert(scopeIndex[i] != scopeIndex[j]);
#endif
    WeightedAllDifferent* decomposableGCF = new WeightedAllDifferent(arity, scopeIndex);
    decomposableGCF->setSemantics(semantics);
    decomposableGCF->setBaseCost(baseCost);
    decomposableGCF->addToCostFunctionNetwork(this);
}

void WCSP::postWSameGcc(int* scopeIndex, int arity, string semantics, Cost baseCost, Value* values, int nbValues, int* lb, int* ub)
{
#ifndef NDEBUG
    for (int i = 0; i < arity / 2; i++)
        for (int j = i + 1; j < arity / 2; j++)
            assert(scopeIndex[i] != scopeIndex[j]);
    for (int i = arity / 2; i < arity; i++)
        for (int j = i + 1; j < arity; j++)
            assert(scopeIndex[i] != scopeIndex[j]);
    for (int i = 0; i < arity / 2; i++)
        for (int j = arity / 2; j < arity; j++)
            assert(scopeIndex[i] != scopeIndex[j]);
#endif
    WeightedSameGcc* decomposableGCF = new WeightedSameGcc(arity, scopeIndex);
    decomposableGCF->setSemantics(semantics);
    decomposableGCF->setBaseCost(baseCost);
    decomposableGCF->setNbValue(nbValues);
    for (int i = 0; i < nbValues; i++) {
        decomposableGCF->setBounds(values[i], lb[i], ub[i]);
    }
    decomposableGCF->addToCostFunctionNetwork(this);
    delete[] values;
    delete[] lb;
    delete[] ub;
}

void WCSP::postWOverlap(int* scopeIndex, int arity, string semantics, Cost baseCost, string comparator, int rightRes)
{
#ifndef NDEBUG
    for (int i = 0; i < arity / 2; i++)
        for (int j = i + 1; j < arity / 2; j++)
            assert(scopeIndex[i] != scopeIndex[j]);
    for (int i = arity / 2; i < arity; i++)
        for (int j = i + 1; j < arity; j++)
            assert(scopeIndex[i] != scopeIndex[j]);
    for (int i = 0; i < arity / 2; i++)
        for (int j = arity / 2; j < arity; j++)
            assert(scopeIndex[i] != scopeIndex[j]);
#endif
    WeightedOverlap* decomposableGCF = new WeightedOverlap(arity, scopeIndex);
    decomposableGCF->setSemantics(semantics);
    decomposableGCF->setBaseCost(baseCost);
    decomposableGCF->setComparator(comparator);
    decomposableGCF->setRightRes(rightRes);
    decomposableGCF->addToCostFunctionNetwork(this);
}

/// \brief create a monolithic global cost function in intension with a particular semantic
/// \param scopeIndex array of enumerated variable indexes (as returned by makeEnumeratedVariable)
/// \param arity size of scopeIndex
/// \param gcname specific \e keyword name of the global cost function (\e eg salldiff, sgcc, sregular, ssame)
/// \param file problem file (\see \ref wcspformat)
/// \deprecated should use postWXXX methods
int WCSP::postGlobalConstraint(int* scopeIndex, int arity, const string& gcname, istream& file, int* constrcounter, bool mult)
{
    if (gcname == "salldiffdp") {
        string semantics;
        Cost baseCost;
        file >> semantics >> baseCost;
        if (mult)
            baseCost *= ToulBar2::costMultiplier;
        postWAllDiff(scopeIndex, arity, semantics, "DAG", baseCost);
        return -1;
    } else if (gcname == "sgccdp") {
        string semantics;
        Cost baseCost;
        int nvalues;
        vector<BoundedObjValue> values;
        file >> semantics >> baseCost >> nvalues;
        for (int i = 0; i < nvalues; i++) {
            int d, high, low;
            file >> d >> low >> high;
            values.push_back(BoundedObjValue(d, high, low));
        }
        if (mult)
            baseCost *= ToulBar2::costMultiplier;
        postWGcc(scopeIndex, arity, semantics, "DAG", baseCost, values);
        return -1;
    }else {
        if(gcname =="knapsack"){
        postKnapsackConstraint(scopeIndex,arity,file);
        return -1;}
    }

    GlobalConstraint* gc = postGlobalCostFunction(scopeIndex, arity, gcname, constrcounter);
    if (gc == NULL)
        return -1;
    if (file) {
        gc->read(file, mult);
    }
    gc->init();
    return gc->wcspIndex;
}

GlobalConstraint* WCSP::postGlobalCostFunction(int* scopeIndex, int arity, const string& gcname, int* constrcounter)
{
#ifndef NDEBUG
    for (int i = 0; i < arity; i++)
        for (int j = i + 1; j < arity; j++)
            assert(scopeIndex[i] != scopeIndex[j]);
#endif
    GlobalConstraint* gc = NULL;

    vector<EnumeratedVariable*> scopeVarsV(arity);
    for (int i = 0; i < arity; i++)
        scopeVarsV[i] = (EnumeratedVariable*)vars[scopeIndex[i]];
    auto scopeVars = scopeVarsV.data();

    if (gcname == "salldiff") {
        gc = new AllDiffConstraint(this, scopeVars, arity);
    } else if (gcname == "sgcc") {
        gc = new GlobalCardinalityConstraint(this, scopeVars, arity);
    } else if (gcname == "ssame") {
        gc = new SameConstraint(this, scopeVars, arity);
    } else if (gcname == "sregular") {
        gc = new RegularFlowConstraint(this, scopeVars, arity);
#ifdef ILOGCPLEX
    } else if (gcname == "slinear") {
        gc = new LPSConstraint(this, scopeVars, arity, constrcounter);
#endif
    } else if (gcname == "samong" || gcname == "samongdp") {
        gc = new AmongConstraint(this, scopeVars, arity);
    } else if (gcname == "sregulardp") {
        gc = new RegularDPConstraint(this, scopeVars, arity);
    } else if (gcname == "sgrammar" || gcname == "sgrammardp") {
        gc = new GrammarConstraint(this, scopeVars, arity);
    } else if (gcname == "MST" || gcname == "smstdp") {
        gc = new TreeConstraint(this, scopeVars, arity);
    } else if (gcname == "max" || gcname == "smaxdp") {
        gc = new MaxConstraint(this, scopeVars, arity);
    } else {
        cout << gcname << " undefined" << endl;
        exit(1);
    }

    if (gc != NULL)
        globalconstrs.push_back(gc);

    return gc;
}

int WCSP::postCliqueConstraint(int* scopeIndex, int arity, istream& file)
{
#ifndef NDEBUG
    for (int i = 0; i < arity; i++)
        for (int j = i + 1; j < arity; j++)
            assert(scopeIndex[i] != scopeIndex[j]);
#endif
    vector<EnumeratedVariable*> scopeVars(arity);
    for (int i = 0; i < arity; i++)
        scopeVars[i] = (EnumeratedVariable*)vars[scopeIndex[i]];
    auto cc = new CliqueConstraint(this, scopeVars.data(), arity);
    cc->read(file);
    if (isDelayedNaryCtr)
        delayedNaryCtr.push_back(cc->wcspIndex);
    else {
        BinaryConstraint* bctr;
        TernaryConstraint* tctr = new TernaryConstraint(this);
        elimTernConstrs.push_back(tctr);
        for (int j = 0; j < 3; j++) {
            if (!ToulBar2::vac)
                bctr = new BinaryConstraint(this);
            else
                bctr = new VACBinaryConstraint(this);
            elimBinConstrs.push_back(bctr);
        }
        cc->propagate();
    }
    return cc->wcspIndex;
}
int WCSP::postKnapsackConstraint(int* scopeIndex, int arity, istream& file)
{
#ifndef NDEBUG
    for (int i = 0; i < arity; i++)
        for (int j = i + 1; j < arity; j++)
            assert(scopeIndex[i] != scopeIndex[j]);
#endif
        //Eliminate variable with weigth 0.
    vector<int> Weightzero;
    vector<Long> weights;
    int nbzeros=0;
    Long readw;
    Long capacity;
    Long MaxWeight=0;
    Long NegCapacity=0;

    file >> capacity;

    bool isclause = (capacity <= 1 && arity > 3);
    Long clausecapacity = 1;
    vector<tValue> clausetuple(arity, 0);

    for (int i = 0; i != arity; ++i) {
        file >> readw;
        isclause &= (readw == 1 || readw == -1);
        if (readw == -1) {
            clausetuple[i] = 1;
            clausecapacity--;
        }
        weights.push_back(readw);
        if(weights[i]==0){
            Weightzero.push_back(i);
            nbzeros++;}
        else if(weights[i]<0)
            NegCapacity-=weights[i];
        else
            MaxWeight+=weights[i];
    }
    isclause &= (clausecapacity == capacity);
    int k=0;
    for(int i=0; i<nbzeros;i++){
        weights.erase(weights.begin() + Weightzero[i]-k);
        for(int j=Weightzero[i]-k;j<arity;j++)
            *(scopeIndex+j)=*(scopeIndex+j+1);
        k++;
        arity--;
    }
    vector<EnumeratedVariable*> scopeVars(arity);
    for (int i = 0; i < arity; i++)
        scopeVars[i] = (EnumeratedVariable*)vars[scopeIndex[i]];
    AbstractNaryConstraint* cc = NULL;
    if (isclause) {
        assert(arity == (int)clausetuple.size());
        if (ToulBar2::verbose >= 3) cout << "Knapsack constraint of arity " << arity << " transformed into clause!" << endl;
        cc = new WeightedClause(this, scopeVars.data(), arity, getUb(), clausetuple);
    } else {
        cc = new KnapsackConstraint(this, scopeVars.data(), arity, capacity, weights, MaxWeight,NegCapacity);
    }
    if (isDelayedNaryCtr)
        delayedNaryCtr.push_back(cc->wcspIndex);
    else {
        BinaryConstraint* bctr;
        TernaryConstraint* tctr = new TernaryConstraint(this);
        elimTernConstrs.push_back(tctr);
        for (int j = 0; j < 3; j++) {
            if (!ToulBar2::vac)
                bctr = new BinaryConstraint(this);
            else
                bctr = new VACBinaryConstraint(this);
            elimBinConstrs.push_back(bctr);
        }
        cc->propagate();
    }
    return cc->wcspIndex;
}

// only DAG-based or network-based propagator
int WCSP::postWAmong(int* scopeIndex, int arity, const string& semantics, const string& propagator, Cost baseCost,
    const vector<Value>& values, int lb, int ub)
{
#ifndef NDEBUG
    for (int i = 0; i < arity; i++)
        for (int j = i + 1; j < arity; j++)
            assert(scopeIndex[i] != scopeIndex[j]);
#endif
    if (propagator == "network") {
        string semantics_ = semantics;
        int* values_ = (int*)&values[0];
        postWAmong(scopeIndex, arity, semantics_, baseCost, values_, values.size(), lb, ub);
        return INT_MIN;
    }

    AmongConstraint* gc = (AmongConstraint*)postGlobalCostFunction(scopeIndex, arity, "samong");
    if (gc == NULL)
        return -1;

    gc->setSemantics(semantics);
    gc->setBaseCost(baseCost);
    gc->setUpperBound(ub);
    gc->setLowerBound(lb);
    for (unsigned int i = 0; i < values.size(); i++)
        gc->addBoundingValue(values[i]);
    gc->init();
    return gc->wcspIndex;
}

int WCSP::postWRegular(int* scopeIndex, int arity, const string& semantics, const string& propagator, Cost baseCost,
    int nbStates,
    const vector<WeightedObjInt>& initial_States,
    const vector<WeightedObjInt>& accepting_States,
    const vector<DFATransition>& Wtransitions)
{
#ifndef NDEBUG
    for (int i = 0; i < arity; i++)
        for (int j = i + 1; j < arity; j++)
            assert(scopeIndex[i] != scopeIndex[j]);
#endif
    if (propagator == "network") { // Warning! semantics not used
        vector<pair<int, Cost>> initial_States_;
        for (unsigned int i = 0; i < initial_States.size(); i++) {
            initial_States_.push_back(pair<int, Cost>(initial_States[i].val, initial_States[i].weight));
        }
        vector<pair<int, Cost>> accepting_States_;
        for (unsigned int i = 0; i < accepting_States.size(); i++) {
            accepting_States_.push_back(pair<int, Cost>(accepting_States[i].val, accepting_States[i].weight));
        }
        vector<Cost> transitionsCosts;
        vector<int*> transitions;
        for (unsigned int i = 0; i < Wtransitions.size(); i++) {
            int* transition = new int[3];
            transition[0] = Wtransitions[i].start;
            transition[1] = Wtransitions[i].end;
            transition[2] = Wtransitions[i].symbol;
            transitions.push_back(transition);
            transitionsCosts.push_back(Wtransitions[i].weight);
        }
        postWRegular(scopeIndex, arity, nbStates, initial_States_, accepting_States_, &transitions[0], transitionsCosts);
        for (unsigned int i = 0; i < Wtransitions.size(); i++) {
            delete[] transitions[i];
        }

        return INT_MIN;
    }

    WeightedAutomaton* wfa = NULL;
    int constrIndex = -1;

    if (propagator == "flow") {
        RegularFlowConstraint* gc = (RegularFlowConstraint*)postGlobalCostFunction(scopeIndex, arity, "sregular");
        if (gc != NULL) {
            constrIndex = gc->wcspIndex;
            gc->setSemantics(semantics);
            gc->setBaseCost(baseCost);
            wfa = gc->getWeightedAutomaton();
        }
    } else {
        RegularDPConstraint* gc = (RegularDPConstraint*)postGlobalCostFunction(scopeIndex, arity, "sregulardp");
        if (gc != NULL) {
            constrIndex = gc->wcspIndex;
            gc->setSemantics(semantics);
            gc->setBaseCost(baseCost);
            wfa = gc->getWeightedAutomaton();
        }
    }

    if (wfa != NULL) {
        wfa->setNumStates(nbStates);
        for (unsigned int i = 0; i < initial_States.size(); i++) {
            wfa->addInitialState(initial_States[i].val);
        }
        for (unsigned int i = 0; i < accepting_States.size(); i++) {
            wfa->addFinalState(accepting_States[i].val);
        }
        for (unsigned int i = 0; i < Wtransitions.size(); i++) {
            wfa->addTransition(Wtransitions[i].start, Wtransitions[i].symbol,
                Wtransitions[i].end, Wtransitions[i].weight);
        }
    }
    if (constrIndex >= 0)
        ((GlobalConstraint*)getCtr(constrIndex))->init();
    return constrIndex;
}

int WCSP::postWGcc(int* scopeIndex, int arity, const string& semantics, const string& propagator, Cost baseCost,
    const vector<BoundedObjValue>& values)
{
#ifndef NDEBUG
    for (int i = 0; i < arity; i++)
        for (int j = i + 1; j < arity; j++)
            assert(scopeIndex[i] != scopeIndex[j]);
#endif
    if (propagator == "network") {
        string semantics_ = semantics;
        int nbValues = values.size();
        Value* values_ = new Value[nbValues];
        int* lb = new int[nbValues];
        int* ub = new int[nbValues];
        for (unsigned int i = 0; i < values.size(); i++) {
            values_[i] = values[i].val;
            lb[i] = values[i].lower;
            ub[i] = values[i].upper;
        }
        postWGcc(scopeIndex, arity, semantics, baseCost, values_, nbValues, lb, ub);
        return INT_MIN;
    }

    if (propagator == "flow") {
        GlobalCardinalityConstraint* gc = (GlobalCardinalityConstraint*)postGlobalCostFunction(scopeIndex, arity, "sgcc");
        if (gc == NULL)
            return -1;

        gc->setSemantics(semantics);
        gc->setBaseCost(baseCost);
        for (unsigned int i = 0; i < values.size(); i++) {
            gc->addValueAndBounds(values[i].val, values[i].lower, values[i].upper);
        }
        gc->init();
        return gc->wcspIndex;
    } else { // DAG-based propagator
        for (unsigned int i = 0; i < values.size(); i++) {
            //Adding a wamong
            vector<Value> values_;
            values_.push_back(values[i].val);
            postWAmong(scopeIndex, arity, semantics, "DAG", baseCost, values_, values[i].lower, values[i].upper);
        }
        return INT_MIN;
    }
}

int WCSP::postWSame(int* scopeIndexG1, int arityG1, int* scopeIndexG2, int arityG2, const string& semantics, const string& propagator, Cost baseCost)
{
    assert(arityG1 >= 2); // does not work for binary or ternary cost functions!!!
    assert(arityG1 == arityG2);
#ifndef NDEBUG
    for (int i = 0; i < arityG1; i++)
        for (int j = i + 1; j < arityG1; j++)
            assert(scopeIndexG1[i] != scopeIndexG1[j]);
    for (int i = 0; i < arityG2; i++)
        for (int j = i + 1; j < arityG2; j++)
            assert(scopeIndexG2[i] != scopeIndexG2[j]);
    for (int i = 0; i < arityG1; i++)
        for (int j = 1; j < arityG2; j++)
            assert(scopeIndexG1[i] != scopeIndexG2[j]);
#endif

    if (propagator == "network") {
        string semantics_ = semantics;
        vector<int> scopeIndex;
        for (int i = 0; i < arityG1; i++)
            scopeIndex.push_back(scopeIndexG1[i]);
        for (int i = 0; i < arityG2; i++)
            scopeIndex.push_back(scopeIndexG2[i]);
        int arity = arityG1 + arityG2;
        postWSame(&scopeIndex[0], arity, semantics_, baseCost);
        return INT_MIN;
    }

    int* scopeIndex = new int[arityG1 + arityG2];
    int arity = arityG1 + arityG2;
    for (int i = 0; i < arityG1; i++) {
        scopeIndex[i] = scopeIndexG1[i];
    }
    for (int i = 0; i < arityG2; i++) {
        scopeIndex[i + arityG1] = scopeIndexG2[i];
    }

    SameConstraint* gc = (SameConstraint*)postGlobalCostFunction(scopeIndex, arity, "ssame");
    if (gc == NULL)
        return -1;

    gc->setSemantics(semantics);
    gc->setBaseCost(baseCost);
    for (int i = 0; i < arityG1; i++)
        gc->addVariablesToGroup((EnumeratedVariable*)vars[scopeIndexG1[i]], 0);
    for (int i = 0; i < arityG2; i++)
        gc->addVariablesToGroup((EnumeratedVariable*)vars[scopeIndexG2[i]], 1);

    delete[] scopeIndex;

    gc->init();
    return gc->wcspIndex;
}

int WCSP::postWAllDiff(int* scopeIndex, int arity, const string& semantics, const string& propagator, Cost baseCost)
{
#ifndef NDEBUG
    for (int i = 0; i < arity; i++)
        for (int j = i + 1; j < arity; j++)
            assert(scopeIndex[i] != scopeIndex[j]);
#endif
    if (propagator == "network") {
        string semantics_ = semantics;
        postWAllDiff(scopeIndex, arity, semantics_, baseCost);
        return INT_MIN;
    }

    if (propagator == "flow") {
        GlobalConstraint* gc = postGlobalCostFunction(scopeIndex, arity, "salldiff");

        if (gc == NULL)
            return -1;

        gc->setSemantics(semantics);
        gc->setBaseCost(baseCost);
        gc->init();
        return gc->wcspIndex;
    } else { // DAG-based propagation using a decomposition into multiple among cost functions
        // Counting the number of value
        int inf = ((EnumeratedVariable*)getVar(scopeIndex[0]))->getInf();
        int sup = ((EnumeratedVariable*)getVar(scopeIndex[0]))->getSup();
        for (int variable = 0; variable < arity; ++variable) {
            int tinf = ((EnumeratedVariable*)getVar(scopeIndex[variable]))->getInf();
            int tsup = ((EnumeratedVariable*)getVar(scopeIndex[variable]))->getSup();
            inf = min(inf, tinf);
            sup = max(sup, tsup);
        }

        // Adding WeightedAmong over each variable
        for (int value = inf; value <= sup; value++) {
            vector<Value> values;
            values.push_back(value);
            postWAmong(scopeIndex, arity, semantics, "DAG", baseCost, values, 0, 1);
        }
        return INT_MIN;
    }
}

int WCSP::postWGrammarCNF(int* scopeIndex, int arity, const string& semantics, const string& propagator, Cost baseCost,
    int nbNonTerminal,
    int startSymbol,
    const vector<CFGProductionRule> WRuleToTerminal)
{

#ifndef NDEBUG
    for (int i = 0; i < arity; i++)
        for (int j = i + 1; j < arity; j++)
            assert(scopeIndex[i] != scopeIndex[j]);
#endif

    GrammarConstraint* gc = (GrammarConstraint*)postGlobalCostFunction(scopeIndex, arity, "sgrammar");

    if (gc == NULL)
        return -1;

    gc->setSemantics(semantics);
    gc->setBaseCost(baseCost);

    WeightedCNFCFG* cnf = gc->getGrammar();

    cnf->setNumNonTerminals(nbNonTerminal);
    cnf->setStartSymbol(startSymbol);

    for (unsigned int i = 0; i < WRuleToTerminal.size(); i++) {
        if (WRuleToTerminal[i].order == 1) {
            cnf->addProduction(WRuleToTerminal[i].from, WRuleToTerminal[i].to[0], 0);
        } else if (WRuleToTerminal[i].order == 2) {
            cnf->addProduction(WRuleToTerminal[i].from, WRuleToTerminal[i].to[0], WRuleToTerminal[i].to[1], 0);
        } else {
            cout << "Warning: either A->v or A->BC is allowed!" << endl;
        }
    }

    gc->init();
    return gc->wcspIndex;
}

int WCSP::postMST(int* scopeIndex, int arity, const string& semantics, const string& propagator, Cost baseCost)
{
#ifndef NDEBUG
    for (int i = 0; i < arity; i++)
        for (int j = i + 1; j < arity; j++)
            assert(scopeIndex[i] != scopeIndex[j]);
#endif
    TreeConstraint* gc = (TreeConstraint*)postGlobalCostFunction(scopeIndex, arity, "MST");
    if (gc == NULL)
        return -1;
    gc->setSemantics(semantics);
    gc->setBaseCost(baseCost);
    gc->init();
    return gc->wcspIndex;
}

int WCSP::postMaxWeight(int* scopeIndex, int arity, const string& semantics, const string& propagator, Cost baseCost,
    const vector<WeightedVarValPair> weightFunction)
{
#ifndef NDEBUG
    for (int i = 0; i < arity; i++)
        for (int j = i + 1; j < arity; j++)
            assert(scopeIndex[i] != scopeIndex[j]);
#endif

    MaxConstraint* gc = (MaxConstraint*)postGlobalCostFunction(scopeIndex, arity, "max");

    if (gc == NULL)
        return -1;

    gc->setSemantics(semantics);
    gc->setBaseCost(baseCost);

    for (unsigned int i = 0; i < weightFunction.size(); i++) {
        gc->setAssignmentWeight((EnumeratedVariable*)vars[weightFunction[i].varIndex],
            weightFunction[i].val,
            weightFunction[i].weight);
    }

    gc->init();
    return gc->wcspIndex;
}

/// \brief add a constant cost to problem lower bound (can be positive or negative)
void WCSP::postNullaryConstraint(Cost cost)
{
    if (cost >= MIN_COST) {
        increaseLb(cost);
    } else {
        decreaseLb(-cost);
    }
}

void WCSP::postNullaryConstraint(Double cost)
{
    postNullaryConstraint((Cost)(round(cost * pow(10, ToulBar2::decimalPoint))));
}

/// \brief add unary costs to enumerated variable \e xIndex
/// \note a unary cost function associated to an enumerated variable is not a Constraint object, it is directly managed inside the EnumeratedVariable class, this is why this function does not return any Constraint index. By doing so, unary costs are better shared inside the cost function network.
void WCSP::postUnary(int xIndex, vector<Cost>& costs)
{
    assert(vars[xIndex]->enumerated());
    EnumeratedVariable* x = (EnumeratedVariable*)vars[xIndex];

    if (ToulBar2::vac) {
        for (unsigned int a = 0; a < x->getDomainInitSize(); a++) {
            Cost c = costs[a];
            histogram(c);
        }
    }
    for (unsigned int a = 0; a < x->getDomainInitSize(); a++) {
        if (costs[a] > MIN_COST)
            x->project(x->toValue(a), costs[a], true);
    }
    x->findSupport();
    x->queueNC();
}

/// \brief add unary costs to interval variable \e xIndex
int WCSP::postUnary(int xIndex, Value* d, int dsize, Cost penalty)
{
    assert(!vars[xIndex]->enumerated());
    Unary* ctr = new Unary(this, (IntervalVariable*)vars[xIndex], d, dsize, penalty);
    return ctr->wcspIndex;
}

/// \brief add unary floating point costs approximated to the current ToulBar2::decimalPoint precision to enumerated variable \e xIndex
/// \note  The ToulBar2::decimapPoint precision must have been set previously for the rest of the CFN existence.
void WCSP::postUnaryConstraint(int xIndex, vector<Double>& dcosts, bool incremental)
{
    assert(vars[xIndex]->enumerated());
    EnumeratedVariable* x = (EnumeratedVariable*)vars[xIndex];

    // normalize the cost function to make it positive
    Double minCost = std::numeric_limits<Double>::infinity();
    for (unsigned int a = 0; a < x->getDomainInitSize(); a++) {
        minCost = min(minCost, dcosts[a]);
    }
    vector<Cost> icosts;
    icosts.reserve(dcosts.size());
    for (unsigned int a = 0; a < x->getDomainInitSize(); a++) {
        icosts[a] = (Cost)(round((dcosts[a] - minCost) * pow(10, ToulBar2::decimalPoint)));
    }
    negCost -= (Cost)(round(minCost * pow(10, ToulBar2::decimalPoint)));
    if (incremental) {
        postIncrementalUnaryConstraint(xIndex, icosts);
    } else {
        postUnaryConstraint(xIndex, icosts);
    }
}

/// \brief create a soft constraint \f$x \geq y + cst\f$ with the associated cost function \f$max( (y + cst - x \leq delta)?(y + cst - x):top , 0 )\f$
int WCSP::postSupxyc(int xIndex, int yIndex, Value cst, Value delta)
{
    assert(xIndex != yIndex);
    if (!vars[xIndex]->enumerated() && !vars[yIndex]->enumerated()) {
        Supxyc* ctr = new Supxyc(this, (IntervalVariable*)vars[xIndex], (IntervalVariable*)vars[yIndex], cst, delta);
        return ctr->wcspIndex;
    } else if (vars[xIndex]->enumerated() && vars[yIndex]->enumerated()) {
        EnumeratedVariable* x = (EnumeratedVariable*)vars[xIndex];
        EnumeratedVariable* y = (EnumeratedVariable*)vars[yIndex];
        vector<Cost> costs;
        for (unsigned int a = 0; a < x->getDomainInitSize(); a++) {
            for (unsigned int b = 0; b < y->getDomainInitSize(); b++) {
                costs.push_back(max((y->toValue(b) + cst - x->toValue(a) <= delta) ? ((Cost)((y->toValue(b) + cst - x->toValue(a)) * ToulBar2::costMultiplier)) : getUb(), MIN_COST));
            }
        }
        return postBinaryConstraint(xIndex, yIndex, costs);
    } else {
        cerr << "Cannot mix variables with interval and enumerated domains!!!" << endl;
        exit(EXIT_FAILURE);
    }
}

/// \brief soft binary disjunctive constraint \f$x \geq y + csty \vee y \geq x + cstx\f$ with associated cost function \f$(x \geq y + csty \vee y \geq x + cstx)?0:penalty\f$
int WCSP::postDisjunction(int xIndex, int yIndex, Value cstx, Value csty, Cost penalty)
{
    assert(xIndex != yIndex);
    if (!vars[xIndex]->enumerated() && !vars[yIndex]->enumerated()) {
        Disjunction* ctr = new Disjunction(this, (IntervalVariable*)vars[xIndex], (IntervalVariable*)vars[yIndex], cstx, csty, penalty);
        return ctr->wcspIndex;
    } else if (vars[xIndex]->enumerated() && vars[yIndex]->enumerated()) {
        EnumeratedVariable* x = (EnumeratedVariable*)vars[xIndex];
        EnumeratedVariable* y = (EnumeratedVariable*)vars[yIndex];
        vector<Cost> costs;
        for (unsigned int a = 0; a < x->getDomainInitSize(); a++) {
            for (unsigned int b = 0; b < y->getDomainInitSize(); b++) {
                costs.push_back(((x->toValue(a) >= y->toValue(b) + csty) || (y->toValue(b) >= x->toValue(a) + cstx)) ? MIN_COST : penalty);
            }
        }
        return postBinaryConstraint(xIndex, yIndex, costs);
    } else {
        cerr << "Cannot mix variables with interval and enumerated domains!!!" << endl;
        exit(EXIT_FAILURE);
    }
}

/// \brief special disjunctive constraint with three implicit hard constraints \f$x \leq xinfty\f$ and \f$y \leq yinfty\f$ and \f$x < xinfty \wedge y < yinfty \Rightarrow (x \geq y + csty \vee y \geq x + cstx)\f$ and an additional cost function \f$((x = xinfty)?costx:0) + ((y= yinfty)?costy:0)\f$
int WCSP::postSpecialDisjunction(int xIndex, int yIndex, Value cstx, Value csty, Value xinfty, Value yinfty, Cost costx, Cost costy)
{
    assert(xIndex != yIndex);
    if (!vars[xIndex]->enumerated() && !vars[yIndex]->enumerated()) {
        SpecialDisjunction* ctr = new SpecialDisjunction(this, (IntervalVariable*)vars[xIndex], (IntervalVariable*)vars[yIndex], cstx, csty, xinfty, yinfty, costx, costy);
        return ctr->wcspIndex;
    } else if (vars[xIndex]->enumerated() && vars[yIndex]->enumerated()) {
        EnumeratedVariable* x = (EnumeratedVariable*)vars[xIndex];
        EnumeratedVariable* y = (EnumeratedVariable*)vars[yIndex];
        vector<Cost> costs;
        for (unsigned int a = 0; a < x->getDomainInitSize(); a++) {
            for (unsigned int b = 0; b < y->getDomainInitSize(); b++) {
                costs.push_back((x->toValue(a) <= xinfty && y->toValue(b) <= yinfty && (x->toValue(a) == xinfty || y->toValue(b) == yinfty || (x->toValue(a) >= y->toValue(b) + csty || y->toValue(b) >= x->toValue(a) + cstx))) ? (((x->toValue(a) == xinfty) ? costx : MIN_COST) + ((y->toValue(b) == yinfty) ? costy : MIN_COST)) : getUb());
            }
        }
        return postBinaryConstraint(xIndex, yIndex, costs);
    } else {
        cerr << "Cannot mix variables with interval and enumerated domains!!!" << endl;
        exit(EXIT_FAILURE);
    }
}

void WCSP::sortConstraints()
{
    for (vector<int>::iterator idctr = delayedNaryCtr.begin(); idctr != delayedNaryCtr.end(); ++idctr) {
        BinaryConstraint* bctr;
        TernaryConstraint* tctr = new TernaryConstraint(this);
        elimTernConstrs.push_back(tctr);
        for (int j = 0; j < 3; j++) {
            if (!ToulBar2::vac)
                bctr = new BinaryConstraint(this);
            else
                bctr = new VACBinaryConstraint(this);
            elimBinConstrs.push_back(bctr);
        }
    }
    for (vector<int>::iterator idctr = delayedNaryCtr.begin(); idctr != delayedNaryCtr.end(); ++idctr) {
        getCtr(*idctr)->propagate();
    }
    delayedNaryCtr.clear();
    isDelayedNaryCtr = false;

    // replaces small arity global cost functions in intension by their equivalent cost functions in extension
    unsigned int i = 0;
    while (i < globalconstrs.size()) {
        if (globalconstrs[i]->connected()) {
            if (globalconstrs[i]->arity() <= 3) {
                globalconstrs[i]->projectNaryBeforeSearch(); // deconnect the current element
                if (i < globalconstrs.size() - 1) {
                    globalconstrs[i] = globalconstrs[globalconstrs.size() - 1]; // and replace it by the last element
                }
                globalconstrs.pop_back(); // decrease vector size by one
            } else {
                i++;
            }
        } else {
            if (i < globalconstrs.size() - 1) {
                globalconstrs[i] = globalconstrs[globalconstrs.size() - 1]; // replace the current element by the last element
            }
            globalconstrs.pop_back(); // decrease vector size by one
        }
    }

    if (ToulBar2::Berge_Dec > 0) {
        // flag pour indiquer si une variable a deja ete visitee initialement a faux
        vector<bool> marked(numberOfVariables(), false);
        vector<int> revdac;
        // nouvel ordre DAC inverse
        //	for (int i = numberOfVariables()-1 ; i >= 0; i--) { if (!marked[i]){ visit(i,revdac,marked,listofsuccessors); }}
        //	for (unsigned int i = 0; i <  numberOfVariables(); i++) { if (!marked[i]){ visit(i,revdac,marked,listofsuccessors); }}

        //	for (unsigned int i = 0; i< numberOfVariables(); i++) {
        //		cout << "listofsuccessors(" << i << "):";
        //		for (unsigned int a = 0; a < listofsuccessors[i].size(); a++) {
        //			cout << " " << listofsuccessors[i][a];
        //		}
        //		cout << endl;
        //	}

        //Mark native variable
        for (int i = ((ToulBar2::nbDecisionVars > 0) ? ToulBar2::nbDecisionVars : numberOfVariables()) - 1; i >= 0; i--) {
            if (!marked[i] && getName(i)[0] != IMPLICIT_VAR_TAG[0]) {
                visit(i, revdac, marked, listofsuccessors);
            }
        }
        //Mark q variable only
        for (int i = numberOfVariables() - 1; i >= 0; i--) {
            if (!marked[i]) {
                visit(i, revdac, marked, listofsuccessors);
            }
        }

        // listofsuccessors.clear(); // appel a la methode clear de l'objet vector

        if (ToulBar2::verbose >= 1) {
            cout << "BERGE DAC reverse order:";
            for (unsigned int i = 0; i < numberOfVariables(); i++) {
                cout << " " << revdac[i];
            }
            cout << endl;
        }

        assert(revdac.size() == numberOfVariables());

        setDACOrder(revdac);
    }
    // postpone costly variable elimination heuristics if too many variables
    if (ToulBar2::varOrder && (numberOfVariables() < LARGE_NB_VARS || reinterpret_cast<uintptr_t>(ToulBar2::varOrder) == MAX_CARD || reinterpret_cast<uintptr_t>(ToulBar2::varOrder) >= ELIM_MAX)) {
        vector<int> order;
        if (isAlreadyTreeDec(ToulBar2::varOrder))
            treeDecFile2Vector(ToulBar2::varOrder, order);
        else
            elimOrderFile2Vector(ToulBar2::varOrder, order);
        setDACOrder(order);
    }
    for (unsigned int i = 0; i < vars.size(); i++) {
        vars[i]->sortConstraints();
    }
    AC.sort(false); // sort in decreasing order to get the smallest DAC index first when doing pop() on this queue
    DAC.sort(true); // sort in increasing order to get the largest DAC index first when doing pop() on this queue
    EAC1.sort(true); // sort in increasing order to get the smallest DAC index first when doing pop() on the EAC2 queue
}

void WCSP::updateCurrentVarsId()
{
    int pos = 0;
    for (unsigned int i = 0; i < vars.size(); i++) {
        if (vars[i]->unassigned()) {
            vars[i]->setCurrentVarId(pos);
            pos++;
        }
    }
}

bool cmpTernaryConstraint(TernaryConstraint* c1, TernaryConstraint* c2)
{
    int v1 = c1->getVar(c1->getDACScopeIndex())->getDACOrder();
    int v2 = c2->getVar(c2->getDACScopeIndex())->getDACOrder();
    if (v1 < v2)
        return true;
    else if (v1 == v2) {
        v1 = min(c1->getVar((c1->getDACScopeIndex() + 1) % 3)->getDACOrder(), c1->getVar((c1->getDACScopeIndex() + 2) % 3)->getDACOrder());
        v2 = min(c2->getVar((c2->getDACScopeIndex() + 1) % 3)->getDACOrder(), c2->getVar((c2->getDACScopeIndex() + 2) % 3)->getDACOrder());
        if (v1 < v2)
            return true;
        else if (v1 == v2) {
            v1 = max(c1->getVar((c1->getDACScopeIndex() + 1) % 3)->getDACOrder(), c1->getVar((c1->getDACScopeIndex() + 2) % 3)->getDACOrder());
            v2 = max(c2->getVar((c2->getDACScopeIndex() + 1) % 3)->getDACOrder(), c2->getVar((c2->getDACScopeIndex() + 2) % 3)->getDACOrder());
            return (v1 < v2);
        }
    }
    return false;
}

void WCSP::processTernary()
{
    // double maxtight = 0;
    // Variable* var;
    // TernaryConstraint *tctr1max = NULL, *tctr2max = NULL;

    // if (ToulBar2::preprocessTernaryHeuristic) {

    //     for (unsigned int i=0; i<vars.size(); i++) {
    // 		TernaryConstraint *tctr1, *tctr2;
    //         double tight = vars[i]->strongLinkedby(var,tctr1,tctr2);
    //         if(tight > maxtight) { maxtight = tight; tctr1max = tctr1; tctr2max = tctr2; }
    // 		if(tctr1 && tctr2 && (tctr1 != tctr2)) {
    // 			tctr1->extendTernary();
    // 			tctr2->extendTernary();
    // 			BinaryConstraint* b = tctr1->commonBinary(tctr2);
    // 			if(b->connected())
    // 			{
    //     			tctr1->projectTernaryBinary(b);
    //     			tctr2->projectTernaryBinary(b);
    //     			b->propagate();
    // 			}
    // 		}
    //     }

    //     if(ToulBar2::verbose > 0) {
    //         cout << "Strongest part has mean cost: " << maxtight;
    //         if(var) cout << "  Variable: " << var->wcspIndex;
    //         if(tctr1max) cout << ", 1. ternary with tight: " << tctr1max->getTightness();
    //         if(tctr2max) cout << ", 2. ternary with tight: " << tctr2max->getTightness();
    //         cout << endl;
    //     }
    // }

    vector<TernaryConstraint*> ternaries;
    for (unsigned int i = 0; i < constrs.size(); i++)
        if (constrs[i]->connected() && !constrs[i]->isSep() && constrs[i]->isTernary()) {
            TernaryConstraint* t = (TernaryConstraint*)constrs[i];
            ternaries.push_back(t);
        }
    for (int i = 0; i < elimTernOrder; i++)
        if (elimTernConstrs[i]->connected()) {
            TernaryConstraint* t = (TernaryConstraint*)elimTernConstrs[i];
            ternaries.push_back(t);
        }
    sort(ternaries.begin(), ternaries.end(), cmpTernaryConstraint);
    for (int i = ternaries.size() - 1; i >= 0; i--) {
        TernaryConstraint* t = ternaries[i];
        //		cout << "PROJECT&SUBTRACT tern(" << t->getVar(0)->getName() << "," << t->getVar(1)->getName() << "," << t->getVar(2)->getName() << ")" << endl;
        t->extendTernary();
        if (ToulBar2::costfuncSeparate)
            t->decompose();
        if (t->connected())
            t->projectTernary();
    }
}

/// \defgroup preprocessing Preprocessing techniques
/// Depending on toulbar2 options, the sequence of preprocessing techniques applied before the search is:
/// -# \e i-bounded variable elimination with user-defined \e i bound
/// -# pairwise decomposition of cost functions (binary cost functions are implicitly decomposed by soft AC and empty cost function removals)
/// -# MinSumDiffusion propagation (see VAC)
/// -# projects\&substracts n-ary cost functions in extension on all the binary cost functions inside their scope (3 < n < max, see toulbar2 options)
/// -# functional variable elimination (see \ref varelim)
/// -# projects\&substracts ternary cost functions in extension on their three binary cost functions inside their scope (before that, extends the existing binary cost functions to the ternary cost function and applies pairwise decomposition)
/// -# creates new ternary cost functions for all triangles (\e ie occurences of three binary cost functions \e xy, \e yz, \e zx)
/// -# removes empty cost functions while repeating #1 and #2 until no new cost functions can be removed
///
/// \note the propagation loop is called after each preprocessing technique (see \ref WCSP::propagate)

void WCSP::preprocessing()
{
    Cost previouslb = getLb();

    Eliminate.clear();
    if (ToulBar2::elimDegree_preprocessing <= -3) {
        int deg = medianDegree();
        Double domsize = medianDomainSize();
        Double size = (Double)numberOfUnassignedVariables() * (sizeof(tValue) * deg + sizeof(Cost)) * Pow(domsize, deg + 1);
        if (ToulBar2::debug >= 2)
            cout << "MAX ESTIMATED ELIM SIZE: " << size << endl;
        assert(ToulBar2::elimSpaceMaxMB > 0);
        if (deg >= 3 && deg <= -ToulBar2::elimDegree_preprocessing && size < (Double)ToulBar2::elimSpaceMaxMB * 1024. * 1024.) {
            ToulBar2::elimDegree_preprocessing = deg;
            if (ToulBar2::verbose >= 0)
                cout << "Generic variable elimination of degree " << deg << endl;
        } else {
            ToulBar2::elimDegree_preprocessing = -1;
            if (ToulBar2::verbose >= 0)
                cout << "Generic variable elimination disabled." << endl;
        }
    }
    if (ToulBar2::elimDegree >= 0 || ToulBar2::elimDegree_preprocessing >= 0 || ToulBar2::preprocessFunctional > 0) {
        initElimConstrs();
        if (ToulBar2::elimDegree_preprocessing >= 0) {
            if (ToulBar2::verbose >= 1)
                cout << "Variable elimination in preprocessing of true degree <= "
                     << ToulBar2::elimDegree_preprocessing << endl;
            ToulBar2::elimDegree_preprocessing_ = ToulBar2::elimDegree_preprocessing;
            maxDegree = -1;
            propagate();
            if (ToulBar2::verbose >= 0)
                cout << "Maximum degree of generic variable elimination: " << maxDegree << endl;
        } else if (ToulBar2::elimDegree >= 0) {
            ToulBar2::elimDegree_ = ToulBar2::elimDegree;
        }
    }

    int posConstrs = 0;
    int posElimTernConstrs = 0;
    if (ToulBar2::costfuncSeparate) {
        double time = cpuTime();
        for (unsigned int i = 0; i < constrs.size(); i++) {
            if (constrs[i]->connected() && !constrs[i]->isSep()) {
                constrs[i]->decompose();
            }
        }
        for (int i = 0; i < elimTernOrder; i++) {
            if (elimTernConstrs[i]->connected() && !elimTernConstrs[i]->isSep()) {
                elimTernConstrs[i]->decompose();
            }
        }
        posConstrs = constrs.size();
        posElimTernConstrs = elimTernOrder;
        if (ToulBar2::verbose >= 0)
            cout << "Cost function decomposition time : " << cpuTime() - time << " seconds.\n";
    }

    propagate();

    // recompute current DAC order and its reverse
    if (ToulBar2::varOrder && numberOfVariables() >= LARGE_NB_VARS && numberOfUnassignedVariables() < LARGE_NB_VARS && reinterpret_cast<uintptr_t>(ToulBar2::varOrder) >= MIN_DEGREE && reinterpret_cast<uintptr_t>(ToulBar2::varOrder) <= APPROX_MIN_DEGREE) {
        vector<int> order;
        if (isAlreadyTreeDec(ToulBar2::varOrder))
            treeDecFile2Vector(ToulBar2::varOrder, order);
        else
            elimOrderFile2Vector(ToulBar2::varOrder, order);
        setDACOrder(order);
    } else {
#ifdef BOOST
        if (ToulBar2::MSTDAC) {
            vector<int> order;
            spanningTreeOrderingBGL(order);
            setDACOrder(order);
        }
#endif
    }
    vector<int> elimorder(numberOfVariables(), -1);
    vector<int> revelimorder(numberOfVariables(), -1);
    for (unsigned int i = 0; i < numberOfVariables(); i++) {
        revelimorder[getVar(i)->getDACOrder()] = i;
        elimorder[numberOfVariables() - getVar(i)->getDACOrder() - 1] = i;
    }
    //	cout << "DAC:";
    //	for (int i = 0; i < numberOfVariables(); i++) {
    //		cout << " " << revelimorder[i];
    //	}
    //	cout << endl;
    //	cout << "REVDAC:";
    //	for (int i = 0; i < numberOfVariables(); i++) {
    //		cout << " " << elimorder[i];
    //	}
    //	cout << endl;
    do {
        previouslb = getLb();
        setDACOrder(revelimorder);
        setDACOrder(elimorder);
        if (ToulBar2::verbose >= 0 && getLb() > previouslb) {
            if (ToulBar2::uai)
                cout << "Reverse DAC dual bound: " << std::fixed << std::setprecision(ToulBar2::decimalPoint) << getDDualBound() << std::setprecision(DECIMAL_POINT) << " energy: " << -(Cost2LogProb(getLb()) + ToulBar2::markov_log) << " (+" << 100. * (getLb() - previouslb) / getLb() << "%)" << endl;
            else
                cout << "Reverse DAC dual bound: " << std::fixed << std::setprecision(ToulBar2::decimalPoint) << getDDualBound() << std::setprecision(DECIMAL_POINT) << " (+" << 100. * (getLb() - previouslb) / getLb() << "%)" << endl;
        }
    } while (getLb() > previouslb && 100. * (getLb() - previouslb) / getLb() > 0.5);

    if (ToulBar2::trwsAccuracy >= 0)
        propagateTRWS();

    if (ToulBar2::preprocessNary > 0) {
        for (unsigned int i = 0; i < constrs.size(); i++) {
            if (constrs[i]->connected() && !constrs[i]->isSep() && constrs[i]->isNary() && constrs[i]->arity() >= 3 && constrs[i]->arity() <= ToulBar2::preprocessNary) {
                NaryConstraint* nary = (NaryConstraint*)constrs[i];
                Long nbtuples = nary->getDomainSizeProduct();
                if ((nbtuples < MAX_NB_TUPLES || nary->size() >= nbtuples) && (nary->size() >= 2 || nary->getDefCost() > MIN_COST)) {
                    nary->keepAllowedTuples(getUb()); // can be very slow!
                    nary->preprojectall2();
                    //			if (nary->connected() && nary->size() >= 4) nary->preproject3();
                }
            }
        }
        //		processTernary();
        propagate();
    }

    // Merge functional (binary bijection only) variables in decision variables
    bool merged = (ToulBar2::preprocessFunctional > 0);
    while (merged) {
        merged = false;
        for (unsigned int i = 0; i < constrs.size(); i++)
            if (constrs[i]->connected() && !constrs[i]->isSep() && constrs[i]->isBinary()
                && (ToulBar2::preprocessFunctional == 1 || constrs[i]->getVar(0)->getDomainSize() == constrs[i]->getVar(1)->getDomainSize())) {
                BinaryConstraint* xy = (BinaryConstraint*)constrs[i];
                EnumeratedVariable* x = (EnumeratedVariable*)xy->getVar(0);
                EnumeratedVariable* y = (EnumeratedVariable*)xy->getVar(1);
                map<Value, Value> functional;
                if (xy->isFunctional(x, y, functional) && y->canbeMerged(x)) {
                    y->mergeTo(xy, functional);
                    merged = true;
                    propagate();
                } else if (xy->isFunctional(y, x, functional) && x->canbeMerged(y)) {
                    x->mergeTo(xy, functional);
                    merged = true;
                    propagate();
                }
            }
        for (int i = 0; i < elimBinOrder; i++)
            if (elimBinConstrs[i]->connected() && !elimBinConstrs[i]->isSep()
                && elimBinConstrs[i]->getVar(0)->getDomainSize() == elimBinConstrs[i]->getVar(1)->getDomainSize()) {
                assert(elimBinConstrs[i]->isBinary());
                BinaryConstraint* xy = (BinaryConstraint*)elimBinConstrs[i];
                EnumeratedVariable* x = (EnumeratedVariable*)xy->getVar(0);
                EnumeratedVariable* y = (EnumeratedVariable*)xy->getVar(1);
                map<Value, Value> functional;
                if (xy->isFunctional(x, y, functional) && y->canbeMerged(x)) {
                    y->mergeTo(xy, functional);
                    merged = true;
                    propagate();
                } else if (xy->isFunctional(y, x, functional) && x->canbeMerged(y)) {
                    x->mergeTo(xy, functional);
                    merged = true;
                    propagate();
                }
            }
    }

    if (ToulBar2::preprocessTernaryRPC) {
        do {
            previouslb = getLb();
            ternaryCompletion();
            setDACOrder(revelimorder);
            processTernary();
            propagate();
            setDACOrder(elimorder);
            processTernary();
            propagate();
            if (ToulBar2::verbose >= 0 && getLb() > previouslb) {
                if (ToulBar2::uai)
                    cout << "PIC dual bound: " << std::fixed << std::setprecision(ToulBar2::decimalPoint) << getDDualBound() << std::setprecision(DECIMAL_POINT) << " energy: " << -(Cost2LogProb(getLb()) + ToulBar2::markov_log) << " (+" << 100. * (getLb() - previouslb) / getLb() << "%, " << numberOfConstraints() << " cost functions)" << endl;
                else
                    cout << "PIC dual bound: " << std::fixed << std::setprecision(ToulBar2::decimalPoint) << getDDualBound() << std::setprecision(DECIMAL_POINT) << " (+" << 100. * (getLb() - previouslb) / getLb() << "%, " << numberOfConstraints() << " cost functions)" << endl;
            }
        } while (getLb() > previouslb && 100. * (getLb() - previouslb) / getLb() > 0.5);
    } else if (ToulBar2::preprocessNary > 0) {
        processTernary();
        propagate();
    }

    if (ToulBar2::minsumDiffusion && ToulBar2::vac)
        vac->minsumDiffusion();
    if (ToulBar2::vac) {
        if (ToulBar2::verbose >= 1)
            cout << "Preprocessing ";
        vac->printStat(true);
        for (unsigned int i = 0; i < vars.size(); i++)
            vars[i]->queueEliminate();
        propagate();
    }

    // Deconnect empty cost functions
    unsigned int nbunvar;
    do {
        nbunvar = numberOfUnassignedVariables();
        for (unsigned int i = 0; i < constrs.size(); i++)
            if (constrs[i]->connected() && !constrs[i]->isSep() && constrs[i]->universal()) {
                if (ToulBar2::verbose >= 3)
                    cout << "deconnect empty cost function: " << *constrs[i];
                constrs[i]->deconnect(true);
            }
        for (int i = 0; i < elimBinOrder; i++)
            if (elimBinConstrs[i]->connected() && !elimBinConstrs[i]->isSep() && elimBinConstrs[i]->universal()) {
                if (ToulBar2::verbose >= 3)
                    cout << "deconnect empty cost function: " << *elimBinConstrs[i];
                elimBinConstrs[i]->deconnect(true);
            }
        for (int i = 0; i < elimTernOrder; i++)
            if (elimTernConstrs[i]->connected() && !elimTernConstrs[i]->isSep() && elimTernConstrs[i]->universal()) {
                if (ToulBar2::verbose >= 3)
                    cout << "deconnect empty cost function: " << *elimTernConstrs[i];
                elimTernConstrs[i]->deconnect(true);
            }
        if (ToulBar2::costfuncSeparate) {
            for (unsigned int i = posConstrs; i < constrs.size(); i++) {
                if (constrs[i]->connected() && !constrs[i]->isSep()) {
                    constrs[i]->decompose();
                }
            }
            for (int i = posElimTernConstrs; i < elimTernOrder; i++) {
                if (elimTernConstrs[i]->connected() && !elimTernConstrs[i]->isSep()) {
                    elimTernConstrs[i]->decompose();
                }
            }
            posConstrs = constrs.size();
            posElimTernConstrs = elimTernOrder;
        }
        propagate();
    } while (numberOfUnassignedVariables() < nbunvar);

    if (ToulBar2::elimDegree >= 0 || ToulBar2::elimDegree_preprocessing >= 0 || ToulBar2::preprocessFunctional > 0) {
        ToulBar2::elimDegree_preprocessing_ = -1;
        if (ToulBar2::elimDegree >= 0) {
            ToulBar2::elimDegree_ = ToulBar2::elimDegree;
            if (ToulBar2::elimDegree_preprocessing < min(2, ToulBar2::elimDegree)) {
                for (int i = numberOfVariables() - 1; i >= 0; i--)
                    vars[i]->queueEliminate();
                propagate();
            }
            if (ToulBar2::verbose >= 1)
                cout << "Variable elimination during search of degree <= "
                     << ToulBar2::elimDegree << endl;
        }
    }

#ifdef BOOST
    if (ToulBar2::MSTDAC) {
        vector<int> order;
        spanningTreeOrderingBGL(order);
        setDACOrder(order);
    }
#endif
    if ((ToulBar2::vac && ToulBar2::useRASPS) || (ToulBar2::vac && ToulBar2::VACthreshold)) {
        ToulBar2::RASPSsaveitThresholds = true;
        ToulBar2::RASPSitThresholds.clear();
        propagate();
        ToulBar2::RASPSsaveitThresholds = false;
        ToulBar2::RASPSlastitThreshold = vac->RASPSFindItThreshold();
        if (ToulBar2::VACthreshold || ToulBar2::RASPSangle < 0) {
            if (ToulBar2::verbose >= 0)
                cout << "RASPS/VAC threshold: " << ToulBar2::RASPSlastitThreshold << endl;
            ToulBar2::costThreshold = ToulBar2::RASPSlastitThreshold;
        } else {
            if (ToulBar2::verbose >= 0)
                cout << "RASPS threshold: " << ToulBar2::RASPSlastitThreshold << endl;
        }
    }
    if (ToulBar2::FullEAC) {
        for (unsigned int i = 0; i < vars.size(); i++) {
            if (vars[i]->unassigned())
                vars[i]->isEAC(); // update EAC supports using values from best known solution if provided
        }
        for (unsigned int i = 0; i < vars.size(); i++) {
            if (vars[i]->unassigned())
                vars[i]->isEAC(); // update Full EAC status using previously updated EAC supports (also used to check all constraints without queueing in FEAC queue)
        }
    }
}

Cost WCSP::finiteUb() const
{
    Cost summaxcost = getLb() + UNIT_COST;
    for (unsigned int i = 0; i < constrs.size(); i++) {
        if (constrs[i]->connected() && !constrs[i]->isSep()) {
            summaxcost += constrs[i]->getMaxFiniteCost();
            if (summaxcost >= MAX_COST)
                return MAX_COST;
        }
    }
    for (int i = 0; i < elimBinOrder; i++) {
        if (elimBinConstrs[i]->connected() && !elimBinConstrs[i]->isSep()) {
            summaxcost += elimBinConstrs[i]->getMaxFiniteCost();
            if (summaxcost >= MAX_COST)
                return MAX_COST;
        }
    }
    for (int i = 0; i < elimTernOrder; i++) {
        if (elimTernConstrs[i]->connected() && !elimTernConstrs[i]->isSep()) {
            summaxcost += elimTernConstrs[i]->getMaxFiniteCost();
            if (summaxcost >= MAX_COST)
                return MAX_COST;
        }
    }
    for (unsigned int i = 0; i < vars.size(); i++) {
        if (NC.empty()) {
            summaxcost += vars[i]->getMaxCost();
        } else {
            if (vars[i]->enumerated()) {
                Cost maxcost = MIN_COST;
                EnumeratedVariable* var = (EnumeratedVariable*)vars[i];
                for (EnumeratedVariable::iterator iter = var->begin(); iter != var->end(); ++iter) {
                    if (var->getCost(*iter) > maxcost)
                        maxcost = var->getCost(*iter);
                }
                summaxcost += maxcost;
            } else {
                summaxcost += max(vars[i]->getInfCost(), vars[i]->getSupCost());
            }
        }
        if (summaxcost >= MAX_COST)
            return MAX_COST;
    }
    return summaxcost;
}

void WCSP::setInfiniteCost()
{
    assert(Store::getDepth() == 0);
    Cost ub = getUb() - getLb();
    assert(ub > 0);
    if (ToulBar2::verbose >= 1)
        cout << "Set infinite cost to " << ub << endl;
    for (unsigned int i = 0; i < constrs.size(); i++) {
        if (constrs[i]->connected() && !constrs[i]->isSep()) {
            constrs[i]->setInfiniteCost(ub);
        }
    }
    for (int i = 0; i < elimBinOrder; i++) {
        if (elimBinConstrs[i]->connected() && !elimBinConstrs[i]->isSep()) {
            elimBinConstrs[i]->setInfiniteCost(ub);
        }
    }
    for (int i = 0; i < elimTernOrder; i++) {
        if (elimTernConstrs[i]->connected() && !elimTernConstrs[i]->isSep()) {
            elimTernConstrs[i]->setInfiniteCost(ub);
        }
    }
}

int WCSP::getMaxCurrentDomainSize() const
{
    int max = (vars.size() > 0) ? 1 : 0;
    for (unsigned int i = 0; i < vars.size(); i++) {
        if (vars[i]->unassigned()) {
            int sz = vars[i]->getDomainSize();
            if (sz > max)
                max = sz;
        }
    }
    return max;
}

unsigned int WCSP::getDomainSizeSum() const
{
    unsigned int sum = 0;
    for (unsigned int i = 0; i < vars.size(); i++) {
        if (vars[i]->unassigned())
            sum += vars[i]->getDomainSize();
    }
    return sum;
}

bool WCSP::getEnumDomain(int varIndex, Value* array)
{
    if (EnumeratedVariable* var = dynamic_cast<EnumeratedVariable*>(vars[varIndex])) {
        var->getDomain(array);
        return true;
    } else
        return false;
}

bool WCSP::getEnumDomainAndCost(int varIndex, ValueCost* array)
{
    if (EnumeratedVariable* var = dynamic_cast<EnumeratedVariable*>(vars[varIndex])) {
        var->getDomainAndCost(array);
        return true;
    } else
        return false;
}

unsigned int WCSP::numberOfConnectedConstraints() const
{
    int res = 0;
    for (unsigned int i = 0; i < constrs.size(); i++)
        if (constrs[i]->connected() && !constrs[i]->isSep())
            res++;
    for (int i = 0; i < elimBinOrder; i++)
        if (elimBinConstrs[i]->connected() && !elimBinConstrs[i]->isSep())
            res++;
    for (int i = 0; i < elimTernOrder; i++)
        if (elimTernConstrs[i]->connected() && !elimTernConstrs[i]->isSep())
            res++;
    return res;
}

unsigned int WCSP::medianArity() const
{
    unsigned int nb = numberOfConnectedConstraints();
    if (nb == 0)
        return 0;
    int arity[nb];
    unsigned int pos = 0;
    for (unsigned int i = 0; i < constrs.size(); i++)
        if (constrs[i]->connected() && !constrs[i]->isSep()) {
            arity[pos] = constrs[i]->arity();
            pos++;
        }
    for (int i = 0; i < elimBinOrder; i++)
        if (elimBinConstrs[i]->connected() && !elimBinConstrs[i]->isSep()) {
            arity[pos] = 2;
            pos++;
        }
    for (int i = 0; i < elimTernOrder; i++)
        if (elimTernConstrs[i]->connected() && !elimTernConstrs[i]->isSep()) {
            arity[pos] = 3;
            pos++;
        }
    assert(pos == numberOfConnectedConstraints() && pos == nb);
    return stochastic_selection<int>(arity, 0, nb - 1, nb / 2);
}

unsigned int WCSP::numberOfConnectedBinaryConstraints() const
{
    unsigned int res = 0;
    for (unsigned int i = 0; i < constrs.size(); i++)
        if (constrs[i]->connected() && constrs[i]->arity() == 2 && !constrs[i]->isSep())
            res++;
    for (int i = 0; i < elimBinOrder; i++)
        if (elimBinConstrs[i]->connected() && !elimBinConstrs[i]->isSep())
            res++;
    return res;
}

unsigned int WCSP::medianDomainSize() const
{
    unsigned int nbunvars = numberOfUnassignedVariables();
    if (nbunvars == 0)
        return 0;
    unsigned int domain[nbunvars];
    unsigned int pos = 0;
    for (unsigned int i = 0; i < vars.size(); i++)
        if (unassigned(i)) {
            domain[pos] = getDomainSize(i);
            pos++;
        }
    assert(pos == numberOfUnassignedVariables() && pos == nbunvars);
    return stochastic_selection<unsigned int>(domain, 0, nbunvars - 1, nbunvars / 2);
}

unsigned int WCSP::medianDegree() const
{
    unsigned int nbunvars = numberOfUnassignedVariables();
    if (nbunvars == 0)
        return 0;
    int degree[nbunvars];
    unsigned int pos = 0;
    for (unsigned int i = 0; i < vars.size(); i++)
        if (unassigned(i)) {
            degree[pos] = getTrueDegree(i);
            pos++;
        }
    assert(pos == numberOfUnassignedVariables() && pos == nbunvars);
    return stochastic_selection<int>(degree, 0, nbunvars - 1, nbunvars / 2);
}

void WCSP::printNCBuckets()
{
    int lastbucket = -1;
    for (int bucket = 0; bucket < NCBucketSize; bucket++) {
        if (NCBuckets[bucket].begin() != NCBuckets[bucket].end())
            lastbucket = bucket;
    }

    for (int bucket = 0; bucket <= lastbucket; bucket++) {
        cout << "NC " << bucket << ":";
        for (VariableList::iterator iter = NCBuckets[bucket].begin(); iter != NCBuckets[bucket].end(); ++iter) {
            cout << " " << (*iter)->getName() << "," << (*iter)->getMaxCostValue() << "," << (*iter)->getMaxCost();

            assert((*iter)->canbe((*iter)->getMaxCostValue()));
            assert((*iter)->getCost((*iter)->getMaxCostValue()) == (*iter)->getMaxCost() || !LUBTEST((*iter)->getMaxCost(), (*iter)->getCost((*iter)->getMaxCostValue())));
            assert((bucket && !PARTIALORDER) ? (to_double((*iter)->getMaxCost()) >= (Long)powl(2., bucket)) : ((*iter)->getMaxCost() > MIN_COST));
            assert(PARTIALORDER || to_double((*iter)->getMaxCost()) < (Long)powl(2., bucket + 1));
        }
        cout << endl;
    }
}

/** \defgroup verbosity Output messages, verbosity options and debugging
 *
 * Depending on verbosity level given as option "-v=level", \p toulbar2 will output:
 * - (level=0, no verbosity) default output mode: shows version number, number of variables and cost functions read in the problem file,
 *  number of unassigned variables and cost functions after preprocessing,
 *  problem upper and lower bounds after preprocessing.
 *  Outputs current best solution cost found,
 *  ends by giving the optimum or "No solution".
 *  Last output line should always be: "end."
 * - (level=-1, no verbosity) restricted output mode: do not print current best solution cost found
 * -# (level=1) shows also search choices ("["\e search_depth \e problem_lower_bound \e problem_upper_bound \e sum_of_current_domain_sizes"] Try" \e variable_index \e operator \e value)
 *  with \e operator being assignment ("=="), value removal ("!="), domain splitting ("<=" or ">=", also showing EAC value in parenthesis)
 * -# (level=2) shows also current domains (\e variable_index \e list_of_current_domain_values "/" \e number_of_cost_functions (see approximate degree in \ref varelim) "/" \e weighted_degree \e list_of_unary_costs "s:" \e support_value) before each search choice
 *  and reports problem lower bound increases, NC bucket sort data (see \ref ncbucket), and basic operations on domains of variables
 * -# (level=3) reports also basic arc EPT operations on cost functions (see \ref softac)
 * -# (level=4) shows also current list of cost functions for each variable and reports more details on arc EPT operations (showing all changes in cost functions)
 * -# (level=5) reports more details on cost functions defined in extension giving their content (cost table by first increasing values in the current domain of the last variable in the scope)
 *
 * For debugging purposes, another option "-Z=level" allows one to monitor the search:
 * -# (level 1) shows current search depth (number of search choices from the root of the search tree) and reports statistics on nogoods for BTD-like methods
 * -# (level 2) idem
 * -# (level 3) also saves current problem into a file before each search choice
 *
 * \note \p toulbar2, compiled in debug mode, can be more verbose and it checks a lot of assertions (pre/post conditions in the code)
 *
 * \note \p toulbar2 will output an help message giving available options if run without any parameters
 *
 */

void WCSP::print(ostream& os)
{
    //    os << "Objective: [" << std::fixed << std::setprecision(ToulBar2::decimalPoint) << getDLb() << "," << getDUb() << "]" << std::setprecision(DECIMAL_POINT) << endl;
    os << "Objective: [" << getLb() << "," << getUb() << "]" << endl;
    os << "Variables:" << endl;
    for (unsigned int i = 0; i < vars.size(); i++)
        os << *vars[i] << endl;
    if (ToulBar2::verbose >= 4) {
        os << "Constraints:" << endl;
        for (unsigned int i = 0; i < constrs.size(); i++)
            if (constrs[i]->connected())
                os << *constrs[i];
        for (int i = 0; i < elimBinOrder; i++)
            if (elimBinConstrs[i]->connected())
                os << *elimBinConstrs[i];
        for (int i = 0; i < elimTernOrder; i++)
            if (elimTernConstrs[i]->connected())
                os << *elimTernConstrs[i];
    }
}

void printClique(ostream& os, int arity, Constraint* ctr)
{
    assert(arity >= 2);
    if (arity > MAX_ARITY / 10) {
        cerr << "warning! cost function arity is too large for primal graph representation." << endl;
        return;
    }
    for (int i = 0; i < arity - 1; i++) {
        for (int j = i + 1; j < arity; j++) {
//            os << "v" << ctr->getVar(i)->wcspIndex + 1 << " -- v" << ctr->getVar(j)->wcspIndex + 1 << " [len= " << ctr->getTightness() << "];" << endl;
            if (ctr->getVar(i)->wcspIndex < ctr->getVar(j)->wcspIndex) {
                os << ctr->getVar(i)->getName() << " -- " << ctr->getVar(j)->getName() << " [len= " << ctr->getTightness() << "];" << endl;
            } else {
                os << ctr->getVar(j)->getName() << " -- " << ctr->getVar(i)->getName() << " [len= " << ctr->getTightness() << "];" << endl;
            }
        }
    }
}

// Warning! make the assumption that all initial domains start at zero!!!
void WCSP::dump(ostream& os, bool original)
{
    Value maxdomsize = 0;
    unsigned int maxdomsizeUI = 0;
    Value xcosts = 0;
    // dump filename
    char Pb_basename[512];
    char Pb_graph[512];
    char Pb_degree[512];

    strcpy(Pb_basename, ToulBar2::problemsaved_filename.c_str());
    strcpy(Pb_graph, Pb_basename);
    strcpy(Pb_degree, Pb_basename);

    if (getLb() > MIN_COST)
        xcosts++;
    for (unsigned int i = 0; i < vars.size(); i++) {
        if (original && vars[i]->getInf() < 0) {
            cerr << "Cannot save domain of variable " << vars[i]->getName() << " with negative values!!!" << endl;
            exit(EXIT_FAILURE);
        }
        if (original) {
            int domsize = (vars[i]->enumerated() ? ((EnumeratedVariable*)vars[i])->getDomainInitSize() : (vars[i]->getSup() + 1));
            if (domsize > maxdomsize)
                maxdomsize = domsize;
        } else {
            if (vars[i]->unassigned() && vars[i]->getDomainSize() > maxdomsizeUI)
                maxdomsizeUI
                    = vars[i]->getDomainSize();
        }
        if (vars[i]->enumerated() && (original || vars[i]->unassigned()))
            xcosts++;
        //          else if (vars[i]->getInfCost() > MIN_COST || vars[i]->getSupCost() > MIN_COST) {
        //              cerr << "Cannot save interval variable " << vars[i]->getName() << " with bound unary costs!!!" << endl;
        //              exit(EXIT_FAILURE);
        //          }
    }
    os << "wcsp " << ((original) ? numberOfVariables() : numberOfUnassignedVariables()) << " "
       << ((original) ? maxdomsize : maxdomsizeUI) << " " << numberOfConnectedConstraints() + xcosts << " "
       << getUb() << endl;
    unsigned int nbvar = 0;
    for (unsigned int i = 0; i < vars.size(); i++) {
        if (original) {
            if (!vars[i]->enumerated())
                os << "-";
            int domsize = (vars[i]->enumerated() ? ((EnumeratedVariable*)vars[i])->getDomainInitSize() : (vars[i]->getSup() + 1));
            os << domsize;
            if (i < vars.size() - 1)
                os << " ";
        } else if (vars[i]->unassigned()) {
            nbvar++;
            if (!vars[i]->enumerated())
                os << "-";
            os << vars[i]->getDomainSize();
            if (nbvar < numberOfUnassignedVariables())
                os << " ";
        }
    }
    if (((original) ? numberOfVariables() : numberOfUnassignedVariables()) > 0)
        os << endl;
    for (unsigned int i = 0; i < constrs.size(); i++)
        if (constrs[i]->connected() && !constrs[i]->isSep())
            constrs[i]->dump(os, original);
    for (int i = 0; i < elimBinOrder; i++)
        if (elimBinConstrs[i]->connected() && !elimBinConstrs[i]->isSep())
            elimBinConstrs[i]->dump(os, original);
    for (int i = 0; i < elimTernOrder; i++)
        if (elimTernConstrs[i]->connected() && !elimTernConstrs[i]->isSep())
            elimTernConstrs[i]->dump(os, original);
    for (unsigned int i = 0; i < vars.size(); i++) {
        if (vars[i]->enumerated() && (original || vars[i]->unassigned())) {
            int size = vars[i]->getDomainSize();
            ValueCost domcost[size]; // replace size by MAX_DOMAIN_SIZE in case of compilation problem
            getEnumDomainAndCost(i, domcost);
            os << "1 " << ((original) ? i : vars[i]->getCurrentVarId()) << " " << getUb() << " " << size << endl;
            for (int v = 0; v < size; v++) {
                os << ((original) ? (domcost[v].value) : v) << " "
                   << ((original) ? domcost[v].cost : min(getUb(), domcost[v].cost)) << endl;
            }
        }
    }
    if (getLb() > MIN_COST)
        os << "0 " << getLb() << " 0" << endl;

    if (!ToulBar2::uaieval && ToulBar2::verbose >= 0) {
        //####################" dump dot file ###############################""
        strcat(Pb_graph, ".dot");
        cout << " Graph structure saved in problem.dot " << endl;
        ofstream pb(Pb_graph);
        pb << " graph graphname {\n " << endl;
        int res = 0;
        for (unsigned int i = 0; i < constrs.size(); i++)
            if (constrs[i]->connected())
                res += (constrs[i]->arity() * (constrs[i]->arity() - 1) / 2);
        for (int i = 0; i < elimBinOrder; i++)
            if (elimBinConstrs[i]->connected())
                res += (elimBinConstrs[i]->arity() * (elimBinConstrs[i]->arity() - 1) / 2);
        for (int i = 0; i < elimTernOrder; i++)
            if (elimTernConstrs[i]->connected())
                res += (elimTernConstrs[i]->arity() * (elimTernConstrs[i]->arity() - 1)
                    / 2);
        pb << "// number of constraint = " << res << " number of variable=  " << numberOfVariables() << endl;
        for (unsigned int i = 0; i < constrs.size(); i++)
            if (!constrs[i]->isSep() && constrs[i]->connected()) {
                //            pb << constrs[i]->getVar(0)->wcspIndex + 1;
                //            for (int j=1; j<constrs[i]->arity(); j++) {
                //                pb << " " << constrs[i]->getVar(j)->wcspIndex + 1;
                //            }
                //            pb << endl;
                printClique(pb, constrs[i]->arity(), constrs[i]);
            }
        for (int i = 0; i < elimBinOrder; i++)
            if (elimBinConstrs[i]->connected()) {
                //            pb << elimBinConstrs[i]->getVar(0)->wcspIndex + 1;
                //            for (int j=1; j<elimBinConstrs[i]->arity(); j++) {
                //                pb << " " << elimBinConstrs[i]->getVar(j)->wcspIndex + 1;
                //            }
                //            pb << endl;
                printClique(pb, elimBinConstrs[i]->arity(), elimBinConstrs[i]);
            }
        for (int i = 0; i < elimTernOrder; i++)
            if (elimTernConstrs[i]->connected()) {
                //            pb << elimTernConstrs[i]->getVar(0)->wcspIndex + 1;
                //            for (int j=1; j<elimTernConstrs[i]->arity(); j++) {
                //                pb << " " << elimTernConstrs[i]->getVar(j)->wcspIndex + 1;
                //            }
                //            pb << endl;
                printClique(pb, elimTernConstrs[i]->arity(), elimTernConstrs[i]);
            }
        pb << "}" << endl;

//####################" end dump dot file ###############################""
//#######################dump degree distribution ###################
#ifdef BOOST
        cout << "Connected components: " << connectedComponents() << endl;
        cout << "Biconnected components: " << biConnectedComponents() << endl;
        cout << "Diameter : " << diameter() << endl;
#endif

        int* degDistrib = new int[vars.size()];

        for (unsigned int i = 0; i < vars.size(); i++)
            degDistrib[i] = 0;
        for (unsigned int i = 0; i < vars.size(); i++)
            if (unassigned(i))
                degDistrib[getTrueDegree(i)]++;

        unsigned int lastnonzero = 0;
        for (unsigned int i = 0; i < vars.size(); i++)
            if (degDistrib[i])
                lastnonzero = i;

        strcat(Pb_degree, ".degree"); // after preprocessing
        ofstream file(Pb_degree);
        for (unsigned int i = 0; i <= lastnonzero; i++)
            if (degDistrib[i])
                file << i << " " << degDistrib[i] << endl;
        delete[] degDistrib;

        //#######################dump degree distribution ###################
    }

    if (ToulBar2::pedigree) {
        string problemname = ToulBar2::problemsaved_filename;
        if (problemname.rfind(".wcsp") != string::npos)
            problemname.replace(problemname.rfind(".wcsp"), 5, ".pre");
        ToulBar2::pedigree->save((problemname.rfind("problem.pre") == string::npos) ? problemname.c_str() : "problem_corrected.pre", this, false, true);
    }
}

void WCSP::dump_CFN(ostream& os, bool original)
{
    bool printed = false;
    std::ios_base::fmtflags f(os.flags());
    // dump filename in ToulBar2::problemsaved_filename

    for (unsigned int i = 0; i < vars.size(); i++) {
        if (vars[i]->getInf() < 0 || !vars[i]->enumerated()) {
            cerr << "Cannot save domain of variable " << vars[i]->getName() << " (negative values or not enumerated)" << endl;
            exit(EXIT_FAILURE);
        }
    }
    // Header
    os << "{\"problem\":{\"name\":\"" << getName() << "\",\"mustbe\":\"" << ((ToulBar2::costMultiplier < 0) ? ">" : "<");
    os << fixed << setprecision(ToulBar2::decimalPoint);
    os << getDPrimalBound() << "\"},\n";

    // Domain variables
    os << "\"variables\":{\n";
    for (unsigned int i = 0; i < vars.size(); i++) {
        assert(enumerated(i));
        EnumeratedVariable* s = static_cast<EnumeratedVariable*>(vars[i]);
        if (original) {
            os << "\"" << s->getName() << "\":";
            os << "[";
            printed = false;
            for (size_t p = 0; p < s->getDomainInitSize(); p++) {
                if (printed)
                    os << ",";
                os << "\"" << ((s->isValueNames()) ? s->getValueName(p) : ("v" + std::to_string(s->toValue(p)))) << "\"";
                printed = true;
            }
        } else if (s->unassigned()) {
            os << "\"" << s->getName() << "\":";
            int domsize = s->getDomainSize();
            Value* values = new Value[domsize];
            s->getDomain(values);
            os << "[";
            printed = false;
            for (int p = 0; p < domsize; p++) {
                if (printed)
                    os << ",";
                os << "\"" << ((s->isValueNames()) ? s->getValueName(s->toIndex(values[p])) : ("v" + std::to_string(values[p]))) << "\"";
                printed = true;
            }
        }
        os << "]";
        if (i < vars.size() - 1)
            os << ",";
        os << "\n";
    }

    os << "},\n\"functions\": {\n";
    for (unsigned int i = 0; i < constrs.size(); i++)
        if (constrs[i]->connected() && !constrs[i]->isSep())
            constrs[i]->dump_CFN(os, original);
    for (int i = 0; i < elimBinOrder; i++)
        if (elimBinConstrs[i]->connected() && !elimBinConstrs[i]->isSep())
            elimBinConstrs[i]->dump_CFN(os, original);
    for (int i = 0; i < elimTernOrder; i++)
        if (elimTernConstrs[i]->connected() && !elimTernConstrs[i]->isSep())
            elimTernConstrs[i]->dump_CFN(os, original);
    for (unsigned int i = 0; i < vars.size(); i++) {
        if (vars[i]->enumerated() && (original || vars[i]->unassigned())) {
            int size = vars[i]->getDomainSize();
            ValueCost domcost[size]; // replace size by MAX_DOMAIN_SIZE in case of compilation problem
            getEnumDomainAndCost(i, domcost);
            os << "\"F_" << ((original) ? i : vars[i]->getCurrentVarId()) << "\":{\"scope\":[";
            os << vars[i]->getName() << "],\"defaultcost\":" << getDPrimalBound()-negCost << ",\n";
            os << "\"costs\":[";
            for (int v = 0; v < size; v++) {
                os << ((original) ? (((EnumeratedVariable *) vars[i])->toIndex(domcost[v].value)) : v) << ","
                   << ((original) ? Cost2RDCost(domcost[v].cost) : min(getDPrimalBound()-negCost, Cost2RDCost(domcost[v].cost)));
                if (v != (size - 1)) {
                    os << ",";
                }
            }
            os << "]},\n";
        }
    }
    os << "\"F\":{\"scope\":[],\"costs\":[" << getDDualBound() << "]}\n}\n}" << endl;
    os.flags(f);
}
ostream& operator<<(ostream& os, WCSP& wcsp)
{
    wcsp.print(os);
    return os;
}

ostream& operator<<(ostream& os, WeightedCSP& wcsp)
{
    wcsp.print(os);
    return os;
}

/*
 * WCSP propagation methods
 *
 */

bool WCSP::verify()
{
    for (unsigned int i = 0; i < vars.size(); i++) {
        if (vars[i]->unassigned()) {
            if (td) {
                if (td->isActiveAndInCurrentClusterSubTree(vars[i]->getCluster())) {
                    if (!vars[i]->verifyNC())
                        return false;
#ifdef DEECOMPLETE
                    if (ToulBar2::DEE_ && !vars[i]->verifyDEE())
                        return false;
#endif
                }
            } else {
                if (!vars[i]->verifyNC())
                    return false;
#ifdef DEECOMPLETE
                if (ToulBar2::DEE_ && !vars[i]->verifyDEE())
                    return false;
#endif
            }
        }
        // Warning! in the CSP case, EDAC is no equivalent to GAC on ternary constraints due to the combination with binary constraints
        // Warning bis! isEAC() may change the current support for variables and constraints during verify (when supports are not valid due to VAC epsilon heuristic for instance)
        bool old_fulleac = false;
        if (vars[i]->enumerated())
            old_fulleac = vars[i]->isFullEAC();
        if (ToulBar2::LcLevel == LC_EDAC && vars[i]->enumerated() && vars[i]->unassigned() && !CSP(getLb(), getUb()) && !((EnumeratedVariable*)vars[i])->isEAC(vars[i]->getSupport())) {
            if (ToulBar2::verbose >= 4)
                cout << endl
                     << *this;
            cout << "warning! support of variable " << vars[i]->getName() << " not EAC!" << endl;
            if (!ToulBar2::vacValueHeuristic)
                return false;
        }
        if (ToulBar2::FullEAC && vars[i]->unassigned() && old_fulleac && old_fulleac != vars[i]->isFullEAC()) {
            if (ToulBar2::verbose >= 4)
                cout << endl
                     << *this;
            if (Store::getDepth() >= 1 || ToulBar2::setvalue != NULL) { // do not report error before preprocessing is done
                cout << endl
                     << "check:" << ((EnumeratedVariable*)vars[i])->checkEACGreedySolution() << endl;
                cout << "warning! support " << vars[i]->getSupport() << " of variable " << vars[i]->getName() << " has wrong FullEAC status!" << endl;
            }
            if (Store::getDepth() >= max(1, abs(ToulBar2::vac)))
                return false;
        }
    }
    if (ToulBar2::LcLevel >= LC_AC) {
        for (unsigned int i = 0; i < constrs.size(); i++) {
            if (constrs[i]->connected() && !constrs[i]->verify())
                return false;
        }
        for (int i = 0; i < elimBinOrder; i++) {
            if (elimBinConstrs[i]->connected() && !elimBinConstrs[i]->verify())
                return false;
        }
        for (int i = 0; i < elimTernOrder; i++) {
            if (elimTernConstrs[i]->connected() && !elimTernConstrs[i]->verify())
                return false;
        }
    }
    return true;
}

void WCSP::whenContradiction()
{
    NC.clear();
    IncDec.clear();
    AC.clear();
    DAC.clear();
    EAC1.clear();
    EAC2.clear();
    Eliminate.clear();
    DEE.clear();
    FEAC.clear();
    objectiveChanged = false;
    nbNodes++;
}

///\defgroup ncbucket NC bucket sort
/// maintains a sorted list of variables having non-zero unary costs in order to make NC propagation incremental.\n
/// - variables are sorted into buckets
/// - each bucket is associated to a single interval of non-zero costs (using a power-of-two scaling, first bucket interval is [1,2[, second interval is [2,4[, etc.)
/// - each variable is inserted into the bucket corresponding to its largest unary cost in its domain
/// - variables having all unary costs equal to zero do not belong to any bucket
///
/// NC propagation will revise only variables in the buckets associated to costs sufficiently large wrt current objective bounds.

void WCSP::propagateNC()
{
    if (ToulBar2::verbose >= 2)
        cout << "NCQueue size: " << NC.getSize() << " (" << NCBucketSize << " buckets maxi)" << endl;
    while (!NC.empty()) {
        Variable* x = NC.pop();
        if (x->unassigned())
            x->propagateNC();
    }
    if (ToulBar2::verbose >= 3) {
        for (unsigned int i = 0; i < vars.size(); i++)
            cout << *vars[i] << endl;
    }
    if (ToulBar2::verbose >= 2)
        printNCBuckets();

    if (objectiveChanged) {
        objectiveChanged = false;
        int bucket = min(cost2log2glb(getUb() - (getLb() + rounding(UNIT_COST) - UNIT_COST)), NCBucketSize - 1);
        if (bucket < 0)
            bucket = 0;
        for (; bucket < NCBucketSize; bucket++) {
            for (VariableList::iterator iter = NCBuckets[bucket].begin(); iter != NCBuckets[bucket].end();) {
                Variable* x = *iter;
                ++iter; // Warning! the iterator could be moved to another place by propagateNC
                if (x->unassigned() && CUT(x->getMaxCost() + getLb(), getUb())) {
                    if (td) {
                        if (td->isActiveAndInCurrentClusterSubTree(x->getCluster()))
                            x->propagateNC();
                    } else
                        x->propagateNC();
                }
            }
        }
    }
    if (objectiveChanged || !NC.empty())
        propagateNC();
}

void WCSP::propagateIncDec()
{
    if (ToulBar2::verbose >= 2)
        cout << "IncDecQueue size: " << IncDec.getSize() << endl;
    while (!IncDec.empty()) {
        int incdec;
        Variable* x = IncDec.pop(&incdec);
        if (x->unassigned())
            x->propagateIncDec(incdec);
    }
}

void WCSP::propagateAC()
{
    if (ToulBar2::verbose >= 2)
        cout << "ACQueue size: " << AC.getSize() << endl;
    if (Store::getDepth() == 0)
        AC.sort(false);
    while (!AC.empty()) {
        EnumeratedVariable* x = (EnumeratedVariable*)((ToulBar2::QueueComplexity) ? AC.pop_min() : AC.pop());
        if (x->unassigned())
            x->propagateAC();
        // Warning! propagateIncDec() necessary to transform inc/dec event into remove event
        propagateIncDec(); // always examine inc/dec events before remove events
    }
}

void WCSP::propagateDAC()
{
    if (ToulBar2::verbose >= 2)
        cout << "DACQueue size: " << DAC.getSize() << endl;
    if (Store::getDepth() == 0)
        DAC.sort(true);
    while (!DAC.empty()) {
        if (ToulBar2::interrupted)
            throw TimeOut();
        EnumeratedVariable* x = (EnumeratedVariable*)((ToulBar2::QueueComplexity) ? DAC.pop_max() : DAC.pop());
        if (x->unassigned())
            x->propagateDAC();
        propagateIncDec(); // always examine inc/dec events before projectFromZero events
    }
}

void WCSP::propagateTRWS()
{
    bool forwardPass = true;
    unsigned int nIteration = 0;
    unsigned int nIterationNoChange = 0;
    Cost previousEbound = MIN_COST;
    Cost ebound = MIN_COST;
    Cost bestUb = getUb();
    vector<int> orders[2] = { vector<int>(numberOfVariables()), vector<int>(numberOfVariables()) };
    vector<unsigned int> ranks[2] = { vector<unsigned int>(numberOfVariables()), vector<unsigned int>(numberOfVariables()) };
    vector<Cost> tmpM(getMaxDomainSize(), numeric_limits<Cost>::max());
    vector<Value> bestPrimalVal(numberOfVariables(), 0);
    vector<int> bestPrimalVar(numberOfVariables(), 0);
    for (unsigned int i = 0; i < numberOfVariables(); ++i) {
        bestPrimalVar[i] = i;
        bestPrimalVal[i] = getSupport(i);
    }

    assert(!td); // warning! tree decomposition must be done after TRW-S

    // Preprocessing: compute order
    if (ToulBar2::trwsOrder) {
        // Use monotonic chains
        unsigned int nVariableUsed = 0;
        vector<bool> variableUsed(numberOfVariables(), false);
        unsigned int firstVariableId = 0;
        while (nVariableUsed != numberOfVariables()) {
            while (variableUsed[getVar(firstVariableId)->wcspIndex])
                ++firstVariableId;
            Variable* firstVariable = getVar(firstVariableId);
            variableUsed[firstVariable->wcspIndex] = true;
            ranks[0][firstVariable->wcspIndex] = nVariableUsed;
            orders[0][nVariableUsed++] = firstVariable->wcspIndex;
            if (assigned(firstVariableId) || !enumerated(firstVariableId))
                continue;
            bool stillPath;
            do {
                stillPath = false;
                for (ConstraintList::iterator iter = firstVariable->getConstrs()->begin(); (iter != firstVariable->getConstrs()->end()) && (!stillPath); ++iter) {
                    Constraint* constraint = (*iter).constr;
                    if (constraint->isBinary()) {
                        BinaryConstraint* binaryConstraint = static_cast<BinaryConstraint*>(constraint);
                        EnumeratedVariable* otherVariable = static_cast<EnumeratedVariable*>(binaryConstraint->getVarDiffFrom(firstVariable));
                        if (!variableUsed[otherVariable->wcspIndex]) {
                            firstVariable = otherVariable;
                            variableUsed[firstVariable->wcspIndex] = true;
                            ranks[0][firstVariable->wcspIndex] = nVariableUsed;
                            orders[0][nVariableUsed++] = firstVariable->wcspIndex;
                            stillPath = true;
                        }
                    }
                }
            } while (stillPath);
        }
    } else {
        // Use DAC order
        for (unsigned int i = 0; i < numberOfVariables(); ++i) {
            Variable* var = getVar(i);
            orders[0][var->getDACOrder()] = var->wcspIndex;
            ranks[0][var->wcspIndex] = var->getDACOrder();
        }
    }

    orders[1] = orders[0];
    reverse(orders[1].begin(), orders[1].end());
    for (unsigned int i = 0; i < numberOfVariables(); ++i) {
        ranks[1][i] = numberOfVariables() - ranks[0][i] - 1;
    }

    // DAC compatible order is the opposite
    if (*(orders[0].begin()) < *(orders[0].rbegin())) {
        swap(orders[0], orders[1]);
        swap(ranks[0], ranks[1]);
    }

    if (ToulBar2::trwsOrder) {
        setDACOrder(orders[0]);
    }

    // Preprocessing: compute gammas
    for (unsigned int i = 0; i < numberOfVariables(); ++i)
        if (unassigned(i) && enumerated(i)) {
            EnumeratedVariable* s = static_cast<EnumeratedVariable*>(getVar(i));
            unsigned int r = ranks[0][s->wcspIndex];
            int nCtrIn = 0, nCtrOut = 0;
            for (ConstraintList::iterator iter = s->getConstrs()->begin(); iter != s->getConstrs()->end(); ++iter) {
                Constraint* constraint = (*iter).constr;
                if (constraint->isBinary()) {
                    BinaryConstraint* binctr = static_cast<BinaryConstraint*>(constraint);
                    if (ranks[0][binctr->getVarDiffFrom(s)->wcspIndex] < r) {
                        ++nCtrIn;
                    } else {
                        ++nCtrOut;
                    }
                }
            }
            s->setTRWSGamma(1.0 / (max<unsigned int>(nCtrIn, nCtrOut) + 0.1));
        }
    // Preprocessing: reset trwsM
    for (unsigned int i = 0; i < numberOfConstraints(); i++) {
        Constraint* ctr = getCtr(i);
        if (ctr->connected() && ctr->isBinary()) {
            BinaryConstraint* binctr = static_cast<BinaryConstraint*>(ctr);
            binctr->trwsM.clear();
        }
    }
    for (int i = 0; i < elimBinOrder; i++) {
        Constraint* ctr = elimBinConstrs[i];
        if (ctr->connected()) {
            BinaryConstraint* binctr = static_cast<BinaryConstraint*>(ctr);
            binctr->trwsM.clear();
        }
    }

    do {
        vector<int>& order = (forwardPass) ? orders[0] : orders[1];
        vector<unsigned int>& rank = (forwardPass) ? ranks[0] : ranks[1];
        ebound = MIN_COST;
        for (unsigned int i = 0; i < numberOfVariables(); ++i) {
            if (unassigned(order[i]) && enumerated(order[i])) {
                if (ToulBar2::interrupted)
                    throw TimeOut();
                EnumeratedVariable* s = static_cast<EnumeratedVariable*>(getVar(order[i]));
                // step 1: normalize unary costs
                vector<Cost> thetaHat(s->getDomainInitSize(), MIN_COST);
                Cost delta = numeric_limits<Cost>::max();
                for (EnumeratedVariable::iterator sIter = s->begin(); sIter != s->end(); ++sIter) {
                    unsigned int j = s->toIndex(*sIter);
                    thetaHat[j] = s->getCost(*sIter);
                    for (ConstraintList::iterator iter = s->getConstrs()->begin(); iter != s->getConstrs()->end(); ++iter) {
                        Constraint* constraint = (*iter).constr;
                        if (constraint->isBinary()) {
                            BinaryConstraint* binctr = static_cast<BinaryConstraint*>(constraint);
                            thetaHat[j] += binctr->trwsM[j];
                        }
                    }
                    delta = min<Cost>(delta, thetaHat[j]);
                }
                if (delta != MIN_COST) {
                    for (EnumeratedVariable::iterator sIter = s->begin(); sIter != s->end(); ++sIter) {
                        unsigned int j = s->toIndex(*sIter);
                        thetaHat[j] -= delta;
                    }
                    ebound += delta;
                }
                // step 2: message update
                for (ConstraintList::iterator iter = s->getConstrs()->begin(); iter != s->getConstrs()->end(); ++iter) {
                    Constraint* constraint = (*iter).constr;
                    if (constraint->isBinary()) {
                        BinaryConstraint* binctr = static_cast<BinaryConstraint*>(constraint);
                        EnumeratedVariable* t = static_cast<EnumeratedVariable*>(binctr->getVarDiffFrom(s));
                        if (rank[s->wcspIndex] < rank[t->wcspIndex]) {
                            std::function<Cost(unsigned int, unsigned int)> getCost1 = [binctr](unsigned int x, unsigned int y) { return binctr->getCost(x, y); };
                            std::function<Cost(unsigned int, unsigned int)> getCost2 = [binctr](unsigned int y, unsigned int x) { return binctr->getCost(x, y); };
                            auto getCost = (s == binctr->getVar(0)) ? getCost1 : getCost2;
                            delta = numeric_limits<Cost>::max();
                            tmpM.resize(max(s->getDomainInitSize(), t->getDomainInitSize()), numeric_limits<Cost>::max());
                            for (EnumeratedVariable::iterator tIter = t->begin(); tIter != t->end(); ++tIter) {
                                unsigned int k = t->toIndex(*tIter);
                                tmpM[k] = numeric_limits<Cost>::max();
                                for (EnumeratedVariable::iterator sIter = s->begin(); sIter != s->end(); ++sIter) {
                                    unsigned int j = s->toIndex(*sIter);
                                    tmpM[k] = min<Cost>(tmpM[k], static_cast<Cost>(trunc(s->getTRWSGamma() * thetaHat[j])) - binctr->trwsM[j] + getCost(j, k));
                                }
                                delta = min<Cost>(delta, tmpM[k]);
                            }
                            binctr->trwsM = tmpM;
                            if (delta != MIN_COST) {
                                for (EnumeratedVariable::iterator tIter = t->begin(); tIter != t->end(); ++tIter) {
                                    unsigned int k = t->toIndex(*tIter);
                                    binctr->trwsM[k] -= delta;
                                }
                                ebound += delta;
                            }
                        }
                    }
                }
            }
        }
        // step 3: compute ub
        if ((!forwardPass) && (ToulBar2::trwsNIterComputeUb > 0) && (nIteration > 0) && (nIteration % ToulBar2::trwsNIterComputeUb == 0)) {
            for (unsigned int i = 0; i < numberOfVariables(); ++i) {
                if (unassigned(order[i]) && enumerated(order[i])) {
                    EnumeratedVariable* s = static_cast<EnumeratedVariable*>(getVar(order[i]));
                    Cost bestCost = numeric_limits<Cost>::max();
                    for (EnumeratedVariable::iterator sIter = s->begin(); sIter != s->end(); ++sIter) {
                        Cost cost = s->getCost(*sIter);
                        unsigned int j = s->toIndex(*sIter);
                        for (ConstraintList::iterator iter = s->getConstrs()->begin(); iter != s->getConstrs()->end(); ++iter) {
                            Constraint* constraint = (*iter).constr;
                            if (constraint->isBinary()) {
                                BinaryConstraint* binctr = static_cast<BinaryConstraint*>(constraint);
                                EnumeratedVariable* t = static_cast<EnumeratedVariable*>(binctr->getVarDiffFrom(s));
                                if (rank[s->wcspIndex] < rank[t->wcspIndex]) {
                                    cost += binctr->trwsM[j];
                                } else {
                                    cost += binctr->getCost(s, t, *sIter, bestPrimalVal[t->wcspIndex]);
                                }
                            }
                        }
                        if (cost < bestCost || (cost == bestCost && s->getSupport() == *sIter)) {
                            bestCost = cost;
                            bestPrimalVal[s->wcspIndex] = *sIter;
                        }
                    }
                }
            }
            int depth = Store::getDepth();
            try {
                Store::store();
                assignLS(bestPrimalVar, bestPrimalVal);
                assert(numberOfUnassignedVariables() == 0);
                ((Solver*)getSolver())->Solver::newSolution();
                bestUb = min<Cost>(getUb(), bestUb);
            } catch (const Contradiction&) {
                whenContradiction();
            }
            Store::restore(depth);
            enforceUb();
        }
        // step 4: reverse ordering
        double change = 0.0;
        if (!forwardPass) {
            ++nIteration;
            change = (ebound == previousEbound) ? 0.0 : (double)(ebound - previousEbound + 1) / (getLb() + ebound + 1);
            nIterationNoChange = (change > ToulBar2::trwsAccuracy) ? 0 : nIterationNoChange + 1;
            if (ToulBar2::verbose >= 0 && !nIterationNoChange) {
                Double Dglb = (ToulBar2::costMultiplier >= 0 ? Cost2ADCost(getLb() + ebound) : Cost2ADCost(bestUb));
                Double Dgub = (ToulBar2::costMultiplier >= 0 ? Cost2ADCost(bestUb) : Cost2ADCost(getLb() + ebound));
                if (ToulBar2::uai)
                    cout << "TRW-S dual bound: " << std::fixed << std::setprecision(ToulBar2::decimalPoint) << Cost2ADCost(getLb() + ebound) << std::setprecision(DECIMAL_POINT) << " energy: " << -(Cost2LogProb(getLb() + ebound) + ToulBar2::markov_log) << " -- primal bound: " << std::setprecision(ToulBar2::decimalPoint) << Cost2ADCost(bestUb) << std::setprecision(DECIMAL_POINT) << " energy: " << -(Cost2LogProb(bestUb) + ToulBar2::markov_log) << " (+" << (100 * change) << "%) (accuracy: " << (100.0 * (bestUb - ebound - getLb()) / (bestUb + 1)) << "%) (iter:" << nIteration << ")" << endl;
                else
                    cout << "TRW-S bounds: [" << std::fixed << std::setprecision(ToulBar2::decimalPoint) << Dglb << ", " << Dgub << std::setprecision(DECIMAL_POINT) << "] (+" << (100 * change) << "%) (accuracy: " << (100.0 * (bestUb - ebound - getLb()) / (bestUb + 1)) << "%) (iter:" << nIteration << ")" << endl;
            }
            previousEbound = ebound;
        }
        forwardPass = !forwardPass;
    } while ((nIteration < ToulBar2::trwsNIter) && (nIterationNoChange < ToulBar2::trwsNIterNoChange));

    // step 5: move to WCSP
    ToulBar2::trwsAccuracy = -1; // stop TRW-S such that VAC can be done
    if (ebound <= MIN_COST) {
        if (ToulBar2::verbose >= 1) {
            cout << "TRW-S did not improve the lower bound." << endl;
        }
    } else {
        Cost delta = MIN_COST;
        for (unsigned int i = 0; i < numberOfVariables(); ++i) {
            if (unassigned(orders[0][i]) && enumerated(orders[0][i])) {
                if (ToulBar2::interrupted)
                    throw TimeOut();
                EnumeratedVariable* s = static_cast<EnumeratedVariable*>(getVar(orders[0][i]));
                for (ConstraintList::iterator iter = s->getConstrs()->begin(); iter != s->getConstrs()->end(); ++iter) {
                    Constraint* constraint = (*iter).constr;
                    if (constraint->isBinary()) {
                        BinaryConstraint* binctr = static_cast<BinaryConstraint*>(constraint);
                        EnumeratedVariable* t = static_cast<EnumeratedVariable*>(binctr->getVarDiffFrom(s));
                        if (ranks[0][s->wcspIndex] < ranks[0][t->wcspIndex]) {
                            for (EnumeratedVariable::iterator sIter = s->begin(); sIter != s->end(); ++sIter) {
                                unsigned int j = s->toIndex(*sIter);
                                binctr->projectTRWS(s, *sIter, binctr->trwsM[j]);
                            }
                        }
                    }
                }
                Cost c = s->normalizeTRWS();
                delta += c;
                for (EnumeratedVariable::iterator sIter = s->begin(); sIter != s->end(); ++sIter) {
                    Cost availableCost = static_cast<Cost>(trunc(s->getTRWSGamma() * s->getCost(*sIter)));
                    for (ConstraintList::iterator iter = s->getConstrs()->begin(); iter != s->getConstrs()->end(); ++iter) {
                        Constraint* constraint = (*iter).constr;
                        if (constraint->isBinary()) {
                            BinaryConstraint* binctr = static_cast<BinaryConstraint*>(constraint);
                            EnumeratedVariable* t = static_cast<EnumeratedVariable*>(binctr->getVarDiffFrom(s));
                            if (ranks[0][s->wcspIndex] < ranks[0][t->wcspIndex]) {
                                binctr->extend(binctr->getIndex(s), *sIter, availableCost);
                            }
                        }
                    }
                }
                for (ConstraintList::iterator iter = s->getConstrs()->begin(); iter != s->getConstrs()->end(); ++iter) {
                    Constraint* constraint = (*iter).constr;
                    if (constraint->isBinary()) {
                        BinaryConstraint* binctr = static_cast<BinaryConstraint*>(constraint);
                        EnumeratedVariable* t = static_cast<EnumeratedVariable*>(binctr->getVarDiffFrom(s));
                        if (ranks[0][s->wcspIndex] < ranks[0][t->wcspIndex]) {
                            std::function<Cost(unsigned int, unsigned int)> getCost1 = [binctr](unsigned int x, unsigned int y) { return binctr->getCostTRWS(x, y); };
                            std::function<Cost(unsigned int, unsigned int)> getCost2 = [binctr](unsigned int y, unsigned int x) { return binctr->getCostTRWS(x, y); };
                            auto getCost = (s == binctr->getVar(0)) ? getCost1 : getCost2;
                            Cost minCost = numeric_limits<Cost>::max();
                            tmpM.resize(max(s->getDomainInitSize(), t->getDomainInitSize()), numeric_limits<Cost>::max());
                            for (EnumeratedVariable::iterator tIter = t->begin(); tIter != t->end(); ++tIter) {
                                unsigned int k = t->toIndex(*tIter);
                                tmpM[k] = numeric_limits<Cost>::max();
                                for (EnumeratedVariable::iterator sIter = s->begin(); sIter != s->end(); ++sIter) {
                                    tmpM[k] = min<Cost>(tmpM[k], getCost(*sIter, *tIter));
                                }
                                minCost = min<Cost>(minCost, tmpM[k]);
                            }
                            for (EnumeratedVariable::iterator tIter = t->begin(); tIter != t->end(); ++tIter) {
                                unsigned int k = t->toIndex(*tIter);
                                tmpM[k] -= minCost;
                            }
                            binctr->trwsM = tmpM;
                            for (EnumeratedVariable::iterator tIter = t->begin(); tIter != t->end(); ++tIter) {
                                binctr->projectTRWS(t, *tIter, binctr->trwsM[t->toValue(*tIter)]);
                            }
                            minCost = binctr->normalizeTRWS();
                            delta += minCost;
                        }
                    }
                }
            }
        }
        increaseLb(delta);
        propagate(); // propagate again without TRWS and possibly with VAC
    }
}

void WCSP::fillEAC2()
{
    assert(EAC2.empty());
    if (ToulBar2::verbose >= 2)
        cout << "EAC1Queue size: " << EAC1.getSize() << endl;
    while (!EAC1.empty()) {
        EnumeratedVariable* x = (EnumeratedVariable*)((ToulBar2::QueueComplexity) ? EAC1.pop_min() : EAC1.pop());
        if (x->unassigned())
            x->fillEAC2(true);
    }
}

void WCSP::propagateEAC()
{
    fillEAC2();
    if (ToulBar2::verbose >= 2)
        cout << "EAC2Queue size: " << EAC2.getSize() << endl;
    if (Store::getDepth() == 0)
        EAC2.sort(false);
    while (!EAC2.empty()) {
        if (ToulBar2::interrupted)
            throw TimeOut();
        EnumeratedVariable* x = (EnumeratedVariable*)((ToulBar2::QueueComplexity) ? EAC2.pop_min() : EAC2.pop());
        if (x->unassigned())
            x->propagateEAC();
        propagateIncDec(); // always examine inc/dec events before projectFromZero events
    }
}

void WCSP::propagateSeparator()
{
    if (ToulBar2::verbose >= 2)
        cout << "PendingSeparator size: " << PendingSeparator.getSize() << endl;
    for (SeparatorList::iterator iter = PendingSeparator.begin(); iter != PendingSeparator.end(); ++iter) {
        (*iter)->propagate();
    }
}

void WCSP::propagateDEE()
{
    if (ToulBar2::verbose >= 2)
        cout << "DEEQueue size: " << DEE.getSize() << endl;
    assert(NC.empty());
    while (!DEE.empty()) {
        if (ToulBar2::interrupted)
            throw TimeOut();
        EnumeratedVariable* x = (EnumeratedVariable*)DEE.pop();
        if (x->unassigned() && !((ToulBar2::divNbSol > 1) && Store::getDepth() == 0)) {
            if (ToulBar2::DEE_ >= 3 || (ToulBar2::DEE_ == 2 && Store::getDepth() == 0)) {
                for (EnumeratedVariable::iterator itera = x->begin(); itera != x->end(); ++itera) {
                    for (EnumeratedVariable::iterator iterb = x->lower_bound(*itera + 1); iterb != x->end(); ++iterb) {
                        assert(x->canbe(*itera));
                        assert(x->canbe(*iterb));
                        assert(*itera != *iterb);
                        x->propagateDEE(*itera, *iterb, false);
                        if (x->cannotbe(*itera))
                            break;
                    }
                }
            } else {
                Value a = x->getSupport();
                Value b = x->getMaxCostValue();
                assert(x->canbe(a) && x->getCost(a) == MIN_COST);
                assert(x->canbe(b) && x->getCost(b) == x->getMaxCost());
                if (a == b) {
                    assert(x->getMaxCost() == MIN_COST);
                    if (a != x->getSup())
                        b = x->getSup();
                    else
                        b = x->getInf();
                }
                assert(a != b);
                x->propagateDEE(a, b);
            }
            propagateNC(); // DEE assumes NC already done
        }
    }
}

void WCSP::propagateFEAC()
{
    if (ToulBar2::verbose >= 2)
        cout << "FEACQueue size: " << FEAC.getSize() << endl;
    while (!FEAC.empty()) {
        if (ToulBar2::interrupted)
            throw TimeOut();
        EnumeratedVariable* x = (EnumeratedVariable*)FEAC.pop();
        if (x->unassigned())
            x->reviseEACGreedySolution();
    }
}

/// \defgroup varelim Variable elimination
/// - \e i-bounded variable elimination eliminates all variables with a degree less than or equal to \e i.
///		It can be done with arbitrary i-bound in preprocessing only and iff all their cost functions are in extension.
/// - \e i-bounded variable elimination with i-bound less than or equal to two can be done during the search.
/// - functional variable elimination eliminates all variables which have a bijective or functional binary hard constraint (\e ie ensuring a one-to-one or several-to-one value mapping) and iff all their cost functions are in extension.
///		It can be done without limit on their degree, in preprocessing only.
/// \note Variable elimination order used in preprocessing is either lexicographic or given by an external file *.order (see toulbar2 options)
/// \note 2-bounded variable elimination during search is optimal in the sense that any elimination order should result in the same final graph
/// \warning It is not possible to display/save solutions when bounded variable elimination is applied in preprocessing
/// \warning toulbar2 maintains a list of current cost functions for each variable.
///		It uses the size of these lists as an approximation of variable degrees.
///		During the search, if variable \e x has three cost functions \e xy, \e xz, \e xyz, its true degree is two but its approximate degree is three.
///		In toulbar2 options, it is the approximate degree which is given by the user for variable elimination during the search (thus, a value at most three).
///		But it is the true degree which is given by the user for variable elimination in preprocessing.

void WCSP::eliminate()
{
    while (!Eliminate.empty()) {
        if (ToulBar2::interrupted)
            throw TimeOut();
        EnumeratedVariable* x = (EnumeratedVariable*)Eliminate.pop();
        if (x->unassigned()) {
            if (td) {
                if (td->isInCurrentClusterSubTree(x->getCluster()))
                    x->eliminate();
            } else
                x->eliminate();
        }
    }
}

/// \defgroup softac Soft arc consistency and problem reformulation
/// Soft arc consistency is an incremental lower bound technique for optimization problems.
/// Its goal is to move costs from high-order (typically arity two or three) cost functions towards the problem lower bound and unary cost functions.
/// This is achieved by applying iteratively local equivalence-preserving problem transformations (EPTs) until some terminating conditions are met.
/// \note \e eg an EPT can move costs between a binary cost function and a unary cost function such that the sum of the two functions remains the same for any complete assignment.
/// \see <em> Arc consistency for Soft Constraints. </em> T. Schiex. Proc. of CP'2000. Singapour, 2000.
/// \note Soft Arc Consistency in toulbar2 is limited to binary and ternary and some global cost functions (\e eg alldifferent, gcc, regular, same).
///		Other n-ary cost functions are delayed for propagation until their number of unassigned variables is three or less.
/// \see <em> Towards Efficient Consistency Enforcement for Global Constraints in Weighted Constraint Satisfaction. </em> Jimmy Ho-Man Lee, Ka Lun Leung. Proc. of IJCAI 2009, pages 559-565. Pasadena, USA, 2009.

/// \defgroup propagation Propagation loop
/// Propagates soft local consistencies and bounded variable elimination until all the propagation queues are empty or a contradiction occurs.\n
/// While (queues are not empty or current objective bounds have changed):
/// -# queue for bounded variable elimination of degree at most two (except at preprocessing)
/// -# BAC queue
/// -# EAC queue
/// -# DAC queue
/// -# AC queue
/// -# monolithic (flow-based and DAG-based) global cost function propagation (partly incremental)
/// -# NC queue
/// -# returns to #1 until all the previous queues are empty
/// -# DEE queue
/// -# returns to #1 until all the previous queues are empty
/// -# VAC propagation (not incremental)
/// -# returns to #1 until all the previous queues are empty (and problem is VAC if enable)
/// -# exploits goods in pending separators for BTD-like methods
///
/// Queues are first-in / first-out lists of variables (avoiding multiple insertions).
/// In case of a contradiction, queues are explicitly emptied by WCSP::whenContradiction

void WCSP::propagate()
{
    if (ToulBar2::interrupted)
        throw TimeOut();
    revise(NULL);
    if (ToulBar2::vac)
        vac->iniThreshold();

    for (vector<GlobalConstraint*>::iterator it = globalconstrs.begin(); it != globalconstrs.end(); it++) {
        (*(it))->init();
    }
    if (isGlobal() && ToulBar2::LcLevel >= LC_EDAC) {
        for (unsigned int i = 0; i < vars.size(); i++) {
            EnumeratedVariable* x = (EnumeratedVariable*)vars[i];
            if (x->unassigned()) {
                x->setCostProvidingPartition(); // For EAC propagation
            }
        }
    }

    do {
        do {
            do {
                do {
                    eliminate();
                    int eac_iter = 0;
                    while (objectiveChanged || !NC.empty() || !IncDec.empty() || ((ToulBar2::LcLevel == LC_AC || ToulBar2::LcLevel >= LC_FDAC) && !AC.empty())
                        || (ToulBar2::LcLevel >= LC_DAC
                               && !DAC.empty())
                        || (ToulBar2::LcLevel == LC_EDAC && !CSP(getLb(), getUb()) && !EAC1.empty())) {
                        eac_iter++;
                        propagateIncDec();
                        if (ToulBar2::LcLevel == LC_EDAC && !CSP(getLb(), getUb()))
                            propagateEAC();
                        assert(IncDec.empty());
                        if (ToulBar2::LcLevel >= LC_DAC)
                            propagateDAC();
                        assert(IncDec.empty());
                        if (ToulBar2::LcLevel == LC_AC || ToulBar2::LcLevel >= LC_FDAC)
                            propagateAC();
                        assert(IncDec.empty());

                        Cost oldLb = getLb();
                        bool cont = true;
                        while (cont) {
                            oldLb = getLb();
                            cont = false;
                            for (vector<GlobalConstraint*>::iterator it = globalconstrs.begin(); it != globalconstrs.end(); it++) {
                                if (ToulBar2::interrupted)
                                    throw TimeOut();
                                (*(it))->propagate();
                                if (ToulBar2::LcLevel == LC_SNIC)
                                    if (!IncDec.empty())
                                        cont = true; //For detecting value removal during SNIC enforcement
                                propagateIncDec();
                            }
                            if (ToulBar2::LcLevel == LC_SNIC)
                                if (!NC.empty() || objectiveChanged)
                                    cont = true; //For detecting value removal and upper bound change
                            propagateNC();
                            if (ToulBar2::LcLevel == LC_SNIC)
                                if (oldLb != getLb() || !AC.empty()) {
                                    cont = true;
                                    AC.clear(); //For detecting value removal and lower bound change
                                }
                        }
                        propagateNC();
                        if (ToulBar2::DEE_) {
                            propagateDEE(); // DEE requires NC and can break soft AC but not VAC
                        }
                        if (ToulBar2::LcLevel == LC_EDAC && eac_iter > MAX_EAC_ITER) {
                            EAC1.clear();
                            cout << "c automatically switch from EDAC to FDAC." << endl;
                            ToulBar2::LcLevel = LC_FDAC;
                            break;
                        } // avoids pathological cases with too many very slow lower bound increase by EAC
                    }
                } while (!Eliminate.empty());

                if (ToulBar2::LcLevel < LC_EDAC || CSP(getLb(), getUb()))
                    EAC1.clear();
                if (ToulBar2::vac && (ToulBar2::trwsAccuracy < 0) && !CSP(getLb(), getUb())) {
                    vac->propagate(); // VAC requires soft AC
                }
            } while (ToulBar2::vac && !CSP(getLb(), getUb()) && !vac->isVAC());
        } while (objectiveChanged || !NC.empty() || !IncDec.empty()
            || ((ToulBar2::LcLevel == LC_AC || ToulBar2::LcLevel >= LC_FDAC) && !AC.empty())
            || (ToulBar2::LcLevel >= LC_DAC && !DAC.empty())
            || (ToulBar2::LcLevel == LC_EDAC && !CSP(getLb(), getUb()) && !EAC1.empty())
            || !Eliminate.empty()
            || (ToulBar2::vac && !CSP(getLb(), getUb()) && !vac->isVAC()));
        // TO BE DONE AFTER NORMAL PROPAGATION
        if (td)
            propagateSeparator();
    } while (objectiveChanged);
    propagateFEAC();
    revise(NULL);

    for (vector<GlobalConstraint*>::iterator it = globalconstrs.begin(); it != globalconstrs.end(); it++) {
        (*(it))->end();
    }
    assert(verify());
    assert(!objectiveChanged);
    assert(NC.empty());
    assert(IncDec.empty());
    if (ToulBar2::LcLevel == LC_AC || ToulBar2::LcLevel >= LC_FDAC)
        assert(AC.empty());
    else
        AC.clear();
    if (ToulBar2::LcLevel >= LC_DAC)
        assert(DAC.empty());
    else
        DAC.clear();
    assert(EAC1.empty());
    assert(EAC2.empty());
    assert(Eliminate.empty());
    DEE.clear(); // DEE might not be empty if verify() has modified supports
    assert(FEAC.empty());
    nbNodes++;
}

void WCSP::restoreSolution(Cluster* c)
{
    static Tuple tctr;
    int elimo = getElimOrder();
    for (int i = elimo - 1; i >= 0; i--) {
        elimInfo ei = elimInfos[i];
        EnumeratedVariable* x = (EnumeratedVariable*)ei.x;
        EnumeratedVariable* y = (EnumeratedVariable*)ei.y;
        EnumeratedVariable* z = (EnumeratedVariable*)ei.z;
        assert(x);
        assert(x->assigned());
        if (c && !c->isVar(x->wcspIndex))
            continue;
        if (y && y->unassigned())
            continue;
        if (z && z->unassigned())
            continue;
        BinaryConstraint* xy = ei.xy;
        BinaryConstraint* xz = ei.xz;
        TernaryConstraint* xyz = ei.xyz;
        Constraint* ctr = ei.ctr;
        Value vy = -1;
        Value vz = -1;
        Cost cctr;
        int xctrindex = -1;
        if (y)
            vy = getValue(y->wcspIndex);
        if (z)
            vz = getValue(z->wcspIndex);
        if (ctr) {
            ctr->firstlex();
            ctr->nextlex(tctr, cctr);
            xctrindex = ctr->getIndex(x);
            assert(xctrindex >= 0);
        }

        Value minv = WRONG_VAL;
        Cost mincost = MAX_COST;
        for (unsigned int vxi = 0; vxi < x->getDomainInitSize(); vxi++) {
            Value vx = x->toValue(vxi);
            if (!x->canbeAfterElim(vx))
                continue;
            Cost cxy = MIN_COST;
            Cost cxz = MIN_COST;
            Cost cxyz = MIN_COST;
            cctr = MIN_COST;
            if (xy) {
                if (xy->getIndex(y) == 0)
                    cxy = xy->getCost(vy, vx);
                else
                    cxy = xy->getCost(vx, vy);
            }
            if (xz) {
                if (xz->getIndex(z) == 0)
                    cxz = xz->getCost(vz, vx);
                else
                    cxz = xz->getCost(vx, vz);
            }
            if (xyz)
                cxyz = xyz->getCost(x, y, z, vx, vy, vz);
            if (ctr) {
                tctr[xctrindex] = vxi;
                cctr = ctr->evalsubstr(tctr, ctr);
            }
            Cost loc = x->getCost(vx) + cxy + cxz + cxyz + cctr;
            //cout << "test " << vx << "," << x->getCost(vx) << "," << cxy << "," << cxz << "," << cxyz << endl;
            if (loc < mincost) {
                mincost = loc;
                minv = vx;
            }
        }
        //cout << i << ": elim " << x->getName() << "_" << minv << ", y= " << ((y)?y->getName():"-") << "_" << vy << ", z= " << ((z)?z->getName():"-") << "_" << vz << endl;
        assert(minv != WRONG_VAL);
        x->assignWhenEliminated(minv);
    }
}

// -----------------------------------------------------------
// Methods for Variable Elimination

// Creates n fake empty constraints and puts them in the pool 'elimBinConstrs'

void WCSP::initElimConstr()
{
    BinaryConstraint* xy = NULL;
    if (!ToulBar2::vac)
        xy = new BinaryConstraint(this);
    else
        xy = new VACBinaryConstraint(this);
    elimBinConstrs.push_back(xy);
    elimInfo ei = { NULL, NULL, NULL, NULL, NULL, NULL, NULL };
    elimInfos.push_back(ei);
}

void WCSP::initElimConstrs()
{
    for (unsigned int i = 0; i < vars.size(); i++)
        initElimConstr();

    vector<int> order;
    if (isAlreadyTreeDec(ToulBar2::varOrder))
        treeDecFile2Vector(ToulBar2::varOrder, order);
    else
        elimOrderFile2Vector(ToulBar2::varOrder, order);
    for (int i = vars.size() - 1; i >= 0; --i)
        vars[order[i]]->queueEliminate();
    elimSpace = 0;
}

// Function that adds a new binary constraint from the pool of fake constraints
BinaryConstraint* WCSP::newBinaryConstr(EnumeratedVariable* x, EnumeratedVariable* y, Constraint* from1, Constraint* from2)
{
    unsigned int newIndex = (int)elimBinOrder;
    assert(newIndex < elimBinConstrs.size());
    BinaryConstraint* ctr = (BinaryConstraint*)elimBinConstrs[newIndex];
    ctr->fillElimConstr(x, y, from1, from2);
    if (ToulBar2::vac)
        ((VACBinaryConstraint*)ctr)->VACfillElimConstr();
    ctr->isDuplicate_ = false;
    return ctr;
}

// warning! Do not propagate this new binary cost function
BinaryConstraint* WCSP::newBinaryConstr(EnumeratedVariable* x, EnumeratedVariable* y, vector<Cost>& costs)
{
    if (!ToulBar2::vac) {
        return new BinaryConstraint(this, x, y, costs);
    } else {
        return new VACBinaryConstraint(this, x, y, costs);
    }
}

// warning! you must create beforehand three binary constraints in fake pool (elimBinConstrs)
// if they do not exist in the main pool (constrs)
TernaryConstraint* WCSP::newTernaryConstr(EnumeratedVariable* x, EnumeratedVariable* y, EnumeratedVariable* z, Constraint* from1)
{
    unsigned int newIndex = (int)elimTernOrder;
    assert(newIndex < elimTernConstrs.size());
    TernaryConstraint* ctr = (TernaryConstraint*)elimTernConstrs[newIndex];
    ctr->fillElimConstr(x, y, z, from1);
    ctr->isDuplicate_ = false;
    return ctr;
}

// warning! Do not propagate this new ternary cost function
TernaryConstraint* WCSP::newTernaryConstr(EnumeratedVariable* x, EnumeratedVariable* y, EnumeratedVariable* z, vector<Cost>& costs)
{
    unsigned int a, b;
    vector<Cost> zerocostsxy;
    vector<Cost> zerocostsxz;
    vector<Cost> zerocostsyz;

    for (a = 0; a < x->getDomainInitSize(); a++) {
        for (b = 0; b < y->getDomainInitSize(); b++) {
            zerocostsxy.push_back(MIN_COST);
        }
    }
    for (a = 0; a < x->getDomainInitSize(); a++) {
        for (b = 0; b < z->getDomainInitSize(); b++) {
            zerocostsxz.push_back(MIN_COST);
        }
    }
    for (a = 0; a < y->getDomainInitSize(); a++) {
        for (b = 0; b < z->getDomainInitSize(); b++) {
            zerocostsyz.push_back(MIN_COST);
        }
    }

    BinaryConstraint* xy = x->getConstr(y);
    BinaryConstraint* xz = x->getConstr(z);
    BinaryConstraint* yz = y->getConstr(z);

    if (!ToulBar2::vac) {
        if (!xy) {
            xy = new BinaryConstraint(this, x, y, zerocostsxy);
            xy->deconnect(true);
        }
        if (!xz) {
            xz = new BinaryConstraint(this, x, z, zerocostsxz);
            xz->deconnect(true);
        }
        if (!yz) {
            yz = new BinaryConstraint(this, y, z, zerocostsyz);
            yz->deconnect(true);
        }
    } else {
        if (!xy) {
            xy = new VACBinaryConstraint(this, x, y, zerocostsxy);
            xy->deconnect(true);
        }
        if (!xz) {
            xz = new VACBinaryConstraint(this, x, z, zerocostsxz);
            xz->deconnect(true);
        }
        if (!yz) {
            yz = new VACBinaryConstraint(this, y, z, zerocostsyz);
            yz->deconnect(true);
        }
    }

    TernaryConstraint* ctr = new TernaryConstraint(this, x, y, z, xy, xz, yz, costs);
    return ctr;
}

Constraint* WCSP::sum(Constraint* ctr1, Constraint* ctr2)
{
    assert(ctr1 != ctr2);
    if (ctr1->order(ctr2) < 0) {
        Constraint* ctraux = ctr1;
        ctr1 = ctr2;
        ctr2 = ctraux;
    }
    if (ToulBar2::verbose >= 1)
        cout << endl
             << "Sum of constraints: " << *ctr1 << " " << *ctr2 << endl;
    assert(ctr1->connected());
    assert(ctr2->connected());
    ctr1->deconnect();
    ctr2->deconnect(true);

    TSCOPE scopeUinv;
    TSCOPE scopeIinv;
    ctr1->scopeUnion(scopeUinv, ctr2);
    ctr1->scopeCommon(scopeIinv, ctr2);
    int arityU = scopeUinv.size();
    int arityI = scopeIinv.size();

    if (arityU == ctr2->arity() && ctr2->extension()) {
        ctr2->sumScopeIncluded(ctr1);
        ctr2->reconnect();
        ctr2->propagate();
        if (ToulBar2::verbose >= 1)
            cout << endl
                 << "Scopes Included.  Has result: " << *ctr2 << endl;
        return ctr2;
    }

    EnumeratedVariable** scopeU = new EnumeratedVariable*[arityU];
    EnumeratedVariable** scopeI = new EnumeratedVariable*[arityI];
    int* scopeUi = new int[arityU];

    int i = 0;
    TSCOPE::iterator it = scopeIinv.begin();
    while (it != scopeIinv.end()) {
        int xi = it->first;
        assert(xi == vars[xi]->wcspIndex);
        scopeU[i] = (EnumeratedVariable*)vars[xi];
        scopeI[i] = scopeU[i];
        scopeUi[i] = vars[xi]->wcspIndex;
        it++;
        i++;
        scopeUinv.erase(xi);
    }
    it = scopeUinv.begin();
    while (it != scopeUinv.end()) {
        int xi = it->first;
        scopeU[i] = (EnumeratedVariable*)vars[xi];
        scopeUi[i] = vars[xi]->wcspIndex;
        it++;
        i++;
    }

    EnumeratedVariable* x = scopeU[0];
    EnumeratedVariable* y = scopeU[1];

    Cost Top = getUb();
    unsigned int vxi, vyi, vzi;
    Tuple tuple, tuple1, tuple2;
    Cost cost, cost1, cost2;
    int ctrIndex = -INT_MAX;
    Constraint* ctr = NULL;
    vector<Cost> costs;

    if (arityU > NARYPROJECTIONSIZE) { // || isGlobal()) {
        ctrIndex = postNaryConstraintBegin(scopeUi, arityU, Top, ctr1->size() * ctr2->size(), true); //TODO: improve estimated number of tuples
        ctr = getCtr(ctrIndex);
        assert(ctr->isNary());
        NaryConstraint* nary = (NaryConstraint*)ctr;

        nary->fillFilters();

        bool tupleXtuple = (ctr1->getDefCost() >= Top) && (ctr2->getDefCost() >= Top);

        if (tupleXtuple) {
            ctr1->first();
            while (ctr1->next(tuple1, cost1)) {
                ctr2->first();
                while (ctr2->next(tuple2, cost2)) {
                    nary->insertSum(tuple1, cost1, ctr1, tuple2, cost2, ctr2, true);
                }
            }
        } else {
            nary->firstlex();
            while (nary->nextlex(tuple, cost)) {
                cost1 = ctr1->evalsubstr(tuple, nary);
                cost2 = ctr2->evalsubstr(tuple, nary);
                if (cost1 + cost2 < Top)
                    nary->setTuple(tuple, cost1 + cost2);
            }
        }
    } else if (arityU == 3) {
        EnumeratedVariable* z = scopeU[2];
        EnumeratedVariable* scopeTernary[3] = { x, y, z };
        Tuple t(3, 0);
        for (vxi = 0; vxi < x->getDomainInitSize(); vxi++)
            for (vyi = 0; vyi < y->getDomainInitSize(); vyi++)
                for (vzi = 0; vzi < z->getDomainInitSize(); vzi++) {
                    Value vx = x->toValue(vxi);
                    Value vy = y->toValue(vyi);
                    Value vz = z->toValue(vzi);
                    Cost costsum = Top;
                    if (x->canbe(vx) && y->canbe(vy) && z->canbe(vz)) {
                        costsum = MIN_COST;
                        if (arityI == 1) {
                            assert(ctr1->isBinary());
                            assert(ctr2->isBinary());
                            costsum += ((BinaryConstraint*)ctr1)->getCost(x, y, vx, vy)
                                + ((BinaryConstraint*)ctr2)->getCost(x, z, vx, vz);
                        } else if (arityI == 2) {
                            assert(ctr1->isBinary());
                            Cost c = MIN_COST;
                            if (ctr2->isTernary()) {
                                c = ((TernaryConstraint*)ctr2)->getCost(x, y, z, vx, vy, vz);
                            } else {
                                assert(ctr2->isNary());
                                t[0] = vxi;
                                t[1] = vyi;
                                t[2] = vzi;
                                c = ((NaryConstraint*)ctr2)->eval(t, scopeTernary);
                            }
                            costsum += ((BinaryConstraint*)ctr1)->getCost(x, y, vx, vy) + c;
                        } else if (arityI == 3) {
                            assert(ctr1->isTernary());
                            assert(ctr2->isTernary());
                            Cost c1 = MIN_COST;
                            Cost c2 = MIN_COST;
                            if (ctr1->isTernary()) {
                                c1 = ((TernaryConstraint*)ctr1)->getCost(x, y, z, vx, vy, vz);
                            } else {
                                assert(ctr1->isNary());
                                t[0] = vxi;
                                t[1] = vyi;
                                t[2] = vzi;
                                c1 = ((NaryConstraint*)ctr1)->eval(t, scopeTernary);
                            }
                            if (ctr2->isTernary()) {
                                c2 = ((TernaryConstraint*)ctr2)->getCost(x, y, z, vx, vy, vz);
                            } else {
                                assert(ctr2->isNary());
                                t[0] = vxi;
                                t[1] = vyi;
                                t[2] = vzi;
                                c2 = ((NaryConstraint*)ctr2)->eval(t, scopeTernary);
                            }
                            costsum += c1 + c2;
                        } else {
                            assert(false);
                        }
                        if (costsum > Top)
                            costsum = Top;
                    }
                    costs.push_back(costsum);
                }
        ctrIndex = postTernaryConstraint(x->wcspIndex, y->wcspIndex, z->wcspIndex, costs);
    } else if (arityU == 2) {
        BinaryConstraint* bctr1 = (BinaryConstraint*)ctr1;
        BinaryConstraint* bctr2 = (BinaryConstraint*)ctr2;
        for (vxi = 0; vxi < x->getDomainInitSize(); vxi++)
            for (vyi = 0; vyi < y->getDomainInitSize(); vyi++) {
                Value vx = x->toValue(vxi);
                Value vy = y->toValue(vyi);
                Cost costsum = Top;
                if (x->canbe(vx) && y->canbe(vy)) {
                    costsum = bctr1->getCost(x, y, vx, vy) + bctr2->getCost(x, y, vx, vy);
                    if (costsum > Top)
                        costsum = Top;
                }
                costs.push_back(costsum);
            }
        ctrIndex = postBinaryConstraint(x->wcspIndex, y->wcspIndex, costs);
    }
    assert(ctrIndex > -INT_MAX);
    delete[] scopeU;
    delete[] scopeUi;
    delete[] scopeI;
    ctr = getCtr(ctrIndex);
    ctr->propagate();
    if (ToulBar2::verbose >= 1)
        cout << endl
             << "Has result: " << *ctr << endl;
    return ctr;
}

void WCSP::project(Constraint*& ctr_inout, EnumeratedVariable* var, Constraint* ctr_copy)
{
    if (ctr_inout->getIndex(var) < 0)
        return;
    unsigned int vxi, vyi, vzi;
    int arity = ctr_inout->arity();
    int truearity = arity;
    if (arity >= 5) { // n-ary (n>=5) with only 4 unassigned variables is projected on a ternary!
        truearity = 0;
        for (int i = 0; i < arity; i++)
            if (ctr_inout->getVar(i)->unassigned())
                truearity++;
        assert(truearity >= 4);
    }

    if (ToulBar2::verbose >= 1)
        cout << endl
             << "Projection of var " << var->wcspIndex << " in ctr: " << *ctr_inout
             << endl;

    if (truearity - 1 > 3) {
        if (!ctr_inout->isNary()) {
            assert(ctr_copy && ctr_copy->isNary());
            ctr_inout->deconnect();
            ctr_inout = ctr_copy->copy();
            ctr_inout->reconnect();
        }
        assert(ctr_inout->isNary());
        ((NaryConstraint*)ctr_inout)->project(var);
        ctr_inout->propagate();
        if (ToulBar2::verbose >= 1)
            cout << endl
                 << "   has result*: " << *ctr_inout << endl;
        return;
    }
    ctr_inout->deconnect();

    int i, j;
    int ivars[3];
    EnumeratedVariable* evars[3];

    j = 0;
    for (i = 0; i < arity; i++) {
        EnumeratedVariable* v = (EnumeratedVariable*)ctr_inout->getVar(i);
        if (v != var && (arity <= 4 || v->unassigned())) {
            ivars[j] = v->wcspIndex;
            evars[j] = v;
            j++;
        }
    }

    Constraint* ctr = NULL;
    Cost Top = MAX_COST; // getUb();
    int ctrIndex;
    Tuple t(arity, 0);
    vector<Cost> costs;
    Cost negcost = 0;

    switch (truearity - 1) {
    case 3: {
        bool isnary = ctr_inout->isNary();
        if (truearity == 4 && arity >= 5) {
            for (i = 0; i < arity; i++) {
                if (ctr_inout->getVar(i)->assigned())
                    t[ctr_inout->getIndex(ctr_inout->getVar(i))] = ((EnumeratedVariable*)ctr_inout->getVar(i))->toIndex(ctr_inout->getVar(i)->getValue());
            }
        }
        for (vxi = 0; vxi < evars[0]->getDomainInitSize(); vxi++) {
            for (vyi = 0; vyi < evars[1]->getDomainInitSize(); vyi++) {
                for (vzi = 0; vzi < evars[2]->getDomainInitSize(); vzi++) {
                    Value v0 = evars[0]->toValue(vxi);
                    Value v1 = evars[1]->toValue(vyi);
                    Value v2 = evars[2]->toValue(vzi);
                    Cost mincost = Top;
                    if (evars[0]->canbe(v0) && evars[1]->canbe(v1) && evars[2]->canbe(v2)) {
                        t[ctr_inout->getIndex(evars[0])] = vxi;
                        t[ctr_inout->getIndex(evars[1])] = vyi;
                        t[ctr_inout->getIndex(evars[2])] = vzi;

                        for (EnumeratedVariable::iterator itv = var->begin(); itv != var->end(); ++itv) {
                            t[ctr_inout->getIndex(var)] = var->toIndex(*itv);
                            Cost c = ((isnary) ? ((NaryConstraint*)ctr_inout)->eval(t) : ctr_inout->evalsubstr(t, ctr_inout)) + var->getCost(*itv);
                            if (ToulBar2::isZ)
                                mincost = LogSumExp(mincost, c);
                            else if (c < mincost)
                                mincost = c;
                        }
                    }
                    if (ToulBar2::isZ && mincost < negcost)
                        negcost = mincost;
                    costs.push_back(mincost);
                }
            }
        }
        assert(negcost <= 0);
        if (negcost < 0) {
            for (vxi = 0; vxi < evars[0]->getDomainInitSize(); vxi++) {
                for (vyi = 0; vyi < evars[1]->getDomainInitSize(); vyi++) {
                    for (vzi = 0; vzi < evars[2]->getDomainInitSize(); vzi++) {
                        costs[vxi * evars[1]->getDomainInitSize() * evars[2]->getDomainInitSize() + vyi * evars[2]->getDomainInitSize() + vzi] -= negcost;
                    }
                }
            }
            decreaseLb(negcost);
        }
        ctrIndex = postTernaryConstraint(ivars[0], ivars[1], ivars[2], costs);
        ctr = getCtr(ctrIndex);
    } break;
    case 2: {
        TernaryConstraint* tctr = (TernaryConstraint*)ctr_inout;
        for (vxi = 0; vxi < evars[0]->getDomainInitSize(); vxi++)
            for (vyi = 0; vyi < evars[1]->getDomainInitSize(); vyi++) {
                Value v0 = evars[0]->toValue(vxi);
                Value v1 = evars[1]->toValue(vyi);
                Cost mincost = Top;
                if (evars[0]->canbe(v0) && evars[1]->canbe(v1)) {
                    for (EnumeratedVariable::iterator itv = var->begin(); itv != var->end(); ++itv) {
                        Cost c = tctr->getCost(evars[0], evars[1], var, v0, v1, *itv) + var->getCost(*itv);
                        if (ToulBar2::isZ)
                            mincost = LogSumExp(mincost, c);
                        else if (c < mincost)
                            mincost = c;
                    }
                }
                if (ToulBar2::isZ && mincost < negcost)
                    negcost = mincost;
                costs.push_back(mincost);
            }
        assert(negcost <= 0);
        if (negcost < 0) {
            for (vxi = 0; vxi < evars[0]->getDomainInitSize(); vxi++) {
                for (vyi = 0; vyi < evars[1]->getDomainInitSize(); vyi++) {
                    costs[vxi * evars[1]->getDomainInitSize() + vyi] -= negcost;
                }
            }
            decreaseLb(negcost);
        }
        ctrIndex = postBinaryConstraint(ivars[0], ivars[1], costs);
        ctr = getCtr(ctrIndex);
    } break;
    case 1: {
        BinaryConstraint* bctr = ((BinaryConstraint*)ctr_inout);
        for (vxi = 0; vxi < evars[0]->getDomainInitSize(); vxi++) {
            Value v0 = evars[0]->toValue(vxi);
            Cost mincost = Top;
            if (evars[0]->canbe(v0)) {
                for (EnumeratedVariable::iterator itv = var->begin(); itv != var->end(); ++itv) {
                    Cost c = bctr->getCost(evars[0], var, v0, *itv) + var->getCost(*itv);
                    if (ToulBar2::isZ)
                        mincost = LogSumExp(mincost, c);
                    else if (c < mincost)
                        mincost = c;
                }
            }
            if (ToulBar2::isZ && mincost < negcost)
                negcost = mincost;
            costs.push_back(mincost);
        }
        assert(negcost <= 0);
        for (EnumeratedVariable::iterator itv0 = evars[0]->begin(); itv0 != evars[0]->end(); ++itv0) {
            vxi = evars[0]->toIndex(*itv0);
            if (costs[vxi] - negcost > MIN_COST)
                evars[0]->project(*itv0, costs[vxi] - negcost);
        }
        evars[0]->findSupport();
        if (negcost < 0)
            decreaseLb(negcost);
    } break;
    default: {
        cerr << "Bad resulting cost function arity during generic variable elimination!";
        exit(EXIT_FAILURE);
    }
    }
    ctr_inout = ctr;
    if (ctr) {
        ctr->propagate();
        if (ToulBar2::verbose >= 1)
            cout << endl
                 << "   has result: " << *ctr_inout << endl;
    }
}

void WCSP::variableElimination(EnumeratedVariable* var)
{
    int degree = var->getTrueDegree();
    if (ToulBar2::verbose >= 1)
        cout << endl
             << "Generic variable elimination of " << var->getName() << "    degree: "
             << var->getDegree() << " true degree: " << degree << " max elim size: " << var->getMaxElimSize() << endl;
    if (degree > maxDegree)
        maxDegree = degree;

    if (var->getDegree() > 0) {

        ConstraintList::iterator it1 = var->getConstrs()->begin();
        ConstraintList::iterator it2;
        Constraint* c1 = (*it1).constr;
        Constraint* c2 = NULL;
        Constraint* csum = c1;
        Constraint* csumcopy = NULL;

        while (var->getDegree() > 1) {
            it1 = var->getConstrs()->begin();
            it2 = var->getConstrs()->rbegin();
            c1 = (*it1).constr;
            c2 = (*it2).constr;
            csum = sum(c1, c2);

            if (getTreeDec()) {
                csum->setCluster(var->getCluster());
            }
        }

        assert(csum->getIndex(var) >= 0);
        csumcopy = csum->copy();
        assert(csumcopy != NULL);
        elimInfo ei = { var, NULL, NULL, NULL, NULL, NULL, csumcopy };
        elimInfos[getElimOrder()] = ei;
        elimOrderInc();
        elimSpace += csumcopy->space();
        project(csum, var, csumcopy);
    } else {
        if (ToulBar2::isZ) { // add all unary loglike into lowerbound or negCost
            Cost clogz = MAX_COST;
            for (EnumeratedVariable::iterator itv = var->begin(); itv != var->end(); ++itv) {
                clogz = LogSumExp(clogz, var->getCost(*itv));
            }
            if (clogz < 0)
                decreaseLb(clogz);
            else
                increaseLb(clogz);
        }
    }
    assert(var->getDegree() == 0);

    var->assign(var->getSupport()); // warning! dummy assigned value
}

bool WCSP::kconsistency(int xIndex, int yIndex, int zIndex, BinaryConstraint* xy, BinaryConstraint* yz, BinaryConstraint* xz)
{
    if ((xIndex == yIndex) || (xIndex == zIndex) || (yIndex == zIndex))
        return false;
    EnumeratedVariable* x = (EnumeratedVariable*)vars[xIndex];
    EnumeratedVariable* y = (EnumeratedVariable*)vars[yIndex];
    EnumeratedVariable* z = (EnumeratedVariable*)vars[zIndex];
    TernaryConstraint* tctr = x->getConstr(y, z);
    if (tctr)
        return false;

    bool added = false;
    vector<Cost> costs;
    Cost ub = getUb();
    Cost minc = ub;
    for (unsigned int a = 0; a < x->getDomainInitSize(); a++) {
        Value va = x->toValue(a);
        Cost costa = x->getCost(va);
        for (unsigned int b = 0; b < y->getDomainInitSize(); b++) {
            Value vb = y->toValue(b);
            Cost costb = y->getCost(vb);
            for (unsigned int c = 0; c < z->getDomainInitSize(); c++) {
                Value vc = z->toValue(c);
                Cost costc = z->getCost(vc);
                Cost ctuple = ub;
                if (x->canbe(va) && y->canbe(vb) && z->canbe(vc)) {
                    ctuple = costa + costb + costc + xy->getCost(x, y, va, vb) + xz->getCost(x, z, va, vc)
                        + yz->getCost(y, z, vb, vc);
                }
                if (ctuple < minc)
                    minc = ctuple;
                costs.push_back(ctuple);
            }
        }
    }

    if (minc > MIN_COST) {
        tctr = new TernaryConstraint(this);
        elimTernConstrs.push_back(tctr);
        tctr = newTernaryConstr(x, y, z);
        tctr->fillElimConstrBinaries();
        tctr->reconnect();
        elimTernOrderInc();

        vector<Cost>::iterator it = costs.begin();
        for (unsigned int a = 0; a < x->getDomainInitSize(); a++) {
            Value va = x->toValue(a);
            Cost costa = x->getCost(va);
            if (x->canbe(va))
                x->extend(va, costa);
            for (unsigned int b = 0; b < y->getDomainInitSize(); b++) {
                Value vb = y->toValue(b);
                Cost costb = y->getCost(vb);
                if (y->canbe(vb))
                    y->extend(vb, costb);
                for (unsigned int c = 0; c < z->getDomainInitSize(); c++) {
                    Value vc = z->toValue(c);
                    Cost costc = z->getCost(vc);
                    if (z->canbe(vc))
                        z->extend(vc, costc);
                    if (x->canbe(va) && y->canbe(vb)) {
                        Cost costab = xy->getCost(x, y, va, vb);
                        xy->addcost(x, y, va, vb, -costab);
                    }
                    if (y->canbe(vb) && z->canbe(vc)) {
                        Cost costbc = yz->getCost(y, z, vb, vc);
                        yz->addcost(y, z, vb, vc, -costbc);
                    }
                    if (x->canbe(va) && z->canbe(vc)) {
                        Cost costac = xz->getCost(x, z, va, vc);
                        xz->addcost(x, z, va, vc, -costac);
                    }
                    tctr->setcost(x, y, z, va, vb, vc, *it - minc);
                    ++it;
                }
            }
        }
        tctr->projectTernary();
        increaseLb(minc);
        if (ToulBar2::verbose >= 1)
            cout << "new ternary(" << x->wcspIndex << "," << y->wcspIndex << ","
                 << z->wcspIndex << ")  newDualBound: " << std::fixed << std::setprecision(ToulBar2::decimalPoint) << getDDualBound() << std::setprecision(DECIMAL_POINT) << endl;
        added = true;
    }
    return added;
}

void WCSP::ternaryCompletion()
{
    if (numberOfUnassignedVariables() < 3)
        return;

    Double nbunvars = numberOfUnassignedVariables();
    Double connectivity = 2. * numberOfConnectedBinaryConstraints() / (nbunvars * (nbunvars - 1));
    Double domsize = medianDomainSize();
    Double size = domsize;
    size = sizeof(StoreCost) * size * size * size * nbunvars * (nbunvars - 1) * (nbunvars - 2) * connectivity * connectivity * connectivity / 6;

    //if (ToulBar2::debug>=2) cout << "MAX ESTIMATED RPC SIZE: " << size << " (" << nbunvars << "," << connectivity <<")" << endl;
    //if (size > 1024. * 1024. * ToulBar2::preprocessTernaryRPC) {
    //cout << "Restricted path consistency disabled (" << size/1024./1024. << " >= " << ToulBar2::preprocessTernaryRPC << " MB)" << endl;
    //ToulBar2::preprocessTernaryRPC = 0;
    //return;
    //}

    vector<TripleVarCostSize> triplelist;
    for (unsigned int i = 0; i < vars.size(); i++) {
        EnumeratedVariable* x = (EnumeratedVariable*)vars[i];
        for (ConstraintList::iterator it = x->getConstrs()->begin(); it != x->getConstrs()->end(); ++it) {
            Constraint* ctr = (*it).constr;
            if (ctr->isBinary()) {
                BinaryConstraint* bctr = (BinaryConstraint*)ctr;
                EnumeratedVariable* y = (EnumeratedVariable*)bctr->getVarDiffFrom(x);
                if (y->wcspIndex > x->wcspIndex)
                    for (ConstraintList::iterator it2 = y->getConstrs()->begin(); it2 != y->getConstrs()->end(); ++it2) {
                        Constraint* ctr2 = (*it2).constr;
                        if (ctr != ctr2 && ctr2->isBinary()) {
                            BinaryConstraint* bctr2 = (BinaryConstraint*)ctr2;
                            EnumeratedVariable* z = (EnumeratedVariable*)bctr2->getVarDiffFrom(y);
                            if (z->wcspIndex > y->wcspIndex)
                                for (ConstraintList::iterator it3 = z->getConstrs()->begin(); it3 != z->getConstrs()->end(); ++it3) {
                                    Constraint* ctr3 = (*it3).constr;
                                    if (ctr2 != ctr3 && ctr3->isBinary()) {
                                        BinaryConstraint* bctr3 = (BinaryConstraint*)ctr3;
                                        EnumeratedVariable* xx = (EnumeratedVariable*)bctr3->getVarDiffFrom(z);
                                        if (x == xx) {
                                            // bool added = kconsistency(x->wcspIndex,
                                            // 		y->wcspIndex, z->wcspIndex, bctr,
                                            // 		bctr2, bctr3);
                                            // if (added)
                                            float tight = bctr->computeTightness() + bctr2->computeTightness() + bctr3->computeTightness();
                                            long unsigned xsize = x->getDomainInitSize();
                                            long unsigned ysize = y->getDomainInitSize();
                                            long unsigned zsize = z->getDomainInitSize();
                                            TripleVarCostSize tvcs = { x, y, z, tight, xsize * ysize * zsize };
                                            triplelist.push_back(tvcs);
                                        }
                                    }
                                }
                        }
                    }
            }
        }
    }

    Double totalsize = 0.;
    Double maxsize = 1024. * 1024. * ToulBar2::preprocessTernaryRPC;
    int ntern = 0;

    sort(triplelist.begin(), triplelist.end());
    for (vector<TripleVarCostSize>::iterator it = triplelist.begin(); it != triplelist.end(); ++it) {
        if (totalsize + (Double)sizeof(StoreCost) * it->size <= maxsize) {
            totalsize += (Double)sizeof(StoreCost) * it->size;
            vector<Cost> costs(it->size, MIN_COST);
            postTernaryConstraint(it->x->wcspIndex, it->y->wcspIndex, it->z->wcspIndex, costs);
            ntern++;
        }
    }
    if (ToulBar2::verbose >= 0)
        cout << "Added " << ntern << "/" << triplelist.size() << " zero-cost ternary cost functions." << endl;
}

// -----------------------------------------------------------
// Methods for Virtual Arc Consistency

void WCSP::histogram(Cost c)
{
    if (vac)
        vac->histogram(c);
}
void WCSP::iniSingleton()
{
    if (vac)
        vac->iniSingleton();
}
void WCSP::updateSingleton()
{
    if (vac)
        vac->updateSingleton();
}
void WCSP::removeSingleton()
{
    if (vac)
        vac->removeSingleton();
}
void WCSP::printVACStat()
{
    if (vac)
        vac->printStat();
}

// -----------------------------------------------------------
// Methods for Cluster Tree Decomposition

bool WCSP::isAlreadyTreeDec(char* filename)
{
    if (filename == NULL || (reinterpret_cast<uintptr_t>(filename) > ELIM_NONE && reinterpret_cast<uintptr_t>(filename) < ELIM_MAX))
        return false;
    ifstream file;
    file.open(filename);
    if (!file)
        return false;
    int clusterid = 0;
    int parentid = 0;
    file >> clusterid;
    if (!file)
        return false;
    file >> parentid;
    file.close();
    if (parentid == -1)
        return true;
    return false;
}

void WCSP::buildTreeDecomposition()
{
    td = new TreeDecomposition(this);
    double time = cpuTime();
    if (isAlreadyTreeDec(ToulBar2::varOrder))
        td->buildFromCovering(ToulBar2::varOrder);
    else if (ToulBar2::approximateCountingBTD)
        td->buildFromOrderForApprox();
    else
        td->buildFromOrder();
    if (ToulBar2::verbose >= 0)
        cout << "Tree decomposition time: " << cpuTime() - time << " seconds." << endl;
    if (!ToulBar2::approximateCountingBTD) {
        vector<int> order;
        td->getElimVarOrder(order);
        // allows propagation to operate on the whole problem without modifying tree decomposition local lower bounds and delta costs
        // it is important for RDS-BTD which assumes zero cluster lower bounds and no delta cost moves
        TreeDecomposition* tmptd = td;
        td = NULL;
        setDACOrder(order);
        td = tmptd;
        // new constraints may be produced by variable elimination that must be correctly assigned to a cluster
        for (unsigned int i = 0; i < numberOfConstraints(); i++)
            if (constrs[i]->getCluster() == -1)
                constrs[i]->assignCluster();
        for (int i = 0; i < elimBinOrder; i++)
            if (elimBinConstrs[i]->connected() && elimBinConstrs[i]->getCluster() == -1)
                elimBinConstrs[i]->assignCluster();
        for (int i = 0; i < elimTernOrder; i++)
            if (elimTernConstrs[i]->connected() && elimTernConstrs[i]->getCluster() == -1)
                elimTernConstrs[i]->assignCluster();
        // check if ternary constraint cluster assignments are valid and do corrections if needed
        for (unsigned int i = 0; i < numberOfConstraints(); i++) {
            Constraint* ctr = getCtr(i);
            if (ctr->connected() && !ctr->isSep()) {
                if (ctr->isTernary()) {
                    TernaryConstraint* tctr = (TernaryConstraint*)ctr;
                    tctr->setDuplicates();
                    assert(tctr->xy->getCluster() == tctr->getCluster() && tctr->xz->getCluster() == tctr->getCluster() && tctr->yz->getCluster() == tctr->getCluster());
                }
            }
        }
        for (int i = 0; i < elimTernOrder; i++)
            if (elimTernConstrs[i]->connected()) {
                Constraint* ctr = elimTernConstrs[i];
                if (ctr->connected() && !ctr->isSep()) {
                    if (ctr->isTernary()) {
                        TernaryConstraint* tctr = (TernaryConstraint*)ctr;
                        tctr->setDuplicates();
                        assert(tctr->xy->getCluster() == tctr->getCluster() && tctr->xz->getCluster() == tctr->getCluster() && tctr->yz->getCluster() == tctr->getCluster());
                    }
                }
            }
    }
}

void WCSP::treeDecFile2Vector(char* filename, vector<int>& order)
{
    assert(order.size() == 0);
    map<int, int> clusterIds;
    int nbclusters = 0;

    set<int> usedvars;

    ifstream file(filename, ios::in);
    string fstr;
    while (getline(file, fstr)) {
        istringstream file(fstr);
        int num;
        file >> num;
        if (!file)
            break;

        clusterIds[num] = nbclusters;
        nbclusters++;

        int num_parent;
        file >> num_parent;
        assert((num_parent == -1) || (clusterIds.find(num_parent) != clusterIds.end()));

        int v = -1;
        while (file >> v) {
            if (usedvars.find(v) == usedvars.end()) {
                order.push_back(v);
                usedvars.insert(v);
            }
        }
    }
    file.close();

    reverse(order.begin(), order.end()); // must return an elimination order, the reverse of a topological order

    if (order.size() != numberOfVariables()) {
        cerr << "Tree decomposition file does not cover all the variables." << endl;
        exit(EXIT_FAILURE);
    }
}

void WCSP::elimOrderFile2Vector(char* elimVarOrder, vector<int>& order)
{
#ifdef BOOST
    if (reinterpret_cast<uintptr_t>(elimVarOrder) > ELIM_NONE && reinterpret_cast<uintptr_t>(elimVarOrder) < ELIM_MAX) {
        switch (reinterpret_cast<uintptr_t>(elimVarOrder)) {
        case MAX_CARD:
            maximumCardinalitySearch(order);
            break;
        case MIN_DEGREE:
            minimumDegreeOrdering(order);
            break;
        case MIN_FILL:
            minimumFillInOrdering(order);
            break;
        case ELIM_MST:
            spanningTreeOrderingBGL(order);
            break;
        case CUTHILL_MCKEE:
            reverseCuthillMcKeeOrderingBGL(order);
            break;
        case APPROX_MIN_DEGREE:
            minimumDegreeOrderingBGL(order);
            break;
        case ELIM_FILE_ORDER:
            order.clear();
            order.reserve(vars.size());
            for (int i = numberOfVariables() - 1; i >= 0; i--)
                order.push_back(i);
            break;
        default: {
            cerr << "Variable elimination order " << reinterpret_cast<uintptr_t>(elimVarOrder) << " not implemented yet!" << endl;
            exit(EXIT_FAILURE);
        }
        }
    } else {
#endif
        ifstream file;
        if (elimVarOrder)
            file.open(elimVarOrder);
        if (!elimVarOrder || !file) {
            if (ToulBar2::verbose >= 1) {
                cout << "Variable elimination order file missing or unreadable... takes reverse problem file variable index order." << endl;
            }
            //		for(unsigned int i=0;i<numberOfVariables();i++) order.push_back(i);
            for (int i = numberOfVariables() - 1; i >= 0; i--)
                order.push_back(i);
        } else {
            while (file) {
                int ix;
                file >> ix;
                if (file)
                    order.push_back(ix);
            }
        }
#ifdef BOOST
    }
#endif
    if (order.size() != numberOfVariables()) {
        cerr << "Variable elimination order file has incorrect number of variables." << endl;
        exit(EXIT_FAILURE);
    }
}

/// \param order variable elimination order (reverse of DAC order)
void WCSP::setDACOrder(vector<int>& order)
{
    if (order.size() != numberOfVariables()) {
        cerr << "DAC order has incorrect number of variables." << endl;
        exit(EXIT_FAILURE);
    }

    // set DAC order to the inverse of the elimination variable ordering
    if (ToulBar2::verbose >= 1)
        cout << "DAC order:";
    for (int i = order.size() - 1; i >= 0; i--) {
        if (ToulBar2::verbose >= 1)
            cout << " " << getVar(order[i])->getName();
        getVar(order[i])->setDACOrder(order.size() - 1 - i);
        if (ToulBar2::DEE >= 2)
            getVar(order[i])->queueDEE();
    }
    if (ToulBar2::verbose >= 1)
        cout << endl;

    sort(divVariables.begin(), divVariables.end(),
        [](const Variable* v1, const Variable* v2) -> bool {
            return (v1->getDACOrder() > v2->getDACOrder());
        });

    for (unsigned int i = 0; i < numberOfConstraints(); i++) {
        Constraint* ctr = getCtr(i);
        ctr->setDACScopeIndex();
        // Postpone global constraint propagation at the end (call to WCSP::propagate())
        if (ctr->connected() && !ctr->isGlobal() && !ctr->isSep())
            ctr->propagate();
    }
    for (int i = 0; i < elimBinOrder; i++) {
        Constraint* ctr = elimBinConstrs[i];
        ctr->setDACScopeIndex();
        if (ctr->connected())
            ctr->propagate();
    }
    for (int i = 0; i < elimTernOrder; i++) {
        Constraint* ctr = elimTernConstrs[i];
        ctr->setDACScopeIndex();
        if (ctr->connected())
            ctr->propagate();
    }
    propagate();
    // recompute all tightness: too slow???
    for (unsigned int i = 0; i < constrs.size(); i++)
        if (constrs[i]->connected())
            constrs[i]->computeTightness();
    for (int i = 0; i < elimBinOrder; i++)
        if (elimBinConstrs[i]->connected())
            elimBinConstrs[i]->computeTightness();
    for (int i = 0; i < elimTernOrder; i++)
        if (elimTernConstrs[i]->connected())
            elimTernConstrs[i]->computeTightness();
}

// -----------------------------------------------------------
// Functions for dealing with probabilities
// Warning: ToulBar2::NormFactor has to be initialized

// Converts the decimal Token to a cost and yells if unfeasible.
// the conversion uses the upper bound precision and ToulBar2::costMultiplier for scaling
// but does not shift cost using negCost. The input string should not contain white chars.

Cost WCSP::decimalToCost(const string& decimalToken, const unsigned int lineNumber) const
{
    size_t dotFound = decimalToken.find('.');
    size_t readIdx;
    if (dotFound == std::string::npos) {
        try {
            Cost cost = (Cost)std::stoll(decimalToken, &readIdx) * ToulBar2::costMultiplier * (Cost)powl(10, ToulBar2::decimalPoint);
            if (decimalToken[readIdx])
                throw std::invalid_argument("Not a cost");
            return cost;
        } catch (const std::invalid_argument&) {
            cerr << "Error: invalid cost '" << decimalToken;
            if (lineNumber)
                cerr << "' at line " << lineNumber << endl;
            else
                cerr << "' in command line" << endl;
            exit(1);
        }
    }

    bool negative = (decimalToken[0] == '-');
    string integerPart = (negative ? decimalToken.substr(1, dotFound - 1) : decimalToken.substr(0, dotFound));
    string decimalPart = decimalToken.substr(dotFound + 1);
    int shift = ToulBar2::decimalPoint - decimalPart.size();

    Cost cost;
    try {
        cost = (std::stoll(integerPart, &readIdx) * (Cost)powl(10, ToulBar2::decimalPoint) * ToulBar2::costMultiplier);
        if (integerPart[readIdx])
            throw std::invalid_argument("Not a cost");
        if (decimalPart.size()) {
            cost += std::stoll(decimalPart, &readIdx) * (Cost)powl(10, shift) * ToulBar2::costMultiplier;
            if (decimalPart[readIdx])
                throw std::invalid_argument("Not a cost");
        }
    } catch (const std::invalid_argument&) {
        cerr << "Error: invalid cost '" << decimalToken;
        if (lineNumber)
            cerr << "' at line " << lineNumber << endl;
        else
            cerr << "' in command line" << endl;
        exit(1);
    }
    if (negative)
        cost = -cost;
    return cost;
}

Cost WCSP::Prob2Cost(TProb p) const
{
    if (p == 0.0)
        return (MAX_COST - UNIT_COST) / MEDIUM_COST / MEDIUM_COST / MEDIUM_COST / MEDIUM_COST;
    TLogProb res = -Log(p) * ToulBar2::NormFactor;
    if (res > to_double(MAX_COST)) {
        cerr << "Overflow when converting probability to cost." << endl;
        exit(EXIT_FAILURE);
    }
    Cost c = (Cost)res;
    if (c > MAX_COST / 2)
        return (MAX_COST - UNIT_COST) / MEDIUM_COST / MEDIUM_COST / MEDIUM_COST / MEDIUM_COST;
    return c;
}

Cost WCSP::LogProb2Cost(TLogProb p) const
{
    TLogProb res = -p * ToulBar2::NormFactor;
    Cost c;
    if (res > to_double(MAX_COST / 2)) {
        c = (MAX_COST - UNIT_COST) / MEDIUM_COST / MEDIUM_COST / MEDIUM_COST / MEDIUM_COST;
        //        if (ToulBar2::verbose >= 2) {
        //            cout << "Warning: converting energy " << -p << " to Top\n";
        //        }
    } else
        c = (Cost)res;
    return c;
}

TLogProb WCSP::Cost2LogProb(Cost c) const
{
    return -to_double(c) / ToulBar2::NormFactor;
}

TProb WCSP::Cost2Prob(Cost c) const
{
    return Exp(-to_double(c) / ToulBar2::NormFactor);
}

Cost WCSP::LogSumExp(Cost c1, Cost c2) const // log[exp(c1) + exp(c2)]
{
    if (c1 >= MAX_COST / 2)
        return c2;
    else if (c2 >= MAX_COST / 2)
        return c1;
    else if (c1 == c2)
        return c1 + LogProb2Cost(Log(2.));
    else {
        if (c1 < c2)
            return c1 + LogProb2Cost(Log1p(Exp(Cost2LogProb(c2 - c1))));
        else
            return c2 + LogProb2Cost(Log1p(Exp(Cost2LogProb(c1 - c2))));
    }
}
TLogProb WCSP::LogSumExp(TLogProb logc1, Cost c2) const // log[exp(c1) + exp(c2)]
{
    TLogProb logc2 = Cost2LogProb(c2);
    if (logc1 == -numeric_limits<TLogProb>::infinity())
        return logc2;
    else if (c2 >= MAX_COST / 2)
        return logc1;
    else {
        if (logc1 >= logc2)
            return logc1 + (Log1p(Exp(logc2 - logc1)));
        else
            return logc2 + (Log1p(Exp(logc1 - logc2)));
    }
}

TLogProb WCSP::LogSumExp(TLogProb logc1, TLogProb logc2) const // log[exp(c1) + exp(c2)]
{
    if (logc1 == -numeric_limits<TLogProb>::infinity())
        return logc2;
    else if (logc2 == -numeric_limits<TLogProb>::infinity())
        return logc1;
    else {
        if (logc1 >= logc2)
            return logc1 + (Log1p(Exp(logc2 - logc1)));
        else
            return logc2 + (Log1p(Exp(logc1 - logc2)));
    }
}

//----------------------------------------
//procedure when berge acycl constant are present in the problem
// toulbar2::Berge_Dec has to be initialized > 0;

void WCSP::visit(int i, vector<int>& revdac, vector<bool>& marked, const vector<vector<int>>& listofsuccessors)
{
    marked[i] = true;
    for (unsigned int j = 0; j < listofsuccessors[i].size(); j++) {
        //   for (int  j = listofsuccessors[i].size()-1 ; j >= 0 ; j--) {
        if (!marked[listofsuccessors[i][j]])
            visit(listofsuccessors[i][j], revdac, marked, listofsuccessors);
    }
    revdac.push_back(i);
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

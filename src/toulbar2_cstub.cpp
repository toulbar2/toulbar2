#include <toulbar2lib.hpp>

extern "C" {

WeightedCSPSolver* tb2Create()
{
    return WeightedCSPSolver::makeWeightedCSPSolver(MAX_COST);
}

void tb2Destroy(WeightedCSPSolver* solver)
{
    delete solver;
}

void tb2Initialize(bool enable_vac)
{
    tb2init();
    ToulBar2::hbfs = false;
    ToulBar2::decimalPoint = 6;
    ToulBar2::vac = enable_vac;
    ToulBar2::vacValueHeuristic = enable_vac;
    ToulBar2::weightedTightness = true;
}

int tb2AddVariable(WeightedCSPSolver* solver, int size, char* name, char** valueNames)
{
    const static bool debug = false;
    if (debug)
        cout << "Entering tb2AddVariable\n";

    auto* prob = solver->getWCSP();
    int vIndex = prob->makeEnumeratedVariable(name, 0, size - 1);

    for (int i = 0; i < size; ++i) {
        string valName(valueNames[i]);
        prob->addValueName(vIndex, valName);
    }
    if (debug)
        cout << "Exiting tb2AddVariable\n";

    return vIndex;
}

void tb2AddFunction(WeightedCSPSolver* solver, int arity, const int* scope, const long double* costs)
{
    const static bool debug = false;
    if (debug)
        cout << "Entering tb2AddFunction\n";

    auto* prob = solver->getWCSP();
    assert(arity >= 0);
    size_t total = 1;

    for (int i = 0; i < arity; i++) {
        total *= prob->getDomainInitSize(scope[i]);
    }
    if (debug)
        cout << "Expected " << total << " values" << endl;

    std::vector<long double> cost_vector(costs, costs + total);

    switch (arity) {
    case 1:
        if (debug) {
            cout << "Posting unary on " << scope[0] << endl;
            for (long double c : cost_vector)
                cout << c << " ";
            cout << endl;
        }
        prob->postUnaryConstraint(scope[0], cost_vector);
        break;

    case 2:
        if (debug) {
            cout << "Posting binary on " << scope[0] << " " << scope[1] << endl;
            for (long double c : cost_vector)
                cout << c << " ";
            cout << endl;
        }
        prob->postBinaryConstraint(scope[0], scope[1], cost_vector);
        break;

    case 3:
        cerr << "Ternary cost functions unimplemented yet!\n";
        break;

    default:
        cerr << "N-ary cost functions unimplemented yet!\n";
        exit(1);
    }
    if (debug)
        cout << "Exiting tb2AddFunction\n";
}

bool tb2Solve(WeightedCSPSolver* solver)
{
    auto* prob = solver->getWCSP();
    prob->sortConstraints(); // Needs to be called before search.
    return solver->solve();
}

long double tb2GetSolution(WeightedCSPSolver* solver, Value* solution)
{
    long double energy = solver->getWCSP()->getDPrimalBound();
    const std::vector<int>& tb2solution = solver->getWCSP()->getSolution();
    for (int i : tb2solution)
        cout << i << " ";
    cout << endl;
    std::copy(tb2solution.begin(), tb2solution.end(), solution);
    return energy;
}

void tb2Debug(const bool debug)
{
    ToulBar2::debug = debug;
}

void tb2Dump(const int level, const char* problem)
{
    ToulBar2::dumpWCSP = level;
    ToulBar2::problemsaved_filename = to_string(problem);
}

void tb2UpdateUb(WeightedCSPSolver* solver, const long double newUb)
{
    auto* prob = solver->getWCSP();
    Cost ctUb = prob->DoubletoCost(newUb);
    prob->updateUb(ctUb);
}

void tb2NoPre()
{
    ToulBar2::elimDegree = -1;
    ToulBar2::elimDegree_preprocessing = -1;
    ToulBar2::preprocessTernaryRPC = 0;
    ToulBar2::preprocessFunctional = 0;
    ToulBar2::preprocessNary = 0;
    ToulBar2::costfuncSeparate = false;
    ToulBar2::MSTDAC = false;
    ToulBar2::DEE = 0;
}

void tb2BTDMode(const int mode)
{
    ToulBar2::btdMode = mode;
}

void tb2BTDRootCluster(const int rCluster)
{
    ToulBar2::btdRootCluster = rCluster;
}

void tb2BTDSubtree(const int sTree)
{
    ToulBar2::btdSubTree = sTree;
}

void tb2SplitClusterMaxSize(const int size)
{
    ToulBar2::splitClusterMaxSize = size;
}

void tb2BoostingBTD(const bool boost)
{
    ToulBar2::boostingBTD = boost;
}
void tb2MaxSeparatorSize(const int size)
{
    ToulBar2::maxSeparatorSize = size;
}

void tb2MinProperVarSize(const int size)
{
    ToulBar2::minProperVarSize = size;
}

void tb2WriteSolution(const char* write)
{
    char* filename = new char[strlen(write) + 1];
    strcpy(filename, write);
    ToulBar2::writeSolution = filename;
}

void tb2ShowSolutions(const bool show)
{
    ToulBar2::showSolutions = show;
}

void tb2AllSolutions(const long int sol)
{
    ToulBar2::allSolutions = sol;
}

void tb2ApproximateCountingBTD(const bool aproxim)
{
    ToulBar2::approximateCountingBTD = aproxim;
}

void tb2BinaryBranching(const bool boost)
{
    ToulBar2::binaryBranching = boost;
}

void tb2StaticVariableOrdering(const bool staticOrdering)
{
    ToulBar2::Static_variable_ordering = staticOrdering;
}

void tb2LastConflict(const bool last)
{
    ToulBar2::lastConflict = last;
}

void tb2DichotomicBranching(const int dicho)
{
    ToulBar2::dichotomicBranching = dicho;
}

void tb2SortDomains(const bool sort)
{
    ToulBar2::sortDomains = sort;
}

void tb2WeightedDegree(const int wDegree)
{
    ToulBar2::weightedDegree = wDegree;
}

void tb2WeightedTightness(const int wTight)
{
    ToulBar2::weightedTightness = wTight;
}

void tb2NbDecisionVars(const int nbDecision)
{
    ToulBar2::nbDecisionVars = nbDecision;
}

void tb2ElimDegree(const int degree)
{
    ToulBar2::elimDegree = degree;
}

void tb2ElimDegree_Preprocessing(const int degree_prepoc)
{
    ToulBar2::elimDegree_preprocessing = degree_prepoc;
}

void tb2ElimSpaceMaxMB(const int size)
{
    ToulBar2::elimSpaceMaxMB = size;
}

void tb2CostFuncSeparate(const bool separate)
{
    ToulBar2::costfuncSeparate = separate;
}

void tb2PartialAssign(WeightedCSPSolver* solver, const char* certificate)
{
    solver->parse_solution(certificate);
}

void tb2DeadEndElimination(const int level)
{
    ToulBar2::DEE = level;
}

void tb2VAC(const int depth)
{
    ToulBar2::vac = depth;
}

void tb2MinSumDiffusion(const int min)
{
    ToulBar2::minsumDiffusion = min;
}

void tb2CostThreshold(WeightedCSPSolver* solver, const long double cost)
{
    auto* prob = solver->getWCSP();
    ToulBar2::costThreshold = prob->DoubletoCost(cost);
}

void tb2CostThresholdPre(WeightedCSPSolver* solver, const long double cost)
{
    auto* prob = solver->getWCSP();
    ToulBar2::costThresholdPre = prob->DoubletoCost(cost);
}

void tb2CostMultiplier(const long double cost)
{
    ToulBar2::costMultiplier = cost;
}

void tb2SingletonConsistency(const bool singleCons)
{
    ToulBar2::singletonConsistency = singleCons;
}

void tb2VACValueHeuristic(const bool vacVal)
{
    ToulBar2::vacValueHeuristic = vacVal;
}

void tb2PreprocessTernaryRPC(const int size)
{
    ToulBar2::preprocessTernaryRPC = size;
}

void tb2PreprocessFunctional(const int func)
{
    ToulBar2::preprocessFunctional = func;
}

void tb2PreprocessNary(const int maxnary)
{
    ToulBar2::preprocessNary = maxnary;
}

void tb2QueueComplexity(const bool queue)
{
    ToulBar2::QueueComplexity = queue;
}

void tb2LDS(const int maxlds)
{
    ToulBar2::lds = maxlds;
}

void tb2Restart(const long maxrestarts)
{
    ToulBar2::restart = maxrestarts;
}

void tb2LCLevel(const int level)
{
    LcLevelType lclevel = (LcLevelType)level;
    ToulBar2::LcLevel = lclevel;
}

void tb2HBFS(const long hbfsgloballimit)
{
    ToulBar2::hbfs = 1;
    ToulBar2::hbfsGlobalLimit = hbfsgloballimit;
}

void tb2HBFSAlpha(const long hbfsalpha)
{
    ToulBar2::hbfsAlpha = hbfsalpha;
}

void tb2HBFSBeta(const long hbfsbeta)
{
    ToulBar2::hbfsBeta = hbfsbeta;
}

void tb2HBFSOpenNodeLimit(const long openlimit)
{
    ToulBar2::hbfs = 1;
    if (ToulBar2::hbfsGlobalLimit == 0)
        ToulBar2::hbfsGlobalLimit = 10000;
    ToulBar2::hbfsOpenNodeLimit = openlimit;
}

void tb2VariableEliminationOrder(const int order)
{
    ToulBar2::varOrder = reinterpret_cast<char*>(abs(order));
}

void tb2Incop(const char* cmd)
{
    if (cmd == NULL || strlen(cmd) == 0) {
        ToulBar2::incop_cmd = "0 1 3 idwa 100000 cv v 0 200 1 0 0";
    } else {
        char* cmd_ = new char[strlen(cmd) + 1];
        strcpy(cmd_, cmd);
        ToulBar2::incop_cmd = cmd_;
    }
}
}
/*
    solver = _swig_property(_Toulbar2.Toulbar2Solver_solver_get, _Toulbar2.Toulbar2Solver_solver_set)
    wcsp = _swig_property(_Toulbar2.Toulbar2Solver_wcsp_get, _Toulbar2.Toulbar2Solver_wcsp_set)
    upperbound = _swig_property(_Toulbar2.Toulbar2Solver_upperbound_get, _Toulbar2.Toulbar2Solver_upperbound_set)
    optimum = _swig_property(_Toulbar2.Toulbar2Solver_optimum_get, _Toulbar2.Toulbar2Solver_optimum_set)
    costshift = _swig_property(_Toulbar2.Toulbar2Solver_costshift_get, _Toulbar2.Toulbar2Solver_costshift_set)
    unsatisfiable = _swig_property(_Toulbar2.Toulbar2Solver_unsatisfiable_get, _Toulbar2.Toulbar2Solver_unsatisfiable_set)
    interrupted = _swig_property(_Toulbar2.Toulbar2Solver_interrupted_get, _Toulbar2.Toulbar2Solver_interrupted_set)
    solution = _swig_property(_Toulbar2.Toulbar2Solver_solution_get, _Toulbar2.Toulbar2Solver_solution_set)

    def add(self, *args): return _Toulbar2.Toulbar2Solver_add(self, *args)
    def initialise(self, *args): return _Toulbar2.Toulbar2Solver_initialise(self, *args)
    def propagate(self): return _Toulbar2.Toulbar2Solver_propagate(self)
    def solve(self): return _Toulbar2.Toulbar2Solver_solve(self)
    def solveAndRestart(self, *args): return _Toulbar2.Toulbar2Solver_solveAndRestart(self, *args)
    def store_solution(self): return _Toulbar2.Toulbar2Solver_store_solution(self)
    def is_opt(self): return _Toulbar2.Toulbar2Solver_is_opt(self)
    def is_sat(self): return _Toulbar2.Toulbar2Solver_is_sat(self)
    def is_unsat(self): return _Toulbar2.Toulbar2Solver_is_unsat(self)
    def getOptimum(self): return _Toulbar2.Toulbar2Solver_getOptimum(self)
    def getBacktracks(self): return _Toulbar2.Toulbar2Solver_getBacktracks(self)
    def getNodes(self): return _Toulbar2.Toulbar2Solver_getNodes(self)
    def getFailures(self): return _Toulbar2.Toulbar2Solver_getFailures(self)
    def getPropags(self): return _Toulbar2.Toulbar2Solver_getPropags(self)
    def getTime(self): return _Toulbar2.Toulbar2Solver_getTime(self)
    def getChecks(self): return _Toulbar2.Toulbar2Solver_getChecks(self)
    def printStatistics(self): return _Toulbar2.Toulbar2Solver_printStatistics(self)
    def getNumVariables(self): return _Toulbar2.Toulbar2Solver_getNumVariables(self)
    def getNumConstraints(self): return _Toulbar2.Toulbar2Solver_getNumConstraints(self)
    def printPython(self): return _Toulbar2.Toulbar2Solver_printPython(self)
    def setHeuristic(self, *args): return _Toulbar2.Toulbar2Solver_setHeuristic(self, *args)
    def setFailureLimit(self, *args): return _Toulbar2.Toulbar2Solver_setFailureLimit(self, *args)
    def setNodeLimit(self, *args): return _Toulbar2.Toulbar2Solver_setNodeLimit(self, *args)
    def setTimeLimit(self, *args): return _Toulbar2.Toulbar2Solver_setTimeLimit(self, *args)
    def setRandomized(self, *args): return _Toulbar2.Toulbar2Solver_setRandomized(self, *args)
    def setRandomSeed(self, *args): return _Toulbar2.Toulbar2Solver_setRandomSeed(self, *args)
    def setVerbosity(self, *args): return _Toulbar2.Toulbar2Solver_setVerbosity(self, *args)
    def debug(self, *args): return _Toulbar2.Toulbar2Solver_debug(self, *args)
    def writeSolution(self, *args): return _Toulbar2.Toulbar2Solver_writeSolution(self, *args)
    def showSolutions(self, *args): return _Toulbar2.Toulbar2Solver_showSolutions(self, *args)
    def dumpWCSP(self, *args): return _Toulbar2.Toulbar2Solver_dumpWCSP(self, *args)
    def nopre(self): return _Toulbar2.Toulbar2Solver_nopre(self)
    def updateUb(self, *args): return _Toulbar2.Toulbar2Solver_updateUb(self, *args)
    def lds(self, *args): return _Toulbar2.Toulbar2Solver_lds(self, *args)
    def restart(self, *args): return _Toulbar2.Toulbar2Solver_restart(self, *args)
    def hbfs(self, *args): return _Toulbar2.Toulbar2Solver_hbfs(self, *args)
    def hbfsAlpha(self, *args): return _Toulbar2.Toulbar2Solver_hbfsAlpha(self, *args)
    def hbfsBeta(self, *args): return _Toulbar2.Toulbar2Solver_hbfsBeta(self, *args)
    def hbfsOpenNodeLimit(self, *args): return _Toulbar2.Toulbar2Solver_hbfsOpenNodeLimit(self, *args)
    def lcLevel(self, *args): return _Toulbar2.Toulbar2Solver_lcLevel(self, *args)
    def QueueComplexity(self, *args): return _Toulbar2.Toulbar2Solver_QueueComplexity(self, *args)
    def allSolutions(self, *args): return _Toulbar2.Toulbar2Solver_allSolutions(self, *args)
    def approximateCountingBTD(self, *args): return _Toulbar2.Toulbar2Solver_approximateCountingBTD(self, *args)
    def binaryBranching(self, *args): return _Toulbar2.Toulbar2Solver_binaryBranching(self, *args)
    def staticVariableOrdering(self, *args): return _Toulbar2.Toulbar2Solver_staticVariableOrdering(self, *args)
    def lastConflict(self, *args): return _Toulbar2.Toulbar2Solver_lastConflict(self, *args)
    def dichotomicBranching(self, *args): return _Toulbar2.Toulbar2Solver_dichotomicBranching(self, *args)
    def sortDomains(self, *args): return _Toulbar2.Toulbar2Solver_sortDomains(self, *args)
    def weightedDegree(self, *args): return _Toulbar2.Toulbar2Solver_weightedDegree(self, *args)
    def weightedTightness(self, *args): return _Toulbar2.Toulbar2Solver_weightedTightness(self, *args)
    def variableEliminationOrdering(self, *args): return _Toulbar2.Toulbar2Solver_variableEliminationOrdering(self, *args)
    def nbDecisionVars(self, *args): return _Toulbar2.Toulbar2Solver_nbDecisionVars(self, *args)
    def partialAssign(self, *args): return _Toulbar2.Toulbar2Solver_partialAssign(self, *args)
    def elimDegree(self, *args): return _Toulbar2.Toulbar2Solver_elimDegree(self, *args)
    def elimDegree_preprocessing(self, *args): return _Toulbar2.Toulbar2Solver_elimDegree_preprocessing(self, *args)
    def elimSpaceMaxMB(self, *args): return _Toulbar2.Toulbar2Solver_elimSpaceMaxMB(self, *args)
    def costfuncSeparate(self, *args): return _Toulbar2.Toulbar2Solver_costfuncSeparate(self, *args)
    def deadEndElimination(self, *args): return _Toulbar2.Toulbar2Solver_deadEndElimination(self, *args)
    def preprocessTernaryRPC(self, *args): return _Toulbar2.Toulbar2Solver_preprocessTernaryRPC(self, *args)
    def preprocessFunctional(self, *args): return _Toulbar2.Toulbar2Solver_preprocessFunctional(self, *args)
    def preprocessNary(self, *args): return _Toulbar2.Toulbar2Solver_preprocessNary(self, *args)
    def btdMode(self, *args): return _Toulbar2.Toulbar2Solver_btdMode(self, *args)
    def splitClusterMaxSize(self, *args): return _Toulbar2.Toulbar2Solver_splitClusterMaxSize(self, *args)
    def maxSeparatorSize(self, *args): return _Toulbar2.Toulbar2Solver_maxSeparatorSize(self, *args)
    def minProperVarSize(self, *args): return _Toulbar2.Toulbar2Solver_minProperVarSize(self, *args)
    def boostingBTD(self, *args): return _Toulbar2.Toulbar2Solver_boostingBTD(self, *args)
    def btdRootCluster(self, *args): return _Toulbar2.Toulbar2Solver_btdRootCluster(self, *args)
    def btdSubTree(self, *args): return _Toulbar2.Toulbar2Solver_btdSubTree(self, *args)
    def vac(self, *args): return _Toulbar2.Toulbar2Solver_vac(self, *args)
    def vacValueHeuristic(self, *args): return _Toulbar2.Toulbar2Solver_vacValueHeuristic(self, *args)
    def costThreshold(self, *args): return _Toulbar2.Toulbar2Solver_costThreshold(self, *args)
    def costThresholdPre(self, *args): return _Toulbar2.Toulbar2Solver_costThresholdPre(self, *args)
    def costMultiplier(self, *args): return _Toulbar2.Toulbar2Solver_costMultiplier(self, *args)
    def singletonConsistency(self, *args): return _Toulbar2.Toulbar2Solver_singletonConsistency(self, *args)
    def minsumDiffusion(self, *args): return _Toulbar2.Toulbar2Solver_minsumDiffusion(self, *args)
    def incop(self, *args): return _Toulbar2.Toulbar2Solver_incop(self, *args)

*/

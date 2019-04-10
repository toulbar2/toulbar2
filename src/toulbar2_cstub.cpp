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
        if (debug)
            cout << "Expected " << total << " values" << endl;
    }
    std::vector<long double> cost_vector(costs, costs + total);

    switch (arity) {
    case 1:
        prob->postUnaryConstraint(scope[0], cost_vector);
        break;

    case 2:
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

void tb2GetSolution(WeightedCSPSolver* solver, Value* solution)
{
    const std::vector<int>& tb2solution = solver->getWCSP()->getSolution();
    std::copy(tb2solution.begin(), tb2solution.end(), solution);
}
}

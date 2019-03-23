#include <toulbar2lib.hpp>

extern "C" {
void toulbar2_initialize(bool enable_vac)
{
    tb2init();
    ToulBar2::hbfs = false;
    ToulBar2::decimalPoint = 6;
    ToulBar2::vac = enable_vac;
    ToulBar2::vacValueHeuristic = enable_vac;
    ToulBar2::weightedTightness = true;
}

WeightedCSPSolver* tb2Create()
{
    return WeightedCSPSolver::makeWeightedCSPSolver(MAX_COST);
}

void tb2Destroy(WeightedCSPSolver* solver)
{
    delete solver;
}

int tb2AddVariable(WeightedCSPSolver* solver, int size, char* name, char** valueNames)
{
    auto* prob = solver->getWCSP();
    int vIndex = prob->makeEnumeratedVariable(name, 0, size - 1);

    for (int i = 0; i < size; ++i) {
        string valName(valueNames[i]);
        prob->addValueName(vIndex, valName);
    }

    return vIndex;
}

void tb2AddFunction(WeightedCSPSolver* solver, int arity, const int* variables, const long double* costs)
{
    auto* prob = solver->getWCSP();
    assert(arity >= 0);
    size_t total = 1;

    for (int i = 0; i < arity; i++) {
        total *= prob->getDomainInitSize(i);
    }
    std::vector<long double> cost_vector(costs, costs + total);

    switch (arity) {
    case 1:
        prob->postUnaryConstraint(variables[0], cost_vector);
        break;

    case 2:
        prob->postBinaryConstraint(variables[0], variables[1], cost_vector);
        break;

    case 3:
        cerr << "Ternary cost functions unimplemented yet!\n";
        break;

    default:
        cerr << "N-ary cost functions unimplemented yet!\n";
        exit(1);
    }
}

bool toulbar2_solve(WeightedCSPSolver* solver)
{
    auto* prob = solver->getWCSP();
    prob->sortConstraints(); // Needs to be called before search.

    return solver->solve();
}

void toulbar2_get_solution(WeightedCSPSolver* solver, Value* solution)
{
    const std::vector<int>& tb2solution = solver->getWCSP()->getSolution();
    std::copy(tb2solution.begin(), tb2solution.end(), solution);
}

void toulbar2_test(int n, char** test)
{
    for (int i = 0; i < n; i++)
        cout << test[i] << endl;
}
}

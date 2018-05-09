/** \file iloglue.cpp
 *  \brief Link with Ilog Solver, adding a global soft constraint representing a weighted CSP and propagated by ToulBar2.
 * 
 */

// these includes are needed if compiled on new g++ versions (>4.0?)
#include <climits>
#include <cstdlib>
#include <cstring>

#include <ilsolver/ilosolver.h>
#include <ilsolver/ilctrace.h>
#include "toulbar2lib.hpp"
#include "tb2domain.hpp"
extern ostream& operator<<(ostream& os, WCSP& wcsp);

// TO BE DONE
// * avoid dedicated search goal:
// * - adds objective goal to the search
// * - encapsulates store/restore in IlcWeightedCSPI demons based on getDepth (store) and numberOfFails (restore) (take care of restore is done before wcsp method calls in ilog search goals!)
// * global wcsp constraint extended by using table constraints
// * why domain-delta getMinDelta and getMaxDelta methods do not work properly?
// * separate increase/decrease events from remove events sent towards ToulBar2
// * test two objective functions on spot5 problems
// * sport scheduling experiments

// use STL namespace
ILOSTLBEGIN

IlcIntVar Objective;
int ProblemSize = 0;
IlcIntVarArray ProblemVars;
int UpperBound = MAX_COST; // best solution cost or initial global upper bound
int* BestSol = NULL; // best solution found during the search

// current IloSolver instance used by libtb2.so to generate a failure
IloSolver IlogSolver;

// current WeightedCSP instance used by value and variable ordering heuristics
WeightedCSP* CurrentWeightedCSP = NULL;

// global weighted csp constraint exploiting toulbar2 propagation
class IlcWeightedCSPI : public IlcConstraintI {
public:
    static vector<IlcWeightedCSPI*> AllIlcWeightedCSPI;
    static int wcspCounter;

    IlcIntVar obj; // objective function
    int size; // |vars|
    IlcIntVarArray vars; // all Ilog variables involved in the WCSP network
    WeightedCSP* wcsp; // WCSP network managed by ToulBar2
    Domain* unassignedVars; // a WCSP domain containing the list of unassigned variables
    IlcInt currentNumberOfFails; // counter of search failures to inform ToulBar2 to reset its propagation queues and update its timestamp
    IlcRevBool synchronized; // if IlcTrue then force a complete synchronization between Ilog and ToulBar2 variable domains and objective

    // unique IlcWeightedCSPI constructor
    // creates an empty WCSP and add soft constraints from a file if available
    IlcWeightedCSPI(IloSolver solver,
        IlcIntVar objective, IlcIntVarArray variables,
        const char* fileName = NULL)
        : IlcConstraintI(solver)
        , obj(objective)
        , size(variables.getSize())
        , vars(variables)
        , wcsp(NULL)
        , unassignedVars(NULL)
        , currentNumberOfFails(0)
        , synchronized(solver, IlcTrue)
    {
        // creates a WCSP object
        wcsp = WeightedCSP::makeWeightedCSP(MAX_COST);
        // load WCSP problem from a file if available
        if (fileName) {
            wcsp->read_wcsp(fileName);
            assert((unsigned int)size == wcsp->numberOfVariables());
        }
        // specific data to check if all variables have been assigned
        unassignedVars = new Domain(0, size - 1);
        // memorizes all WeightedCSP instances
        assert(wcsp->getIndex() == wcspCounter);
        AllIlcWeightedCSPI.push_back(this);
        wcspCounter++;
        CurrentWeightedCSP = wcsp;
    }

    // destructor
    ~IlcWeightedCSPI()
    {
        AllIlcWeightedCSPI[wcsp->getIndex()] = NULL;
        delete wcsp;
        delete unassignedVars;
    }

    // domain synchronization between obj&vars (Ilog) and wcsp (ToulBar2)
    void synchronize()
    {
        if (ToulBar2::verbose >= 2)
            cout << "Domain synchronization between IlogSolver and Toulbar2!" << endl;
        for (int i = 0; i < size; i++) {
            if (ToulBar2::verbose >= 2)
                cout << vars[i] << " (" << wcsp->getInf(i) << "," << wcsp->getSup(i) << ")" << endl;
            vars[i].setMin(wcsp->getInf(i));
            vars[i].setMax(wcsp->getSup(i));
            for (int d = wcsp->getInf(i); d <= wcsp->getSup(i); d++) {
                if (wcsp->cannotbe(i, d)) {
                    vars[i].removeValue(d);
                }
            }
            wcsp->increase(i, vars[i].getMin());
            wcsp->decrease(i, vars[i].getMax());
            for (int d = vars[i].getMin(); d <= vars[i].getMax(); d++) {
                if (!vars[i].isInDomain(d)) {
                    wcsp->remove(i, d);
                }
            }
        }
        obj.setMin(wcsp->getLb());
        obj.setMax(wcsp->getUb() - 1);
        wcsp->decreaseUb(obj.getMax() + 1);
        UpperBound = wcsp->getUb();
    }

    // if a search node failure has just occured then informs ToulBar2 to reset its propagation queues and update its timestamp
    void checkFailure()
    {
        if (getSolver().getNumberOfFails() != currentNumberOfFails) {
            currentNumberOfFails = getSolver().getNumberOfFails();
            wcsp->whenContradiction();
        }
    }

    // links the WCSP variables to the ILOG variables
    void post();

    // global propagation using WCSP propagation queues
    void propagate()
    {
        checkFailure();
        if (synchronized) {
            synchronized.setValue(getSolver(), IlcFalse);
            synchronize();
        }
        if (ToulBar2::verbose >= 2)
            cout << "ILOG: propagate wcsp index " << wcsp->getIndex() << endl;
        wcsp->decreaseUb(obj.getMax() + 1);
        wcsp->propagate();
    }

    // variable varIndex has been assigned
    void whenValue(const IlcInt varIndex)
    {
        checkFailure();
        if (ToulBar2::verbose >= 2)
            cout << "ILOG: " << vars[varIndex].getName() << " = " << vars[varIndex].getValue() << endl;
        wcsp->assign(varIndex, vars[varIndex].getValue());
        if (unassignedVars->canbe(varIndex)) {
            unassignedVars->erase(varIndex);
            if (unassignedVars->empty()) {
                assert(wcsp->verify());
                obj.setValue(wcsp->getLb());
            }
        }
        push(); // global propagation done after local propagation
    }

    // check only modifications on the objective variable
    void whenRange()
    {
        checkFailure();
        wcsp->enforceUb(); // fail if lower bound >= upper bound and enforce NC*
        if (obj.getMax() + 1 < wcsp->getUb()) {
            wcsp->decreaseUb(obj.getMax() + 1);
        }
        push(); // global propagation done after local propagation
    }

    // variable varIndex has its domain reduced
    void whenDomain(const IlcInt varIndex)
    {
        checkFailure();
        if (!vars[varIndex].isBound()) {
            for (IlcIntVarDeltaIterator iter(vars[varIndex]); iter.ok(); ++iter) {
                IlcInt val = *iter;
                if (ToulBar2::verbose >= 2)
                    cout << "ILOG: " << vars[varIndex].getName() << " != " << val << endl;
                wcsp->remove(varIndex, val);
            }
            push(); // global propagation done after local propagation
        }
    }

    // Not implemented!
    //   IlcBool isViolated() const;
    //   void metaPostDemon(IlcDemonI*);
    //   IlcConstraintI* makeOpposite() const;
};

vector<IlcWeightedCSPI*> IlcWeightedCSPI::AllIlcWeightedCSPI;
int IlcWeightedCSPI::wcspCounter = 0;

void tb2setvalue(int wcspId, int varIndex, Value value, void* solver)
{
    assert(wcspId < IlcWeightedCSPI::wcspCounter);
    assert(IlcWeightedCSPI::AllIlcWeightedCSPI[wcspId] != NULL);
    assert(varIndex < IlcWeightedCSPI::AllIlcWeightedCSPI[wcspId]->size);
    if (ToulBar2::verbose >= 2)
        cout << "TOULBAR2: x" << varIndex << " = " << value << endl;
    IlcWeightedCSPI::AllIlcWeightedCSPI[wcspId]->vars[varIndex].setValue(value);
}

void tb2removevalue(int wcspId, int varIndex, Value value, void* solver)
{
    assert(wcspId < IlcWeightedCSPI::wcspCounter);
    assert(IlcWeightedCSPI::AllIlcWeightedCSPI[wcspId] != NULL);
    assert(varIndex < IlcWeightedCSPI::AllIlcWeightedCSPI[wcspId]->size);
    if (ToulBar2::verbose >= 2)
        cout << "TOULBAR2: x" << varIndex << " != " << value << endl;
    IlcWeightedCSPI::AllIlcWeightedCSPI[wcspId]->vars[varIndex].removeValue(value);
}

void tb2setmin(int wcspId, int varIndex, Value value, void* solver)
{
    assert(wcspId < IlcWeightedCSPI::wcspCounter);
    assert(IlcWeightedCSPI::AllIlcWeightedCSPI[wcspId] != NULL);
    assert(varIndex < IlcWeightedCSPI::AllIlcWeightedCSPI[wcspId]->size);
    if (ToulBar2::verbose >= 2)
        cout << "TOULBAR2: x" << varIndex << " >= " << value << endl;
    IlcWeightedCSPI::AllIlcWeightedCSPI[wcspId]->vars[varIndex].setMin(value);
}

void tb2setmax(int wcspId, int varIndex, Value value, void* solver)
{
    assert(wcspId < IlcWeightedCSPI::wcspCounter);
    assert(IlcWeightedCSPI::AllIlcWeightedCSPI[wcspId] != NULL);
    assert(varIndex < IlcWeightedCSPI::AllIlcWeightedCSPI[wcspId]->size);
    if (ToulBar2::verbose >= 2)
        cout << "TOULBAR2: x" << varIndex << " <= " << value << endl;
    IlcWeightedCSPI::AllIlcWeightedCSPI[wcspId]->vars[varIndex].setMax(value);
}

void tb2setminobj(int wcspId, int varIndex, Value value, void* solver)
{
    assert(wcspId < IlcWeightedCSPI::wcspCounter);
    assert(IlcWeightedCSPI::AllIlcWeightedCSPI[wcspId] != NULL);
    assert(varIndex == -1);
    if (ToulBar2::verbose >= 2)
        cout << "TOULBAR2: obj"
             << " >= " << value << endl;
    IlcWeightedCSPI::AllIlcWeightedCSPI[wcspId]->obj.setMin(value);
}

ILCCTDEMON1(IlcWeightedCSPWhenValueDemon, IlcWeightedCSPI, whenValue, IlcInt, varIndex);
ILCCTDEMON1(IlcWeightedCSPWhenDomainDemon, IlcWeightedCSPI, whenDomain, IlcInt, varIndex);
ILCCTDEMON0(IlcWeightedCSPWhenRangeDemon, IlcWeightedCSPI, whenRange);

void IlcWeightedCSPI::post()
{
    ToulBar2::setvalue = ::tb2setvalue;
    ToulBar2::removevalue = ::tb2removevalue;
    ToulBar2::setmin = ::tb2setmin;
    ToulBar2::setmax = ::tb2setmax;
    ToulBar2::setminobj = ::tb2setminobj;
    for (int i = 0; i < size; i++) {
        vars[i].whenValue(IlcWeightedCSPWhenValueDemon(getSolver(), this, i));
        vars[i].whenDomain(IlcWeightedCSPWhenDomainDemon(getSolver(), this, i));
    }
    obj.whenRange(IlcWeightedCSPWhenRangeDemon(getSolver(), this));
}

IlcConstraint IlcWeightedCSP(IlcIntVar objective, IlcIntVarArray variables, const char* filename)
{
    IloSolver solver = objective.getSolver();
    return IlcConstraint(new (solver.getHeap())
            IlcWeightedCSPI(solver, objective, variables, filename));
}
ILOCPCONSTRAINTWRAPPER3(IloWeightedCSP, solver, IloIntVar, obj, IloIntVarArray, vars,
    const char*, filename)
{
    use(solver, obj);
    use(solver, vars);
    return IlcWeightedCSP(solver.getIntVar(obj), solver.getIntVarArray(vars), filename);
}

ILCGOAL2(IlcGuess, IlcIntVar, var, IlcInt, value)
{
    Store::store();
    if (ToulBar2::verbose >= 1)
        cout << "[" << getSolver().getSearchNode().getDepth() << "," << Objective.getMin() << "," << UpperBound << "] Try " << var << " = " << value << endl;
    var.setValue(value);
    return 0;
}

ILCGOAL3(IlcRefute, IlcIntVar, var, IlcInt, value, IlcInt, depth)
{
    Store::restore(depth);
    Store::store(); // => store.getDepth() == getSolver().getSearchNode().getDepth()
    Objective.setMax(UpperBound - 1);
    if (ToulBar2::verbose >= 1)
        cout << "[" << getSolver().getSearchNode().getDepth() << "," << Objective.getMin() << "," << UpperBound << "] Refute " << var << " != " << value << endl;
    var.removeValue(value);
    return 0;
}

ILCGOAL0(IlcNewSolution)
{
    if (getSolver().isInSearch() && !getSolver().isInRecomputeMode()) {
        if (ToulBar2::verbose >= 0)
            cout << "New solution: " << Objective.getMin() << " (" << getSolver().getNumberOfFails() << " fails)" << endl;
        if (ToulBar2::showSolutions) {
            for (int i = 0; i < ProblemSize; i++) {
                cout << ProblemVars[i] << endl;
            }
        }
    }
    CurrentWeightedCSP->updateUb(Objective.getMin());
    UpperBound = Objective.getMin();
    for (int i = 0; i < ProblemSize; i++) {
        BestSol[i] = ProblemVars[i].getValue();
    }
    fail();
    return 0;
}
ILOCPGOALWRAPPER0(IloNewSolution, solver)
{
    return IlcNewSolution(solver);
}

ILCGOAL2(IlcInstantiateVar, IlcIntVar, var, IlcInt, varIndex)
{
    IlcInt value;
    int depth = Store::getDepth();

    if (var.isBound())
        return 0;

    // value ordering heuristic: try the unary support first
    value = CurrentWeightedCSP->getSupport(varIndex);
    if (!var.isInDomain(value)) {
        value = var.getMin();
    }
    return IlcOr(IlcGuess(getSolver(), var, value),
        IlcAnd(IlcRefute(getSolver(), var, value, depth), this));
}

// variable ordering heuristic: selects the first unassigned variable with the smallest ratio current domain size divided by actual current degree in the WCSP network
IlcInt IlcChooseMinSizeIntDivMaxDegree(const IlcIntVarArray vars)
{
    int varIndex = -1;
    ;
    double best = MAX_VAL - MIN_VAL;

    for (int i = 0; i < vars.getSize(); i++)
        if (!vars[i].isBound()) {
            // remove following "+1" when isolated variables are automatically assigned
            double heuristic = (double)CurrentWeightedCSP->getDomainSize(i) / (CurrentWeightedCSP->getDegree(i) + 1);
            if (varIndex < 0 || heuristic < best - 1. / 100001.) {
                best = heuristic;
                varIndex = i;
            }
        }
    return varIndex;
}

// For Queens problem
// variable ordering heuristic: among the unassigned variables with the smallest current domain size, selects the first one with the smallest domain value
IlcChooseIndex2(IlcChooseMinSizeMin, var.getSize(), var.getMin(), IlcIntVar)

    ILCGOAL1(IlcGenerateVars, IlcIntVarArray, vars)
{
    IlcInt index = IlcChooseMinSizeIntDivMaxDegree(vars);
    //  IlcInt index = IlcChooseMinSizeMin(vars);
    if (index == -1)
        return 0;
    return IlcAnd(IlcInstantiateVar(getSolver(), vars[index], index), this);
}
ILOCPGOALWRAPPER1(IloGenerateVars, solver, IloIntVarArray, vars)
{
    return IlcGenerateVars(solver, solver.getIntVarArray(vars));
}

void alldiff(IloEnv& env, IloModel& model, IloIntVarArray& vars)
{
    cout << "Add 1 hard AllDiff constraint on all the (permutation) problem variables." << endl;
    model.add(IloAllDiff(env, vars));
}

void zebra(IloEnv& env, IloModel& model, IloIntVarArray& vars)
{
    cout << "Add 5 hard AllDiff constraints for the Zebra problem." << endl;
    for (int i = 0; i < 5; i++) {
        int pos = i * 5;
        IloIntVarArray vars1(env, 5);
        vars1[0] = vars[0 + pos];
        vars1[1] = vars[1 + pos];
        vars1[2] = vars[2 + pos];
        vars1[3] = vars[3 + pos];
        vars1[4] = vars[4 + pos];
        model.add(IloAllDiff(env, vars1));
    }
}

void wqueens(IloEnv& env, IloModel& model, IloIntVarArray& vars)
{
    cout << "Add 3 hard AllDiff constraints for the Queens problem." << endl;

    model.add(vars); // ensure vars are the main decision variables

    int nqueen = vars.getSize();

    IloIntVarArray vars1(env, nqueen, -2 * nqueen, 2 * nqueen);
    IloIntVarArray vars2(env, nqueen, -2 * nqueen, 2 * nqueen);

    for (IloInt i = 0; i < nqueen; i++) {
        model.add(vars1[i] == vars[i] + i);
        model.add(vars2[i] == vars[i] - i);
    }

    model.add(IloAllDiff(env, vars));
    model.add(IloAllDiff(env, vars1));
    model.add(IloAllDiff(env, vars2));
}

void quasi(IloEnv& env, IloModel& model, IloIntVarArray& vars)
{
    int n = sqrt((double)vars.getSize());
    cout << "Add " << n * 2 << " hard AllDiff constraints for the \"homogeneous\" QuasiGroup problem." << endl;
    for (int i = 0; i < n; i++) {
        int pos = i * n;
        IloIntVarArray vars1(env, n);
        for (int j = 0; j < n; j++) {
            vars1[j] = vars[pos];
            pos++;
        }
        model.add(IloAllDiff(env, vars1));
    }
    for (int j = 0; j < n; j++) {
        int pos = j;
        IloIntVarArray vars1(env, n);
        for (int i = 0; i < n; i++) {
            vars1[i] = vars[pos];
            pos += n;
        }
        model.add(IloAllDiff(env, vars1));
    }
}

// Usage: iloglue problem_name.wcsp [verbosity]
int main(int argc, char** argv)
{
    string pbname;
    int nbvar, nbval, nbconstr;
    IloEnv env;
    IloTimer timer(env);

    if (argc >= 3)
        ToulBar2::verbose = atoi(argv[2]);

    try {
        IloModel model(env);

        // open the file
        ifstream file(argv[1]);
        if (!file) {
            cerr << "Could not open file " << argv[1] << endl;
            exit(EXIT_FAILURE);
        }

        // reads problem name and sizes
        file >> pbname;
        file >> nbvar;
        file >> nbval;
        file >> nbconstr;

        // creates the objective function
        IloIntVar obj(env, 0, MAX_COST, "objective");

        // creates the problem variables
        IloIntVarArray vars(env, nbvar, 0, MAX_DOMAIN_SIZE - 1);
        model.add(vars);
        for (int i = 0; i < nbvar; i++) {
            char* name = new char[16];
            sprintf(name, "x%d", i);
            vars[i].setName(name);
        }

        // creates a global weighted CSP constraint
        model.add(IloWeightedCSP(env, obj, vars, argv[1]));

        if (strstr(argv[1], "zebra"))
            zebra(env, model, vars);
        if (strstr(argv[1], "wqueens"))
            wqueens(env, model, vars);
        if (strstr(argv[1], "quasi"))
            quasi(env, model, vars);
        if (strstr(argv[1], "alldiff"))
            alldiff(env, model, vars);

        //     model.add(IloMinimize(env, obj)); DOES NOT WORK???

        IloSolver solver(model);
        IlogSolver = solver;

        // creates a goal to store the best solution found DOES NOT WORK???
        //     IloSolution solution(env);
        //     solution.add(obj);
        //     solution.add(vars);
        //     IloGoal storeSolution = IloStoreSolution(env, solution);

        Objective = solver.getIntVar(obj);
        ProblemSize = nbvar;
        ProblemVars = solver.getIntVarArray(vars);
        BestSol = new int[nbvar];
        BestSol[0] = -1;

        timer.start();

        //    IlcPrintTrace trace(solver, IlcTraceConstraint);
        //    solver.setTraceMode(IlcTrue);
        //    solver.setTrace(trace);

        // chooses high propagation level for AllDiff
        solver.setDefaultFilterLevel(IloAllDiffCt, IloExtendedLevel);

        // finds an optimal solution
        solver.solve(IloGenerateVars(env, vars) && IloNewSolution(env));
        Store::restore(0);
        // restores the best solution and shows it
        //     solver.solve(IloRestoreSolution(env,solution));
        //     cout << solver.getStatus() << " Solution" << endl;
        //     cout << "Optimum: " << solver.getValue(obj) << " in " << solver.getNumberOfFails() << " fails and " << solver.getTime() << " seconds." << endl;
        cout << "Optimum: " << UpperBound << " in " << solver.getNumberOfFails() << " fails and " << solver.getTime() << " seconds." << endl;
        if (ToulBar2::verbose >= 0 && BestSol[0] != -1) {
            cout << "Optimal solution:";
            for (int i = 0; i < nbvar; i++) {
                cout << " " << BestSol[i];
            }
            cout << endl;
        }
        solver.printInformation();
    } catch (IloException& ex) {
        cout << "Error: " << ex << endl;
    }
    env.end();
    return 0;
}

// Create the tuple set
// IloIntTupleSet(IloEnv env, const int arity);
// tuple.add(IloIntArray(env, 3, i, j, m));
// Create a table constraint = constraint in extension
// IloConstraint IloTableConstraint(const IloEnv env,
//                                  const IloNumVarArray vars,
//                                  const IloNumTupleSet set,
//                                  IloBool compatible);

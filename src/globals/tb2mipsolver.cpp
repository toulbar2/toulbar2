#include "tb2mipsolver.hpp"

#ifdef ILOGCPLEX

MIP* MIP::makeMIP()
{
    return new IlogMIP();
}

IlogMIP::IlogMIP()
{
    model = new IloModel(env);

    var = new IloNumVarArray(env);

    obj = new IloObjective(env, 0, IloObjective::Minimize);

    con = new IloRangeArray(env);

    cplex = new IloCplex(env);

    sols = new IloNumArray(env);

    buObjExpr = new IloNumExprArg();

    cols.clear();
    rowCount = 0;
    colCount = 0;
    objValue = 0;

    cplex->setOut(env.getNullStream());
    cplex->setWarning(env.getNullStream());
    cplex->setError(env.getNullStream());
    called = 0;
}

void IlogMIP::clear()
{

    cplex->clear();
    model->end();
    *model = IloModel(env);
    con->endElements();
    *con = IloRangeArray(env);
    obj->end();
    *obj = IloObjective(env, 0, IloObjective::Minimize);
    var->endElements();
    *var = IloNumVarArray(env);
    sols->clear();
    *sols = IloNumArray(env);

    cols.clear();
    rowCount = 0;
    colCount = 0;
    objValue = 0;
}

void IlogMIP::end()
{
    model->add(*var);
    model->add(*obj);
    model->add(*con);
    cplex->extract(*model);
}

void IlogMIP::addRows(int n)
{
    for (int i = 0; i < n; i++) {
        con->add(IloRange(env, -IloInfinity, IloInfinity));
        rowCount++;
    }
} // add n new inequalities to the linear program

void IlogMIP::addInt(int n)
{
    for (int i = 0; i < n; i++) {
        var->add(IloNumVar(env, 0.0, IloInfinity));
        cols.push_back(0);
        colCount++;
    }
} // add n integer variables to the linear program

void IlogMIP::addBool(int n)
{
    for (int i = 0; i < n; i++) {
        var->add(IloNumVar(env, 0.0, 0.0));
        cols.push_back(0);
        colCount++;
    }
} // add n boolean variables to the linear program

void IlogMIP::addCols(int n)
{
    for (int i = 0; i < n; i++) {
        var->add(IloNumVar(env, 0.0, IloInfinity));
        cols.push_back(0);
        colCount++;
    }
} // add n numeric variables to the linear program

void IlogMIP::rowBound(int n, int lower, int upper)
{

    (*con)[n].setBounds(lower, upper);
} // set the bounds of a variable

void IlogMIP::rowLowerBound(int n, int lower)
{
    (*con)[n].setLB(lower);
} // set the lower bound of a variable

void IlogMIP::rowUpperBound(int n, int upper)
{
    (*con)[n].setUB(upper);
} // set the upper bound of a variable

void IlogMIP::rowCoeff(int n, int count, int indexes[], double values[])
{
    for (int i = 0; i < count; i++) {
        (*con)[n].setLinearCoef((*var)[indexes[i]], values[i]);
    }
} // set the coefficients of the variables of a row

int IlogMIP::sol(int var1)
{

    if (ceil((*sols)[var1]) - (*sols)[var1] < 0.000001) {
        return ceil((*sols)[var1]);
    } else {
        return floor((*sols)[var1]);
    }

} // return the value of a variable (rounded down)

int IlogMIP::solValue()
{
    return objValue;
} // return the minimal of the current linear program

int IlogMIP::colUpperBound(int var1)
{
    return (*var)[var1].getUB();
} // return the lower bound of a value

void IlogMIP::colUpperBound(int var1, int i)
{
    (*var)[var1].setUB(i);
} // set the upper bound of a value

int IlogMIP::augment(int var1)
{
    (*var)[var1].setLB(1);
    solve();
    int cost = solValue();
    assert(sols->getSize() == 0 || sol(var1) == 1);
    (*var)[var1].setLB(0);
    return cost;
} // compute the minimal when a value is used

int IlogMIP::objCoeff(int var1)
{
    return cols[var1];
} // get the coefficient of a variable in the objective function

void IlogMIP::objCoeff(int var1, int i)
{
    cols[var1] = i;
    obj->setLinearCoef((*var)[var1], i);
} // set the coefficient of a variable in the objective function

int IlogMIP::solve()
{

    unsigned t0 = clock();
    cplex->solve();
    called += clock() - t0;
    if (cplex->getStatus() == IloAlgorithm::Infeasible) {
        throw Contradiction(); // cannot call THROWCONTRADICTION because no conflict function for weighted degree heuristic
    }
    if (cplex->getStatus() != IloAlgorithm::Optimal) {
        cerr << "Solution status = " << cplex->getStatus() << std::endl;
        cerr << "IlogMIP solver error." << endl;
        cerr << *con << endl;
        cerr << *obj << endl;
        cerr << *sols << endl;
        cerr << "cost = " << objValue << endl;
        throw InternalError();
    }

    if (cplex->getObjValue() - floor(cplex->getObjValue()) < 0.000001) {
        objValue = floor(cplex->getObjValue());
    } else {
        objValue = ceil(cplex->getObjValue());
    }

    if (objValue < 0) {
        objValue = 0;
    }

    cplex->getValues(*sols, *var);

    return objValue;
} // solve the current linear program for the minimal and store the values of the variables in the solution

int IlogMIP::sol(int varindex, int value)
{
    return sol(mapvar[varindex][value]);
} // check if the current solution is using this value

void IlogMIP::removeValue(int varindex, int value)
{
    colUpperBound(mapvar[varindex][value], 0);
} // remove the value from the domain of the variable

int IlogMIP::augment(int varindex, int value)
{
    if (sols->getSize() > 0 && sol(varindex, value) == 1) {
        return solValue();
    } else {
        return augment(mapvar[varindex][value]);
    }
} // return the minimum cost when this value is used

int IlogMIP::coeff(int varindex, int value)
{
    return objCoeff(mapvar[varindex][value]);
} // return the cost projected on this value

void IlogMIP::increaseCoeff(int varindex, int value, int newCoeff)
{
    int var1 = mapvar[varindex][value];
    cols[var1] += newCoeff;
    obj->setLinearCoef((*var)[var1], cols[var1]);
} // increase the cost projected on this value

void IlogMIP::backup()
{
    *buObjExpr = obj->getExpr();
} // backup the current solution (used before extensions)

int IlogMIP::restore()
{
    obj->setExpr(*buObjExpr);
    return solve();
} // restore the solution to the saved one

#endif

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

/** \file tb2mipsolver.hpp
 *  \brief Wrapper interface for MIP solver CPLEX
 */

#ifndef TB2MIPSOLVER_HPP_
#define TB2MIPSOLVER_HPP_

#include "core/tb2types.hpp"

class MIP { //Wrapper Interface for MIP solver

private:
    MIP* solver;

public:
    MIP();

    virtual ~MIP();

    virtual unsigned called_time()
    {
        if (solver)
            return solver->called_time();
        return 0;
    }

    virtual void clear()
    {
        if (solver)
            solver->clear();
    }

    virtual void end()
    {
        if (solver)
            solver->end();
    }

    virtual void addRows(int n)
    {
        if (solver)
            solver->addRows(n);
    }

    virtual void addInt(int n)
    {
        if (solver)
            solver->addInt(n);
    }

    virtual void addBool(int n)
    {
        if (solver)
            solver->addBool(n);
    }

    virtual void addCols(int n)
    {
        if (solver)
            solver->addCols(n);
    }

    virtual void rowBound(int n, int upper, int lower)
    {
        if (solver)
            solver->rowBound(n, upper, lower);
    }

    virtual void rowLowerBound(int n, int lower)
    {
        if (solver)
            solver->rowLowerBound(n, lower);
    }

    virtual void rowUpperBound(int n, int upper)
    {
        if (solver)
            solver->rowUpperBound(n, upper);
    }

    virtual void rowCoeff(int n, int count, int indexes[], double values[])
    {
        if (solver)
            solver->rowCoeff(n, count, indexes, values);
    }

    virtual int solValue()
    {
        if (solver)
            return solver->solValue();
        return 0;
    }

    virtual int sol(int var1)
    {
        if (solver)
            return solver->sol(var1);
        return 0;
    }

    virtual int colUpperBound(int var1)
    {
        if (solver)
            return solver->colUpperBound(var1);
        return 0;
    }

    virtual void colUpperBound(int var1, int i)
    {
        if (solver)
            return solver->colUpperBound(var1, i);
    }

    virtual int augment(int var1)
    {
        if (solver)
            return solver->augment(var1);
        return 0;
    }

    virtual int objCoeff(int var1)
    { //Get the coefficient of the variable from the MIP
        if (solver)
            return solver->objCoeff(var1);
        return 0;
    }

    virtual void objCoeff(int var1, int i)
    { // Set the coefficient of the variable
        if (solver)
            solver->objCoeff(var1, i);
    }

    virtual int solve()
    { //return the optimal value from the MIP
        if (solver)
            solver->solve();
        return 0;
    }

    virtual int sol(int varindex, int value)
    {
        if (solver)
            solver->sol(varindex, value);
        return 0;
    }

    virtual void removeValue(int varindex, int value)
    {
        if (solver)
            solver->removeValue(varindex, value);
    }

    virtual int augment(int varindex, int value)
    {
        if (solver)
            return solver->augment(varindex, value);
        return 0;
    }
    virtual int coeff(int varindex, int value)
    {
        if (solver)
            return solver->coeff(varindex, value);
        return 0;
    }

    virtual void increaseCoeff(int varindex, int value, int newCoeff)
    {
        if (solver)
            return solver->increaseCoeff(varindex, value, newCoeff);
    }

    virtual void getDomain(int varindex, vector<int>& domain)
    {
        if (solver)
            return solver->getDomain(varindex, domain);
    }

    virtual void backup()
    {
        if (solver)
            solver->backup();
    }

    virtual int restore()
    {
        if (solver)
            return solver->restore();
        return 0;
    }
};

#endif /*TB2MIPSOLVER_HPP_*/

#ifdef ILOGCPLEX
#ifndef TB2ILOGMIPSOLVER_HPP_
#define TB2ILOGMIPSOLVER_HPP_

#include "tb2types.hpp"
#include <ilcplex/ilocplex.h>

class IlogMIP : public MIP {
private:
    IloEnv env;
    IloModel* model;
    IloCplex* cplex;
    IloNumVarArray* var;
    IloObjective* obj;
    IloNumExprArg* buObjExpr;
    IloRangeArray* con;
    IloNumArray* sols;
    vector<int> cols;
    int rowCount;
    int colCount;
    int objValue;
    unsigned called;

public:
    map<Value, int>* mapvar;

    IlogMIP();

    virtual ~IlogMIP()
    {
        env.end();
    }

    unsigned called_time()
    {
        return called;
    }

    void clear();

    void end();

    void addRows(int n);

    void addInt(int n);

    void addBool(int n);

    void addCols(int n);

    void rowBound(int n, int upper, int lower);

    void rowLowerBound(int n, int lower);

    void rowUpperBound(int n, int upper);

    void rowCoeff(int n, int count, int indexes[], double values[]);

    int solValue();

    int sol(int var1);

    int colUpperBound(int var1);

    void colUpperBound(int var1, int i);

    int augment(int var1);

    int objCoeff(int var1);

    void objCoeff(int var1, int i);

    int solve();

    virtual int sol(int varindex, int value); // if the current solution using this value
    virtual void removeValue(int varindex, int value); // remove the value from the domain of the structure
    virtual int augment(int varindex, int value); // return the minimum cost when this value is used
    virtual int coeff(int varindex, int value); // return the projected cost on this value
    virtual void increaseCoeff(int varindex, int value, int newCoeff); // increase the projected cost on this value

    virtual void getDomain(int varindex, vector<int>& domain)
    {
        domain.clear();
        for (map<Value, int>::iterator v = mapvar[varindex].begin(); v != mapvar[varindex].end(); v++) {
            if (colUpperBound(v->second) == 1) {
                domain.push_back(v->first);
            }
        }
    } // return the corresponding domain in the linear program of a variable in the WCSP

    virtual void backup();
    virtual int restore();
};

#endif /*TB2ILOGMIPSOLVER_HPP_*/

#endif

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

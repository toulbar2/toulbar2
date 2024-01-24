/** \file tb2mipsolver.hpp
 *  \brief Wrapper interface for MIP solver CPLEX
 */

#ifndef TB2MIPSOLVER_HPP_
#define TB2MIPSOLVER_HPP_

#include "core/tb2types.hpp"

class MIP { // Wrapper Interface for MIP solver
public:
#ifdef ILOGCPLEX
    static MIP* makeMIP();
#endif
    virtual ~MIP()
    {
    }

    virtual unsigned called_time() = 0;

    virtual void clear() = 0;

    virtual void end() = 0;

    virtual void addRows(int n) = 0;

    virtual void addInt(int n) = 0;

    virtual void addBool(int n) = 0;

    virtual void addCols(int n) = 0;

    virtual void rowBound(int n, int upper, int lower) = 0;

    virtual void rowLowerBound(int n, int lower) = 0;

    virtual void rowUpperBound(int n, int upper) = 0;

    virtual void rowCoeff(int n, int count, int indexes[], double values[]) = 0;

    virtual int solValue() = 0;

    virtual int sol(int var1) = 0;

    virtual int colUpperBound(int var1) = 0;

    virtual void colUpperBound(int var1, int i) = 0;

    virtual int augment(int var1) = 0;

    virtual int objCoeff(int var1) = 0;

    virtual void objCoeff(int var1, int i) = 0;

    virtual int solve() = 0;

    virtual int sol(int varindex, int value) = 0;

    virtual void removeValue(int varindex, int value) = 0;

    virtual int augment(int varindex, int value) = 0;

    virtual int coeff(int varindex, int value) = 0;

    virtual void increaseCoeff(int varindex, int value, int newCoeff) = 0;

    virtual void getDomain(int varindex, vector<int>& domain) = 0;

    virtual void backup() = 0;

    virtual int restore() = 0;
};

#ifdef ILOGCPLEX
#include <ilcplex/ilocplex.h>

class IlogMIP FINAL : public MIP {
private:
    IloEnv env;
    IloModel* model;
    IloCplex* cplex;
    IloNumVarArray* var;
    IloObjective* obj;
    IloNumExprArg* buObjExpr;
    IloRangeArray* con;
    vector<int> cols;
    int rowCount;
    int colCount;
    int objValue;
    unsigned called;

public:
    map<Value, int>* mapvar;
    IloNumArray* sols;

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

#endif

#endif /*TB2MIPSOLVER_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

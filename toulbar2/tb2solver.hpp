/** \file tb2solver.hpp
 *  \brief Generic solver.
 * 
 */
 
#ifndef TB2SOLVER_HPP_
#define TB2SOLVER_HPP_

#include "tb2variable.hpp"
#include "tb2wcsp.hpp"

class Solver
{
    Store store;
    long long nbNodes;
    long long nbBacktracks;
    vector<Variable *> vars;
    VariableList unassignedVars;
    
    // Monocriterion: only one WCSP for the moment
    Cost upperBound;
    Variable objective;
    WCSP wcsp;

    // Heuristics and search methods
    Variable *getVarMinDomainDivMaxDegree();
    Variable *getNextUnassignedVar();
    void binaryChoicePoint(Variable *x, Value value);
    void narySortedChoicePoint(Variable *x);
    void naryChoicePoint(Variable *x);
    void recursiveSolve();
    void whenContradiction();
    
public:
    Solver(int storeSize, Cost initUpperBound);
    
    ~Solver();
    
    Store *getStore() {return &store;}
    vector<Variable *> *getVars() {return &vars;}
    VariableList *getUnassignedVars() {return &unassignedVars;}
    Variable *getObjective() {return &objective;}
    
    void read_wcsp(const char *fileName);
    
    bool solve();
};

#endif /*TB2SOLVER_HPP_*/

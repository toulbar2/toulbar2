/** \file tb2solver.hpp
 *  \brief Generic solver.
 * 
 */
 
#ifndef TB2SOLVER_HPP_
#define TB2SOLVER_HPP_

#include "toulbar2.hpp"
template <class T> struct DLink;
template <class T> class BTList;

class Solver
{
    static Solver *currentSolver;

    Store *store;
    Long nbNodes;
    Long nbBacktracks;
    WeightedCSP *wcsp;
    DLink<Value> *allVars;
    BTList<Value> *unassignedVars;
    int lastConflictVar;
        
    // Heuristics and search methods
	int getVarMinDomainDivMaxWeightedDegreeLastConflict();
	int getVarMinDomainDivMaxWeightedDegree();
    int getVarMinDomainDivMaxDegreeLastConflict();
    int getVarMinDomainDivMaxDegree();
    int getNextUnassignedVar();
    int getVarVACHeuristic();
    void increase(int varIndex, Value value);
    void decrease(int varIndex, Value value);
    void assign(int varIndex, Value value);
    void remove(int varIndex, Value value);
   
    void singletonConsistency();
   
    void binaryChoicePoint(int xIndex, Value value);
    void binaryChoicePointLDS(int xIndex, Value value, int discrepancy);
    void narySortedChoicePoint(int xIndex);
    void narySortedChoicePointLDS(int xIndex, int discrepancy);
    void newSolution();
    void recursiveSolve();
    void recursiveSolveLDS(int discrepancy);
    
public:
    Solver(int storeSize, Cost initUpperBound);
    ~Solver();
    
    void read_wcsp(const char *fileName);
    void read_random(int n, int m, vector<int>& p, int seed);
    
    bool solve();
    
    void dump_wcsp(const char *fileName);
    void read_solution(const char *fileName);
    
    friend void setvalue(int wcspId, int varIndex, Value value);

    WeightedCSP* getWCSP() { return wcsp; }
};

#endif /*TB2SOLVER_HPP_*/

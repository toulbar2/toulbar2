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
    Long nbBacktracksLimit;
    WeightedCSP *wcsp;
    DLink<Value> *allVars;
    BTList<Value> *unassignedVars;
    int lastConflictVar;

    BigInteger nbSol;
    Long nbSGoods;				//number of #good which created
    Long nbSGoodsUse;			//number of #good which used
    map<int,BigInteger > ubSol;	// upper bound of solution number
    double timeDeconnect;		// time for the disconnection

    // Heuristics and search methods
    /// \warning hidden feature: do not branch on variable indexes from ToulBar2::nbDecisionVars to the last variable
    void initVarHeuristic();
    int getVarMinDomainDivMaxWeightedDegreeLastConflictRandomized();
	int getVarMinDomainDivMaxWeightedDegreeLastConflict();
	int getVarMinDomainDivMaxWeightedDegreeRandomized();
	int getVarMinDomainDivMaxWeightedDegree();
    int getVarMinDomainDivMaxDegreeLastConflictRandomized();
    int getVarMinDomainDivMaxDegreeLastConflict();
    int getVarMinDomainDivMaxDegreeRandomized();
    int getVarMinDomainDivMaxDegree();
    int getNextUnassignedVar();
    int getMostUrgent();
    void increase(int varIndex, Value value);
    void decrease(int varIndex, Value value);
    void assign(int varIndex, Value value);
    void remove(int varIndex, Value value);
    void conflict() {}
    void enforceUb();
    void singletonConsistency();

    void binaryChoicePoint(int xIndex, Value value);
    void binaryChoicePointLDS(int xIndex, Value value, int discrepancy);
    void narySortedChoicePoint(int xIndex);
    void narySortedChoicePointLDS(int xIndex, int discrepancy);
    void newSolution();
    void recursiveSolve();
    void recursiveSolveLDS(int discrepancy);
    Value postponeRule(int varIndex);
    void scheduleOrPostpone(int varIndex);

    int getVarMinDomainDivMaxWeightedDegreeLastConflictRandomized(Cluster *cluster);
	int getVarMinDomainDivMaxWeightedDegreeLastConflict(Cluster *cluster);
	int getVarMinDomainDivMaxWeightedDegreeRandomized(Cluster *cluster);
	int getVarMinDomainDivMaxWeightedDegree(Cluster *cluster);
    int getVarMinDomainDivMaxDegreeLastConflictRandomized(Cluster *cluster);
    int getVarMinDomainDivMaxDegreeLastConflict(Cluster *cluster);
    int getVarMinDomainDivMaxDegreeRandomized(Cluster *cluster);
    int getVarMinDomainDivMaxDegree(Cluster *cluster);
    int getNextUnassignedVar(Cluster *cluster);

    Cost binaryChoicePoint(Cluster *cluster, Cost lbgood, Cost cub, int varIndex, Value value);
    Cost recursiveSolve(Cluster *cluster, Cost lbgood, Cost cub);
    void russianDollSearch(Cluster *c, Cost cub);

    BigInteger binaryChoicePointSBTD(Cluster *cluster, int varIndex, Value value);
    BigInteger sharpBTD(Cluster *cluster);
    void approximate(BigInteger& nbsol, TreeDecomposition* td);

public:
    Solver(int storeSize, Cost initUpperBound);
    ~Solver();

    void read_wcsp(const char *fileName);
    void read_random(int n, int m, vector<int>& p, int seed, bool forceSubModular = false );

    bool solve();

    bool solve_symmax2sat(int n, int m, int *posx, int *posy, double *cost, int *sol);

    void dump_wcsp(const char *fileName, bool original = true);
    void read_solution(const char *fileName);
    void parse_solution(const char *certificate);

    friend void setvalue(int wcspId, int varIndex, Value value);

    WeightedCSP* getWCSP() { return wcsp; }
};

class NbBacktracksOut
{
public:
    NbBacktracksOut() {if (ToulBar2::verbose >= 2) cout << "... limit on the number of backtracks reached!" << endl;}
};

int solveSymMax2SAT(int n, int m, int *posx, int *posy, double *cost, int *sol);
extern "C" int solvesymmax2sat_(int *n, int *m, int *posx, int *posy, double *cost, int *sol);

#endif /*TB2SOLVER_HPP_*/

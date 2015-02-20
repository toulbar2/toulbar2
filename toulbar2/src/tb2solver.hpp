/** \file tb2solver.hpp
 *  \brief Generic solver.
 *
 */

#ifndef TB2SOLVER_HPP_
#define TB2SOLVER_HPP_

#include "toulbar2lib.hpp"
#include "tb2store.hpp"

template <class T> struct DLink;
template <class T> class BTList;

const double epsilon = 1e-6; // 1./100001.

const int MAX_BRANCH_SIZE = 1000000;
const int CHOICE_POINT_LIMIT = INT_MAX - 2 * MAX_BRANCH_SIZE;
const int OPEN_NODE_LIMIT = INT_MAX;

class Solver : public WeightedCSPSolver
{
public:
    class OpenNode
    {
    public:
        Cost cost;      // lower bound associated to the open node
        int first;      // first position in the list of choice points corresponding to a branch in order to reconstruct the open node
        int last;       // last position (excluded) in the list of choice points corresponding to a branch in order to reconstruct the open node

        OpenNode(Cost cost_, int first_, int last_) : cost(cost_), first(first_), last(last_) {}
        bool operator<(const OpenNode& right) const {return (cost > right.cost) || (cost == right.cost && ((last-first) < (right.last-right.first)));} // reverse order to get the open node with the smallest lower bound first and deepest depth next
    };

    class OpenList : public priority_queue<OpenNode>
    {
    public:
        Cost glb;   // current cluster lower bound built from closed nodes
        Cost cub;   // current cluster upper bound

        OpenList(Cost lb, Cost ub) : glb(lb), cub(ub) {}
        OpenList() : glb(MAX_COST), cub(MAX_COST) {}

        void clear() {glb=MAX_COST; cub=MAX_COST; while (!empty()) pop();}
    };

    typedef enum {
        CP_ASSIGN = 0, CP_REMOVE = 1, CP_INCREASE = 2, CP_DECREASE = 3, CP_REMOVE_RANGE = 4, CP_MAX
    } ChoicePointOp;
    static const string CPOperation[CP_MAX]; // for pretty print

    struct ChoicePoint {
        ChoicePointOp op;   // choice point operation
        int varIndex;       // variable wcsp's index
        Value value;        // variable's value
        bool reverse;       // true if the choice point corresponds to the last right branch of an open node

        ChoicePoint(ChoicePointOp op_, int var_, Value val_, bool rev_) : op(op_), varIndex(var_), value(val_), reverse(rev_) {}
    };

    class CPStore : public vector<ChoicePoint>
    {
        public:
            int start;       // beginning of the current branch
            int stop;        // deepest saved branch end (should be free at this position)
            StoreInt index;  // current branch depth (should be free at this position)

            CPStore(Store *s) : start(0), stop(0), index(0, &s->storeInt) {}

            void addChoicePoint(ChoicePointOp op, int varIndex, Value value, bool reverse);
            void store() {start = stop; index = start;}
    };

    void addChoicePoint(ChoicePointOp op, int varIndex, Value value, bool reverse);
    void addOpenNode(CPStore &cp, OpenList &open, Cost lb);
    void restore(CPStore &cp, OpenNode node);

protected:
    Store *store;
    Long nbNodes;
    Long nbBacktracks;
    Long nbBacktracksLimit;
    WeightedCSP *wcsp;
    DLink<Value> *allVars;
    BTList<Value> *unassignedVars;
    int lastConflictVar;
    void *searchSize;

    BigInteger nbSol;
    Long nbSGoods;				//number of #good which created
    Long nbSGoodsUse;			//number of #good which used
    map<int,BigInteger > ubSol;	// upper bound of solution number
    double timeDeconnect;		// time for the disconnection

    CPStore *cp;                // choice point cache for open nodes (except BTD)
    OpenList *open;             // list of open nodes (except BTD)
    Long hybridBFSLimit;        // limit on number of backtracks for hybrid search (except BTD)
    Long nbHybrid;
    Long nbHybridContinue;

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
    void increase(int varIndex, Value value, bool reverse = false);
    void decrease(int varIndex, Value value, bool reverse = false);
    void assign(int varIndex, Value value, bool reverse = false);
    void remove(int varIndex, Value value, bool reverse = false);
    void remove(int varIndex, ValueCost *array, int first, int last, bool reverse = false);
    void conflict() {}
    void enforceUb();
    void singletonConsistency();

    void binaryChoicePoint(int xIndex, Value value, Cost lb = MIN_COST);
    void binaryChoicePointLDS(int xIndex, Value value, int discrepancy);
    void narySortedChoicePoint(int xIndex, Cost lb = MIN_COST);
    void narySortedChoicePointLDS(int xIndex, int discrepancy);
    void newSolution();
    void recursiveSolve(Cost lb = MIN_COST);
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

    pair<Cost, Cost> binaryChoicePoint(Cluster *cluster, Cost lbgood, Cost cub, int varIndex, Value value);
    pair<Cost, Cost> recursiveSolve(Cluster *cluster, Cost lbgood, Cost cub);
    pair<Cost,Cost> hybridSolve(Cluster *root, Cost clb, Cost cub);
    pair<Cost,Cost> hybridSolve() {return hybridSolve(NULL,  wcsp->getLb(), wcsp->getUb());}
    void russianDollSearch(Cluster *c, Cost cub);

    BigInteger binaryChoicePointSBTD(Cluster *cluster, int varIndex, Value value);
    BigInteger sharpBTD(Cluster *cluster);
    void approximate(BigInteger& nbsol, TreeDecomposition* td);

public:
    Solver(int storeSize, Cost initUpperBound);
    ~Solver();

    void read_wcsp(const char *fileName);
    void read_random(int n, int m, vector<int>& p, int seed, bool forceSubModular = false );

    Long getNbNodes() const {return nbNodes;}
    Long getNbBacktracks() const {return nbBacktracks;}
    set<int> getUnassignedVars() const;

    bool solve();

    Cost narycsp(string cmd, vector<Value> &solution);

    bool solve_symmax2sat(int n, int m, int *posx, int *posy, double *cost, int *sol);

    void dump_wcsp(const char *fileName, bool original = true);
    void read_solution(const char *fileName);
    void parse_solution(const char *certificate);

    Cost getSolution(vector<Value>& solution);

    friend void setvalue(int wcspId, int varIndex, Value value, void *solver);

    WeightedCSP* getWCSP() { return wcsp; }
};

class NbBacktracksOut
{
public:
    NbBacktracksOut() {if (ToulBar2::verbose >= 2) cout << "... limit on the number of backtracks reached!" << endl;}
};

class TimeOut
{
public:
    TimeOut() {if (ToulBar2::verbose >= 2) cout << "... time limit reached!" << endl;}
};

int solveSymMax2SAT(int n, int m, int *posx, int *posy, double *cost, int *sol);
extern "C" int solvesymmax2sat_(int *n, int *m, int *posx, int *posy, double *cost, int *sol);

#endif /*TB2SOLVER_HPP_*/

/** \file tb2solver.hpp
 *  \brief Generic solver.
 *
 */

#ifndef TB2SOLVER_HPP_
#define TB2SOLVER_HPP_

#include "toulbar2lib.hpp"
#include "utils/tb2store.hpp"

template <class T>
class DLink;
template <class T>
class BTList;

class NeighborhoodStructure;
class RandomNeighborhoodChoice;
class ClustersNeighborhoodStructure;
class RandomClusterChoice;
class ParallelRandomClusterChoice;

const double epsilon = 1e-6; // 1./100001.

class Solver : public WeightedCSPSolver {
public:
    class OpenNode {
        Cost cost; // global lower bound associated to the open node
    public:
        ptrdiff_t first; // first position in the list of choice points corresponding to a branch in order to reconstruct the open node
        ptrdiff_t last; // last position (excluded) in the list of choice points corresponding to a branch in order to reconstruct the open node

        OpenNode(Cost cost_, ptrdiff_t first_, ptrdiff_t last_)
            : cost(cost_)
            , first(first_)
            , last(last_)
        {
        }
        bool operator<(const OpenNode& right) const { return (cost > right.cost) || (cost == right.cost && ((last - first) < (right.last - right.first) || ((last - first) == (right.last - right.first) && last >= right.last))); } // reverse order to get the open node with first, the smallest lower bound, and next, the deepest depth, and next, the oldest time-stamp

        Cost getCost(Cost delta = MIN_COST) const { return MAX(MIN_COST, cost - delta); }
    };

    class CPStore;
    class OpenList FINAL : public priority_queue<OpenNode> {
        Cost clb; // current cluster lower bound built from closed nodes (independent of any soft arc consistency cost moves)
        Cost cub; // current cluster upper bound (independent of any soft arc consistency cost moves)
    public:
        OpenList(Cost lb, Cost ub)
            : clb(lb)
            , cub(ub)
        {
        }
        OpenList()
            : clb(MAX_COST)
            , cub(MAX_COST)
        {
        } /// \warning use also this method to clear an open list

        bool finished() const
        {
            assert(clb <= cub);
            return (empty() || CUT(top().getCost(), clb));
        }
        Cost getLb(Cost delta = MIN_COST) const { return MIN(MAX(MIN_COST, clb - delta), (empty() ? MAX_COST : top().getCost(delta))); }

        Cost getClosedNodesLb(Cost delta = MIN_COST) const { return MAX(MIN_COST, clb - delta); }
        void setClosedNodesLb(Cost lb, Cost delta = MIN_COST)
        {
            clb = MAX(MIN_COST, lb + delta);
            assert(clb <= cub);
        }
        void updateClosedNodesLb(Cost lb, Cost delta = MIN_COST) { clb = MIN(clb, MAX(MIN_COST, lb + delta)); }

        Cost getUb(Cost delta = MIN_COST) const { return MAX(MIN_COST, cub - delta); }
        void setUb(Cost ub, Cost delta = MIN_COST) { cub = MAX(MIN_COST, ub + delta); }
        void updateUb(Cost ub, Cost delta = MIN_COST)
        {
            Cost tmpub = MAX(MIN_COST, ub + delta);
            cub = MIN(cub, tmpub);
            clb = MIN(clb, tmpub);
        }

        size_type capacity() const { return c.capacity(); }
        priority_queue::container_type::iterator begin() { return c.begin(); }
        priority_queue::container_type::iterator end() { return c.end(); }
    };

    class SolutionTrie {
    public:
        class TrieNode {
        public:
            TrieNode(size_t w = 0);
            ~TrieNode();
            vector<vector<TrieNode*>> insertSolution(const vector<Value>& sol, unsigned int pos, vector<vector<TrieNode*>> nodesAtPos);
            vector<TrieNode*> sons;
            vector<vector<TrieNode*>> insertNode(Value v, unsigned int pos, vector<vector<TrieNode*>> nodesAtPos);
            bool present(Value v);
            void printTrie(vector<Value>& sol);
            static size_t nbSolutions;
            static vector<size_t> widths;
        };

        SolutionTrie(){};
        ~SolutionTrie(){};
        void init(const vector<Variable*>& vv);
        void insertSolution(const vector<Value>& sol);
        void printTrie();
        size_t getNbSolutions() { return root.nbSolutions; };
        vector<vector<TrieNode*>> getNodesAtPos() { return nodesAtPos; };

    private:
        TrieNode root;
        vector<vector<TrieNode*>> nodesAtPos;
    };

    Mdd computeMDD(SolutionTrie* solTrie, Cost cost);
    ostream& printLayers(ostream& os, Mdd mdd);

    typedef enum {
        CP_ASSIGN = 0,
        CP_REMOVE = 1,
        CP_INCREASE = 2,
        CP_DECREASE = 3,
        CP_REMOVE_RANGE = 4,
        CP_MAX
    } ChoicePointOp;
    static const string CPOperation[CP_MAX]; // for pretty print

    struct ChoicePoint {
        ChoicePointOp op; // choice point operation
        int varIndex; // variable wcsp's index
        Value value; // variable's value
        bool reverse; // true if the choice point corresponds to the last right branch of an open node

        ChoicePoint(ChoicePointOp op_, int var_, Value val_, bool rev_)
            : op(op_)
            , varIndex(var_)
            , value(val_)
            , reverse(rev_)
        {
        }
    };

    class CPStore FINAL : public vector<ChoicePoint> {
    public:
        ptrdiff_t start; // beginning of the current branch
        ptrdiff_t stop; // deepest saved branch end (should be free at this position)
        StoreCost index; // current branch depth (should be free at this position)

        CPStore()
            : start(0)
            , stop(0)
            , index(0)
        {
        }

        void addChoicePoint(ChoicePointOp op, int varIndex, Value value, bool reverse);
        void store()
        {
            start = stop;
            index = start;
        }
    };

    void addChoicePoint(ChoicePointOp op, int varIndex, Value value, bool reverse);
    void addOpenNode(CPStore& cp, OpenList& open, Cost lb, Cost delta = MIN_COST); ///< \param delta cost moved out from the cluster by soft arc consistency
    void restore(CPStore& cp, OpenNode node);

protected:
    friend class NeighborhoodStructure;
    friend class RandomNeighborhoodChoice;
    friend class ClustersNeighborhoodStructure;
    friend class RandomClusterChoice;
    friend class ParallelRandomClusterChoice;
    friend class VACExtension;

    Long nbNodes;
    Long nbBacktracks;
    Long nbBacktracksLimit;
    WeightedCSP* wcsp;
    DLink<Value>* allVars;
    BTList<Value>* unassignedVars;
    int lastConflictVar;
    void* searchSize;

    BigInteger nbSol;
    Long nbSGoods; //number of #good which created
    Long nbSGoodsUse; //number of #good which used
    map<int, BigInteger> ubSol; // upper bound of solution number
    double timeDeconnect; // time for the disconnection
    int tailleSep;

    CPStore* cp; // choice point cache for open nodes (except BTD)
    OpenList* open; // list of open nodes (except BTD)
    Long hbfsLimit; // limit on number of backtracks for hybrid search (except BTD)
    Long nbHybrid;
    Long nbHybridContinue;
    Long nbHybridNew;
    Long nbRecomputationNodes;

    //only for pretty print of optimality gap information
    Cost initialLowerBound;
    Cost globalLowerBound;
    Cost globalUpperBound;
    int initialDepth;
    void initGap(Cost newlb, Cost newub);
    void showGap(Cost newlb, Cost newub);

    Cost prevDivSolutionCost;
    SolutionTrie solTrie;
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
    void remove(int varIndex, ValueCost* array, int first, int last, bool reverse = false);
    void conflict() {}
    void enforceUb();
    void singletonConsistency();
    void binaryChoicePoint(int xIndex, Value value, Cost lb = MIN_COST);
    void binaryChoicePointLDS(int xIndex, Value value, int discrepancy);
    void narySortedChoicePoint(int xIndex, Cost lb = MIN_COST);
    void narySortedChoicePointLDS(int xIndex, int discrepancy);
    Value postponeRule(int varIndex);
    void scheduleOrPostpone(int varIndex);

    int getVarMinDomainDivMaxWeightedDegreeLastConflictRandomized(Cluster* cluster);
    int getVarMinDomainDivMaxWeightedDegreeLastConflict(Cluster* cluster);
    int getVarMinDomainDivMaxWeightedDegreeRandomized(Cluster* cluster);
    int getVarMinDomainDivMaxWeightedDegree(Cluster* cluster);
    int getVarMinDomainDivMaxDegreeLastConflictRandomized(Cluster* cluster);
    int getVarMinDomainDivMaxDegreeLastConflict(Cluster* cluster);
    int getVarMinDomainDivMaxDegreeRandomized(Cluster* cluster);
    int getVarMinDomainDivMaxDegree(Cluster* cluster);
    int getNextUnassignedVar(Cluster* cluster);

    pair<Cost, Cost> binaryChoicePoint(Cluster* cluster, Cost lbgood, Cost cub, int varIndex, Value value);
    pair<Cost, Cost> recursiveSolve(Cluster* cluster, Cost lbgood, Cost cub);
    pair<Cost, Cost> hybridSolve(Cluster* root, Cost clb, Cost cub);
    pair<Cost, Cost> russianDollSearch(Cluster* c, Cost cub);

    BigInteger binaryChoicePointSBTD(Cluster* cluster, int varIndex, Value value);
    BigInteger sharpBTD(Cluster* cluster);
    void approximate(BigInteger& nbsol, TreeDecomposition* td);

public:
    Solver(Cost initUpperBound);
    ~Solver();

    Cost read_wcsp(const char* fileName);
    void read_random(int n, int m, vector<int>& p, int seed, bool forceSubModular = false, string globalname = "");

    Long getNbNodes() const FINAL { return nbNodes; }
    Long getNbBacktracks() const FINAL { return nbBacktracks; }
    set<int> getUnassignedVars() const;
    unsigned int numberOfUnassignedVariables() const; // faster than its WCSP linear-time counterpart, but it is valid only during search

    virtual bool solve(bool first = true);

    // internal methods called by solve, for advanced programmers only!!!
    void beginSolve(Cost ub);
    Cost preprocessing(Cost ub);
    void recursiveSolve(Cost lb = MIN_COST);
    void recursiveSolveLDS(int discrepancy);
    pair<Cost, Cost> hybridSolve() { return hybridSolve(NULL, wcsp->getLb(), wcsp->getUb()); }
    void endSolve(bool isSolution, Cost cost, bool isComplete);
    // end of internal solve methods

    Cost narycsp(string cmd, vector<Value>& solution);

    bool solve_symmax2sat(int n, int m, int* posx, int* posy, double* cost, int* sol);

    void dump_wcsp(const char* fileName, bool original = true, ProblemFormat format = WCSP_FORMAT);
    void read_solution(const char* fileName, bool updateValueHeuristic = true);
    void parse_solution(const char* certificate, bool updateValueHeuristic = true);

    virtual void newSolution();
    const vector<Value> getSolution() { return wcsp->getSolution(); }
    Double getSolutionValue() const { return wcsp->getSolutionValue(); }
    Cost getSolutionCost() const { return wcsp->getSolutionCost(); }
    Cost getSolution(vector<Value>& solution) const
    {
        Cost cost = MAX_COST;
        solution = wcsp->getSolution(&cost);
        return cost;
    }
    vector<pair<Double, vector<Value>>> getSolutions() const { return wcsp->getSolutions(); }

    friend void setvalue(int wcspId, int varIndex, Value value, void* solver);

    WeightedCSP* getWCSP() FINAL { return wcsp; }
};

class SolverOut : public std::exception {
public:
    SolverOut()
    {
        ToulBar2::limited = true;
        if (ToulBar2::verbose >= 2)
            cout << SolverOut::what() << endl;
    }
    virtual const char* what() const throw() { return "... some solver limit was reached!"; }
};

class NbBacktracksOut : public SolverOut {
public:
    const char* what() const throw() FINAL { return "... limit on the number of backtracks reached!"; }
};

class NbSolutionsOut : public SolverOut {
public:
    const char* what() const throw() FINAL { return "... limit on the number of solutions reached!"; }
};

class DivSolutionOut : public SolverOut {
public:
    const char* what() const throw() FINAL { return "... optimal diverse solution found!"; }
};

class TimeOut : public SolverOut {
public:
    const char* what() const throw() FINAL { return "... time limit reached!"; }
};

int solveSymMax2SAT(int n, int m, int* posx, int* posy, double* cost, int* sol);
extern "C" int solvesymmax2sat_(int* n, int* m, int* posx, int* posy, double* cost, int* sol);

#endif /*TB2SOLVER_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

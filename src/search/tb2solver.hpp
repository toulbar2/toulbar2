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
    private:
        friend class serialization::access;
        template <class Archive>
        void serialize(Archive& ar, const unsigned int version)
        {
            ar& cost; // node lower bound
            ar& first; // pointer of type intptr_t = ptrdiff_t based on signed integer type
            ar& last; // means the "last" choice point in CPStore = vector<ChoicePoint> is at the adr (last-1)
        }

        Cost cost; // global lower bound associated to the open node
    public:
        ptrdiff_t first; // first position in the list of choice points corresponding to a branch in order to reconstruct the open node
        ptrdiff_t last; // last position (excluded) in the list of choice points corresponding to a branch in order to reconstruct the open node

        OpenNode() {} // default constructor added to avoid boost/serialization/access.hpp:130:9: error
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
    class OpenList FINAL : public std::priority_queue<OpenNode> {
    private:
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
        void init()
        {
            clb = MAX_COST;
            cub = MAX_COST;
        }

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
    private:
        friend class serialization::access;
        template <class Archive>
        void serialize(Archive& ar, const unsigned int version)
        {
            ar& op;
            ar& varIndex;
            ar& value;
            ar& reverse;
        }
    public:
        ChoicePointOp op; // choice point operation
        int varIndex; // variable wcsp's index
        Value value; // variable's value
        bool reverse; // true if the choice point corresponds to the last right branch of an open node

        ChoicePoint() {} // default constructor added to avoid boost/serialization/access.hpp:130:9: error
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

#ifdef OPENMPI
    /**
     * \brief class to send/receive work to/from workers in the form of an object i.e. a message in MPI's semantic //
     *
     */
    class Work {
    private:
        friend class serialization::access;
        template <class Archive>
        void serialize(Archive& ar, const unsigned int version)
        {
            ar& hbfs;
            ar& nbNodes;
            ar& nbBacktracks;
            ar& nbDEE;
            ar& lb;
            ar& ub;
            ar& sol;
            ar& open;
            ar& cp;
        }

    public:
        bool hbfs = true; // false if master cp/open memory is out
        Long nbNodes = 0;
        Long nbBacktracks = 0;
        Long nbDEE = 0;
        Cost lb; // best lower bound known by the master or found by the worker
        Cost ub; // best current solution a.k.a incumbent solution
        vector<Value> sol;
        vector<OpenNode> open; // priority queue which will contain the node(s) to eXchange between the processes
        vector<ChoicePoint> cp; // vector of choice points

        Work() {} // dummy message when stopping

        /**
         * @brief constructor used by the master
         */
        Work(const CPStore& cpMaster_, OpenList& openMaster_, const Cost lbMaster_, const Cost ubMaster_, vector<Value>& sol_)
        : lb(lbMaster_)
        , ub(ubMaster_)
        {

            open.push_back(openMaster_.top());

            openMaster_.pop(); // pop up directly the open queue!!

            for (ptrdiff_t i = open[0].first; i < open[0].last; i++) { // create a sequence of decisions in the vector vec.  node.last: index in CPStore of the past-the-end element
                cp.push_back(cpMaster_[i]);
            }

            sol.swap(sol_); //after the swap sol_ in the worker is an empty vector
        }

        /**
         * @brief constructor used by the worker
         */
        Work(CPStore& cpWorker_, OpenList& openWorker_, Long nbNodes_, Long nbBacktracks_, Long nbDEE_, const Cost lbWorker_, const Cost ubWorker_, vector<Value>& sol_)
        : nbNodes(nbNodes_)
        , nbBacktracks(nbBacktracks_)
        , nbDEE(nbDEE_)
        , lb(lbWorker_)
        , ub(ubWorker_)
        {

            while (!openWorker_.empty()) // init of vector of OpenNode nodeX
            {
                OpenNode node = openWorker_.top();

                open.push_back(node);

                openWorker_.pop(); // pop up directly the open queue!!
            }
            openWorker_.init(); // clb=cub=MAX_COST  method added to init openList attributes
            // because emptying an openList is not enough we have to initialize its attributes

            if (open.size() > 0) {
                for (ptrdiff_t i = 0; i < cpWorker_.stop; i++) {
                    cp.push_back(cpWorker_[i]);
                }
            }

            cpWorker_.clear(); // size = 0  added to put new cp out of the while(1)
//            cpWorker_.shrink_to_fit(); // to win space we shrink the vector: capacity=0 //TODO: test if it speeds-up things or not

            sol.swap(sol_); //after the swap sol_ in the worker is an empty vector
        }

        friend ostream& operator<<(ostream& os, Work& work)
        {
            os << "[lb: " << work.lb << ", ub: " << work.ub << ", nbvars: " << work.sol.size() << ", nodes: " << work.open.size() << ", choices: " << work.cp.size() << "]";
            return os;
        }
    };
#endif

    void addChoicePoint(ChoicePointOp op, int varIndex, Value value, bool reverse);
    void addOpenNode(CPStore& cp, OpenList& open, Cost lb, Cost delta = MIN_COST); ///< \param delta cost moved out from the cluster by soft arc consistency
    void restore(CPStore& cp, OpenNode node);
    char opSymbol(const CPStore& cp, const ptrdiff_t idx, OpenNode nd);

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
    vector<DLink<Value>*> allVars;
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

    void Manage_Freedom(Cluster* cluster);
    double nbChoices;
    double nbForcedChoices;
    double nbForcedChoiceChange;
    double nbChoiceChange;
    double nbReadOnly;
    int solveDepth;

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

#ifdef OPENMPI
    mpi::communicator world;
    queue<int> idleQ; //MASTER ONLY container with the rank of the free workers
    std::unordered_map<int, Cost> activeWork; //MASTER ONLY map the rank i of a worker with the lb cost of an open node
    std::unordered_map<int, Cost> bestsolWork; //MASTER ONLY map the rank i of a worker with the cost of the best solution sent by the master

    inline bool MPI_interrupted()
    {
        if (world.iprobe(mpi::any_source, DIETAG)) return true;
        else return false;
    }
    pair<Cost, Cost> hybridSolveMaster(Cluster* root, Cost clb, Cost cub);
    pair<Cost, Cost> hybridSolveWorker(Cluster* root, Cost clb, Cost cub);
#endif

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
    unsigned int numberOfUnassignedVariables() const; // faster than its WCSP linear-time counterpart, but it is valid only during search (otherwise returns -1)
    void updateVarHeuristic(); /// \brief to be called if DAC order has been changed after preprocessing (initVarHeuristic call)

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

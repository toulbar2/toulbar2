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

typedef int (Solver::*intFunctionCall_t)();

class Solver : public WeightedCSPSolver {
public:
    static Solver* CurrentSolver; // Current solver used by open node heuristics

    class OpenNode {
    private:
#ifdef OPENMPI
        friend class serialization::access;
        template <class Archive>
        void serialize(Archive& ar, const unsigned int version)
        {
            ar& cost; // node lower bound
            ar& first; // pointer of type intptr_t = ptrdiff_t based on signed integer type
            ar& last; // means the "last" choice point in CPStore = vector<ChoicePoint> is at the adr (last-1)
        }
#endif
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
        bool operator<(const OpenNode& right) const { int res = 0; if (ToulBar2::sortBFS) return ((cost > right.cost) || (cost == right.cost && (res = Solver::recHeuristicCmp(*this, first, right, right.first)) == 2) || (cost == right.cost && res == 1 && last >= right.last)); else return (cost > right.cost) || (cost == right.cost && ((last - first) < (right.last - right.first) || ((last - first) == (right.last - right.first) && last >= right.last))); } // reverse order to get the open node with first, the smallest lower bound, and next, the deepest depth, and next, the oldest time-stamp

        Cost getCost(Cost delta = MIN_COST) const { return MAX(MIN_COST, cost - delta); }
    };

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

    typedef enum : unsigned char {
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
#ifdef OPENMPI
        friend class serialization::access;
        template <class Archive>
        void serialize(Archive& ar, const unsigned int version)
        {
            ar& varIndex;
            ar& value;
            ar& op;
            ar& reverse;
        }
#endif
    public:
        int varIndex; // variable wcsp's index
        Value value; // variable's value
        ChoicePointOp op; // choice point operation
        bool reverse; // true if the choice point corresponds to the last right branch of an open node

        ChoicePoint() {} // default constructor added to avoid boost/serialization/access.hpp:130:9: error
        ChoicePoint(ChoicePointOp op_, int var_, Value val_, bool rev_)
            : varIndex(var_)
            , value(val_)
            , op(op_)
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
            ar& nbRecomputationNodes;
            ar& lb;
            ar& ub;
            ar& open;
            ar& cp;
            ar& sol;
        }

    public:
        bool hbfs; // false if master cp/open memory is out
        Long nbNodes;
        Long nbBacktracks;
        Long nbDEE;
        Long nbRecomputationNodes;
        Cost lb; // best lower bound known by the master or found by the worker
        Cost ub; // best current solution a.k.a incumbent solution
        vector<OpenNode> open; // priority queue which will contain the node(s) to eXchange between the processes
        vector<ChoicePoint> cp; // vector of choice points
        vector<Value> sol;

        Work()
            : hbfs(true)
            , nbNodes(0)
            , nbBacktracks(0)
            , nbDEE(0)
            , nbRecomputationNodes(0)
            , lb(MIN_COST)
            , ub(MAX_COST)
        {
        } // dummy message when stopping

        /**
         * @brief constructor used by the master
         */
        Work(const CPStore& cpMaster_, OpenList& openMaster_, const Cost lbMaster_, const Cost ubMaster_, vector<Value>& sol_)
            : hbfs(true)
            , nbNodes(0)
            , nbBacktracks(0)
            , nbDEE(0)
            , nbRecomputationNodes(0)
            , lb(lbMaster_)
            , ub(ubMaster_)
        {

            open.push_back(openMaster_.top());

            openMaster_.pop(); // pop up directly the open queue!!

            for (ptrdiff_t i = open[0].first; i < open[0].last; i++) { // create a sequence of decisions
                assert(i >= 0 && (size_t)i < cpMaster_.size());
                assert(i < cpMaster_.stop);
                cp.push_back(cpMaster_[i]);
                assert(cp.back().op >= 0 && cp.back().op < CP_MAX);
            }

            if (sol_.size() > 0) {
                //            sol.swap(sol_); //after the swap sol_ in the worker is an empty vector
                sol = sol_;
            }
        }

        /**
         * @brief constructor used by the worker
         */
        Work(CPStore& cpWorker_, OpenList& openWorker_, const Long nbNodes_, const Long nbBacktracks_, const Long nbDEE_, const Long nbRecomputationNodes_, const Cost lbWorker_, const Cost ubWorker_, vector<Value>& sol_)
            : hbfs(true)
            , nbNodes(nbNodes_)
            , nbBacktracks(nbBacktracks_)
            , nbDEE(nbDEE_)
            , nbRecomputationNodes(nbRecomputationNodes_)
            , lb(lbWorker_)
            , ub(ubWorker_)
        {

            ptrdiff_t cpmax = 0; // stop collecting choice points if unused. however, we must start collecting choice points at position 0 to be consistent with stored node.first and node.last even if they are not used.
            while (!openWorker_.empty()) // fill-in local vector of open nodes
            {
                OpenNode node = openWorker_.top();
                assert(node.first >= 0);
                assert(node.first <= node.last);
                assert(node.last <= cpWorker_.stop);
                cpmax = MAX(cpmax, node.last);

                open.push_back(node);

                openWorker_.pop(); // pop up directly the open queue!!
            }

            if (!ToulBar2::burst) {
                openWorker_.init(); // clb=cub=MAX_COST  method added to init openList attributes
                // because emptying an openList is not enough we have to initialize its attributes
            }

            if (open.size() > 0) {
                assert(cpmax <= cpWorker_.stop);
                for (ptrdiff_t i = 0; i < cpmax; i++) {
                    assert((size_t)i < cpWorker_.size());
                    cp.push_back(cpWorker_[i]);
                    assert(cp.back().op >= 0 && cp.back().op < CP_MAX);
                }
            }

            if (!ToulBar2::burst)
                cpWorker_.clear(); // size = 0  added to put new cp out of the while(1)
            //            cpWorker_.shrink_to_fit(); // to win space we shrink the vector: capacity=0 //TODO: test if it speeds-up things or not

            if (sol_.size() > 0) {
                //            sol.swap(sol_); //after the swap sol_ in the worker is an empty vector
                sol = sol_;
            }
        }

        friend ostream& operator<<(ostream& os, Work& work)
        {
            os << "[lb: " << work.lb << ", ub: " << work.ub << ", nbvars: " << work.sol.size() << ", nodes: " << work.open.size() << ", choices: " << work.cp.size() << "]";
            return os;
        }
    };
#endif

    // returns 0 if left preferred or 1 if equal or 2 if right preferred
    // prefer open nodes with best (smallest dom/max_at_any_depth(wdeg+1)) variable heuristic values first (lexicographic order)
    // in case of prefix equality, prefer the shortest (included) open node to favor discrepancy first
    static int recHeuristicCmp(const OpenNode& left, ptrdiff_t curLeft, const OpenNode& right, ptrdiff_t curRight) {
        if (curLeft < left.last) {
            if (curRight < right.last) {
                int varLeft = (*Solver::CurrentSolver->cp)[curLeft].varIndex;
                int varRight  = (*Solver::CurrentSolver->cp)[curRight].varIndex;
                double heurLeft = (double) Solver::CurrentSolver->getWCSP()->getDomainSize(varLeft) / (double)(Solver::CurrentSolver->heuristics[varLeft] + 1);
                double heurRight = (double) Solver::CurrentSolver->getWCSP()->getDomainSize(varRight) / (double)(Solver::CurrentSolver->heuristics[varRight] + 1);
                if (heurLeft < heurRight - (double)ToulBar2::epsilon * heurRight) {
                    return 0;
                } else if (heurLeft < heurRight + (double)ToulBar2::epsilon * heurRight) {
                    Cost costLeft = Solver::CurrentSolver->getWCSP()->getMaxUnaryCost(varLeft);
                    Cost costRight = Solver::CurrentSolver->getWCSP()->getMaxUnaryCost(varRight);
                    if (costLeft > costRight) {
                        return 0;
                    } else if (costLeft < costRight) {
                        return 2;
                    } else {
                        return recHeuristicCmp(left, curLeft+1, right, curRight+1);
                    }
                } else {
                    return 2;
                }
            } else { // right=last
                return 2;
            }
        } else { // cur=last
            if (curRight < right.last) {
                return 0;
            } else { //right=last
                return 1;
            }
        }
    }

    void addChoicePoint(ChoicePointOp op, int varIndex, Value value, bool reverse);
    void addOpenNode(CPStore& cp, OpenList& open, Cost lb, Cost delta = MIN_COST); ///< \param delta cost moved out from the cluster by soft arc consistency
    void restore(CPStore& cp, OpenNode node);
    char opSymbol(const CPStore& cp, const ptrdiff_t idx, OpenNode nd);
    void epsDumpSubProblems(CPStore& cp, OpenList& open);

protected:
    friend class NeighborhoodStructure;
    friend class RandomNeighborhoodChoice;
    friend class ClustersNeighborhoodStructure;
    friend class RandomClusterChoice;
    friend class ParallelRandomClusterChoice;
    friend class VACExtension;

    bool self_wcsp; // true if the wcsp has been created inside the solver object, false otherwise
    Long nbNodes;
    Long nbBacktracks;
    Long nbBacktracksLimit;
    WeightedCSP* wcsp;
    vector<DLink<Value>*> allVars;
    BTList<Value>* unassignedVars;
    int lastConflictVar;
    void* searchSize;
    vector<Long> heuristics; // precomputed heuristic value for each variable

    BigInteger nbSol;
    Long nbSGoods; // number of #good which created
    Long nbSGoodsUse; // number of #good which used
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

    // only for pretty print of optimality gap information
    Cost initialLowerBound;
    Cost globalLowerBound;
    Cost globalUpperBound;
    int initialDepth;
    void initGap(Cost newlb, Cost newub);
    void showGap(Cost newlb, Cost newub);
    Double getDDualBound() const { return wcsp->Cost2ADCost(globalLowerBound); }

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
    int greedy(intFunctionCall_t func);

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
    void singletonConsistency(int restricted = INT_MAX);
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
    double hbfsWaitingTime;
    Long initWorkerNbNodes;
    Long initWorkerNbBacktracks;
    Long initWorkerNbDEE;
    Long initWorkerNbRecomputationNodes;
    mpi::communicator world;
    queue<int> idleQ; // MASTER ONLY container with the rank of the free workers
    std::unordered_map<int, OpenNode> activeWork; // MASTER ONLY map the rank i of a worker with an open node currently explored by the worker
    std::unordered_map<int, Cost> bestsolWork; // MASTER ONLY map the rank i of a worker with the cost of the best solution sent by the master

    inline bool MPI_interrupted()
    {
        if (world.iprobe(mpi::any_source, DIETAG))
            return true;
        else
            return false;
    }
    pair<Cost, Cost> hybridSolveMaster(Cluster* root, Cost clb, Cost cub);
    pair<Cost, Cost> hybridSolveWorker(Cluster* root, Cost clb, Cost cub);
#endif

    pair<Cost, Cost> russianDollSearch(Cluster* c, Cost cub);

    BigInteger binaryChoicePointSBTD(Cluster* cluster, int varIndex, Value value);
    BigInteger sharpBTD(Cluster* cluster);
    void approximate(BigInteger& nbsol, TreeDecomposition* td);

    Cost logZCurrentEstimate();

public:
    Solver(Cost initUpperBound, WeightedCSP* wcsp = NULL);

    virtual ~Solver();

    Cost read_wcsp(const char* fileName);
    void read_random(int n, int m, vector<int>& p, int seed, bool forceSubModular = false, string globalname = "");

    Long getNbNodes() const FINAL { return nbNodes; }
    Long getNbBacktracks() const FINAL { return nbBacktracks; }
    set<int> getUnassignedVars() const;
    int numberOfUnassignedVariables() const; // faster than its WCSP linear-time counterpart, but it is valid only during search (otherwise returns -1)
    void updateVarHeuristic(); ///< \brief to be called if DAC order has been changed after preprocessing (initVarHeuristic call)

    virtual bool solve(bool first = true);

    // internal methods called by solve, for advanced programmers only!!!
    void beginSolve(Cost ub);
    Cost preprocessing(Cost ub);
    void recursiveSolve(Cost lb = MIN_COST);
    void recursiveSolveLDS(int discrepancy);
    pair<Cost, Cost> hybridSolve() { return hybridSolve(NULL, wcsp->getLb(), wcsp->getUb()); }
    void endSolve(bool isSolution, Cost cost, bool isComplete);
    // end of internal solve methods

    // INCOP local search
    Cost narycsp(string cmd, vector<Value>& solution);

    // PILS local search
    /// \brief solves the current problem using PILS local search @ Francois Beuvin, David Simoncini, Sebastien Verel
    /// \param solution initial starting solution for the first run only
    /// \param nbruns number of runs
    /// \param perturb_id perturbation mode: 0 means static perturbation, 1: random perturbation, 2: adaptive perturbation
    /// \param perturb_s size of the static perturbation (proportional to sequence length)
    /// \param flatMax number of iterations allowed without improvement
    /// \param nEvalHC maximum number of evaluations during each steepest descent
    /// \param nEvalMax maximum total number of evaluations
    /// \param strengthMin minimum size of the random/adaptive perturbation (proportional to sequence length)
    /// \param strengthMax maximum size of the random/adaptive perturbation (proportional to sequence length)
    /// \param incrFactor increasing factor of perturbation size
    /// \param decrFactor decreasing factor of perturbation size
    /// \return best solution cost found (and updates upper bound and best solution found as a side-effect)
    /// \warning cannot solve problems with non-binary cost functions (it ignores them!)
    Cost pils(vector<Value>& solution, int nbruns = 3, int perturb_id = 0, double perturb_s = 0.333, unsigned long long flatMax = 100, unsigned long long nEvalHC = 500, unsigned long long nEvalMax = 10000, double strengthMin = 0.1, double strengthMax = 0.5, double incrFactor = 0.1, double decrFactor = 0.1);
    Cost pils(string cmd, vector<Value>& solution);

    // LR-BCD local search
    Cost lrBCD(size_t maxiter, int k, size_t nbR, vector<Value>& solution);
    Cost lrBCD(string cmd, vector<Value>& solution);

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
    friend void newsolution(int wcspId, void* solver);

    WeightedCSP* getWCSP() FINAL { return wcsp; }
};

#ifdef OPENMPI
BOOST_IS_MPI_DATATYPE(Solver::OpenNode)
BOOST_IS_BITWISE_SERIALIZABLE(Solver::OpenNode) // only for simple datatypes with no pointer or any varying size objects (vectors, etc.)
BOOST_CLASS_IMPLEMENTATION(Solver::OpenNode, serialization::object_serializable)
BOOST_CLASS_TRACKING(Solver::OpenNode, serialization::track_never)

BOOST_IS_MPI_DATATYPE(Solver::ChoicePoint)
BOOST_IS_BITWISE_SERIALIZABLE(Solver::ChoicePoint)
BOOST_CLASS_IMPLEMENTATION(Solver::ChoicePoint, serialization::object_serializable)
BOOST_CLASS_TRACKING(Solver::ChoicePoint, serialization::track_never)

BOOST_CLASS_IMPLEMENTATION(Solver::Work, serialization::object_serializable)
BOOST_CLASS_TRACKING(Solver::Work, serialization::track_never)
#endif

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

class BestSolFound : public SolverOut {
public:
    const char* what() const throw() FINAL { return "... best known solution found!"; }
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

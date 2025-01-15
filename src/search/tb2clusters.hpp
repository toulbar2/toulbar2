/** \file tb2clusters.hpp
 *  \brief Cluster Tree Decomposition data-structures.
 *
 */

#ifndef TB2CLUSTERS_HPP_
#define TB2CLUSTERS_HPP_

#include "core/tb2wcsp.hpp"
#include "tb2solver.hpp"
#include "core/tb2enumvar.hpp"
#include "core/tb2naryconstr.hpp"

class Cluster;

typedef set<int> TVars;
typedef ConstraintSet TCtrs;
// typedef map<int,Value>     TAssign;

// sort variables by their DAC order
struct CmpVarStruct {
    static WCSP* wcsp;
    bool operator()(const int lhs, const int rhs) const;
};
typedef set<int, CmpVarStruct> TVarsSorted;

// sort clusters by their id if non-negative else by pointer addresses (warning! stochastic behavior!!)
struct CmpClusterStructBasic {
    bool operator()(const Cluster* lhs, const Cluster* rhs) const;
};
typedef set<Cluster*, CmpClusterStructBasic> TClusters;

// data structure for connected components
typedef set<TClusters> component;
// cluster visited or not
typedef map<Cluster*, bool> cluster_visited;
// sort cluster sons by mean separator size first and by number of variables in their subtree next
struct CmpClusterStruct {
    bool operator()(const Cluster* lhs, const Cluster* rhs) const;
};
typedef set<Cluster*, CmpClusterStruct> TClustersSorted;

typedef triplet<Cost, Cost, Solver::OpenList> TPairNG;
typedef pair<Cost, Tuple> TPairSol;

typedef map<Tuple, TPairNG> TNoGoods;
typedef map<Tuple, TPairSol> TSols;

// for solution counting :
typedef pair<Cost, BigInteger> TPairSG;
typedef map<Tuple, TPairSG> TSGoods;

typedef map<Tuple, bool> TFrees;
typedef map<Tuple, bool> TFreesSols;
typedef map<Tuple, int> TFreesLimit;

class Separator : public AbstractNaryConstraint {
private:
    Cluster* cluster;
    TVars vars;
    vector<vector<StoreCost>> delta; // structure to record the costs that leave the cluster
    StoreInt isUsed;
    StoreCost lbPrevious;
    StoreInt optPrevious;

    TNoGoods nogoods;
    TSGoods sgoods; // for solution counting
    TSols solutions;
    DLink<Separator*> linkSep; // link to insert the separator in PendingSeparator list

    TFrees frees; // remember freedom status of a given separator assignment
    TFreesLimit freesLimit; // and the corresponding number of failed searches in the subproblem associated to the given separator assignment
    TFreesSols freesSol; // freedom status used when recording a solution for a given separator assignment (the status may change later, thus freesSol <> frees)

    Tuple t; // temporary buffer for a separator tuple
    Tuple s; // temporary buffer for a solution tuple

public:
    Separator(WCSP* wcsp, EnumeratedVariable** scope_in, int arity_in);
    Separator(WCSP* wcsp);

    void setup(Cluster* cluster_in);

    TVars& getVars() { return vars; }
    int getNbVars() { return vars.size(); }
    bool is(int i) { return vars.find(i) != vars.end(); }
    TVars::iterator begin() { return vars.begin(); }
    TVars::iterator end() { return vars.end(); }

    void set(Cost clb, Cost cub, Solver::OpenList** open = NULL);
    bool get(Cost& clb, Cost& cub, Solver::OpenList** open = NULL);

    void setF(bool free); ///< \brief update freedom status of a given separator assignment
    bool getF(bool& free); ///< \brief return true if the corresponding freedom status is found (and update free status), false otherwise
    bool setFInc(); ///< \brief initialize freedom counter if needed (freesLimit) and return true if it is below a given limit (\see ToulBar2::heuristicFreedomLimit)
    void freeIncS(); ///< \brief increase freedom counter

    void setSg(Cost c, BigInteger nb);
    BigInteger getSg(Cost& res, BigInteger& nb);

    void solRec(Cost ub);
    bool solGet(TAssign& a, Tuple& sol, bool& free);

    void resetLb();
    void resetUb();

    void queueSep() { wcsp->queueSeparator(&linkSep); }
    void unqueueSep() { wcsp->unqueueSeparator(&linkSep); }
    bool isInQueueSep() { return wcsp->isInQueueSeparator(&linkSep); }

    void addDelta(unsigned int posvar, Value value, Cost cost)
    {
        assert(posvar < vars.size());
        delta[posvar][wcsp->toIndex(posvar, value)] += cost;
    }
    Cost getCurrentDeltaUb(); // separator variables may be unassigned, returns sum of maximum delta per variable
    Cost getCurrentDeltaLb(); // separator variables may be unassigned, returns sum of minimum delta per variable

    bool used() { return isUsed; }

    void assign(int varIndex);
    void propagate();

    // NaryConstraint methods not used by Separator class
    double computeTightness() { return 0; }
    bool verify() { return true; }
    void increase(int index) {}
    void decrease(int index) {}
    void remove(int index) {}
    // ConstraintSet subConstraint(){TCtrs s; return s;}
    void print(ostream& os);
};

class Cluster {
private:
    static int clusterCounter; ///< count the number of instances of Cluster class
    int instance; ///< instance number
    TreeDecomposition* td;
    WCSP* wcsp;
    int id; // corresponding to the vector index of the cluster in the tree decomposition
    TVars vars; // contains all variables inside a cluster including separator variables
    TVarsSorted sortedVars; // contains all cluster's variables sorted by DAC order after makeRooted is done
    TCtrs ctrs; // intermediate usage by bucket elimination (DO NOT USE)
    TClusters edges; // adjacent clusters (includes parent cluster before makeRooted is done)
    TClustersSorted sortedEdges; // cluster sons are sorted after makeRooted is done

    Cluster* parent; // parent cluster
    TClusters descendants; // set of cluster descendants (including itself)
    TVars varsTree; // set of variables in cluster descendants (including itself)
    TVarsSorted sortedVarsTree; // set of cluster's tree variables sorted by DAC order after makeRooted is done
    vector<bool> quickdescendants;
    map<Cluster*, TVars> quickIntersections; // set of variables corresponding to the intersection between a neighbor cluster and itself

    Separator* sep; // associated separator with parent cluster
    StoreCost lb; // current cluster lower bound deduced by propagation
    Cost ub; // current cluster best known solution cost
    Cost lbRDS; // global cluster lower bound found by RDS
    StoreInt active; // unactive if a nogood including this cluster has been used by propagation

    StoreBigInteger countElimVars;
    int num_part; // for approximation: number of the corresponding partition

    StoreInt freedom_on; // current freedom status for the cluster
    StoreInt isCurrentlyInTD; // true if the cluster is part of the current tree decomposition
    int depth; // current depth of the cluster in the rooted tree decomposition

public:
    Cluster(TreeDecomposition* tdin);
    ~Cluster();

    void setup();

    int getIndex() const { return instance; } ///< \brief instantiation occurrence number of current Cluster object
    int getId() const { return id; } ///< \brief temporary/final index of the cluster in the current tree decomposition
    void setId(int iid) { id = iid; }

    WCSP* getWCSP() const { return wcsp; }
    TreeDecomposition* getTreeDec() const { return td; }

    Separator* getSep() const { return sep; }
    void setSep(Separator* sepin) { sep = sepin; }
    int sepSize() const
    {
        if (sep)
            return sep->arity();
        else
            return 0;
    }
    void deconnectSep(); // deconnect all the constraints on separator variables and assigns separator variables to their support value
    void deconnectDiff(TCtrs& listCtrsTot, TCtrs& listCtrs);
    bool isSepVar(int i)
    {
        if (!sep)
            return false;
        return sep->is(i);
    }

    Solver::CPStore* cp; // choice point cache for open nodes related to this cluster
    Solver::OpenList* open; // list of open nodes related to this cluster
    Long hbfsGlobalLimit; // global limit on number of backtracks for hybrid search on the subproblem rooted to this cluster
    Long hbfsLimit; // local limit on number of backtracks for hybrid search on this cluster only
    Long nbBacktracks; // current number of backtracks related to this cluster
    Long getNbBacktracksClusterTree() const
    {
        Long res = nbBacktracks;
        for (TClusters::const_iterator iter = beginEdges(); iter != endEdges(); ++iter)
            res += (*iter)->getNbBacktracksClusterTree();
        return res;
    }
    vector<Cluster*> sons; // copy of edges allowing sorting

    bool isVar(int i)
    {
        TVars::iterator it = vars.find(i);
        return it != vars.end();
    }
    int getNbVars() const { return vars.size(); }
    TVars& getVars() { return vars; }
    int getNbVarsTree() const { return varsTree.size(); }
    TVars& getVarsTree() { return varsTree; }
    TCtrs getCtrsTree();
    void addVars(TVars& vars);
    void addVar(Variable* x);
    void removeVar(Variable* x);

    void setParent(Cluster* p) { parent = p; }
    Cluster* getParent() const { return parent; }
    TClusters& getEdges() { return edges; }
    void addEdges(TClusters& cls);
    void addEdge(Cluster* c);
    void removeEdge(Cluster* c);
    TClusters::iterator removeEdge(TClusters::iterator it);

    TClusters& getDescendants() { return descendants; }
    bool isEdge(Cluster* c);
    void accelerateDescendants();
    bool isDescendant(Cluster* c) { return quickdescendants[c->getId()]; } ///< test if c is a descendant of this

    void accelerateIntersections();
    void quickIntersection(Cluster* cj, TVars& cjsep);

    TCtrs& getCtrs() { return ctrs; }
    void addCtrs(TCtrs& ctrsin);
    void addCtr(Constraint* c);
    void removeCtr(Constraint* c);
    void clearCtrs();
    void sum(TCtrs& c1, TCtrs& c2, TCtrs& ctout);

    bool isActive() const
    {
        int a = active;
        return a == 1;
    }
    void deactivate();
    void reactivate();

    Cost getLb() { return lb; }
    void setLb(Cost c) { lb = c; }
    void increaseLb(Cost addToLb) { lb += addToLb; }
    Cost getUb() const { return ub; }
    void setUb(Cost c) { ub = c; }
    Cost getLbRDS()
    {
        Cost delta = getCurrentDeltaUb();
        return MAX(lbRDS - delta, MIN_COST);
    }
    void setLbRDS(Cost c)
    {
        assert(!sep || sep->getCurrentDeltaUb() == MIN_COST);
        lbRDS = c;
    }
    Cost getLbRec() const;
    Cost getLbRecRDS();

    void addDelta(int posvar, Value value, Cost cost)
    {
        assert(sep);
        sep->addDelta(posvar, value, cost);
    }
    Cost getCurrentDeltaUb() { return (sep) ? sep->getCurrentDeltaUb() : MIN_COST; }
    Cost getCurrentDeltaLb() { return (sep) ? sep->getCurrentDeltaLb() : MIN_COST; }

    void nogoodRec(Cost clb, Cost cub, Solver::OpenList** open = NULL)
    {
        if (sep)
            sep->set(clb, cub, open);
    }
    bool nogoodGet(Cost& clb, Cost& cub, Solver::OpenList** open = NULL) { return sep->get(clb, cub, open); }

    void resetLbRec();
    void resetUbRec(Cluster* rootCluster);

    void sgoodRec(Cost c, BigInteger nb)
    {
        if (sep)
            sep->setSg(c, nb);
    }
    BigInteger sgoodGet()
    {
        Cost c = MIN_COST;
        BigInteger nb;
        sep->getSg(c, nb);
        return nb;
    }
    BigInteger getCount() { return countElimVars; }
    void multCount(unsigned int s) { countElimVars = countElimVars * s; }
    void cartProduct(BigInteger& cartProd);
    int getPart() { return num_part; }
    void setPart(int num) { num_part = num; }

    // set freedom status for the cluster and the separator
    void freeRec(bool free)
    {
        if (sep)
            sep->setF(free);
        setFreedom(free);
    }
    // the returned boolean says if the information (freedom status attached to the separator assignment) is found or not
    // return true if found else false. update freedom status of the cluster
    bool freeGet()
    {
        if (!sep)
            return false;
        bool free, found;
        found = sep->getF(free);
        if (found)
            setFreedom(free);
        return found;
    }
    // initialize freesLimit if new else checks if it does not overcome the freedom limit and update freedom status
    void freeRecInc()
    {
        bool ok;
        ok = sep && sep->setFInc();
        if (ok)
            setFreedom(true);
        else
            setFreedom(false);
    }
    // increase freesLimit by one (every time we did not improve the search)
    void freeInc()
    {
        if (sep)
            sep->freeIncS();
    }

    void solutionRec(Cost c)
    {
        setUb(c);
        if (sep)
            sep->solRec(c);
    }
    void getSolution(TAssign& sol); // updates sol by the recorded solution found for a separator assignment also given in sol

    void setWCSP2Cluster(); // sets the WCSP to the cluster problem, deconnecting the rest
    void getElimVarOrder(vector<int>& elimVarOrder);

    // freedom of variables choice
    bool getFreedom() { return freedom_on; }
    void setFreedom(bool f) { freedom_on = f; } // { if (ToulBar2::verbose >= 1 && freedom_on != f) cout << " status of cluster " << getId() << " changed from " << freedom_on << " to " << f << endl;  freedom_on = f; }
    bool isVarTree(int i)
    {
        TVars::iterator it = varsTree.find(i);
        return it != varsTree.end();
    }
    bool getIsCurrInTD() { return isCurrentlyInTD; }
    void setIsCurrInTD(bool f) { isCurrentlyInTD = f; }
    int getDepth() { return depth; }
    void setDepth(int d) { depth = d; }

    TVars::iterator beginVars() { return vars.begin(); }
    TVars::iterator endVars() { return vars.end(); }
    TVars::iterator beginVarsTree() { return varsTree.begin(); }
    TVars::iterator endVarsTree() { return varsTree.end(); }
    TVars::iterator beginSep() { return sep->begin(); }
    TVars::iterator endSep() { return sep->end(); }
    TCtrs::iterator beginCtrs() { return ctrs.begin(); }
    TCtrs::iterator endCtrs() { return ctrs.end(); }
    TClusters::iterator beginEdges() const { return edges.begin(); }
    TClusters::iterator endEdges() const { return edges.end(); }
    TClusters::reverse_iterator rbeginEdges() const { return edges.rbegin(); }
    TClusters::reverse_iterator rendEdges() const { return edges.rend(); }
    TClusters::iterator beginDescendants() { return descendants.begin(); }
    TClusters::iterator endDescendants() { return descendants.end(); }

    void sortEdgesRec()
    {
        for (TClusters::iterator iter = beginEdges(); iter != endEdges(); ++iter)
            (*iter)->sortEdgesRec();
        TClustersSorted tmpset(edges.begin(), edges.end(), CmpClusterStruct());
        sortedEdges = tmpset;
    }
    TClusters::iterator beginSortedEdges() const { return sortedEdges.begin(); }
    TClusters::iterator endSortedEdges() const { return sortedEdges.end(); }

    void sortVarsRec()
    {
        for (TClusters::iterator iter = beginEdges(); iter != endEdges(); ++iter)
            (*iter)->sortVarsRec();
        TVarsSorted tmpset(vars.begin(), vars.end(), CmpVarStruct());
        sortedVars = tmpset;
    }
    TVarsSorted::iterator beginSortedVars() { return sortedVars.begin(); }
    TVarsSorted::iterator endSortedVars() { return sortedVars.end(); }

    void sortVarsTreeRec()
    {
        for (TClusters::iterator iter = beginEdges(); iter != endEdges(); ++iter)
            (*iter)->sortVarsTreeRec();
        TVarsSorted tmpset(varsTree.begin(), varsTree.end(), CmpVarStruct());
        sortedVarsTree = tmpset;
    }
    TVarsSorted::iterator beginSortedVarsTree() { return sortedVarsTree.begin(); }
    TVarsSorted::iterator endSortedVarsTree() { return sortedVarsTree.end(); }

    void print();
    void dump();
    void printStats()
    {
        if (!sep)
            return;
        sep->print(cout);
    }

    void printStatsRec()
    {
        TClusters::iterator it = beginSortedEdges();
        while (it != endSortedEdges()) {
            (*it)->sep->print(cout);
            (*it)->printStatsRec();
            ++it;
        }
    }

    bool isLeaf() { return sortedEdges.begin() == sortedEdges.end(); }
};

class TreeDecomposition {
private:
    WCSP* wcsp;
    vector<Cluster*> clusters;
    std::list<Cluster*> roots; // intermediate list used by makeRooted method, only one root at the end

    Cluster* rootRDS; // root cluster of the current RDS iteration

    StoreInt currentCluster; // used to restrict local propagation (NC) and boosting by variable elimination to the current cluster's subtree
    vector<StoreInt> deltaModified; // accelerator avoiding unnecessary checks to delta structure if it is empty (Boolean value)

    int max_depth;

public:
    // connected component of the tree
    component tree_component;
    // temporary components
    TClusters comp;

    TreeDecomposition(WCSP* wcsp_in);

    ~TreeDecomposition();

    WCSP* getWCSP() { return wcsp; }

    int getNbOfClusters() { return clusters.size(); }
    Cluster* getCluster(int i)
    {
        assert(0 <= i && i < (int)clusters.size());
        return clusters[i];
    }

    void setCurrentCluster(Cluster* c) { currentCluster = c->getId(); }
    Cluster* getCurrentCluster() { return getCluster(currentCluster); }

    bool isInCurrentClusterSubTree(int idc);
    bool isActiveAndInCurrentClusterSubTree(int idc);

    // main function to build a cluster tree/path decomposition:
    //  - builds a bucket for each variable following a given variable elimination order or directly from a tree decomposition file
    //  - builds a tree/path decomposition from the buckets
    //  - associate constraints to clusters, with special treatment for ternary constraints (duplicate flag)
    void buildFromCovering(string filename);
    void buildFromOrder();
    void buildFromOrderNext(vector<int>& order);
    void getElimVarOrder(vector<int>& elimVarOrder);
    void treeFusions(); // merges all redundant clusters
    bool treeFusion(); // one fusion step
    void pathFusions(vector<int>& order); // builds a path decomposition of clusters from a given order

    void buildFromOrderForApprox(); // builds a decomposition for approximation solution counting
    void maxchord(int sizepart, vector<int>& order, ConstraintSet& totalusedctrs, TVars& inusedvars, TVars& currentusedvars, vector<Variable*>& currentRevElimOrder, ConstraintSet& currentusedctrs);
    void insert(int sizepart, vector<Variable*> currentRevElimOrder, ConstraintSet currentusedctrs);

    void fusion(Cluster* ci, Cluster* cj);
    void reduceHeight(Cluster* c, vector<Cluster*> path);
    int getNextUnassignedVar(TVars* vars);
    int getVarMinDomainDivMaxWeightedDegree(TVars* vars);
    void splitClusterRec(Cluster* c, Cluster* father, unsigned int maxsize, TClusters& unvisited);
    TVars boostingVarElimRec(Cluster* c, Cluster* father, Cluster* grandfather, unsigned int maxsize, TClusters& unvisited);
    void mergeClusterRec(Cluster* c, Cluster* father, unsigned int maxsepsize, unsigned int minpropervar, TClusters& unvisited);
    void heuristicFusionRec(Cluster* c, Cluster* noc);

    void makeDescendants(Cluster* c);
    bool isDescendant(Variable* x, Variable* y) { return getCluster(x->getCluster())->isDescendant(getCluster(y->getCluster())); }

    int makeRooted(); // defines a rooted cluster tree decomposition from an undirected one

    // reachable clusters from a cluster
    void DFSUtil(Cluster* c, cluster_visited& c_visited);

    // all cluster connected components of the tree decomposition
    int connectedComponents();

    void makeRootedRec(Cluster* c, Cluster* father, TClusters& unvisited);
    Cluster* getBiggerCluster(TClusters& unvisited);
    Cluster* getCluster_height_rootsize_min(TClusters& unvisited);
    Cluster* getCluster_height_rootsize_max(TClusters& unvisited);
    Cluster* getClusterMinHeight(TClusters& unvisited);
    Cluster* getRoot() { return roots.front(); }
    Cluster* getRootRDS() { return rootRDS; }
    void setRootRDS(Cluster* rdsroot) { rootRDS = rdsroot; }
    void setDuplicates(bool init = false); // deal with two or more ternary constraints having the same included binary constraint belonging to different clusters

    int height(Cluster* r);
    int height(Cluster* r, Cluster* father);

    void intersection(TVars& v1, TVars& v2, TVars& vout);
    void intersection(Cluster* c, Cluster* cj, TVars& vout);
    void difference(TVars& v1, TVars& v2, TVars& vout);
    void sum(TVars v1, TVars v2, TVars& vout); // copy inputs to avoid any overlap issues with output
    void sum(TVars& v1, TVars& v2); // it assumes vout = v1
    bool included(TVars& v1, TVars& v2); // true if v1 is included in v2
    void clusterSum(TClusters v1, TClusters v2, TClusters& vout);
    void clusterSum(TClusters& v1, TClusters& v2);
    void ctrSum(TCtrs v1, TCtrs v2, TCtrs& vout);
    void ctrSum(TCtrs& v1, TCtrs& v2);

    bool isDeltaModified(int varIndex) { return deltaModified[varIndex]; }
    Cost getLbRecRDS()
    {
        Cluster* c = getCluster(currentCluster);
        Cost res = c->getLbRecRDS();
        return MAX(res, c->getLbRDS());
    }
    void addDelta(int c, EnumeratedVariable* x, Value value, Cost cost);
    void newSolution(Cost lb);

    bool verify();

    void print(Cluster* c = NULL, int recnum = 0);
    void printStats(Cluster* c = NULL);
    void dump(Cluster* c = NULL);

    // manage freedom
    void updateInTD(Cluster* c); ///< \brief deconnect cluster subtree and its separators according to the current tree decomposition

    void computeDepths(Cluster* c, int parent_depth); ///< \brief recursively set depth of all clusters in a given cluster subtree
    int getMaxDepth() { return max_depth; }
    Cluster* lowestCommonAncestor(Cluster* c1, Cluster* c2); ///< \brief compute the lowest common ancestor cluster of two clusters in a rooted tree decomposition
    bool isSameCluster(Cluster* c1, Cluster* c2); ///< \brief return true if both clusters are the same or they have been merged together by adaptive BTD
    bool isSameCluster(int c1, int c2) { return c1 == c2 || isSameCluster(getCluster(c1), getCluster(c2)); }
};

#endif

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

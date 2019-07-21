/** \file tb2solver.hpp
 *  \brief Generic solver.
 *
 */
#ifndef TB2SOLVER_HPP_
#define TB2SOLVER_HPP_

#include "toulbar2lib.hpp"
#include "utils/tb2store.hpp"
//kad
#include <utils/tb2files_kad.hpp>

#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>

/*
 #include <mpi.h>
 #include <boost/serialization/array.hpp>
 #include <boost/serialization/string.hpp>
 #include <boost/mpi/datatype.hpp> // for optimization during send if objects contain only PODs: int,float, ...
 #include <boost/serialization/priority_queue.hpp>
 #include <boost/serialization/stack.hpp>
 #include <boost/serialization/queue.hpp>
 #include <boost/serialization/deque.hpp>
 #include <boost/serialization/list.hpp>
 */

template<class T>
class DLink;
template<class T>
class BTList;

class NeighborhoodStructure;
class RandomNeighborhoodChoice;
class ClustersNeighborhoodStructure;
class RandomClusterChoice;
class ParallelRandomClusterChoice;

const double epsilon = 1e-6; // 1./100001.

class Solver: public WeightedCSPSolver {
public:
	class OpenNode {
	private:
		friend class boost::serialization::access;
		template<class Archive>
		void serialize(Archive &ar, const unsigned int version) {
			ar & cost;  // node lower bound
			ar & first; // pointer of type intptr_t = ptrdiff_t based on signed integer type
			ar & last; // means the "last" choice point in CPStore = vector<ChoicePoint> is at the adr (last-1)
		}
	protected: // \Warning cost changed from private to protected to be accessed directly by derived class OpenNode2
		Cost cost; // global lower bound associated with the open node
	public:
		ptrdiff_t first; // first position in the list of choice points corresponding to a branch in order to reconstruct the open node
		ptrdiff_t last; // last position (excluded) in the list of choice points corresponding to a branch in order to reconstruct the open node

		OpenNode(Cost cost_, ptrdiff_t first_, ptrdiff_t last_) :
				cost(cost_), first(first_), last(last_) {
		}
		OpenNode() {
		}
		bool operator<(const OpenNode &right) const {
			return (cost > right.cost)
					|| (cost == right.cost
							&& ((last - first) < (right.last - right.first)
									|| ((last - first)
											== (right.last - right.first)
											&& last >= right.last)));
		} // reverse order to get the open node with first, the smallest lower bound, and next, the deepest depth, and next, the oldest time-stamp

		Cost getCost(Cost delta = MIN_COST) const {
			return MAX(MIN_COST, cost - delta);
		}

	};

	class CPStore;
	class OpenList FINAL : public priority_queue<OpenNode> {
	private:
		Cost clb; // current cluster lower bound built from closed nodes (independent of any soft arc consistency cost moves)
		Cost cub;// current cluster upper bound (independent of any soft arc consistency cost moves)
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
		Cost getLb(Cost delta = MIN_COST) const {return MIN(MAX(MIN_COST, clb - delta), (empty() ? MAX_COST : top().getCost(delta)));}

		Cost getClosedNodesLb(Cost delta = MIN_COST) const {return MAX(MIN_COST, clb - delta);}
		void setClosedNodesLb(Cost lb, Cost delta = MIN_COST)
		{
			clb = MAX(MIN_COST, lb + delta);
			assert(clb <= cub);
		}
		void updateClosedNodesLb(Cost lb, Cost delta = MIN_COST) {clb = MIN(clb, MAX(MIN_COST, lb + delta));}

		Cost getUb(Cost delta = MIN_COST) const {return MAX(MIN_COST, cub - delta);}
		void setUb(Cost ub, Cost delta = MIN_COST) {cub = MAX(MIN_COST, ub + delta);}
		void updateUb(Cost ub, Cost delta = MIN_COST)
		{
			Cost tmpub = MAX(MIN_COST, ub + delta);
			cub = MIN(cub, tmpub);
			clb = MIN(clb, tmpub);
		}

		size_type capacity() const {return c.capacity();}
	};

	typedef enum {
		CP_ASSIGN = 0,
		CP_REMOVE = 1,
		CP_INCREASE = 2,
		CP_DECREASE = 3,
		CP_REMOVE_RANGE = 4,
		CP_MAX
	}ChoicePointOp;

	static const string CPOperation[CP_MAX]; // for pretty print

	struct ChoicePoint {
		//kad
	private:
		friend class boost::serialization::access;
		template<class Archive>
		void serialize(Archive & ar, const unsigned int version) {
			ar & op;
			ar & varIndex;
			ar & value;
			ar & reverse;
		}
		//kad
	public:
		ChoicePointOp op;// choice point operation
		int varIndex;// variable wcsp's index
		Value value;// variable's value
		bool reverse;// true if the choice point corresponds to the last right branch of an open node
		ChoicePoint() {}; // kad : default ctor added to avoid boost/serialization/access.hpp:130:9: error
		ChoicePoint(ChoicePointOp op_, int var_, Value val_, bool rev_)
		: op(op_)
		, varIndex(var_)
		, value(val_)
		, reverse(rev_)
		{
		}

	};

//kad
	/**
	 * \brief class to send work to workers in the form of an object i.e. a message in MPI's semantic
	 *
	 */

	class Work
	{
	private:
		friend class boost::serialization::access;
		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			ar & vecCp;
			ar & nodeX; // TODO : Is it necessary to transmit a hole node: maybe nodes.cost and the value (last-first) would suffice ?
			ar & ub;
			ar & sender;// sender rank can probably be taken from mpi status object
			//ar & subProblemId; // to be used only if we want to memorize pair (i,j)
		}

	public:
		// TODO : See if std arrays or built in arrays are usable here and faster than std vectors

		vector<ChoicePoint> vecCp;// vector of choice points

		vector<OpenNode> nodeX;// priority queue which will contain the node(s) to eXchange between the processes

		Cost ub;//  Best current solution a.k.a incumbent solution

		int sender;// rank of the process that send the msg. nb: sender rank can be taken from mpi status object but in non blocking mode we have to use probe() function to get the sender or the tag etc.

		// TODO: Do we have to transmit the number of backtracks Z ?

		//Long subProblemId; // to be used only if we want to memorize pair (i,j). the idea is to improve the resilience in case of a worker becomes out of order during computation
		// Work(const CPStore &cp, const vector<OpenNode> & openVec_, const int ub_, const int sender_ = 0, Long subProblemId=0)

		/**
		 * @brief constructor used by the master with 3 arguments that pop up directly the open queue of the master or the worker.
		 * If it is the master, we pop up only one node ( boolean oneNode = true), if it is a worker all the nodes of open are popped up (oneNOde=false)
		 * @param cpMaster_ : CPStore, i.e. vector of choice points, of either the master or the worker
		 * @param openMaster_ : OpenList open (priority queue of OpenNode) of either the master or the worker. open_ is modified by this constructor
		 * @param ubMaster_ : current best Solution of either the master or the worker
		 * @param sender : rank of the sender. By default, if the sender is the master rank = 0 by convention.
		 */
		Work(const CPStore& cpMaster_, OpenList& openMaster_, const int ubMaster_)
		: ub(ubMaster_)
		, sender(0)
		{
			//OpenNode node = openMaster_.top();
			nodeX.push_back(openMaster_.top());
			nodeX.shrink_to_fit(); // to transmit a vector with size = capacity . to be tested if this improve speedup or not. Indeed, boost.MPI may take care correctly of this problem.
			openMaster_.pop(); // pop up directly the queue open_ !!

			// TODO: See if no need to pass a whole node and if transmission of only node's cost + ub + vecCp  is sufficient. sender=0 not necessary and vector nodeX contains first and last which are
			// necessary only in the master to compute vecCp.

			for(ptrdiff_t i = nodeX[0].first; i < nodeX[0].last; i++)  // create a sequence of decisions in the vector vec.  node.last: index in CPStore of the past-the-end element
				vecCp.push_back(cpMaster_[i]);
			vecCp.shrink_to_fit(); // to be tested

		}


		Work(const CPStore & cpWorker_, OpenList & openWorker_, const int ubWorker_,  const int sender_) // ctor used by the workers with 4 arguments. All the cp
		: ub(ubWorker_)
		, sender(sender_)
		{

			while(!openWorker_.empty()) // init of vector of OpenNode nodeX
			{
				OpenNode node = openWorker_.top();
				nodeX.push_back(node);
				openWorker_.pop(); // pop up directly the queue open !!
			}

			nodeX.shrink_to_fit(); // to transmit a vector with size = capacity

			for(size_t i = 0; i < cpWorker_.size(); i++) // init of vector of ChoicePoint with ALL THE CPs in cpWorker_
				vecCp.push_back(cpWorker_[i]);

			vecCp.shrink_to_fit();

		}

		Work() {}

	};

	/*
	 * To optimize when it is a class of plain old objects POD : int, long, ...
	 * \warning possible errors seg fault with pointers attributes and non POD attributes
	 */
	//BOOST_IS_MPI_DATATYPE(Work);
	//kad

	class CPStore FINAL : public vector<ChoicePoint> {
	public:
		ptrdiff_t start;	// beginning of the current branch
		ptrdiff_t stop;// deepest saved branch end (should be free at this position)
		StoreCost index;// current branch depth (should be free at this position)

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

	Long nbNodes;
	Long nbBacktracks;
	Long nbBacktracksLimit;
	WeightedCSP* wcsp;
	DLink<Value>* allVars;
	BTList<Value>* unassignedVars;
	int lastConflictVar;
	void* searchSize;

	BigInteger nbSol;
	Long nbSGoods;//number of #good which created
	Long nbSGoodsUse;//number of #good which used
	map<int, BigInteger> ubSol;// upper bound of solution number
	double timeDeconnect;// time for the disconnection
	int tailleSep;

	CPStore* cp;// choice point cache for open nodes (except BTD)
	OpenList* open;// list of open nodes (except BTD)
	Long hbfsLimit;// limit on number of backtracks for hybrid search (except BTD)
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

	//kad
	/**
	 * @brief create the string of partial Assignments to put in subProblems.txt i.e. lines like -x=",...."
	 * @param cp
	 * @param open
	 * @param nbCores
	 * @return string of sub problems
	 */
	string epsSubProblems(const CPStore& cp, OpenList& open, const int nbCores);
	//string epsCommand(const CPStore& cp, OpenList& open, const int nbCores);
	string opSymbol(const CPStore& cp, const ptrdiff_t idx, OpenNode nd);
	//kad

	Cost beginSolve(Cost ub);
	Cost preprocessing(Cost ub);
	void endSolve(bool isSolution, Cost cost, bool isComplete);
	void binaryChoicePoint(int xIndex, Value value, Cost lb = MIN_COST);
	void binaryChoicePointLDS(int xIndex, Value value, int discrepancy);
	void narySortedChoicePoint(int xIndex, Cost lb = MIN_COST);
	void narySortedChoicePointLDS(int xIndex, int discrepancy);
	void recursiveSolve(Cost lb = MIN_COST);
	void recursiveSolveLDS(int discrepancy);
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
	pair<Cost, Cost> hybridSolve() {return hybridSolve(NULL, wcsp->getLb(), wcsp->getUb());}
	//kad
	/**
	 * @brief hybrid B&B parallelized with MPI with paradigm master-worker. Activated by -para option
	 * The master send in non blocking mode the first node to worker 1 with clb=w0 and cub=k (see notations of wcsp formalism)
	 * and wait for (ie receive in blocking mode) the worker to send its open nodes and UB.
	 * Then the master update its open queue and cp, the vector of choice points aka decisions,
	 * and distribute these nodes to worker 2 3 .... 7 1 (case of 7 workers, 1 master).
	 * It is the number of backtrack authorized that give the signal to the worker to send its nodes.
	 *
	 * TODO : if possible, to avoid master bottleneck because it centralizes all the communications, maybe a better scalabilty
	 * could be reach by transferring the initiative to return nodes to the master. This way the workers won't be idle
	 * because the master would be occupied to transfer nodes and communicate with a large number of workers, say 100 000 !
	 * In this case, the workers will pursue their DFS waiting for the master to tell them to return their open nodes.
	 */
	pair<Cost, Cost> hybridSolvePara(Cost clb, Cost cub);
	pair<Cost, Cost> hybridSolveParaBck(Cost clb, Cost cub);//kad temporary backup. do nothing .has to be deleted at some point
	pair<Cost, Cost> hybridSolveParaBck2(Cost clb, Cost cub);//kad temporary backup.
	pair<Cost, Cost> hybridSolvePara() {return hybridSolvePara(wcsp->getLb(), wcsp->getUb());}
	//kad
	pair<Cost, Cost> russianDollSearch(Cluster* c, Cost cub);

	BigInteger binaryChoicePointSBTD(Cluster* cluster, int varIndex, Value value);
	BigInteger sharpBTD(Cluster* cluster);
	void approximate(BigInteger& nbsol, TreeDecomposition* td);

public:
	Solver(Cost initUpperBound);
	~Solver();

	Cost read_wcsp(const char* fileName);
	void read_random(int n, int m, vector<int>& p, int seed, bool forceSubModular = false, string globalname = "");

	Long getNbNodes() const FINAL {return nbNodes;}
	Long getNbBacktracks() const FINAL {return nbBacktracks;}
	set<int> getUnassignedVars() const;

	virtual bool solve();

	Cost narycsp(string cmd, vector<Value>& solution);

	bool solve_symmax2sat(int n, int m, int* posx, int* posy, double* cost, int* sol);

	void dump_wcsp(const char* fileName, bool original = true);
	void read_solution(const char* fileName, bool updateValueHeuristic = true);
	void parse_solution(const char* certificate);

	virtual void newSolution();
	Cost getSolution(vector<Value>& solution);

	friend void setvalue(int wcspId, int varIndex, Value value, void* solver);

	WeightedCSP* getWCSP() FINAL {return wcsp;}
};

class NbBacktracksOut {
public:
	NbBacktracksOut() {
		ToulBar2::limited = true;
		if (ToulBar2::verbose >= 2)
			cout << "... limit on the number of backtracks reached!" << endl;
	}
};

class NbSolutionsOut {
public:
	NbSolutionsOut() {
		ToulBar2::limited = true;
		if (ToulBar2::verbose >= 2)
			cout << "... limit on the number of solutions reached!" << endl;
	}
};

class TimeOut {
public:
	TimeOut() {
		ToulBar2::limited = true;
		if (ToulBar2::verbose >= 2)
			cout << "... time limit reached!" << endl;
	}
};

int solveSymMax2SAT(int n, int m, int *posx, int *posy, double *cost, int *sol);
extern "C" int solvesymmax2sat_(int *n, int *m, int *posx, int *posy,
		double *cost, int *sol);

#endif /*TB2SOLVER_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

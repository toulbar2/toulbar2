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
#include <mpi.h>
#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
#include <boost/mpi/datatype.hpp> // for optimization during send if objects contain only PODs: int,float, ...
#include <boost/serialization/vector.hpp>
#include <boost/serialization/priority_queue.hpp>
#include <boost/serialization/stack.hpp>
/*
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
		friend class boost::serialization::access;

		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			ar & boost::serialization::base_object<priority_queue>(*this);
			ar & clb;
			ar & cub;
		}
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
		void serialize(Archive & ar, const unsigned int version)
		{
			ar & op;
			ar & varIndex;
			ar & varIndex;
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

	// OpenNode2 without derivation from OpenNode class
	/*
	 class OpenNode2 {
	 private:
	 friend class boost::serialization::access;

	 template<class Archive>
	 void serialize(Archive &ar, const unsigned int version) {
	 ar & cost;  // node lower bound
	 ar & first; // pointer of type intptr_t = ptrdiff_t based on signed integer type
	 ar & last; // means the "last" choice point in CPStore = vector<ChoicePoint> is at the index (last-1)
	 }
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

	 */

	// OpenNode using derivation
	class OpenNode2 : public OpenNode {
	private:
		friend class boost::serialization::access;

		template<class Archive>
		void serialize(Archive &ar, const unsigned int version) {
			ar & boost::serialization::base_object<OpenNode>(*this);
			ar & decisions;
		}
	public:
		vector<ChoicePoint> decisions;

		OpenNode2() {
		}
		/*
		 OpenNode2(Cost cost_, ptrdiff_t first_, ptrdiff_t last_)
		 : OpenNode(cost_, first_, last_)
		 {
		 }
		 */
		OpenNode2(Cost cost_, const vector<ChoicePoint> & decisions_)
		{
			decisions = decisions_;
			first = 0;
			last = (ptrdiff_t )decisions.size();  // bad cast?!
			cost = cost_;
		}

	};

	class OpenList2 FINAL : public priority_queue<OpenNode2> {
	private:
		friend class boost::serialization::access;

		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			ar & boost::serialization::base_object<priority_queue>(*this);
			ar & clb;
			ar & cub;
		}
		Cost clb; // current cluster lower bound built from closed nodes (independent of any soft arc consistency cost moves)
		Cost cub;// current cluster upper bound (independent of any soft arc consistency cost moves)
	public:
		OpenList2(Cost lb, Cost ub)
		: clb(lb)
		, cub(ub)
		{
		}
		OpenList2()
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

			ar & nodeVec;
			ar & ub;
			ar & vecDecisions;
			ar & sender; // sender rank can probably be taken from mpi status object
			//ar & subProblemId; // to be used only if we want to memorize pair (i,j)
		}

	public:
		vector<OpenNode> nodeVec; // priority queue which will contain the node(s) to send to other processes

		Cost ub; //  Best current solution a.k.a incumbent solution

		int sender; // rank of the process that send the msg. nb: sender rank can probably be taken from mpi status object

		vector<vector<ChoicePoint>> vecDecisions;

		// TODO: Do we have to transmit the number of backtracks Z ?

		//Long subProblemId; // to be used only if we want to memorize pair (i,j)

		// Work(const CPStore &cp, const vector<OpenNode> & openVec_, const int ub_, const int sender_ = 0, Long subProblemId=0)


		Work(const CPStore &cp, const vector<OpenNode> & nodeVec_, const int ub_, const int sender_=0)
		: nodeVec(nodeVec_)
		, ub(ub_)
		, sender(sender_)
		{
			for (size_t j = 0; j< nodeVec.size();j++)
			{
				OpenNode node = nodeVec[j];
				vector<ChoicePoint> vec;
				for(ptrdiff_t i = node.first; i < node.last; i++)
				{
					vec.push_back(cp[i]);
				}
				vecDecisions.push_back(vec);
				vec.clear();
			}
		}
/*
		 void vector2CPStore(CPStore *cp)
		 {
		 for(size_t i = 0; i<=decisions.size()-1; i++)
		 cp->push_back(decisions[i]);
		 }
*/
		Work() {}

	};

	/*
	 * To optimize when it is a class of plain old objects POD : int, long, ...
	 * \warning possible errors seg fault with pointers attributes and non POD attributes
	 */
	//BOOST_IS_MPI_DATATYPE(MasterToWorker);
	//kad
	class CPStore FINAL : public vector<ChoicePoint> {
		//kad
		/* private:
		 friend class boost::serialization::access;

		 template<class Archive>
		 void serialize(Archive & ar, const unsigned int version)
		 {
		 ar & boost::serialization::base_object<vector>(*this);
		 ar & start;
		 ar & stop;  // attribute
		 ar & index;
		 }
		 */
		//kad
	public:
		ptrdiff_t start;// beginning of the current branch
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
	string epsSubProblems(const CPStore& cp, OpenList& open, const int nbCores); //kad
	//string epsCommand(const CPStore& cp, OpenList& open, const int nbCores); //kad
	string opSymbol(const CPStore& cp, const ptrdiff_t idx, OpenNode nd );//kad
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
	pair<Cost, Cost> hybridSolvePara(Cost clb, Cost cub);
	pair<Cost, Cost> hybridSolveParaBck(Cost clb, Cost cub);//kad temporary backup. do nothing .has to be deleted at some point
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

/** \file tb2treeconstr.hpp
 *  \brief Dynamic programming based global constraint : tree
 */

#ifndef TB2TREECONSTR_HPP_
#define TB2TREECONSTR_HPP_

#include "tb2dpglobalconstr.hpp"
#include "tb2rangeminquery.hpp"

#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <string>

using namespace std;

class TreeConstraint : public DPGlobalConstraint
{
	private:

		int curTreeCost;

		struct Edge {
			int u;
			int v;
			Cost weight;
			Edge(int u, int v, Cost w): u(u), v(v), weight(w) {}
			bool operator< (const Edge &e) const {return weight < e.weight;}
		};

		int minTreeEdgeCost;
		int maxTreeEdgeCost;
		set<pair<int, int> > treeEdge;
                
                struct CCTreeNode;  // Forward declaration                
                vector<CCTreeNode> nodeStore;
                //typedef vector<CCTreeNode>::iterator CCTreeNodePtr;                
                typedef CCTreeNode* CCTreeNodePtr;                
                                
		struct CCTreeNode {		
			int nodeIndex;
			int u;
			int v;
			Cost weight;
			int height;
			CCTreeNodePtr parent;
			CCTreeNodePtr left;
			CCTreeNodePtr right;
			CCTreeNode():nodeIndex(0), u(-1), v(-1), weight(MIN_COST), height(0), parent(NULL), left(NULL), right(NULL) {}
		};                                

		vector<CCTreeNodePtr> ccTree;
		vector<CCTreeNodePtr> inorder;
		vector<CCTreeNodePtr> inorderNodeHeight;
		vector<int> pos;	
		CCTreeNodePtr ccTreeRoot;
		RangeMinQuery<int> RMQ;

                //CCTreeNodePtr PtrNULL() {return nodeStore.end();}
                CCTreeNodePtr PtrNULL() {return NULL;}
                CCTreeNodePtr createNewNode();
                
		void joinCCTrees(int u, int v, Cost weight);
		CCTreeNodePtr findRoot(CCTreeNodePtr node);
		void InorderTransveral(CCTreeNodePtr root);

		// disjoint data set	
		vector<int> p;	
		int findParent(int index, vector<int>& p);
		void unionSet(int u, int v, vector<int>& p);

		map<int, int> val2VarIndex;	

		int recomputeCurMST();				
		int recomputeMST(vector<Edge> &edgeList);				

	protected:

		Cost minCostOriginal();
		Cost minCostOriginal(int var, Value val, bool changed);
		Result minCost(int var, Value val, bool changed);
                
                // This is a hard constraint. SNIC and D(G)AC* are equivalent to AC
                
                void propagateStrongNIC() {
                    propagateAC();
                }
		
		void propagateDAC() {
                   if (ToulBar2::LcLevel == LC_DAC) propagateAC();
                }                                

		// No need to run anything for (weak) ED(G)AC*
		bool isEAC(int var, Value val) {return true;}
		void findFullSupportEAC(int var) {}

	public:
		TreeConstraint(WCSP * wcsp, EnumeratedVariable ** scope, int arity);
		virtual ~TreeConstraint();

		Cost eval(String s);

		void read(istream & file) {} //No parameter needed
                         void initMemoization();
		string getName(){return "MST";}
};

#endif

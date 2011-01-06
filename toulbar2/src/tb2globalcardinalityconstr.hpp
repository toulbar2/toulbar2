//#include "glpk.h"
#include "stddef.h"
#include "tb2flowbasedconstr.hpp"

#define upper_bound first
#define lower_bound second

class GlobalCardinalityConstraint : public FlowBasedGlobalConstraint 
{
	private:
		map<Value, pair<int, int> > bound;
		void buildIndex();
		void buildGraph(Graph &g);
		Cost constructFlow(Graph &g);
		pair<int,int> mapto(int varindex, Value val) {
			return make_pair(varindex+1, mapval[val]);
		}
	public:
		static const int VALUE = 1;
		static const int VAR = 0;
		GlobalCardinalityConstraint(WCSP *wcsp, EnumeratedVariable** scope_in, int
		arity_in);

		~GlobalCardinalityConstraint() {
			/*if (ToulBar2::consistencyLevel != FINE_IC) {
				cout << "no. of GAC propagation = " << count << endl;
				cout << "no. of FDAC propagation = " << count_fdac << endl;
				cout << "no. of error = " << error << endl;
			}*/ 
		}

		string getName() {return "GCC constraint";}
		Cost evalOriginal (String s);
		void read(ifstream &file);

        void print(ostream& os);
        void dump(ostream& os, bool original = true);
};



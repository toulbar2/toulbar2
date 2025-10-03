/** \file tb2globalcardinalityconstr.hpp
 *  \brief Flow based global cost function : sgcc_flow
 */

#ifndef TB2GLOBALCARDCONSTR_HPP_
#define TB2GLOBALCARDCONSTR_HPP_

// #include "glpk.h"
// #include "stddef.h"
#include "tb2flowbasedconstr.hpp"

class SoftGlobalCardinalityConstraint : public FlowBasedGlobalConstraint {
private:
    map<Value, pair<int, int>> bound;
    void buildIndex();
    size_t GetGraphAllocatedSize();
    void buildGraph(Graph& g);
    Cost constructFlow(Graph& g);
    pair<int, int> mapto(int varindex, Value val)
    {
        return make_pair(varindex + 1, mapval[val]);
    }
    // JP Start// This array stores the respective weight of each bound
    map<Value, pair<Cost, Cost>> weights;
    int nDistinctDomainValue;
    // JP End//
public:
    // JP Start// New type
    static const int EMPTY = -1;
    static const int WVALUE = 2;
    // JP End//
    static const int VALUE = 1;
    static const int VAR = 0;
    SoftGlobalCardinalityConstraint(WCSP* wcsp, EnumeratedVariable** scope_in, int arity_in);

    ~SoftGlobalCardinalityConstraint()
    {
        /*if (ToulBar2::consistencyLevel != FINE_IC) {
                                cout << "no. of GAC propagation = " << count << endl;
                                cout << "no. of FDAC propagation = " << count_fdac << endl;
                                cout << "no. of error = " << error << endl;
                        }*/
    }

    string getName();
    Cost evalOriginal(const Tuple& s);
    void read(istream& file, bool mult = true);

    // GlobalCostFunctionParameters* getParameters() {return this;}
    void addValueAndBounds(Value value, int upper = -1, int lower = 0)
    {
        if (upper == -1)
            upper = arity();
        bound[value] = make_pair(lower, upper);
    }
    void addValueAndWeights(Value value, Cost wexcess = -1, Cost wshortage = -1)
    {
        if (wexcess == -1)
            wexcess = def;
        if (wshortage == -1)
            wshortage = def;
        weights[value] = make_pair(wshortage, wexcess);
    }
    void organizeConfig();

    //    void print(ostream& os);
    void dump(ostream& os, bool original = true);
};

#endif /*TB2GLOBALCARDCONSTR_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

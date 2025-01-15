/*
 * \file tb2localsearch.hpp
 * \brief abstract class for local search solvers
 *
 *  Created on: 3 mars 2015
 *      Author: Abdelkader Ouali
 *      PhD Student: LITIO, University of Oran ; GREYC, University of Caen.
 */

#ifndef TB2LOCALSEARCH_HPP_
#define TB2LOCALSEARCH_HPP_

#include "search/tb2solver.hpp"

class LocalSearch : public Solver {
protected:
    map<int, Value> bestSolution;
    map<int, Value> lastSolution;
    Cost bestUb;
    Cost lastUb;

public:
    LocalSearch(Cost initUpperBound);
    virtual ~LocalSearch();

    Cost generateInitSolution(VNSSolutionInitMethod mode, map<int, Value>& solutionInit, bool& complete);
    Cost evaluate_partialInstantiation(vector<int>& variables, vector<Value>& values);
    Cost evaluate_partialInstantiation(map<int, Value>& solution)
    {
        vector<int> variables;
        vector<Value> values;
        for (map<int, Value>::iterator it = solution.begin(); it != solution.end(); ++it) {
            variables.push_back((*it).first);
            values.push_back((*it).second);
        }
        return evaluate_partialInstantiation(variables, values);
    }
    bool repair_recursiveSolve(int discrepancy, vector<int>& variables, vector<Value>& values, Cost ls_ub = MAX_COST); /// \warning if discrepancy>=0 then explores with LDS else with a complete search
    bool repair_recursiveSolve(vector<int>& variables, vector<Value>& values, Cost ls_ub = MAX_COST) { return repair_recursiveSolve(-1, variables, values, ls_ub); } ///  explores with a complete search

    virtual void newSolution();
};

#endif /* TB2LOCALSEARCH_HPP_ */

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

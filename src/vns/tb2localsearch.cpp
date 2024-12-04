/*
 * tb2localsearch.cpp
 *
 *  Created on: 3 mars 2015
 *      Author: Abdelkader Ouali
 *      Phd. Student : LITIO, University of Oran. GREYC, University of Caen.
 */

#include "tb2localsearch.hpp"
#include "core/tb2wcsp.hpp"

LocalSearch::LocalSearch(Cost initUpperBound)
    : Solver(initUpperBound)
    , bestUb(MAX_COST)
    , lastUb(MAX_COST)
{
}

LocalSearch::~LocalSearch()
{
}

void LocalSearch::newSolution()
{
    Solver::newSolution();
    //	static_cast<WCSP*>(wcsp)->registerConflicts();
    lastSolution.clear();
    for (unsigned int i = 0; i < wcsp->numberOfVariables(); ++i) {
        lastSolution[i] = wcsp->getValue(i);
    }
    lastUb = wcsp->getLb();
    if (ToulBar2::lds)
        throw TimeOut(); // force VNS to restart with smallest search parameters at each solution
}

/// This function generates the initial solution
/// \param mode : the generation method
/// \param solutionInit : a map to store the initial solution
Cost LocalSearch::generateInitSolution(VNSSolutionInitMethod mode, map<int, Value>& solutionInit, bool& complete)
{
    if (lastUb < MAX_COST && lastSolution.size() == wcsp->numberOfVariables()) { // reuse INCOP solution or any solution found
        for (map<int, Value>::iterator it = lastSolution.begin(); it != lastSolution.end(); ++it) {
            solutionInit[(*it).first] = (*it).second;
        }
        complete = (lastUb == wcsp->getLb());
        return lastUb;
    }

    Cost cost = MAX_COST;
    complete = false;
    vector<int> dumvariables;
    vector<Value> dumvalues;
    int lds = ToulBar2::lds;
    switch (mode) {
    case LS_INIT_RANDOM:
        if (ToulBar2::verbose >= 1)
            cout << "solution init random" << endl;
        for (unsigned int i = 0; i < wcsp->numberOfVariables(); ++i) {
            int res;
            Value* val = new Value[wcsp->getDomainSize(i)];
            wcsp->getEnumDomain(i, val);
            res = (myrand() % wcsp->getDomainSize(i));
            solutionInit[i] = *(val + res);
            delete[] val;
        }
        cost = evaluate_partialInstantiation(solutionInit);
        break;
    case LS_INIT_INF:
        if (ToulBar2::verbose >= 1)
            cout << "solution init inf" << endl;
        for (unsigned int i = 0; i < wcsp->numberOfVariables(); ++i) {
            solutionInit[i] = wcsp->getInf(i);
        }
        cost = evaluate_partialInstantiation(solutionInit);
        break;
    case LS_INIT_SUP:
        if (ToulBar2::verbose >= 1)
            cout << "solution init sup" << endl;
        for (unsigned int i = 0; i < wcsp->numberOfVariables(); ++i) {
            solutionInit[i] = wcsp->getSup(i);
        }
        cost = evaluate_partialInstantiation(solutionInit);
        break;
    case LS_INIT_DFBB:
        if (ToulBar2::verbose >= 1)
            cout << "solution init DFBB" << endl;
        ToulBar2::lds = 1; // ensures DFBB will stop after the first solution is found if any exists
        complete = repair_recursiveSolve(dumvariables, dumvalues, wcsp->getUb()); // forbidden assignments are NOT allowed!
        assert(complete || lastUb < MAX_COST);
        for (map<int, Value>::iterator it = lastSolution.begin(); it != lastSolution.end(); ++it) {
            solutionInit[(*it).first] = (*it).second;
        }
        cost = lastUb;
        ToulBar2::lds = lds;
        break;
    case LS_INIT_LDS0:
    default: // search using LDS with 0 or more discrepancies
        if (ToulBar2::verbose >= 1)
            cout << "solution init LDS " << mode << endl;
        ToulBar2::lds = 0; // ensures LDS will explore without stopping at the first solution
        complete = repair_recursiveSolve(abs(mode), dumvariables, dumvalues, wcsp->getUb()); // first, forbidden assignments are NOT allowed!
        if (!complete && lastUb == MAX_COST) {
            complete = repair_recursiveSolve(abs(mode), dumvariables, dumvalues, MAX_COST); // if nothing found, forbidden assignments are allowed!
        }
        assert(complete || lastUb < MAX_COST);
        for (map<int, Value>::iterator it = lastSolution.begin(); it != lastSolution.end(); ++it) {
            solutionInit[(*it).first] = (*it).second;
        }
        cost = lastUb;
        ToulBar2::lds = lds;
        break;
    }
    return cost;
}

Cost LocalSearch::evaluate_partialInstantiation(
    vector<int>& variables, vector<Value>& values)
{
    Cost cost = MAX_COST;
    int storedepth = Store::getDepth();
    Store::store();
    try {
        wcsp->setUb(MAX_COST);
        wcsp->assignLS(variables, values);
        cost = wcsp->getLb();
    } catch (const Contradiction&) {
        wcsp->whenContradiction();
    }
    Store::restore(storedepth);
    return cost;
}

bool LocalSearch::repair_recursiveSolve(int discrepancy, vector<int>& variables, vector<Value>& values, Cost ls_ub)
{
    lastUb = MAX_COST;
    lastSolution.clear();
    ToulBar2::limited = false;
    if (ToulBar2::backtrackLimit < LONGLONG_MAX) {
        nbBacktracksLimit = nbBacktracks + ToulBar2::backtrackLimit;
    }
    Long hbfs_ = ToulBar2::hbfs;
    ToulBar2::hbfs = 0; // HBFS not compatible with LDS
    bool solutionBasedPhaseSaving_ = ToulBar2::solutionBasedPhaseSaving;
    ToulBar2::solutionBasedPhaseSaving = false; // VNS prefers randomization inside neighborhood search
    int storedepth = Store::getDepth();
    Cost lb = wcsp->getLb();
    Store::store();
    try {
        wcsp->setUb(ls_ub);
        wcsp->enforceUb();
        wcsp->propagate();
        lb = wcsp->getLb();
        int nbvar = unassignedVars->getSize();
        ToulBar2::limited = true;
        wcsp->assignLS(variables, values);
        if (unassignedVars->getSize() == nbvar)
            ToulBar2::limited = false;
        if (ToulBar2::DEE == 4)
            ToulBar2::DEE_ = 0; // only PSNS in preprocessing
        try {
            try {
                if (discrepancy >= 0)
                    recursiveSolveLDS(discrepancy);
                else
                    recursiveSolve();
            } catch (const NbBacktracksOut&) {
                assert(ToulBar2::limited == true);
            }
        } catch (const TimeOut&) {
            assert(ToulBar2::limited == true);
            if (ToulBar2::interrupted)
                throw TimeOut();
        }
    } catch (const Contradiction&) {
        wcsp->whenContradiction();
    }
    Store::restore(storedepth);
    ToulBar2::hbfs = hbfs_;
    ToulBar2::solutionBasedPhaseSaving = solutionBasedPhaseSaving_;
    return (!ToulBar2::limited || lastUb == lb);
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

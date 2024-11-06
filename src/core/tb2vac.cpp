/** \file tb2vac.cpp
 *  \brief VAC implementation for binary cost functions.
 *
 *  \defgroup VAC Virtual Arc Consistency enforcing
 *  The three phases of VAC are enforced in three different "Pass".
 *  Bool(P) is never built. Instead specific functions (getVACCost) booleanize the WCSP on the fly.
 *  The domain variables of Bool(P) are the original variable domains (saved and restored using trailing at each iteration).
 *  All the counter data-structures (k) are timestamped to avoid clearing them at each iteration.
 *
 *  Note: Simultaneously AC (and potentially DAC, EAC) are maintained by proper queuing.
 *
 *  Note: usual domain events such that assign/remove should not be called during VAC phase 1, use removeVAC instead.
 *
 *  See : <em> Soft Arc Consistency Revisited. </em> Cooper et al. Artificial Intelligence. 2010.
 */

#include <math.h> /* atan2 */
#include "search/tb2clusters.hpp"
#include "tb2vacutils.hpp"
#include "tb2knapsack.hpp"

class tVACStat {
public:
    int var;
    Cost sumlb;
    Long nlb;

    tVACStat(int varin)
    {
        var = varin;
        sumlb = MIN_COST;
        nlb = 0;
    }
};

bool cmp_function(tVACStat* v1, tVACStat* v2)
{
    return v1->sumlb > v2->sumlb;
}

VACExtension::VACExtension(WCSP* w)
    : wcsp(w)
    , nbIterations(0)
    , inconsistentVariable(-1)
    , PBconflict(-1)
    , prevItThreshold(MIN_COST)
    , itThreshold(MIN_COST)
    , breakCycles(0)
    , minlambda(MAX_COST)
    , sumlb(MIN_COST)
    , nlb(0)
    , sumvars(0)
    , sumk(0)
    , theMaxK(0)
    , bneckVar(-1)
    , bneckCF(NULL)
    , bneckCost(MIN_COST)
{
    queueP = new stack<pair<int, Value>>;
    queueR = new stack<pair<int, Value>>;
}

VACExtension::~VACExtension()
{
    delete queueP;
    delete queueR;
}

void VACExtension::init()
{
    VACVariable* xi;
    for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
        xi = (VACVariable*)wcsp->getVar(i);
        xi->setThreshold(MIN_COST);
    }
    iniThreshold();

    // for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
    //    tVACStat* vacinfo = new tVACStat(i);
    //    heapAccess[i] = vacinfo;
    //    heap.insert(heap.end(), vacinfo);
    // }
}

void VACExtension::histogram(Cost c)
{
    if (c != MIN_COST) {
        tScale::iterator it = scaleCost.find(c);
        if (it == scaleCost.end())
            scaleCost[c] = 0;
        else
            it->second++;
    }
}

void VACExtension::histogram()
{
    int cumulus = 0;
    int packetsize = 50;
    bool toomany = true;
    while (toomany) {
        scaleVAC.clear();
        tScale::iterator it = scaleCost.begin();
        while (it != scaleCost.end()) {
            cumulus += it->second;
            if (cumulus > packetsize) {
                scaleVAC.push_front(it->first);
                cumulus = 0;
            }
            ++it;
        }
        toomany = scaleVAC.size() > 20;
        if (toomany)
            packetsize *= 2;
    }

    if (ToulBar2::verbose >= 1) {
        cout << "Reverse Costs Range and Scale: " << scaleCost.rbegin()->first;
        list<Cost>::iterator itl = scaleVAC.begin();
        while (itl != scaleVAC.end()) {
            cout << " " << *itl;
            ++itl;
        }
        cout << " " << scaleCost.begin()->first << endl;
    }
}

void VACExtension::iniThreshold()
{
    if (!ToulBar2::RASPS && Store::getDepth() >= abs(ToulBar2::vac)) {
        return;
    }
    if (scaleCost.size() > 0 && scaleVAC.size() == 0)
        histogram();
    Cost c = MAX_COST;
    bool done = false;
    list<Cost>::iterator it = scaleVAC.begin();
    while ((it != scaleVAC.end()) && !done) {
        c = *it;
        done = (c < wcsp->getUb());
        ++it;
    }
    if (!done) {
        c = max(UNIT_COST, wcsp->getUb() - 1);
    }
    itThreshold = c;
}

void VACExtension::iniThreshold(Cost threshold)
{
    itThreshold = threshold;
}

void VACExtension::nextScaleCost()
{
    Cost c = MAX_COST;
    bool done = false;
    list<Cost>::iterator it = scaleVAC.begin();
    while ((it != scaleVAC.end()) && !done) {
        c = *it;
        done = c < itThreshold;
        ++it;
    }
    if (!done)
        c = itThreshold / (UNIT_COST + UNIT_COST);

    if (Store::getDepth() <= 1) {
        if (c < ToulBar2::costThresholdPre)
            c = MIN_COST;
    } else if (c < ToulBar2::costThreshold)
        c = MIN_COST;

    itThreshold = c;
}

void VACExtension::resetSupports()
{ // TODO: reset TernaryConstraint supports
    for (unsigned int i = 0; i < wcsp->numberOfConstraints(); i++) {
        Constraint* ctr = wcsp->getCtr(i);
        if (ctr->isBinary() && ctr->connected()) {
            ((BinaryConstraint*)ctr)->resetSupports();
        }
    }
    for (int i = 0; i < wcsp->elimBinOrder; i++) {
        Constraint* ctr = wcsp->elimBinConstrs[i];
        if (ctr->isBinary() && ctr->connected()) {
            ((BinaryConstraint*)ctr)->resetSupports();
        }
    }
}

// do not need to revise all variables because we assume soft AC already done
bool VACExtension::enqueueVAC(Cost threshold, Cost previousThreshold)
{
    wcsp->revise(NULL);
    assert(VAC.empty());
    assert(previousThreshold == -1 || previousThreshold > threshold);

#ifdef INCREMENTALVAC
    for (Queue::iterator iter = VAC2.begin(); iter != VAC2.end(); ++iter) {
        VACVariable* x = (VACVariable*)iter.getElt()->content.var;
        x->queueVAC();
    }
#endif
    VACVariable* x;
    TreeDecomposition* td = wcsp->getTreeDec();
    int bucket2 = cost2log2gub(previousThreshold);
    if (bucket2 < 0 || bucket2 >= wcsp->getNCBucketSize())
        bucket2 = wcsp->getNCBucketSize() - 1;
    int bucket = cost2log2glb(threshold);
    if (bucket < 0)
        bucket = 0;
    for (; bucket <= bucket2; bucket++) {
        VariableList* varlist = wcsp->getNCBucket(bucket);
        for (VariableList::iterator iter = varlist->begin(); iter != varlist->end();) {
            x = (VACVariable*)*iter;
            if (!(x->isNull(x->getMaxCost()))) { // SdG: do not check if x is assigned because it might contain a single value in its domain with a nonzero cost (smaller than the previous threshold)
                if (td) {
                    if (td->isActiveAndInCurrentClusterSubTree(x->getCluster())) {
                        x->queueVAC();
#ifdef INCREMENTALVAC
                        x->queueVAC2();
#endif
                    }
                } else {
                    x->queueVAC();
#ifdef INCREMENTALVAC
                    x->queueVAC2();
#endif
                }
            }
            ++iter;
        }
    }
    return !VAC.empty();
}

bool VACExtension::propagate()
{
    if (!ToulBar2::RASPS && Store::getDepth() >= abs(ToulBar2::vac)) {
        inconsistentVariable = -1;
        return false;
    }

    if (firstTime()) {
        init();
        if (ToulBar2::verbose >= 1)
            cout << "Dual bound before VAC: " << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->getDDualBound() << std::setprecision(DECIMAL_POINT) << endl;
    }

    bool util = false;

    breakCycles = 0;

    acSupport.clear();
    bool acSupportOK = false;
    if (ToulBar2::VAClin) {
        for (KnapsackList::iterator iter = wcsp->knapsackList.begin(); iter != wcsp->knapsackList.end(); ++iter) {
            if ((*iter)->connected()) {
                (*iter)->InitVac();
            }
        }
    }
    while (!util && itThreshold != MIN_COST) {
        int storedepth = Store::getDepth();
        bool isvac = true;
#ifdef AC2001
        resetSupports();
#endif
#ifdef INCREMENTALVAC
        prevItThreshold = -1;
        clear();
        if (!CSP(wcsp->getLb(), wcsp->getUb())) {
            Store::store();
        }
#endif
        while (isvac && itThreshold != MIN_COST) {
#ifndef INCREMENTALVAC
            prevItThreshold = -1;
            clear();
            if (!CSP(wcsp->getLb(), wcsp->getUb())) {
                Store::store();
            }
            setNbAssigned(nbvars - nbunassigned);
#endif
            minlambda = wcsp->getUb() - wcsp->getLb();
            inconsistentVariable = -1;
            PBconflict = -1;
            nbIterations++;
            assert(VAC.empty());
            assert(wcsp->NC.empty());
            // skip this current itThreshold if there are no new value removals do be done compared to previous itThreshold
            while (!enqueueVAC(itThreshold, prevItThreshold) && itThreshold != MIN_COST) {
                prevItThreshold = itThreshold;
                nextScaleCost();
            }
            if (itThreshold == MIN_COST) {
                Store::restore(storedepth);
                inconsistentVariable = -1;
                return false;
            }
            enforcePass1();
            isvac = isVAC();
            if (PBconflict != -1)
                isvac = false;
            // assert(!isvac || checkPass1());
            if (!isvac && CSP(wcsp->getLb(), wcsp->getUb())) {
                if (ToulBar2::weightedDegree)
                    wcsp->conflict();
                Store::restore(storedepth);
                throw Contradiction();
            }
            if (ToulBar2::vacValueHeuristic && isvac) {
                acSupportOK = true;
                acSupport.clear();
                // fill SeekSupport with ALL variables if in preprocessing (i.e. before the search)
                if (Store::getDepth() <= 1 || ToulBar2::debug) {
                    for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
                        ((VACVariable*)wcsp->getVar(i))->queueSeekSupport();
                    }
                }
                // remember first arc consistent domain values in Bool(P) before restoring domains and check VAC-integrality
                while (!SeekSupport.empty()) {
                    VACVariable* x = (VACVariable*)SeekSupport.pop();
                    bool vacintegral = ToulBar2::FullEAC && x->getDomainSize() == 1;
                    if (vacintegral && (x->getInf() != x->getSupport() || !x->isFullEAC())) {
                        acSupport.push_back(make_tuple(x, x->getInf(), true));
                    } else if (!vacintegral && x->cannotbe(x->getSupport())) {
                        Value bestValue = wcsp->getBestValue(x->wcspIndex);
                        if (x->canbe(bestValue) && x->getCost(bestValue) == MIN_COST) { // favor solution-based phase saving
                            acSupport.push_back(make_tuple(x, bestValue, false));
                        } else {
                            for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
                                if (x->getCost(*iterX) == MIN_COST) {
                                    acSupport.push_back(make_tuple(x, *iterX, false));
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            if (isvac) {
                ToulBar2::nbTimesIsVAC++;
                if (itThreshold > 1)
                    ToulBar2::nbTimesIsVACitThresholdMoreThanOne++;

                if (ToulBar2::verbose >= 7) {
                    cout << "End of VAC pass 1..." << endl;
                    cout << *wcsp;
                }

                if (ToulBar2::RASPSsaveitThresholds) {
                    ToulBar2::RASPSnbStrictACVariables = wcsp->numberOfVariables() - wcsp->numberOfUnassignedVariables();
                    double ratio = (ToulBar2::RASPSnbStrictACVariables == 0) ? 0.0000000001 : (((double)ToulBar2::RASPSnbStrictACVariables / (double)wcsp->numberOfVariables()) / (double)itThreshold);
                    if (ToulBar2::verbose >= 0) {
                        cout << std::fixed << std::setprecision(7);
                        cout << "Threshold: " << itThreshold << " NbAssignedVar: " << ToulBar2::RASPSnbStrictACVariables << " Ratio: " << ratio << " SumOfDomainSize: " << wcsp->getDomainSizeSum() << endl;
                        cout << std::fixed << std::setprecision(DECIMAL_POINT);
                    }
                    ToulBar2::RASPSitThresholds.push_back(std::make_pair(itThreshold, ratio));
                }

                if (ToulBar2::RASPS) {
                    ToulBar2::RASPS = false;
                    assert(ToulBar2::restart < 1);
                    Store::store();

                    Solver* solver = (Solver*)wcsp->getSolver();
                    Long storehbfs = ToulBar2::hbfs;
                    Long storehbfsGlobalLimit = ToulBar2::hbfsGlobalLimit;
                    Long storehbfsLimit = solver->hbfsLimit;
                    Long storenbBacktracksLimit = solver->nbBacktracksLimit;
                    Long storenbBacktracks = solver->nbBacktracks;
                    Long storerestart = ToulBar2::restart;
                    bool storeLimited = ToulBar2::limited;
                    int storeVac = ToulBar2::vac;
                    bool storeFullEAC = ToulBar2::FullEAC;

                    ToulBar2::hbfs = 0;
                    ToulBar2::hbfsGlobalLimit = 0;
                    solver->hbfsLimit = LONGLONG_MAX;
                    solver->nbBacktracksLimit = solver->nbBacktracks + ToulBar2::RASPSnbBacktracks;
                    if (ToulBar2::useRASPS <= 1) {
                        ToulBar2::restart = 1; // random variable selection for breaking ties in DFS
                    } else {
                        ToulBar2::restart = -1; // no randomness for LDS
                    }
                    ToulBar2::limited = false;
                    ToulBar2::vac = 0;
                    ToulBar2::FullEAC = false;

                    Cost lastUB = wcsp->getUb();
                    try {
                        try {
                            // Current WCSP is AC(Bool(P))
                            vector<int> variables;
                            vector<Value> values;

                            for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
                                if (wcsp->getVar(i)->enumerated()) {
                                    EnumeratedVariable* xi = (EnumeratedVariable*)wcsp->getVar(i);
                                    if (xi->getInf() == xi->getSup()) {
                                        variables.push_back(i);
                                        values.push_back(xi->getInf());
                                    } else {
                                        xi->queueDAC();
                                        xi->queueEAC1();
                                        xi->queueAC();
                                    }
                                    xi->findSupport(); // update support values
                                    xi->propagateNC(); // and update maxcost values before propagate done in assignLS
                                }
                            }
                            if (variables.size() > 0)
                                wcsp->assignLS(variables, values, true); // option true: make sure already assigned variables are removed from Solver::unassignedVars
                            wcsp->propagate();
                            if (ToulBar2::useRASPS <= 1) {
                                // cout << "call to recursiveSolve from VAC" << endl;
                                solver->recursiveSolve(wcsp->getLb()); // if a new solution is found, UB will be updated automatically
                            } else {
                                // cout << "call to recursiveSolveLDS from VAC" << endl;
                                solver->recursiveSolveLDS(ToulBar2::useRASPS - 1);
                            }
                        } catch (const Contradiction&) {
                            wcsp->whenContradiction();
                        }
                    } catch (const NbBacktracksOut&) {
                    }

                    ToulBar2::hbfs = storehbfs;
                    ToulBar2::hbfsGlobalLimit = storehbfsGlobalLimit;
                    solver->hbfsLimit = storehbfsLimit;
                    solver->nbBacktracksLimit = storenbBacktracksLimit;
                    if (solver->nbBacktracksLimit < LONGLONG_MAX)
                        solver->nbBacktracksLimit += solver->nbBacktracks - storenbBacktracks;
                    ToulBar2::restart = storerestart;
                    ToulBar2::limited = storeLimited; // still a complete search
                    ToulBar2::vac = storeVac;
                    ToulBar2::FullEAC = storeFullEAC;

                    Store::restore(storedepth);
                    inconsistentVariable = -1;
                    return (wcsp->getUb() < lastUB);
                }
            }

            if (isvac) {
                prevItThreshold = itThreshold;
                nextScaleCost();
#ifdef INCREMENTALVAC
                if (itThreshold != MIN_COST) {
                    while (!wcsp->NC.empty()) { // update maxCost and maxCostValue after value removals by AC2001
                        EnumeratedVariable* x = (EnumeratedVariable*)wcsp->NC.pop();
                        Cost maxcost = MIN_COST;
                        Value maxcostvalue = x->getSup() + 1;
                        // Warning! the first value must be visited because it may be removed
                        for (EnumeratedVariable::iterator iter = x->begin(); iter != x->end(); ++iter) {
                            Cost cost = x->getCost(*iter);
                            if (LUB(&maxcost, cost) || x->cannotbe(maxcostvalue)) {
                                maxcostvalue = *iter;
                            }
                        }
                        x->setMaxCostValue(maxcostvalue);
                        if (maxcost != x->getMaxCost()) {
                            assert(maxcost < x->getMaxCost());
                            x->setMaxCost(maxcost);
                            int newbucket = min(cost2log2gub(maxcost), wcsp->getNCBucketSize() - 1);
                            x->changeNCBucket(newbucket);
                        }
                    }
                }
#endif
                //                if (ToulBar2::preprocessTernaryRPC && Store::getDepth() <= 1) {
                //                    double mintight = 1e20;
                //                    Variable *var1 = NULL;
                //                    Variable *var2 = NULL;
                //                    for (unsigned int i = 0; i < wcsp->numberOfConstraints(); i++) {
                //                        if (wcsp->getCtr(i)->connected() && wcsp->getCtr(i)->isBinary()) {
                //                            wcsp->getCtr(i)->resetTightness();
                //                            double tight = wcsp->getCtr(i)->getTightness();
                //                            if (mintight > tight) {
                //                                mintight = tight;
                //                                var1 = wcsp->getCtr(i)->getVar(0);
                //                                var2 = wcsp->getCtr(i)->getVar(1);
                //                            }
                //                        }
                //                    }
                //                    for (int i = 0; i < wcsp->getElimBinOrder(); i++) {
                //                        if (wcsp->getCtr(-i -1)->connected()) {
                //                            wcsp->getCtr(-i -1)->resetTightness();
                //                            double tight = wcsp->getCtr(-i -1)->getTightness();
                //                            if (mintight > tight) {
                //                                mintight = tight;
                //                                var1 = wcsp->getCtr(-i -1)->getVar(0);
                //                                var2 = wcsp->getCtr(-i -1)->getVar(1);
                //                            }
                //                            mintight = MIN(mintight, tight);
                //                        }
                //                    }
                //                    if (ToulBar2::verbose >= 0 && var1 && var2) {
                //                        cout << "Minimum binary cost function tightness after VAC: " << mintight << " on variables " << var1->getName() << " and " << var2->getName() << endl;
                //                    }
                //                }
            }
#ifndef INCREMENTALVAC
            Store::restore(storedepth);

#endif
        }
#ifdef INCREMENTALVAC
        Store::restore(storedepth);
#endif
        if (!isvac) {
            if (ToulBar2::VAClin) {
                for (KnapsackList::iterator iter = wcsp->knapsackList.begin(); iter != wcsp->knapsackList.end(); ++iter) {
                    if ((*iter)->connected()) {
                        (*iter)->ResetVACLastValTested();
                    }
                }
            }
            enforcePass2();
            if (ToulBar2::verbose > 0)
                cout << "VAC dual bound: " << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->getDDualBound() << std::setprecision(DECIMAL_POINT) << "    incvar: " << inconsistentVariable << "    minlambda: " << minlambda << "      itThreshold: " << itThreshold << endl;
            if (CUT(wcsp->getLb() + minlambda, wcsp->getUb())) {
                if (ToulBar2::weightedDegree)
                    wcsp->conflict();
                throw Contradiction();
            }
            util = enforcePass3();
            for (KnapsackList::iterator iter = wcsp->knapsackList.begin(); iter != wcsp->knapsackList.end(); ++iter) {
                if ((*iter)->connected()) {
                    (*iter)->InitVac();
                }
            }
        } else {
            nextScaleCost();
        }
    }

    if (ToulBar2::vacValueHeuristic && acSupportOK && isVAC()) {
        // update current unary support if possible && needed
        for (vector<tuple<VACVariable*, Value, bool>>::iterator iter = acSupport.begin(); iter != acSupport.end(); ++iter) {
            VACVariable* x = std::get<0>(*iter);
            Value val = std::get<1>(*iter);
            bool vacintegral = std::get<2>(*iter);
            if (x->canbe(val)) {
                if (x->getCost(val) == MIN_COST) {
                    if (ToulBar2::verbose > 0 && (x->getSupport() != val || (vacintegral && !x->isFullEAC())))
                        cout << "CHANGE SUPPORT " << x->getName() << " from " << x->getSupport() << ((x->isFullEAC()) ? "!" : "") << " to " << val << ((vacintegral) ? "!" : "") << endl;
                    if (vacintegral && !x->isFullEAC()) {
                        x->setFullEAC(); // TODO: is it better to set VAC-integrality to true even if current unary cost is not zero?
#ifndef NDEBUG
                        x->queueFEAC(); // TODO: is it better to set VAC-integrality to true even if some cost functions are violated?
#endif
                    }
                    x->setSupport(val);
                }
            } else {
                if (ToulBar2::verbose > 0)
                    cout << "WARNING: BAD EAC SUPPORT " << val << " FOR VARIABLE " << *x << endl;
            }
        }
    }

    return util;
}

bool VACExtension::enforcePass1(VACVariable* xj, VACBinaryConstraint* cij)
{
    bool wipeout = false;
    VACVariable* xi;
    xi = (VACVariable*)cij->getVarDiffFrom(xj);
    for (EnumeratedVariable::iterator it = xi->begin(); it != xi->end(); ++it) {
        Value v = *it;
        if (xi->getVACCost(v) > MIN_COST) {
            xi->removeVAC(v);
            xi->setPBkillers(v,{});
        } else if (cij->revise(xi, v)) {
            wipeout = xi->removeVAC(v);
            xi->setKiller(v, xj->wcspIndex);
            xi->setPBkillers(v, {});
            queueP->push(pair<int, Value>(xi->wcspIndex, v));
            if (wipeout) {
                inconsistentVariable = xi->wcspIndex;
                return true;
            }
            xi->queueVAC();
#ifdef INCREMENTALVAC
            xi->queueVAC2();
#endif
            if (ToulBar2::vacValueHeuristic)
                xi->queueSeekSupport();
        }
    }
    return false;
}

void VACExtension::enforcePass1()
{
    bool firstit = true;
    if (ToulBar2::verbose >= 4)
        cout << "enforcePass1 itThreshold: " << itThreshold << endl;
    PBkillersctr.clear();
    while (!VAC.empty()) {
        if (ToulBar2::interrupted && !ToulBar2::isZ)
            throw TimeOut();
        VACVariable* xj = (VACVariable*)VAC.pop_first();
        for (EnumeratedVariable::iterator it = xj->begin(); it != xj->end(); ++it) {
            if (xj->getVACCost(*it) > MIN_COST){
                bool wipeout = xj->removeVAC(*it);
                xj->setPBkillers(*it,{});
                if (wipeout) {
                    inconsistentVariable = xj->wcspIndex;
                    return;
                }
            }else if(xj->canbe(*it)){
                xj->setPBkillers(*it,{});
            }
        }
        for (ConstraintList::iterator itc = xj->getConstrs()->begin();
             itc != xj->getConstrs()->end(); ++itc) {
            Constraint* c = (*itc).constr;
            if (c->isBinary() && !c->isDuplicate()) {
                VACBinaryConstraint* cij = (VACBinaryConstraint*)c;
                if (enforcePass1(xj, cij))
                    return;
            }
        }
        //Apply phase 1 on all the connected knapsack constraints after propagating all binary cost functions
        if (VAC.empty() && ToulBar2::VAClin) {
            for (KnapsackList::iterator iter = wcsp->knapsackList.begin(); iter != wcsp->knapsackList.end(); ++iter) {
                KnapsackConstraint* k = (*iter);
                if (k->connected() && (k->VACNeedPropagate() || firstit)) {
                    bool wipeout;
                    VACVariable* xi;
                    Cost OPT = k->VACPass1(killers, killed, wcsp->getUb(), itThreshold);
                    //OPT=-1 Only if the constraint is infeasible
                    if (OPT == -1) {
                        PBconflict = killers[0].first;
                        //assert(killers.size() > 1);
                        PBkillersctr = {killers.begin() + 1, killers.end()};
                        return;
                    }
                    if (OPT > MIN_COST) {
                        assert(OPT>=itThreshold);
                        //assert(killers.size() > 1);
                        PBconflict = killers[0].first;
                        PBkillersctr = {killers.begin() + 1, killers.end()};
                        minlambda = min(minlambda, OPT);
                        k->IncreasekConstraintVAC(-k->getkConstraintVAC()); // set to zero counter kConstraintVAC
                        k->IncreasekConstraintVAC(1);
                        return;
                    }
                    for (int j = 0; j < (int)killed.size(); ++j) {
                        xi = (VACVariable*)wcsp->getVar(killed[j].first);
                        xi->setPBkillers(killed[j].second, killers);
                        wipeout = xi->removeVAC(killed[j].second);
                        queueP->push(pair<int, Value>(xi->wcspIndex, killed[j].second));
                        if (wipeout) {
                            inconsistentVariable = xi->wcspIndex;
                            return;
                        }
                        xi->queueVAC();
#ifdef INCREMENTALVAC
                        xi->queueVAC2();
#endif
                        if (ToulBar2::vacValueHeuristic)
                            xi->queueSeekSupport();
                    }
                }
            }
            firstit=false;
        }
    }
    inconsistentVariable = -1;
}

bool VACExtension::checkPass1() const
{
    VACBinaryConstraint* cij;
    VACVariable *xi, *xj;
    bool supportFound;
    TreeDecomposition* td = wcsp->getTreeDec();

    for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
        xi = (VACVariable*)wcsp->getVar(i);
        if (td && !td->isActiveAndInCurrentClusterSubTree(xi->getCluster()))
            continue;
        for (ConstraintList::iterator iter = xi->getConstrs()->begin();
             iter != xj->getConstrs()->end(); ++iter) {
            Constraint* c = (*iter).constr;
            if (c->isBinary() && !c->isDuplicate()) {
                cij = (VACBinaryConstraint*)c;
                xj = (VACVariable*)cij->getVarDiffFrom(xi);
                for (EnumeratedVariable::iterator iti = xi->begin(); iti != xi->end(); ++iti) {
                    Value v = *iti;
                    supportFound = false;
                    for (EnumeratedVariable::iterator itj = xj->begin(); itj != xj->end();
                         ++itj) {
                        Value w = *itj;
                        if ((xi->getVACCost(v) == MIN_COST)
                            && (xj->getVACCost(w) == MIN_COST)
                            && (cij->getVACCost(xi, xj, v, w) == MIN_COST)) {
                            supportFound = true;
                            break;
                        }
                    }
                    if (!supportFound) {
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

void VACExtension::enforcePass2()
{
    int i0 = inconsistentVariable;
    int i, j;
    VACVariable *xi0, *xi, *xj;
    VACBinaryConstraint* cij;
    Cost tmplambda;
    Value v;
    int maxk = 1;
    // if (ToulBar2::verbose > 0)  cout << "VAC Enforce Pass 2" << endl;
    if (PBconflict == -1) {
        assert(i0 >= 0);
        xi0 = (VACVariable*)wcsp->getVar(i0);
        for (EnumeratedVariable::iterator iti0 = xi0->begin();
             iti0 != xi0->end(); ++iti0) {
            v = *iti0;
            xi0->addToK(v, 1, nbIterations);
            Cost cost = xi0->getVACCost(v);
            if (cost > MIN_COST) {
                if (cost < minlambda) {
                    minlambda = cost; // NB: we don't need to check for bottleneck here as k=1 necessarily
                }
            } else {
                xi0->setMark(v, nbIterations);
            }
        }
    } else {
        //Mark all the values in the explanation
        for (int k = 0; k < (int)PBkillersctr.size(); ++k) {
            xi = (VACVariable*)wcsp->getVar(PBkillersctr[k].first);
            if(xi->canbe(PBkillersctr[k].second)){
                xi->addToK(PBkillersctr[k].second, 1, nbIterations);
                Cost cost = xi->getVACCost(PBkillersctr[k].second);
                if (cost > MIN_COST) {
                    if (cost < minlambda) {
                        minlambda = cost; // NB: we don't need to check for bottleneck here as k=1 necessarily
                    }
                } else {
                    assert(xi->canbe(PBkillersctr[k].second));
                    xi->setMark(PBkillersctr[k].second, nbIterations);
                }
            }
        }
    }
    while (!queueP->empty()) {
        i = queueP->top().first;
        v = queueP->top().second;
        queueP->pop();
        xi = (VACVariable*)wcsp->getVar(i);
        assert(xi->canbe(v));
        if (xi->isMarked(v, nbIterations) && minlambda >= UNIT_COST) {
            if (xi->getPBkillers(v).empty()) {
                j = xi->getKiller(v);
                xj = (VACVariable*)wcsp->getVar(j);
                queueR->push(pair<int, Value>(i, v));
                cij = (VACBinaryConstraint*)xi->getConstrNotDuplicate(xj);
                assert(cij);
                assert(!cij->isDuplicate());
                // if (ToulBar2::verbose > 6) cout << "x" << xi->wcspIndex << "," << v << "   killer: " << xj->wcspIndex << endl;

                for (EnumeratedVariable::iterator itj = xj->begin(); itj != xj->end(); ++itj) {
                    Value w = *itj;
                    Cost costij = cij->getVACCost(xi, xj, v, w);
                    if (costij > MIN_COST) {
                        int tmpK = xi->getK(v, nbIterations);
                        if (xj->getKiller(w) == i && xj->isMarked(w, nbIterations))
                            tmpK += xj->getK(w, nbIterations);
                        if (!CUT(wcsp->getLb() + costij,
                                wcsp->getUb())) {
                            if ((costij / tmpK) < minlambda) {
                                minlambda = costij / tmpK;
                                if (minlambda < UNIT_COST) { // A cost function bottleneck here !
                                    bneckCost = costij;
                                    bneckCF = cij;
                                    bneckVar = -1;
                                }
                            }
                        } else {
                            if ((costij / tmpK) < minlambda) { // costij should be made infinite to avoid to decrease minlambda
                                Cost cost = (((Double)tmpK * minlambda < MAX_COST)?(tmpK * minlambda):MAX_COST) - costij;
                                assert(cost > MIN_COST);
                                assert(ToulBar2::verbose < 1 || ((cout << "inflate(C" << cij->getVar(0)->getName() << "," << cij->getVar(1)->getName() << ", (" << ((xi == cij->getVar(0)) ? v : w) << "," << ((xi == cij->getVar(0)) ? w : v) << "), " << cost << ")" << endl), true));
                                cij->addcost(xi, xj, v, w, cost);
                            }
                        }
                    } else {
                        int tmpK = xi->getK(v, nbIterations) - cij->getK(xj, w, nbIterations);
                        if (tmpK > 0) {
                            xj->addToK(w, tmpK, nbIterations);
                            if (xj->getK(w, nbIterations) > maxk)
                                maxk = xj->getK(w, nbIterations);
                            cij->setK(xj, w, xi->getK(v, nbIterations), nbIterations);
                            Cost cost = xj->getVACCost(w);
                            if (cost == MIN_COST) {
                                xj->setMark(w, nbIterations);
                            } else { // we assume NC has been done before
                                assert(!CUT(wcsp->getLb() + cost, wcsp->getUb()));
                                tmplambda = cost / xj->getK(w, nbIterations);
                                xj->setPBkillers(w,{});
                                if (tmplambda < minlambda) {
                                    minlambda = tmplambda;
                                    if (minlambda < UNIT_COST) { // A unary cost bottleneck here
                                        bneckVar = j;
                                        bneckCF = NULL;
                                        bneckCost = cost;
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                pair<int, Value> PBKILLERSxi0 = xi->getPBkillers(v)[0];
                queueR->push(pair<int, Value>(i, v));
                assert(wcsp->constrs[PBKILLERSxi0.first]->isKnapsack());
                KnapsackConstraint* knap = (KnapsackConstraint*)wcsp->constrs[PBKILLERSxi0.first];
                int usek = 0;
                // alreadysendk is only useful for the values in NotVarVal. It captures the maximal number of quantum requested by the already processed values in NotVarVal.
                int alreadysendk=0; 
                Cost OPT = knap->VACPass2(PBKILLERSxi0.second, { i, v }, killers, wcsp->getUb(), xi->getK(v, nbIterations));
                vector<Value> wasLastVal = knap->waslastValue(i,v);
                usek = xi->getK(v, nbIterations);
                if (wasLastVal.size()>1) {
                    for (int k = 0; k < (int)wasLastVal.size(); ++k) {
                        if (wasLastVal[k]!=v && knap->getVACLastVALChecked(xi->wcspIndex,wasLastVal[k]) && alreadysendk < xi->getK(wasLastVal[k], nbIterations)) {
                            alreadysendk = xi->getK(wasLastVal[k], nbIterations);
                        }
                    }
                    usek=max(0,usek-alreadysendk);
                }
                knap->IncreasekConstraintVAC(usek);
                xi->setPBkillers(v, killers);
                assert(!killers.empty());
                Cost tmplambda = OPT / max(1, knap->getkConstraintVAC());
                if (tmplambda < minlambda)
                    minlambda = tmplambda;
                if (minlambda < UNIT_COST) { // A unary cost bottleneck here
                    bneckVar = i;
                    bneckCF = NULL;
                    bneckCost = OPT;
                } else if(usek>0) {
                    knap->RestVACGroupExt();
                    for (int k = 1; k < (int)killers.size(); ++k) {
                        xj = (VACVariable*)wcsp->getVar(killers[k].first);
                        if (xj->canbe(killers[k].second)) {
                            Value w = killers[k].second;
                            xj->addToK(w, usek, nbIterations);
                            if (xj->getK(w, nbIterations) > maxk)
                                maxk = xj->getK(w, nbIterations);
                            knap->setkVAC(xj, w, knap->getkVAC(xj, w, nbIterations) + usek, nbIterations);
                            Cost cost = xj->getVACCost(w);
                            if (cost == MIN_COST)
                                xj->setMark(w, nbIterations);
                            else { // we assume NC has been done before
                                assert(!CUT(wcsp->getLb() + cost, wcsp->getUb()));
                                tmplambda = cost / xj->getK(w, nbIterations);
                                xj->setPBkillers(w,{});
                                if (tmplambda < minlambda) {
                                    minlambda = tmplambda;
                                    if (minlambda < UNIT_COST) { // A unary cost bottleneck here
                                        bneckVar = killers[k].first;
                                        bneckCF = NULL;
                                        bneckCost = cost;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    // if (maxK == 0) {
    //   maxK = wcsp->getUb() - wcsp->getLb();
    // }
    if (ToulBar2::verbose > 1)
        cout << "minLambda: " << minlambda << "\t\t (dualb = " << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->getDDualBound() << ", primalb = " << wcsp->getDPrimalBound() << ")" << std::setprecision(DECIMAL_POINT) << endl;
}

bool VACExtension::enforcePass3()
{
    bool util = (minlambda >= UNIT_COST);
    /*if(util) {
           Cost ub = wcsp->getUb();
           Cost lb = wcsp->getLb();
           util = ( (ub - lb)/(10000*ToulBar2::costMultiplier) ) < minlambda;
           } */
    Cost lambda = minlambda;

    // if (ToulBar2::verbose > 2) cout << "VAC Enforce Pass 3.   minlambda " << minlambda << " , var: " << inconsistentVariable << endl;

    VACVariable* xi0 = (PBconflict == -1) ? (VACVariable*)wcsp->getVar(inconsistentVariable) : NULL;

    int maxk = 0;

    if (!util) { // Empty R ?
        assert(bneckVar != -1 || bneckCF != NULL);
        if (bneckVar != -1) {
            VACVariable *xi = (VACVariable*)wcsp->getVar(bneckVar);
            xi->setThreshold(bneckCost + UNIT_COST);
        } else if (bneckCF != NULL) {
            bneckCF->setThreshold(bneckCost + UNIT_COST);
        }
        breakCycles++;
        if (ToulBar2::verbose > 1) {
            cout << "BreakCycle: bneckCost=" << bneckCost + UNIT_COST << ", bneckVar=" << bneckVar << ", bneckCF=";
            if (bneckCF)
                cout << *bneckCF;
            else
                cout << bneckCF;
            cout << endl;
        }
        if (breakCycles >= 5) {
            if (ToulBar2::verbose > 1)
                cout << "BreakCycle stops!" << endl;
            inconsistentVariable = -1;
            itThreshold = MIN_COST;
        }
        return false;
    }
    // update general stats
    nlb++;
    sumlb += lambda;
    sumvars += queueR->size();
    maxk = 0;
    tempvars.clear();
    while (!queueR->empty()) {
        int j = queueR->top().first;
        Value w = queueR->top().second;
        queueR->pop();
        VACVariable *xj = (VACVariable*)wcsp->getVar(j);
        if (xj->getPBkillers(w).empty()) {
            int i = xj->getKiller(w);
            VACVariable *xi = (VACVariable*)wcsp->getVar(i);
            VACBinaryConstraint *cij = (VACBinaryConstraint*)xi->getConstrNotDuplicate(xj);
            assert(cij);
            assert(!cij->isDuplicate());
            int xjk = xj->getK(w, nbIterations);
            if (maxk < xjk)
                maxk = xjk;
            for (EnumeratedVariable::iterator iti = xi->begin(); iti != xi->end(); ++iti) {
                Value v = *iti;
                if (cij->getK(xi, v, nbIterations) != 0) {
                    int tmpK = cij->getK(xi, v, nbIterations);
                    Cost ecost = (((Double)tmpK * lambda < MAX_COST)?(tmpK * lambda):MAX_COST);
                    cij->setK(xi, v, 0, nbIterations);
                    cij->VACextend(xi, v, ecost);
                    // extention from unary to binary cost function may break soft AC/DAC in both directions due to isNull/itThreshold
                    if (ToulBar2::LcLevel == LC_AC) {
                        xi->queueAC();
                        xj->queueAC();
                    } else {
                        if (cij->getDACScopeIndex() == cij->getIndex(xi)) {
                            xi->queueAC();
                            xi->queueEAC1();
                            xj->queueDAC();
                        } else {
                            xi->queueDAC();
                            xj->queueAC();
                            xj->queueEAC1();
                        }
                    }
                }
            }
            int tmpK = xj->getK(w, nbIterations);
            cij->VACproject(xj, w, ((Double)tmpK * lambda < MAX_COST)?(tmpK * lambda):MAX_COST);
            tempvars.insert(xj->wcspIndex);
        } else {
            const vector<pair<int, Value>>& PBKILLERSxj = xj->getPBkillers(w);
            assert(wcsp->constrs[PBKILLERSxj[0].first]->isKnapsack());
            KnapsackConstraint* knap = (KnapsackConstraint *)wcsp->constrs[PBKILLERSxj[0].first];
            int xjk = xj->getK(w, nbIterations);
            if (maxk < xjk)
                maxk = xjk;
            knap->RestVACGroupExt();
            for (int k = 1; k < (int)PBKILLERSxj.size(); ++k) {
                VACVariable *xi = (VACVariable*)wcsp->getVar(PBKILLERSxj[k].first);
                if (xi->canbe(PBKILLERSxj[k].second) && knap->getkVAC(xi, PBKILLERSxj[k].second, nbIterations) != 0) {
                    int tmpK = knap->getkVAC(xi, PBKILLERSxj[k].second, nbIterations);
                    Cost ecost = ((Double)tmpK * lambda < MAX_COST)?(tmpK * lambda):MAX_COST;
                    knap->VACextend(xi, PBKILLERSxj[k].second, ecost);
                    if (ToulBar2::LcLevel == LC_AC) { // SdG: to be compatible with verify for knapsack constraints
                        xi->queueAC();
                    }
                    vector<Value> LastValues = knap->waslastValue(xi->wcspIndex, PBKILLERSxj[k].second);
                    for (int l = 0; l < (int)LastValues.size(); ++l) {
                        knap->setkVAC(xi, LastValues[l], 0, nbIterations);
                        assert(xi->getCost(LastValues[l])>= MIN_COST);
                    }
                    assert(knap->getkVAC(xi, PBKILLERSxj[k].second, nbIterations)>=0);
                }
            }
            vector<Value> LastValues = knap->waslastValue(xj->wcspIndex,w);
            for (int k = 0; k < (int)LastValues.size(); ++k) {
                if(!xj->getPBkillers(LastValues[k]).empty() && xj->getPBkillers(LastValues[k])[0].first==PBKILLERSxj[0].first){
                    assert(xj->canbe(LastValues[k]));
                    if(xjk < xj->getK(LastValues[k],nbIterations)){
                        xjk = xj->getK(LastValues[k],nbIterations);
                    }
                    xj->setK(LastValues[k], 0, nbIterations);
                }
            }
            knap->VACproject(xj, w, ((Double)xjk * lambda < MAX_COST)?(xjk * lambda):MAX_COST); // SdG: VACproject may project more than what is needed after, possibly breaking (EAC) supports and DAC
            tempvars.insert(xj->wcspIndex);
        }
    }
    sumk += maxk;
    if (maxk > theMaxK)
        theMaxK = maxk;
    if (PBconflict == -1) {
        xi0->extendAll(lambda);
        xi0->projectLB(lambda);
    } else {
        assert(wcsp->constrs[PBconflict]->isKnapsack());
        TreeDecomposition* td = wcsp->getTreeDec();
        KnapsackConstraint* k2 = (KnapsackConstraint *)wcsp->constrs[PBconflict];
        k2->RestVACGroupExt();
        k2->VACPass3(EPT, minlambda, PBkillersctr);
#ifndef NDEBUG
        sort(EPT.begin(), EPT.end(), [&](auto& x, auto& y) {return x.second < y.second;}); // SdG: performs all projections before extensions in order to keep positive unary costs inside the sequence of EPTs
#endif
        for (int l = 0; l < (int)EPT.size(); ++l) {
            VACVariable *xi = (VACVariable*)wcsp->getVar(EPT[l].first.first);
            if (EPT[l].second > MIN_COST) {
                if (td && xi->canbe(EPT[l].first.second))
                    td->addDelta(k2->getCluster(), xi, EPT[l].first.second, -EPT[l].second);
                xi->extend(EPT[l].first.second, EPT[l].second);
                assert(xi->getCost(EPT[l].first.second)>= MIN_COST);
            } else {
                if (td && xi->canbe(EPT[l].first.second))
                    td->addDelta(k2->getCluster(), xi, EPT[l].first.second, -EPT[l].second);
                xi->project(EPT[l].first.second, -EPT[l].second);
                tempvars.insert(xi->wcspIndex);
            }
        }
        if (k2->connected()) {
            k2->VACObjConsistency(tempvars);
        }
    }
    for (int i : tempvars) {
        ((EnumeratedVariable *)wcsp->getVar(i))->findSupport();
    }

    return true;
}

// void VACExtension::updateStat(Cost lambda)
// {
//   //tVACStat* v = heapAccess[inconsistentVariable];
//   if(varAssign >= 0) {
//    tVACStat* v = heapAccess[varAssign];
//    v->sumlb += lambda;
//    v->nlb++;
//   }
// }

// Cost VACExtension::getVarCostStat( int i )
// {
//   tVACStat* v = heap[i];
//   return v->sumlb;
// }

// Long VACExtension::getVarTimesStat( int i )
// {
//   if(i < 0) return 0;
//   tVACStat* v = heap[i];
//   return v->nlb;
// }

void VACExtension::iniSingleton()
{
    singletonI.clear();
    for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
        int size = wcsp->getDomainSize(i);
        for (int a = 0; a < size; a++)
            singletonI.insert(wcsp->getMaxDomainSize() * i + a);
    }
}

void VACExtension::updateSingleton()
{
    set<Long>& s1 = singleton;
    set<Long> s2(singletonI);
    singletonI.clear();
    set_intersection(s1.begin(), s1.end(),
        s2.begin(), s2.end(),
        inserter(singletonI, singletonI.begin()));
    singleton.clear();
}

void VACExtension::removeSingleton()
{
    set<Long>& s = singletonI;
    set<Long>::iterator it = s.begin();
    while (it != s.end()) {
        int ivar = *it / wcsp->getMaxDomainSize();
        Value a = *it % wcsp->getMaxDomainSize();
        Variable* var = wcsp->getVar(ivar);
        var->remove(a);
        var->queueNC();
        ++it;
    }
    wcsp->propagate();
}

void VACExtension::clear()
{
    while (!queueP->empty())
        queueP->pop();
    while (!queueR->empty())
        queueR->pop();
    // Cannot use  VAC.clear() as it will not reset timeStamps which are based on the number of wcsp propagate calls and not VAC iterations
    while (!VAC.empty())
        VAC.pop();
#ifdef INCREMENTALVAC
    while (!wcsp->NC.empty())
        wcsp->NC.pop();
    while (!VAC2.empty())
        VAC2.pop();
#endif
    if (ToulBar2::vacValueHeuristic)
        while (!SeekSupport.empty())
            SeekSupport.pop();
}

void VACExtension::queueVAC(DLink<VariableWithTimeStamp>* link)
{
    assert(ToulBar2::vac);
    VAC.push(link, wcsp->getNbNodes());
}
#ifdef INCREMENTALVAC
void VACExtension::queueVAC2(DLink<VariableWithTimeStamp>* link)
{
    assert(ToulBar2::vac);
    VAC2.push(link, wcsp->getNbNodes());
}
#endif

void VACExtension::queueSeekSupport(DLink<VariableWithTimeStamp>* link)
{
    assert(ToulBar2::vac);
    SeekSupport.push(link, wcsp->getNbNodes());
}

void VACExtension::printStat(bool ini)
{
    if (ToulBar2::verbose >= 1) {
        Double mean = to_double(sumlb) / (Double)nlb;
        cout << "VAC mean lb/incr: " << std::setprecision(DECIMAL_POINT) << mean << "     total increments: " << nlb
             << "     cyclesize: " << (double)sumvars / (double)nlb << "     k: " << (double)sumk / (double)nlb << " (mean), " << theMaxK << " (max)" << endl;
    }
    if (ini && nlb > 0 && sumlb > MIN_COST && ToulBar2::verbose >= 1) {
        if (ToulBar2::uai)
            cout << "VAC dual bound: " << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->getDDualBound() << std::setprecision(DECIMAL_POINT) << " energy: " << -(wcsp->Cost2LogProb(wcsp->getLb()) + ToulBar2::markov_log) << " (iter:" << nlb << ")" << endl;
        else
            cout << "VAC dual bound: " << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->getDDualBound() << std::setprecision(DECIMAL_POINT) << " (iter:" << nlb << ")" << endl;
    }
    if (ToulBar2::verbose >= 0) {
        cout << "Number of VAC iterations: " << nbIterations << endl;
        cout << "Number of times is VAC: " << ToulBar2::nbTimesIsVAC << " Number of times isvac and itThreshold > 1: " << ToulBar2::nbTimesIsVACitThresholdMoreThanOne << endl;
    }
    // sort(heap.begin(), heap.end(), cmp_function);
    /*cout << "Vars: ";
           vector<tVACStat*>::iterator it = heap.begin();
           while(it != heap.end()) {
           tVACStat* v = *it;
           if(v->sumlb != MIN_COST) cout << "(" << v->var << "," << v->sumlb << ") ";
           ++it;
           }
           cout << endl; */

    sumk = 0;
    theMaxK = 0;
    sumvars = 0;
    sumlb = MIN_COST;
    nlb = 0;
}

void VACExtension::printTightMatrix()
{
    ofstream ofs("problem.dat");

    Cost Top = wcsp->getUb();
    for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
        for (unsigned int j = 0; j < wcsp->numberOfVariables(); j++) {
            if (i != j) {
                EnumeratedVariable* x = (EnumeratedVariable*)wcsp->getVar(i);
                EnumeratedVariable* y = (EnumeratedVariable*)wcsp->getVar(j);
                Constraint* bctr = x->getConstr(y);
                double t = 0;
                if (bctr)
                    t = bctr->getTightness();
                if (t > to_double(Top))
                    t = to_double(Top);
                t = t * 256.0 / to_double(Top);
                ofs << t << " ";
            } else
                ofs << 0 << " ";
        }
        ofs << endl;
    }
}

/* Min-Sum diffusion algorithm */
void VACExtension::minsumDiffusion()
{
    for (int times = 0; times < 2; times++) {
        bool change = true;
        int maxit = ToulBar2::minsumDiffusion;
        if (ToulBar2::verbose >= 0) {
            cout << "MinSumDiffusion: " << endl;
            cout << "   max iterations " << maxit << endl;
            cout << "   dual bound = " << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->getDDualBound() << std::setprecision(DECIMAL_POINT) << endl;
        }
        int ntimes = 0;
        while (change && (ntimes < maxit)) {
            change = false;
            int nchanged = 0;
            for (unsigned int i = 0; i < wcsp->numberOfVariables();
                 i++)
                if (wcsp->unassigned(i)) {
                    VACVariable* evar = (VACVariable*)wcsp->getVar(i);
                    if (evar->averaging()) {
                        change = true;
                        nchanged++;
                        evar->findSupport();
                    }
                }
            ntimes++;
            // cout << "it " << ntimes << "   changed: " << nchanged << endl;
        }
        if (ToulBar2::verbose >= 0)
            cout << "   done iterations: " << ntimes << endl;
        for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++)
            if (wcsp->unassigned(i)) {
                EnumeratedVariable* evar = (EnumeratedVariable*)wcsp->getVar(i);
                evar->findSupport();
            }
        for (unsigned int i = 0; i < wcsp->numberOfConstraints(); i++)
            if (wcsp->getCtr(i)->connected())
                wcsp->getCtr(i)->propagate();
        for (int i = 0; i < wcsp->getElimBinOrder(); i++)
            if (wcsp->getElimBinCtr(i)->connected()
                && !wcsp->getElimBinCtr(i)->isSep())
                wcsp->getElimBinCtr(i)->propagate();
        for (int i = 0; i < wcsp->getElimTernOrder(); i++)
            if (wcsp->getElimTernCtr(i)->connected()
                && !wcsp->getElimTernCtr(i)->isSep())
                wcsp->getElimTernCtr(i)->propagate();
        wcsp->propagate();
        if (ToulBar2::verbose >= 0)
            cout << "   dual bound = " << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->getDDualBound() << std::setprecision(DECIMAL_POINT) << endl;
        //    printTightMatrix();
    }
}

Cost VACExtension::RASPSFindItThreshold()
{
    Cost result = ToulBar2::costThreshold;
    unsigned int size = ToulBar2::RASPSitThresholds.size();

    if (size > 0) {
        if (size < 3) {
            return ToulBar2::RASPSitThresholds[size - 1].first;
        }

        double lastRatio = ToulBar2::RASPSitThresholds[size - 1].second;
        for (unsigned int i = 0; i < size; i++) {
            ToulBar2::RASPSitThresholds[i].second = ToulBar2::RASPSitThresholds[i].second / lastRatio;
        }

        unsigned int i = 1;
        double stepSize = 2.0 / (double)size;
        while (i < size - 1 && atan2(ToulBar2::RASPSitThresholds[i + 1].second - ToulBar2::RASPSitThresholds[i - 1].second, stepSize) * 180.0 / PI < (double)abs(ToulBar2::RASPSangle)) {
            result = ToulBar2::RASPSitThresholds[i].first;
            i++;
        }
    }
    return result;
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

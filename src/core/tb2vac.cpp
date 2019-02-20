/** \file tb2vac.cpp
 *  \brief VAC implementation for binary cost functions.
 *
 *      \defgroup VAC Virtual Arc Consistency enforcing
 *  The three phases of VAC are enforced in three different "Pass".
 *  Bool(P) is never built. Instead specific functions (getVACCost) booleanize the WCSP on the fly.
 *  The domain variables of Bool(P) are the original variable domains (saved and restored using trailing at each iteration) 
 *  All the counter data-structures (k) are timestamped to avoid clearing them at each iteration.
 *  \note Simultaneously AC (and potentially DAC, EAC) are maintained by proper queuing.
 *  \see <em> Soft Arc Consistency Revisited. </em> Cooper et al. Artificial Intelligence. 2010.
 */

#include "tb2vacutils.hpp"
#include "search/tb2clusters.hpp"
#include <math.h> /* atan2 */
#define PI 3.14159265

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
{
    queueP = new stack<pair<int, int>>;
    queueR = new stack<pair<int, int>>;
    minlambda = MAX_COST;
    sumlb = MIN_COST;
    sumvars = 0;
    sumk = 0;
    theMaxK = 0;
    nlb = 0;
    // varAssign = -1;
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
    if (!done)
        c = UNIT_COST;
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

    if (Store::getDepth() == 0) {
        if (c < ToulBar2::costThresholdPre)
            c = MIN_COST;
    } else if (c < ToulBar2::costThreshold)
        c = MIN_COST;

    itThreshold = c;
}

void VACExtension::resetSupports()
{
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
    //    if (ToulBar2::verbose >= 8) cout << "previousThreshold: " << previousThreshold << " (" <<  bucket2 << "), threshold " << threshold  << " (" <<  bucket << ")" << endl;
    for (; bucket <= bucket2; bucket++) {
        VariableList* varlist = wcsp->getNCBucket(bucket);
        for (VariableList::iterator iter = varlist->begin(); iter != varlist->end();) {
            x = (VACVariable*)*iter;
            if (x->unassigned() && !(x->isNull(x->getMaxCost()))) {
                //                if (ToulBar2::RINS) {
                //                    for (EnumeratedVariable::iterator it = x->begin(); it != x->end(); ++it) {
                //                        if (x->getVACCost(*it) > MIN_COST) x->removeVAC(*it);
                //                    }
                //                }
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
    //cout << "CALL to VAC::propagate()" << endl;
    if (Store::getDepth() >= ToulBar2::vac) {
        inconsistentVariable = -1;
        return false;
    }

    assert(wcsp->verify());
    if (firstTime()) {
        init();
        if (ToulBar2::verbose >= 1)
            cout << "Dual bound before VAC: " << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->getDDualBound() << std::setprecision(DECIMAL_POINT) << endl;
    }

    bool isvac = true;
    bool util = false;

    breakCycles = 0;

    static vector<pair<VACVariable*, Value>> acSupport; /// \warning NOT SAFE FOR MULTITHREADING!!!
    bool acSupportOK = false;

    for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
        if (wcsp->getVar(i)->enumerated()) {
            EnumeratedVariable* xi = (EnumeratedVariable*)wcsp->getVar(i);
            xi->domSizeInBoolOfP = xi->getDomainSize();
        }
    }

    while (!util && itThreshold != MIN_COST) {
        int storedepth = Store::getDepth();
#ifdef AC2001
        resetSupports();
#endif
#ifdef INCREMENTALVAC
        prevItThreshold = -1;
        clear();
        Store::store();
#endif
        while (isvac && itThreshold != MIN_COST) {
#ifndef INCREMENTALVAC
            prevItThreshold = -1;
            clear();
            Store::store();
#endif
            minlambda = wcsp->getUb() - wcsp->getLb();
            inconsistentVariable = -1;
            nbIterations++;
            //assert(wcsp->verify()); // it modifies binary supports???
            assert(wcsp->NC.empty());
            assert(VAC.empty());
            while (!enqueueVAC(itThreshold, prevItThreshold) && itThreshold != MIN_COST) {
                prevItThreshold = itThreshold;
                nextScaleCost();
            }
            if (itThreshold == MIN_COST) {
                Store::restore(storedepth);
                inconsistentVariable = -1;
                return false;
            }
            if (ToulBar2::verbose >= 4)
                cout << "VAC itThreshold: " << itThreshold << " before enforcePass1 (prevItThreshold: " << prevItThreshold << ")" << endl;
            if (ToulBar2::verbose >= 8)
                cout << *wcsp << endl;
            enforcePass1();
            isvac = isVAC();
            //        cout << "and after enforcePass1" << endl << *wcsp << endl;
            assert(!isvac || checkPass1());
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
                // remember first arc consistent domain values in Bool(P) before restoring domains
                int nbassigned = 0;
                int nbassignedzero = 0;
                while (!SeekSupport.empty()) {
                    VACVariable* x = (VACVariable*)SeekSupport.pop();
                    if (x->assigned())
                        nbassigned++;
                    pair<VACVariable*, Value> p;
                    p.first = x;
                    p.second = x->getSup() + 1;
                    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
                        if (x->getCost(*iterX) == MIN_COST) {
                            if (x->assigned())
                                nbassignedzero++;
                            p.second = *iterX;
                            break;
                        }
                    }
                    if (x->canbe(p.second))
                        acSupport.push_back(p);
                }
                if (ToulBar2::debug && nbassignedzero > 0)
                    cout << "[" << Store::getDepth() << "] " << nbassignedzero << "/" << nbassigned - nbassignedzero << "/" << wcsp->numberOfUnassignedVariables() << " fixed/singletonnonzerocost/unassigned" << endl;
            }

            if (isvac) {

                ToulBar2::RINS_nbStrictACVariables = 0;

                for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
                    if (wcsp->getVar(i)->enumerated()) {
                        EnumeratedVariable* xi = (EnumeratedVariable*)wcsp->getVar(i);
                        xi->moreThanOne = false;
                        xi->strictACValue = xi->getSupport();
                        xi->domSizeInBoolOfP = 0;
                        int currentDomSize = xi->getDomainSize();
                        int initDomSize = xi->getDomainInitSize();
                        ValueCost domcost[currentDomSize];
                        wcsp->getEnumDomainAndCost(i, domcost);
                        xi->RINS_valuesToBeRemoved.clear();
                        int currentIndex = 0;

                        for (int v = 0; v < initDomSize; v++) {
                            if (currentIndex < currentDomSize && domcost[currentIndex].value == v) {
                                if (((VACVariable*)xi)->getVACCost(domcost[currentIndex].value) == MIN_COST) {
                                    xi->domSizeInBoolOfP += 1;
                                    xi->strictACValue = domcost[currentIndex].value;
                                }
                                currentIndex++;
                            } else {
                                xi->RINS_valuesToBeRemoved.push_back(v);
                            }
                        }
                        xi->moreThanOne = (xi->domSizeInBoolOfP > 1) ? true : false;
                        if (ToulBar2::strictAC == 2 && xi->domSizeInBoolOfP == 1) {
                            for (ConstraintList::iterator itc = xi->getConstrs()->begin();
                                 itc != xi->getConstrs()->end(); ++itc) {
                                Constraint* c = (*itc).constr;
                                if (c->isBinary()) {
                                    EnumeratedVariable* xj = (EnumeratedVariable*)(((BinaryConstraint*)c)->getVarDiffFrom(xi));
                                    if (xj->domSizeInBoolOfP > 1) // TODO: xj may be not refreshed yet
                                        xi->moreThanOne = true;
                                }
                            }
                        }
                        if (ToulBar2::strictAC == 3 && xi->domSizeInBoolOfP > 1 && xi->getCost(xi->strictACValue) == MIN_COST) {
                            xi->moreThanOne = false;
                            for (ConstraintList::iterator itc = xi->getConstrs()->begin();
                                 !xi->moreThanOne && itc != xi->getConstrs()->end(); ++itc) {
                                Constraint* c = (*itc).constr;
                                if (c->isBinary()) {
                                    BinaryConstraint* cij = ((BinaryConstraint*)c);
                                    EnumeratedVariable* xj = (EnumeratedVariable*)cij->getVarDiffFrom(xi);
                                    if (xj->wcspIndex < xi->wcspIndex && !xj->moreThanOne && cij->getCost(xi, xj, xi->strictACValue, xj->strictACValue) > MIN_COST) {
                                        xi->moreThanOne = true;
                                    }
                                }
                            }
                        }

                        if (xi->domSizeInBoolOfP == 1) {
                            ToulBar2::RINS_nbStrictACVariables++;
                        }
                    }
                }

                ToulBar2::nbTimesIsVAC++;
                if (itThreshold > 1)
                    ToulBar2::nbTimesIsVACitThresholdMoreThanOne++;

                if (ToulBar2::RINS_saveitThresholds) {
                    double ratio = (ToulBar2::RINS_nbStrictACVariables == 0) ? 0.0000000001 : (((double)ToulBar2::RINS_nbStrictACVariables / (double)ToulBar2::nbvar) / (double)itThreshold);
                    //cout << "itThreshold: " << itThreshold << " strictAC: " << ToulBar2::RINS_nbStrictACVariables << " ratio: " << ratio << endl;
                    ToulBar2::RINS_itThresholds.push_back(std::make_pair(itThreshold, ratio));
                }

                //cout << "Nb Variables With BoolDomSize Zero: " << nbDomSizeZero << " One: " << ToulBar2::RINS_nbStrictACVariables << " More Than One: " << nbDomSizeMore << endl;
                //cout << ((double)ToulBar2::RINS_nbStrictACVariables / (double)ToulBar2::nbvar) / (double) itThreshold << endl;
                //cout << "Nb Variables That Changed State: " << nbVariablesChanged << endl;

                if (ToulBar2::RINS) {
                    ToulBar2::RINS = false;
                    if (ToulBar2::RINS_HBFSnodes > 0) {
                        Store::restore(storedepth);
                        Store::store();
                        /*string fileName = "beforeAssignment_";
                        fileName += std::to_string(ToulBar2::RINS_HBFSnodes);
                        fileName += ".wcsp";
                        ofstream pb(fileName.c_str());
                        wcsp->dump(pb, true);*/
                    }
                    Solver* solver = (Solver*)wcsp->getSolver();
                    //cout << "[" << Store::getDepth() << "," << wcsp->getNbNodes() << "]" << " VAC Propagate RINS = true" << endl;
                    //cout << "itThreshold: " << itThreshold << " nbStrictACVariables / nbVariables: " << (double)ToulBar2::RINS_nbStrictACVariables / (double)ToulBar2::nbvar << endl;
                    int storehbfs = ToulBar2::hbfs;
                    int storehbfsGlobalLimit = ToulBar2::hbfsGlobalLimit;
                    int storehbfsLimit = solver->hbfsLimit;
                    int storeVac = ToulBar2::vac;
                    int storenbBacktracksLimit = solver->nbBacktracksLimit;
                    int storerestart = ToulBar2::restart;
                    bool storeLimited = ToulBar2::limited;
                    int storenbBacktracks = solver->nbBacktracks;
                    int storeStrictAC = ToulBar2::strictAC;
                    assert(storerestart < 1);

                    ToulBar2::strictAC = 0;
                    ToulBar2::vac = 0;
                    ToulBar2::hbfs = 0;
                    ToulBar2::hbfsGlobalLimit = 0;
                    ToulBar2::limited = false;
                    if (ToulBar2::useRINS <= 1) {
                        ToulBar2::restart = 1; // random variable selection for breaking ties in DFS
                    } else {
                        ToulBar2::restart = -1; // no randomness for LDS
                    }
                    solver->hbfsLimit = LONGLONG_MAX;
                    solver->nbBacktracksLimit = solver->nbBacktracks + 1000;

                    Cost lastUB = wcsp->getUb();
                    try {
                        try {
                            // Current WCSP is AC(Bool(P))
                            vector<int> variables;
                            vector<Value> values;

                            for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
                                if (wcsp->getVar(i)->enumerated()) {
                                    EnumeratedVariable* xi = (EnumeratedVariable*)wcsp->getVar(i);
                                    //if(xi->strictACValue != xi->RINS_lastValue)
                                    //cout << i << " strictAC val: " << xi->strictACValue << " incumbentValue: " << xi->RINS_lastValue << endl;
                                    int nbValues = xi->RINS_valuesToBeRemoved.size();
                                    for (int j = 0; j < nbValues; j++) {
                                        if (ToulBar2::RINS_HBFSnodes == 0) {
                                            xi->remove(xi->RINS_valuesToBeRemoved[j]);
                                        }
                                        if (ToulBar2::RINS_HBFSnodes > 0 && ToulBar2::useRINS == -1 && xi->RINS_valuesToBeRemoved[j] != xi->RINS_lastValue) {
                                            xi->remove(xi->RINS_valuesToBeRemoved[j]);
                                        }
                                    }
                                    if (xi->domSizeInBoolOfP == 1) {
                                        if (ToulBar2::RINS_HBFSnodes == 0 || xi->RINS_lastValue == -1) {
                                            //cout << "add " << i << " to variables" << endl;
                                            variables.push_back(i);
                                            values.push_back(xi->strictACValue);
                                        } else if (xi->strictACValue == xi->RINS_lastValue) {
                                            //cout << *iter << " ";
                                            variables.push_back(i);
                                            values.push_back(xi->strictACValue);
                                        }
                                        //else update domain
                                    }
                                    if (xi->cannotbe(xi->getSupport()))
                                        xi->findSupport(); // update support values
                                    xi->propagateNC(); // and update maxcost values before propagate done in assignLS
                                }
                            }
                            if (variables.size() > 0)
                                wcsp->assignLS(variables, values, true); // option true: make sure already assigned variables are removed from Solver::unassignedVars
                            /*string fileName = "afterAssignment_";
                            fileName += std::to_string(ToulBar2::RINS_HBFSnodes);
                            fileName += ".wcsp";
                        ofstream pb(fileName.c_str());
                        wcsp->dump(pb, true);*/
                            if (ToulBar2::useRINS <= 1) {
                                cout << "call to recursiveSolve from VAC" << endl;
                                solver->recursiveSolve(wcsp->getLb()); // look at its search tree (if a new solution is found, UB should be updated automatically)
                            } else {
                                cout << "call to recursiveSolveLDS from VAC" << endl;
                                solver->recursiveSolveLDS(ToulBar2::useRINS - 1);
                            }
                        } catch (Contradiction) {
                            wcsp->whenContradiction();
                        }
                    } catch (NbBacktracksOut) {
                    }

                    ToulBar2::limited = storeLimited; // still a complete search
                    ToulBar2::hbfs = storehbfs;
                    ToulBar2::hbfsGlobalLimit = storehbfsGlobalLimit;
                    ToulBar2::restart = storerestart;
                    solver->hbfsLimit = storehbfsLimit;
                    solver->nbBacktracksLimit = storenbBacktracksLimit + solver->nbBacktracks - storenbBacktracks;
                    ToulBar2::vac = storeVac;
                    ToulBar2::strictAC = storeStrictAC;

                    //cout << "[" << Store::getDepth() << "," << wcsp->getNbNodes() << "]" << " VAC Propagate RINS = false" << endl;
                    //ToulBar2::RINS = false;
                    Store::restore(storedepth);
                    inconsistentVariable = -1;
                    return (wcsp->getUb() < lastUB);
                }

                /*string fileName = "problem_";
			fileName += to_string(wcsp->getNbNodes());
			fileName += "_";
			fileName += to_string(nbIterations);
			fileName += ".wcsp";
			ofstream pb(fileName.c_str());
			wcsp->dump_strictAC(pb, true);*/
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
            }
#ifndef INCREMENTALVAC
            Store::restore(storedepth);

#endif
        }
#ifdef INCREMENTALVAC
        Store::restore(storedepth);
#endif
        if (!isvac) {
            enforcePass2();
            if (ToulBar2::verbose > 0)
                cout << "VAC dual bound: " << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->getDDualBound() << std::setprecision(DECIMAL_POINT) << "    incvar: " << inconsistentVariable << "    minlambda: " << minlambda << "      itThreshold: " << itThreshold << endl;
            if (CUT(wcsp->getLb() + minlambda, wcsp->getUb())) {
                if (ToulBar2::weightedDegree)
                    wcsp->conflict();
                throw Contradiction();
            }
            util = enforcePass3();
        } else {
            nextScaleCost();
        }
    }

    if (ToulBar2::vacValueHeuristic && acSupportOK) {
        // update current unary support if possible && needed
        for (vector<pair<VACVariable*, Value>>::iterator iter = acSupport.begin(); iter != acSupport.end(); ++iter) {
            VACVariable* x = iter->first;
            if (x->canbe(iter->second)) {
                assert(x->getCost(iter->second) == MIN_COST);
                if (ToulBar2::verbose > 0
                    && x->getSupport() != iter->second)
                    cout << "CHANGE SUPPORT " << x->getName() << " from " << x->getSupport() << " to " << iter->second << endl;
                x->setSupport(iter->second);
            } else {
                if (ToulBar2::verbose > 0)
                    cout << "WARNING: BAD AC SUPPORT " << iter->second << " FOR VARIABLE " << *x << endl;
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
    //    if (ToulBar2::verbose >= 8) cout << "revise " << *((EnumeratedVariable *) xi) << endl;
    for (EnumeratedVariable::iterator it = xi->begin(); it != xi->end(); ++it) {
        Value v = *it;
        //        if (ToulBar2::verbose >= 8) cout << "check variable " << xi->getName() << " with value " << v << " and cost " << xi->getVACCost(v) << endl;
        if (xi->getVACCost(v) > MIN_COST) {
            xi->removeVAC(v);
        } else if (cij->revise(xi, v)) {
            wipeout = xi->removeVAC(v);
            xi->setKiller(v, xj->wcspIndex);
            queueP->push(pair<int, int>(xi->wcspIndex, v));
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
    VACVariable* xj;
    VACBinaryConstraint* cij;

    while (!VAC.empty()) {
        xj = (VACVariable*)VAC.pop_first();
        for (EnumeratedVariable::iterator it = xj->begin(); it != xj->end(); ++it) {
            if (xj->getVACCost(*it) > MIN_COST)
                xj->removeVAC(*it);
        }
        for (ConstraintList::iterator itc = xj->getConstrs()->begin();
             itc != xj->getConstrs()->end(); ++itc) {
            Constraint* c = (*itc).constr;
            if (c->isBinary()) {
                cij = (VACBinaryConstraint*)c;
                if (enforcePass1(xj, cij))
                    return;
            }
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
            if (c->isBinary()) {
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

    //if (ToulBar2::verbose > 0)  cout << "VAC Enforce Pass 2" << endl;

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
        } else
            xi0->setMark(v, nbIterations);
    }

    while (!queueP->empty()) {
        i = queueP->top().first;
        v = queueP->top().second;
        queueP->pop();
        xi = (VACVariable*)wcsp->getVar(i);
        if (xi->isMarked(v, nbIterations)) {
            j = xi->getKiller(v);
            xj = (VACVariable*)wcsp->getVar(j);
            queueR->push(pair<int, int>(i, v));
            cij = (VACBinaryConstraint*)xi->getConstr(xj);
            assert(cij);
            //if (ToulBar2::verbose > 6) cout << "x" << xi->wcspIndex << "," << v << "   killer: " << xj->wcspIndex << endl;

            for (EnumeratedVariable::iterator itj = xj->begin();
                 itj != xj->end(); ++itj) {
                Value w = *itj;
                Cost costij = cij->getVACCost(xi, xj, v, w);
                if (costij > MIN_COST) {
                    int tmpK = xi->getK(v, nbIterations);
                    if (xj->getKiller(w) == i
                        && xj->isMarked(w, nbIterations))
                        tmpK += xj->getK(w, nbIterations);
                    if (!CUT(wcsp->getLb() + costij,
                            wcsp->getUb())) {
                        if ((costij / tmpK) < minlambda) {
                            minlambda = costij / tmpK;
                            if (minlambda < UNIT_COST) { //A cost function bottleneck here !
                                bneckCost = costij;
                                bneckCF = cij;
                                bneckVar = -1;
                            }
                        }
                    } else {
                        if ((costij / tmpK) < minlambda) { // costij should be made infinite to avoid to decrease minlambda
                            Cost cost = tmpK * minlambda - costij;
                            assert(cost > MIN_COST);
                            assert(ToulBar2::verbose < 1 || ((cout << "inflate(C" << cij->getVar(0)->getName() << "," << cij->getVar(1)->getName() << ", (" << ((xi == cij->getVar(0)) ? v : w) << "," << ((xi == cij->getVar(0)) ? w : v) << "), " << cost << ")" << endl), true));
                            cij->addcost(xi, xj, v, w, cost);
                        }
                    }
                } else {
                    int tmpK = xi->getK(v,
                                   nbIterations)
                        - cij->getK(xj, w, nbIterations);
                    if (tmpK > 0) {
                        xj->addToK(w, tmpK,
                            nbIterations);
                        cij->setK(xj, w,
                            xi->getK(v,
                                nbIterations),
                            nbIterations);
                        Cost cost = xj->getVACCost(w);
                        if (cost == MIN_COST)
                            xj->setMark(w,
                                nbIterations);
                        else { // we assume NC has been done before
                            assert(!CUT(wcsp->getLb() + cost, wcsp->getUb()));
                            tmplambda = cost / xj->getK(w, nbIterations);
                            if (tmplambda < minlambda) {
                                minlambda = tmplambda;
                                if (minlambda < UNIT_COST) { // A unary cost bottleneck here
                                    bneckVar
                                        = j;
                                    bneckCF
                                        = NULL;
                                    bneckCost
                                        = cost;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    //if (maxK == 0) {
    //  maxK = wcsp->getUb() - wcsp->getLb();
    //}
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

    //if (ToulBar2::verbose > 2) cout << "VAC Enforce Pass 3.   minlambda " << minlambda << " , var: " << inconsistentVariable << endl;
    int i, j;
    VACVariable *xi, *xj;
    VACBinaryConstraint* cij;
    int i0 = inconsistentVariable;
    VACVariable* xi0 = (VACVariable*)wcsp->getVar(i0);
    Value w;

    int maxk = 0;

    if (!util) { // Empty R ?
        assert(bneckVar != -1 || bneckCF != NULL);
        if (bneckVar != -1) {
            xi = (VACVariable*)wcsp->getVar(bneckVar);
            xi->setThreshold(bneckCost + UNIT_COST);
        } else {
            bneckCF->setThreshold(bneckCost + UNIT_COST);
        }
        breakCycles++;
        //if (ToulBar2::verbose > 1) cout << "BreakCycle   (var: " << minvar << ", val= " << minval << ")   thr: " <<  xi->getThreshold() << endl;
        if (breakCycles > 5) {
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

    while (!queueR->empty()) {
        j = queueR->top().first;
        w = queueR->top().second;
        queueR->pop();
        xj = (VACVariable*)wcsp->getVar(j);
        i = xj->getKiller(w);
        xi = (VACVariable*)wcsp->getVar(i);
        cij = (VACBinaryConstraint*)xi->getConstr(xj);
        assert(cij);

        int xjk = xj->getK(w, nbIterations);
        if (maxk < xjk)
            maxk = xjk;

        for (EnumeratedVariable::iterator iti = xi->begin();
             iti != xi->end(); ++iti) {
            Value v = *iti;
            if (cij->getK(xi, v, nbIterations) != 0) {
                Cost ecost = lambda * cij->getK(xi, v, nbIterations);
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
        cij->VACproject(xj, w, lambda * xj->getK(w, nbIterations));
    }
    sumk += maxk;
    if (maxk > theMaxK)
        theMaxK = maxk;

    xi0->extendAll(lambda);
    xi0->projectLB(lambda);
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
            singletonI.insert(MAX_DOMAIN_SIZE * i + a);
    }
}

void VACExtension::updateSingleton()
{
    set<int>& s1 = singleton;
    set<int> s2(singletonI);
    singletonI.clear();
    set_intersection(s1.begin(), s1.end(),
        s2.begin(), s2.end(),
        inserter(singletonI, singletonI.begin()));
    singleton.clear();
}

void VACExtension::removeSingleton()
{
    set<int>& s = singletonI;
    set<int>::iterator it = s.begin();
    while (it != s.end()) {
        int ivar = *it / MAX_DOMAIN_SIZE;
        Value a = *it % MAX_DOMAIN_SIZE;
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
        long double mean = to_double(sumlb) / (long double)nlb;
        cout << "VAC mean lb/incr: " << std::setprecision(DECIMAL_POINT) << mean << "     total increments: " << nlb
             << "     cyclesize: " << (double)sumvars / (double)nlb << "     k: " << (double)sumk / (double)nlb << " (mean), " << theMaxK << " (max)" << endl;
    }
    if (ini && nlb > 0 && sumlb > MIN_COST && ToulBar2::verbose >= 0) {
        if (ToulBar2::uai)
            cout << "VAC dual bound: " << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->getDDualBound() << std::setprecision(DECIMAL_POINT) << " energy: " << -(wcsp->Cost2LogProb(wcsp->getLb()) + ToulBar2::markov_log) << " (iter:" << nlb << ")" << endl;
        else
            cout << "VAC dual bound: " << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->getDDualBound() << std::setprecision(DECIMAL_POINT) << " (iter:" << nlb << ")" << endl;
    }

    cout << "Number of VAC iterations: " << nbIterations << endl;
    cout << "Number of times is VAC: " << ToulBar2::nbTimesIsVAC << " Number of times isvac and itThreshold > 1: " << ToulBar2::nbTimesIsVACitThresholdMoreThanOne << endl;
    //sort(heap.begin(), heap.end(), cmp_function);
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
            //cout << "it " << ntimes << "   changed: " << nchanged << endl;
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

void VACExtension::RINS_finditThreshold()
{

    cout << "call to find itThreshold with size " << ToulBar2::RINS_itThresholds.size() << endl;

    if (ToulBar2::RINS_itThresholds.size() == 0) {
        ToulBar2::RINS_lastitThreshold = ToulBar2::costThresholdPre; // ToulBar2::costThreshold
        return;
    }

    if (ToulBar2::RINS_itThresholds.size() < 3) {
        ToulBar2::RINS_lastitThreshold = ToulBar2::RINS_itThresholds[ToulBar2::RINS_itThresholds.size() - 1].first;
        return;
    }

    double lastRatio = ToulBar2::RINS_itThresholds[ToulBar2::RINS_itThresholds.size() - 1].second;
    //cout << std::fixed << std::setprecision(7);
    //cout << lastRatio << endl;

    for (unsigned int i = 0; i < ToulBar2::RINS_itThresholds.size(); i++) {
        //cout << ToulBar2::RINS_itThresholds[i].second << " ";
        ToulBar2::RINS_itThresholds[i].second = ToulBar2::RINS_itThresholds[i].second / lastRatio;
        //cout << ToulBar2::RINS_itThresholds[i].second << endl;
    }

    unsigned int i = 1;
    double stepSize = 2.0 / (double)ToulBar2::RINS_itThresholds.size();

    //cout << "getUb: " << wcsp->getUb() << endl;

    while (i < ToulBar2::RINS_itThresholds.size() && atan2(ToulBar2::RINS_itThresholds[i + 1].second - ToulBar2::RINS_itThresholds[i - 1].second, stepSize) * 180.0 / PI < (double)ToulBar2::RINS_angle) {
        //cout << ToulBar2::RINS_itThresholds[i].second << endl;
        ToulBar2::RINS_lastitThreshold = ToulBar2::RINS_itThresholds[i].first;
        i++;
    }

    //ToulBar2::costThreshold = ToulBar2::RINS_lastitThreshold;

    //cout << std::fixed << std::setprecision(DECIMAL_POINT);
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

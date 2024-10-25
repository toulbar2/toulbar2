/*
 * **************** Data-structure to manage Cluster Tree Decomposition *******************
 *
 */

#include "tb2clusters.hpp"
#include "core/tb2naryconstr.hpp"
// #include "applis/tb2pedigree.hpp"
// #include "applis/tb2haplotype.hpp"

/*
 * Comparison between cluster sons
 *
 */

int Cluster::clusterCounter = 0;

bool CmpClusterStructBasic::operator()(const Cluster* lhs, const Cluster* rhs) const
{
    return lhs && rhs && (lhs->getIndex() < rhs->getIndex());
}
bool CmpClusterStruct::operator()(const Cluster* lhs, const Cluster* rhs) const
{
    if (ToulBar2::bilevel && lhs->getParent() == lhs->getTreeDec()->getRoot()) {
        // do not sort clusters by size if bilevel at depth 1 (keep id ordering)
        return lhs && rhs && (lhs->getIndex() < rhs->getIndex());
    } else {
        return lhs && rhs && (lhs->sepSize() < rhs->sepSize() || (lhs->sepSize() == rhs->sepSize() && (lhs->getNbVarsTree() < rhs->getNbVarsTree() || (lhs->getNbVarsTree() == rhs->getNbVarsTree() && lhs->getIndex() < rhs->getIndex()))));
    }
}

/*
 * Comparison between variables (used for tie breaking)
 *
 */

WCSP* CmpVarStruct::wcsp = NULL;

bool CmpVarStruct::operator()(const int lhs, const int rhs) const
{
    return CmpVarStruct::wcsp->getDACOrder(lhs) < CmpVarStruct::wcsp->getDACOrder(rhs);
}

/*
 * Separator class derived from NaryConstraint
 *
 */

Separator::Separator(WCSP* wcsp, EnumeratedVariable** scope_in, int arity_in)
    : AbstractNaryConstraint(wcsp, scope_in, arity_in)
    , cluster(NULL)
    , isUsed(false)
    , lbPrevious(MIN_COST)
    , optPrevious(false)
{
    for (int i = 0; i < arity_in; i++) {
        unsigned int domsize = scope_in[i]->getDomainInitSize();
        vars.insert(scope_in[i]->wcspIndex);
        if (domsize > (unsigned int)std::numeric_limits<tValue>::max()) {
            cerr << "Nary constraints overflow. Extend tValue type range." << endl;
            throw BadConfiguration();
        }
    }
    t = Tuple(arity_in, 0);

    linkSep.content = this;

    // initial "delayed" propagation
    if (arity_ == 0) {
        queueSep();
    } else {
        for (int i = 0; i < arity_; i++) {
            if (getVar(i)->assigned())
                assign(i);
        }
    }
}

Separator::Separator(WCSP* wcsp)
    : AbstractNaryConstraint(wcsp)
    , isUsed(false)
    , lbPrevious(MIN_COST)
    , optPrevious(false)
{
}

void Separator::setup(Cluster* cluster_in)
{
    cluster = cluster_in;
    AbstractNaryConstraint::cluster = cluster_in->getParent()->getId();
    delta.clear();
    TVars::iterator it = vars.begin();
    while (it != vars.end()) {
        EnumeratedVariable* var = (EnumeratedVariable*)cluster->getWCSP()->getVar(*it);
        delta.push_back(vector<StoreCost>(var->getDomainInitSize(), StoreCost(MIN_COST)));
        ++it;
    }

    // take into account the fact that the cluster may be used as its descendants
    int nvars = cluster->getNbVarsTree();
    if (!nvars)
        return;

    s = Tuple(cluster->getNbVarsTree(), 0);
}

void Separator::assign(int varIndex)
{
    if (connected(varIndex)) {
        deconnect(varIndex);
        assert(getNonAssigned() >= 0);
        if (getNonAssigned() == 0) {
            if (ToulBar2::bilevel && (!cluster || cluster->getParent() == wcsp->getTreeDec()->getRoot()))
                return; // TODO: how to reuse Problem2 nogood if it exists? (but should never collect NegProblem2 separator)
            assert(!cluster || cluster->isActive());
            queueSep();
        }
    }
}

void Separator::propagate()
{
    if (ToulBar2::verbose >= 3)
        cout << this << " propagate C" << cluster->getId() << " " << getNonAssigned() << " " << cluster->getParent()->getId() << " " << connected() << endl;
    for (int i = 0; connected() && i < arity_; i++) {
        if (getVar(i)->assigned())
            assign(i);
    }
    if (cluster->getIsCurrInTD() && getNonAssigned() == 0 && wcsp->getTreeDec()->isInCurrentClusterSubTree(cluster->getParent()->getId())) {
        wcsp->revise(this);
        if (ToulBar2::allSolutions) {
            Cost res = MIN_COST;
            BigInteger nb = 0.;
            getSg(res, nb);
            if (nb == 0.) {
                if (ToulBar2::verbose >= 1)
                    cout << "use #good " << this << endl;
                THROWCONTRADICTION;
            } else if (nb > 0.)
                unqueueSep();
        } else {
            Cost clb = MIN_COST;
            Cost cub = MAX_COST;
            get(clb, cub, &cluster->open);
            bool opt = (clb == cub);
            if (cluster->isActive()) {
                Cost lbpropa = cluster->getLbRec();
                Cost lb = clb - lbpropa;
                if (opt || lb > MIN_COST) {
                    if (ToulBar2::verbose >= 1)
                        cout << "nogood C" << cluster->getId() << " used in advance (lbpropa=" << lbpropa << " ,lb+=" << lb << ")" << endl;
                    assert(lb >= MIN_COST);
                    if (opt)
                        unqueueSep();
                    // atomic operations:
                    isUsed = true;
                    cluster->deactivate();
                    assert(cluster->getParent()->getId() == Constraint::cluster);
                    cluster->getParent()->increaseLb(lbpropa);
                    if (lb > MIN_COST)
                        projectLB(lb); // project into global lb and into parent cluster
                    lbPrevious = clb;
                    optPrevious = opt;
                    // end of atomic operations.
                }
            } else if (isUsed && cluster->getParent()->isActive()) {
                if (clb > lbPrevious || (opt == true && optPrevious == false)) {
                    if (ToulBar2::verbose >= 1)
                        cout << "nogood C" << cluster->getId() << " used in advance (lbPrevious=" << lbPrevious << " ,optPrevious=" << optPrevious << " ,clb=" << clb << " ,opt=" << opt << ")" << endl;
                    if (opt)
                        unqueueSep();
                    // atomic operations:
                    if (clb > lbPrevious)
                        projectLB(clb - lbPrevious); // project into global lb and into parent cluster
                    lbPrevious = clb;
                    optPrevious = opt;
                    // end of atomic operations.
                }
            }
        }
    }
}

void Separator::set(Cost clb, Cost cub, Solver::OpenList** open)
{
    assert(ToulBar2::bilevel || clb <= cub);
    int i = 0;
    WCSP* wcsp = cluster->getWCSP();
    Cost deltares = MIN_COST;
    if (ToulBar2::verbose >= 1)
        cout << "( ";
    TVars::iterator it = vars.begin();
    while (it != vars.end()) {
        assert(wcsp->assigned(*it));
        tValue val = wcsp->toIndex(*it, wcsp->getValue(*it));
        if (ToulBar2::verbose >= 1)
            cout << "(" << *it << "," << val << ") ";
        t[i] = val;
        deltares += delta[i][val];
        ++it;
        i++;
    }
    if (ToulBar2::verbose >= 1)
        cout << ")";
    assert(clb < cub || clb + deltares >= MIN_COST);
    // assert(nogoods.find(t) == nogoods.end() || nogoods[Tuple(t)].second <= MAX(MIN_COST, c + deltares));
    TNoGoods::iterator itng = nogoods.find(t);
    if (ToulBar2::verbose >= 3) {
        cout << " <C" << cluster->getId() << ",";
        cout << t;
        cout << "," << MAX(MIN_COST, clb + deltares) << "," << MAX(MIN_COST, cub + deltares) << ">" << endl;
    }
    if (open && ToulBar2::hbfs) {
        if (*open) {
            // open node list already found => the corresponding nogood has been created before
            assert(itng != nogoods.end());
            assert(*open == &itng->second.third);
            itng->second.first = MAX(itng->second.first, clb + deltares);
            itng->second.second = MIN(itng->second.second, MAX(MIN_COST, cub + ((cub < MAX_COST) ? deltares : MIN_COST)));
            if (ToulBar2::verbose >= 1)
                cout << " Learn nogood " << itng->second.first << ", cub= " << itng->second.second << ", delta= " << deltares << " on cluster " << cluster->getId() << endl;
        } else {
            assert(itng == nogoods.end());
            nogoods[t] = make_triplet(MAX(MIN_COST, clb + deltares), MAX(MIN_COST, cub + ((cub < MAX_COST) ? deltares : MIN_COST)), Solver::OpenList());
            if (ToulBar2::verbose >= 1)
                cout << " Learn nogood " << nogoods[t].first << ", cub= " << nogoods[t].second << ", delta= " << deltares << " on cluster " << cluster->getId() << endl;
            *open = &nogoods[t].third;
        }
    } else {
        if (itng == nogoods.end()) {
            nogoods[t] = make_triplet(MAX(MIN_COST, clb + deltares), MAX(MIN_COST, cub + ((cub < MAX_COST) ? deltares : MIN_COST)), Solver::OpenList());
            if (ToulBar2::verbose >= 1)
                cout << " Learn nogood " << nogoods[t].first << ", cub= " << nogoods[t].second << ", delta= " << deltares << " on cluster " << cluster->getId() << endl;
        } else {
            if (ToulBar2::bilevel && ToulBar2::verbose >= 0)
                cout << "Warning! nogood already solved on cluster " << cluster->getId() << " !?!" << endl;
            itng->second.first = MAX(itng->second.first, clb + deltares);
            itng->second.second = MIN(itng->second.second, MAX(MIN_COST, cub + ((cub < MAX_COST) ? deltares : MIN_COST)));
            if (ToulBar2::verbose >= 1)
                cout << " Learn nogood " << itng->second.first << ", cub= " << itng->second.second << ", delta= " << deltares << " on cluster " << cluster->getId() << endl;
        }
    }
}

void Separator::setSg(Cost c, BigInteger nb)
{
    int i = 0;
    WCSP* wcsp = cluster->getWCSP();
    Cost deltares = MIN_COST;
    if (ToulBar2::verbose >= 1)
        cout << "( ";
    TVars::iterator it = vars.begin();
    while (it != vars.end()) {
        assert(wcsp->assigned(*it));
        tValue val = wcsp->toIndex(*it, wcsp->getValue(*it));
        if (ToulBar2::verbose >= 1)
            cout << "(" << *it << "," << val << ") ";
        t[i] = val;
        deltares += delta[i][val];
        ++it;
        i++;
    }
    assert(c + deltares >= MIN_COST);
    if (ToulBar2::verbose >= 1)
        cout << ") Learn #good with " << nb << " solutions" << endl; // /" << cluster->getVarsTree().size() << endl;
    sgoods[t] = TPairSG(MAX(MIN_COST, c + deltares), nb);
}

void Separator::setF(bool free)
{
    int i = 0;
    WCSP* wcsp = cluster->getWCSP();
    if (ToulBar2::verbose >= 1)
        cout << "( ";
    TVars::iterator it = vars.begin();
    while (it != vars.end()) {
        assert(wcsp->assigned(*it));
        tValue val = wcsp->toIndex(*it, wcsp->getValue(*it));
        if (ToulBar2::verbose >= 1)
            cout << "(" << *it << "," << val << ") ";
        t[i] = val;
        ++it;
        i++;
    }
    if (ToulBar2::verbose >= 1)
        cout << ") Learn from heuristic of freedom with " << free << endl;
    frees[t] = free;
}

bool Separator::setFInc()
{
    int i = 0;
    WCSP* wcsp = cluster->getWCSP();
    TVars::iterator it = vars.begin();
    while (it != vars.end()) {
        tValue val = wcsp->toIndex(*it, wcsp->getValue(*it));
        t[i] = val;
        ++it;
        i++;
    }

    int nb;
    TFreesLimit::iterator itsg = freesLimit.find(t);
    if (itsg != freesLimit.end()) {
        nb = itsg->second;
        if (nb < ToulBar2::heuristicFreedomLimit) {
            assert(frees.find(t) != frees.end() && frees[t] == true);
            return true;
        } else {
            if (ToulBar2::verbose >= 1) {
                cout << " limit of " << nb << " reached for cluster " << cluster->getId() << " with separator assignment " << t << endl;
            }
            frees[t] = false;
            return false;
        }
    } else {
        freesLimit[t] = 0;
        frees[t] = true;
        return true;
    }
}

void Separator::freeIncS()
{
    int i = 0;
    WCSP* wcsp = cluster->getWCSP();
    TVars::iterator it = vars.begin();
    while (it != vars.end()) {
        tValue val = wcsp->toIndex(*it, wcsp->getValue(*it));
        t[i] = val;
        ++it;
        i++;
    }
    freesLimit[t] += 1;
    if (ToulBar2::verbose >= 1) {
        cout << " hybridSolve ends without any improvement for cluster " << cluster->getId() << " (separator limit: " << freesLimit[t] << " for assignment " << t << endl;
    }
}

Cost Separator::getCurrentDeltaUb()
{
    int i = 0;
    WCSP* wcsp = cluster->getWCSP();
    Cost sumdelta = MIN_COST;
    TVars::iterator it = vars.begin();
    while (it != vars.end()) {
        if (wcsp->assigned(*it)) {
            tValue val = wcsp->toIndex(*it, wcsp->getValue(*it));
            sumdelta += delta[i][val];
        } else {
            EnumeratedVariable* x = (EnumeratedVariable*)wcsp->getVar(*it);
            if (wcsp->td->isDeltaModified(x->wcspIndex)) {
                Cost del = -MAX_COST;
                for (EnumeratedVariable::iterator itx = x->begin(); itx != x->end(); ++itx) {
                    tValue val = x->toIndex(*itx);
                    // Cost unaryc = x->getCost(val);
                    // Could use delta[i][val]-unaryc for pure RDS with only one separator per variable
                    if (del < delta[i][val])
                        del = delta[i][val];
                }
                assert(del > -MAX_COST);
                sumdelta += del;
            }
        }
        ++it;
        i++;
    }
    return sumdelta;
}

Cost Separator::getCurrentDeltaLb()
{
    int i = 0;
    WCSP* wcsp = cluster->getWCSP();
    Cost sumdelta = MIN_COST;
    TVars::iterator it = vars.begin();
    while (it != vars.end()) {
        if (wcsp->assigned(*it)) {
            tValue val = wcsp->toIndex(*it, wcsp->getValue(*it));
            sumdelta += delta[i][val];
        } else {
            EnumeratedVariable* x = (EnumeratedVariable*)wcsp->getVar(*it);
            if (wcsp->td->isDeltaModified(x->wcspIndex)) {
                Cost del = MAX_COST;
                for (EnumeratedVariable::iterator itx = x->begin(); itx != x->end(); ++itx) {
                    tValue val = x->toIndex(*itx);
                    // Cost unaryc = x->getCost(val);
                    // Could use delta[i][val]-unaryc for pure RDS with only one separator per variable
                    if (del > delta[i][val])
                        del = delta[i][val];
                }
                assert(del < MAX_COST);
                sumdelta += del;
            }
        }
        ++it;
        i++;
    }
    return sumdelta;
}

bool Separator::get(Cost& clb, Cost& cub, Solver::OpenList** open)
{
    int i = 0;
    clb = MIN_COST;
    cub = MIN_COST;

    if (ToulBar2::verbose >= 1)
        cout << "( ";
    TVars::iterator it = vars.begin();
    while (it != vars.end()) {
        assert(wcsp->assigned(*it));
        tValue val = wcsp->toIndex(*it, wcsp->getValue(*it));
        if (ToulBar2::verbose >= 1)
            cout << "(" << *it << "," << val << ") ";
        t[i] = val; // build the tuple
        clb -= delta[i][val]; // delta structure
        cub -= delta[i][val]; // delta structure
        ++it;
        i++;
    }
    TNoGoods::iterator itng = nogoods.find(t);
    if (itng != nogoods.end()) {
        TPairNG& p = itng->second; // it is crucial here to get a reference to the data triplet object instead of a copy, otherwise open node list would be copied
        if (ToulBar2::verbose >= 1)
            cout << ") Use nogood " << p.first << ", delta=" << clb << " (cub=" << p.second << ") on cluster " << cluster->getId() << " (active=" << cluster->isActive() << ")" << endl;
        assert(p.first < p.second || clb + p.first >= MIN_COST);
        clb += p.first;
        cub += p.second;
        cub = MAX(MIN_COST, cub);
        cluster->setUb(cub);
        if (open)
            *open = &p.third;
        if (ToulBar2::btdMode >= 2) {
            Cost lbrds = cluster->getLbRDS();
            assert(clb < cub || clb >= lbrds);
            clb = MAX(lbrds, clb);
        } else {
            clb = MAX(MIN_COST, clb);
        }
        return true;
    } else {
        clb = (ToulBar2::btdMode >= 2) ? cluster->getLbRDS() : MIN_COST;
        cub = MAX_COST;
        cluster->setUb(MAX_COST);
        if (open)
            *open = NULL;
        if (ToulBar2::verbose >= 1)
            cout << ") NOT FOUND for cluster " << cluster->getId() << endl;
        return false;
    }
}

BigInteger Separator::getSg(Cost& res, BigInteger& nb)
{
    int i = 0;
    res = MIN_COST;
    if (ToulBar2::verbose >= 1)
        cout << "( ";
    TVars::iterator it = vars.begin();
    while (it != vars.end()) {
        assert(wcsp->assigned(*it));
        tValue val = wcsp->toIndex(*it, wcsp->getValue(*it));
        if (ToulBar2::verbose >= 1)
            cout << "(" << *it << "," << val << ") ";
        t[i] = val; // build the tuple
        res -= delta[i][val]; // delta structure
        ++it;
        i++;
    }
    TSGoods::iterator itsg = sgoods.find(t);
    if (itsg != sgoods.end()) {
        TPairSG p = itsg->second;
        if (ToulBar2::verbose >= 1)
            cout << ") Use #good  with nb = " << p.second << "solutions on cluster " << cluster->getId() << endl;
        /*		assert(res + p.first >= MIN_COST);
                res += p.first;*/
        nb = p.second;
        /*res = MAX(MIN_COST,res);*/
        return nb;
    } else {
        /*res = MIN_COST;*/
        if (ToulBar2::verbose >= 1)
            cout << ") NOT FOUND for cluster " << cluster->getId() << endl;
        return nb = -1;
    }
}

bool Separator::getF(bool& free)
{
    int i = 0;
    if (ToulBar2::verbose >= 1)
        cout << "( ";
    TVars::iterator it = vars.begin();
    while (it != vars.end()) {
        assert(wcsp->assigned(*it));
        tValue val = wcsp->toIndex(*it, wcsp->getValue(*it));
        if (ToulBar2::verbose >= 1)
            cout << "(" << *it << "," << val << ") ";
        t[i] = val; // build the tuple
        ++it;
        i++;
    }
    TFrees::iterator itsg = frees.find(t);
    if (itsg != frees.end()) {
        if (ToulBar2::verbose >= 1)
            cout << ") Use freedom with value = " << itsg->second << " on cluster " << cluster->getId() << endl;
        free = itsg->second;
        return true;
    } else {
        if (ToulBar2::verbose >= 1)
            cout << ") freedom NOT FOUND for cluster " << cluster->getId() << endl;
        return false;
    }
}

bool Separator::solGet(TAssign& a, Tuple& sol, bool& free)
{
    int i = 0;
    TVars::iterator it = vars.begin();
    while (it != vars.end()) {
        tValue val = wcsp->toIndex(*it, a[*it]);
        t[i] = val; // build the tuple
        ++it;
        i++;
    }
    TPairSol p;
    TSols::iterator itsol = solutions.find(t);
    if (itsol != solutions.end()) {
        p = itsol->second;
        sol = p.second;

        if (ToulBar2::verbose >= 1) {
            cout << "asking  solution  sep:";
            cout << t;
            cout << "  cost: " << p.first << endl;
            cout << "  sol: " << sol << endl;
        }

        // update the freedom status
        TFrees::iterator itsg = freesSol.find(t);
        assert(itsg != freesSol.end());
        free = itsg->second;

        return true;
    }
    return false;
}

void Separator::solRec(Cost ub)
{
    WCSP* wcsp = cluster->getWCSP();

    Cost deltares = MIN_COST;
    int i = 0;
    TVars::iterator it = vars.begin();
    while (it != vars.end()) {
        assert(wcsp->assigned(*it));
        tValue val = wcsp->toIndex(*it, wcsp->getValue(*it));
        t[i] = val; // build the tuple
        deltares += delta[i][val];
        ++it;
        i++;
    }

    //  	TPairSol p;
    //  	TSols::iterator itsol = solutions.find(t);
    //  	if(itsol != solutions.end()) {
    //  		p = itsol->second;
    //  	    assert(p.first > ub + deltares); // previous known solution must be worse
    //  	}

    wcsp->restoreSolution(cluster);

    TVars::iterator iter_begin, iter_end;

    if (cluster->getFreedom()) {
        iter_begin = cluster->beginVarsTree();
        iter_end = cluster->endVarsTree();
    } else {
        iter_begin = cluster->beginVars();
        iter_end = cluster->endVars();
    }

    i = 0;
    it = iter_begin;
    while (it != iter_end) {
        assert(wcsp->assigned(*it));
        if (!cluster->isSepVar(*it)) {
            tValue val = wcsp->toIndex(*it, wcsp->getValue(*it));
            s[i] = val;
            i++;
        }
        ++it;
    }

    solutions[t] = TPairSol(ub + deltares, Tuple(s.begin(), s.begin() + i)); // remember only proper variables
    freesSol[t] = cluster->getFreedom();

    if (ToulBar2::verbose >= 1) {
        cout << "recording solution  "
             << " cost: " << ub << " + delta: " << deltares;
        cout << " sol: " << s << " sep: " << t << endl;
    }
}

void Separator::resetLb()
{
    TNoGoods::iterator it = nogoods.begin();
    while (it != nogoods.end()) {
        (it->second).first = MIN_COST;
        (it->second).third = Solver::OpenList();
        ++it;
    }
}

void Separator::resetUb()
{
    TNoGoods::iterator it = nogoods.begin();
    while (it != nogoods.end()) {
        (it->second).second = MAX_COST;
        (it->second).third = Solver::OpenList();
        ++it;
    }
}

void Separator::print(ostream& os)
{
    os << this << " nogoods(";
    Double totaltuples = 1;
    for (int i = 0; i < arity_; i++) {
        os << scope[i]->getName();
        if (i < arity_ - 1)
            os << ",";
        totaltuples = totaltuples * scope[i]->getDomainInitSize();
    }
    os << ")    ";
    os << " |nogoods| = " << nogoods.size() << " / " << totaltuples << " min:" << ((nogoods.size() == 0) ? MIN_COST : min_element(nogoods.begin(), nogoods.end(), nogoods.value_comp())->second.first) << " (" << cluster->getNbBacktracksClusterTree() << " bt)";
    if (ToulBar2::verbose >= 4) {
        os << "nogoods: {";
        TNoGoods::iterator it = nogoods.begin();
        while (it != nogoods.end()) {
            TPairNG p = it->second;
            os << "<";
            for (unsigned int i = 0; i < it->first.size(); i++) {
                os << it->first[i];
                if (i < it->first.size() - 1)
                    os << " ";
            }
            os << "," << p.first << ">";
            if (it != nogoods.end())
                os << " ";
            it++;
        }
        os << "} " << endl;
    }
    os << endl;
}

/*
 * Cluster class
 *
 */

Cluster::Cluster(TreeDecomposition* tdin)
    : td(tdin)
    , wcsp(tdin->getWCSP())
    , id(-1)
    , parent(NULL)
    , sep(NULL)
    , lb(MIN_COST)
    , ub(MAX_COST)
    , lbRDS(MIN_COST)
    , active(true)
    , countElimVars(1)
    , num_part(-1)
    , freedom_on(false)
    , isCurrentlyInTD(true)
    , depth(-1)
    , cp(NULL)
    , open(NULL)
    , hbfsGlobalLimit(LONGLONG_MAX)
    , hbfsLimit(LONGLONG_MAX)
    , nbBacktracks(0)
{
    instance = clusterCounter++;
}

Cluster::~Cluster()
{
    id = -1;
    if (sep) {
        sep->deconnect();
        if (sep->isInQueueSep()) {
            sep->unqueueSep();
        }
        // Notice: separators will be deleted inside ~WCSP with the other constraints
    }
    if (cp) {
        delete cp;
    }
    // Notice: do not need to delete open because it is done in ~Solver for the root cluster only (the other clusters have a pointer to nogoods'open field, an object not created by new!)
}

void Cluster::addVar(Variable* x) { vars.insert(x->wcspIndex); }
void Cluster::removeVar(Variable* x) { vars.erase(x->wcspIndex); }

void Cluster::addVars(TVars& morevars)
{
    vars.insert(morevars.begin(), morevars.end());
}

void Cluster::addCtr(Constraint* c) { ctrs.insert(c); }
void Cluster::removeCtr(Constraint* c) { ctrs.erase(c); }
void Cluster::clearCtrs() { ctrs.clear(); }

void Cluster::addEdge(Cluster* c) { edges.insert(c); }

TClusters::iterator Cluster::removeEdge(TClusters::iterator it)
{
    TClusters::iterator itaux = it;
    ++it;
    edges.erase(itaux);
    return it;
}

void Cluster::removeEdge(Cluster* c)
{
    TClusters::iterator it = edges.find(c);
    if (it != edges.end())
        edges.erase(it);
}

void Cluster::addEdges(TClusters& cls)
{
    edges.insert(cls.begin(), cls.end());
}

void Cluster::addCtrs(TCtrs& ctrsin)
{
    ctrs.insert(ctrsin.begin(), ctrsin.end());
}

TCtrs Cluster::getCtrsTree()
{
    TCtrs ctrsTree;
    for (TClusters::iterator it = descendants.begin(); it != descendants.end(); ++it) {
        Cluster* c = *it;
        ctrsTree.insert(c->getCtrs().begin(), c->getCtrs().end());
    }
    return ctrsTree;
}

bool compareCtrs(const Constraint* lhs, const Constraint* rhs)
{
    assert(lhs);
    assert(rhs);
    return lhs != rhs;
}

void Cluster::deconnectDiff(TCtrs& listCtrsTot, TCtrs& listCtrs)
{
    TCtrs listDiff;
    set_difference(listCtrsTot.begin(), listCtrsTot.end(), listCtrs.begin(), listCtrs.end(), inserter(listDiff, listDiff.begin()), compareCtrs);
    for (TCtrs::iterator itctr = listDiff.begin(); itctr != listDiff.end(); ++itctr) {
        Constraint* ctr = *itctr;
        ctr->deconnect();
    }
}

void Cluster::deconnectSep()
{
    if (!sep)
        return;
    TVars::iterator its = beginSep();
    while (its != endSep()) {
        Variable* x = wcsp->getVar(*its);
        ConstraintList* xctrs = x->getConstrs();
        for (ConstraintList::iterator it = xctrs->begin(); it != xctrs->end(); ++it) {
            Constraint* ctr = (*it).constr;
            Cluster* ctrc = td->getCluster(ctr->getCluster()); // warning! return parent cluster if sep
            if (!(ctr->isSep() && isDescendant(ctrc))) {
                // keep descendant separators connected
                //			    if (ctr->isSep()) cout << "deconnect separator parent " << ctr->cluster << " " << *ctr << endl;
                ctr->deconnect();
            }
        }
        x->assign(x->getSupport());
        ++its;
    }
}

void Cluster::resetLbRec()
{
    if (sepSize() > 0)
        sep->resetLb();
    if (this != td->getRoot())
        open = NULL;
    for (TClusters::iterator iter = beginEdges(); iter != endEdges(); ++iter) {
        (*iter)->resetLbRec();
    }
}

void Cluster::resetUbRec(Cluster* root)
{
    if (!sep || sepSize() == 0 || !root->sep || root->sepSize() == 0)
        return;
    TVars inter;
    td->intersection(sep->getVars(), root->sep->getVars(), inter);
    if (inter.size() > 0)
        sep->resetUb();
    if (this != td->getRoot())
        open = NULL;
    for (TClusters::iterator iter = beginEdges(); iter != endEdges(); ++iter) {
        (*iter)->resetUbRec(root);
    }
}

Cost Cluster::getLbRec() const
{
    assert(isActive());
    Cost res = lb;
    for (TClusters::iterator iter = beginEdges(); iter != endEdges(); ++iter) {
        if ((*iter)->isActive())
            res += (*iter)->getLbRec();
    }
    return res;
}

Cost Cluster::getLbRecRDS()
{
    assert(isActive());
    Cost res = lb;
    for (TClusters::iterator iter = beginEdges(); iter != endEdges(); ++iter) {
        if ((*iter)->isActive()) {
            Cost propa = (*iter)->getLbRecRDS();
            Cost rds = (*iter)->getLbRDS();
            res += MAX(propa, rds);
        }
    }
    return res;
}

void Cluster::reactivate()
{
    assert(!isActive());
    active = true;
    for (TClusters::iterator iter = beginEdges(); iter != endEdges(); ++iter) {
        assert(!(*iter)->isActive());
        if (!(*iter)->sep->used())
            (*iter)->reactivate();
    }
}

void Cluster::deactivate()
{
    if (isActive()) {
        if (ToulBar2::verbose >= 1)
            cout << "deactive cluster " << getId() << endl;
        active = false;
        for (TClusters::iterator iter = beginEdges(); iter != endEdges(); ++iter) {
            (*iter)->deactivate();
        }
    }
}

void Cluster::setWCSP2Cluster()
{
    for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
        if (!isVar(i)) {
            EnumeratedVariable* x = (EnumeratedVariable*)wcsp->getVar(i);
            for (ConstraintList::iterator it = x->getConstrs()->begin(); it != x->getConstrs()->end(); ++it) {
                Constraint* ctr = (*it).constr;
                ctr->deconnect();
            }
        }
    }
    for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
        if (!isVar(i)) {
            EnumeratedVariable* x = (EnumeratedVariable*)wcsp->getVar(i);
            x->assign(x->getSupport());
        }
    }
}

void Cluster::getElimVarOrder(vector<int>& elimVarOrder)
{
    for (TClustersSorted::reverse_iterator iter = sortedEdges.rbegin(); iter != sortedEdges.rend(); ++iter) {
        Cluster* cluster = *iter;
        cluster->getElimVarOrder(elimVarOrder);
    }
    for (TVarsSorted::reverse_iterator itp = sortedVars.rbegin(); itp != sortedVars.rend(); ++itp) {
        if (!isSepVar(*itp)) {
            elimVarOrder.push_back(*itp);
        }
    }
}

// side-effect: remember last solution
void Cluster::getSolution(TAssign& sol)
{
    static Tuple s; // FIXME: unsafe???

    TVars::iterator it, iter_begin, iter_end;

    bool free = getFreedom();

    if (getParent() == NULL || this == td->getRootRDS()) {
        if (vars.size() == 0) {
            if (free) {
                iter_begin = beginVarsTree();
                iter_end = endVarsTree();

                it = iter_begin;
                while (it != iter_end) {
                    assert(wcsp->assigned(*it));
                    sol[*it] = wcsp->getValue(*it);
                    ++it;
                }
            }
        } else {
            if (free) {
                iter_begin = beginVarsTree();
                iter_end = endVarsTree();
            } else {
                iter_begin = beginVars();
                iter_end = endVars();
            }
            it = iter_begin;
            while (it != iter_end) {
                assert(wcsp->assigned(*it));
                sol[*it] = wcsp->getValue(*it);
                ++it;
            }
        }
    }
    if (sep) {
#ifndef NDEBUG
        bool found = sep->solGet(sol, s, free);
        assert(found);
#else
        sep->solGet(sol, s, free);
#endif

        if (free) {
            iter_begin = beginVarsTree();
            iter_end = endVarsTree();
        } else {
            iter_begin = beginVars();
            iter_end = endVars();
        }

        int i = 0;
        it = iter_begin;
        while (it != iter_end) {
            if (!isSepVar(*it)) {
                sol[*it] = wcsp->toValue(*it, s[i]);
                //				cout << *it << " := " << sol[*it] << endl;
                if (!ToulBar2::verifyOpt && ToulBar2::solutionBasedPhaseSaving)
                    wcsp->setBestValue(*it, sol[*it]);
                i++;
            }
            ++it;
        }
    }

    if (!free) {
        for (TClusters::iterator iter = beginEdges(); iter != endEdges(); ++iter) {
            Cluster* cluster = *iter;
            if (ToulBar2::bilevel && td->getRoot() == this && cluster == *rbeginEdges())
                break; // Do not reconstruct solution for NegProblem2
            cluster->getSolution(sol);
        }
    }
}

bool Cluster::isEdge(Cluster* c)
{
    for (TClusters::iterator iter = beginEdges(); iter != endEdges(); ++iter) {
        Cluster* cluster = *iter;
        if (c == cluster)
            return true;
    }
    return false;
}

void Cluster::setup()
{
    if (sep)
        sep->setup(this);
    if (ToulBar2::hbfs) {
        if (cp)
            delete cp;
        cp = new Solver::CPStore();
    }
}

void Cluster::accelerateDescendants()
{
    quickdescendants = vector<bool>(td->getNbOfClusters(), false);
    for (auto itc = beginDescendants(); itc != endDescendants(); ++itc) {
        Cluster* c = *itc;
        quickdescendants[c->getId()] = true;
    }
}

void Cluster::accelerateIntersections()
{
    quickIntersections.clear();
    for (TClusters::iterator itc = beginEdges(); itc != endEdges(); ++itc) {
        Cluster* cj = *itc;
        TVars cjsep;
        td->intersection(getVars(), cj->getVars(), cjsep);
        quickIntersections[cj] = cjsep;
    }
}

void Cluster::quickIntersection(Cluster* cj, TVars& cjsep)
{
    map<Cluster*, TVars>::iterator itcj = quickIntersections.find(cj);
    if (itcj != quickIntersections.end()) {
        cjsep = itcj->second;
    } else {
        TVars cjsep2;
        td->intersection(getVars(), cj->getVars(), cjsep2);
        quickIntersections[cj] = cjsep2; // warning! copy constructor
        cjsep = quickIntersections[cj]; // warning! copy by reference only
    }
}

void Cluster::print()
{
    // cout << "(" << id << ",n:" << getNbVars() << ",lb:" << getLb() << ") ";

    cout << "cluster " << getId();

    cout << " vars {";
    TVars::iterator itp = beginVars();
    while (itp != endVars()) {
        if (!isSepVar(*itp)) {
            cout << wcsp->getVar(*itp)->getName() << ",";
            // cout << *itp << "C" << wcsp->getVar(*itp)->getCluster() << ",";
            assert(wcsp->getVar(*itp)->getCluster() == -1 || wcsp->getVar(*itp)->getCluster() == getId());
        }
        ++itp;
    }
    cout << "\b}";

    if (sep) {
        cout << " U sep {";
        TVars::iterator its = beginSep();
        while (its != endSep()) {
            cout << wcsp->getVar(*its)->getName();
            ++its;
            if (its != endSep())
                cout << ",";
        }
        cout << "}";
    }

    if (!edges.empty()) {
        cout << " sons {";
        if (sortedEdges.size() == edges.size()) {
            TClusters::iterator itc = beginSortedEdges();
            while (itc != endSortedEdges()) {
                cout << (*itc)->getId();
                ++itc;
                if (itc != endSortedEdges())
                    cout << ",";
            }
        } else {
            TClusters::iterator itc = beginEdges();
            while (itc != endEdges()) {
                cout << (*itc)->getId();
                ++itc;
                if (itc != endEdges())
                    cout << ",";
            }
        }
        cout << "}";
    }

    /*
        cout << " ctrs {";
        TCtrs::iterator itctr = beginCtrs();
        while(itctr != endCtrs()) {
          Constraint* ctr = *itctr;
          cout << "( ";
          for(int i=0;i<ctr->arity();i++) cout << ctr->getVar(i)->wcspIndex << " ";
          cout << ">C" << ctr->getCluster() << ")";
          ++itctr;
        }
        cout << "}";
        cout << " descendants {";
        TClusters::iterator itd = beginDescendants();
        while(itd != endDescendants()) {
      cout << (*itd)->getId();
      ++itd;
      if(itd != endDescendants()) cout << ",";
        }
        cout << "}";
     */

    cout << endl;
}

void Cluster::dump()
{
    // cout << "(" << id << ",n:" << getNbVars() << ",lb:" << getLb() << ") ";

    char clusterVarsFilename[128];
    char sepVarsFilename[128];
    char sonsFilename[128];
    char fatherFilename[128];
    char sepSizeFilename[128];

    sprintf(clusterVarsFilename, "%s.info/%d.vars", getWCSP()->getName().c_str(), getId());
    sprintf(sepVarsFilename, "%s.info/%d.sep", getWCSP()->getName().c_str(), getId());
    sprintf(sonsFilename, "%s.info/%d.sons", getWCSP()->getName().c_str(), getId());
    sprintf(fatherFilename, "%s.info/%d.father", getWCSP()->getName().c_str(), getId());
    sprintf(sepSizeFilename, "%s.info/%d.sepsize", getWCSP()->getName().c_str(), getId());

    ofstream clusterVarsFile(clusterVarsFilename);
    ofstream sepVarsFile(sepVarsFilename);
    ofstream sonsFile(sonsFilename);
    ofstream fatherFile(fatherFilename);
    ofstream sepSizeFile(sepSizeFilename);

    if (parent) {
        fatherFile << parent->getId();
    } else {
        fatherFile << "-1";
    }
    fatherFile.close();

    long double separatorSize = 1.0;

    if (sep) {
        TVars::iterator its = beginSep();
        while (its != endSep()) {
            clusterVarsFile << wcsp->getVar(*its)->getName() << " ";
            sepVarsFile << wcsp->getVar(*its)->getName() << " ";
            separatorSize *= (1.0 + wcsp->getVar(*its)->getSup() - wcsp->getVar(*its)->getInf());
            ++its;
        }
    }
    sepSizeFile << separatorSize;

    TVars::iterator itp = beginVars();
    while (itp != endVars()) {
        if (!isSepVar(*itp)) {
            clusterVarsFile << wcsp->getVar(*itp)->getName() << " ";
            assert(wcsp->getVar(*itp)->getCluster() == -1 || wcsp->getVar(*itp)->getCluster() == getId());
        }
        ++itp;
    }

    if (getNbVars() == 0)
        clusterVarsFile << " ";

    if (!edges.empty()) {
        TClusters::iterator itc = beginEdges();
        while (itc != endEdges()) {
            sonsFile << (*itc)->getId();
            ++itc;
            if (itc != endEdges())
                sonsFile << " ";
        }
    }

    clusterVarsFile.close();
    sepVarsFile.close();
    sonsFile.close();
    sepSizeFile.close();
}

void Cluster::cartProduct(BigInteger& prodCart)
{
    for (TVars::iterator it = varsTree.begin(); it != varsTree.end(); it++) {
        Variable* x = (Variable*)wcsp->getVar(*it);
        prodCart *= x->getDomainSize();
    }
}

/*
 * Tree Decomposition class
 *
 */

TreeDecomposition::TreeDecomposition(WCSP* wcsp_in)
    : wcsp(wcsp_in)
    , rootRDS(NULL)
    , currentCluster(-1)
    , deltaModified(vector<StoreInt>(wcsp_in->numberOfVariables(), StoreInt(false)))
    , max_depth(-1)
{
}

TreeDecomposition::~TreeDecomposition()
{
    for(auto& c: clusters) {
        if (c) {
            delete c;
        }
    }
}

bool TreeDecomposition::isInCurrentClusterSubTree(int idc)
{
    if (idc < 0)
        return false;
    Cluster* ci = getCurrentCluster();
    Cluster* cj = getCluster(idc);
    assert(ci->isActive());
    return ci->isDescendant(cj);
}

bool TreeDecomposition::isActiveAndInCurrentClusterSubTree(int idc)
{
    if (idc < 0)
        return false;
    Cluster* ci = getCurrentCluster();
    Cluster* cj = getCluster(idc);
    assert(ci->isActive());
    if (!cj->isActive())
        return false;
    else
        return ci->isDescendant(cj);
}

void TreeDecomposition::fusion(Cluster* ci, Cluster* cj)
{
    if (!ci)
        return;
    if (!cj)
        return;

    if (ToulBar2::verbose >= 1)
        cout << "fusion: " << ci->getId() << " " << cj->getId() << endl;

    ci->addVars(cj->getVars());
    ci->addCtrs(cj->getCtrs());
    ci->addEdges(cj->getEdges());
    TClusters::iterator itk = cj->beginEdges();
    while (itk != cj->endEdges()) {
        Cluster* ck = *itk;
        ++itk;
        ck->removeEdge(cj);
        ck->addEdge(ci);
    }
    ci->removeEdge(ci);
    clusters[cj->getId()] = NULL;
    if (ToulBar2::verbose >= 1) {
        cout << "fusion ci " << ci->getId() << ",  cj " << cj->getId() << endl;
        ci->print();
    }
    delete cj;
}

bool TreeDecomposition::treeFusion()
{
    bool done = false;
    for (int j = clusters.size() - 1; j >= 0; j--) {
        if (!clusters[j])
            continue;
        Cluster* cj = clusters[j];
        if (ToulBar2::verbose >= 3) {
            cout << "fusion testing ";
            cj->print();
        }

        TClusters::iterator it = cj->beginEdges();
        while (it != cj->endEdges()) {
            Cluster* c = *it;
            assert(c == clusters[c->getId()]);

            if ((c->getId() < cj->getId()) && (included(c->getVars(), cj->getVars()) || included(cj->getVars(), c->getVars()))) {
                c->addVars(cj->getVars());
                c->addCtrs(cj->getCtrs());
                c->addEdges(cj->getEdges());
                TClusters::iterator itk = cj->beginEdges();
                while (itk != cj->endEdges()) {
                    Cluster* ck = *itk;
                    ck->removeEdge(cj);
                    ck->addEdge(c);
                    ++itk;
                }
                c->removeEdge(c);
                clusters[cj->getId()] = NULL;
                if (ToulBar2::verbose >= 1) {
                    cout << "fusion ci " << c->getId() << ",  cj " << cj->getId() << endl;
                    c->print();
                }
                delete cj;
                //					done = true;
                break;
            }
            ++it;
        }
    }
    return done;
}

void TreeDecomposition::treeFusions()
{
    while (treeFusion())
        ;

    TClusters visited;
    //	Cluster* croot = getBiggerCluster(visited);
    //	heuristicFusionRec(croot, croot);

    int treewidth = 0;
    TClusters sclu;
    for (unsigned int i = 0; i < clusters.size(); i++) {
        if (clusters[i]) {
            Cluster* c = clusters[i];
            sclu.insert(c);
            if (c->getNbVars() > treewidth)
                treewidth = c->getNbVars();
        }
    }
    int i = 0;
    clusters.clear();
    TClusters::iterator it = sclu.begin();
    while (it != sclu.end()) {
        Cluster* c = *it;
        c->setId(i++);
        clusters.push_back(*it);
        ++it;
    }
    if (ToulBar2::verbose >= 2)
        cout << "Tree decomposition width  : " << treewidth - 1 << endl;
}

void TreeDecomposition::pathFusions(vector<int>& order)
{
    vector<Cluster*> rds;
    int size = clusters.size();
    vector<bool> connected;

    // detect singleton variables
    for (int i = 0; i < size; i++) {
        bool isconnected = (clusters[i]->getNbVars() > 1);
        for (int j = 0; j < i; j++) {
            if (clusters[j]->isVar(order[i])) {
                isconnected = true;
                break;
            }
        }
        connected.push_back(isconnected);
    }

    for (int i = 0; i < size; i++) {
        assert(clusters[i] && clusters[i]->isVar(order[i]));
        Cluster* c = new Cluster(this);
        c->addVar(wcsp->getVar(order[i]));
        if (connected[i]) {
            for (int j = 0; j < i; j++) {
                if (clusters[j]->isVar(order[i])) {
                    for (int l = j + 1; l < i; l++) {
                        if (connected[l])
                            rds[l]->addVar(wcsp->getVar(order[j]));
                    }
                    c->addVar(wcsp->getVar(order[j]));
                }
            }
            if (rds.size() > 0) {
                int last = rds.size() - 1;
                while (last >= 0 && !connected[last])
                    last--;
                if (last >= 0) {
                    c->addEdge(rds[last]);
                    rds[last]->addEdge(c);
                }
            }
        }
        rds.push_back(c);
    }

    // fusion on included clusters
    for (int i = 0; i < size; i++) {
        if (i < size - 1 && included(rds[i]->getVars(), rds[i + 1]->getVars())) {
            rds[i + 1]->removeEdge(rds[i]);
            rds[i]->removeEdge(rds[i + 1]);
            if (rds[i]->getEdges().size() > 0) {
                assert(rds[i]->getEdges().size() == 1);
                TClusters::iterator it = rds[i]->beginEdges();
                rds[i + 1]->addEdge(*it);
                (*it)->removeEdge(rds[i]);
                (*it)->addEdge(rds[i + 1]);
            }
            delete rds[i];
            rds[i] = NULL;
        }
    }
    for(auto&c : clusters) {
        delete c;
    }
    clusters.clear();
    for (int i = 0; i < size; i++) {
        Cluster* c = rds[i];
        if (c) {
            c->setId(clusters.size());
            clusters.push_back(c);
        }
    }
}

// Minimize tree height when the separators are included
// path contains the ordered list of clusters from the root to the parent cluster of c
// time complexity in O(n) using DFS
void TreeDecomposition::reduceHeight(Cluster* c, vector<Cluster*> path)
{
    Cluster* cparent = NULL;
    if (path.size() > 0)
        cparent = path.back();
    assert(c != cparent);
    TClusters::iterator itj;
    itj = c->beginEdges();
    while (itj != c->endEdges()) {
        Cluster* cj = *itj;
        ++itj; // warning! done before removing the corresponding edge (c -> cj) so that the iterator remains valid
        if (cj != cparent) { // warning! the current tree decomposition is undirected
            TVars cjsep;
            intersection(c, cj, cjsep);
            if (cparent && included(cjsep, cparent->getVars())) {
                // reconnect the child cluster closest to the root of the tree decomposition following path information
                int pos = path.size() - 1;
                assert(pos >= 0 && path[pos] == cparent);
                while (pos >= 1 && included(cjsep, path[pos - 1]->getVars())) {
                    pos--;
                }
                assert(pos >= 0 && path[pos] != c);
                //                cout << "move " << cj->getId() << " from " << c->getId() << " to " << path[pos]->getId() << endl;
                //                cout << "with path:";
                //                for (int i=0; i<path.size(); i++) cout << " " << path[i]->getId();
                //                cout << endl;
                c->removeEdge(cj);
                path[pos]->addEdge(cj);
                cj->removeEdge(c);
                cj->addEdge(path[pos]);
                reduceHeight(cj, vector<Cluster*>(path.begin(), path.begin() + pos + 1)); // continue recursively on cj with an updated path to the root
            } else if (!cparent && cjsep.size() == 0) { // warning! it is done before a meta-root is created with connected components as its children
                c->removeEdge(cj);
                cj->removeEdge(c);
                reduceHeight(cj, vector<Cluster*>()); // cj becomes a root
            } else {
                vector<Cluster*> newpath = path;
                newpath.push_back(c); // add current cluster c in the path as father of cj
                reduceHeight(cj, newpath); // continue recursively on cj
            }
        }
    }
}

int TreeDecomposition::getNextUnassignedVar(TVars* vars)
{
    return *(vars->begin());
}

int TreeDecomposition::getVarMinDomainDivMaxWeightedDegree(TVars* vars)
{
    int varIndex = -1;
    Cost worstUnaryCost = MIN_COST;
    double best = MAX_VAL - MIN_VAL;

    for (TVars::iterator iter = vars->begin(); iter != vars->end(); ++iter) {
        Cost unarymediancost = MIN_COST;
        int domsize = wcsp->getDomainSize(*iter);
        if (ToulBar2::weightedTightness) {
            ValueCost array[domsize];
            wcsp->getEnumDomainAndCost(*iter, array);
            unarymediancost = stochastic_selection<ValueCost>(array, 0, domsize - 1, domsize / 2).cost;
        }
        double heuristic = (double)domsize / (double)(wcsp->getWeightedDegree(*iter) + 1 + unarymediancost);
        if (varIndex < 0 || heuristic < best - (double)ToulBar2::epsilon * best
            || (heuristic < best + (double)ToulBar2::epsilon * best && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
            best = heuristic;
            varIndex = *iter;
            worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
        }
    }
    return varIndex;
}

void TreeDecomposition::splitClusterRec(Cluster* c, Cluster* father, unsigned int maxsize, TClusters& unvisited)
{
    TVars cvars = c->getVars();
    //  cout << c->getId() << " " << cvars.size() << endl;
    TVars csep;
    if (father) {
        intersection(father, c, csep);
    }
    TVars cproper;
    difference(cvars, csep, cproper);
    if (cproper.size() > maxsize && (!father || c->getEdges().size() != 1)) {
        Cluster* cprev = NULL;
        TClusters cedges = c->getEdges();
        if (father)
            cedges.erase(father);
        while (cproper.size() > 0) {
            TVars cnewvars;
            int varIndex = ((ToulBar2::Static_variable_ordering) ? getNextUnassignedVar(&cproper) : getVarMinDomainDivMaxWeightedDegree(&cproper));
            for (unsigned int i = 0; i < maxsize && varIndex >= 0; i++) {
                cnewvars.insert(varIndex);
                cproper.erase(varIndex);
                varIndex = ((ToulBar2::Static_variable_ordering) ? getNextUnassignedVar(&cproper) : getVarMinDomainDivMaxWeightedDegree(&cproper));
            }
            if (!cprev) {
                c->getVars().clear();
                c->addVars(csep);
                c->addVars(cnewvars);
                c->getEdges().clear();
                if (father)
                    c->addEdge(father);
                cprev = c;
            } else {
                Cluster* cnew = new Cluster(this);
                cnew->setId(clusters.size());
                clusters.push_back(cnew);
                unvisited.insert(cnew);
                cnew->addVars(cprev->getVars());
                cnew->addVars(cnewvars);
                cnew->addEdge(cprev);
                cprev->addEdge(cnew);
                cprev = cnew;
            }
        }
        assert(cprev->getEdges().size() == 1);
        father = *(cprev->beginEdges());
        cprev->addEdges(cedges);
        TClusters::iterator itj = cprev->beginEdges();
        while (itj != cprev->endEdges()) {
            Cluster* cj = *itj;
            if (cj != father) {
                cj->removeEdge(c);
                cj->addEdge(cprev);
            }
            ++itj;
        }
        c = cprev;
    }
    TClusters::iterator itj = c->beginEdges();
    while (itj != c->endEdges()) {
        Cluster* cj = *itj;
        if (cj != father)
            splitClusterRec(cj, c, maxsize, unvisited);
        ++itj;
    }
}

TVars TreeDecomposition::boostingVarElimRec(Cluster* c, Cluster* father, Cluster* grandfather, unsigned int maxsize, TClusters& unvisited)
{
    TVars addedVarBySons;
    TClusters::iterator itj = c->beginEdges();
    while (itj != c->endEdges()) {
        Cluster* cj = *itj;
        ++itj; // warning! must be done before going inside boostingVarElimRec as it can delete current cluster/iterator by the following removeEdge(c) operation
        if (cj != father) {
            TVars cjaddedvars;
            cjaddedvars = boostingVarElimRec(cj, c, father, maxsize, unvisited);
            sum(addedVarBySons, cjaddedvars);
        }
    }
    if (father && c->getEdges().size() == 1) {
        TVars fathersep;
        if (grandfather) {
            intersection(father, grandfather, fathersep);
        }
        TVars cvars;
        sum(fathersep, addedVarBySons);
        difference(c->getVars(), fathersep, cvars);
        if (cvars.size() <= maxsize) {
            //	  	  cout << c->getId() << " which has " << cvars.size() << " vars (except whose from " << ((grandfather)?grandfather->getId():-1) << ") is merged into " << father->getId() << endl;
            TVars csep;
            intersection(c, father, csep);
            TVars cproper;
            difference(c->getVars(), csep, cproper);
            father->addVars(cproper);
            father->removeEdge(c);
            father->accelerateIntersections();
            unvisited.erase(c);
            clusters.back()->setId(c->getId());
            clusters[c->getId()] = clusters.back();
            clusters.pop_back();
            sum(addedVarBySons, cproper);
            delete c;
        }
    }
    return addedVarBySons;
}

void TreeDecomposition::mergeClusterRec(Cluster* c, Cluster* father, unsigned int maxsepsize, unsigned int minpropervar, TClusters& unvisited)
{
    TClusters::iterator itj = c->beginEdges();
    while (itj != c->endEdges()) {
        Cluster* cj = *itj;
        ++itj; // warning! must be done before going inside boostingVarElimRec as it can delete current cluster/iterator by the following removeEdge(c) operation
        if (cj != father) {
            mergeClusterRec(cj, c, maxsepsize, minpropervar, unvisited);
        }
    }
    if (father) {
        TVars csep;
        intersection(c, father, csep);
        assert(csep.size() > 0);
        if ((csep.size() > maxsepsize) || (c->getVars().size() - csep.size() < minpropervar)) {
            father->addVars(c->getVars());
            father->addEdges(c->getEdges());
            TClusters::iterator itk = c->beginEdges();
            while (itk != c->endEdges()) {
                Cluster* ck = *itk;
                ++itk;
                ck->removeEdge(c);
                ck->addEdge(father);
            }
            father->removeEdge(father);
            father->removeEdge(c);
            father->accelerateIntersections();
            unvisited.erase(c);
            clusters.back()->setId(c->getId());
            clusters[c->getId()] = clusters.back();
            clusters.pop_back();
            delete c;
        }
    }
}

// Specific code for cluster fusion based on heuristic criteria
void TreeDecomposition::heuristicFusionRec(Cluster* c, Cluster* noc)
{
    TClusters::iterator it = c->beginEdges();
    while (it != c->endEdges()) {
        Cluster* cj = *it;
        ++it;
        if (cj == c)
            continue;
        if (cj == noc)
            continue;
        heuristicFusionRec(cj, c);
    }

    it = c->beginEdges();
    while (it != c->endEdges()) {
        Cluster* cj = *it;
        ++it;
        if (cj == c)
            continue;
        if (cj == noc)
            continue;
        TVars varsum;
        TVars varinter;
        sum(c->getVars(), cj->getVars(), varsum);
        intersection(c->getVars(), cj->getVars(), varinter);

        int dif = 2;
        bool bf1 = (varinter.size() > 2) && (varsum.size() <= (unsigned int)c->getNbVars() + dif);
        bool bf2 = (varinter.size() > 2) && (varsum.size() <= (unsigned int)cj->getNbVars() + dif);
        bool bf3 = (varinter.size() > 100);
        if (bf1 || bf2 || bf3) {
            fusion(c, cj);
        }
    }
}

Cluster* TreeDecomposition::getBiggerCluster(TClusters& unvisited)
{
    Cluster* cmax = NULL;
    int maxsize = 0;
    for (TClusters::iterator itc = unvisited.begin(); itc != unvisited.end(); ++itc) {
        Cluster* c = *itc;
        assert(c);
        if (c->getNbVars() > maxsize) {
            maxsize = c->getNbVars();
            cmax = c;
            if (ToulBar2::btdMode == 3)
                break;
        }
    }
    return cmax;
}

Cluster* TreeDecomposition::getClusterMinHeight(TClusters& unvisited)
{

    TClusters::iterator itc = unvisited.begin();
    Cluster* c_start = *itc;
    assert(c_start);
    int minheight = height(c_start);
    Cluster* cmin_height = c_start;
    ++itc;
    while (itc != unvisited.end()) {
        Cluster* c = *itc;
        assert(c);
        if (ToulBar2::reduceHeight) {
            reduceHeight(c, vector<Cluster*>());
        }
        if (height(c) < minheight) {
            minheight = height(c);
            cmin_height = c;
            if (ToulBar2::btdMode == 3)
                break;
        }
        ++itc;
    }
    return cmin_height;
}

Cluster* TreeDecomposition::getCluster_height_rootsize_max(TClusters& unvisited)
{
    Cluster* cmax_ratio = NULL;
    float maxratio = 0;
    for (TClusters::iterator itc = unvisited.begin(); itc != unvisited.end(); ++itc) {
        Cluster* c = *itc;
        assert(c);
        if (ToulBar2::reduceHeight) {
            reduceHeight(c, vector<Cluster*>());
        }
        if (float(c->getNbVars()) / float((height(c) - c->getNbVars())) >= maxratio) {
            maxratio = float(c->getNbVars()) / float(height(c) - c->getNbVars());
            cmax_ratio = c;
            if (ToulBar2::btdMode == 3)
                break;
        }
    }
    return cmax_ratio;
}

Cluster* TreeDecomposition::getCluster_height_rootsize_min(TClusters& unvisited)
{
    TClusters::iterator itc = unvisited.begin();
    Cluster* c_start = *itc;
    assert(c_start);
    float minratio = float(c_start->getNbVars()) / float(height(c_start) - c_start->getNbVars());
    Cluster* cmin_ratio = c_start;
    ++itc;
    while (itc != unvisited.end()) {
        Cluster* c = *itc;
        assert(c);
        if (ToulBar2::reduceHeight) {
            reduceHeight(c, vector<Cluster*>());
        }
        if (float(c->getNbVars()) / float((height(c) - c->getNbVars())) < minratio) {
            minratio = float(c->getNbVars()) / float(height(c) - c->getNbVars());
            cmin_ratio = c;
            if (ToulBar2::btdMode == 3)
                break;
        }
        ++itc;
    }
    return cmin_ratio;
}

int TreeDecomposition::height(Cluster* r, Cluster* father)
{
    int maxh = 0;
    TClusters::iterator it = r->beginEdges();
    while (it != r->endEdges()) {
        Cluster* adjr = *it;
        if (adjr != father) {
            int h = height(adjr, r);
            if (h > maxh)
                maxh = h;
        }
        ++it;
    }
    TVars rsep;
    intersection(r, father, rsep);
    return maxh + r->getNbVars() - rsep.size();
}

int TreeDecomposition::height(Cluster* r)
{
    int maxh = 0;
    TClusters::iterator it = r->beginEdges();
    while (it != r->endEdges()) {
        int h = height(*it, r);
        if (h > maxh)
            maxh = h;
        ++it;
    }
    return maxh + r->getNbVars();
}

void TreeDecomposition::makeDescendants(Cluster* c)
{
    c->getDescendants().insert(c);
    sum(c->getVarsTree(), c->getVars());
    TClusters::iterator itj = c->beginEdges();
    while (itj != c->endEdges()) {
        Cluster* cj = *itj;
        ++itj;
        makeDescendants(cj);
        clusterSum(c->getDescendants(), cj->getDescendants());
        sum(c->getVarsTree(), cj->getVarsTree());
    }
}

void TreeDecomposition::makeRootedRec(Cluster* c, Cluster* father, TClusters& unvisited)
{
    // cout << "makeRootedRec " << c << " (" << c->getId() << ") with father cluster " << father << " (" << ((father)?father->getId():NULL) << ")" << endl;
    TClusters::iterator itj = c->beginEdges();
    while (itj != c->endEdges()) {
        Cluster* cj = *itj;
        assert(cj != c);
        assert(cj != father);
        assert(cj->getId() >= 0 && clusters[cj->getId()] == cj); // must be a valid cluster inside the list of clusters
        // cout << "makeRootedRec " << c << " (" << c->getId() << ")" << " on child cluster " << cj << " (" << cj->getId() << ")" << endl;
        ++itj;
        cj->removeEdge(c);
        cj->setParent(c);
        unvisited.erase(cj);
        if (ToulBar2::searchMethod == DFBB) {
            TVars cjsep;
            intersection(c, cj, cjsep);

            //------- Add the constraint separator
            int i = 0;
            int arity = cjsep.size();
            EnumeratedVariable** scopeVars = new EnumeratedVariable*[arity]; // warning! it might allocate too much memory in the stack, it is safer to use the heap!
            TVars::iterator it = cjsep.begin();
            while (it != cjsep.end()) {
                scopeVars[i] = (EnumeratedVariable*)wcsp->getVar(*it);
                ++it;
                i++;
            }
            cj->setSep(new Separator(wcsp, scopeVars, arity));
            if (ToulBar2::approximateCountingBTD)
                cj->addCtr(cj->getSep());
            delete[] scopeVars;
            //-------
        }

        makeRootedRec(cj, c, unvisited);
    }
}

void TreeDecomposition::computeDepths(Cluster* c, int parent_depth)
{
    int c_depth = parent_depth + 1;
    if (c_depth > max_depth)
        max_depth = c_depth;
    c->setDepth(c_depth);
    for (TClusters::iterator iter = c->beginSortedEdges(); iter != c->endSortedEdges(); ++iter) {
        computeDepths(*iter, c_depth);
    }
}

void TreeDecomposition::DFSUtil(Cluster* c, cluster_visited& c_visited)
{
    c_visited[c] = true;
    comp.insert(c);
    for (TClusters::iterator itc = c->getEdges().begin(); itc != c->getEdges().end(); ++itc) {
        if (!c_visited[*itc]) {
            DFSUtil(*itc, c_visited);
        }
    }
}

int TreeDecomposition::connectedComponents()
{
    cluster_visited c_visited;
    // Mark all the clusters as not visited
    for (vector<Cluster*>::iterator itc = clusters.begin(); itc != clusters.end(); ++itc) {
        c_visited[*itc] = false;
    }

    for (vector<Cluster*>::iterator it = clusters.begin(); it != clusters.end(); ++it) {
        if (c_visited[*it] == false) {
            comp.clear();
            DFSUtil(*it, c_visited);
            tree_component.insert(comp);
        }
    }
    if (ToulBar2::verbose >= 1) {
        cout << "Number of connect components : " << tree_component.size() << endl;
    }
    return tree_component.size();
}

int TreeDecomposition::makeRooted()

{
    bool isalreadyrooted = (roots.size() > 0);
    Cluster* root = NULL;
    list<Cluster*> temproots;
    if (isalreadyrooted) {
        temproots = roots;
    }

    connectedComponents();

    for (unsigned int i = 0; i < clusters.size(); i++) {
        Cluster* c = clusters[i];
        c->accelerateIntersections();
    }

    for (auto it = tree_component.begin(); it != tree_component.end(); ++it) {
        TClusters unvisited(*it);
        //        cout << "clusters component:";
        //        for (auto ccit=unvisited.begin(); ccit != unvisited.end(); ++ccit) cout << ' ' << *ccit;
        //        cout << endl;
        bool selected = false;
        while (unvisited.size() > 0) {
            if (isalreadyrooted) {
                if (temproots.size() > 0) {
                    root = temproots.front();
                    temproots.pop_front();
                } else {
                    // Error, some clusters are missing in the decomposition tree
                    cerr << "Input tree decomposition file is not valid! (may-be cycles within cluster parents)" << endl;
                    throw WrongFileFormat();
                }
            } else {
                if (!selected && ToulBar2::btdRootCluster >= 0 && ToulBar2::btdRootCluster < (int)(*it).size() && unvisited.find(getCluster(ToulBar2::btdRootCluster)) != unvisited.end()) {
                    root = getCluster(ToulBar2::btdRootCluster);
                    selected = true;
                } else {
                    switch (ToulBar2::rootHeuristic) {
                    case 0:
                        root = getBiggerCluster(unvisited);
                        if (ToulBar2::verbose >= 0 && ToulBar2::searchMethod == DFBB && (*it).size() > 1)
                            cout << "Get root cluster C" << root->getId() << " with max. size: " << root->getNbVars() << endl;
                        break;
                    case 1:
                        root = getCluster_height_rootsize_max(unvisited);
                        if (ToulBar2::verbose >= 0 && ToulBar2::searchMethod == DFBB && (*it).size() > 1)
                            cout << "Get root cluster C" << root->getId() << " with max. ratio size/(height-size): " << float(root->getNbVars()) / float(height(root) - root->getNbVars()) << endl;
                        break;
                    case 2:
                        root = getCluster_height_rootsize_min(unvisited);
                        if (ToulBar2::verbose >= 0 && ToulBar2::searchMethod == DFBB && (*it).size() > 1)
                            cout << "Get root cluster C" << root->getId() << " with min. ratio size/(height-size): " << float(root->getNbVars()) / float(height(root) - root->getNbVars()) << endl;
                        break;
                    case 3:
                        root = getClusterMinHeight(unvisited);
                        if (ToulBar2::verbose >= 0 && ToulBar2::searchMethod == DFBB && (*it).size() > 1)
                            cout << "Get root cluster C" << root->getId() << " with min. height: " << height(root) << endl;
                        break;
                    default:
                        cerr << "Unknown root cluster heuristic " << ToulBar2::rootHeuristic << endl;
                        throw BadConfiguration();
                    }
                }
                roots.push_back(root);
                reduceHeight(root, vector<Cluster*>());
                if (ToulBar2::splitClusterMaxSize >= 1)
                    splitClusterRec(root, NULL, ToulBar2::splitClusterMaxSize, unvisited);
                if (ToulBar2::maxSeparatorSize >= 0 || ToulBar2::minProperVarSize >= 2)
                    mergeClusterRec(root, NULL, ToulBar2::maxSeparatorSize, ToulBar2::minProperVarSize, unvisited);
                if (ToulBar2::boostingBTD > 0. && ToulBar2::elimDegree >= 1)
                    boostingVarElimRec(root, NULL, NULL, ToulBar2::elimDegree, unvisited);
                reduceHeight(root, vector<Cluster*>());
            }
            unvisited.erase(root);
            makeRootedRec(root, NULL, unvisited);
            makeDescendants(root);
        }
        assert(unvisited.size() == 0);
    }

    if (ToulBar2::searchMethod != DFBB)
        return 0;

    // if it is a forest then create a unique meta-root cluster with empty separators with its children
    // if it is not a forest but adaptive BTD is used then always create a meta-root cluster
    if (roots.size() > 1 || (ToulBar2::btdMode == 1 && ToulBar2::heuristicFreedom)) {
        root = new Cluster(this);
        root->setId(clusters.size());
        clusters.push_back(root);

        for (list<Cluster*>::iterator iter = roots.begin(); iter != roots.end(); ++iter) {
            Cluster* oneroot = *iter;
            assert(oneroot->getSep() == NULL);

            EnumeratedVariable** scopeVars = new EnumeratedVariable*[1];
            oneroot->setSep(new Separator(wcsp, scopeVars, 0));
            delete[] scopeVars;

            if (oneroot->getNbVars() <= 1 && oneroot->getDescendants().size() == 1) {
                oneroot->getSep()->unqueueSep();
            }
            root->addEdge(oneroot);
            oneroot->setParent(root);
            root->getDescendants().insert(root);

            clusterSum(root->getDescendants(), oneroot->getDescendants());
            sum(root->getVarsTree(), oneroot->getVarsTree()); // the variables of varsTree should also be updated
        }
        roots.clear();
        roots.push_back(root);
    }

    for (unsigned int i = 0; i < clusters.size(); i++) {
        Cluster* c = clusters[i];
        assert(c->getId() == (int)i);
        c->accelerateDescendants();
        if (c->getSep())
            c->getSep()->setSep();
        if (!ToulBar2::approximateCountingBTD) {
            int posx = 0;
            TVars::iterator itv = c->beginVars();
            while (itv != c->endVars()) {
                Variable* var = wcsp->getVar(*itv);
                if (!c->isSepVar(var->wcspIndex))
                    var->setCluster(c->getId());
                else {
                    var->setSep();
                    var->addCluster(c->getId(), posx++); // we add the cluster and also the position of the variable for the delta structure
                }
                ++itv;
            }
        }
        c->setup();
    }
    rootRDS = NULL;
    assert(getRoot() == root);
    root->sortEdgesRec();
    root->sortVarsRec();
    root->sortVarsTreeRec();

    int treewidth = 0;
    for (unsigned int i = 0; i < clusters.size(); i++) {
        Cluster* c = clusters[i];
        if (c->getNbVars() > treewidth)
            treewidth = c->getNbVars();
    }
    if (ToulBar2::verbose >= 0)
        cout << "Tree decomposition width  : " << treewidth - 1 << endl;

    // compute the depth of each cluster
    computeDepths(getRoot(), -1);

    int h = height(root);

    return h;
}

void TreeDecomposition::setDuplicates(bool init)
{
    if (ToulBar2::approximateCountingBTD)
        return;

    static unsigned int curCtr = 0;
    static int curElimBin = 0;
    static int curElimTern = 0;

    if (init) {
        curCtr = 0;
        curElimBin = 0;
        curElimTern = 0;
        // clear all ctrs in clusters
        for (Cluster* c : clusters) {
            c->clearCtrs();
        }
    }

    // assign constraints to clusters and check for duplicate ternary constraints
    for (; curCtr < wcsp->numberOfConstraints(); curCtr++) {
        Constraint* ctr = wcsp->getCtr(curCtr);
        ctr->assignCluster();
    }
    for (; curElimBin < wcsp->elimBinOrder; curElimBin++)
        if (wcsp->elimBinConstrs[curElimBin]->connected()) {
            Constraint* ctr = wcsp->elimBinConstrs[curElimBin];
            ctr->assignCluster();
        }
    for (; curElimTern < wcsp->elimTernOrder; curElimTern++)
        if (wcsp->elimTernConstrs[curElimTern]->connected()) {
            Constraint* ctr = wcsp->elimTernConstrs[curElimTern];
            ctr->assignCluster();
        }

    // check if ternary constraint cluster assignments are valid and do corrections if needed
    for (unsigned int i = 0; i < wcsp->numberOfConstraints(); i++) {
        Constraint* ctr = wcsp->getCtr(i);
        if (ctr->connected() && !ctr->isSep()) {
            if (ctr->isTernary()) {
                TernaryConstraint* tctr = (TernaryConstraint*)ctr;
                tctr->setDuplicates();
            }
        }
    }
    for (int i = 0; i < wcsp->elimTernOrder; i++)
        if (wcsp->elimTernConstrs[i]->connected()) {
            Constraint* ctr = wcsp->elimTernConstrs[i];
            if (ctr->connected() && !ctr->isSep()) {
                assert(ctr->isTernary());
                TernaryConstraint* tctr = (TernaryConstraint*)ctr;
                tctr->setDuplicates();
            }
        }
    // setDuplicates may add binary cost functions
    curElimBin = wcsp->elimBinOrder;
    assert(curElimTern == wcsp->elimTernOrder);
}

void TreeDecomposition::buildFromCovering(string filename)
{
    if (clusters.size() > 0) {
        for (unsigned int i = 0; i < clusters.size(); i++) {
            Cluster* c = clusters[i];
            if (c)
                delete c;
        }
    }
    clusters.clear();
    roots.clear();

    map<int, int> clusterIds;
    int nbclusters = 0;

    ConstraintSet usedctrs;
    //    vector<int> order;

    istringstream sfile(filename.c_str());
    ifstream ffile(filename.c_str(), std::ios::in);
    istream& file = (ToulBar2::bilevel) ? reinterpret_cast<istream&>(sfile) : reinterpret_cast<istream&>(ffile); // use filename as a string containing an explicit covering when doing bilevel optimization
    string fstr;
    while (getline(file, fstr)) {
        std::istringstream file(fstr);
        int num;
        file >> num;
        if (!file)
            break;

        Cluster* C = new Cluster(this);
        clusterIds[num] = nbclusters;
        C->setId(nbclusters);
        clusters.push_back(C);
        nbclusters++;

        int num_parent;
        file >> num_parent;
        assert((num_parent == -1) || (clusterIds.find(num_parent) != clusterIds.end()));

        int v = -1;
        while (file >> v) {
            if (!C->isVar(v)) {
                C->addVar(wcsp->getVar(v));
                //                if ((num_parent == -1) || (!clusters[clusterIds[num_parent]]->isVar(v))) {
                //                  order.push_back(v);
                //                }
            }
        }

        if (num_parent >= 0) {
            C->addEdge(clusters[clusterIds[num_parent]]);
            clusters[clusterIds[num_parent]]->addEdge(C);
        } else {
            roots.push_back(C);
        }

        for (TVars::iterator iter = C->getVars().begin(); iter != C->getVars().end(); iter++) {
            ConstraintList* xctrs = wcsp->getVar(*iter)->getConstrs();
            for (ConstraintList::iterator it = xctrs->begin(); it != xctrs->end(); ++it) {
                Constraint* ctr = (*it).constr;
                bool used = usedctrs.find(ctr) != usedctrs.end();
                if (!used) {
                    int k = 0;
                    while ((k < ctr->arity()) && (C->isVar(ctr->getVar(k)->wcspIndex)))
                        k++;

                    if (k == ctr->arity()) {
                        usedctrs.insert(ctr);
                        C->addCtr(ctr);
                    }
                }
            }
        }
    }
    if (!ToulBar2::bilevel) {
        ffile.close();
    }
    //    reverse(order.begin(), order.end()); // must return an elimination order, the reverse of a topological order

    // buildFromOrderNext(order);
    // if (ToulBar2::btdMode == 3) pathFusions(order); // bug: it assumes cluster.size() == order.size()
    // else treeFusions(); // we assume there is no cluster separator included into another cluster separator in the input tree decomposition

    for (unsigned int i = 0; i < clusters.size(); i++) {
        Cluster* c = clusters[i];
        if (c)
            c->getDescendants().clear();
    }

    int h = makeRooted();

    if (ToulBar2::verbose >= 0)
        cout << "Tree decomposition height : " << h << endl;
    setDuplicates(true);
    if (ToulBar2::verbose >= 0)
        cout << "Number of clusters         : " << clusters.size() << endl;
    if (ToulBar2::debug >= 1 || ToulBar2::verbose >= 1)
        print();
    if (ToulBar2::dumpWCSP)
        dump();
    assert(verify());
}

void TreeDecomposition::buildFromOrder()
{
    vector<int> order;
    assert(!((WCSP*)wcsp)->isAlreadyTreeDec(ToulBar2::varOrder));
    ((WCSP*)wcsp)->elimOrderFile2Vector(ToulBar2::varOrder, order);
    if (!ToulBar2::varOrder)
        reverse(order.begin(), order.end());

    if (clusters.size() > 0) {
        for (unsigned int i = 0; i < clusters.size(); i++) {
            Cluster* c = clusters[i];
            if (c)
                delete c;
        }
    }
    clusters.clear();

    for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
        Cluster* c = new Cluster(this);
        c->setId(i);
        c->addVar(wcsp->getVar(order[i]));
        clusters.push_back(c);
    }
    ConstraintSet usedctrs;

    for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
        Variable* x = wcsp->getVar(order[i]);
        Cluster* c = clusters[i];

        ConstraintList* xctrs = x->getConstrs();
        for (ConstraintList::iterator it = xctrs->begin(); it != xctrs->end(); ++it) {
            Constraint* ctr = (*it).constr;
            bool used = usedctrs.find(ctr) != usedctrs.end();
            if (!used) {
                usedctrs.insert(ctr);
                c->addCtr(ctr);
                for (int k = 0; k < ctr->arity(); k++)
                    if (ctr->getVar(k)->unassigned())
                        c->addVar(ctr->getVar(k));
            }
        }

        for (unsigned int j = i + 1; j < wcsp->numberOfVariables(); j++) {
            if (c->isVar(order[j])) {
                Cluster* cj = clusters[j];
                TVars::iterator it = c->beginVars();
                while (it != c->endVars()) {
                    cj->addVar(wcsp->getVar(*it));
                    ++it;
                }
                cj->removeVar(x);
                c->addEdge(cj);
                cj->addEdge(c);
                break;
            }
        }
    }
    buildFromOrderNext(order);
}

void TreeDecomposition::buildFromOrderForApprox()
{

    vector<int> order;
    bool firstComponent = true;
    int sizepart = 0; // number of parts in the built partition
    ConstraintSet totalusedctrs; // constraints already in a part
    vector<int> degreeinusedctr; // number of constraints not adding for each variable
    //	int nbcstr = 0;					//
    double time;

    assert(!((WCSP*)wcsp)->isAlreadyTreeDec(ToulBar2::varOrder));
    ((WCSP*)wcsp)->elimOrderFile2Vector(ToulBar2::varOrder, order);
    if (!ToulBar2::varOrder)
        reverse(order.begin(), order.end());

    if (clusters.size() > 0) {
        for (unsigned int i = 0; i < clusters.size(); i++) {
            Cluster* c = clusters[i];
            if (c)
                delete c;
        }
    }
    clusters.clear();
    for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
        Variable* x = wcsp->getVar(i);
        //		degree.push_back(x->getTrueDegree());
        degreeinusedctr.push_back(x->getDegree());
    }
    time = cpuTime();
    while (totalusedctrs.size() < wcsp->numberOfConnectedConstraints()) //&& nbparties<4)
    {
        ConstraintSet currentusedctrs; // liste des contraintes contenues dans la partie courante
        TVars currentusedvars; // liste des variables contenues dans la partie courante
        TVars inusedvars; // liste des variables qui n'ont pas encore ete etudiee dans la partie courante
        vector<Variable*> currentRevElimOrder; // liste des variables dans l'ordre inverse construit
        sizepart++;
        for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
            if (wcsp->unassigned(i)) {
                if (wcsp->getDegree(i) == 0) {
                    if (firstComponent) {
                        currentRevElimOrder.push_back(wcsp->getVar(i));
                        currentusedvars.insert(i);
                    }
                } else {
                    if (degreeinusedctr[i] > 0)
                        inusedvars.insert(i);
                }
            }
        }

        maxchord(sizepart, order, totalusedctrs, inusedvars, currentusedvars, currentRevElimOrder, currentusedctrs);

        // insert into tree decomposition

        // supprime les variables qui n'ont aucune contraintes dans la partition
        for (vector<Variable*>::iterator it = currentRevElimOrder.begin(); it != currentRevElimOrder.end();) {
            if (currentusedvars.find((*it)->wcspIndex) == currentusedvars.end()) {
                it = currentRevElimOrder.erase(it);
            } else
                ++it;
        }

        if (sizepart == 1)
            cout << endl;
        cout << "part " << sizepart << " : " << currentRevElimOrder.size() << " variables and " << currentusedctrs.size() << " constraints (really added)\n";
        if (ToulBar2::debug >= 1 || ToulBar2::verbose >= 3) // affichage
        {
            cout << "\tVariables : ";
            for (vector<Variable*>::iterator it = currentRevElimOrder.begin(); it != currentRevElimOrder.end(); it++) {
                cout << (*it)->wcspIndex << " ";
            }
            cout << endl;
            cout << "\tContraintes : ";
            for (ConstraintSet::iterator it = currentusedctrs.begin(); it != currentusedctrs.end(); it++) {
                cout << "[";
                for (int k = 0; k < (*it)->arity(); k++) {
                    cout << (*it)->getVar(k)->wcspIndex;
                    if (k != (*it)->arity() - 1)
                        cout << " ";
                }
                cout << "] ";
            }
            cout << endl;
        }

        insert(sizepart, currentRevElimOrder, currentusedctrs);
        firstComponent = false;
    }
    time = cpuTime() - time;
    cout << "--> number of parts : " << sizepart << endl;
    cout << "--> time : " << time << " seconds. " << endl
         << endl;
    buildFromOrderNext(order);
}

void TreeDecomposition::buildFromOrderNext(vector<int>& order)
{

    if (ToulBar2::verbose >= 2) {
        cout << "----- Before fusions process: " << endl;
        for (unsigned int i = 0; i < clusters.size(); i++) {
            if (!clusters[i])
                continue;
            Cluster* c = clusters[i];
            c->print();
        }
        cout << "----- fusions process starting... " << endl;
    }

    if (ToulBar2::btdMode == 3)
        pathFusions(order);
    else
        treeFusions();

    if (ToulBar2::verbose >= 2) {
        cout << "----- After fusions process: " << endl;
        for (unsigned int i = 0; i < clusters.size(); i++) {
            if (!clusters[i])
                continue;
            Cluster* c = clusters[i];
            c->print();
        }
        cout << "----- fusions process ended... " << endl;
    }

    for (unsigned int i = 0; i < clusters.size(); i++) {
        Cluster* c = clusters[i];
        c->getDescendants().clear();
    }

    roots.clear();
    int h = makeRooted();
    if (ToulBar2::searchMethod != DFBB)
        return;
    if (ToulBar2::verbose >= 0)
        cout << "Tree decomposition height : " << h << endl;
    setDuplicates(true);
    if (ToulBar2::verbose >= 0)
        cout << "Number of clusters        : " << clusters.size() << endl;
    if (ToulBar2::debug >= 1 || ToulBar2::verbose >= 1)
        print();
    if (ToulBar2::dumpWCSP)
        dump();
    assert(verify());
}

void TreeDecomposition::maxchord(int sizepart, vector<int>& order, ConstraintSet& totalusedctrs, TVars& inusedvars, TVars& currentusedvars, vector<Variable*>& currentRevElimOrder, ConstraintSet& currentusedctrs)
{
    vector<TVars> listeVars(wcsp->numberOfVariables()); // liste des voisins d'ordre superieur de chaque variable
    int nbcstr = 0;
    double time, timetot = 0;
    while (inusedvars.size() > 0) {
        int maxsize = -1;
        Variable* maxvar = NULL; /* next variable */

        // Choose the nex variable
        for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
            Variable* x = wcsp->getVar(order[i]);
            if (inusedvars.find(x->wcspIndex) != inusedvars.end()) {
                int size = listeVars[x->wcspIndex].size();
                if (size > maxsize) {
                    maxsize = size;
                    maxvar = x;
                }
            }
        }

        if (maxvar) {
            //				cout << "Variable choisie: " << maxvar->wcspIndex << " ";
            ConstraintList* xctrs = maxvar->getConstrs();
            for (ConstraintList::iterator it = xctrs->begin(); it != xctrs->end(); ++it) {
                Constraint* ctr = (*it).constr;
                bool used = totalusedctrs.find(ctr) != totalusedctrs.end();
                if (!used) {
                    TVars scopectr;
                    TVars sc;
                    for (int k = 0; k < ctr->arity(); k++) {
                        Variable* x = ctr->getVar(k);
                        if (x->wcspIndex != maxvar->wcspIndex && wcsp->unassigned(x->wcspIndex)) {
                            sc.insert(x->wcspIndex);
                            if (inusedvars.find(x->wcspIndex) != inusedvars.end())
                                scopectr.insert(x->wcspIndex);
                        }
                    }
                    if (scopectr.size() == 0) { // all edges of the ctr are in the sub graph => the cstr is added in this current part
                        if (included(sc, listeVars[maxvar->wcspIndex])) {
                            ConstraintSet subctr;
                            nbcstr++;
                            currentusedctrs.insert(ctr);
                            totalusedctrs.insert(ctr);
                            time = cpuTime();
                            subctr = ctr->subConstraint();
                            ctrSum(totalusedctrs, subctr);
                            ctrSum(currentusedctrs, subctr);
                            time = time - cpuTime();
                            timetot += time;
                            sum(currentusedvars, sc);
                            currentusedvars.insert(maxvar->wcspIndex);
                        }
                    }

                    for (TVars::iterator i = scopectr.begin(); i != scopectr.end(); ++i) {
                        int vari = wcsp->getVar(*i)->wcspIndex;
                        int varj = maxvar->wcspIndex;
                        if (included(listeVars[vari], listeVars[varj])) {
                            listeVars[(*i)].insert(varj);
                            //--degree[(*i)];
                        }
                    }
                }
            }
            currentRevElimOrder.push_back(maxvar);
            inusedvars.erase(maxvar->wcspIndex);
        }
    }
}

void TreeDecomposition::insert(int sizepart, vector<Variable*> currentRevElimOrder, ConstraintSet currentusedctrs)
{
    int firstCluster = clusters.size();
    for (unsigned int i = 0; i < currentRevElimOrder.size(); i++) {
        Cluster* c = new Cluster(this);
        c->setId(clusters.size());
        c->addVar(currentRevElimOrder[currentRevElimOrder.size() - i - 1]);
        clusters.push_back(c);
    }
    ConstraintSet usedctrs;

    for (unsigned int i = 0; i < currentRevElimOrder.size(); i++) {
        Cluster* c = clusters[firstCluster + i];
        c->setPart(sizepart);
        Variable* x = currentRevElimOrder[currentRevElimOrder.size() - i - 1];

        ConstraintList* xctrs = x->getConstrs();
        for (ConstraintList::iterator it = xctrs->begin(); it != xctrs->end(); ++it) {
            Constraint* ctr = (*it).constr;
            bool used = usedctrs.find(ctr) != usedctrs.end();
            if (!used) {
                if (currentusedctrs.find(ctr) != currentusedctrs.end()) {
                    usedctrs.insert(ctr);
                    c->addCtr(ctr);
                    for (int k = 0; k < ctr->arity(); k++) {
                        if (ctr->getVar(k)->unassigned()) {
                            // assert(currentusedvars.find(ctr->getVar(k)->wcspIndex) != currentusedvars.end());
                            c->addVar(ctr->getVar(k));
                        }
                    }
                }
            }
        }

        for (unsigned int j = i + 1; j < currentRevElimOrder.size(); j++) {
            if (c->isVar(currentRevElimOrder[currentRevElimOrder.size() - j - 1]->wcspIndex)) {
                Cluster* cj = clusters[firstCluster + j];
                TVars::iterator it = c->beginVars();
                while (it != c->endVars()) {
                    cj->addVar(wcsp->getVar(*it));
                    ++it;
                }
                cj->removeVar(x);
                c->addEdge(cj);
                cj->addEdge(c);
                break;
            }
        }
    }
}

void TreeDecomposition::getElimVarOrder(vector<int>& elimVarOrder)
{
    getRoot()->getElimVarOrder(elimVarOrder);
}

void TreeDecomposition::addDelta(int cyid, EnumeratedVariable* x, Value value, Cost cost)
{
    Cluster* cy = getCluster(cyid);
    Cluster* cx = getCluster(x->getCluster());
    if (!cy->isDescendant(cx) && !isSameCluster(cy, cx)) {
        int ckid, posx;
        assert(x->clusters.size() > 0);
        if (cost != MIN_COST && !deltaModified[x->wcspIndex])
            deltaModified[x->wcspIndex] = true;
        x->beginCluster();
        while (x->nextCluster(ckid, posx)) {
            Cluster* ck = getCluster(ckid);
            if (ck->isDescendant(cy) || isSameCluster(ck, cy)) {
                if (ToulBar2::verbose >= 2)
                    cout << "add delta " << cost << " to var " << x->wcspIndex << " (cluster " << cx->getId() << ") value " << value << " from subtree " << ck->getId() << " (cluster " << cyid << ")" << endl;
                ck->addDelta(posx, value, cost);
            }
        }
    }
}

// warning! variables are not assigned to the current new solution
// use assignment "a" instead
void TreeDecomposition::newSolution(Cost lb)
{
    ToulBar2::deltaUb = max(ToulBar2::deltaUbAbsolute, (Cost)(ToulBar2::deltaUbRelativeGap * (Double)lb));
    wcsp->setUb(lb);

    TAssign a;

    Cluster* root = getRoot();
    wcsp->restoreSolution(root);
    root->getSolution(a);

    if ((ToulBar2::elimDegree > 0 || ToulBar2::elimDegree_preprocessing > 0 || ToulBar2::preprocessFunctional > 0) && root->getNbVars() == 0) {
        // recorded solutions in clusters containing a single variable eliminated in preprocessing may be wrong due to variable elimination in preprocessing; must be recovered after complete assignment and restoreSolution
        for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
            if (wcsp->enumerated(i)) {
                EnumeratedVariable* x = (EnumeratedVariable*)wcsp->getVar(i);
                x->assignWhenEliminated(a[i]);
            }
        }
        wcsp->restoreSolution();
        for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
            if (wcsp->enumerated(i)) {
                a[i] = wcsp->getValue(i);
            }
        }
    }
    if (!ToulBar2::isZ)
        wcsp->setSolution(lb, &a);

    if (ToulBar2::showSolutions) {
        wcsp->printSolution();
        cout << endl;
    }

    if (!ToulBar2::uaieval && ToulBar2::writeSolution && ToulBar2::solutionFile != NULL) {
        if (!ToulBar2::allSolutions)
            fseek(ToulBar2::solutionFile, ToulBar2::solutionFileRewindPos, SEEK_SET);
        wcsp->printSolution(ToulBar2::solutionFile);
        fprintf(ToulBar2::solutionFile, "\n");
    }

    if (ToulBar2::xmlflag) {
        cout << "o " << std::fixed << std::setprecision(0) << wcsp->Cost2ADCost(lb) << std::setprecision(DECIMAL_POINT) << endl; //" ";
        ((WCSP*)wcsp)->solution_XML(false);
    }
    if (ToulBar2::maxsateval) {
        cout << "o " << lb << endl;
    }
    if (ToulBar2::uaieval && !ToulBar2::isZ) {
        wcsp->solution_UAI(lb);
    }
    // warning: cannot read solution from variable assignments
    // else if(ToulBar2::pedigree){
    // 	ToulBar2::pedigree->printSol(wcsp);
    // }
    // else if(ToulBar2::haplotype){
    //   ToulBar2::haplotype->printSol(wcsp);
    // }

    if (ToulBar2::newsolution)
        (*ToulBar2::newsolution)(wcsp->getIndex(), wcsp->getSolver());
}

void TreeDecomposition::intersection(TVars& v1, TVars& v2, TVars& vout)
{
    assert(&vout != &v1);
    assert(&vout != &v2);
    set_intersection(v1.begin(), v1.end(),
        v2.begin(), v2.end(),
        inserter(vout, vout.begin()));
}

// returns precomputed intersection if available or compute it
void TreeDecomposition::intersection(Cluster* c, Cluster* cj, TVars& vout)
{
    c->quickIntersection(cj, vout);
#ifndef NDEBUG
    TVars vcompare;
    intersection(c->getVars(), cj->getVars(), vcompare);
    assert(vcompare == vout);
#endif
}

void TreeDecomposition::difference(TVars& v1, TVars& v2, TVars& vout)
{
    assert(&vout != &v1);
    assert(&vout != &v2);
    set_difference(v1.begin(), v1.end(),
        v2.begin(), v2.end(),
        inserter(vout, vout.begin()));
}

void TreeDecomposition::sum(TVars v1, TVars v2, TVars& vout)
{
    set_union(v1.begin(), v1.end(),
        v2.begin(), v2.end(),
        inserter(vout, vout.begin()));
}

void TreeDecomposition::sum(TVars& v1, TVars& v2)
{
    v1.insert(v2.begin(), v2.end());
}

bool TreeDecomposition::included(TVars& v1, TVars& v2)
{
    TVars vout;
    set_union(v1.begin(), v1.end(),
        v2.begin(), v2.end(),
        inserter(vout, vout.begin()));
    return vout.size() == v2.size();
}

void TreeDecomposition::clusterSum(TClusters v1, TClusters v2, TClusters& vout)
{
    set_union(v1.begin(), v1.end(),
        v2.begin(), v2.end(),
        inserter(vout, vout.begin()));
}

void TreeDecomposition::clusterSum(TClusters& v1, TClusters& v2)
{
    v1.insert(v2.begin(), v2.end());
}

void TreeDecomposition::ctrSum(TCtrs v1, TCtrs v2, TCtrs& vout)
{
    set_union(v1.begin(), v1.end(),
        v2.begin(), v2.end(),
        inserter(vout, vout.begin()));
}

void TreeDecomposition::ctrSum(TCtrs& v1, TCtrs& v2)
{
    v1.insert(v2.begin(), v2.end());
}

bool TreeDecomposition::verify()
{
    if (!ToulBar2::approximateCountingBTD) {
        for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
            Variable* x = wcsp->getVar(i);
            if (x->assigned())
                continue;

            Cluster* ci = clusters[x->getCluster()];
            if (!ci->isVar(x->wcspIndex) || ci->isSepVar(x->wcspIndex)) {
                cout << "cluster: " << ci->getId() << " , var " << x->wcspIndex << endl;
                return false;
            }
            //  	    ConstraintList* xctrs = x->getConstrs();
            //  	    for (ConstraintList::iterator it=xctrs->begin(); it != xctrs->end(); ++it) {
            //              Constraint* ctr = (*it).constr;
            //  			Cluster* cj  = clusters[ctr->getCluster()];
            //              int arity = ctr->arity();
            //              for(i=0;i<arity;i++) {
            //          		Variable* x = ctr->getVar(i);

            //              }
            //  	    }
        }
    }
    return true;
}

void TreeDecomposition::printStats(Cluster* c)
{
    if (!c)
        return;
    c->printStats();
    TClusters::iterator itj = c->beginEdges();
    while (itj != c->endEdges()) {
        Cluster* cj = *itj;
        ++itj;
        printStats(cj);
    }
}

void TreeDecomposition::print(Cluster* c, int recnum)
{
    if (!c) {
        //  		for(unsigned int i=0;i<wcsp->numberOfVariables();i++) {
        //  			Variable* x = wcsp->getVar(i);
        //  			x->beginCluster();
        //  			int c,posx;
        //  			cout << x->wcspIndex << " appears in sep {";
        //  			while(x->nextCluster(c,posx)) {
        //  				cout << c << " ";
        //  			}
        //  			cout << "}" << endl;
        //  		}
        if (roots.empty())
            return;
        c = *roots.begin();
    }

    for (int i = 0; i < recnum; i++)
        cout << "  ";
    c->print();

    TClusters::iterator ita = c->beginSortedEdges();
    while (ita != c->endSortedEdges()) {
        print(*ita, recnum + 1);
        ++ita;
    }
}

void TreeDecomposition::dump(Cluster* c)
{
    if (!c) {
        char tmpName[256];
        sprintf(tmpName, "%s.info", getWCSP()->getName().c_str());
#ifdef __WIN32__
        mkdir(tmpName);
#else
        mkdir(tmpName, 0777);
#endif

        sprintf(tmpName, "%s.info/root", getWCSP()->getName().c_str());

        ofstream rootFile(tmpName);
        if (roots.empty()) {
            rootFile.close();
            return;
        }
        c = *roots.begin();
        rootFile << c->getId();
        rootFile.close();
    }

    c->dump();

    TClusters::iterator ita = c->beginEdges();
    while (ita != c->endEdges()) {
        dump(*ita);
        ++ita;
    }
}

// deconnect cluster subtree rooted at c and the associated separators
void TreeDecomposition::updateInTD(Cluster* c)
{
    if (c->getFreedom()) {
        // all the descendants clusters are discarded from the TD
        for (TClusters::iterator iter = c->beginDescendants(); iter != c->endDescendants(); ++iter)
            if ((*iter)->getId() != c->getId()) {
                (*iter)->setIsCurrInTD(false);
                (*iter)->getSep()->deconnect();
            }
    }
}

Cluster* TreeDecomposition::lowestCommonAncestor(Cluster* c1, Cluster* c2)
{
    if (c1->getDepth() < c2->getDepth()) {
        while (!c1->isDescendant(c2)) {
            c1 = c1->getParent();
            assert(c1 != NULL);
        }
        return c1;
    } else {
        while (!c2->isDescendant(c1)) {
            c2 = c2->getParent();
            assert(c2 != NULL);
        }
        return c2;
    }
}

bool TreeDecomposition::isSameCluster(Cluster* c1, Cluster* c2)
{
    if (c1 == c2)
        return true;
    if (ToulBar2::heuristicFreedom) {
        if (c1->getIsCurrInTD() && c2->getIsCurrInTD())
            return false;
        if ((c1->getFreedom() || !c1->getIsCurrInTD()) && c1->isDescendant(c2))
            return true;
        if ((c2->getFreedom() || !c2->getIsCurrInTD()) && c2->isDescendant(c1))
            return true;
        Cluster* lca = lowestCommonAncestor(c1, c2);
        assert(getRoot() != lca || (!lca->getFreedom() && lca->getIsCurrInTD()));
        if (lca->getFreedom() || !lca->getIsCurrInTD())
            return true;
    }
    return false;
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

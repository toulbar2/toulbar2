
#include "tb2naryconstr.hpp"
#include "search/tb2clusters.hpp"

NaryConstraint::NaryConstraint(WCSP* wcsp, EnumeratedVariable** scope_in, int arity_in, Cost defval, Long nbtuples)
    : AbstractNaryConstraint(wcsp, scope_in, arity_in)
    , pf(NULL)
    , costs(NULL)
    , costSize(0)
    , default_cost(defval)
    , filters(NULL)
{
    for (int i = 0; i < arity_in; i++) {
        conflictWeights.push_back(0);
        unsigned int domsize = scope_in[i]->getDomainInitSize();
        if (domsize > (unsigned int)std::numeric_limits<tValue>::max()) {
            cerr << "Nary constraints overflow. Extend tValue type range." << endl;
            throw BadConfiguration();
        }
    }

    Cost Top = wcsp->getUb();
    if (default_cost > Top)
        default_cost = Top;

    pf = new TUPLES;
    if (nbtuples != 0 && expandtodo(min(nbtuples, getDomainSizeProduct())))
        expand();

    // Cannot propagate here because cost tuples are not known yet
}

NaryConstraint::NaryConstraint(WCSP* wcsp)
    : AbstractNaryConstraint(wcsp)
    , pf(NULL)
    , costs(NULL)
    , costSize(0)
    , default_cost(wcsp->getUb())
    , filters(NULL)
{
    pf = new TUPLES;
}

NaryConstraint::~NaryConstraint()
{
    if (pf)
        delete pf;
    if (costs)
        delete[] costs;
    if (filters)
        delete filters;
}

// USED ONLY DURING SEARCH
void NaryConstraint::assign(int varIndex)
{
    if (connected(varIndex)) {
        deconnect(varIndex);

        if (size() <= NARYDECONNECTSIZE && universal()) { // check if it is satisfied (clause)
            //	  cout << "nary cost function satisfied: " << this << endl;
            deconnect();
            return;
        }

        if (getNonAssigned() <= NARYPROJECTIONSIZE && (getNonAssigned() < 3 || maxInitDomSize <= NARYPROJECTION3MAXDOMSIZE || prodInitDomSize <= NARYPROJECTION3PRODDOMSIZE)) {
            //	  cout << "Assign var " << *getVar(varIndex) << "  in  " << *this;
            deconnect();
            projectNary();
        } else {
            if (ToulBar2::FullEAC)
                reviseEACGreedySolution();
        }
    }
}

// NOT USED
// void NaryConstraint::projectFromZero(int index)
//{
//    int i;
//    int nsup = 0;
//    Cost minc = MAX_COST;
//    Cost c;
//    EnumeratedVariable* var = (EnumeratedVariable*) getVar(index);
//    for (EnumeratedVariable::iterator iter = var->begin(); iter != var->end(); ++iter) {
//        c = var->getCost(*iter);
//        if(c == MIN_COST) { if(nsup) return; nsup++; }
//        else if(c < minc) minc = c;
//
//    }
//    for(i=0;i<arity_;i++) {
//        if(i != index) {
//            var = (EnumeratedVariable*) getVar(i);
//            if(getVar(i)->unassigned()) {
//                nsup = 0;
//                for (EnumeratedVariable::iterator iter = var->begin(); iter != var->end(); ++iter) {
//                    c = var->getCost(*iter);
//                    if(c == MIN_COST) { if(nsup) return; nsup++; }
//                    else if(c < minc) minc = c;
//                }
//            }
//        }
//    }
//    int a = arity_;
//    Tuple cht(a);
//    for(i=0;i<a;i++) {
//        var = (EnumeratedVariable*) getVar(i);
//        cht[i] = var->toIndex(var->getSup());
//    }
//    c = eval(cht);
//    if(c < minc) minc = c;
//    if(c > MIN_COST) starrule(cht,minc);
//}

// NOT USED
// void NaryConstraint::starrule(const Tuple& t, Cost minc)
//{
//    int i;
//    for(i=0;i<arity_;i++) {
//        EnumeratedVariable* var = (EnumeratedVariable*) getVar(i);
//        if(getVar(i)->unassigned()) {
//            for (EnumeratedVariable::iterator iter = var->begin(); iter != var->end(); ++iter) {
//                Cost c = var->getCost(*iter);
//                if(c > MIN_COST) var->extend(*iter, minc);
//            }
//        }
//
//    }
//    setTuple( t, eval(t) - minc );
//    projectLB(minc);
//}

void NaryConstraint::expand()
{
    if (pf == NULL)
        return;
    assert(costs == NULL);

    ptrdiff_t sz = getDomainInitSizeProduct();
    if (ToulBar2::elimSpaceMaxMB && (Double)sz * sizeof(Cost) > (Double)ToulBar2::elimSpaceMaxMB * 1024. * 1024.)
        return;

    try {
        costs = new Cost[sz];
    } catch (const std::bad_alloc&) {
        if (ToulBar2::verbose >= 1)
            cout << "Warning! nary expand cannot be done! " << sz << endl;
        costs = NULL;
        return;
    }
    costSize = sz;
    std::fill(costs, costs + sz, default_cost);
    for (TUPLES::iterator it = pf->begin(); it != pf->end(); ++it) {
        costs[getCostsIndex(it->first)] = it->second;
    }
    delete pf;
    pf = NULL;
}

void NaryConstraint::resetFilters()
{
    if (filters) {
        delete filters;
        filters = NULL;
    }
}

void NaryConstraint::fillFilters()
{
    if (filters == NULL) {
        filters = new ConstraintSet;
        for (int i = 0; i < arity_; i++) {
            EnumeratedVariable* v = (EnumeratedVariable*)getVar(i);
            for (ConstraintList::iterator iter = v->getConstrs()->begin(); iter != v->getConstrs()->end(); ++iter) {
                Constraint* ctr = (*iter).constr;
                if (scopeIncluded(ctr))
                    filters->insert(ctr);
            }
        }
    }
}

Cost NaryConstraint::eval(const Tuple& tin, EnumeratedVariable** scope_in)
{
    if (scope_in) {
        for (int i = 0; i < arity_; i++) {
            int pos = getIndex(scope_in[i]);
            evalTuple[pos] = tin[i];
        }
        return eval(evalTuple);
    } else
        return eval(tin);
}

Cost NaryConstraint::eval(const Tuple& s)
{
    if (pf) {
        TUPLES& f = *pf;
        TUPLES::iterator it = f.find(s);
        if (it != f.end())
            return it->second;
        else
            return default_cost;
    } else
        return costs[getCostsIndex(s)];
}

bool NaryConstraint::checkEACGreedySolution(int index, Value support)
{
    int a = arity();
    for (int i = 0; i < a; i++) {
        EnumeratedVariable* var = (EnumeratedVariable*)getVar(i);
        evalTuple[i] = var->toIndex((i == index) ? support : var->getSupport());
    }
    return (eval(evalTuple) == MIN_COST);
}

bool NaryConstraint::reviseEACGreedySolution(int index, Value support)
{
    bool result = checkEACGreedySolution(index, support);
    if (!result) {
        if (index >= 0) {
            getVar(index)->unsetFullEAC();
        } else {
            int a = arity();
            for (int i = 0; i < a; i++) {
                getVar(i)->unsetFullEAC();
            }
        }
    }
    return result;
}

/// Set new default cost to df (df <= Top), keep existing costs SMALLER than this default cost in a Map or a table
void NaryConstraint::keepAllowedTuples(Cost df)
{
    static Tuple t;
    assert(CUT(wcsp->getUb(), df));
    if (pf) {
        TUPLES* pfnew = new TUPLES;

        Cost c;
        firstlex();
        while (nextlex(t, c)) {
            if (c < df)
                (*pfnew)[t] = c;
        }

        delete pf;
        pf = pfnew;
    } else {
        int a = arity_;
        vector<unsigned int> t(a, 0);
        for (ptrdiff_t idx = 0; idx < costSize; idx++) {
            bool ok = true;
            for (int i = 0; i < a && ok; i++) {
                EnumeratedVariable* var = (EnumeratedVariable*)getVar(i);
                ok = ok && var->canbe(var->toValue(t[i]));
            }
            if (!ok || costs[idx] > df)
                costs[idx] = df;
            int i = a - 1;
            while (i >= 0 && t[i] == ((EnumeratedVariable*)getVar(i))->getDomainInitSize() - 1) {
                t[i] = 0;
                i--;
            }
            if (i >= 0)
                t[i]++;
        }
    }
    default_cost = df;
}

bool NaryConstraint::consistent(const Tuple& t)
{
    int a = arity_;
    bool ok = true;
    for (int i = 0; i < a && ok; i++) {
        EnumeratedVariable* var = (EnumeratedVariable*)getVar(i);
        ok = ok && var->canbe(var->toValue(t[i]));
    }
    return ok;
}

void NaryConstraint::first()
{
    if (pf)
        tuple_it = pf->begin();
    else
        firstlex();
}

bool NaryConstraint::next(Tuple& t, Cost& c)
{
    bool ok = false;
    if (pf) {
        while (!ok && (tuple_it != pf->end())) {
            t = tuple_it->first;
            c = tuple_it->second;
            ok = (c != default_cost) && consistent(t);
            tuple_it++;
        }
    } else {
        while (!ok && nextlex(t, c)) {
            ok = (c != default_cost);
        }
    }
    return ok;
}

void NaryConstraint::first(EnumeratedVariable* vx, EnumeratedVariable* vz)
{
    it_values.clear();
    EnumeratedVariable* var;
    it_values.push_back(vx->begin());
    for (int i = 0; i < arity_; i++) {
        var = (EnumeratedVariable*)getVar(i);
        if (var != vx && var != vz)
            it_values.push_back(var->begin());
    }
    it_values.push_back(vz->begin());
}

bool NaryConstraint::separability(EnumeratedVariable* vx, EnumeratedVariable* vz)
{
    static Tuple t1, t;
    int a = arity_, i = a - 1, k;
    bool finished = false;
    bool neweq = true;
    bool sev = true; // false if vx and vz are not separable
    Cost diff = 0;
    Cost c1, c;
    vector<int> ordre(a, -1);
    EnumeratedVariable* var;
    EnumeratedVariable* scope_in[a];
    t1.resize(a);
    t.resize(a);

    EnumeratedVariable::iterator itvzfirst = vz->begin();
    EnumeratedVariable::iterator itvznext = itvzfirst;
    if (itvznext != vz->end())
        ++itvznext;
    first(vx, vz);

    scope_in[0] = vx; //.push_back( alpha );
    ordre[a - 1] = 0;
    k = 1;
    if (ToulBar2::verbose >= 1)
        cout << "[ " << vx->getName() << " ";
    for (int j = 0; j < a; j++) {
        var = (EnumeratedVariable*)getVar(j);
        if (var != vx && var != vz) {
            scope_in[k] = var;
            ordre[k - 1] = k;
            if (ToulBar2::verbose >= 1)
                cout << var->getName() << " ";
            ++k;
        }
    }

    scope_in[a - 1] = vz; //.push_back( beta );
    ordre[a - 2] = a - 1;
    // if(it_values[ordre[a-2]] != vz->end()) ++it_values[ordre[a-2]];
    if (ToulBar2::verbose >= 1)
        cout << vz->getName() << " ] ? "; // << endl;

    while (sev && itvznext != vz->end()) {
        first(vx, vz);
        while (sev && it_values[ordre[0]] != scope_in[ordre[0]]->end()) {
            if (it_values[ordre[a - 2]] == itvzfirst /*vz->begin()*/)
                ++it_values[ordre[a - 2]];
            // if(it_values[ordre[a-2]] != vz->end()){
            for (int j = 0; j < a; j++) {
                t[j] = scope_in[j]->toIndex(*it_values[j]); // PROBLEME CAR it_values[ordre[a-2]] = vz->end()
                t1[j] = scope_in[j]->toIndex(*it_values[j]);
            }
            t1[a - 1] = scope_in[a - 1]->toIndex(*itvzfirst /*vz->begin()*/);

            // etude de la difference
            c1 = eval(t1, scope_in);
            c = eval(t, scope_in);
            if (!universe(c, c1, wcsp->getUb())) {

                if (ToulBar2::verbose >= 3) {
                    if (!neweq)
                        cout << "= \n";
                    else
                        cout << endl;
                }
                if (ToulBar2::verbose >= 3) {
                    cout << "C";
                    for (unsigned int j = 0; j < t.size(); j++) {
                        cout << t[j] << ".";
                    }
                    cout << " - C";
                    for (unsigned int j = 0; j < t1.size(); j++) {
                        cout << t1[j] << ".";
                    }
                    cout << " = " << c << " - " << c1;
                }
                if (neweq) {
                    diff = squareminus(c, c1, wcsp->getUb());
                    neweq = false;
                } else
                    sev = (diff == squareminus(c, c1, wcsp->getUb()));
                if (ToulBar2::verbose >= 3)
                    cout << " = " << squareminus(c, c1, wcsp->getUb()) << endl;
            } else if (ToulBar2::verbose >= 3)
                cout << "universe\n";
            finished = false;
            i = a - 1;
            while (sev && !finished) {
                var = scope_in[ordre[i]];
                ++it_values[ordre[i]];
                finished = it_values[ordre[i]] != var->end();
                if (!finished) {
                    // it_values[ordre[i]] = var->begin();
                    if (i > 0) {
                        it_values[ordre[i]] = var->begin();
                        if (i == a - 1)
                            neweq = true;
                        i--;
                    } else
                        finished = true;
                }
            }
            if (it_values[ordre[a - 2]] == itvzfirst /*vz->begin()*/ && it_values[ordre[a - 2]] != vz->end())
                ++it_values[ordre[a - 2]];
            //}
        }
        ++itvzfirst;
        ++itvznext;
        if (ToulBar2::verbose >= 3)
            cout << "--\n";
    }
    return sev;
}

void NaryConstraint::separate(EnumeratedVariable* vx, EnumeratedVariable* vz)
{
    static Tuple t, tX, tZ;
    Cost cost, minCost = MAX_COST;
    int a = arity_;
    t.resize(a);
    tX.resize(a - 1);
    tZ.resize(a - 1);
    int index, k;
    EnumeratedVariable* var = NULL;
    Constraint* existX = NULL;
    Constraint* existZ = NULL;
    Constraint* naryx = NULL;
    Constraint* naryz = NULL;
    EnumeratedVariable* scope_in[a];
    EnumeratedVariable* subscopeX[a - 1];
    EnumeratedVariable* subscopeZ[a - 1];
    int scopeX[a - 1];
    int scopeZ[a - 1];

    first(vx, vz);

    // initialisation des scopes
    scope_in[0] = vx;
    scopeX[0] = vx->wcspIndex;
    subscopeX[0] = vx;

    k = 1;

    for (int j = 0; j < a; j++) {
        var = (EnumeratedVariable*)getVar(j);
        if (var != vx && var != vz) {
            scope_in[k] = var;

            subscopeX[k] = var;
            scopeX[k] = var->wcspIndex;

            subscopeZ[k - 1] = var;
            scopeZ[k - 1] = var->wcspIndex;

            ++k;
        }
    }

    scope_in[a - 1] = vz;

    scopeZ[a - 2] = vz->wcspIndex;
    subscopeZ[a - 2] = vz;

    // creation de la nouvelle contrainte
    if (a == 3) {
        int size = subscopeX[0]->getDomainInitSize() * subscopeX[1]->getDomainInitSize();
        vector<Cost> costs(size, 0);
        existX = subscopeX[0]->getConstr(subscopeX[1]);
        naryx = wcsp->newBinaryConstr(subscopeX[0], subscopeX[1], costs);
    } else if (a == 4) {
        int size = subscopeX[0]->getDomainInitSize() * subscopeX[1]->getDomainInitSize() * subscopeX[2]->getDomainInitSize();
        vector<Cost> costs(size, 0);
        existX = subscopeX[0]->getConstr(subscopeX[1], subscopeX[2]);
        naryx = wcsp->newTernaryConstr(subscopeX[0], subscopeX[1], subscopeX[2], costs);
    } else {
        index = wcsp->postNaryConstraintBegin(scopeX, a - 1, wcsp->getUb(), size() / vz->getDomainInitSize(), true);
        naryx = wcsp->getCtr(index);
    }

    minCost = MAX_COST;
    int i;
    if (ToulBar2::verbose >= 3) {
        cout << "\n[ ";
        for (int j = 0; j < a - 1; j++)
            cout << scopeX[j] << " ";
        cout << " ]" << endl;
    }
    while (it_values[0] != scope_in[0]->end()) {
        for (int j = 0; j < a - 1; j++) {
            t[j] = scope_in[j]->toIndex(*it_values[j]);
            tX[j] = scope_in[j]->toIndex(*it_values[j]);
        }
        t[a - 1] = scope_in[a - 1]->toIndex(*it_values[a - 1]);

        cost = eval(t, scope_in);
        if (cost < minCost)
            minCost = cost;
        int finished = false;
        i = a - 1;
        while (!finished) {

            ++it_values[i];
            finished = it_values[i] != scope_in[i]->end();
            if (!finished) {
                if (i == a - 1) {
                    if (ToulBar2::verbose >= 3) {
                        for (unsigned int j = 0; j < tX.size(); j++) {
                            cout << tX[j] << " ";
                        }
                        cout << "  " << minCost << endl;
                    }
                    naryx->setTuple(tX, minCost);
                    minCost = MAX_COST;
                }
                if (i > 0) {
                    it_values[i] = scope_in[i]->begin();
                    i--;
                } else
                    finished = true;
            }
        }
    }

    if (ToulBar2::verbose >= 3)
        cout << "-----------------------------------" << endl;
    if (a == 3) {
        int size = subscopeZ[0]->getDomainInitSize() * subscopeZ[1]->getDomainInitSize();
        vector<Cost> costs(size, 0);
        existZ = subscopeZ[0]->getConstr(subscopeZ[1]);
        naryz = wcsp->newBinaryConstr(subscopeZ[0], subscopeZ[1], costs);
    } else if (a == 4) {

        int size = subscopeZ[0]->getDomainInitSize() * subscopeZ[1]->getDomainInitSize() * subscopeZ[2]->getDomainInitSize();
        vector<Cost> costs(size, 0);
        existZ = subscopeZ[0]->getConstr(subscopeZ[1], subscopeZ[2]);
        naryz = wcsp->newTernaryConstr(subscopeZ[0], subscopeZ[1], subscopeZ[2], costs);
    } else {
        index = wcsp->postNaryConstraintBegin(scopeZ, a - 1, wcsp->getUb(), size() / vx->getDomainInitSize(), true);
        naryz = wcsp->getCtr(index);
    }

    if (ToulBar2::verbose >= 3) {
        cout << "[ ";
        for (int j = 0; j < a - 1; j++)
            cout << scopeZ[j] << " ";
        cout << " ]" << endl;
    }
    first(vx, vz);

    Cost diffcost;
    Cost cX;
    while (it_values[1] != scope_in[1]->end()) {
        for (int j = 1; j < a; j++) {
            t[j] = scope_in[j]->toIndex(*it_values[j]);
            tZ[j - 1] = scope_in[j]->toIndex(*it_values[j]);
        }
        // t[0] = scope_in[0]->toIndex(*it_values[0]);
        EnumeratedVariable::iterator it = scope_in[0]->begin();
        do {
            t[0] = scope_in[0]->toIndex(*it);
            cost = eval(t, scope_in);
            cX = naryx->evalsubstr(t, naryx);
            diffcost = squareminus(cost, cX, wcsp->getUb());
            ++it;
            // cout << t[0] << endl;
        } while (it != scope_in[0]->end() && cost >= wcsp->getUb() && cX >= wcsp->getUb());

        cost = eval(t, scope_in);
        Cost cX = naryx->evalsubstr(t, naryx);
        diffcost = squareminus(cost, cX, wcsp->getUb());
        if (ToulBar2::verbose >= 3) {
            for (unsigned int j = 0; j < tZ.size(); j++) {
                cout << tZ[j] << " ";
            }

            cout << "  " << diffcost << endl;
        }
        assert(diffcost >= MIN_COST);
        naryz->setTuple(tZ, diffcost);
        int finished = false;
        i = a - 1;
        while (!finished) {

            ++it_values[i];
            finished = it_values[i] != scope_in[i]->end();
            if (!finished) {
                if (i > 1) {
                    it_values[i] = scope_in[i]->begin();
                    i--;
                } else
                    finished = true;
            }
        }
    }

    if (ToulBar2::verbose >= 3)
        cout << "-----------------------------------" << endl;
    assert(verifySeparate(naryx, naryz));

    if (!naryx->universal()) {
        if (ToulBar2::verbose == 1) {
            cout << "\n[ ";
            for (int j = 0; j < a - 1; j++)
                cout << scopeX[j] << " ";
            cout << " ]" << endl;
        }
        if (a == 3) {
            assert(naryx->isBinary());
            if (existX && !existX->universal()) {
                assert(existX->isBinary());
                ((BinaryConstraint*)naryx)->addCosts((BinaryConstraint*)existX);
                existX->deconnect();
            }
        } else if (a == 4) {
            assert(naryx->isTernary());
            if (existX && !existX->universal()) {
                assert(existX->isTernary());
                ((TernaryConstraint*)naryx)->addCosts((TernaryConstraint*)existX);
                existX->deconnect();
            }
        }
        naryx->reconnect();
        naryx->propagate();
    } else {
        naryx->deconnect();
    }
    if (!naryz->universal()) {
        if (ToulBar2::verbose == 1) {
            cout << "[ ";
            for (int j = 0; j < a - 1; j++)
                cout << scopeZ[j] << " ";
            cout << " ]" << endl;
        }
        if (a == 3) {
            assert(naryz->isBinary());
            if (existZ && !existZ->universal()) {
                assert(existZ->isBinary());
                ((BinaryConstraint*)naryz)->addCosts((BinaryConstraint*)existZ);
                existZ->deconnect();
            }
        } else if (a == 4) {
            assert(naryz->isTernary());
            if (existZ && !existZ->universal()) {
                assert(existZ->isTernary());
                ((TernaryConstraint*)naryz)->addCosts((TernaryConstraint*)existZ);
                existZ->deconnect();
            }
        }
        naryz->reconnect();
        naryz->propagate();
    } else {
        naryz->deconnect();
    }
    deconnect();
}

// THIS CODE IS NEVER USED!!!
// void NaryConstraint::permute( EnumeratedVariable** scope_in )
//{
//    TUPLES* pf_old = pf;
//    pf = new TUPLES;
//
//    TUPLES::iterator it = pf_old->begin();
//    while(it != pf_old->end()) {
//        Tuple s(it->first);
//        setTuple(s, it->second, scope_in );
//        it++;
//    }
//
//    delete pf_old;
//    for(int i=0; i<arity_; i++)
//    {
//        map<int,int>::iterator it_pos =  scope_inv.find(scope_in[i]->wcspIndex);
//        int i_old = it_pos->second;
//        it_pos->second = i;
//
//        scope_inv[ scope_in[i]->wcspIndex ] = i;
//        scope[i] = scope_in[i];
//
//        DLink<ConstraintLink>* l = links[i];
//        links[i] = links[i_old];
//        links[i_old] = l;
//
//        links[i]->content.scopeIndex = i;
//        l->content.scopeIndex = i_old;
//    }
//}

// for adding a tuple in f
// scope_in contains the order of the values in Tuple tin
// void NaryConstraint::setTuple( const Tuple& tin, Cost c, EnumeratedVariable** scope_in )
//{
//    assert(scope_in);
//    Tuple t(tin);
//    for(int i = 0; i < arity_; i++) {
//        int pos = getIndex(scope_in[i]);
//        t[pos] = tin[i];
//    }
//    if (pf) (*pf)[t] = c;
//    else costs[getCostsIndex(t)] = c;
//}

// void NaryConstraint::addtoTuple( const Tuple& tin, Cost c, EnumeratedVariable** scope_in )
//{
//     assert(scope_in);
//     Tuple t(tin);
//     for(int i = 0; i < arity_; i++) {
//         int pos = getIndex(scope_in[i]);
//         t[pos] = tin[i];
//     }
//     if (pf) (*pf)[t] += c;
//     else costs[getCostsIndex(t)] += c;
// }

void NaryConstraint::addtoTuples(EnumeratedVariable* x, Value v, Cost costi)
{
    static Tuple tuple;
    Cost Top = wcsp->getUb();
    int tindex = getIndex(x);
    assert(tindex >= 0);
    Cost cost;
    if (getDefCost() < Top) {
        firstlex();
        while (nextlex(tuple, cost)) {
            if (x->toValue(tuple[tindex]) == v) {
                if (cost < Top)
                    setTuple(tuple, cost + costi);
            }
        }
    } else {
        first();
        while (next(tuple, cost)) {
            if (x->toValue(tuple[tindex]) == v) {
                if (cost < Top)
                    setTuple(tuple, cost + costi);
            }
        }
    }
}

void NaryConstraint::addtoTuples(Cost costi)
{
    assert(costi >= MIN_COST || getMinCost() >= -costi);
    static Tuple tuple;
    Cost Top = wcsp->getUb();
    Cost cost;
    Cost olddf = default_cost;
    default_cost = Top;
    if (pf) {
        first();
        while (next(tuple, cost)) {
            if (cost < Top)
                setTuple(tuple, cost + costi); // authorizes resulting costs to be greater than Top
        }
        if (olddf < Top && olddf + costi >= MIN_COST)
            default_cost = olddf + costi;
    } else {
        firstlex();
        while (nextlex(tuple, cost)) {
            if (cost < Top)
                setTuple(tuple, cost + costi);
        }
    }
}

void NaryConstraint::setInfiniteCost(Cost ub)
{
    Cost mult_ub = ((wcsp->getUb() < (MAX_COST / MEDIUM_COST)) ? (max(LARGE_COST, wcsp->getUb() * MEDIUM_COST)) : wcsp->getUb());
    if (pf) {
        for (TUPLES::iterator it = pf->begin(); it != pf->end(); ++it) {
            Cost c = it->second; // TODO: add corresponding unary costs
            if (CUT(c, ub))
                it->second = mult_ub;
        }
    } else {
        for (ptrdiff_t idx = 0; idx < costSize; idx++) {
            Cost c = costs[idx]; // TODO: add corresponding unary costs
            if (CUT(c, ub))
                costs[idx] = mult_ub;
        }
    }
    if (CUT(default_cost, ub))
        default_cost = mult_ub;
}

void NaryConstraint::insertSum(const Tuple& t1, Cost c1, Constraint* ctr1, const Tuple& t2, Cost c2, Constraint* ctr2, bool bFilters)
{
    Cost Top = wcsp->getUb();
    if (c1 >= Top)
        return;
    if (c2 >= Top)
        return;
    Cost csum = c1 + c2;

    for (int i = 0; i < arity_; i++) {
        EnumeratedVariable* v = scope[i];
        int pos = i;
        int pos1 = ctr1->getIndex(v);
        int pos2 = ctr2->getIndex(v);

        if ((pos1 >= 0) && (pos2 >= 0)) {
            if (t1[pos1] != t2[pos2])
                return;
            evalTuple[pos] = t1[pos1];
        } else if (pos1 >= 0)
            evalTuple[pos] = t1[pos1];
        else if (pos2 >= 0)
            evalTuple[pos] = t2[pos2];
        else
            evalTuple[pos] = 0;

        Cost unaryc = v->getCost(v->toValue(evalTuple[pos]));
        if (unaryc >= Top)
            return;
        csum += unaryc;
        if (csum >= Top)
            return;
    }

    if (bFilters && filters && (default_cost >= Top)) {
        ConstraintSet::iterator it = filters->begin();
        while (it != filters->end()) {
            Constraint* ctr = *it;
            if (ctr->connected()) {
                assert(ctr != this);
                Cost c = ctr->evalsubstr(evalTuple, this);
                if (c >= Top)
                    return;
                csum += c;
            }
            if (csum >= Top)
                return;
            ++it;
        }
    }

    if (pf)
        (*pf)[evalTuple] = c1 + c2;
    else
        costs[getCostsIndex(evalTuple)] = c1 + c2;
}

// THIS CODE IS NEVER USED!!!
// TODO: make it compatible with table representation rather than Map
// void NaryConstraint::sum( NaryConstraint* nary )
//{
//    deconnect(true);
//
//    map<int,int> snew;
//    set_union( scope_inv.begin(), scope_inv.end(),
//            nary->scope_inv.begin(), nary->scope_inv.end(),
//            inserter(snew, snew.begin()) );
//
//    arity_ = snew.size();
//    EnumeratedVariable** scope1 = scope;
//    DLink<ConstraintLink>** links1 = links;
//    scope = new EnumeratedVariable* [arity_];
//    links = new DLink<ConstraintLink>* [arity_];
//
//    int i = 0;
//    map<int,int>::iterator its = snew.begin();
//    while(its != snew.end()) {
//        EnumeratedVariable* var = (EnumeratedVariable*) wcsp->getVar(its->first);
//        its->second = i;
//        scope[i] =  var;
//        int index1 = getIndex(var);
//        if(index1 >= 0) {
//            links[i] = links1[index1];
//            ConstraintLink e = {this, i};
//            links[i]->content = e;
//        }
//        else links[i] = nary->links[ nary->getIndex(var) ];
//
//        i++;
//        its++;
//    }
//
//    TUPLES& f1 = *pf;
//    TUPLES& f2 = *nary->pf;
//    TUPLES::iterator  it1 = f1.begin();
//    TUPLES& f = * new TUPLES;
//    pf = &f;
//
//    Tuple t1,t2;
//    Cost c1,c2;
//    while(it1 != f1.end()) {
//        t1 = it1->first;
//        c1 =  it1->second;
//        TUPLES::iterator  it2 = f2.begin();
//        while(it2 != f2.end()) {
//            t2 = it2->first;
//            c2 =  it2->second;
//            insertSum(t1, c1, this, t2, c2, nary);
//            it2++;
//        }
//        it1++;
//    }
//
//    scope_inv = snew;
//    delete [] scope1;
//    delete [] links1;
//
//    reconnect();
//}

// Projection of variable x of the nary constraint
// complexity O(2|f|)
// this function is independent of the search
void NaryConstraint::project(EnumeratedVariable* x)
{
    int xindex = getIndex(x);
    if (xindex < 0)
        return;
    assert(x->getDegree() == 1);

    Cost Top = wcsp->getUb();
    Cost negcost = 0;

    if (pf) {
        Tuple t, tnext, tproj, tswap;
        Cost c;
        Value val;
        TUPLES fproj;
        TUPLES::iterator it;
        // First part of the projection: complexity O(|f|) we swap positions between the projected variable and the last variable
        while (!pf->empty()) {
            it = pf->begin();
            t = it->first;
            c = it->second;
            assert(x->getDegree() == 1);
            val = x->toValue(t[xindex]);
            c += ((x->canbe(val)) ? (x->getCost(val)) : Top);
            if (c > Top)
                c = Top;
            tswap = t;
            tValue a = tswap[arity_ - 1];
            tswap[arity_ - 1] = tswap[xindex];
            tswap[xindex] = a;
            fproj[tswap] = c;
            pf->erase(it);
        }

        // Second part of the projection: complexity O(|f|) as the projected variable is in the last position,
        // it is sufficient to look for tuples with the same arity-1 prefix. If there are less than d (domain of
        // the projected variable) tuples, we have also to perform the minimum with default_cost
        // this is only true when the tuples are LEXICOGRAPHICALY ordered
        set<Value> markValue;
        it = fproj.begin();
        if (it != fproj.end()) {
            t = it->first;
            c = it->second;
            val = x->toValue(t[arity_ - 1]);
            t.resize(arity_ - 1);
            markValue.clear();
            if (x->canbe(val)) {
                markValue.insert(val);
            }
            bool end = false;
            while (!end) {
                it++;
                end = (it == fproj.end());
                bool sameprefix = false;

                Cost cnext = MAX_COST;
                Value valnext = -1;
                if (!end) {
                    tnext = it->first;
                    cnext = it->second;
                    valnext = x->toValue(tnext[arity_ - 1]);
                    tnext.resize(arity_ - 1);
                    sameprefix = (t == tnext);
                }
                if (!sameprefix) {
                    if (markValue.size() < x->getDomainSize()) {
                        for (EnumeratedVariable::iterator itv = x->begin(); itv != x->end(); ++itv) {
                            if (markValue.find(*itv) == markValue.end()) {
                                Cost udefcost = default_cost + x->getCost(*itv);
                                if (ToulBar2::isZ) {
                                    c = wcsp->LogSumExp(c, udefcost);
                                } else if (udefcost < c)
                                    c = udefcost;
                            }
                        }
                    }
                    if (ToulBar2::isZ && c < negcost)
                        negcost = c;
                    if (c != default_cost || ToulBar2::isZ) {
                        (*pf)[t] = c;
                    }
                    t = tnext;
                    c = cnext;
                    val = valnext;
                    markValue.clear();
                    if (x->canbe(val)) {
                        markValue.insert(val);
                    }
                } else {
                    if (x->canbe(valnext)) {
                        markValue.insert(valnext);
                    }
                    if (ToulBar2::isZ) {
                        c = wcsp->LogSumExp(c, cnext);
                    } else if (cnext < c)
                        c = cnext;
                }
            }
        }

    } else {
        ptrdiff_t sz = costSize / x->getDomainInitSize();
        Cost* costs_ = new Cost[sz];
        std::fill(costs_, costs_ + sz, Top);
        int a = arity_;
        vector<unsigned int> t(a, 0);
        for (ptrdiff_t idx = 0; idx < costSize; idx++) {
            bool ok = true;
            for (int i = 0; i < a && ok; i++) {
                EnumeratedVariable* var = (EnumeratedVariable*)getVar(i);
                ok = ok && var->canbe(var->toValue(t[i]));
            }
            Cost c = ((ok) ? costs[idx] : Top);
            if (ok) {
                c += x->getCost(x->toValue(t[xindex]));
                if (c > Top)
                    c = Top;
            }

            vector<unsigned int> tswap(t);
            tswap[xindex] = tswap[a - 1];
            ptrdiff_t idx_ = 0;
            ptrdiff_t base_ = 1;
            for (int i = a - 2; i >= 0; --i) {
                idx_ += tswap[i] * base_;
                base_ *= ((EnumeratedVariable*)getVar((i == xindex) ? (a - 1) : i))->getDomainInitSize();
            }
            assert(base_ == sz);
            assert(idx_ < sz);
            assert(idx_ >= 0);
            if (ToulBar2::isZ) {
                c = wcsp->LogSumExp(costs_[idx_], c);
                if (c < negcost)
                    negcost = c;
                costs_[idx_] = c;
            } else
                costs_[idx_] = min(costs_[idx_], c);

            int i = a - 1;
            while (i >= 0 && t[i] == ((EnumeratedVariable*)getVar(i))->getDomainInitSize() - 1) {
                t[i] = 0;
                i--;
            }
            if (i >= 0)
                t[i]++;
        }
        costSize = sz;
        delete[] costs;
        costs = costs_;
    }
    assert(negcost <= 0);
    if (negcost < 0) {
        if (pf) {
            for (TUPLES::iterator it = pf->begin(); it != pf->end(); ++it) {
                (*pf)[(*it).first] = (*it).second - negcost;
            }
        } else {
            for (ptrdiff_t idx = 0; idx < costSize; idx++) {
                costs[idx] -= negcost;
            }
        }
        default_cost -= negcost;
        assert(default_cost >= MIN_COST); // THIS TEST COULD BE REMOVED BUT NEGATIVE DEFAULT COST EFFECTS ARE UNSPECIFIED
        wcsp->decreaseLb(negcost);
    }

    // reduce nary by replacing variable x by the last variable in the scope (if needed) and then remove it from the scope by assuming it will be assigned
    if (xindex < arity_ - 1) {
        // swap of links
        DLink<ConstraintLink>* linkx = links[xindex];
        links[xindex] = links[arity_ - 1];
        links[arity_ - 1] = linkx;
        // swap of scope array
        scope[xindex] = scope[arity_ - 1];
        scope[arity_ - 1] = x;
        // update of links indexs
        links[xindex]->content.scopeIndex = xindex;
        links[arity_ - 1]->content.scopeIndex = arity_ - 1;
        scope_inv[scope[xindex]->wcspIndex] = xindex;
    }

    // shift left scope_dac array starting at position of variable x in scope_dac
    int dacindex = 0;
    while (dacindex < arity_ - 1 && scope_dac[dacindex] != x) {
        dacindex++;
    }
    assert(dacindex == arity_ - 1 || scope_dac[dacindex] == x);
    for (; dacindex < arity_ - 1; dacindex++) {
        scope_dac[dacindex] = scope_dac[dacindex + 1];
    }

    if (x->unassigned()) {
        deconnect(arity_ - 1);
    }
    scope_inv.erase(x->wcspIndex);
    arity_--;
    evalTuple.resize(arity_);
}

// Projects out all variables except x,y,z
// and gives the result at fproj
void NaryConstraint::projectxyz(EnumeratedVariable* x,
    EnumeratedVariable* y,
    EnumeratedVariable* z,
    TUPLES& fproj)
{
    assert(CUT(default_cost, wcsp->getUb()));
    static Tuple t;
    Tuple txyz(3, 0);
    Cost c;
    TUPLES::iterator itproj;

    // compute in one pass of all tuples the projection
    first();
    while (next(t, c)) {
        if (c == default_cost)
            continue;
        txyz[0] = t[getIndex(x)];
        txyz[1] = t[getIndex(y)];
        txyz[2] = t[getIndex(z)];
        itproj = fproj.find(txyz);
        if (itproj != fproj.end()) {
            if (c < itproj->second)
                fproj[txyz] = c;
        } else {
            fproj[txyz] = c;
        }
    }

    // finally we subtract the projection from the initial function
    first();
    while (next(t, c)) {
        if (c == default_cost)
            continue;
        txyz[0] = t[getIndex(x)];
        txyz[1] = t[getIndex(y)];
        txyz[2] = t[getIndex(z)];
        itproj = fproj.find(txyz);
        if (itproj != fproj.end()) {
            assert(CUT(c, itproj->second));
            //        if (!CUT(c + wcsp->getLb(), wcsp->getUb())) {
            if (pf)
                (*pf)[t] -= itproj->second;
            else
                costs[getCostsIndex(t)] -= itproj->second;
            //        }
        } else
            assert(false);
    }
}

// Projects out all variables except x,y
// and gives the result at fproj
void NaryConstraint::projectxy(EnumeratedVariable* x,
    EnumeratedVariable* y,
    TUPLES& fproj)
{
    assert(CUT(default_cost, wcsp->getUb()));
    static Tuple t;
    Tuple txy(2, 0);
    Cost c;
    TUPLES::iterator itproj;

    // compute in one pass of all tuples the projection
    first();
    while (next(t, c)) {
        if (c == default_cost)
            continue;
        txy[0] = t[getIndex(x)];
        txy[1] = t[getIndex(y)];
        itproj = fproj.find(txy);
        if (itproj != fproj.end()) {
            if (c < itproj->second)
                fproj[txy] = c;
        } else {
            fproj[txy] = c;
        }
    }

    // finally we subtract the projection from the initial function
    first();
    while (next(t, c)) {
        if (c == default_cost)
            continue;
        txy[0] = t[getIndex(x)];
        txy[1] = t[getIndex(y)];
        itproj = fproj.find(txy);
        if (itproj != fproj.end()) {
            assert(CUT(c, itproj->second));
            //		  if (!CUT(c + wcsp->getLb(), wcsp->getUb())) {
            if (pf)
                (*pf)[t] -= itproj->second;
            else
                costs[getCostsIndex(t)] -= itproj->second;
            //		  }
        } else
            assert(false);
    }
}

void NaryConstraint::preproject3(TernaryConstraint* ctr)
{
    EnumeratedVariable* x = (EnumeratedVariable*)ctr->getVar(0);
    assert(getIndex(x) >= 0);
    EnumeratedVariable* y = (EnumeratedVariable*)ctr->getVar(1);
    assert(getIndex(y) >= 0);
    EnumeratedVariable* z = (EnumeratedVariable*)ctr->getVar(2);
    assert(getIndex(z) >= 0);

    TUPLES fproj;
    projectxyz(x, y, z, fproj);

    unsigned int sizex = x->getDomainInitSize();
    unsigned int sizey = y->getDomainInitSize();
    unsigned int sizez = z->getDomainInitSize();

    vector<Cost> xyz;
    for (unsigned int a = 0; a < sizex; a++)
        for (unsigned int b = 0; b < sizey; b++)
            for (unsigned int c = 0; c < sizez; c++)
                xyz.push_back(default_cost);

    TUPLES::iterator it = fproj.begin();
    while (it != fproj.end()) {
        unsigned int a = it->first[0];
        unsigned int b = it->first[1];
        unsigned int c = it->first[2];
        xyz[a * sizey * sizez + b * sizez + c] = it->second;
        it++;
    }
    if (fproj.size() > 0 || default_cost > MIN_COST) {
        if (ToulBar2::verbose >= 3) {
            cout << "[" << Store::getDepth() << ",W" << wcsp->getIndex() << "] preproject nary " << this << " to ternary (" << x->getName() << "," << y->getName() << "," << z->getName() << ")" << endl;
        }
        wcsp->postTernaryConstraint(x->wcspIndex, y->wcspIndex, z->wcspIndex, xyz);
    }
}

inline bool cmp_pairvars(pair<EnumeratedVariable*, EnumeratedVariable*> pv1, pair<EnumeratedVariable*, EnumeratedVariable*> pv2)
{
    return (pv1.first->wcspIndex < pv2.first->wcspIndex || (pv1.first->wcspIndex == pv2.first->wcspIndex && pv1.second->wcspIndex < pv2.second->wcspIndex));
}
void NaryConstraint::preprojectall2()
{
    Tuple t;
    assert(connected());
    assert(CUT(default_cost, wcsp->getUb()));

    for (int i = 0; i < arity_; i++) {
        for (int j = i + 1; j < arity_; j++) {
            EnumeratedVariable* x = scope_dac[i]; // (EnumeratedVariable*)wcsp->getVar((*it1).first);
            EnumeratedVariable* y = scope_dac[j]; // (EnumeratedVariable*)wcsp->getVar((*it2).first);

            TUPLES fproj;
            projectxy(x, y, fproj);

            vector<Cost> xy;
            unsigned int a, b;
            unsigned int sizex = x->getDomainInitSize();
            unsigned int sizey = y->getDomainInitSize();

            for (a = 0; a < sizex; a++)
                for (b = 0; b < sizey; b++)
                    xy.push_back(default_cost);

            TUPLES::iterator it = fproj.begin();
            while (it != fproj.end()) {
                t = it->first;
                a = t[0];
                b = t[1];
                xy[a * sizey + b] = it->second;
                it++;
            }
            if (fproj.size() > 0 || default_cost > MIN_COST) {
                if (ToulBar2::verbose >= 3) {
                    cout << "[" << Store::getDepth() << ",W" << wcsp->getIndex() << "] preproject nary " << this << " to binary (" << x->getName() << "," << y->getName() << ")" << endl;
                }
                wcsp->postBinaryConstraint(x->wcspIndex, y->wcspIndex, xy);
            }
            if (deconnected())
                return;
        }
    }
}

double NaryConstraint::computeTightness()
{
    int count = 0;
    double sum = 0;
    double tight = -1;
    Cost* costs_ = new Cost[size()];
    if (pf) {
        TUPLES::iterator it = pf->begin();
        while (it != pf->end()) {
            Cost c = it->second;
            sum += to_double(min(wcsp->getUb(), c));
            costs_[count] = min(wcsp->getUb(), c);
            count++;
            it++;
        }
    } else {
        for (ptrdiff_t idx = 0; idx < costSize; idx++) {
            Cost c = costs[idx];
            sum += to_double(min(wcsp->getUb(), c));
            costs_[count] = min(wcsp->getUb(), c);
            count++;
        }
    }
    Long psize = getDomainSizeProduct();
    if (psize >= LONGLONG_MAX) {
        tight = to_double(min(wcsp->getUb(), default_cost));
    } else {
        if (ToulBar2::weightedTightness == 2) {
            tight = to_double(stochastic_selection<Cost>(costs_, 0, count - 1, count / 2)); // TO BE IMPPROVED!!!
        } else {
            if (psize > count) {
                tight = ((double)default_cost * (psize - count) + sum) / (double)psize;
            } else {
                tight = sum / (double)count;
            }
        }
    }
    delete[] costs_;
    return tight;
}

void NaryConstraint::print(ostream& os)
{
    os << endl
       << this << " f(";

    int unassigned_ = 0;
    Long totaltuples = 1;
    for (int i = 0; i < arity_; i++) {
        if (scope[i]->unassigned())
            unassigned_++;
        os << wcsp->getName(scope[i]->wcspIndex);
        if (i < arity_ - 1)
            os << ",";
        totaltuples = totaltuples * scope[i]->getDomainInitSize();
    }
    os << ")    ";
    if (ToulBar2::weightedDegree) {
        os << "/" << getConflictWeight();
        for (int i = 0; i < arity_; i++) {
            os << "," << conflictWeights[i];
        }
    }
    if (pf)
        os << " |f| = " << size() << " / " << totaltuples;
    else
        os << " |f| = " << size();
    os << "   default_cost: " << default_cost;
    os << "   arity: " << arity_;
    os << "   unassigned: " << getNonAssigned() << "/" << unassigned_ << "         ";

    //	assert(getNonAssigned() == unassigned_); // not valid when used with assignLS

    /*TSCOPE::iterator it = scope_inv.begin();
        while(it != scope_inv.end()) {
                os << "(" << it->first << ",idx: " << it->second << ") ";
                ++it;
        }*/
    os << endl;

    if (ToulBar2::verbose >= 4) {
        if (pf) {
            os << "tuples: {";
            TUPLES::iterator it = pf->begin();
            while (it != pf->end()) {
                Tuple t = it->first;
                Cost c = it->second;
                it++;
                os << "<";
                for (unsigned int i = 0; i < t.size(); i++) {
                    os << t[i];
                    if (i < t.size() - 1)
                        os << " ";
                }
                os << "," << c << ">";
                if (it != pf->end())
                    os << " ";
            }
            os << "} " << endl;
        } else {
            os << "costs: [";
            for (ptrdiff_t idx = 0; idx < costSize; idx++) {
                if (ToulBar2::verbose >= 8) {
                    Long tvalues = idx;
                    Long product = costSize;
                    os << "(" << idx << ",[";
                    for (int i = 0; i < arity_; i++) {
                        product /= scope[i]->getDomainInitSize();
                        os << scope[i]->toValue(tvalues / product);
                        tvalues %= product;
                        if (i < arity_ - 1)
                            os << ",";
                    }
                    assert(product == 1);
                    os << "]," << costs[idx] << ")";
                    if (idx < costSize - 1)
                        os << ",";
                } else {
                    os << " " << costs[idx];
                }
            }
            os << "] " << endl;
        }
    }
}

void NaryConstraint::dump(ostream& os, bool original)
{
    if (original) {
        os << arity_;
        for (int i = 0; i < arity_; i++)
            os << " " << scope[i]->wcspIndex;
        os << " " << default_cost << " " << size() << endl;
        if (pf) {
            TUPLES::iterator it = pf->begin();
            while (it != pf->end()) {
                Tuple t = it->first;
                Cost c = it->second;
                it++;
                for (unsigned int i = 0; i < t.size(); i++) {
                    os << scope[i]->toValue(t[i]) << " ";
                }
                os << c << endl;
            }
        } else {
            int a = arity_;
            vector<unsigned int> t(a, 0);
            for (ptrdiff_t idx = 0; idx < costSize; idx++) {
                for (int i = 0; i < a; i++) {
                    os << scope[i]->toValue(t[i]) << " ";
                }
                os << costs[idx] << endl;
                int i = a - 1;
                while (i >= 0 && t[i] == scope[i]->getDomainInitSize() - 1) {
                    t[i] = 0;
                    i--;
                }
                if (i >= 0)
                    t[i]++;
            }
        }
    } else {
        os << getNonAssigned();
        for (int i = 0; i < arity_; i++)
            if (scope[i]->unassigned())
                os << " " << scope[i]->getCurrentVarId();
        Tuple tuple;
        Cost cost;
        Long nbtuples = 0;
        first();
        while (next(tuple, cost)) {
            nbtuples++;
        }
        os << " " << default_cost << " " << nbtuples << endl;
        first();
        while (next(tuple, cost)) {
            for (unsigned int i = 0; i < tuple.size(); i++) {
                if (scope[i]->unassigned())
                    os << scope[i]->toCurrentIndex(scope[i]->toValue(tuple[i])) << " ";
            }
            os << ((original) ? cost : min(wcsp->getUb(), cost)) << endl;
        }
    }
}

void NaryConstraint::dump_CFN(ostream& os, bool original)
{

    bool printed = false;
    os << "\"F_";

    if (original) {
        printed = false;
        for (int i = 0; i < arity_; i++) {
            if (printed)
                os << "_";
            os << scope[i]->wcspIndex;
            printed = true;
        }

        os << "\":{\"scope\":[";
        printed = false;
        for (int i = 0; i < arity_; i++) {
            if (printed)
                os << ",";
            os << "\"" << name2cfn(scope[i]->getName()) << "\"";
            printed = true;
        }
        os << "],\"defaultcost\":" << wcsp->DCost2Decimal(wcsp->Cost2RDCost(default_cost)) << ",\n\"costs\":[";
        if (pf) {
            TUPLES::iterator it = pf->begin();
            printed = false;
            while (it != pf->end()) {
                Tuple t = it->first;
                Cost c = it->second;
                it++;
                os << endl;
                for (unsigned int i = 0; i < t.size(); i++) {
                    if (printed)
                        os << ",";
                    os << ((scope[i]->isValueNames()) ? name2cfn(scope[i]->getValueName(t[i])) : to_string(t[i]));
                    printed = true;
                }
                os << "," << wcsp->DCost2Decimal(wcsp->Cost2RDCost(c));
            }
        } else {
            int a = arity_;
            vector<unsigned int> t(a, 0);
            printed = false;
            for (ptrdiff_t idx = 0; idx < costSize; idx++) {
                os << endl;
                for (int i = 0; i < a; i++) {
                    if (printed)
                        os << ",";
                    os << ((scope[i]->isValueNames()) ? name2cfn(scope[i]->getValueName(t[i])) : to_string(t[i]));
                    printed = true;
                }
                os << "," << wcsp->DCost2Decimal(wcsp->Cost2RDCost(costs[idx]));
                int i = a - 1;
                while (i >= 0 && t[i] == scope[i]->getDomainInitSize() - 1) {
                    t[i] = 0;
                    i--;
                }
                if (i >= 0)
                    t[i]++;
            }
        }
    } else {
        for (int i = 0; i < arity_; i++)
            if (scope[i]->unassigned()) {
                if (printed)
                    os << "_";
                os << scope[i]->getCurrentVarId();
                printed = true;
            }
        os << "\":{\"scope\":[";
        printed = false;
        for (int i = 0; i < arity_; i++)
            if (scope[i]->unassigned()) {
                if (printed)
                    os << ",";
                os << "\"" << name2cfn(scope[i]->getName()) << "\"";
                printed = true;
            }
        os << "],\"defaultcost\":" << ((wcsp->getUb() > default_cost) ? wcsp->DCost2Decimal(wcsp->Cost2RDCost(default_cost)) : "inf") << ",\n\"costs\":[";

        Tuple t;
        Cost cost;
        first();
        while (next(t, cost)) {
            os << endl;
            printed = false;
            for (unsigned int i = 0; i < t.size(); i++) {
                if (scope[i]->unassigned()) {
                    if (printed)
                        os << ",";
                    os << scope[i]->toCurrentIndex(scope[i]->toValue(t[i]));
                    printed = true;
                }
            }
            os << "," << ((wcsp->getUb() > cost) ? wcsp->DCost2Decimal(wcsp->Cost2RDCost(cost)) : "inf");
        }
    }
    os << "]},\n";
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

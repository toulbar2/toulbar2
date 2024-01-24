#include "tb2linearconstr.hpp"
#include "core/tb2wcsp.hpp"

// #define verify

LinearConstraint::LinearConstraint(WCSP* wcsp, EnumeratedVariable** scope_in, int arity_in)
    : GlobalConstraint(wcsp, scope_in, arity_in, 0)
    , initTest(false)
    , mip(NULL)
{
#ifdef ILOGCPLEX
    mip = MIP::makeMIP();
#endif
}

Cost LinearConstraint::solveMIP(MIP& mip)
{
    return mip.solve();
}

void LinearConstraint::initStructure()
{

    cost = buildMIP();

    buObj = new int[count];

    propagate();
}

void LinearConstraint::end()
{
    mip->called_time();
    if (deconnected())
        return;
}

void LinearConstraint::checkRemoved(MIP& mip, Cost& cost, vector<int>& rmv)
{ // get current domain from mip and prune other values

    pair<Cost, bool> result;
    vector<int> cDomain;
    // bool deleted = false;
    bool flag = false;
    for (int i = 0; i < arity_; i++) {
        cDomain.clear();
        getDomainFromMIP(mip, i, cDomain); // getDomain
        sort(cDomain.begin(), cDomain.end());

        EnumeratedVariable* y = (EnumeratedVariable*)getVar(i);
        for (EnumeratedVariable::iterator v = y->begin(); v != y->end(); ++v) {

            vector<int>::iterator it = find(cDomain.begin(), cDomain.end(), *v);
            if (it == cDomain.end()) {
                cerr << "non exist a value ?" << endl;
                for (vector<int>::iterator v = cDomain.begin(); v != cDomain.end(); v++) {
                    cerr << *v << " ";
                }
                cerr << endl;
                for (EnumeratedVariable::iterator v = y->begin(); v != y->end(); ++v) {
                    cerr << *v << " ";
                }
                cerr << endl;
                throw InternalError();
            }
            cDomain.erase(it);
            // deleted = true;
        }

        if (!cDomain.empty()) {
            rmv.push_back(i);
            for (vector<int>::iterator v = cDomain.begin(); v != cDomain.end(); v++) {
                int var1 = mapvar[i][*v];
#ifdef ILOGCPLEX
                if (((IlogMIP&)mip).sols->getSize() > 0 && mip.sol(var1) == 1) { // checking if this value is being used
                    flag = true;
                }
#endif
                mip.colUpperBound(var1, 0); // removeDomain
            }
            // deleted = true;
        }
    }
    if (flag) {
        cost = solveMIP(); // solve
    }
}

void LinearConstraint::findProjection(MIP& mip, Cost& cost, int varindex, map<Value, Cost>& delta)
{

    pair<Cost, bool> result;
    delta.clear();
    EnumeratedVariable* x = (EnumeratedVariable*)getVar(varindex);
    for (EnumeratedVariable::iterator j = x->begin(); j != x->end(); ++j) {

        int var1 = mapvar[varindex][*j];
        int tmp = cost;
        cost = tmp = mip.augment(var1); // make sure this value is used...

        assert(tmp >= 0);
        delta[*j] = tmp;
    }
}

void LinearConstraint::augmentStructure(MIP& mip, Cost& cost, int varindex, map<Value, Cost>& delta)
{

    for (map<Value, Cost>::iterator i = delta.begin(); i != delta.end(); i++) {

        int var1 = mapvar[varindex][i->first];
        mip.objCoeff(var1, mip.objCoeff(var1) - i->second); // update unary cost
#ifdef ILOGCPLEX
        if (((IlogMIP&)mip).sols->getSize() > 0 && mip.sol(var1) == 1) { // using this value?
            cost -= i->second;
        }
#endif
    }
}

void LinearConstraint::changeAfterExtend(vector<int>& supports, vector<map<Value, Cost>>& deltas)
{

    bucost = cost;
    if (buObj == NULL) {
        buObj = new int[count];
    }
    for (int i = 0; i < count; i++) {
        buObj[i] = mip->objCoeff(i); // retrieve unary cost
    }
    for (unsigned int i = 0; i < supports.size(); i++) {
        for (map<Value, Cost>::iterator v = deltas[i].begin(); v != deltas[i].end(); v++)
            v->second *= -1;
        augmentStructure(*mip, cost, supports[i], deltas[i]);
        for (map<Value, Cost>::iterator v = deltas[i].begin(); v != deltas[i].end(); v++)
            v->second *= -1;
    }
    cost = solveMIP(); // solve
}

void LinearConstraint::changeAfterProject(vector<int>& supports, vector<map<Value, Cost>>& deltas)
{
    for (unsigned int i = 0; i < supports.size(); i++) {
        augmentStructure(*mip, cost, supports[i], deltas[i]);
    }
}

void LinearConstraint::getDomainFromMIP(MIP& mip, int varindex, vector<int>& domain)
{

    domain.clear();
    for (map<Value, int>::iterator v = mapvar[varindex].begin(); v != mapvar[varindex].end(); v++) {
        if (mip.colUpperBound(v->second) == 1) {
            domain.push_back(v->first);
        }
    }
}

unsigned LinearConstraint::called_time()
{
    return mip->called_time();
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

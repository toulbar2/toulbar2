/** \file tb2linearconstr.hpp
 *  \brief Global cost functions using linear programming for propagation
 *
 */

#ifndef TB2LINEARCONSTR_HPP_
#define TB2LINEARCONSTR_HPP_

#include "tb2globalconstr.hpp"
#include "tb2mipsolver.hpp"

class LinearConstraint : public GlobalConstraint {
protected:
    bool initTest;
    MIP* mip;

    int* buObj;

    Cost cost, bucost;

    int domainSize;

    int count; // total number of domain values in the scope of the global cost function

    map<Value, int>* mapvar;

    // compute the projection from the linear program. store the projected cost in
    // the map delta
    virtual void findProjection(MIP& mip, Cost& cost, int varindex, map<Value, Cost>& delta);
    void findProjection(int varindex, map<Value, Cost>& delta)
    {
        findProjection(*mip, cost, varindex, delta);
    }

    // check whether the linear program corresponding to the current domains
    // remove any edge which is corresponded to an infeasible assignment
    virtual void checkRemoved(MIP& mip, Cost& cost, vector<int>& rmv);
    void checkRemoved(vector<int>& rmv)
    {
        checkRemoved(*mip, cost, rmv);
    }

    virtual void changeAfterExtend(vector<int>& supports, vector<map<Value, Cost>>& deltas);
    virtual void changeAfterProject(vector<int>& supports, vector<map<Value, Cost>>& deltas);
    virtual void undoExtend()
    {
        cost = bucost;
        for (int i = 0; i < count; i++) {
            mip->objCoeff(i, buObj[i]);
        }
    }

    // construct the linear program
    virtual Cost buildMIP(MIP& mip) { return 0; }
    inline Cost buildMIP() { return buildMIP(*mip); }

    inline void augmentMIP(int varindex, map<Value, Cost>& delta)
    {
        augmentStructure(*mip, cost, varindex, delta);
    }

    // compute the minimal of the linear program
    virtual Cost solveMIP(MIP& mip);
    inline Cost solveMIP() { return solveMIP(*mip); }

    // compute the domains of a variable from the linear program
    virtual void getDomainFromMIP(MIP& mip, int varindex, vector<int>& domain);

    // augment the cost to the linear program
    virtual void augmentStructure(MIP& mip, Cost& cost, int varindex, map<Value, Cost>& delta);

    // compute the cost according to the original cost structure
    virtual Cost evalOriginal(const Tuple& s) { return MIN_COST; }
    virtual Cost getMinCost()
    {
        return cost;
    }

public:
    LinearConstraint(WCSP* wcsp, EnumeratedVariable** scope_in, int arity_in);

    ~LinearConstraint() {}

    virtual void read(istream& file, bool mult = true) {}
    virtual void initStructure();
    virtual void end();

    unsigned called_time();
};

#endif /*TB2LINEARCONSTR_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

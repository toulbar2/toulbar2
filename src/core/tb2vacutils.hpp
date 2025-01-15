/** \file tb2vacutils.hpp
 *  \brief Set of useful classes to enforce VAC
 */

#ifndef TB2VACUTILS_HPP_
#define TB2VACUTILS_HPP_

#include "tb2enumvar.hpp"
#include "tb2binconstr.hpp"
//#include "tb2ternaryconstr.hpp"
#include "tb2vac.hpp"

class VACVariable : public EnumeratedVariable {

public:
private:
    VACExtension* vac; /**< Ref. to the VAC data-structures */

    vector<Long> mark; /**< The boolean M used to mark values whose deletion is needed to wipe-out */
    vector<Long> k_timeStamp; /**< timestamp for the counter k (one per value) */
    vector<int> k; /**< Number of cost requests per value for all cost functions */
    vector<int> killer; /**< The killer of each value : the other variable index (binary case)*/
    vector<vector<pair<int, Value>>> PBkillers;
    int maxk; /**< The Max number of cost requests seen on this variable, used for stats */
    Long maxk_timeStamp; /**< timestamp for maxk */

    StoreCost myThreshold; /** The local threshold used to break loops */

    DLink<VariableWithTimeStamp> linkVACQueue;
#ifdef INCREMENTALVAC
    DLink<VariableWithTimeStamp> linkVAC2Queue;
#endif
    DLink<VariableWithTimeStamp> linkSeekSupport;

    void init();

public:
    VACVariable(WCSP* wcsp, string n, Value iinf, Value isup);
    VACVariable(WCSP* wcsp, string n, vector<Value>& dom);
    ~VACVariable();

    bool removeVAC(Value v)
    {
        if (v == inf) {
            if (v == sup)
                return true;
            inf = domain.increase(v + 1);
        } else if (v == sup) {
            if (v == inf)
                return true;
            sup = domain.decrease(v - 1);
        } else if (canbe(v)) {
            domain.erase(v);
        }
#ifdef INCREMENTALVAC
        if (v == maxCostValue || PARTIALORDER) {
            queueNC();
        }
#endif
        return false;
    }

    int getMaxK(Long timeStamp)
    {
        if (maxk_timeStamp < timeStamp)
            return 0;
        else
            return maxk;
    }

    int getK(Value v, Long timeStamp)
    {
        if (k_timeStamp[toIndex(v)] < timeStamp)
            return 0;
        else
            return k[toIndex(v)];
    }
    void setK(Value v, int c, Long timeStamp)
    {
        k[toIndex(v)] = c;
        k_timeStamp[toIndex(v)] = timeStamp;
        if (maxk_timeStamp < timeStamp) {
            maxk = 0;
            maxk_timeStamp = timeStamp;
        }
    }

    void addToK(Value v, int c, Long timeStamp)
    {
        if (k_timeStamp[toIndex(v)] < timeStamp)
            k[toIndex(v)] = c;
        else
            k[toIndex(v)] += c;
        if (maxk_timeStamp < timeStamp)
            maxk = k[toIndex(v)];
        else if (maxk < k[toIndex(v)])
            maxk = k[toIndex(v)];
        maxk_timeStamp = timeStamp;
        k_timeStamp[toIndex(v)] = timeStamp;
    }

    bool isMarked(Value v, Long timeStamp) { return (mark[toIndex(v)] >= timeStamp); }
    void setMark(Value v, Long timeStamp) { mark[toIndex(v)] = timeStamp; }

    int getKiller(Value v) { return killer[toIndex(v)]; }
    void setKiller(Value v, int i) { killer[toIndex(v)] = i; }
    const vector<pair<int, Value>>& getPBkillers(Value v) const { return PBkillers[toIndex(v)]; }
    void setPBkillers(Value v, const vector<pair<int, Value>> & i) { PBkillers[toIndex(v)] = i; }

    Cost getVACCost(Value v)
    {
        Cost c = getCost(v);
        if (isNull(c))
            return MIN_COST;
        else
            return c;
    }

    void setThreshold(Cost c) { myThreshold = c; }
    Cost getThreshold() { return myThreshold; }

    bool isSimplyNull(Cost c) { return (vac->isNull(c)); }
    bool isNull(Cost c) { return (vac->isNull(c) || (c < myThreshold)); }

    void queueVAC() { wcsp->vac->queueVAC(&linkVACQueue); }
#ifdef INCREMENTALVAC
    void queueVAC2()
    {
        wcsp->vac->queueVAC2(&linkVAC2Queue);
    }
#endif
    void queueSeekSupport()
    {
        wcsp->vac->queueSeekSupport(&linkSeekSupport);
    }

    void VACproject(Value v, Cost c) /**< Increases unary cost and may queue for NC enforcing (maintaining maxCost and maxCostValue during VAC-lin which may create unary costs without consuming them) */
    {
        assert(ToulBar2::verbose < 4 || ((cout << "[" << Store::getDepth() << ",W" << wcsp->getIndex() << "] project " << getName() << " (" << v << ") += " << c << endl), true));
        assert(c > MIN_COST);
        costs[toIndex(v)] += c;
        if (ToulBar2::VAClin && !wcsp->knapsackList.empty() && (v == maxCostValue || LUBTEST(maxCost, getCost(v))))
            queueNC();
    }
    void VACextend(Value v, Cost c) /**< Decreases unary cost and may queue for NC enforcing (maintaining maxCost and maxCostValue) */
    {
        assert(ToulBar2::verbose < 4 || ((cout << "[" << Store::getDepth() << ",W" << wcsp->getIndex() << "] extend " << getName() << " (" << v << ") -= " << c << endl), true));
        assert(c > MIN_COST);
        costs[toIndex(v)] -= c;
        if (v == maxCostValue || PARTIALORDER) {
            queueNC();
        }
    }

    bool averaging(); /**< For Min-Sum diffusion */

    void print(ostream& os)
    {
        EnumeratedVariable::print(os);
        cout << " Threshold: " << myThreshold;
    }

    void KilledOne();
};

/**
 * A class that stores information about a binary cost function
 */
class VACBinaryConstraint : public BinaryConstraint {

private:
    vector<int> kX; /**< The k_XY(X,v) counters: nb. of cost request on X,v by this cost function */
    vector<int> kY; /**< The k_XY(Y,v) counters: nb. of cost request on Y,v by this cost function */
    vector<Long> kX_timeStamp;
    vector<Long> kY_timeStamp;

    Cost myThreshold; /** The local threshold used to break loops */

public:
    VACBinaryConstraint(WCSP* wcsp, EnumeratedVariable* xx, EnumeratedVariable* yy, vector<Cost>& tab);
    VACBinaryConstraint(WCSP* wcsp);
    ~VACBinaryConstraint();

    void VACfillElimConstr()
    {
        for (unsigned int a = kX.size(); a < sizeX; a++) {
            kX.push_back(0);
            kX_timeStamp.push_back(0);
        }
        for (unsigned int b = kY.size(); b < sizeY; b++) {
            kY.push_back(0);
            kY_timeStamp.push_back(0);
        }
    }

    int getK(VACVariable* var, Value v, Long timeStamp)
    {
        if (var == (VACVariable*)getVar(0)) {
            if (kX_timeStamp[var->toIndex(v)] < timeStamp)
                return 0;
            else
                return kX[var->toIndex(v)];
        } else {
            if (kY_timeStamp[var->toIndex(v)] < timeStamp)
                return 0;
            else
                return kY[var->toIndex(v)];
        }
    }

    void setK(VACVariable* var, Value v, int c, Long timeStamp)
    {
        if (var == getVar(0)) {
            kX[var->toIndex(v)] = c;
            kX_timeStamp[var->toIndex(v)] = timeStamp;
        } else {
            kY[var->toIndex(v)] = c;
            kY_timeStamp[var->toIndex(v)] = timeStamp;
        }
    }

    void setThreshold(Cost c) { myThreshold = c; }
    Cost getThreshold() { return myThreshold; }

    bool isNull(Cost c)
    {
        VACVariable* xi = (VACVariable*)getVar(0);
        return (xi->isSimplyNull(c) || (c < myThreshold));
    }

    Cost getVACCost(VACVariable* xx, VACVariable* yy, Value v, Value w)
    {
        Cost c = getCost(xx, yy, v, w);
        if (isNull(c))
            return MIN_COST;
        else
            return c;
    }
    void VACproject(VACVariable* x, Value v, Cost c); /**< Modifies Delta counters, then VAC projects on value */
    void VACextend(VACVariable* x, Value v, Cost c); /**< Modifies Delta counters, then VAC extends from value */

    bool revise(VACVariable* var, Value v); /**< AC2001 based Revise for Pass1 : Revise value wrt this cost function */

    void print(ostream& os)
    {
        BinaryConstraint::print(os);
        cout << "Threshold: " << myThreshold << endl;
    }
};

/**
 * A class that stores information about a ternary cost function
 */
// class VACTernaryConstraint : public TernaryConstraint {
//
// private:
//     vector<int> kX; /**< The k_XYZ(X,v) counters: nb. of cost request on X,v by this cost function */
//     vector<int> kY; /**< The k_XYZ(Y,v) counters: nb. of cost request on Y,v by this cost function */
//     vector<int> kZ; /**< The k_XYZ(Z,v) counters: nb. of cost request on Z,v by this cost function */
//     vector<Long> kX_timeStamp;
//     vector<Long> kY_timeStamp;
//     vector<Long> kZ_timeStamp;
//
//     StoreCost myThreshold; /** The local thresold used to break loops */
//
// public:
//     VACTernaryConstraint(WCSP* wcsp, EnumeratedVariable* xx, EnumeratedVariable* yy, EnumeratedVariable* zz, BinaryConstraint* xy, BinaryConstraint* xz, BinaryConstraint* yz, vector<Cost>& tab);
//     VACTernaryConstraint(WCSP* wcsp);
//     ~VACTernaryConstraint();
//
//     int getK(VACVariable* var, Value v, Long timeStamp);
//     void setK(VACVariable* var, Value v, int c, Long timeStamp);
//
//     void setThreshold(Cost c) { myThreshold = c; }
//     Cost getThreshold() { return myThreshold; }
//
//     bool isNull(Cost c);
//
//     Cost getVACCost(VACVariable* xx, VACVariable* yy, VACVariable* zz, Value u, Value v, Value w)
//     {
//         Cost c = getCost(xx, yy, zz, u, v, w);
//         if (isNull(c))
//             return MIN_COST;
//         else
//             return c;
//     }
//     void VACproject(VACVariable* x, Value v, Cost c); /**< Modifies Delta counters, then VAC projects on value */
//     void VACextend(VACVariable* x, Value v, Cost c); /**< Modifies Delta counters, then VAC extends from value */
//
//     bool revise(VACVariable* var, Value v); /**< AC2001 based Revise for Pass1 : Revise value wrt this cost function */
//
//     friend ostream& operator<<(ostream& os, VACTernaryConstraint& c)
//     {
//         return os;
//     }
// };

#endif /*TB2VACUTILS_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

/** \file tb2binconstr.hpp
 *  \brief Binary constraint applied on variables with enumerated domains.
 *  \warning EAC is not applied on duplicated binary constraints
 */

#ifndef TB2BINCONSTR_HPP_
#define TB2BINCONSTR_HPP_

#include "tb2abstractconstr.hpp"
#include "tb2enumvar.hpp"
#include "tb2wcsp.hpp"

struct Functor_getCost {
    BinaryConstraint& obj;
    inline Functor_getCost(BinaryConstraint& in)
        : obj(in)
    {
    }
    inline Cost operator()(EnumeratedVariable* xx, EnumeratedVariable* yy, Value vx, Value vy) const;
};
struct Functor_getCostReverse {
    BinaryConstraint& obj;
    inline Functor_getCostReverse(BinaryConstraint& in)
        : obj(in)
    {
    }
    inline Cost operator()(EnumeratedVariable* xx, EnumeratedVariable* yy, Value vy, Value vx) const;
};

class BinaryConstraint : public AbstractBinaryConstraint<EnumeratedVariable, EnumeratedVariable> {
protected:
    unsigned int sizeX;
    unsigned int sizeY;
    vector<StoreCost> deltaCostsX;
    vector<StoreCost> deltaCostsY;
#ifdef NO_STORE_BINARY_COSTS
    vector<Cost> costs;
#else
    vector<StoreCost> costs;
#endif
    vector<Value> supportX;
    vector<Value> supportY;

    template <typename T>
    void findSupport(T getCost, EnumeratedVariable* x, EnumeratedVariable* y,
        vector<Value>& supportX, vector<StoreCost>& deltaCostsX);
    template <typename T>
    void findFullSupport(T getCost, EnumeratedVariable* x, EnumeratedVariable* y,
        vector<Value>& supportX, vector<StoreCost>& deltaCostsX,
        vector<Value>& supportY, vector<StoreCost>& deltaCostsY);
    template <typename T>
    void projection(T getCost, EnumeratedVariable* x, EnumeratedVariable* y, Value valueY, vector<StoreCost>& deltaCostsX);
    template <typename T>
    bool verify(T getCost, EnumeratedVariable* x, EnumeratedVariable* y);

    // return true if unary support of x is broken
    bool project(EnumeratedVariable* x, Value value, Cost cost, vector<StoreCost>& deltaCostsX);
    void extend(EnumeratedVariable* x, Value value, Cost cost, vector<StoreCost>& deltaCostsX);

    void findSupportX() { findSupport(Functor_getCost(*this), x, y, supportX, deltaCostsX); }
    void findSupportY() { findSupport(Functor_getCostReverse(*this), y, x, supportY, deltaCostsY); }
    void findFullSupportX() { findFullSupport(Functor_getCost(*this), x, y, supportX, deltaCostsX, supportY, deltaCostsY); }
    void findFullSupportY() { findFullSupport(Functor_getCostReverse(*this), y, x, supportY, deltaCostsY, supportX, deltaCostsX); }
    void projectX() { projection(Functor_getCost(*this), x, y, y->getValue(), deltaCostsX); }
    void projectY() { projection(Functor_getCostReverse(*this), y, x, x->getValue(), deltaCostsY); }
    bool verifyX() { return verify(Functor_getCost(*this), x, y); }
    bool verifyY() { return verify(Functor_getCostReverse(*this), y, x); }
    bool projectX(Value value, Cost cost) { return project(x, value, cost, deltaCostsX); }
    bool projectY(Value value, Cost cost) { return project(y, value, cost, deltaCostsY); }
    void extendX(Value value, Cost cost) { extend(x, value, cost, deltaCostsX); }
    void extendY(Value value, Cost cost) { extend(y, value, cost, deltaCostsY); }

public:
    vector<Cost> trwsM;

    BinaryConstraint(WCSP* wcsp, EnumeratedVariable* xx, EnumeratedVariable* yy, vector<Cost>& tab);

    BinaryConstraint(WCSP* wcsp);

    virtual ~BinaryConstraint() {}

    bool extension() const FINAL { return true; }
    bool isBinary() const FINAL { return true; }

    Cost getCost(Value vx, Value vy) const
    {
        unsigned int ix = x->toIndex(vx);
        unsigned int iy = y->toIndex(vy);
        Cost res = costs[ix * sizeY + iy];
        // if (res >= wcsp->getUb() || res - deltaCostsX[ix] - deltaCostsY[iy] + wcsp->getLb() >= wcsp->getUb()) return wcsp->getUb();
        res -= deltaCostsX[ix] + deltaCostsY[iy];
        assert(res >= MIN_COST);
        return res;
    }
    Cost getCostTRWS(Value vx, Value vy) const
    {
        unsigned int ix = x->toIndex(vx);
        unsigned int iy = y->toIndex(vy);
        return costs[ix * sizeY + iy] - (deltaCostsX[ix] + deltaCostsY[iy]);
    }
    void projectTRWS(EnumeratedVariable* var, Value value, Cost cost)
    {
        vector<StoreCost>& deltaCosts = (var == x) ? deltaCostsX : deltaCostsY;
        deltaCosts[var->toIndex(value)] += cost;
        var->project(value, cost, true);
    }
    Cost normalizeTRWS()
    {
        Cost minCost = numeric_limits<Cost>::max();
        for (EnumeratedVariable::iterator xIter = x->begin(); xIter != x->end(); ++xIter) {
            unsigned int ix = x->toIndex(*xIter);
            for (EnumeratedVariable::iterator yIter = y->begin(); yIter != y->end(); ++yIter) {
                unsigned int iy = y->toIndex(*yIter);
                minCost = min<Cost>(minCost, costs[ix * sizeY + iy] - (deltaCostsX[ix] + deltaCostsY[iy]));
            }
        }
        if (minCost != MIN_COST) {
            for (EnumeratedVariable::iterator xIter = x->begin(); xIter != x->end(); ++xIter) {
                unsigned int ix = x->toIndex(*xIter);
                for (EnumeratedVariable::iterator yIter = y->begin(); yIter != y->end(); ++yIter) {
                    unsigned int iy = y->toIndex(*yIter);
                    costs[ix * sizeY + iy] -= minCost;
                }
            }
        }
        return minCost;
    }

    vector<Cost> getDeltaCostsX()
    {
        vector<Cost> deltas;
        for (Cost cost : deltaCostsX) {
            deltas.push_back(cost);
        }
        return deltas;
    }
    vector<Cost> getDeltaCostsY()
    {
        vector<Cost> deltas;
        for (Cost cost : deltaCostsY) {
            deltas.push_back(cost);
        }
        return deltas;
    }
    void setDeltaCostsX(const vector<Cost>& update)
    {
        assert(update.size() == sizeX);
        for (unsigned int a = 0; a < sizeX; a++) {
            deltaCostsX[a] = update[a];
        }
    }
    void setDeltaCostsY(const vector<Cost>& update)
    {
        assert(update.size() == sizeY);
        for (unsigned int b = 0; b < sizeY; b++) {
            deltaCostsY[b] = update[b];
        }
    }
    Cost getCost(EnumeratedVariable* xx, EnumeratedVariable* yy, Value vx, Value vy)
    {
        unsigned int vindex[2];
        vindex[getIndex(xx)] = xx->toIndex(vx);
        vindex[getIndex(yy)] = yy->toIndex(vy);
        Cost res = costs[vindex[0] * sizeY + vindex[1]];
        // if (res >= wcsp->getUb() || res - deltaCostsX[vindex[0]] - deltaCostsY[vindex[1]] + wcsp->getLb() >= wcsp->getUb()) return wcsp->getUb();
        res -= deltaCostsX[vindex[0]] + deltaCostsY[vindex[1]];
        assert(res >= MIN_COST);
        return res;
    }

    Cost getCost() FINAL
    {
        int vX = x->getValue();
        int vY = y->getValue();
        return getCost(vX, vY);
    }

    void addcost(Value vx, Value vy, Cost mincost)
    {
#ifdef NO_STORE_BINARY_COSTS
        if (Store::getDepth() > 0) {
            cerr << "Cannot modify costs in binary cost functions during search!" << endl;
            throw BadConfiguration();
        }
#endif
        assert(ToulBar2::verbose < 4 || ((cout << "[" << Store::getDepth() << ",W" << wcsp->getIndex() << "] addcost(C" << getVar(0)->getName() << "," << getVar(1)->getName() << "," << vx << "," << vy << "), " << mincost << ")" << endl), true));
        assert(mincost >= MIN_COST || !LUBTEST(getCost(vx, vy), -mincost) || ToulBar2::isZ); // Warning! negative costs can be added temporally by variable elimination on the fly
        unsigned int ix = x->toIndex(vx);
        unsigned int iy = y->toIndex(vy);
        costs[ix * sizeY + iy] += mincost;
    }

    void addcost(EnumeratedVariable* xin, EnumeratedVariable* yin, Value vx, Value vy, Cost mincost)
    {
#ifdef NO_STORE_BINARY_COSTS
        if (Store::getDepth() > 0) {
            cerr << "Cannot modify costs in binary cost functions during search!" << endl;
            throw BadConfiguration();
        }
#endif
        assert(ToulBar2::verbose < 4 || ((cout << "[" << Store::getDepth() << ",W" << wcsp->getIndex() << "] addcost(C" << xin->getName() << "," << yin->getName() << "," << vx << "," << vy << "), " << mincost << ")" << endl), true));
        assert(mincost >= MIN_COST || !LUBTEST(getCost(xin, yin, vx, vy), -mincost));
        if (xin == x) {
            costs[x->toIndex(vx) * sizeY + y->toIndex(vy)] += mincost;
        } else {
            costs[x->toIndex(vy) * sizeY + y->toIndex(vx)] += mincost;
        }
    }

    void setCost(Cost c)
    {
#ifdef NO_STORE_BINARY_COSTS
        if (Store::getDepth() > 0) {
            cerr << "Cannot modify costs in binary cost functions during search!" << endl;
            throw BadConfiguration();
        }
#endif
        assert(std::all_of(deltaCostsX.begin(), deltaCostsX.end(), [](Cost c) { return c == MIN_COST; }));
        assert(std::all_of(deltaCostsY.begin(), deltaCostsY.end(), [](Cost c) { return c == MIN_COST; }));
        for (unsigned int a = 0; a < sizeX; a++)
            for (unsigned int b = 0; b < sizeY; b++)
                costs[a * sizeY + b] = c;
    }

    void setcost(EnumeratedVariable* xin, EnumeratedVariable* yin, Value vx, Value vy, Cost mincost)
    {
#ifdef NO_STORE_BINARY_COSTS
        if (Store::getDepth() > 0) {
            cerr << "Cannot modify costs in binary cost functions during search!" << endl;
            throw BadConfiguration();
        }
#endif
        assert(ToulBar2::verbose < 4 || ((cout << "[" << Store::getDepth() << ",W" << wcsp->getIndex() << "] setcost(C" << xin->getName() << "," << yin->getName() << "," << vx << "," << vy << "), " << mincost << ")" << endl), true));
        assert(std::all_of(deltaCostsX.begin(), deltaCostsX.end(), [](Cost c) { return c == MIN_COST; }));
        assert(std::all_of(deltaCostsY.begin(), deltaCostsY.end(), [](Cost c) { return c == MIN_COST; }));
        if (xin == x)
            costs[x->toIndex(vx) * sizeY + y->toIndex(vy)] = mincost;
        else
            costs[x->toIndex(vy) * sizeY + y->toIndex(vx)] = mincost;
    }

    void setcost(Value vx, Value vy, Cost mincost)
    {
#ifdef NO_STORE_BINARY_COSTS
        if (Store::getDepth() > 0) {
            cerr << "Cannot modify costs in binary cost functions during search!" << endl;
            throw BadConfiguration();
        }
#endif
        assert(ToulBar2::verbose < 4 || ((cout << "[" << Store::getDepth() << ",W" << wcsp->getIndex() << "] setcost(C" << getVar(0)->getName() << "," << getVar(1)->getName() << "," << vx << "," << vy << "), " << mincost << ")" << endl), true));
        assert(std::all_of(deltaCostsX.begin(), deltaCostsX.end(), [](Cost c) { return c == MIN_COST; }));
        assert(std::all_of(deltaCostsY.begin(), deltaCostsY.end(), [](Cost c) { return c == MIN_COST; }));
        costs[x->toIndex(vx) * sizeY + y->toIndex(vy)] = mincost;
    }

    void addCosts(EnumeratedVariable* xin, EnumeratedVariable* yin, vector<Cost>& costsin)
    {
#ifdef NO_STORE_BINARY_COSTS
        if (Store::getDepth() > 0) {
            cerr << "Cannot modify costs in binary cost functions during search!" << endl;
            throw BadConfiguration();
        }
#endif
        assert(ToulBar2::verbose < 4 || ((cout << "[" << Store::getDepth() << ",W" << wcsp->getIndex() << "] add binary cost vector to (C" << getVar(0)->getName() << "," << getVar(1)->getName() << ") " << costsin[0] << "," << costsin[1] << "," << costsin[2] << "," << costsin[3] << " ..." << endl), true));
        assert(costsin.size() <= costs.size());
        unsigned int ix, iy;
        for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
            for (EnumeratedVariable::iterator itery = y->begin(); itery != y->end(); ++itery) {
                ix = x->toIndex(*iterx);
                iy = y->toIndex(*itery);
                if (xin == x)
                    costs[ix * sizeY + iy] += costsin[ix * sizeY + iy];
                else
                    costs[ix * sizeY + iy] += costsin[iy * sizeX + ix];
            }
        }
    }

    void addCosts(BinaryConstraint* xy)
    {
#ifdef NO_STORE_BINARY_COSTS
        if (Store::getDepth() > 0) {
            cerr << "Cannot modify costs in binary cost functions during search!" << endl;
            throw BadConfiguration();
        }
#endif
        assert(ToulBar2::verbose < 4 || ((cout << "[" << Store::getDepth() << ",W" << wcsp->getIndex() << "] add binary cost function to C(" << getVar(0)->getName() << "," << getVar(1)->getName() << ")" << endl), true));
        assert(((x == xy->x) && (y == xy->y)) || ((x == xy->y) && (y == xy->x)));
        unsigned int ix, iy;
        for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
            for (EnumeratedVariable::iterator itery = y->begin(); itery != y->end(); ++itery) {
                ix = x->toIndex(*iterx);
                iy = y->toIndex(*itery);
                Cost c = costs[ix * sizeY + iy];
                // if(costs[ix * sizeY + iy] < wcsp->getUb()) //BUG with BTD: ub is only local, deltaCosts should be considered
                {
                    costs[ix * sizeY + iy] = c + xy->getCost(x, y, *iterx, *itery);
                }
            }
        }
    }

    void clearCosts()
    {
#ifdef NO_STORE_BINARY_COSTS
        if (Store::getDepth() > 0) {
            cerr << "Cannot modify costs in binary cost functions during search!" << endl;
            throw BadConfiguration();
        }
#endif
        assert(ToulBar2::verbose < 4 || ((cout << "[" << Store::getDepth() << ",W" << wcsp->getIndex() << "] clearcosts(C" << getVar(0)->getName() << "," << getVar(1)->getName() << ")" << endl), true));
        for (unsigned int i = 0; i < sizeX; i++)
            deltaCostsX[i] = MIN_COST;
        for (unsigned int j = 0; j < sizeY; j++)
            deltaCostsY[j] = MIN_COST;
        for (unsigned int i = 0; i < sizeX; i++) {
            for (unsigned int j = 0; j < sizeY; j++) {
                costs[i * sizeY + j] = MIN_COST;
            }
        }
    }

    void clearFiniteCosts()
    {
#ifdef NO_STORE_BINARY_COSTS
        if (Store::getDepth() > 0) {
            cerr << "Cannot modify finite costs in binary cost functions during search!" << endl;
            throw BadConfiguration();
        }
#endif
        assert(ToulBar2::verbose < 4 || ((cout << "[" << Store::getDepth() << ",W" << wcsp->getIndex() << "] clearfinitecosts(C" << getVar(0)->getName() << "," << getVar(1)->getName() << ")" << endl), true));
        for (unsigned int i = 0; i < sizeX; i++)
            deltaCostsX[i] = MIN_COST;
        for (unsigned int j = 0; j < sizeY; j++)
            deltaCostsY[j] = MIN_COST;
        for (unsigned int i = 0; i < sizeX; i++) {
            for (unsigned int j = 0; j < sizeY; j++) {
                if (!CUT(costs[i * sizeY + j], wcsp->getUb())) {
                    costs[i * sizeY + j] = MIN_COST;
                }
            }
        }
    }

    Cost getMaxFiniteCost()
    {
        Cost ub = wcsp->getUb();
        Cost maxcost = MIN_COST;
        for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
            for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                Cost cost = getCost(*iterX, *iterY);
                if (cost < ub && cost > maxcost)
                    maxcost = cost;
            }
        }
        return maxcost;
    }

    void setInfiniteCost(Cost ub)
    {
        bool modified = false;
        Cost mult_ub = ((wcsp->getUb() < (MAX_COST / MEDIUM_COST)) ? (max(LARGE_COST, wcsp->getUb() * MEDIUM_COST)) : wcsp->getUb());
        for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
            unsigned int ix = x->toIndex(*iterx);
            for (EnumeratedVariable::iterator itery = y->begin(); itery != y->end(); ++itery) {
                unsigned int iy = y->toIndex(*itery);
                Cost cost = costs[ix * sizeY + iy] + x->getCost(*iterx) + y->getCost(*itery);
                Cost delta = deltaCostsX[ix] + deltaCostsY[iy];
                if (CUT(cost - delta, ub)) {
                    costs[ix * sizeY + iy] = mult_ub + delta;
                    modified = true;
                }
            }
        }
        if (modified) {
            propagate();
        }
    }

    Cost evalsubstr(const Tuple& s, Constraint* ctr) FINAL
    {
        Value vals[2];
        int count = 0;

        for (int i = 0; i < 2; i++) {
            EnumeratedVariable* var = (EnumeratedVariable*)getVar(i);
            int ind = ctr->getIndex(var);
            if (ind >= 0) {
                vals[i] = var->toValue(s[ind]);
                count++;
            }
        }
        if (count == 2)
            return getCost(vals[0], vals[1]);
        else
            return MIN_COST;
    }
    Cost evalsubstr(const Tuple& s, NaryConstraint* ctr) FINAL { return evalsubstr(s, (Constraint*)ctr); } // NaryConstraint class undefined

    void resetSupports()
    {
        supportX.assign(supportX.size(), y->getInf());
        supportY.assign(supportY.size(), x->getInf());
    }

    Value getSupport(EnumeratedVariable* var, Value v)
    {
        if (var == x)
            return supportX[x->toIndex(v)];
        else
            return supportY[y->toIndex(v)];
    }

    void setSupport(EnumeratedVariable* var, Value v, Value s)
    {
        if (var == x)
            supportX[x->toIndex(v)] = s;
        else
            supportY[y->toIndex(v)] = s;
    }

    EnumeratedVariable* xvar;
    EnumeratedVariable* yvar;
    EnumeratedVariable::iterator itvx;
    EnumeratedVariable::iterator itvy;

    void first()
    {
        itvx = x->begin();
        itvy = y->begin();
        xvar = x;
        yvar = y;
    }

    bool next(Tuple& t, Cost& c)
    {
        Tuple tch(2, 0);
        if (itvx != xvar->end()) {
            unsigned int ix = xvar->toIndex(*itvx);
            tch[0] = ix;
            if (itvy != yvar->end()) {
                unsigned int iy = yvar->toIndex(*itvy);
                tch[1] = iy;
                t = tch;
                c = getCost(xvar, yvar, *itvx, *itvy);
                ++itvy;
                return true;
            } else {
                ++itvx;
                itvy = yvar->begin();
                return next(t, c);
            }
        }
        return false;
    }

    void firstlex() { first(); }
    bool nextlex(Tuple& t, Cost& c) { return next(t, c); }

    void setTuple(const Tuple& t, Cost c) FINAL
    {
        Value v0 = x->toValue(t[0]);
        Value v1 = y->toValue(t[1]);
        Cost oldc = getCost(v0, v1);
        addcost(v0, v1, c - oldc);
    }

    void addtoTuple(const Tuple& t, Cost c) FINAL
    {
        Value v0 = x->toValue(t[0]);
        Value v1 = y->toValue(t[1]);
        addcost(v0, v1, c);
    }

    void fillElimConstr(EnumeratedVariable* xin, EnumeratedVariable* yin, Constraint* from1, Constraint* from2)
    {
        x = xin;
        y = yin;
        sizeX = x->getDomainInitSize();
        sizeY = y->getDomainInitSize();
        if (sizeX > deltaCostsX.size())
            deltaCostsX.resize(sizeX, StoreCost(MIN_COST));
        if (sizeY > deltaCostsY.size())
            deltaCostsY.resize(sizeY, StoreCost(MIN_COST));
        if (sizeX > supportX.size())
            supportX.resize(sizeX);
        if (sizeY > supportY.size())
            supportY.resize(sizeY);
        if (max(sizeX, sizeY) > trwsM.size())
            trwsM.resize(max(sizeX, sizeY), MIN_COST);
        if ((size_t)sizeX * (size_t)sizeY > costs.size())
#ifdef NO_STORE_BINARY_COSTS
            costs.resize((size_t)sizeX * (size_t)sizeY, MIN_COST);
#else
            costs.resize((size_t)sizeX * (size_t)sizeY, StoreCost(MIN_COST));
#endif
        linkX->removed = true;
        linkY->removed = true;
        linkX->content.constr = this;
        linkY->content.constr = this;
        linkX->content.scopeIndex = 0;
        linkY->content.scopeIndex = 1;
        setDACScopeIndex();
        resetConflictWeight(); // warning! it is done before having updated costs!!
        elimFrom(from1, from2);
    }

    bool project(int varIndex, Value value, Cost cost)
    {
        if (varIndex == 0)
            return projectX(value, cost);
        else
            return projectY(value, cost);
    }
    void extend(int varIndex, Value value, Cost cost)
    {
        if (varIndex == 0)
            return extendX(value, cost);
        else
            return extendY(value, cost);
    }

    void propagate()
    {
        if (ToulBar2::dumpWCSP % 2) // do not propagate if problem is dumped before preprocessing
            return;
        if (x->assigned()) {
            assign(0);
            return;
        }
        if (y->assigned()) {
            assign(1);
            return;
        }
        if (getDACScopeIndex() == 0) {
            x->queueAC();
            x->queueEAC1();
            if (ToulBar2::LcLevel >= LC_DAC)
                y->queueDAC();
            else
                y->queueAC();
        } else {
            y->queueAC();
            y->queueEAC1();
            if (ToulBar2::LcLevel >= LC_DAC)
                x->queueDAC();
            else
                x->queueAC();
        }
        if (ToulBar2::FullEAC) {
            reviseEACGreedySolution();
            for (ConstraintList::iterator iter = x->getConstrs()->begin(); iter != x->getConstrs()->end(); ++iter) {
                if ((*iter).constr->isTernary() && (*iter).constr->getIndex(y) >= 0) {
                    (*iter).constr->reviseEACGreedySolution(); // all ternary cost functions including this current binary cost function need to be revised
                }
            }
        }
    }
    void remove(int varIndex)
    {
        if (varIndex == 0)
            y->queueDEE();
        else
            x->queueDEE();
        if (ToulBar2::LcLevel == LC_AC) {
            if (varIndex == 0)
                findSupportY();
            else
                findSupportX();
        } else {
            if (getDACScopeIndex() == 0) {
                if (varIndex == 0)
                    findSupportY();
            } else {
                if (varIndex == 1)
                    findSupportX();
            }
        }
    }
    void projectFromZero(int varIndex)
    {
        if (getDACScopeIndex() == 0) {
            if (varIndex == 1)
                findFullSupportX();
        } else {
            if (varIndex == 0)
                findFullSupportY();
        }
    }
    // Trick! instead of doing remove(index) now, let AC queue do the job.
    // So several incdec events on the same constraint can be merged into one AC event
    void increase(int index)
    {
        if (index == 0)
            x->queueAC();
        else
            y->queueAC();
    }
    void decrease(int index)
    {
        if (index == 0)
            x->queueAC();
        else
            y->queueAC();
    } // Trick! instead of remove(index);
    void assign(int varIndex)
    {
        deconnect(); // Warning! deconnection has to be done before the projection
        if (varIndex == 0) {
            projectY();
        } else {
            projectX();
        }
    }

    bool checkEACGreedySolution(int index = -1, Value a = 0) FINAL
    {
        assert(x->canbe((getIndex(x) == index) ? a : x->getSupport()));
        assert(x->getCost((getIndex(x) == index) ? a : x->getSupport()) == MIN_COST);
        assert(y->canbe((getIndex(y) == index) ? a : y->getSupport()));
        assert(y->getCost((getIndex(y) == index) ? a : y->getSupport()) == MIN_COST);
        return (getCost((getIndex(x) == index) ? a : x->getSupport(), (getIndex(y) == index) ? a : y->getSupport()) == MIN_COST);
    }

    bool reviseEACGreedySolution(int index = -1, Value a = 0) FINAL
    {
        bool result = checkEACGreedySolution(index, a);
        if (!result) {
            if (index >= 0)
                getVar(index)->unsetFullEAC();
            else {
                x->unsetFullEAC();
                y->unsetFullEAC();
            }
        }
        return result;
    }

    void fillEAC2(int varIndex)
    {
        assert(!isDuplicate());
        if (getDACScopeIndex() == 0) {
            if (varIndex == 0) {
                assert(y->canbe(y->getSupport()));
                unsigned int yindex = y->toIndex(y->getSupport());
                if (x->cannotbe(supportY[yindex]) || x->getCost(supportY[yindex]) > MIN_COST || getCost(supportY[yindex], y->getSupport()) > MIN_COST || (ToulBar2::vacValueHeuristic && Store::getDepth() < abs(ToulBar2::vac))) {
                    y->queueEAC2();
                }
            }
        } else {
            if (varIndex == 1) {
                assert(x->canbe(x->getSupport()));
                unsigned int xindex = x->toIndex(x->getSupport());
                if (y->cannotbe(supportX[xindex]) || y->getCost(supportX[xindex]) > MIN_COST || getCost(x->getSupport(), supportX[xindex]) > MIN_COST || (ToulBar2::vacValueHeuristic && Store::getDepth() < abs(ToulBar2::vac))) {
                    x->queueEAC2();
                }
            }
        }
    }

    bool isEAC(int varIndex, Value a)
    {
        assert(!isDuplicate());
        if (ToulBar2::QueueComplexity && varIndex == getDACScopeIndex() && !ToulBar2::FullEAC)
            return true;
        if (varIndex == 0) {
            assert(y->canbe(y->getSupport()));
            assert(y->getCost(y->getSupport()) == MIN_COST);
            if (getCost(a, y->getSupport()) > MIN_COST) {
                if (ToulBar2::FullEAC)
                    x->unsetFullEAC();
                unsigned int xindex = x->toIndex(a);
                if (y->cannotbe(supportX[xindex]) || y->getCost(supportX[xindex]) > MIN_COST || getCost(a, supportX[xindex]) > MIN_COST) {
                    for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                        if (y->getCost(*iterY) == MIN_COST && getCost(a, *iterY) == MIN_COST) {
                            supportX[xindex] = *iterY;
                            return true;
                        }
                    }
                    return false;
                }
            }
        } else {
            assert(x->canbe(x->getSupport()));
            assert(x->getCost(x->getSupport()) == MIN_COST);
            if (getCost(x->getSupport(), a) > MIN_COST) {
                if (ToulBar2::FullEAC)
                    y->unsetFullEAC();
                unsigned int yindex = y->toIndex(a);
                if (x->cannotbe(supportY[yindex]) || x->getCost(supportY[yindex]) > MIN_COST || getCost(supportY[yindex], a) > MIN_COST) {
                    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
                        if (x->getCost(*iterX) == MIN_COST && getCost(*iterX, a) == MIN_COST) {
                            supportY[yindex] = *iterX;
                            return true;
                        }
                    }
                    return false;
                }
            }
        }
        return true;
    }

    void findFullSupportEAC(int varIndex)
    {
        assert(!isDuplicate());
        if (ToulBar2::QueueComplexity && varIndex == getDACScopeIndex() && !ToulBar2::FullEAC)
            return;
        if (varIndex == 0)
            findFullSupportX();
        else
            findFullSupportY();
    }

    bool verify()
    {
        if (ToulBar2::LcLevel == LC_DAC) {
            if (getDACScopeIndex() == 0)
                return verifyX();
            else
                return verifyY();
        } else {
            return verifyX() && verifyY();
        }
    }

    pair<pair<Cost, Cost>, pair<Cost, Cost>> getMaxCost(int varIndex, Value a, Value b);

    double computeTightness()
    {
        int count = 0;
        double sum = 0;
        double tight = -1;
        Cost* costs = new Cost[x->getDomainSize() * y->getDomainSize()];
        for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
            for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                Cost c = getCost(*iterX, *iterY);
                sum += to_double(min(wcsp->getUb(), c));
                costs[count] = min(wcsp->getUb(), c);
                count++;
            }
        }
        if (ToulBar2::weightedTightness == 2) {
            tight = to_double(stochastic_selection<Cost>(costs, 0, count - 1, count / 2));
        } else {
            tight = sum / (double)count;
        }
        delete[] costs;
        return tight;
    }

    EnumeratedVariable* commonVar(BinaryConstraint* bctr)
    {
        if (getIndex(bctr->getVar(0)) >= 0)
            return (EnumeratedVariable*)bctr->getVar(0);
        else if (getIndex(bctr->getVar(1)) >= 0)
            return (EnumeratedVariable*)bctr->getVar(1);
        else
            return NULL;
    }

    void permute(EnumeratedVariable* xin, Value a, Value b);

    bool isFunctional(EnumeratedVariable* x, EnumeratedVariable* y, map<Value, Value>& functional);

    virtual void print(ostream& os);
    void dump(ostream& os, bool original = true);
    void dump_CFN(ostream& os, bool original = true);
    Long size() const FINAL { return (Long)sizeX * sizeY; }
#ifdef NO_STORE_BINARY_COSTS
    Long space() const FINAL
    {
        return (Long)sizeof(Cost) * sizeX * sizeY;
    }
#else
    Long space() const FINAL
    {
        return (Long)sizeof(StoreCost) * sizeX * sizeY;
    }
#endif

    friend struct Functor_getCost;
    friend struct Functor_getCostReverse;
};

inline Cost Functor_getCost::operator()(EnumeratedVariable* xx, EnumeratedVariable* yy, Value vx, Value vy) const
{
    assert(xx == obj.x);
    assert(yy == obj.y);
    return obj.getCost(vx, vy);
}
inline Cost Functor_getCostReverse::operator()(EnumeratedVariable* yy, EnumeratedVariable* xx, Value vy, Value vx) const
{
    assert(xx == obj.x);
    assert(yy == obj.y);
    return obj.getCost(vx, vy);
}

template <typename T>
void BinaryConstraint::findSupport(T getCost, EnumeratedVariable* x, EnumeratedVariable* y,
    vector<Value>& supportX, vector<StoreCost>& deltaCostsX)
{
    assert(connected());
    wcsp->revise(this);
    if (ToulBar2::verbose >= 3)
        cout << "[" << Store::getDepth() << ",W" << wcsp->getIndex() << "] findSupport C" << x->getName() << "," << y->getName() << endl;
    bool supportBroken = false;
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        unsigned int xindex = x->toIndex(*iterX);
        Value support = supportX[xindex];
        if (y->cannotbe(support) || getCost(x, y, *iterX, support) > MIN_COST) {
            Value minCostValue = y->getInf();
            Cost minCost = getCost(x, y, *iterX, minCostValue);
            EnumeratedVariable::iterator iterY = y->begin();
            for (++iterY; minCost > MIN_COST && iterY != y->end(); ++iterY) {
                Cost cost = getCost(x, y, *iterX, *iterY);
                if (GLB(&minCost, cost)) {
                    minCostValue = *iterY;
                }
            }
            if (minCost > MIN_COST) {
                supportBroken |= project(x, *iterX, minCost, deltaCostsX);
                if (deconnected())
                    return;
            }
            supportX[xindex] = minCostValue;
        }
    }
    if (supportBroken) {
        x->findSupport();
    }
}

template <typename T>
void BinaryConstraint::findFullSupport(T getCost, EnumeratedVariable* x, EnumeratedVariable* y,
    vector<Value>& supportX, vector<StoreCost>& deltaCostsX,
    vector<Value>& supportY, vector<StoreCost>& deltaCostsY)
{
    assert(connected());
    wcsp->revise(this);
    if (ToulBar2::verbose >= 3)
        cout << "[" << Store::getDepth() << ",W" << wcsp->getIndex() << "] findFullSupport C" << x->getName() << "," << y->getName() << endl;
    bool supportBroken = false;
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        unsigned int xindex = x->toIndex(*iterX);
        Value support = supportX[xindex];
        if (y->cannotbe(support) || getCost(x, y, *iterX, support) + y->getCost(support) > MIN_COST) {
            Value minCostValue = y->getInf();
            Cost minCost = getCost(x, y, *iterX, minCostValue) + y->getCost(minCostValue);
            EnumeratedVariable::iterator iterY = y->begin();
            for (++iterY; minCost > MIN_COST && iterY != y->end(); ++iterY) {
                Cost cost = getCost(x, y, *iterX, *iterY) + y->getCost(*iterY);
                if (GLB(&minCost, cost)) {
                    minCostValue = *iterY;
                }
            }
            if (minCost > MIN_COST) {
                // extend unary to binary
                for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                    Cost cost = getCost(x, y, *iterX, *iterY);
                    if (GLBTEST(minCost, cost)) {
                        extend(y, *iterY, minCost - cost, deltaCostsY);
                        supportY[y->toIndex(*iterY)] = *iterX;
                        //                         if (ToulBar2::vac) {
                        //                             x->queueVAC2();
                        //                             y->queueVAC2();
                        //                         }
                    }
                }
                supportBroken |= project(x, *iterX, minCost, deltaCostsX);
                if (deconnected())
                    return;
            }
            supportX[xindex] = minCostValue;
        }
    }
    if (supportBroken) {
        x->findSupport();
    }
}

template <typename T>
void BinaryConstraint::projection(T getCost, EnumeratedVariable* x, EnumeratedVariable* y, Value valueY, vector<StoreCost>& deltaCostsX)
{
    x->queueDEE();
    bool supportBroken = false;
    wcsp->revise(this);
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        Cost cost = getCost(x, y, *iterX, valueY);
        if (cost > MIN_COST) {
            supportBroken |= project(x, *iterX, cost, deltaCostsX);
        }
    }
    if (supportBroken) {
        x->findSupport();
    }
}

template <typename T>
bool BinaryConstraint::verify(T getCost, EnumeratedVariable* x, EnumeratedVariable* y)
{
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        Cost minCost = getCost(x, y, *iterX, y->getInf());
        if (ToulBar2::LcLevel >= LC_DAC && getDACScopeIndex() == getIndex(x))
            minCost += y->getCost(y->getInf());
        EnumeratedVariable::iterator iterY = y->begin();
        for (++iterY; minCost > MIN_COST && iterY != y->end(); ++iterY) {
            Cost cost = getCost(x, y, *iterX, *iterY);
            if (ToulBar2::LcLevel >= LC_DAC && getDACScopeIndex() == getIndex(x))
                cost += y->getCost(*iterY);
            GLB(&minCost, cost);
        }
        if (minCost > MIN_COST) {
            cout << *this;
            return false;
        }
    }
    return true;
}

#endif /*TB2BINCONSTR_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

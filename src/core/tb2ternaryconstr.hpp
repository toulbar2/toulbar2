/** \file tb2ternaryconstr.hpp
 *  \brief Ternary constraint applied on variables with enumerated domains.
 *  Functional variables are automatically detected and used to reduce time complexity of propagation by a factor of its domain size (isEAC, findSupport and findFullSupport).
 *  If the first variable in the scope is functional wrt to the two others then it is also used to reduce space complexity (quadratic instead of cubic)
 *  \warning EAC is not applied on duplicated ternary constraints
 */

#ifndef TB2TERNARYCONSTR_HPP_
#define TB2TERNARYCONSTR_HPP_

#include "tb2abstractconstr.hpp"
#include "tb2enumvar.hpp"
#include "tb2binconstr.hpp"

struct Functor_getCostXYZ {
    TernaryConstraint& obj;
    inline Functor_getCostXYZ(TernaryConstraint& in)
        : obj(in)
    {
    }
    inline Cost operator()(EnumeratedVariable* xx, EnumeratedVariable* yy, EnumeratedVariable* zz, Value vx, Value vy, Value vz) const;
};
struct Functor_getCostXZY {
    TernaryConstraint& obj;
    inline Functor_getCostXZY(TernaryConstraint& in)
        : obj(in)
    {
    }
    inline Cost operator()(EnumeratedVariable* xx, EnumeratedVariable* yy, EnumeratedVariable* zz, Value vx, Value vy, Value vz) const;
};
struct Functor_getCostYXZ {
    TernaryConstraint& obj;
    inline Functor_getCostYXZ(TernaryConstraint& in)
        : obj(in)
    {
    }
    inline Cost operator()(EnumeratedVariable* xx, EnumeratedVariable* yy, EnumeratedVariable* zz, Value vx, Value vy, Value vz) const;
};
struct Functor_getCostYZX {
    TernaryConstraint& obj;
    inline Functor_getCostYZX(TernaryConstraint& in)
        : obj(in)
    {
    }
    inline Cost operator()(EnumeratedVariable* xx, EnumeratedVariable* yy, EnumeratedVariable* zz, Value vx, Value vy, Value vz) const;
};
struct Functor_getCostZXY {
    TernaryConstraint& obj;
    inline Functor_getCostZXY(TernaryConstraint& in)
        : obj(in)
    {
    }
    inline Cost operator()(EnumeratedVariable* xx, EnumeratedVariable* yy, EnumeratedVariable* zz, Value vx, Value vy, Value vz) const;
};
struct Functor_getCostZYX {
    TernaryConstraint& obj;
    inline Functor_getCostZYX(TernaryConstraint& in)
        : obj(in)
    {
    }
    inline Cost operator()(EnumeratedVariable* xx, EnumeratedVariable* yy, EnumeratedVariable* zz, Value vx, Value vy, Value vz) const;
};

struct Functor_getCostWithBinariesXYZ {
    TernaryConstraint& obj;
    inline Functor_getCostWithBinariesXYZ(TernaryConstraint& in)
        : obj(in)
    {
    }
    inline Cost operator()(EnumeratedVariable* xx, EnumeratedVariable* yy, EnumeratedVariable* zz, Value vx, Value vy, Value vz) const;
};
struct Functor_getCostWithBinariesXZY {
    TernaryConstraint& obj;
    inline Functor_getCostWithBinariesXZY(TernaryConstraint& in)
        : obj(in)
    {
    }
    inline Cost operator()(EnumeratedVariable* xx, EnumeratedVariable* yy, EnumeratedVariable* zz, Value vx, Value vy, Value vz) const;
};
struct Functor_getCostWithBinariesYXZ {
    TernaryConstraint& obj;
    inline Functor_getCostWithBinariesYXZ(TernaryConstraint& in)
        : obj(in)
    {
    }
    inline Cost operator()(EnumeratedVariable* xx, EnumeratedVariable* yy, EnumeratedVariable* zz, Value vx, Value vy, Value vz) const;
};
struct Functor_getCostWithBinariesYZX {
    TernaryConstraint& obj;
    inline Functor_getCostWithBinariesYZX(TernaryConstraint& in)
        : obj(in)
    {
    }
    inline Cost operator()(EnumeratedVariable* xx, EnumeratedVariable* yy, EnumeratedVariable* zz, Value vx, Value vy, Value vz) const;
};
struct Functor_getCostWithBinariesZXY {
    TernaryConstraint& obj;
    inline Functor_getCostWithBinariesZXY(TernaryConstraint& in)
        : obj(in)
    {
    }
    inline Cost operator()(EnumeratedVariable* xx, EnumeratedVariable* yy, EnumeratedVariable* zz, Value vx, Value vy, Value vz) const;
};
struct Functor_getCostWithBinariesZYX {
    TernaryConstraint& obj;
    inline Functor_getCostWithBinariesZYX(TernaryConstraint& in)
        : obj(in)
    {
    }
    inline Cost operator()(EnumeratedVariable* xx, EnumeratedVariable* yy, EnumeratedVariable* zz, Value vx, Value vy, Value vz) const;
};

struct Functor_addCostXYZ {
    TernaryConstraint& obj;
    inline Functor_addCostXYZ(TernaryConstraint& in)
        : obj(in)
    {
    }
    inline void operator()(EnumeratedVariable* xx, EnumeratedVariable* yy, EnumeratedVariable* zz, Value vx, Value vy, Value vz, Cost c);
};
struct Functor_addCostXZY {
    TernaryConstraint& obj;
    inline Functor_addCostXZY(TernaryConstraint& in)
        : obj(in)
    {
    }
    inline void operator()(EnumeratedVariable* xx, EnumeratedVariable* yy, EnumeratedVariable* zz, Value vx, Value vy, Value vz, Cost c);
};
struct Functor_addCostYXZ {
    TernaryConstraint& obj;
    inline Functor_addCostYXZ(TernaryConstraint& in)
        : obj(in)
    {
    }
    inline void operator()(EnumeratedVariable* xx, EnumeratedVariable* yy, EnumeratedVariable* zz, Value vx, Value vy, Value vz, Cost c);
};
struct Functor_addCostYZX {
    TernaryConstraint& obj;
    inline Functor_addCostYZX(TernaryConstraint& in)
        : obj(in)
    {
    }
    inline void operator()(EnumeratedVariable* xx, EnumeratedVariable* yy, EnumeratedVariable* zz, Value vx, Value vy, Value vz, Cost c);
};
struct Functor_addCostZXY {
    TernaryConstraint& obj;
    inline Functor_addCostZXY(TernaryConstraint& in)
        : obj(in)
    {
    }
    inline void operator()(EnumeratedVariable* xx, EnumeratedVariable* yy, EnumeratedVariable* zz, Value vx, Value vy, Value vz, Cost c);
};
struct Functor_addCostZYX {
    TernaryConstraint& obj;
    inline Functor_addCostZYX(TernaryConstraint& in)
        : obj(in)
    {
    }
    inline void operator()(EnumeratedVariable* xx, EnumeratedVariable* yy, EnumeratedVariable* zz, Value vx, Value vy, Value vz, Cost c);
};

struct Functor_getFunctionXYZ {
    TernaryConstraint& obj;
    inline Functor_getFunctionXYZ(TernaryConstraint& in)
        : obj(in)
    {
    }
    inline Value operator()(EnumeratedVariable* xx, EnumeratedVariable* yy, Value vx, Value vy) const;
};
struct Functor_getFunctionXZY {
    TernaryConstraint& obj;
    inline Functor_getFunctionXZY(TernaryConstraint& in)
        : obj(in)
    {
    }
    inline Value operator()(EnumeratedVariable* xx, EnumeratedVariable* yy, Value vx, Value vy) const;
};
struct Functor_getFunctionYXZ {
    TernaryConstraint& obj;
    inline Functor_getFunctionYXZ(TernaryConstraint& in)
        : obj(in)
    {
    }
    inline Value operator()(EnumeratedVariable* xx, EnumeratedVariable* yy, Value vx, Value vy) const;
};
struct Functor_getFunctionYZX {
    TernaryConstraint& obj;
    inline Functor_getFunctionYZX(TernaryConstraint& in)
        : obj(in)
    {
    }
    inline Value operator()(EnumeratedVariable* xx, EnumeratedVariable* yy, Value vx, Value vy) const;
};
struct Functor_getFunctionZXY {
    TernaryConstraint& obj;
    inline Functor_getFunctionZXY(TernaryConstraint& in)
        : obj(in)
    {
    }
    inline Value operator()(EnumeratedVariable* xx, EnumeratedVariable* yy, Value vx, Value vy) const;
};
struct Functor_getFunctionZYX {
    TernaryConstraint& obj;
    inline Functor_getFunctionZYX(TernaryConstraint& in)
        : obj(in)
    {
    }
    inline Value operator()(EnumeratedVariable* xx, EnumeratedVariable* yy, Value vx, Value vy) const;
};

class TernaryConstraint : public AbstractTernaryConstraint<EnumeratedVariable, EnumeratedVariable, EnumeratedVariable> {
protected:
    unsigned int sizeX;
    unsigned int sizeY;
    unsigned int sizeZ;
#ifdef NO_STORE_TERNARY_COSTS
    vector<Cost> costs;
#else
    vector<StoreCost> costs;
#endif
    vector<StoreCost> deltaCostsX;
    vector<StoreCost> deltaCostsY;
    vector<StoreCost> deltaCostsZ;
    Cost top;
    bool functionalX;
    vector<Value> functionX;
    bool functionalY;
    vector<Value> functionY;
    bool functionalZ;
    vector<Value> functionZ;
    vector<pair<Value, Value>> supportX;
    vector<pair<Value, Value>> supportY;
    vector<pair<Value, Value>> supportZ;
#ifdef NO_STORE_TERNARY_COSTS
    vector<Cost> costsYZ;
#else
    vector<StoreCost> costsYZ;
#endif
    inline Value getFunctionX(Value vy, Value vz) const
    {
        return functionX[y->toIndex(vy) * sizeZ + z->toIndex(vz)];
    }
    inline Value getFunctionY(Value vx, Value vz) const { return functionY[x->toIndex(vx) * sizeZ + z->toIndex(vz)]; }
    inline Value getFunctionZ(Value vx, Value vy) const { return functionZ[x->toIndex(vx) * sizeY + y->toIndex(vy)]; }

    // return true if unary support of x is broken
    bool project(EnumeratedVariable* x, Value value, Cost cost, vector<StoreCost>& deltaCostsX);
    void extend(EnumeratedVariable* x, Value value, Cost cost, vector<StoreCost>& deltaCostsX);

    template <typename T1, typename T2, typename T3>
    void project(T1 getCost, T2 addCost, bool functionalZ, T3 getFunctionZ, BinaryConstraint* xy, EnumeratedVariable* x, EnumeratedVariable* y, EnumeratedVariable* z, Value valx, Value valy, Cost cost);
    template <typename T1, typename T2, typename T3>
    void extend(T1 getCost, T2 addCost, bool functionalZ, T3 getFunctionZ, BinaryConstraint* xy, EnumeratedVariable* x, EnumeratedVariable* y, EnumeratedVariable* z, Value valx, Value valy, Cost cost);

    template <typename T1, typename T2, typename T3>
    void findSupport(T1 getCost, bool functionalY, T2 getFunctionY, bool functionalZ, T3 getFunctionZ,
        EnumeratedVariable* x, EnumeratedVariable* y, EnumeratedVariable* z,
        int getIndexX, int getIndexY, int getIndexZ,
        vector<pair<Value, Value>>& supportX, vector<StoreCost>& deltaCostsX,
        vector<pair<Value, Value>>& supportY, vector<pair<Value, Value>>& supportZ);
    template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10>
    void findFullSupport(T1 getCost, T2 getCostXZY, T3 getCostYZX, T4 getCostWithBinaries, T5 addCost, T6 addCostXZY, T7 addCostYZX, bool functionalX, T8 getFunctionX, bool functionalY, T9 getFunctionY, bool functionalZ, T10 getFunctionZ,
        EnumeratedVariable* x, EnumeratedVariable* y, EnumeratedVariable* z,
        int getIndexX, int getIndexY, int getIndexZ,
        vector<pair<Value, Value>>& supportX, vector<StoreCost>& deltaCostsX,
        vector<pair<Value, Value>>& supportY, vector<StoreCost>& deltaCostsY,
        vector<pair<Value, Value>>& supportZ, vector<StoreCost>& deltaCostsZ,
        BinaryConstraint* xy, BinaryConstraint* xz, BinaryConstraint* yz);
    bool verify(EnumeratedVariable* x, EnumeratedVariable* y, EnumeratedVariable* z);

    template <typename T1, typename T2, typename T3>
    bool isEAC(T1 getCostWithBinaries, bool functionalY, T2 getFunctionY, bool functionalZ, T3 getFunctionZ,
        EnumeratedVariable* x, Value a, EnumeratedVariable* y, EnumeratedVariable* z,
        vector<pair<Value, Value>>& supportX);

    void findSupportX() { findSupport(Functor_getCostXYZ(*this), functionalY, Functor_getFunctionYXZ(*this), functionalZ, Functor_getFunctionZXY(*this), x, y, z, 0, 1, 2, supportX, deltaCostsX, supportY, supportZ); }
    void findSupportY() { findSupport(Functor_getCostYXZ(*this), functionalX, Functor_getFunctionXYZ(*this), functionalZ, Functor_getFunctionZYX(*this), y, x, z, 1, 0, 2, supportY, deltaCostsY, supportX, supportZ); }
    void findSupportZ() { findSupport(Functor_getCostZXY(*this), functionalX, Functor_getFunctionXZY(*this), functionalY, Functor_getFunctionYZX(*this), z, x, y, 2, 0, 1, supportZ, deltaCostsZ, supportX, supportY); }
    void findFullSupportX()
    {
#ifdef NO_STORE_TERNARY_COSTS
        if (Store::getDepth() > 0) {
            findSupportX();
            return;
        }
#endif
        if (y->wcspIndex < z->wcspIndex)
            findFullSupport(Functor_getCostXYZ(*this), Functor_getCostXZY(*this), Functor_getCostYZX(*this), Functor_getCostWithBinariesXYZ(*this), Functor_addCostXYZ(*this), Functor_addCostXZY(*this), Functor_addCostYZX(*this), functionalX, Functor_getFunctionXYZ(*this), functionalY, Functor_getFunctionYXZ(*this), functionalZ, Functor_getFunctionZXY(*this), x, y, z, 0, 1, 2, supportX, deltaCostsX, supportY, deltaCostsY, supportZ, deltaCostsZ, xy, xz, yz);
        else
            findFullSupport(Functor_getCostXZY(*this), Functor_getCostXYZ(*this), Functor_getCostZYX(*this), Functor_getCostWithBinariesXZY(*this), Functor_addCostXZY(*this), Functor_addCostXYZ(*this), Functor_addCostZYX(*this), functionalX, Functor_getFunctionXZY(*this), functionalZ, Functor_getFunctionZXY(*this), functionalY, Functor_getFunctionYXZ(*this), x, z, y, 0, 2, 1, supportX, deltaCostsX, supportZ, deltaCostsZ, supportY, deltaCostsY, xz, xy, yz);
    }
    void findFullSupportY()
    {
#ifdef NO_STORE_TERNARY_COSTS
        if (Store::getDepth() > 0) {
            findSupportY();
            return;
        }
#endif
        if (x->wcspIndex < z->wcspIndex)
            findFullSupport(Functor_getCostYXZ(*this), Functor_getCostYZX(*this), Functor_getCostXZY(*this), Functor_getCostWithBinariesYXZ(*this), Functor_addCostYXZ(*this), Functor_addCostYZX(*this), Functor_addCostXZY(*this), functionalY, Functor_getFunctionYXZ(*this), functionalX, Functor_getFunctionXYZ(*this), functionalZ, Functor_getFunctionZYX(*this), y, x, z, 1, 0, 2, supportY, deltaCostsY, supportX, deltaCostsX, supportZ, deltaCostsZ, xy, yz, xz);
        else
            findFullSupport(Functor_getCostYZX(*this), Functor_getCostYXZ(*this), Functor_getCostZXY(*this), Functor_getCostWithBinariesYZX(*this), Functor_addCostYZX(*this), Functor_addCostYXZ(*this), Functor_addCostZXY(*this), functionalY, Functor_getFunctionYZX(*this), functionalZ, Functor_getFunctionZYX(*this), functionalX, Functor_getFunctionXYZ(*this), y, z, x, 1, 2, 0, supportY, deltaCostsY, supportZ, deltaCostsZ, supportX, deltaCostsX, yz, xy, xz);
    }
    void findFullSupportZ()
    {
#ifdef NO_STORE_TERNARY_COSTS
        if (Store::getDepth() > 0) {
            findSupportZ();
            return;
        }
#endif
        if (x->wcspIndex < y->wcspIndex)
            findFullSupport(Functor_getCostZXY(*this), Functor_getCostZYX(*this), Functor_getCostXYZ(*this), Functor_getCostWithBinariesZXY(*this), Functor_addCostZXY(*this), Functor_addCostZYX(*this), Functor_addCostXYZ(*this), functionalZ, Functor_getFunctionZXY(*this), functionalX, Functor_getFunctionXZY(*this), functionalY, Functor_getFunctionYZX(*this), z, x, y, 2, 0, 1, supportZ, deltaCostsZ, supportX, deltaCostsX, supportY, deltaCostsY, xz, yz, xy);
        else
            findFullSupport(Functor_getCostZYX(*this), Functor_getCostZXY(*this), Functor_getCostYXZ(*this), Functor_getCostWithBinariesZYX(*this), Functor_addCostZYX(*this), Functor_addCostZXY(*this), Functor_addCostYXZ(*this), functionalZ, Functor_getFunctionZYX(*this), functionalY, Functor_getFunctionYZX(*this), functionalX, Functor_getFunctionXZY(*this), z, y, x, 2, 1, 0, supportZ, deltaCostsZ, supportY, deltaCostsY, supportX, deltaCostsX, yz, xz, xy);
    }
    bool verifyX() { return verify(x, y, z); }
    bool verifyY() { return verify(y, x, z); }
    bool verifyZ() { return verify(z, x, y); }

public:
    TernaryConstraint(WCSP* wcsp,
        EnumeratedVariable* xx,
        EnumeratedVariable* yy,
        EnumeratedVariable* zz,
        BinaryConstraint* xy,
        BinaryConstraint* xz,
        BinaryConstraint* yz,
        vector<Cost>& tab);

    TernaryConstraint(WCSP* wcsp);

    void setBinaries(BinaryConstraint* xyin, BinaryConstraint* xzin, BinaryConstraint* yzin)
    {
        xy = xyin;
        xz = xzin;
        yz = yzin;
    }

    BinaryConstraint* xy;
    BinaryConstraint* xz;
    BinaryConstraint* yz;

    virtual ~TernaryConstraint() {}

    bool extension() const FINAL { return true; }
    bool isTernary() const FINAL { return true; }

    Cost getCost(Value vx, Value vy, Value vz) const
    {
        unsigned int ix = x->toIndex(vx);
        unsigned int iy = y->toIndex(vy);
        unsigned int iz = z->toIndex(vz);
        Cost res = ((costs.empty()) ? ((vx == functionX[iy * sizeZ + iz]) ? (costsYZ[iy * sizeZ + iz] - deltaCostsX[ix] - deltaCostsY[iy] - deltaCostsZ[iz]) : top) : (costs[(size_t)ix * sizeY * sizeZ + (size_t)iy * sizeZ + iz] - deltaCostsX[ix] - deltaCostsY[iy] - deltaCostsZ[iz]));
        assert(res >= MIN_COST);
        return res;
    }

    Cost getCost(EnumeratedVariable* xx, EnumeratedVariable* yy, EnumeratedVariable* zz, Value vx, Value vy, Value vz) const
    {
        unsigned int vindex[3];
        vindex[getIndex(xx)] = xx->toIndex(vx);
        vindex[getIndex(yy)] = yy->toIndex(vy);
        vindex[getIndex(zz)] = zz->toIndex(vz);
        Cost res = ((costs.empty()) ? ((x->toValue(vindex[0]) == functionX[vindex[1] * sizeZ + vindex[2]]) ? (costsYZ[vindex[1] * sizeZ + vindex[2]] - deltaCostsX[vindex[0]] - deltaCostsY[vindex[1]] - deltaCostsZ[vindex[2]]) : top) : (costs[(size_t)vindex[0] * sizeY * sizeZ + (size_t)vindex[1] * sizeZ + vindex[2]] - deltaCostsX[vindex[0]] - deltaCostsY[vindex[1]] - deltaCostsZ[vindex[2]]));
        assert(res >= MIN_COST);
        return res;
    }

    Cost getCost() FINAL
    {
        Value vX = x->getValue();
        Value vY = y->getValue();
        Value vZ = z->getValue();
        return getCost(vX, vY, vZ);
    }

    Cost getCostWithBinaries(Value vx, Value vy, Value vz) const
    {
        unsigned int ix = x->toIndex(vx);
        unsigned int iy = y->toIndex(vy);
        unsigned int iz = z->toIndex(vz);
        Cost res = ((costs.empty()) ? ((vx == functionX[iy * sizeZ + iz]) ? (costsYZ[iy * sizeZ + iz] - deltaCostsX[ix] - deltaCostsY[iy] - deltaCostsZ[iz]) : top) : (costs[(size_t)ix * sizeY * sizeZ + (size_t)iy * sizeZ + iz] - deltaCostsX[ix] - deltaCostsY[iy] - deltaCostsZ[iz]));
        if (xy->connected())
            res += xy->getCost(x, y, vx, vy);
        if (xz->connected())
            res += xz->getCost(x, z, vx, vz);
        if (yz->connected())
            res += yz->getCost(y, z, vy, vz);
        assert(res >= MIN_COST);
        return res;
    }

    Cost getCostWithBinaries(EnumeratedVariable* xx, EnumeratedVariable* yy, EnumeratedVariable* zz, Value vx, Value vy, Value vz) const
    {
        pair<unsigned int, Value> vindex[3];
        vindex[getIndex(xx)] = pair<unsigned int, Value>(xx->toIndex(vx), vx);
        vindex[getIndex(yy)] = pair<unsigned int, Value>(yy->toIndex(vy), vy);
        vindex[getIndex(zz)] = pair<unsigned int, Value>(zz->toIndex(vz), vz);
        Cost res = ((costs.empty()) ? ((vindex[0].second == functionX[vindex[1].first * sizeZ + vindex[2].first]) ? (costsYZ[vindex[1].first * sizeZ + vindex[2].first] - deltaCostsX[vindex[0].first] - deltaCostsY[vindex[1].first] - deltaCostsZ[vindex[2].first]) : top) : (costs[(size_t)vindex[0].first * sizeY * sizeZ + (size_t)vindex[1].first * sizeZ + vindex[2].first] - deltaCostsX[vindex[0].first] - deltaCostsY[vindex[1].first] - deltaCostsZ[vindex[2].first]));
        if (xy->connected())
            res += xy->getCost(x, y, vindex[0].second, vindex[1].second);
        if (xz->connected())
            res += xz->getCost(x, z, vindex[0].second, vindex[2].second);
        if (yz->connected())
            res += yz->getCost(y, z, vindex[1].second, vindex[2].second);
        assert(res >= MIN_COST);
        return res;
    }

    void addCosts(TernaryConstraint* xyz)
    {
#ifdef NO_STORE_TERNARY_COSTS
        if (Store::getDepth() > 0) {
            cerr << "Cannot modify costs in ternary cost functions during search!" << endl;
            throw BadConfiguration();
        }
#endif
        unsigned int ix, iy, iz;
        for (EnumeratedVariable::iterator itery = y->begin(); itery != y->end(); ++itery) {
            for (EnumeratedVariable::iterator iterz = z->begin(); iterz != z->end(); ++iterz) {
                for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
                    ix = x->toIndex(*iterx);
                    iy = y->toIndex(*itery);
                    iz = z->toIndex(*iterz);
                    // if(costs[ix*sizeY*sizeZ + iy*sizeZ + iz] < wcsp->getUb()) // BUG with BTD (local ub, deltaCosts missing)
                    if (costs.empty()) {
                        if (*iterx == functionX[iy * sizeZ + iz])
                            costsYZ[iy * sizeZ + iz] += xyz->getCost(x, y, z, *iterx, *itery, *iterz);
                    } else
                        costs[(size_t)ix * sizeY * sizeZ + (size_t)iy * sizeZ + iz] += xyz->getCost(x, y, z, *iterx, *itery, *iterz);
                }
            }
        }
    }

    void addCosts(EnumeratedVariable* xin, EnumeratedVariable* yin, EnumeratedVariable* zin, vector<Cost>& costsin)
    {
#ifdef NO_STORE_TERNARY_COSTS
        if (Store::getDepth() > 0) {
            cerr << "Cannot modify costs in ternary cost functions during search!" << endl;
            throw BadConfiguration();
        }
#endif
        assert(costsin.size() <= costs.size() || functionalX);

        unsigned int vindex[3];
        unsigned int sizeYin = yin->getDomainInitSize();
        unsigned int sizeZin = zin->getDomainInitSize();

        for (EnumeratedVariable::iterator itery = yin->begin(); itery != yin->end(); ++itery) {
            for (EnumeratedVariable::iterator iterz = zin->begin(); iterz != zin->end(); ++iterz) {
                for (EnumeratedVariable::iterator iterx = xin->begin(); iterx != xin->end(); ++iterx) {

                    unsigned int vxin = xin->toIndex(*iterx);
                    unsigned int vyin = yin->toIndex(*itery);
                    unsigned int vzin = zin->toIndex(*iterz);

                    vindex[getIndex(xin)] = vxin;
                    vindex[getIndex(yin)] = vyin;
                    vindex[getIndex(zin)] = vzin;

                    // if(costs[vindex[0]*sizeY*sizeZ + vindex[1]*sizeZ + vindex[2]]  < wcsp->getUb()) // BUG with BTD (local ub, deltaCosts missing)
                    if (costs.empty()) {
                        if (x->toValue(vindex[0]) == functionX[vindex[1] * sizeZ + vindex[2]])
                            costsYZ[vindex[1] * sizeZ + vindex[2]] += costsin[vxin * sizeYin * sizeZin + vyin * sizeZin + vzin];
                    } else
                        costs[(size_t)vindex[0] * sizeY * sizeZ + (size_t)vindex[1] * sizeZ + vindex[2]] += costsin[(size_t)vxin * sizeYin * sizeZin + (size_t)vyin * sizeZin + vzin];
                }
            }
        }
    }

    void addCost(Value vxi, Value vyi, Value vzi, Cost c)
    {
#ifdef NO_STORE_TERNARY_COSTS
        if (Store::getDepth() > 0) {
            cerr << "Cannot modify costs in ternary cost functions during search!" << endl;
            throw BadConfiguration();
        }
#endif
        assert(c >= MIN_COST || !LUBTEST(getCost(vxi, vyi, vzi), -c));
        unsigned int vx = x->toIndex(vxi);
        unsigned int vy = y->toIndex(vyi);
        unsigned int vz = z->toIndex(vzi);
        if (costs.empty()) {
            if (vxi == functionX[vy * sizeZ + vz])
                costsYZ[vy * sizeZ + vz] += c;
        } else {
            if (c < MIN_COST && (functionalX || functionalY || functionalZ)) {
                if ((!functionalX || getFunctionX(vyi, vzi) == vxi) && (!functionalY || getFunctionY(vxi, vzi) == vyi) && (!functionalZ || getFunctionZ(vxi, vyi) == vzi)) {
                    costs[(size_t)vx * sizeY * sizeZ + (size_t)vy * sizeZ + vz] += c; // does not subtract infinity if known by a functional constraint
                }
            } else {
                costs[(size_t)vx * sizeY * sizeZ + (size_t)vy * sizeZ + vz] += c;
            }
        }
    }

    void addCost(EnumeratedVariable* xin, EnumeratedVariable* yin, EnumeratedVariable* zin, Value vxi, Value vyi, Value vzi, Cost c)
    {
#ifdef NO_STORE_TERNARY_COSTS
        if (Store::getDepth() > 0) {
            cerr << "Cannot modify costs in ternary cost functions during search!" << endl;
            throw BadConfiguration();
        }
#endif
        assert(c >= MIN_COST || !LUBTEST(getCost(xin, yin, zin, vxi, vyi, vzi), -c));

        unsigned int vindex[3];
        unsigned int vx = xin->toIndex(vxi);
        unsigned int vy = yin->toIndex(vyi);
        unsigned int vz = zin->toIndex(vzi);

        vindex[getIndex(xin)] = vx;
        vindex[getIndex(yin)] = vy;
        vindex[getIndex(zin)] = vz;

        if (costs.empty()) {
            if (x->toValue(vindex[0]) == functionX[vindex[1] * sizeZ + vindex[2]])
                costsYZ[vindex[1] * sizeZ + vindex[2]] += c;
        } else {
            if (c < MIN_COST && (functionalX || functionalY || functionalZ)) {
                Value valxi = x->toValue(vindex[0]);
                Value valyi = y->toValue(vindex[1]);
                Value valzi = z->toValue(vindex[2]);
                if ((!functionalX || getFunctionX(valyi, valzi) == valxi) && (!functionalY || getFunctionY(valxi, valzi) == valyi) && (!functionalZ || getFunctionZ(valxi, valyi) == valzi)) {
                    costs[(size_t)vindex[0] * sizeY * sizeZ + (size_t)vindex[1] * sizeZ + vindex[2]] += c; // does not subtract infinity if known by a functional constraint
                }
            } else {
                costs[(size_t)vindex[0] * sizeY * sizeZ + (size_t)vindex[1] * sizeZ + vindex[2]] += c;
            }
        }
    }

    void setcost(Value vxi, Value vyi, Value vzi, Cost c)
    {
#ifdef NO_STORE_TERNARY_COSTS
        if (Store::getDepth() > 0) {
            cerr << "Cannot modify costs in ternary cost functions during search!" << endl;
            throw BadConfiguration();
        }
#endif
        assert(std::all_of(deltaCostsX.begin(), deltaCostsX.end(), [](Cost c) { return c == MIN_COST; }));
        assert(std::all_of(deltaCostsY.begin(), deltaCostsY.end(), [](Cost c) { return c == MIN_COST; }));
        assert(std::all_of(deltaCostsZ.begin(), deltaCostsZ.end(), [](Cost c) { return c == MIN_COST; }));
        unsigned int vx = x->toIndex(vxi);
        unsigned int vy = y->toIndex(vyi);
        unsigned int vz = z->toIndex(vzi);
        if (costs.empty()) {
            if (vxi == functionX[vy * sizeZ + vz])
                costsYZ[vy * sizeZ + vz] = c;
            else if (!CUT(wcsp->getLb() + c, wcsp->getUb())) {
                cerr << "cannot reset a forbidden tuple in ternary functional cost functions!" << endl;
                throw InternalError();
            }
        } else
            costs[(size_t)vx * sizeY * sizeZ + (size_t)vy * sizeZ + vz] = c;
    }

    void setcost(EnumeratedVariable* xin, EnumeratedVariable* yin, EnumeratedVariable* zin, Value vxi, Value vyi, Value vzi, Cost c)
    {
#ifdef NO_STORE_TERNARY_COSTS
        if (Store::getDepth() > 0) {
            cerr << "Cannot modify costs in ternary cost functions during search!" << endl;
            throw BadConfiguration();
        }
#endif
        assert(std::all_of(deltaCostsX.begin(), deltaCostsX.end(), [](Cost c) { return c == MIN_COST; }));
        assert(std::all_of(deltaCostsY.begin(), deltaCostsY.end(), [](Cost c) { return c == MIN_COST; }));
        assert(std::all_of(deltaCostsZ.begin(), deltaCostsZ.end(), [](Cost c) { return c == MIN_COST; }));
        unsigned int vindex[3];
        unsigned int vx = xin->toIndex(vxi);
        unsigned int vy = yin->toIndex(vyi);
        unsigned int vz = zin->toIndex(vzi);
        vindex[getIndex(xin)] = vx;
        vindex[getIndex(yin)] = vy;
        vindex[getIndex(zin)] = vz;
        if (costs.empty()) {
            if (x->toValue(vindex[0]) == functionX[vindex[1] * sizeZ + vindex[2]])
                costsYZ[vindex[1] * sizeZ + vindex[2]] = c;
            else if (!CUT(c, wcsp->getUb())) {
                cerr << "cannot reset a forbidden tuple in ternary functional cost functions!" << endl;
                throw InternalError();
            }
        } else
            costs[(size_t)vindex[0] * sizeY * sizeZ + (size_t)vindex[1] * sizeZ + vindex[2]] = c;
    }

    /// \warning Cannot remove a forbidden tuple in ternary functional cost functions!
    void clearCosts()
    {
#ifdef NO_STORE_TERNARY_COSTS
        if (Store::getDepth() > 0) {
            cerr << "Cannot modify costs in ternary cost functions during search!" << endl;
            throw BadConfiguration();
        }
#endif
        assert(ToulBar2::verbose < 4 || ((cout << "[" << Store::getDepth() << ",W" << wcsp->getIndex() << "] clearcosts(C" << getVar(0)->getName() << "," << getVar(1)->getName() << "," << getVar(2)->getName() << ")" << endl), true));
        for (unsigned int i = 0; i < sizeX; i++)
            deltaCostsX[i] = MIN_COST;
        for (unsigned int j = 0; j < sizeY; j++)
            deltaCostsY[j] = MIN_COST;
        for (unsigned int k = 0; k < sizeZ; k++)
            deltaCostsZ[k] = MIN_COST;
        if (costs.empty()) {
            cerr << "Cannot remove a forbidden tuple in ternary functional cost functions!" << endl;
            throw InternalError();
        } else {
            for (unsigned int i = 0; i < sizeX; i++) {
                for (unsigned int j = 0; j < sizeY; j++) {
                    for (unsigned int k = 0; k < sizeZ; k++) {
                        costs[i * sizeY * sizeZ + j * sizeZ + k] = MIN_COST;
                    }
                }
            }
        }
    }

    void clearFiniteCosts()
    {
#ifdef NO_STORE_TERNARY_COSTS
        if (Store::getDepth() > 0) {
            cerr << "Cannot modify finite costs in ternary cost functions during search!" << endl;
            throw BadConfiguration();
        }
#endif
        assert(ToulBar2::verbose < 4 || ((cout << "[" << Store::getDepth() << ",W" << wcsp->getIndex() << "] clearfinitecosts(C" << getVar(0)->getName() << "," << getVar(1)->getName() << "," << getVar(2)->getName() << ")" << endl), true));
        for (unsigned int i = 0; i < sizeX; i++)
            deltaCostsX[i] = MIN_COST;
        for (unsigned int j = 0; j < sizeY; j++)
            deltaCostsY[j] = MIN_COST;
        for (unsigned int k = 0; k < sizeZ; k++)
            deltaCostsZ[k] = MIN_COST;
        if (costs.empty()) {
            for (unsigned int j = 0; j < sizeY; j++) {
                for (unsigned int k = 0; k < sizeZ; k++) {
                    assert(!CUT(costsYZ[j * sizeZ + k], top));
                    if (!CUT(costsYZ[j * sizeZ + k], wcsp->getUb())) {
                        costsYZ[j * sizeZ + k] = MIN_COST;
                    }
                }
            }
        } else {
            for (unsigned int i = 0; i < sizeX; i++) {
                for (unsigned int j = 0; j < sizeY; j++) {
                    for (unsigned int k = 0; k < sizeZ; k++) {
                        if (!CUT(costs[i * sizeY * sizeZ + j * sizeZ + k], wcsp->getUb())) {
                            costs[i * sizeY * sizeZ + j * sizeZ + k] = MIN_COST;
                        }
                    }
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
                for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
                    Cost cost = getCost(*iterX, *iterY, *iterZ);
                    if (cost < ub && cost > maxcost)
                        maxcost = cost;
                }
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
                for (EnumeratedVariable::iterator iterz = z->begin(); iterz != z->end(); ++iterz) {
                    unsigned int iz = z->toIndex(*iterz);
                    if (costs.empty()) {
                        if (*iterx == functionX[iy * sizeZ + iz]) {
                            Cost cost = costsYZ[iy * sizeZ + iz] + x->getCost(*iterx) + y->getCost(*itery) + z->getCost(*iterz);
                            Cost delta = deltaCostsX[ix] + deltaCostsY[iy] + deltaCostsZ[iz];
                            if (CUT(cost - delta, ub)) {
                                costsYZ[iy * sizeZ + iz] = mult_ub + delta;
                                modified = true;
                            }
                        }
                    } else {
                        Cost cost = costs[(size_t)ix * sizeY * sizeZ + (size_t)iy * sizeZ + iz] + x->getCost(*iterx) + y->getCost(*itery) + z->getCost(*iterz);
                        Cost delta = deltaCostsX[ix] + deltaCostsY[iy] + deltaCostsZ[iz];
                        if (CUT(cost - delta, ub)) {
                            costs[(size_t)ix * sizeY * sizeZ + (size_t)iy * sizeZ + iz] = mult_ub + delta;
                            modified = true;
                        }
                    }
                }
            }
        }
        if (costs.empty()) {
            if (CUT(top, ub))
                top = mult_ub;
        }
        if (modified) {
            propagate();
        }
    }

    void resetSupports()
    {
        using std::make_pair;
        supportX.assign(supportX.size(), make_pair(y->getInf(), z->getInf()));
        supportY.assign(supportY.size(), make_pair(x->getInf(), z->getInf()));
        supportZ.assign(supportZ.size(), make_pair(x->getInf(), y->getInf()));
    }

    pair<Value, Value> getSupport(EnumeratedVariable* var, Value v)
    {
        if (var == x)
            return supportX[x->toIndex(v)];
        else if (var == y)
            return supportY[y->toIndex(v)];
        else
            return supportZ[z->toIndex(v)];
    }

    void setSupport(EnumeratedVariable* var, Value v, pair<Value, Value> s)
    {
        if (var == x)
            supportX[x->toIndex(v)] = s;
        else if (var == y)
            supportY[y->toIndex(v)] = s;
        else
            supportZ[z->toIndex(v)] = s;
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
        if (z->assigned()) {
            assign(2);
            return;
        }
        x->queueAC();
        y->queueAC();
        z->queueAC();
        x->queueDAC();
        y->queueDAC();
        z->queueDAC();
        x->queueEAC1();
        y->queueEAC1();
        z->queueEAC1();
        if (ToulBar2::FullEAC)
            reviseEACGreedySolution();
    }

    void remove(int varIndex)
    {
        switch (varIndex) {
        case 0:
            y->queueDEE();
            z->queueDEE();
            if (ToulBar2::LcLevel == LC_AC || getDACScopeIndex() != 1)
                findSupportY();
            if (connected() && (ToulBar2::LcLevel == LC_AC || getDACScopeIndex() != 2))
                findSupportZ();
            break;
        case 1:
            x->queueDEE();
            z->queueDEE();
            if (ToulBar2::LcLevel == LC_AC || getDACScopeIndex() != 0)
                findSupportX();
            if (connected() && (ToulBar2::LcLevel == LC_AC || getDACScopeIndex() != 2))
                findSupportZ();
            break;
        case 2:
            x->queueDEE();
            y->queueDEE();
            if (ToulBar2::LcLevel == LC_AC || getDACScopeIndex() != 0)
                findSupportX();
            if (connected() && (ToulBar2::LcLevel == LC_AC || getDACScopeIndex() != 1))
                findSupportY();
            break;
        }
    }

    void projectFromZero(int varIndex)
    {
        switch (varIndex) {
        case 0:
            if (getDACScopeIndex() == 1)
                findFullSupportY();
            else if (getDACScopeIndex() == 2)
                findFullSupportZ();
            break;
        case 1:
            if (getDACScopeIndex() == 0)
                findFullSupportX();
            else if (getDACScopeIndex() == 2)
                findFullSupportZ();
            break;
        case 2:
            if (getDACScopeIndex() == 0)
                findFullSupportX();
            else if (getDACScopeIndex() == 1)
                findFullSupportY();
            break;
        }
    }

    // Trick! instead of doing remove(index) now, let AC queue do the job.
    // So several incdec events on the same constraint can be merged into one AC event
    void increase(int index)
    {
        if (index == 0)
            x->queueAC();
        else if (index == 1)
            y->queueAC();
        else
            z->queueAC();
    }
    void decrease(int index)
    {
        if (index == 0)
            x->queueAC();
        else if (index == 1)
            y->queueAC();
        else
            z->queueAC();
    }

    void assign(int varIndex)
    {
        deconnect();
        switch (varIndex) {
        case 0:
            projectTernaryBinary(Functor_getCostXYZ(*this), Functor_getCostYZX(*this), Functor_addCostYZX(*this), functionalX, Functor_getFunctionXYZ(*this), x, y, z, yz);
            break;
        case 1:
            projectTernaryBinary(Functor_getCostYXZ(*this), Functor_getCostXZY(*this), Functor_addCostXZY(*this), functionalY, Functor_getFunctionYXZ(*this), y, x, z, xz);
            break;
        case 2:
            projectTernaryBinary(Functor_getCostZXY(*this), Functor_getCostXYZ(*this), Functor_addCostXYZ(*this), functionalZ, Functor_getFunctionZXY(*this), z, x, y, xy);
            break;
        }
    }

    bool checkEACGreedySolution(int index = -1, Value a = 0) FINAL
    {
        assert(x->canbe((getIndex(x) == index) ? a : x->getSupport()));
        assert(x->getCost((getIndex(x) == index) ? a : x->getSupport()) == MIN_COST);
        assert(y->canbe((getIndex(y) == index) ? a : y->getSupport()));
        assert(y->getCost((getIndex(y) == index) ? a : y->getSupport()) == MIN_COST);
        assert(z->canbe((getIndex(z) == index) ? a : z->getSupport()));
        assert(z->getCost((getIndex(z) == index) ? a : z->getSupport()) == MIN_COST);
        return (getCostWithBinaries(x, y, z, (getIndex(x) == index) ? a : x->getSupport(), (getIndex(y) == index) ? a : y->getSupport(), (getIndex(z) == index) ? a : z->getSupport()) == MIN_COST);
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
                z->unsetFullEAC();
            }
        }
        return result;
    }

    void fillEAC2(EnumeratedVariable* x, EnumeratedVariable* y, EnumeratedVariable* z,
        vector<pair<Value, Value>>& supportX)
    {
        assert(x->canbe(x->getSupport()));
        assert(getIndex(y) < getIndex(z));
        unsigned int xindex = x->toIndex(x->getSupport());
        Value ysupport = supportX[xindex].first;
        Value zsupport = supportX[xindex].second;
        if (y->cannotbe(ysupport) || z->cannotbe(zsupport) || getCostWithBinaries(x, y, z, x->getSupport(), ysupport, zsupport) + y->getCost(ysupport) + z->getCost(zsupport) > MIN_COST || (ToulBar2::vacValueHeuristic && Store::getDepth() < abs(ToulBar2::vac))) {
            x->queueEAC2();
        }
    }

    void fillEAC2(int varIndex)
    {
        assert(!isDuplicate());
        switch (varIndex) {
        case 0:
            fillEAC2(y, x, z, supportY);
            fillEAC2(z, x, y, supportZ);
            break;
        case 1:
            fillEAC2(x, y, z, supportX);
            fillEAC2(z, x, y, supportZ);
            break;
        case 2:
            fillEAC2(x, y, z, supportX);
            fillEAC2(y, x, z, supportY);
            break;
        }
    }

    bool isEAC(int varIndex, Value a)
    {
        assert(!isDuplicate());
        if (ToulBar2::QueueComplexity && varIndex == getDACScopeIndex() && !ToulBar2::FullEAC)
            return true;
        switch (varIndex) {
        case 0:
            return isEAC(Functor_getCostWithBinariesXYZ(*this), functionalY, Functor_getFunctionYXZ(*this), functionalZ, Functor_getFunctionZXY(*this), x, a, y, z, supportX);
            break;
        case 1:
            return isEAC(Functor_getCostWithBinariesYXZ(*this), functionalX, Functor_getFunctionXYZ(*this), functionalZ, Functor_getFunctionZYX(*this), y, a, x, z, supportY);
            break;
        case 2:
            return isEAC(Functor_getCostWithBinariesZXY(*this), functionalX, Functor_getFunctionXZY(*this), functionalY, Functor_getFunctionYZX(*this), z, a, x, y, supportZ);
            break;
        default:
            throw InternalError();
        }
        return true;
    }

    bool checkTreeDecomposition(); ///< \brief if tree decomposition then xy, xz, xz binary constraints attached to this ternary should all belong to the same cluster

    void findFullSupportEAC(int varIndex)
    {
        assert(!isDuplicate());
        if (ToulBar2::QueueComplexity && varIndex == getDACScopeIndex() && !ToulBar2::FullEAC)
            return;
        assert(checkTreeDecomposition());
        switch (varIndex) {
        case 0:
            findFullSupportX();
            break;
        case 1:
            findFullSupportY();
            break;
        case 2:
            findFullSupportZ();
            break;
        }
    }

    bool verify();

    pair<pair<Cost, Cost>, pair<Cost, Cost>> getMaxCost(int varIndex, Value a, Value b);

    template <typename T1, typename T2, typename T3, typename T4>
    void projectTernaryBinary(T1 getCost, T2 getCostYZX, T3 addCostYZX, bool functionalX, T4 getFunctionX, EnumeratedVariable* xx, EnumeratedVariable* yy, EnumeratedVariable* zz, BinaryConstraint* yzin);

    void projectTernary()
    {
        if (x->getDACOrder() < y->getDACOrder() && y->getDACOrder() < z->getDACOrder()) {
            projectTernaryBinary(Functor_getCostZXY(*this), Functor_getCostXYZ(*this), Functor_addCostXYZ(*this), functionalZ, Functor_getFunctionZXY(*this), z, x, y, xy);
            if (connected())
                projectTernaryBinary(Functor_getCostYXZ(*this), Functor_getCostXZY(*this), Functor_addCostXZY(*this), functionalY, Functor_getFunctionYXZ(*this), y, x, z, xz);
            if (connected())
                projectTernaryBinary(Functor_getCostXYZ(*this), Functor_getCostYZX(*this), Functor_addCostYZX(*this), functionalX, Functor_getFunctionXYZ(*this), x, y, z, yz);
        } else if (x->getDACOrder() < z->getDACOrder() && z->getDACOrder() < y->getDACOrder()) {
            projectTernaryBinary(Functor_getCostYXZ(*this), Functor_getCostXZY(*this), Functor_addCostXZY(*this), functionalY, Functor_getFunctionYXZ(*this), y, x, z, xz);
            if (connected())
                projectTernaryBinary(Functor_getCostZXY(*this), Functor_getCostXYZ(*this), Functor_addCostXYZ(*this), functionalZ, Functor_getFunctionZXY(*this), z, x, y, xy);
            if (connected())
                projectTernaryBinary(Functor_getCostXYZ(*this), Functor_getCostYZX(*this), Functor_addCostYZX(*this), functionalX, Functor_getFunctionXYZ(*this), x, y, z, yz);
        } else if (y->getDACOrder() < x->getDACOrder() && x->getDACOrder() < z->getDACOrder()) {
            projectTernaryBinary(Functor_getCostZXY(*this), Functor_getCostXYZ(*this), Functor_addCostXYZ(*this), functionalZ, Functor_getFunctionZXY(*this), z, x, y, xy);
            if (connected())
                projectTernaryBinary(Functor_getCostXYZ(*this), Functor_getCostYZX(*this), Functor_addCostYZX(*this), functionalX, Functor_getFunctionXYZ(*this), x, y, z, yz);
            if (connected())
                projectTernaryBinary(Functor_getCostYXZ(*this), Functor_getCostXZY(*this), Functor_addCostXZY(*this), functionalY, Functor_getFunctionYXZ(*this), y, x, z, xz);
        } else if (y->getDACOrder() < z->getDACOrder() && z->getDACOrder() < x->getDACOrder()) {
            projectTernaryBinary(Functor_getCostXYZ(*this), Functor_getCostYZX(*this), Functor_addCostYZX(*this), functionalX, Functor_getFunctionXYZ(*this), x, y, z, yz);
            if (connected())
                projectTernaryBinary(Functor_getCostZXY(*this), Functor_getCostXYZ(*this), Functor_addCostXYZ(*this), functionalZ, Functor_getFunctionZXY(*this), z, x, y, xy);
            if (connected())
                projectTernaryBinary(Functor_getCostYXZ(*this), Functor_getCostXZY(*this), Functor_addCostXZY(*this), functionalY, Functor_getFunctionYXZ(*this), y, x, z, xz);
        } else if (z->getDACOrder() < x->getDACOrder() && x->getDACOrder() < y->getDACOrder()) {
            projectTernaryBinary(Functor_getCostYXZ(*this), Functor_getCostXZY(*this), Functor_addCostXZY(*this), functionalY, Functor_getFunctionYXZ(*this), y, x, z, xz);
            if (connected())
                projectTernaryBinary(Functor_getCostXYZ(*this), Functor_getCostYZX(*this), Functor_addCostYZX(*this), functionalX, Functor_getFunctionXYZ(*this), x, y, z, yz);
            if (connected())
                projectTernaryBinary(Functor_getCostZXY(*this), Functor_getCostXYZ(*this), Functor_addCostXYZ(*this), functionalZ, Functor_getFunctionZXY(*this), z, x, y, xy);
        } else if (z->getDACOrder() < y->getDACOrder() && y->getDACOrder() < x->getDACOrder()) {
            projectTernaryBinary(Functor_getCostXYZ(*this), Functor_getCostYZX(*this), Functor_addCostYZX(*this), functionalX, Functor_getFunctionXYZ(*this), x, y, z, yz);
            if (connected())
                projectTernaryBinary(Functor_getCostYXZ(*this), Functor_getCostXZY(*this), Functor_addCostXZY(*this), functionalY, Functor_getFunctionYXZ(*this), y, x, z, xz);
            if (connected())
                projectTernaryBinary(Functor_getCostZXY(*this), Functor_getCostXYZ(*this), Functor_addCostXYZ(*this), functionalZ, Functor_getFunctionZXY(*this), z, x, y, xy);
        } else
            throw InternalError();
    }

    void extendTernary()
    { // extend binary cost functions to the ternary cost function
        Cost c;
        bool isbincost = false;
        for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
            for (EnumeratedVariable::iterator itery = y->begin(); itery != y->end(); ++itery) {
                for (EnumeratedVariable::iterator iterz = z->begin(); iterz != z->end(); ++iterz) {
                    if (xy->connected()) {
                        c = xy->getCost(x, y, *iterx, *itery);
                        addCost(x, y, z, *iterx, *itery, *iterz, c);
                        if (c > MIN_COST)
                            isbincost = true;
                    }
                    if (xz->connected()) {
                        c = xz->getCost(x, z, *iterx, *iterz);
                        addCost(x, y, z, *iterx, *itery, *iterz, c);
                        if (c > MIN_COST)
                            isbincost = true;
                    }
                    if (yz->connected()) {
                        c = yz->getCost(y, z, *itery, *iterz);
                        addCost(x, y, z, *iterx, *itery, *iterz, c);
                        if (c > MIN_COST)
                            isbincost = true;
                    }
                }
            }
        }

        xy->clearCosts();
        xz->clearCosts();
        yz->clearCosts();

        xy->deconnect(true);
        xz->deconnect(true);
        yz->deconnect(true);

        // extend unary costs to the ternary cost function
        if (isbincost && ToulBar2::LcLevel >= LC_FDAC) {
            for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx)
                extend(x, *iterx, x->getCost(*iterx), deltaCostsX);
            for (EnumeratedVariable::iterator itery = y->begin(); itery != y->end(); ++itery)
                extend(y, *itery, y->getCost(*itery), deltaCostsY);
            for (EnumeratedVariable::iterator iterz = z->begin(); iterz != z->end(); ++iterz)
                extend(z, *iterz, z->getCost(*iterz), deltaCostsZ);
        }
    }

    BinaryConstraint* commonBinary(TernaryConstraint* t)
    {
        if ((t->getIndex(xy->getVar(0)) >= 0) && (t->getIndex(xy->getVar(1)) >= 0))
            return xy;
        else if ((t->getIndex(xz->getVar(0)) >= 0) && (t->getIndex(xz->getVar(1)) >= 0))
            return xz;
        else if ((t->getIndex(yz->getVar(0)) >= 0) && (t->getIndex(yz->getVar(1)) >= 0))
            return yz;
        return NULL;
    }

    double computeTightness();

    // add weights from auxilliary binary constraints if they are deconnected (otherwise their weight will be taken into account in the list of active constraints for variable varIndex)
    // use weights minus one due to default initialization of conflictWeight to 1
    Long getConflictWeight(int varIndex) const
    {
        switch (varIndex) {
        case 0:
            return Constraint::getConflictWeight() + ((xy->deconnected()) ? (xy->getConflictWeight(xy->getIndex(x)) - 1) : 0) + ((xz->deconnected()) ? (xz->getConflictWeight(xz->getIndex(x)) - 1) : 0);
            break;
        case 1:
            return Constraint::getConflictWeight() + ((xy->deconnected()) ? (xy->getConflictWeight(xy->getIndex(y)) - 1) : 0) + ((yz->deconnected()) ? (yz->getConflictWeight(yz->getIndex(y)) - 1) : 0);
            break;
        case 2:
            return Constraint::getConflictWeight() + ((xz->deconnected()) ? (xz->getConflictWeight(xz->getIndex(z)) - 1) : 0) + ((yz->deconnected()) ? (yz->getConflictWeight(yz->getIndex(z)) - 1) : 0);
            break;
        }
        return Constraint::getConflictWeight();
    }
    Long getConflictWeight() const
    {
        return Constraint::getConflictWeight() + ((xy && xy->deconnected()) ? (xy->getConflictWeight() - 1) : 0) + ((xz && xz->deconnected()) ? (xz->getConflictWeight() - 1) : 0) + ((yz && yz->deconnected()) ? (yz->getConflictWeight() - 1) : 0);
    }

    EnumeratedVariable* xvar;
    EnumeratedVariable* yvar;
    EnumeratedVariable* zvar;
    EnumeratedVariable::iterator itvx;
    EnumeratedVariable::iterator itvy;
    EnumeratedVariable::iterator itvz;

    void first()
    {
        itvx = x->begin();
        itvy = y->begin();
        itvz = z->begin();
        xvar = x;
        yvar = y;
        zvar = z;
    }

    bool next(Tuple& t, Cost& c)
    {
        Tuple tch(3, 0);
        if (itvx != xvar->end()) {
            unsigned int ix = xvar->toIndex(*itvx);
            tch[0] = ix;
            if (itvy != yvar->end()) {
                unsigned int iy = yvar->toIndex(*itvy);
                tch[1] = iy;
                if (itvz != zvar->end()) {
                    unsigned int iz = zvar->toIndex(*itvz);
                    tch[2] = iz;
                    t = tch;
                    c = getCost(xvar, yvar, zvar, *itvx, *itvy, *itvz);
                    ++itvz;
                    return true;
                } else {
                    ++itvy;
                    itvz = zvar->begin();
                    return next(t, c);
                }
            } else {
                ++itvx;
                itvy = yvar->begin();
                return next(t, c);
            }
        }
        return false;
    }

    void first(EnumeratedVariable* alpha, EnumeratedVariable* beta)
    {
        int pos_alpha = getIndex(alpha);
        int pos_beta = getIndex(beta);
        xvar = NULL;
        yvar = NULL;
        zvar = NULL;
        switch (pos_alpha) {
        case 0:
            itvz = x->begin();
            zvar = x;
            break;
        case 1:
            itvz = y->begin();
            zvar = y;
            break;
        case 2:
            itvz = z->begin();
            zvar = z;
            break;
        }
        switch (pos_beta) {
        case 0:
            itvy = x->begin();
            yvar = x;
            break;
        case 1:
            itvy = y->begin();
            yvar = y;
            break;
        case 2:
            itvy = z->begin();
            yvar = z;
            break;
        }
        switch (3 - pos_alpha - pos_beta) {
        case 0:
            itvx = x->begin();
            xvar = x;
            break;
        case 1:
            itvx = y->begin();
            xvar = y;
            break;
        case 2:
            itvx = z->begin();
            xvar = z;
            break;
        }
        assert(xvar != yvar);
        assert(xvar != zvar);
        assert(yvar != zvar);
    }
    bool separability(EnumeratedVariable* alpha, EnumeratedVariable* beta) FINAL;
    void separate(EnumeratedVariable* a, EnumeratedVariable* c) FINAL;

    void firstlex() { first(); }
    bool nextlex(Tuple& t, Cost& c) { return next(t, c); }

    void setTuple(const Tuple& t, Cost c) FINAL
    {
        Value v0 = x->toValue(t[0]);
        Value v1 = y->toValue(t[1]);
        Value v2 = z->toValue(t[2]);
        Cost oldc = getCost(v0, v1, v2);
        addCost(v0, v1, v2, c - oldc);
    }

    void addtoTuple(const Tuple& t, Cost c) FINAL
    {
        Value v0 = x->toValue(t[0]);
        Value v1 = y->toValue(t[1]);
        Value v2 = z->toValue(t[2]);
        addCost(v0, v1, v2, c);
    }

    Cost evalsubstr(const Tuple& s, Constraint* ctr) FINAL
    {
        Value vals[3];
        int count = 0;

        for (int i = 0; i < 3; i++) {
            EnumeratedVariable* var = (EnumeratedVariable*)getVar(i);
            int ind = ctr->getIndex(var);
            if (ind >= 0) {
                vals[i] = var->toValue(s[ind]);
                count++;
            }
        }
        if (count == 3)
            return getCost(vals[0], vals[1], vals[2]);
        else
            return MIN_COST;
    }
    Cost evalsubstr(const Tuple& s, NaryConstraint* ctr) FINAL { return evalsubstr(s, (Constraint*)ctr); } // NaryConstraint class undefined

    void fillElimConstr(EnumeratedVariable* xin, EnumeratedVariable* yin, EnumeratedVariable* zin, Constraint* from1)
    {
        assert(!functionalX && costsYZ.empty());
        x = xin;
        y = yin;
        z = zin;
        sizeX = x->getDomainInitSize();
        sizeY = y->getDomainInitSize();
        sizeZ = z->getDomainInitSize();
        if (sizeX > deltaCostsX.size())
            deltaCostsX.resize(sizeX, StoreCost(MIN_COST));
        if (sizeY > deltaCostsY.size())
            deltaCostsY.resize(sizeY, StoreCost(MIN_COST));
        if (sizeZ > deltaCostsZ.size())
            deltaCostsZ.resize(sizeZ, StoreCost(MIN_COST));
        if (sizeX > supportX.size())
            supportX.resize(sizeX);
        if (sizeY > supportY.size())
            supportY.resize(sizeY);
        if (sizeZ > supportZ.size())
            supportZ.resize(sizeZ);
        if ((size_t)sizeX * (size_t)sizeY * (size_t)sizeZ > costs.size())
#ifdef NO_STORE_TERNARY_COSTS
            costs.resize((size_t)sizeX * (size_t)sizeY * (size_t)sizeZ, MIN_COST);
#else
            costs.resize((size_t)sizeX * (size_t)sizeY * (size_t)sizeZ, StoreCost(MIN_COST));
#endif
        linkX->removed = true;
        linkY->removed = true;
        linkZ->removed = true;
        linkX->content.constr = this;
        linkY->content.constr = this;
        linkZ->content.constr = this;
        linkX->content.scopeIndex = 0;
        linkY->content.scopeIndex = 1;
        linkZ->content.scopeIndex = 2;
        setDACScopeIndex();
        // resetConflictWeight(); // warning! it is done before having updated costs!! it will be done later in fillElimConstrBinaries
        elimFrom(from1);
    }

    void fillxy();
    void fillxz();
    void fillyz();
    void fillElimConstrBinaries();
    void setDuplicates();

    void print(ostream& os);
    void dump(ostream& os, bool original = true);
    void dump_CFN(ostream& os, bool original = true);
    Long size() const FINAL { return (Long)sizeX * sizeY * sizeZ; }
#ifdef NO_STORE_TERNARY_COSTS
    Long space() const FINAL
    {
        return (Long)sizeof(Cost) * sizeX * sizeY * sizeZ;
    }
#else
    Long space() const FINAL
    {
        return (Long)sizeof(StoreCost) * sizeX * sizeY * sizeZ;
    }
#endif
    friend struct Functor_getCostXYZ;
    friend struct Functor_getCostXZY;
    friend struct Functor_getCostYXZ;
    friend struct Functor_getCostYZX;
    friend struct Functor_getCostZXY;
    friend struct Functor_getCostZYX;

    friend struct Functor_getCostWithBinariesXYZ;
    friend struct Functor_getCostWithBinariesXZY;
    friend struct Functor_getCostWithBinariesYXZ;
    friend struct Functor_getCostWithBinariesYZX;
    friend struct Functor_getCostWithBinariesZXY;
    friend struct Functor_getCostWithBinariesZYX;

    friend struct Functor_addCostXYZ;
    friend struct Functor_addCostXZY;
    friend struct Functor_addCostYXZ;
    friend struct Functor_addCostYZX;
    friend struct Functor_addCostZXY;
    friend struct Functor_addCostZYX;

    friend struct Functor_getFunctionXYZ;
    friend struct Functor_getFunctionXZY;
    friend struct Functor_getFunctionYXZ;
    friend struct Functor_getFunctionYZX;
    friend struct Functor_getFunctionZXY;
    friend struct Functor_getFunctionZYX;
};

inline Cost Functor_getCostXYZ::operator()(EnumeratedVariable* xx, EnumeratedVariable* yy, EnumeratedVariable* zz, Value vx, Value vy, Value vz) const
{
    assert(xx == obj.x);
    assert(yy == obj.y);
    assert(zz == obj.z);
    return obj.getCost(vx, vy, vz);
}
inline Cost Functor_getCostXZY::operator()(EnumeratedVariable* xx, EnumeratedVariable* zz, EnumeratedVariable* yy, Value vx, Value vz, Value vy) const
{
    assert(xx == obj.x);
    assert(yy == obj.y);
    assert(zz == obj.z);
    return obj.getCost(vx, vy, vz);
}
inline Cost Functor_getCostYXZ::operator()(EnumeratedVariable* yy, EnumeratedVariable* xx, EnumeratedVariable* zz, Value vy, Value vx, Value vz) const
{
    assert(xx == obj.x);
    assert(yy == obj.y);
    assert(zz == obj.z);
    return obj.getCost(vx, vy, vz);
}
inline Cost Functor_getCostYZX::operator()(EnumeratedVariable* yy, EnumeratedVariable* zz, EnumeratedVariable* xx, Value vy, Value vz, Value vx) const
{
    assert(xx == obj.x);
    assert(yy == obj.y);
    assert(zz == obj.z);
    return obj.getCost(vx, vy, vz);
}
inline Cost Functor_getCostZXY::operator()(EnumeratedVariable* zz, EnumeratedVariable* xx, EnumeratedVariable* yy, Value vz, Value vx, Value vy) const
{
    assert(xx == obj.x);
    assert(yy == obj.y);
    assert(zz == obj.z);
    return obj.getCost(vx, vy, vz);
}
inline Cost Functor_getCostZYX::operator()(EnumeratedVariable* zz, EnumeratedVariable* yy, EnumeratedVariable* xx, Value vz, Value vy, Value vx) const
{
    assert(xx == obj.x);
    assert(yy == obj.y);
    assert(zz == obj.z);
    return obj.getCost(vx, vy, vz);
}

inline Cost Functor_getCostWithBinariesXYZ::operator()(EnumeratedVariable* xx, EnumeratedVariable* yy, EnumeratedVariable* zz, Value vx, Value vy, Value vz) const
{
    assert(xx == obj.x);
    assert(yy == obj.y);
    assert(zz == obj.z);
    return obj.getCostWithBinaries(vx, vy, vz);
}
inline Cost Functor_getCostWithBinariesXZY::operator()(EnumeratedVariable* xx, EnumeratedVariable* zz, EnumeratedVariable* yy, Value vx, Value vz, Value vy) const
{
    assert(xx == obj.x);
    assert(yy == obj.y);
    assert(zz == obj.z);
    return obj.getCostWithBinaries(vx, vy, vz);
}
inline Cost Functor_getCostWithBinariesYXZ::operator()(EnumeratedVariable* yy, EnumeratedVariable* xx, EnumeratedVariable* zz, Value vy, Value vx, Value vz) const
{
    assert(xx == obj.x);
    assert(yy == obj.y);
    assert(zz == obj.z);
    return obj.getCostWithBinaries(vx, vy, vz);
}
inline Cost Functor_getCostWithBinariesYZX::operator()(EnumeratedVariable* yy, EnumeratedVariable* zz, EnumeratedVariable* xx, Value vy, Value vz, Value vx) const
{
    assert(xx == obj.x);
    assert(yy == obj.y);
    assert(zz == obj.z);
    return obj.getCostWithBinaries(vx, vy, vz);
}
inline Cost Functor_getCostWithBinariesZXY::operator()(EnumeratedVariable* zz, EnumeratedVariable* xx, EnumeratedVariable* yy, Value vz, Value vx, Value vy) const
{
    assert(xx == obj.x);
    assert(yy == obj.y);
    assert(zz == obj.z);
    return obj.getCostWithBinaries(vx, vy, vz);
}
inline Cost Functor_getCostWithBinariesZYX::operator()(EnumeratedVariable* zz, EnumeratedVariable* yy, EnumeratedVariable* xx, Value vz, Value vy, Value vx) const
{
    assert(xx == obj.x);
    assert(yy == obj.y);
    assert(zz == obj.z);
    return obj.getCostWithBinaries(vx, vy, vz);
}

inline void Functor_addCostXYZ::operator()(EnumeratedVariable* xx, EnumeratedVariable* yy, EnumeratedVariable* zz, Value vx, Value vy, Value vz, Cost c)
{
    assert(xx == obj.x);
    assert(yy == obj.y);
    assert(zz == obj.z);
    return obj.addCost(vx, vy, vz, c);
}
inline void Functor_addCostXZY::operator()(EnumeratedVariable* xx, EnumeratedVariable* zz, EnumeratedVariable* yy, Value vx, Value vz, Value vy, Cost c)
{
    assert(xx == obj.x);
    assert(yy == obj.y);
    assert(zz == obj.z);
    return obj.addCost(vx, vy, vz, c);
}
inline void Functor_addCostYXZ::operator()(EnumeratedVariable* yy, EnumeratedVariable* xx, EnumeratedVariable* zz, Value vy, Value vx, Value vz, Cost c)
{
    assert(xx == obj.x);
    assert(yy == obj.y);
    assert(zz == obj.z);
    return obj.addCost(vx, vy, vz, c);
}
inline void Functor_addCostYZX::operator()(EnumeratedVariable* yy, EnumeratedVariable* zz, EnumeratedVariable* xx, Value vy, Value vz, Value vx, Cost c)
{
    assert(xx == obj.x);
    assert(yy == obj.y);
    assert(zz == obj.z);
    return obj.addCost(vx, vy, vz, c);
}
inline void Functor_addCostZXY::operator()(EnumeratedVariable* zz, EnumeratedVariable* xx, EnumeratedVariable* yy, Value vz, Value vx, Value vy, Cost c)
{
    assert(xx == obj.x);
    assert(yy == obj.y);
    assert(zz == obj.z);
    return obj.addCost(vx, vy, vz, c);
}
inline void Functor_addCostZYX::operator()(EnumeratedVariable* zz, EnumeratedVariable* yy, EnumeratedVariable* xx, Value vz, Value vy, Value vx, Cost c)
{
    assert(xx == obj.x);
    assert(yy == obj.y);
    assert(zz == obj.z);
    return obj.addCost(vx, vy, vz, c);
}

inline Value Functor_getFunctionXYZ::operator()(EnumeratedVariable* yy, EnumeratedVariable* zz, Value vy, Value vz) const
{
    assert(yy == obj.y);
    assert(zz == obj.z);
    return obj.getFunctionX(vy, vz);
}
inline Value Functor_getFunctionXZY::operator()(EnumeratedVariable* zz, EnumeratedVariable* yy, Value vz, Value vy) const
{
    assert(yy == obj.y);
    assert(zz == obj.z);
    return obj.getFunctionX(vy, vz);
}
inline Value Functor_getFunctionYXZ::operator()(EnumeratedVariable* xx, EnumeratedVariable* zz, Value vx, Value vz) const
{
    assert(xx == obj.x);
    assert(zz == obj.z);
    return obj.getFunctionY(vx, vz);
}
inline Value Functor_getFunctionYZX::operator()(EnumeratedVariable* zz, EnumeratedVariable* xx, Value vz, Value vx) const
{
    assert(xx == obj.x);
    assert(zz == obj.z);
    return obj.getFunctionY(vx, vz);
}
inline Value Functor_getFunctionZXY::operator()(EnumeratedVariable* xx, EnumeratedVariable* yy, Value vx, Value vy) const
{
    assert(xx == obj.x);
    assert(yy == obj.y);
    return obj.getFunctionZ(vx, vy);
}
inline Value Functor_getFunctionZYX::operator()(EnumeratedVariable* yy, EnumeratedVariable* xx, Value vy, Value vx) const
{
    assert(xx == obj.x);
    assert(yy == obj.y);
    return obj.getFunctionZ(vx, vy);
}

template <typename T1, typename T2, typename T3>
void TernaryConstraint::project(T1 getCost, T2 addCost, bool functionalZ, T3 getFunctionZ, BinaryConstraint* xy, EnumeratedVariable* x, EnumeratedVariable* y, EnumeratedVariable* z, Value valx, Value valy, Cost cost)
{
    assert(ToulBar2::verbose < 4 || ((cout << "[" << Store::getDepth() << ",W" << wcsp->getIndex() << "] project(C" << getVar(0)->getName() << "," << getVar(1)->getName() << "," << getVar(2)->getName() << ", ((" << x->getName() << "," << valx << "),(" << y->getName() << "," << valy << ")), " << cost << ")" << endl), true));
    assert(cost >= MIN_COST);
// BUG!
// if (functionalZ) {
//   Value valz = getFunctionZ(x,y, valx, valy);
//   if (valz != WRONG_VAL && z->canbe(valz)) {
// 	if (!CUT(getCost(x,y,z, valx,valy,valz), wcsp->getUb())) { // keeps forbidden costs into ternaries to get strong GAC3
// 	  addCost(x,y,z,valx,valy,valz,-cost);
// 	}
//   }
// } else {
#ifndef NO_STORE_TERNARY_COSTS
    for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
        if (!CUT(getCost(x, y, z, valx, valy, *iterZ), wcsp->getUb())) { // keeps forbidden costs into ternaries to get strong GAC3
            addCost(x, y, z, valx, valy, *iterZ, -cost);
        }
    }
#endif
    // }
    xy->addcost(x, y, valx, valy, cost);
#ifdef DEECOMPLETE
    getVar(0)->queueDEE();
    getVar(1)->queueDEE();
    getVar(2)->queueDEE();
#endif
}

template <typename T1, typename T2, typename T3>
void TernaryConstraint::extend(T1 getCost, T2 addCost, bool functionalZ, T3 getFunctionZ, BinaryConstraint* xy, EnumeratedVariable* x, EnumeratedVariable* y, EnumeratedVariable* z, Value valx, Value valy, Cost cost)
{
    assert(ToulBar2::verbose < 4 || ((cout << "[" << Store::getDepth() << ",W" << wcsp->getIndex() << "] extend(C" << getVar(0)->getName() << "," << getVar(1)->getName() << "," << getVar(2)->getName() << ", ((" << x->getName() << "," << valx << "),(" << y->getName() << "," << valy << ")), " << cost << ")" << endl), true));
    assert(cost >= MIN_COST);
    // BUG!
    // if (functionalZ) {
    //   Value valz = getFunctionZ(x,y, valx, valy);
    //   if (valz != WRONG_VAL && z->canbe(valz)) {
    //     addCost(x,y,z,valx,valy,valz,cost);
    //   }
    // } else {
    for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
        addCost(x, y, z, valx, valy, *iterZ, cost);
    }
    // }
    assert(xy->connected());
    xy->addcost(x, y, valx, valy, -cost);
}

template <typename T1, typename T2, typename T3>
void TernaryConstraint::findSupport(T1 getCost, bool functionalY, T2 getFunctionY, bool functionalZ, T3 getFunctionZ,
    EnumeratedVariable* x, EnumeratedVariable* y, EnumeratedVariable* z,
    int getIndexX, int getIndexY, int getIndexZ,
    vector<pair<Value, Value>>& supportX, vector<StoreCost>& deltaCostsX,
    vector<pair<Value, Value>>& supportY, vector<pair<Value, Value>>& supportZ)
{
    assert(getIndex(y) < getIndex(z)); // check that support.first/.second is consistent with y/z parameters
    assert(connected());
    wcsp->revise(this);
    if (ToulBar2::verbose >= 3)
        cout << "[" << Store::getDepth() << ",W" << wcsp->getIndex() << "] findSupport C" << x->getName() << "?"
             << "," << y->getName() << ((functionalY) ? "!" : "") << "," << z->getName() << ((functionalZ) ? "!" : "") << endl;
    if (ToulBar2::verbose >= 7)
        cout << *x << endl
             << *y << endl
             << *z << endl
             << *this;
    bool supportBroken = false;
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        unsigned int xindex = x->toIndex(*iterX);
        pair<Value, Value> support = supportX[xindex];
        if (y->cannotbe(support.first) || z->cannotbe(support.second) || getCost(x, y, z, *iterX, support.first, support.second) > MIN_COST) {
            support = std::make_pair(y->getInf(), z->getInf());
            Cost minCost = MAX_COST;
            if (functionalZ) {
                for (EnumeratedVariable::iterator iterY = y->begin(); minCost > MIN_COST && iterY != y->end(); ++iterY) {
                    Value valZ = getFunctionZ(x, y, *iterX, *iterY);
                    if (valZ != WRONG_VAL && z->canbe(valZ)) {
                        Cost cost = getCost(x, y, z, *iterX, *iterY, valZ);
                        if (GLB(&minCost, cost)) {
                            support = std::make_pair(*iterY, valZ);
                        }
                    }
                }
            } else if (functionalY) {
                for (EnumeratedVariable::iterator iterZ = z->begin(); minCost > MIN_COST && iterZ != z->end(); ++iterZ) {
                    Value valY = getFunctionY(x, z, *iterX, *iterZ);
                    if (valY != WRONG_VAL && y->canbe(valY)) {
                        Cost cost = getCost(x, y, z, *iterX, valY, *iterZ);
                        if (GLB(&minCost, cost)) {
                            support = std::make_pair(valY, *iterZ);
                        }
                    }
                }
            } else {
                for (EnumeratedVariable::iterator iterY = y->begin(); minCost > MIN_COST && iterY != y->end(); ++iterY) {
                    for (EnumeratedVariable::iterator iterZ = z->begin(); minCost > MIN_COST && iterZ != z->end(); ++iterZ) {
                        Cost cost = getCost(x, y, z, *iterX, *iterY, *iterZ);
                        if (GLB(&minCost, cost)) {
                            support = std::make_pair(*iterY, *iterZ);
                        }
                    }
                }
            }
            if (minCost > MIN_COST) {
                supportBroken |= project(x, *iterX, minCost, deltaCostsX);
                if (deconnected())
                    return;
            }
            supportX[xindex] = support;
            assert(getIndexY < getIndexZ);
            // warning! do not break DAC support for the variable in the constraint scope having the smallest wcspIndex
            if (getIndexY != getDACScopeIndex()) {
                unsigned int yindex = y->toIndex(support.first);
                if (getIndexX < getIndexZ)
                    supportY[yindex] = std::make_pair(*iterX, support.second);
                else
                    supportY[yindex] = std::make_pair(support.second, *iterX);
            }
            if (getIndexZ != getDACScopeIndex()) {
                unsigned int zindex = z->toIndex(support.second);
                if (getIndexX < getIndexY)
                    supportZ[zindex] = std::make_pair(*iterX, support.first);
                else
                    supportZ[zindex] = std::make_pair(support.first, *iterX);
            }
        }
    }
    if (supportBroken) {
        x->findSupport();
    }
}

// take into account associated binary constraints and perform unary extension to binary instead of ternary
template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10>
void TernaryConstraint::findFullSupport(T1 getCost, T2 getCostXZY, T3 getCostYZX, T4 getCostWithBinaries, T5 addCost, T6 addCostXZY, T7 addCostYZX, bool functionalX, T8 getFunctionX, bool functionalY, T9 getFunctionY, bool functionalZ, T10 getFunctionZ,
    EnumeratedVariable* x, EnumeratedVariable* y, EnumeratedVariable* z,
    int getIndexX, int getIndexY, int getIndexZ,
    vector<pair<Value, Value>>& supportX, vector<StoreCost>& deltaCostsX,
    vector<pair<Value, Value>>& supportY, vector<StoreCost>& deltaCostsY,
    vector<pair<Value, Value>>& supportZ, vector<StoreCost>& deltaCostsZ,
    BinaryConstraint* xy, BinaryConstraint* xz, BinaryConstraint* yz)
{
    assert(connected());
    wcsp->revise(this);
    if (ToulBar2::verbose >= 3)
        cout << "[" << Store::getDepth() << ",W" << wcsp->getIndex() << "] findFullSupport C" << x->getName() << ((functionalX) ? "!" : "") << "," << y->getName() << ((functionalY) ? "!" : "") << "," << z->getName() << ((functionalZ) ? "!" : "") << endl;
    if (ToulBar2::verbose >= 7)
        cout << *x << endl
             << *y << endl
             << *z << endl
             << *this << *xy << *xz << *yz;
    bool supportBroken = false;
    bool supportReversed = (getIndexY > getIndexZ);
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        unsigned int xindex = x->toIndex(*iterX);
        pair<Value, Value> support = (supportReversed) ? std::make_pair(supportX[xindex].second, supportX[xindex].first) : std::make_pair(supportX[xindex].first, supportX[xindex].second);
        if (y->cannotbe(support.first) || z->cannotbe(support.second) || getCostWithBinaries(x, y, z, *iterX, support.first, support.second) + y->getCost(support.first) + z->getCost(support.second) > MIN_COST) {
            support = std::make_pair(y->getInf(), z->getInf());
            Cost minCost = MAX_COST;
            if (functionalZ) {
                for (EnumeratedVariable::iterator iterY = y->begin(); minCost > MIN_COST && iterY != y->end(); ++iterY) {
                    Value valZ = getFunctionZ(x, y, *iterX, *iterY);
                    if (valZ != WRONG_VAL && z->canbe(valZ)) {
                        Cost cost = getCostWithBinaries(x, y, z, *iterX, *iterY, valZ) + y->getCost(*iterY) + z->getCost(valZ);
                        if (GLB(&minCost, cost)) {
                            support = std::make_pair(*iterY, valZ);
                        }
                    }
                }
            } else if (functionalY) {
                for (EnumeratedVariable::iterator iterZ = z->begin(); minCost > MIN_COST && iterZ != z->end(); ++iterZ) {
                    Value valY = getFunctionY(x, z, *iterX, *iterZ);
                    if (valY != WRONG_VAL && y->canbe(valY)) {
                        Cost cost = getCostWithBinaries(x, y, z, *iterX, valY, *iterZ) + y->getCost(valY) + z->getCost(*iterZ);
                        if (GLB(&minCost, cost)) {
                            support = std::make_pair(valY, *iterZ);
                        }
                    }
                }
            } else {
                for (EnumeratedVariable::iterator iterY = y->begin(); minCost > MIN_COST && iterY != y->end(); ++iterY) {
                    for (EnumeratedVariable::iterator iterZ = z->begin(); minCost > MIN_COST && iterZ != z->end(); ++iterZ) {
                        Cost cost = getCostWithBinaries(x, y, z, *iterX, *iterY, *iterZ) + y->getCost(*iterY) + z->getCost(*iterZ);
                        if (GLB(&minCost, cost)) {
                            support = std::make_pair(*iterY, *iterZ);
                        }
                    }
                }
            }
            if (CUT(minCost + wcsp->getLb(), wcsp->getUb())) {
                supportBroken |= project(x, *iterX, minCost, deltaCostsX);
                if (deconnected())
                    return;
                continue;
            }
            assert(minCost < MAX_COST);

            if (minCost > MIN_COST) {
                // extend unary to binary
                for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                    if (y->getCost(*iterY) > MIN_COST) {
                        Cost costfromy = MIN_COST;
                        for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
                            Cost cost = getCost(x, y, z, *iterX, *iterY, *iterZ);
                            if (LUBTEST(cost, minCost)) {
                                Cost zcost = z->getCost(*iterZ);
                                Cost xycost = (xy->connected()) ? xy->getCost(x, y, *iterX, *iterY) : MIN_COST;
                                Cost xzcost = (xz->connected()) ? xz->getCost(x, z, *iterX, *iterZ) : MIN_COST;
                                Cost yzcost = (yz->connected()) ? yz->getCost(y, z, *iterY, *iterZ) : MIN_COST;
                                Cost remain = minCost - (cost + xycost + xzcost + yzcost + zcost);
                                LUB(&costfromy, remain);
                            }
                        }
                        assert(costfromy <= y->getCost(*iterY));
                        if (costfromy > MIN_COST) {
                            assert(x->unassigned() && y->unassigned());
                            xy->reconnect(); // must be done before using the constraint
                            xy->extend(xy->getIndex(y), *iterY, costfromy);
                        }
                    }
                }
                for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
                    if (z->getCost(*iterZ) > MIN_COST) {
                        Cost costfromz = MIN_COST;
                        for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                            Cost cost = getCost(x, y, z, *iterX, *iterY, *iterZ);
                            if (LUBTEST(cost, minCost)) {
                                Cost xycost = (xy->connected()) ? xy->getCost(x, y, *iterX, *iterY) : MIN_COST;
                                Cost xzcost = (xz->connected()) ? xz->getCost(x, z, *iterX, *iterZ) : MIN_COST;
                                Cost yzcost = (yz->connected()) ? yz->getCost(y, z, *iterY, *iterZ) : MIN_COST;
                                Cost remain = minCost - (cost + xycost + xzcost + yzcost);
                                LUB(&costfromz, remain);
                            }
                        }
                        assert(costfromz <= z->getCost(*iterZ));
                        if (costfromz > MIN_COST) {
                            assert(x->unassigned() && z->unassigned());
                            xz->reconnect(); // must be done before using the constraint
                            xz->extend(xz->getIndex(z), *iterZ, costfromz);
                        }
                    }
                }
                // extend binary to ternary
                if (yz->connected()) {
                    for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                        Cost xycost = (xy->connected()) ? xy->getCost(x, y, *iterX, *iterY) : MIN_COST;
                        for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
                            Cost yzcost = yz->getCost(y, z, *iterY, *iterZ);
                            if (yzcost > MIN_COST) {
                                Cost cost = getCost(x, y, z, *iterX, *iterY, *iterZ);
                                if (LUBTEST(cost, minCost)) {
                                    Cost xzcost = (xz->connected()) ? xz->getCost(x, z, *iterX, *iterZ) : MIN_COST;
                                    Cost remain = minCost - (cost + xycost + xzcost);
                                    if (remain > MIN_COST)
                                        extend(getCostYZX, addCostYZX, functionalX, getFunctionX, yz, y, z, x, *iterY, *iterZ, remain);
                                }
                            }
                        }
                    }
                }
                if (xy->connected()) {
                    for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                        Cost xycost = xy->getCost(x, y, *iterX, *iterY);
                        if (xycost > MIN_COST) {
                            Cost costfromxy = MIN_COST;
                            for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
                                Cost cost = getCost(x, y, z, *iterX, *iterY, *iterZ);
                                if (LUBTEST(cost, minCost)) {
                                    Cost xzcost = (xz->connected()) ? xz->getCost(x, z, *iterX, *iterZ) : MIN_COST;
                                    //                                    Cost yzcost = (yz->connected())?yz->getCost(y,z,*iterY,*iterZ):MIN_COST;
                                    //                                    Cost remain = minCost - (cost + xzcost + yzcost);
                                    Cost remain = minCost - (cost + xzcost);
                                    LUB(&costfromxy, remain);
                                }
                            }
                            assert(costfromxy <= xycost);
                            if (costfromxy > MIN_COST) {
                                extend(getCost, addCost, functionalZ, getFunctionZ, xy, x, y, z, *iterX, *iterY, costfromxy);
                            }
                        }
                    }
                }
                if (xz->connected()) {
                    for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
                        Cost xzcost = xz->getCost(x, z, *iterX, *iterZ);
                        if (xzcost > MIN_COST) {
                            Cost costfromxz = MIN_COST;
                            for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                                Cost cost = getCost(x, y, z, *iterX, *iterY, *iterZ);
                                if (LUBTEST(cost, minCost)) {
                                    //                                    Cost yzcost = (yz->connected())?yz->getCost(y,z,*iterY,*iterZ):MIN_COST;
                                    //                                    Cost remain = minCost - (cost + yzcost);
                                    Cost remain = minCost - cost;
                                    LUB(&costfromxz, remain);
                                }
                            }
                            assert(costfromxz <= xzcost);
                            if (costfromxz > MIN_COST) {
                                extend(getCostXZY, addCostXZY, functionalY, getFunctionY, xz, x, z, y, *iterX, *iterZ, costfromxz);
                            }
                        }
                    }
                }
                //                if (yz->connected()) {
                //                    for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                //                        for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
                //                            Cost yzcost = yz->getCost(y,z,*iterY,*iterZ);
                //                            if (yzcost > MIN_COST) {
                //                                Cost cost = getCost(x,y,z,*iterX,*iterY,*iterZ);
                //                                if (LUBTEST(cost, minCost)) {
                //                                    assert(yzcost >= minCost - cost);
                //                                    extend(getCostYZX,addCostYZX,functionalX,getFunctionX,yz,y,z,x,*iterY,*iterZ,minCost - cost);
                //                                }
                //                            }
                //                        }
                //                    }
                //                }
                supportBroken |= project(x, *iterX, minCost, deltaCostsX);
                if (deconnected())
                    return;
            }

            supportX[xindex] = (supportReversed) ? std::make_pair(support.second, support.first) : std::make_pair(support.first, support.second);

            unsigned int yindex = y->toIndex(support.first);
            unsigned int zindex = z->toIndex(support.second);
            if (getIndexX < getIndexZ)
                supportY[yindex] = std::make_pair(*iterX, support.second);
            else
                supportY[yindex] = std::make_pair(support.second, *iterX);
            if (getIndexX < getIndexY)
                supportZ[zindex] = std::make_pair(*iterX, support.first);
            else
                supportZ[zindex] = std::make_pair(support.first, *iterX);
        }
    }
    if (supportBroken) {
        x->findSupport();
    }
}

template <typename T1, typename T2, typename T3, typename T4>
void TernaryConstraint::projectTernaryBinary(T1 getCost, T2 getCostYZX, T3 addCostYZX, bool functionalX, T4 getFunctionX,
    EnumeratedVariable* x, EnumeratedVariable* y, EnumeratedVariable* z, BinaryConstraint* yzin)
{
    // cout << "PROJECT " << *this <<endl;
    // cout << "on " << y->getName() << "," << z->getName() << endl;

    bool flag = false;
    for (EnumeratedVariable::iterator itery = y->begin(); itery != y->end(); ++itery) {
        for (EnumeratedVariable::iterator iterz = z->begin(); iterz != z->end(); ++iterz) {
            Cost mincost = MAX_COST;
            for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
                Cost curcost = getCost(x, y, z, *iterx, *itery, *iterz);
                GLB(&mincost, curcost);
            }
            if (mincost > MIN_COST) {
                flag = true;
                project(getCostYZX, addCostYZX, functionalX, getFunctionX, yzin, y, z, x, *itery, *iterz, mincost);
            }
        }
    }

    assert((yzin == xy || yzin == xz || yzin == yz) && checkTreeDecomposition());

    if (flag) {
        if (y->unassigned() && z->unassigned())
            yzin->reconnect();
        yzin->propagate();
    }
}

template <typename T1, typename T2, typename T3>
bool TernaryConstraint::isEAC(T1 getCostWithBinaries, bool functionalY, T2 getFunctionY, bool functionalZ, T3 getFunctionZ,
    EnumeratedVariable* x, Value a, EnumeratedVariable* y, EnumeratedVariable* z,
    vector<pair<Value, Value>>& supportX)
{
    assert(y->canbe(y->getSupport()));
    assert(y->getCost(y->getSupport()) == MIN_COST);
    assert(z->canbe(z->getSupport()));
    assert(z->getCost(z->getSupport()) == MIN_COST);
    if (getCostWithBinaries(x, y, z, a, y->getSupport(), z->getSupport()) > MIN_COST) {
        if (ToulBar2::FullEAC)
            x->unsetFullEAC();
        unsigned int xindex = x->toIndex(a);
        assert(getIndex(y) < getIndex(z));
        pair<Value, Value> support = supportX[xindex];
        if (y->cannotbe(support.first) || z->cannotbe(support.second) || getCostWithBinaries(x, y, z, a, support.first, support.second) + y->getCost(support.first) + z->getCost(support.second) > MIN_COST) {
            if (functionalZ) {
                for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                    Value valZ = getFunctionZ(x, y, a, *iterY);
                    if (valZ != WRONG_VAL && z->canbe(valZ)) {
                        if (getCostWithBinaries(x, y, z, a, *iterY, valZ) + y->getCost(*iterY) + z->getCost(valZ) == MIN_COST) {
                            supportX[xindex] = std::make_pair(*iterY, valZ);
                            return true;
                        }
                    }
                }
            } else if (functionalY) {
                for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
                    Value valY = getFunctionY(x, z, a, *iterZ);
                    if (valY != WRONG_VAL && y->canbe(valY)) {
                        if (getCostWithBinaries(x, y, z, a, valY, *iterZ) + y->getCost(valY) + z->getCost(*iterZ) == MIN_COST) {
                            supportX[xindex] = std::make_pair(valY, *iterZ);
                            return true;
                        }
                    }
                }
            } else {
                for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                    for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
                        if (getCostWithBinaries(x, y, z, a, *iterY, *iterZ) + y->getCost(*iterY) + z->getCost(*iterZ) == MIN_COST) {
                            supportX[xindex] = std::make_pair(*iterY, *iterZ);
                            return true;
                        }
                    }
                }
            }
            //        cout << x->getName() << " = " << a << " not EAC due to constraint " << *this << endl;
            return false;
        }
    }
    return true;
}

// A triangle of three binary cost functions (maxRPC/PIC)
// class Triangle : public AbstractTernaryConstraint<EnumeratedVariable,EnumeratedVariable,EnumeratedVariable>
//{
//    BinaryConstraint* xy;
//    BinaryConstraint* xz;
//    BinaryConstraint* yz;
//    TernaryConstraint* xyz;
//
// public:
//    Triangle(WCSP *wcsp,
//					  EnumeratedVariable *xx,
//					  EnumeratedVariable *yy,
//					  EnumeratedVariable *zz,
//					  BinaryConstraint* xy,
//					  BinaryConstraint* xz,
//					  BinaryConstraint* yz,
//					  StoreStack<Cost, Cost> *storeCost);
//
//	~Triangle() {}
//
//    bool isTriangle() const {return true;}
//
//    void activateTernary();
//};

#endif /*TB2TERNARYCONSTR_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

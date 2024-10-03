/** \file tb2amongconstr.hpp
 *  \brief Dynamic programming based global cost function : samong_dp
 */

#ifndef TB2AMONGCONSTR_HPP_
#define TB2AMONGCONSTR_HPP_

#include "tb2dpglobalconstr.hpp"

class AmongConstraint : public DPGlobalConstraint {
private:
    template <class Source>
    struct TableCell {
        Cost val;
        Source source;
    };

    typedef TableCell<int> DPTableCell;
    DPTableCell** f;
    DPTableCell** invf;
    DPTableCell** curf;
    Cost top;

    typedef TableCell<Value> UnaryTableCell;
    UnaryTableCell *minBarU, *minU;

    template <class T>
    void resizeTable(T**& table, int width, int heigth)
    {
        assert(width >= arity() + 1);
        table = new T*[width];
        for (int i = 0; i <= arity(); i++) {
            table[i] = new T[heigth];
        }
    }

    template <class T>
    void deleteTable(T**& table)
    {
        for (int i = 0; i <= arity(); i++)
            delete[] table[i];
        delete[] table;
        table = NULL;
    }

    set<Value> V;
    int ub, lb;

    void recomputeTable(DPTableCell** table, DPTableCell** invTable = NULL, int startRow = 0);
    void recompute();

    Cost computeMinU(int var);

    Cost computeMinBarU(int var);

protected:
    Cost minCostOriginal();
    Cost minCostOriginal(int var, Value val, bool changed);
    Result minCost(int var, Value val, bool changed);

public:
    AmongConstraint(WCSP* wcsp, EnumeratedVariable** scope, int arity);
    virtual ~AmongConstraint();

    Cost evalOriginal(const Tuple& s);

    void read(istream& file, bool mult = true);
    void setUpperBound(int upper) { ub = upper; }
    void setLowerBound(int lower) { lb = lower; }
    void addBoundingValue(Value value) { V.insert(value); }
    virtual void initMemoization();

    string getName()
    {
        string name = "samong";
        name += to_string("_") + to_string(lb) + to_string("_") + to_string(ub) + to_string("_") + to_string(V.size());
        for (set<Value>::iterator iter = V.begin(); iter != V.end(); ++iter)
            name += to_string("_") + to_string(*iter);
        return name;
    }
    void dump(ostream& os, bool original = true);
};

#endif /*TB2AMONGCONSTR_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

/** \file tb2globaldecomposable.hpp
 *  \brief Decomposable global cost functions : WeightedRegular, WeightedAmong
 */
#ifndef TB2GLOBALDECOMPOSABLE_HPP_
#define TB2GLOBALDECOMPOSABLE_HPP_

#include "tb2wcsp.hpp"
#include "tb2types.hpp"
#include "tb2enumvar.hpp"
#include "utils/tb2automaton.hpp"

class DecomposableGlobalCostFunction {
protected:
    int arity;
    int* scope;
    string label;

public:
    DecomposableGlobalCostFunction();
    DecomposableGlobalCostFunction(unsigned int _arity, int* _scope);
    virtual ~DecomposableGlobalCostFunction();
    static std::unique_ptr<DecomposableGlobalCostFunction> FactoryDGCF(string type, unsigned int _arity, int* _scope, istream& file, bool mult = true);

    int getArity() { return arity; }
    int getVarIndex(int i) { return scope[i]; }
    void setLabel(string _label) { label = _label; }
    string getLabel() { return label; }

    virtual Cost evaluate(Value* tuple) = 0;
    virtual void addToCostFunctionNetwork(WCSP* wcsp) = 0;
    virtual void display() = 0;

    void color(int);
};

class WeightedAmong : public DecomposableGlobalCostFunction {
private:
    set<Value> values;
    string semantics;
    Cost baseCost;
    unsigned int lb;
    unsigned int ub;

public:
    WeightedAmong();
    WeightedAmong(unsigned int _arity, int* _scope);
    WeightedAmong(unsigned int _arity, int* _scope, istream& file, bool mult = true);
    ~WeightedAmong();

    inline void addValue(Value _value) { values.insert(_value); }
    inline void setSemantics(string _semantics) { semantics = _semantics; }
    inline void setBaseCost(Cost _baseCost) { baseCost = _baseCost; }
    inline void setBounds(unsigned int _lb, unsigned int _ub)
    {
        lb = _lb;
        ub = _ub;
    }

    Cost evaluate(Value* tuple);
    void addToCostFunctionNetwork(WCSP* wcsp);
    void display();
};

class WeightedRegular : public DecomposableGlobalCostFunction {
private:
    WFA* automaton;

public:
    WeightedRegular();
    WeightedRegular(unsigned int _arity, int* _scope);
    WeightedRegular(unsigned int _arity, int* _scope, istream& file, bool mult = true);
    ~WeightedRegular();

    inline void setWFA(WFA* _automaton) { automaton = _automaton; }

    Cost evaluate(Value* tuple)
    {
        cout << "WeightedRegular::evaluate => no yet implemented" << endl;
        return 0;
    }
    void addToCostFunctionNetwork(WCSP* wcsp);
    void display();
};

class WeightedSum : public DecomposableGlobalCostFunction {
private:
    string comparator;
    int rightRes;
    string semantics;
    Cost baseCost;

public:
    WeightedSum();
    WeightedSum(unsigned int _arity, int* _scope);
    WeightedSum(unsigned int _arity, int* _scope, istream& file, bool mult = true);
    ~WeightedSum();

    inline void setSemantics(string _semantics) { semantics = _semantics; }
    inline void setBaseCost(Cost _baseCost) { baseCost = _baseCost; }
    inline void setComparator(string _comparator) { comparator = _comparator; }
    inline void setRightRes(int _rightRes) { rightRes = _rightRes; }

    Cost evaluate(Value* tuple);
    void addToCostFunctionNetwork(WCSP* wcsp);
    void display();
};

class WeightedOverlap : public DecomposableGlobalCostFunction {
private:
    string semantics;
    Cost baseCost;
    string comparator;
    int rightRes;

public:
    WeightedOverlap();
    WeightedOverlap(unsigned int _arity, int* _scope);
    WeightedOverlap(unsigned int _arity, int* _scope, istream& file, bool mult = true);
    ~WeightedOverlap();

    inline void setSemantics(string _semantics) { semantics = _semantics; }
    inline void setBaseCost(Cost _baseCost) { baseCost = _baseCost; }
    inline void setComparator(string _comparator) { comparator = _comparator; }
    inline void setRightRes(int _rightRes) { rightRes = _rightRes; }

    Cost evaluate(Value* tuple);
    void addToCostFunctionNetwork(WCSP* wcsp);
    void display();
};

// WeightedAmong with a Variable as right member
class WeightedVarAmong : public DecomposableGlobalCostFunction {
private:
    set<Value> values;
    string semantics;
    Cost baseCost;
    unsigned int index;

public:
    WeightedVarAmong();
    WeightedVarAmong(unsigned int _arity, int* _scope);
    WeightedVarAmong(unsigned int _arity, int* _scope, istream& file, bool mult = true);
    ~WeightedVarAmong();

    inline void addValue(Value _value) { values.insert(_value); }
    inline void setSemantics(string _semantics) { semantics = _semantics; }
    inline void setBaseCost(Cost _baseCost) { baseCost = _baseCost; }
    inline void setIndex(unsigned int _index) { index = _index; }

    Cost evaluate(Value* tuple)
    {
        cerr << "Not yet implemented" << endl;
        return 0;
    }
    void addToCostFunctionNetwork(WCSP* wcsp);
    void display();
};

class WeightedVarSum : public DecomposableGlobalCostFunction {
private:
    string comparator;
    unsigned int index;
    string semantics;
    Cost baseCost;

public:
    WeightedVarSum();
    WeightedVarSum(unsigned int _arity, int* _scope);
    WeightedVarSum(unsigned int _arity, int* _scope, istream& file, bool mult = true);
    ~WeightedVarSum();

    inline void setSemantics(string _semantics) { semantics = _semantics; }
    inline void setBaseCost(Cost _baseCost) { baseCost = _baseCost; }
    inline void setComparator(string _comparator) { comparator = _comparator; }
    inline void setIndex(unsigned int _index) { index = _index; }

    Cost evaluate(Value* tuple)
    {
        cerr << "Not yet implemented" << endl;
        return 0;
    }
    void addToCostFunctionNetwork(WCSP* wcsp);
    void display();
};

// WeightedAllDifferent using several WeightedAmongs
class WeightedAllDifferent : public DecomposableGlobalCostFunction {
private:
    string semantics;
    Cost baseCost;

public:
    WeightedAllDifferent();
    WeightedAllDifferent(unsigned int _arity, int* _scope);
    WeightedAllDifferent(unsigned int _arity, int* _scope, istream& file, bool mult = true);
    ~WeightedAllDifferent();

    inline void setSemantics(string _semantics) { semantics = _semantics; }
    inline void setBaseCost(Cost _baseCost) { baseCost = _baseCost; }

    Cost evaluate(Value* tuple)
    {
        cerr << "Not yet implemented" << endl;
        return 0;
    }
    void addToCostFunctionNetwork(WCSP* wcsp);
    void display();
};

// WeightedGcc using several WeightedAmongs
class WeightedGcc : public DecomposableGlobalCostFunction {
private:
    string semantics;
    Cost baseCost;
    int nbValueToWatch;
    map<Value, pair<unsigned int, unsigned int>> bounds;

public:
    WeightedGcc();
    WeightedGcc(unsigned int _arity, int* _scope);
    WeightedGcc(unsigned int _arity, int* _scope, istream& file, bool mult = true);
    ~WeightedGcc();

    inline void setSemantics(string _semantics) { semantics = _semantics; }
    inline void setBaseCost(Cost _baseCost) { baseCost = _baseCost; }
    inline void setNbValue(int _nbValue) { nbValueToWatch = _nbValue; }
    void setBounds(Value value, unsigned int lb, unsigned int ub);

    Cost evaluate(Value* tuple)
    {
        cerr << "Not yet implemented" << endl;
        return 0;
    }
    void rec_sum_counters(WCSP* wcsp, int* scope, int arity, int totlb, int totub, int* counters, int* lb, int* ub, int nb, int rec);
    void addToCostFunctionNetwork(WCSP* wcsp);
    void display();
};

// WeightedSame using several WeightedAmongs
class WeightedSame : public DecomposableGlobalCostFunction {
private:
    string semantics;
    Cost baseCost;

public:
    WeightedSame();
    WeightedSame(unsigned int _arity, int* _scope);
    WeightedSame(unsigned int _arity, int* _scope, istream& file, bool mult = true);
    ~WeightedSame();

    inline void setSemantics(string _semantics) { semantics = _semantics; }
    inline void setBaseCost(Cost _baseCost) { baseCost = _baseCost; }

    Cost evaluate(Value* tuple)
    {
        cerr << "Not yet implemented" << endl;
        return 0;
    }
    void addToCostFunctionNetwork(WCSP* wcsp);
    void display();
};

// WeightedSameGcc using several WeightedAmongs and Unary
class WeightedSameGcc : public DecomposableGlobalCostFunction {
private:
    string semantics;
    Cost baseCost;
    int nbValueToWatch;
    map<Value, pair<unsigned int, unsigned int>> bounds;

public:
    WeightedSameGcc();
    WeightedSameGcc(unsigned int _arity, int* _scope);
    WeightedSameGcc(unsigned int _arity, int* _scope, istream& file, bool mult = true);
    ~WeightedSameGcc();

    inline void setSemantics(string _semantics) { semantics = _semantics; }
    inline void setBaseCost(Cost _baseCost) { baseCost = _baseCost; }
    inline void setNbValue(int _nbValue) { nbValueToWatch = _nbValue; }
    void setBounds(Value value, unsigned int lb, unsigned int ub);

    Cost evaluate(Value* tuple)
    {
        cerr << "Not yet implemented" << endl;
        return 0;
    }
    void addToCostFunctionNetwork(WCSP* wcsp);
    void display();
};

class WeightedDiverse : public DecomposableGlobalCostFunction {
private:
    vector<Value> values;
    string semantics;
    Cost baseCost;
    unsigned int distance;
    int method; // Dual=0, Hidden=1, Ternary=2

public:
    WeightedDiverse();
    WeightedDiverse(unsigned int _arity, int* _scope);
    WeightedDiverse(unsigned int _arity, int* _scope, istream& file, int _method, bool mult = true);
    ~WeightedDiverse();

    inline void addValue(Value _value) { values.push_back(_value); }
    inline void setSemantics(string _semantics) { semantics = _semantics; }
    inline void setBaseCost(Cost _baseCost) { baseCost = _baseCost; }
    inline void setDistance(unsigned int dist) { distance = dist; }
    inline void setMethod(int m) { method = m; }

    Cost evaluate(Value* tuple);
    void addToCostFunctionNetwork(WCSP* wcsp);
    void display();
};

#endif

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

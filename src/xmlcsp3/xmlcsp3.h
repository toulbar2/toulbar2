/*
 * xcsp3.h
 *
 *  Created on: May 20th 2022
 *      Author: sdegivry
 */

#ifndef SRC_XCSP3_XMLCSP3_H_
#define SRC_XCSP3_XMLCSP3_H_

#include "XCSP3CoreParser.h"
#include "XCSP3CoreCallbacks.h"
#include "XCSP3Variable.h"

#include "toulbar2lib.hpp"

using namespace XCSP3Core;

#define MAX_COST ((Cost)INT_MAX * 1024)

class MySolverCallbacks : public XCSP3CoreCallbacks {
    public:
        WeightedCSP *problem;
        std::map<string,int> mapping;
        vector<vector<int> > lastTuples;

    virtual ~MySolverCallbacks() {}

    void beginInstance(InstanceType type) {
        XCSP3CoreCallbacks::intensionUsingString = false;
        XCSP3CoreCallbacks::recognizeSpecialIntensionCases = true;
        XCSP3CoreCallbacks::recognizeSpecialCountCases = false;
        XCSP3CoreCallbacks::recognizeNValuesCases = false;
        problem->updateUb(MAX_COST);
    }

    void endInstance() {
        problem->sortConstraints();
    }

    void buildVariableInteger(string id, int minValue, int maxValue) override {
        int v = problem->makeEnumeratedVariable(id, minValue, maxValue);
        mapping[id] = v;
    }

    void buildVariableInteger(string id, vector<int> &values) override {
        int v = problem->makeEnumeratedVariable(id,values);
        mapping[id] = v;
    }

    int getMyVar(string var) {
        assert(mapping.find(var) != mapping.end());
        assert(mapping[var] >= 0 && mapping[var] < (int)problem->numberOfVariables());
        return mapping[var];
    }

    int getMyVar(XVariable *var) {
        return getMyVar(var->id);
    }

    // transforms a vector of XVariable in vector of toulbar2 variable indices and add it to dest
    void toMyVariables(vector<XVariable*> &src, vector<int> &dest) {
        set<int> control;
        for(unsigned int i = 0;i<src.size();i++) {
            dest.push_back(getMyVar(src[i]));
            control.insert(getMyVar(src[i]));
        }
        assert(dest.size() == control.size());
        assert(dest.size() > 0);
    }

    // transforms a vector of string Variable name in vector of toulbar2 variable indices and add it to dest
    void toMyVariables(vector<string> &src, vector<int> &dest) {
        set<int> control;
        for(unsigned int i = 0;i<src.size();i++) {
            dest.push_back(getMyVar(src[i]));
            control.insert(getMyVar(src[i]));
        }
        assert(dest.size() == control.size());
        assert(dest.size() > 0);
    }

    void recursiveExpanse(vector<int> &vars, vector<int> &tuple, unsigned int pos, int value, vector<vector<int>> &tuples, vector<int> expanded) {
        if (pos < vars.size()) {
            if (value == STAR) {
                for (unsigned int a = 0; a < problem->getDomainInitSize(vars[pos]); a++) {
                    recursiveExpanse(vars, tuple, pos, problem->toValue(vars[pos], a), tuples, expanded);
                }
            } else {
                expanded.push_back(value);
                recursiveExpanse(vars, tuple, pos+1, (pos+1 < vars.size())?tuple[pos+1]:STAR, tuples, expanded);
            }
        } else {
            assert(expanded.size() == tuple.size());
            tuples.push_back(expanded);
        }
    }

    void buildConstraintExtension(vector<int> vars, vector<vector<int> > &tuples, bool isSupport, bool hasStar) {
        if(hasStar) {
            vector<vector<int> > newtuples;
            for(auto& tuple:tuples) {
                recursiveExpanse(vars, tuple, 0, tuple[0], newtuples, vector<int>());
            }
            assert(newtuples.size() >= tuples.size());
            tuples = newtuples; // warning it will change the input tuples
        }
        lastTuples = tuples; // copy tuples for future usage by buildConstraintExtensionAs after tuples expansion
        assert(vars.size() > 0);
        if (vars.size()==1) {
            vector<Cost> costs(problem->getDomainInitSize(vars[0]), (isSupport)?MAX_COST:MIN_COST);
            for(auto& tuple:tuples) {
                if (isSupport) {
                    costs[problem->toIndex(vars[0], tuple[0])] = MIN_COST;
                } else {
                    costs[problem->toIndex(vars[0], tuple[0])] = MAX_COST;
                }
            }
            problem->postUnaryConstraint(vars[0], costs);
        } else if (vars.size()==2) {
            vector<Cost> costs(problem->getDomainInitSize(vars[0]) * problem->getDomainInitSize(vars[1]), (isSupport)?MAX_COST:MIN_COST);
            for(auto& tuple:tuples) {
                if (isSupport) {
                    costs[problem->toIndex(vars[0], tuple[0]) * problem->getDomainInitSize(vars[1]) + problem->toIndex(vars[1], tuple[1])] = MIN_COST;
                } else {
                    costs[problem->toIndex(vars[0], tuple[0]) * problem->getDomainInitSize(vars[1]) + problem->toIndex(vars[1], tuple[1])] = MAX_COST;
                }
            }
            problem->postBinaryConstraint(vars[0], vars[1], costs);
        } else if (vars.size()==3) {
            vector<Cost> costs(problem->getDomainInitSize(vars[0]) * problem->getDomainInitSize(vars[1]) * problem->getDomainInitSize(vars[2]), (isSupport)?MAX_COST:MIN_COST);
            for(auto& tuple:tuples) {
                if (isSupport) {
                    costs[problem->toIndex(vars[0], tuple[0]) * problem->getDomainInitSize(vars[1]) * problem->getDomainInitSize(vars[2]) + problem->toIndex(vars[1], tuple[1]) * problem->getDomainInitSize(vars[2]) + problem->toIndex(vars[2], tuple[2])] = MIN_COST;
                } else {
                    costs[problem->toIndex(vars[0], tuple[0]) * problem->getDomainInitSize(vars[1]) * problem->getDomainInitSize(vars[2]) + problem->toIndex(vars[1], tuple[1]) * problem->getDomainInitSize(vars[2]) + problem->toIndex(vars[2], tuple[2])] = MAX_COST;
                }
            }
            problem->postTernaryConstraint(vars[0], vars[1], vars[2], costs);
        } else {
            int ctridx = problem->postNaryConstraintBegin(vars, (isSupport)?MAX_COST:MIN_COST, tuples.size());
            for(auto& tuple:tuples) {
                if (isSupport) {
                    problem->postNaryConstraintTuple(ctridx, tuple, MIN_COST);
                } else {
                    problem->postNaryConstraintTuple(ctridx, tuple, MAX_COST);
                }
            }
            problem->postNaryConstraintEnd(ctridx);
        }
    }

    void buildConstraintExtension(string id, vector<string> list, vector<vector<int> > &tuples, bool isSupport, bool hasStar) {
        vector<int> vars;
        toMyVariables(list,vars);
        buildConstraintExtension(vars, tuples, isSupport, hasStar);
    }

    void buildConstraintExtension(string id, vector<XVariable *> list, vector<vector<int> > &tuples, bool isSupport, bool hasStar) override {
        vector<string> thelist;
        for (XVariable *var : list) {
            thelist.push_back(var->id);
        }
        buildConstraintExtension(id, thelist, tuples, isSupport, hasStar);
    }

    void buildConstraintExtension(int var, vector<int> &tuples, bool isSupport, bool hasStar) {
        // This function is called for unary extensional constraint.
        assert(tuples[0] != STAR || tuples.size()==1);
        lastTuples = vector<vector<int> >();
        lastTuples.push_back(tuples);
        if (tuples[0] == STAR) {
            if (!isSupport) {
                throw Contradiction();
            } else {
                return;
            }
        }
        vector<Cost> costs(problem->getDomainInitSize(var), (isSupport)?MAX_COST:MIN_COST);
        for(auto value:tuples) {
            assert(value != STAR);
            if (isSupport) {
                costs[problem->toIndex(var, value)] = MIN_COST;
            } else {
                costs[problem->toIndex(var, value)] = MAX_COST;
            }
        }
        problem->postUnaryConstraint(var, costs);
    }

    void buildConstraintExtension(string id, string variable, vector<int> &tuples, bool isSupport, bool hasStar) {
        int var = getMyVar(variable);
        buildConstraintExtension(var, tuples, isSupport, hasStar);
    }

    void buildConstraintExtension(string id, XVariable *variable, vector<int> &tuples, bool isSupport, bool hasStar) override {
        buildConstraintExtension(id, variable->id, tuples, isSupport, hasStar);
    }

    // This function is called with group of constraint where the set of tuples is exactly the same
    // than the previous one (then, you can save time/memory using the same set of tuples.
    void buildConstraintExtensionAs(string id, vector<XVariable *> list, bool isSupport, bool hasStar) override {
        buildConstraintExtension(id, list, lastTuples, isSupport, false);
    }

    // returns a WCSP variable index corresponding to the evaluation of Tree expression
    int buildConstraintIntension(Tree *tree) {
        set<int> values;
        vector<string> list = tree->listOfVariables;
        vector<int> vars;
        toMyVariables(list,vars);
        vector<vector<int> > tuples;
        map<string, int> tuplemap;
        vector<int> tuplevec;
        unsigned int pos = 0;
        for (unsigned int i=0; i<vars.size(); i++) {
            Value val = problem->getInf(vars[i]);
            tuplemap[list[i]] = val;
            tuplevec.push_back(val);
        }
        tuplevec.push_back(0); // dummy evaluation result associated to reified variable
        do {
            int res = tree->evaluate(tuplemap);
            values.insert(res);
            tuplevec[vars.size()] = res; // update corresponding value assigned to reified variable
            tuples.push_back(tuplevec);
            while (pos < vars.size() && tuplevec[pos] == problem->getSup(vars[pos])) {
                pos++;
            }
            if (pos < vars.size()) {
                Value nextval = problem->nextValue(vars[pos], tuplevec[pos]);
                assert(nextval != tuplevec[pos]);
                tuplemap[list[pos]] = nextval;
                tuplevec[pos] = nextval;
                for (unsigned int i=0; i<pos; i++) {
                    Value val = problem->getInf(vars[i]);
                    tuplemap[list[i]] = val;
                    tuplevec[i] = val;
                }
                pos = 0;
            }
        } while (pos < vars.size());
        string extravarname = IMPLICIT_VAR_TAG + "expr" + to_string(problem->numberOfVariables());
        vector<int> vvalues(values.begin(), values.end());
        int extravar = problem->makeEnumeratedVariable(extravarname, vvalues);
        mapping[extravarname] = extravar;
        vars.push_back(extravar); // add reified variable in the scope
        buildConstraintExtension(vars, tuples, true, false);
        return extravar;
    }

    void buildConstraintIntension(string id, Tree *tree) override {
        vector<string> list = tree->listOfVariables;
        vector<int> vars;
        toMyVariables(list,vars);
        vector<vector<int> > tuples;
        map<string, int> tuplemap;
        vector<int> tuplevec;
        unsigned int pos = 0;
        for (unsigned int i=0; i<vars.size(); i++) {
            Value val = problem->getInf(vars[i]);
            tuplemap[list[i]] = val;
            tuplevec.push_back(val);
        }
        do {
            if (tree->evaluate(tuplemap)) {
                tuples.push_back(tuplevec);
            }
            while (pos < vars.size() && tuplevec[pos] == problem->getSup(vars[pos])) {
                pos++;
            }
            if (pos < vars.size()) {
                Value nextval = problem->nextValue(vars[pos], tuplevec[pos]);
                assert(nextval != tuplevec[pos]);
                tuplemap[list[pos]] = nextval;
                tuplevec[pos] = nextval;
                for (unsigned int i=0; i<pos; i++) {
                    Value val = problem->getInf(vars[i]);
                    tuplemap[list[i]] = val;
                    tuplevec[i] = val;
                }
                pos = 0;
            }
        } while (pos < vars.size());
        buildConstraintExtension(id, list, tuples, true, false);
    }

    void buildConstraintPrimitive(OrderType op, int varx, int k, int vary) {
        assert(varx != vary);
        vector<Cost> costs;
        for (unsigned int a = 0; a < problem->getDomainInitSize(varx); a++) {
            for (unsigned int b = 0; b < problem->getDomainInitSize(vary); b++) {
                switch (op) {
                case OrderType::LE:
                    costs.push_back((problem->toValue(varx, a) + k <= problem->toValue(vary, b))?MIN_COST:MAX_COST);
                    break;
                case OrderType::LT:
                    costs.push_back((problem->toValue(varx, a) + k < problem->toValue(vary, b))?MIN_COST:MAX_COST);
                    break;
                case OrderType::GE:
                    costs.push_back((problem->toValue(varx, a) + k >= problem->toValue(vary, b))?MIN_COST:MAX_COST);
                    break;
                case OrderType::GT:
                    costs.push_back((problem->toValue(varx, a) + k > problem->toValue(vary, b))?MIN_COST:MAX_COST);
                    break;
                case OrderType::IN:
                case OrderType::EQ:
                    costs.push_back((problem->toValue(varx, a) + k == problem->toValue(vary, b))?MIN_COST:MAX_COST);
                    break;
                case OrderType::NE:
                    costs.push_back((problem->toValue(varx, a) + k != problem->toValue(vary, b))?MIN_COST:MAX_COST);
                    break;
                default:
                    cerr << "Sorry operator " << op << " not implemented!" << endl;
                    throw WrongFileFormat();
                }
            }
        }
        problem->postBinaryConstraint(varx, vary, costs);
    }

    void buildConstraintPrimitiveMult(OrderType op, int varx, int mult, int vary) {
        assert(varx != vary);
        vector<Cost> costs;
        for (unsigned int a = 0; a < problem->getDomainInitSize(varx); a++) {
            for (unsigned int b = 0; b < problem->getDomainInitSize(vary); b++) {
                switch (op) {
                case OrderType::LE:
                    costs.push_back((problem->toValue(varx, a) * mult <= problem->toValue(vary, b))?MIN_COST:MAX_COST);
                    break;
                case OrderType::LT:
                    costs.push_back((problem->toValue(varx, a) * mult < problem->toValue(vary, b))?MIN_COST:MAX_COST);
                    break;
                case OrderType::GE:
                    costs.push_back((problem->toValue(varx, a) * mult >= problem->toValue(vary, b))?MIN_COST:MAX_COST);
                    break;
                case OrderType::GT:
                    costs.push_back((problem->toValue(varx, a) * mult > problem->toValue(vary, b))?MIN_COST:MAX_COST);
                    break;
                case OrderType::IN:
                case OrderType::EQ:
                    costs.push_back((problem->toValue(varx, a) * mult == problem->toValue(vary, b))?MIN_COST:MAX_COST);
                    break;
                case OrderType::NE:
                    costs.push_back((problem->toValue(varx, a) * mult != problem->toValue(vary, b))?MIN_COST:MAX_COST);
                    break;
                default:
                    cerr << "Sorry operator " << op << " not implemented!" << endl;
                    throw WrongFileFormat();
                }
            }
        }
        problem->postBinaryConstraint(varx, vary, costs);
    }

    void buildConstraintPrimitive(string id, OrderType op, string x, int k, string y)  {
        int varx = getMyVar(x);
        int vary = getMyVar(y);
        buildConstraintPrimitive(op, varx, k, vary);
    }

    void buildConstraintPrimitive(string id, OrderType op, XVariable *x, int k, XVariable *y) override {
        buildConstraintPrimitive(id, op, x->id, k, y->id);
    }

    void buildConstraintPrimitive(OrderType op, int varx, int k) {
        vector<Cost> costs;
        for (unsigned int a = 0; a < problem->getDomainInitSize(varx); a++) {
            switch (op) {
            case OrderType::LE:
                costs.push_back((problem->toValue(varx, a) <= k)?MIN_COST:MAX_COST);
                break;
            case OrderType::LT:
                costs.push_back((problem->toValue(varx, a) < k)?MIN_COST:MAX_COST);
                break;
            case OrderType::GE:
                costs.push_back((problem->toValue(varx, a) >= k)?MIN_COST:MAX_COST);
                break;
            case OrderType::GT:
                costs.push_back((problem->toValue(varx, a) > k)?MIN_COST:MAX_COST);
                break;
            case OrderType::IN:
            case OrderType::EQ:
                costs.push_back((problem->toValue(varx, a) == k)?MIN_COST:MAX_COST);
                break;
            case OrderType::NE:
                costs.push_back((problem->toValue(varx, a) != k)?MIN_COST:MAX_COST);
                break;
            default:
                cerr << "Sorry operator " << op << " not implemented!" << endl;
                throw WrongFileFormat();
            }
        }
        problem->postUnaryConstraint(varx, costs);
    }

    void buildConstraintPrimitive(string id, OrderType op, string x, int k) {
        int varx = getMyVar(x);
        buildConstraintPrimitive(op, varx, k);
    }

    void buildConstraintPrimitive(string id, OrderType op, XVariable *x, int k) override {
        buildConstraintPrimitive(id, op, x->id, k);
    }

    void buildConstraintPrimitive(int varx, bool in, int min, int max) {
        vector<Cost> costs;
        for (unsigned int a = 0; a < problem->getDomainInitSize(varx); a++) {
            if (in) {
                costs.push_back((problem->toValue(varx, a) >= min && problem->toValue(varx, a) <= max)?MIN_COST:MAX_COST);
            } else {
                costs.push_back((problem->toValue(varx, a) < min || problem->toValue(varx, a) > max)?MIN_COST:MAX_COST);
            }
        }
        problem->postUnaryConstraint(varx, costs);
    }

    void buildConstraintPrimitive(string id, string x, bool in, int min, int max) {
        int varx = getMyVar(x);
        buildConstraintPrimitive(varx, in, min, max);
    }

    void buildConstraintPrimitive(string id, XVariable *x, bool in, int min, int max) override {
        buildConstraintPrimitive(id, x->id, in, min, max);
    }

    void buildConstraintMult(int varx, int vary, int varz) {
        vector<Cost> costs;
        if (varx != vary) {
            if (varx != varz && vary != varz) {
                for (unsigned int a = 0; a < problem->getDomainInitSize(varx); a++) {
                    for (unsigned int b = 0; b < problem->getDomainInitSize(vary); b++) {
                        for (unsigned int c = 0; c < problem->getDomainInitSize(varz); c++) {
                            costs.push_back((problem->toValue(varx, a) * problem->toValue(vary, b) == problem->toValue(varz, c))?MIN_COST:MAX_COST);
                        }
                    }
                }
                problem->postTernaryConstraint(varx, vary, varz, costs);
            } else if (varx == varz) {
                for (unsigned int a = 0; a < problem->getDomainInitSize(varx); a++) {
                    for (unsigned int b = 0; b < problem->getDomainInitSize(vary); b++) {
                        costs.push_back((problem->toValue(varx, a) * problem->toValue(vary, b) == problem->toValue(varx, a))?MIN_COST:MAX_COST);
                    }
                }
                problem->postBinaryConstraint(varx, vary, costs);
            } else {
                assert(vary == varz);
                for (unsigned int a = 0; a < problem->getDomainInitSize(varx); a++) {
                    for (unsigned int b = 0; b < problem->getDomainInitSize(vary); b++) {
                        costs.push_back((problem->toValue(varx, a) * problem->toValue(vary, b) == problem->toValue(vary, b))?MIN_COST:MAX_COST);
                    }
                }
                problem->postBinaryConstraint(varx, vary, costs);
            }
        } else {
            assert(varx == vary);
            if (varx != varz) {
                for (unsigned int a = 0; a < problem->getDomainInitSize(varx); a++) {
                    for (unsigned int c = 0; c < problem->getDomainInitSize(varz); c++) {
                        costs.push_back((problem->toValue(varx, a) * problem->toValue(varx, a) == problem->toValue(varz, c))?MIN_COST:MAX_COST);
                    }
                }
                problem->postBinaryConstraint(varx, varz, costs);
            } else {
                assert(varx == varz);
                for (unsigned int a = 0; a < problem->getDomainInitSize(varx); a++) {
                    costs.push_back((problem->toValue(varx, a) * problem->toValue(varx, a) == problem->toValue(varx, a))?MIN_COST:MAX_COST);
                }
                problem->postUnaryConstraint(varx, costs);
            }
        }
    }

    void buildConstraintMult(string id, string x, string y, string z) {
        int varx = getMyVar(x);
        int vary = getMyVar(y);
        int varz = getMyVar(z);
        buildConstraintMult(varx, vary, varz);
    }

    void buildConstraintMult(string id, XVariable *x, XVariable *y, XVariable *z) override {
        buildConstraintMult(id, x->id, y->id, z->id);
    }

    void buildConstraintAlldifferent(vector<int> &vars) {
        if (vars.size() <= 50) {
            for (unsigned int i = 0; i < vars.size(); i++) {
                for (unsigned int j = i+1; j < vars.size(); j++) {
                    vector<Cost> costs;
                    for (unsigned int a = 0; a < problem->getDomainInitSize(vars[i]); a++) {
                        for (unsigned int b = 0; b < problem->getDomainInitSize(vars[j]); b++) {
                            costs.push_back((problem->toValue(vars[i], a)!=problem->toValue(vars[j], b))?MIN_COST:MAX_COST);
                        }
                    }
                    problem->postBinaryConstraint(vars[i], vars[j], costs);
                }
            }
        } else {
            problem->postWAllDiff(vars, "hard", "knapsack", MAX_COST);
        }
    }

    void buildConstraintAlldifferent(string id, vector<string> &list) {
        vector<int> vars;
        toMyVariables(list,vars);
        buildConstraintAlldifferent(vars);
    }

    void buildConstraintAlldifferent(string id, vector<XVariable *> &list) override {
        vector<string> slist;
        for (auto& x:list) {
            slist.push_back(x->id);
        }
        buildConstraintAlldifferent(id, slist);
    }

    void buildConstraintAlldifferent(string id, vector<Tree *> &list) override {
        vector<int> vars;
        for (unsigned int i = 0; i < list.size(); i++) {
            vars.push_back(buildConstraintIntension(list[i]));
        }
        assert(vars.size() == list.size());
        buildConstraintAlldifferent(vars);
    }

    void buildConstraintSum(vector<int> &vars, vector<int> &coefs, XCondition &cond) {
        assert(vars.size() == coefs.size());
        vector<int> myvars = vars;
        vector<int> mycoefs = coefs;
        int rightcoef = 0;
        string extravarname = "";
        int extravar = -1;
        string params = "";
        switch (cond.operandType) {
        case OperandType::VARIABLE:
            myvars.push_back(getMyVar(cond.var));
            mycoefs.push_back(-1);
            rightcoef -= cond.val;
        case OperandType::INTEGER:
            rightcoef += cond.val;
            switch (cond.op) {
            case OrderType::LE:
                params = to_string(-rightcoef);
                for (unsigned int i=0; i < myvars.size(); i++) {
                    int domsize = problem->getDomainInitSize(myvars[i]);
                    params += " " + to_string(domsize);
                    for (int idval=0; idval < domsize; idval++) {
                        int value = problem->toValue(myvars[i], idval);
                        params += " " + to_string(value) + " " + to_string(-mycoefs[i] * value);
                    }
                }
                problem->postKnapsackConstraint(myvars, params, false, true, false);
                break;
            case OrderType::LT:
                params = to_string(-rightcoef + 1);
                for (unsigned int i=0; i < myvars.size(); i++) {
                    int domsize = problem->getDomainInitSize(myvars[i]);
                    params += " " + to_string(domsize);
                    for (int idval=0; idval < domsize; idval++) {
                        int value = problem->toValue(myvars[i], idval);
                        params += " " + to_string(value) + " " + to_string(-mycoefs[i] * value);
                    }
                }
                problem->postKnapsackConstraint(myvars, params, false, true, false);
                break;
            case OrderType::GE:
                params = to_string(rightcoef);
                for (unsigned int i=0; i < myvars.size(); i++) {
                    int domsize = problem->getDomainInitSize(myvars[i]);
                    params += " " + to_string(domsize);
                    for (int idval=0; idval < domsize; idval++) {
                        int value = problem->toValue(myvars[i], idval);
                        params += " " + to_string(value) + " " + to_string(mycoefs[i] * value);
                    }
                }
                problem->postKnapsackConstraint(myvars, params, false, true, false);
                break;
            case OrderType::GT:
                params = to_string(rightcoef + 1);
                for (unsigned int i=0; i < myvars.size(); i++) {
                    int domsize = problem->getDomainInitSize(myvars[i]);
                    params += " " + to_string(domsize);
                    for (int idval=0; idval < domsize; idval++) {
                        int value = problem->toValue(myvars[i], idval);
                        params += " " + to_string(value) + " " + to_string(mycoefs[i] * value);
                    }
                }
                problem->postKnapsackConstraint(myvars, params, false, true, false);
                break;
            case OrderType::NE:
                extravarname = IMPLICIT_VAR_TAG + "sum" + to_string(problem->numberOfVariables());
                extravar = problem->makeEnumeratedVariable(extravarname, 0, 1);
                mapping[extravarname] = extravar;
                myvars.push_back(extravar);
                mycoefs.push_back(INT_MAX);
                params = to_string(rightcoef + 1);
                for (unsigned int i=0; i < myvars.size(); i++) {
                    int domsize = problem->getDomainInitSize(myvars[i]);
                    params += " " + to_string(domsize);
                    for (int idval=0; idval < domsize; idval++) {
                        int value = problem->toValue(myvars[i], idval);
                        params += " " + to_string(value) + " " + to_string(mycoefs[i] * value);
                    }
                }
                problem->postKnapsackConstraint(myvars, params, false, true, false);
                params = to_string((Long)-rightcoef + 1 - INT_MAX);
                for (unsigned int i=0; i < myvars.size(); i++) {
                    int domsize = problem->getDomainInitSize(myvars[i]);
                    params += " " + to_string(domsize);
                    for (int idval=0; idval < domsize; idval++) {
                        int value = problem->toValue(myvars[i], idval);
                        params += " " + to_string(value) + " " + to_string(-mycoefs[i] * value);
                    }
                }
                problem->postKnapsackConstraint(myvars, params, false, true, false);
                break;
            case OrderType::EQ:
                params = to_string(rightcoef);
                for (unsigned int i=0; i < myvars.size(); i++) {
                    int domsize = problem->getDomainInitSize(myvars[i]);
                    params += " " + to_string(domsize);
                    for (int idval=0; idval < domsize; idval++) {
                        int value = problem->toValue(myvars[i], idval);
                        params += " " + to_string(value) + " " + to_string(mycoefs[i] * value);
                    }
                }
                problem->postKnapsackConstraint(myvars, params, false, true, false);
                params = to_string(-rightcoef);
                for (unsigned int i=0; i < myvars.size(); i++) {
                    int domsize = problem->getDomainInitSize(myvars[i]);
                    params += " " + to_string(domsize);
                    for (int idval=0; idval < domsize; idval++) {
                        int value = problem->toValue(myvars[i], idval);
                        params += " " + to_string(value) + " " + to_string(-mycoefs[i] * value);
                    }
                }
                problem->postKnapsackConstraint(myvars, params, false, true, false);
                break;
            default:
                cerr << "Sorry operator " << cond.op << " not implemented in sum constraint!" << endl;
                throw WrongFileFormat();
            }
            break;
        case OperandType::INTERVAL:
            assert(cond.op == OrderType::IN);
            // sum >= cond.min
            params = to_string(cond.min);
            for (unsigned int i=0; i < myvars.size(); i++) {
                int domsize = problem->getDomainInitSize(myvars[i]);
                params += " " + to_string(domsize);
                for (int idval=0; idval < domsize; idval++) {
                    int value = problem->toValue(myvars[i], idval);
                    params += " " + to_string(value) + " " + to_string(mycoefs[i] * value);
                }
            }
            problem->postKnapsackConstraint(myvars, params, false, true, false);
            // sum <= cond.max
            params = to_string(-cond.max);
            for (unsigned int i=0; i < myvars.size(); i++) {
                int domsize = problem->getDomainInitSize(myvars[i]);
                params += " " + to_string(domsize);
                for (int idval=0; idval < domsize; idval++) {
                    int value = problem->toValue(myvars[i], idval);
                    params += " " + to_string(value) + " " + to_string(-mycoefs[i] * value);
                }
            }
            problem->postKnapsackConstraint(myvars, params, false, true, false);
            break;
        default:
            cerr << "Sorry operandType " << cond.operandType << " not implemented in sum constraint!" << endl;
            throw WrongFileFormat();
        }
    }

    void buildConstraintSum(string id, vector<XVariable *> &list, vector<int> &coefs, XCondition &cond) override {
        vector<int> vars;
        toMyVariables(list,vars);
        buildConstraintSum(vars, coefs, cond);
    }

    void buildConstraintSum(string id, vector<XVariable *> &list, XCondition &cond) override {
        vector<int> coefs(list.size(), 1);
        buildConstraintSum(id, list, coefs, cond);
    }

    void buildConstraintSum(string id, vector<Tree *> &trees, vector<int> &coefs, XCondition &cond) override {
        vector<int> vars;
        for (unsigned int i = 0; i < trees.size(); i++) {
            vars.push_back(buildConstraintIntension(trees[i]));
        }
        assert(vars.size() == trees.size());
        buildConstraintSum(vars, coefs, cond);
    }

    void buildConstraintSum(string id, vector<Tree *> &trees, XCondition &cond) override {
        vector<int> coefs(trees.size(), 1);
        buildConstraintSum(id, trees, coefs, cond);
    }

    void buildConstraintCount(vector<int> &vars, vector<int> &values, XCondition &cond) {
        vector<int> myvars = vars;
        vector<int> myvalues = values;
        int rightcoef = 0;
        string params = "";
        string params2 = "";
        int domsize;
        int nbval;
        switch (cond.operandType) {
            case OperandType::VARIABLE:
                myvars.push_back(getMyVar(cond.var));
                switch (cond.op) {
                    case OrderType::LE:
                        params = to_string(0);
                        params2="";
                        for (unsigned int i=0; i < myvars.size()-1; i++) {
                            domsize = problem->getDomainInitSize(myvars[i]);
                            //params += " " + to_string(domsize);
                            nbval=0;
                            params2="";
                            for (int idval=0; idval < domsize; idval++) {
                                int value = problem->toValue(myvars[i], idval);
                                auto it = find(values.begin(), values.end(), value);
                                if (it != values.end()) {
                                    nbval+=1;
                                    params2 += " " + to_string(value) + " " + to_string(-1);
                                }
                            }
                            params+=" "+to_string(nbval)+params2;
                        }
                        domsize=problem->getDomainInitSize(myvars.back());
                        params += " " + to_string(domsize);
                        for (int idval=0; idval < domsize; idval++) {
                            int value = problem->toValue(myvars.back(), idval);
                            params += " " + to_string(value) + " " + to_string(value);
                        }
                        problem->postKnapsackConstraint(myvars, params, false, true, false);
                        break;
                    case OrderType::LT:
                        params = to_string(1);
                        params2="";
                        for (unsigned int i=0; i < myvars.size()-1; i++) {
                            domsize = problem->getDomainInitSize(myvars[i]);
                            //params += " " + to_string(domsize);
                            nbval=0;
                            params2="";
                            for (int idval=0; idval < domsize; idval++) {
                                int value = problem->toValue(myvars[i], idval);
                                auto it = find(values.begin(), values.end(), value);
                                if (it != values.end()) {
                                    nbval+=1;
                                    params2 += " " + to_string(value) + " " + to_string(-1);
                                }
                            }
                            params+=" "+to_string(nbval)+params2;
                        }
                        domsize=problem->getDomainInitSize(myvars.back());
                        params += " " + to_string(domsize);
                        for (int idval=0; idval < domsize; idval++) {
                            int value = problem->toValue(myvars.back(), idval);
                            params += " " + to_string(value) + " " + to_string(value);
                        }
                        problem->postKnapsackConstraint(myvars, params, false, true, false);
                        break;
                    case OrderType::GE:
                        params = to_string(0);
                        params2="";
                        for (unsigned int i=0; i < myvars.size()-1; i++) {
                            domsize = problem->getDomainInitSize(myvars[i]);
                            //params += " " + to_string(domsize);
                            nbval=0;
                            params2="";
                            for (int idval=0; idval < domsize; idval++) {
                                int value = problem->toValue(myvars[i], idval);
                                auto it = find(values.begin(), values.end(), value);
                                if (it != values.end()) {
                                    nbval+=1;
                                    params2 += " " + to_string(value) + " " + to_string(1);
                                }
                            }
                            params+=" "+to_string(nbval)+params2;
                        }
                        domsize=problem->getDomainInitSize(myvars.back());
                        params += " " + to_string(domsize);
                        for (int idval=0; idval < domsize; idval++) {
                            int value = problem->toValue(myvars.back(), idval);
                            params += " " + to_string(value) + " " + to_string(-value);
                        }
                        problem->postKnapsackConstraint(myvars, params, false, true, false);
                        break;
                    case OrderType::GT:
                        params = to_string(1);
                        params2="";
                        for (unsigned int i=0; i < myvars.size()-1; i++) {
                            domsize = problem->getDomainInitSize(myvars[i]);
                            //params += " " + to_string(domsize);
                            nbval=0;
                            params2="";
                            for (int idval=0; idval < domsize; idval++) {
                                int value = problem->toValue(myvars[i], idval);
                                auto it = find(values.begin(), values.end(), value);
                                if (it != values.end()) {
                                    nbval+=1;
                                    params2 += " " + to_string(value) + " " + to_string(1);
                                }
                            }
                            params+=" "+to_string(nbval)+params2;
                        }
                        domsize=problem->getDomainInitSize(myvars.back());
                        params += " " + to_string(domsize);
                        for (int idval=0; idval < domsize; idval++) {
                            int value = problem->toValue(myvars.back(), idval);
                            params += " " + to_string(value) + " " + to_string(-value);
                        }
                        problem->postKnapsackConstraint(myvars, params, false, true, false);
                        break;
                    case OrderType::NE:
                        params = to_string(1);
                        params2="";
                        for (unsigned int i=0; i < myvars.size()-1; i++) {
                            domsize = problem->getDomainInitSize(myvars[i]);
                            //params += " " + to_string(domsize);
                            nbval=0;
                            params2="";
                            for (int idval=0; idval < domsize; idval++) {
                                int value = problem->toValue(myvars[i], idval);
                                auto it = find(values.begin(), values.end(), value);
                                if (it != values.end()) {
                                    nbval+=1;
                                    params2 += " " + to_string(value) + " " + to_string(1);
                                }
                            }
                            params+=" "+to_string(nbval)+params2;
                        }
                        domsize=problem->getDomainInitSize(myvars.back());
                        params += " " + to_string(domsize);
                        for (int idval=0; idval < domsize; idval++) {
                            int value = problem->toValue(myvars.back(), idval);
                            params += " " + to_string(value) + " " + to_string(-value);
                        }
                        problem->postKnapsackConstraint(myvars, params, false, true, false);
                        params = to_string(1);
                        params2="";
                        for (unsigned int i=0; i < myvars.size()-1; i++) {
                            domsize = problem->getDomainInitSize(myvars[i]);
                            //params += " " + to_string(domsize);
                            nbval=0;
                            params2="";
                            for (int idval=0; idval < domsize; idval++) {
                                int value = problem->toValue(myvars[i], idval);
                                auto it = find(values.begin(), values.end(), value);
                                if (it != values.end()) {
                                    nbval+=1;
                                    params2 += " " + to_string(value) + " " + to_string(-1);
                                }
                            }
                            params+=" "+to_string(nbval)+params2;
                        }
                        domsize=problem->getDomainInitSize(myvars.back());
                        params += " " + to_string(domsize);
                        for (int idval=0; idval < domsize; idval++) {
                            int value = problem->toValue(myvars.back(), idval);
                            params += " " + to_string(value) + " " + to_string(value);
                        }
                        problem->postKnapsackConstraint(myvars, params, false, true, false);
                        break;
                    case OrderType::EQ:
                        params = to_string(0);
                        params2="";
                        for (unsigned int i=0; i < myvars.size()-1; i++) {
                            domsize = problem->getDomainInitSize(myvars[i]);
                            //params += " " + to_string(domsize);
                            nbval=0;
                            params2="";
                            params2="";
                            for (int idval=0; idval < domsize; idval++) {
                                int value = problem->toValue(myvars[i], idval);
                                auto it = find(values.begin(), values.end(), value);
                                if (it != values.end()) {
                                    nbval+=1;
                                    params2 += " " + to_string(value) + " " + to_string(1);
                                }
                            }
                            params+=" "+to_string(nbval)+params2;
                        }
                        domsize=problem->getDomainInitSize(myvars.back());
                        params += " " + to_string(domsize);
                        for (int idval=0; idval < domsize; idval++) {
                            int value = problem->toValue(myvars.back(), idval);
                            params += " " + to_string(value) + " " + to_string(-value);
                        }
                        problem->postKnapsackConstraint(myvars, params, false, true, false);
                        params = to_string(0);
                        params2="";
                        for (unsigned int i=0; i < myvars.size()-1; i++) {
                            domsize = problem->getDomainInitSize(myvars[i]);
                            //params += " " + to_string(domsize);
                            nbval=0;
                            params2="";
                            for (int idval=0; idval < domsize; idval++) {
                                int value = problem->toValue(myvars[i], idval);
                                auto it = find(values.begin(), values.end(), value);
                                if (it != values.end()) {
                                    nbval+=1;
                                    params2 += " " + to_string(value) + " " + to_string(-1);
                                }
                            }
                            params+=" "+to_string(nbval)+params2;
                        }
                        domsize=problem->getDomainInitSize(myvars.back());
                        params += " " + to_string(domsize);
                        for (int idval=0; idval < domsize; idval++) {
                            int value = problem->toValue(myvars.back(), idval);
                            params += " " + to_string(value) + " " + to_string(value);
                        }
                        problem->postKnapsackConstraint(myvars, params, false, true, false);
                        break;
                    default:
                        cerr << "Sorry operator " << cond.op << " not implemented in sum constraint!" << endl;
                        throw WrongFileFormat();
                }
                break;
            case OperandType::INTEGER:
                rightcoef += cond.val;
                switch (cond.op) {
                    case OrderType::LE:
                        params = to_string(-rightcoef);
                        params2="";
                        for (unsigned int i=0; i < myvars.size(); i++) {
                            domsize = problem->getDomainInitSize(myvars[i]);
                            //params += " " + to_string(domsize);
                            nbval=0;
                            params2="";
                            for (int idval=0; idval < domsize; idval++) {
                                int value = problem->toValue(myvars[i], idval);
                                auto it = find(values.begin(), values.end(), value);
                                if (it != values.end()) {
                                    nbval+=1;
                                    params2 += " " + to_string(value) + " " + to_string(-1);
                                }
                            }
                            params+=" "+to_string(nbval)+params2;
                        }
                        problem->postKnapsackConstraint(myvars, params, false, true, false);
                        break;
                    case OrderType::LT:
                        params = to_string(-rightcoef+1);
                        params2="";
                        for (unsigned int i=0; i < myvars.size(); i++) {
                            domsize = problem->getDomainInitSize(myvars[i]);
                            //params += " " + to_string(domsize);
                            nbval=0;
                            params2="";
                            for (int idval=0; idval < domsize; idval++) {
                                int value = problem->toValue(myvars[i], idval);
                                auto it = find(values.begin(), values.end(), value);
                                if (it != values.end()) {
                                    nbval+=1;
                                    params2 += " " + to_string(value) + " " + to_string(-1);
                                }
                            }
                            params+=" "+to_string(nbval)+params2;
                        }
                        problem->postKnapsackConstraint(myvars, params, false, true, false);
                        break;
                    case OrderType::GE:
                        params = to_string(rightcoef);
                        params2="";
                        for (unsigned int i=0; i < myvars.size(); i++) {
                            domsize = problem->getDomainInitSize(myvars[i]);
                            //params += " " + to_string(domsize);
                            nbval=0;
                            params2="";
                            for (int idval=0; idval < domsize; idval++) {
                                int value = problem->toValue(myvars[i], idval);
                                auto it = find(values.begin(), values.end(), value);
                                if (it != values.end()) {
                                    nbval+=1;
                                    params2 += " " + to_string(value) + " " + to_string(1);
                                }
                            }
                            params+=" "+to_string(nbval)+params2;
                        }
                        problem->postKnapsackConstraint(myvars, params, false, true, false);
                        break;
                    case OrderType::GT:
                        params = to_string(rightcoef+1);
                        params2="";
                        for (unsigned int i=0; i < myvars.size(); i++) {
                            domsize = problem->getDomainInitSize(myvars[i]);
                            //params += " " + to_string(domsize);
                            nbval=0;
                            params2="";
                            for (int idval=0; idval < domsize; idval++) {
                                int value = problem->toValue(myvars[i], idval);
                                auto it = find(values.begin(), values.end(), value);
                                if (it != values.end()) {
                                    nbval+=1;
                                    params2 += " " + to_string(value) + " " + to_string(1);
                                }
                            }
                            params+=" "+to_string(nbval)+params2;
                        }
                        problem->postKnapsackConstraint(myvars, params, false, true, false);
                        break;
                    case OrderType::NE:
                        params = to_string(rightcoef+1);
                        params2="";
                        for (unsigned int i=0; i < myvars.size(); i++) {
                            domsize = problem->getDomainInitSize(myvars[i]);
                            //params += " " + to_string(domsize);
                            nbval=0;
                            params2="";
                            for (int idval=0; idval < domsize; idval++) {
                                int value = problem->toValue(myvars[i], idval);
                                auto it = find(values.begin(), values.end(), value);
                                if (it != values.end()) {
                                    nbval+=1;
                                    params2 += " " + to_string(value) + " " + to_string(1);
                                }
                            }
                            params+=" "+to_string(nbval)+params2;
                        }
                        problem->postKnapsackConstraint(myvars, params, false, true, false);
                        params = to_string(-rightcoef+1);
                        params2="";
                        for (unsigned int i=0; i < myvars.size(); i++) {
                            domsize = problem->getDomainInitSize(myvars[i]);
                            //params += " " + to_string(domsize);
                            nbval=0;
                            params2="";
                            for (int idval=0; idval < domsize; idval++) {
                                int value = problem->toValue(myvars[i], idval);
                                auto it = find(values.begin(), values.end(), value);
                                if (it != values.end()) {
                                    nbval+=1;
                                    params2 += " " + to_string(value) + " " + to_string(-1);
                                }
                            }
                            params+=" "+to_string(nbval)+params2;
                        }
                        problem->postKnapsackConstraint(myvars, params, false, true, false);
                        break;
                    case OrderType::EQ:
                        params = to_string(rightcoef);
                        params2="";
                        for (unsigned int i=0; i < myvars.size(); i++) {
                            domsize = problem->getDomainInitSize(myvars[i]);
                            //params += " " + to_string(domsize);
                            nbval=0;
                            params2="";
                            for (int idval=0; idval < domsize; idval++) {
                                int value = problem->toValue(myvars[i], idval);
                                auto it = find(values.begin(), values.end(), value);
                                if (it != values.end()) {
                                    nbval+=1;
                                    params2 += " " + to_string(value) + " " + to_string(1);
                                }
                            }
                            params+=" "+to_string(nbval)+params2;
                        }
                        problem->postKnapsackConstraint(myvars, params, false, true, false);
                        params = to_string(-rightcoef);
                        params2="";
                        for (unsigned int i=0; i < myvars.size(); i++) {
                            domsize = problem->getDomainInitSize(myvars[i]);
                            //params += " " + to_string(domsize);
                            nbval=0;
                            params2="";
                            for (int idval=0; idval < domsize; idval++) {
                                int value = problem->toValue(myvars[i], idval);
                                auto it = find(values.begin(), values.end(), value);
                                if (it != values.end()) {
                                    nbval+=1;
                                    params2 += " " + to_string(value) + " " + to_string(-1);
                                }
                            }
                            params+=" "+to_string(nbval)+params2;
                        }
                        problem->postKnapsackConstraint(myvars, params, false, true, false);
                        break;
                    default:
                        cerr << "Sorry operator " << cond.op << " not implemented in sum constraint!" << endl;
                        throw WrongFileFormat();
                }
                break;
            case OperandType::INTERVAL:
                assert(cond.op == OrderType::IN);
                // sum >= cond.min
                params = to_string(cond.min);
                params2="";
                for (unsigned int i=0; i < myvars.size(); i++) {
                    domsize = problem->getDomainInitSize(myvars[i]);
                    //params += " " + to_string(domsize);
                    nbval=0;
                    params2="";
                    for (int idval=0; idval < domsize; idval++) {
                        int value = problem->toValue(myvars[i], idval);
                        auto it = find(values.begin(), values.end(), value);
                        if (it != values.end()) {
                            nbval+=1;
                            params2 += " " + to_string(value) + " " + to_string(1);
                        }
                    }
                    params+=" "+to_string(nbval)+params2;
                }
                problem->postKnapsackConstraint(myvars, params, false, true, false);
                // sum <= cond.max
                params = to_string(-cond.max);
                params2="";
                for (unsigned int i=0; i < myvars.size(); i++) {
                    domsize = problem->getDomainInitSize(myvars[i]);
                    //params += " " + to_string(domsize);
                    nbval=0;
                    params2="";
                    for (int idval=0; idval < domsize; idval++) {
                        int value = problem->toValue(myvars[i], idval);
                        auto it = find(values.begin(), values.end(), value);
                        if (it != values.end()) {
                            nbval+=1;
                            params2 += " " + to_string(value) + " " + to_string(-1);
                        }
                    }
                    params+=" "+to_string(nbval)+params2;
                }
                problem->postKnapsackConstraint(myvars, params, false, true, false);
                break;
            default:
                cerr << "Sorry operandType " << cond.operandType << " not implemented in sum constraint!" << endl;
                throw WrongFileFormat();
        }
    }

    void buildConstraintCount(string id, vector<XVariable *> &list, vector<int> &values, XCondition &cond) override {
        vector<int> vars;
        toMyVariables(list,vars);
        buildConstraintCount(vars, values, cond);
    }

    void buildConstraintCount(string id, vector<Tree*> &trees, vector<int> &values, XCondition &cond) override {
        vector<int> vars;
        for (unsigned int i = 0; i < trees.size(); i++) {
            vars.push_back(buildConstraintIntension(trees[i]));
        }
        assert(vars.size() == trees.size());
        buildConstraintCount(vars, values, cond);
    }

    void buildConstraintNValues(vector<int> &vars, vector<int> &except, XCondition &cond) {
        vector<int> Bvar;
        int n=vars.size();
        vector<int> Diffvalue;
        int currentval,domsize;
        vector <Cost> costs;
        bool emptyExcept=false;
        string params;
        int rightcoef;
        if(except.empty())
            emptyExcept= true;
        for (int i = 0; i < n; i++) {
            for (unsigned int a = 0; a < problem->getDomainInitSize(vars[i]); a++) {
                currentval = problem->toValue(vars[i], a);
                if (!emptyExcept && currentval != except[0]){
                    auto it = find(Diffvalue.begin(), Diffvalue.end(), currentval);
                    if (it == Diffvalue.end()) {
                        string extravarname = IMPLICIT_VAR_TAG + "Nval" + to_string(problem->numberOfVariables());
                        int extravar = problem->makeEnumeratedVariable(extravarname, 0, 1);
                        mapping[extravarname] = extravar;
                        Bvar.push_back(extravar);
                        Diffvalue.push_back(currentval);
                        costs.clear();
                        for (int b = 0; b < (int)problem->getDomainInitSize(vars[i]); b++) {
                            if (b == currentval) {
                                costs.push_back(MAX_COST);
                            } else {
                                costs.push_back(MIN_COST);
                            }
                        }
                        for (int b = 0; b < (int)problem->getDomainInitSize(vars[i]); b++) {
                            if (b == currentval) {
                                costs.push_back(MIN_COST);
                            } else {
                                costs.push_back(MAX_COST);
                            }
                        }
                        problem->postBinaryConstraint(Bvar.back(), vars[i], costs);
                    }else{
                        costs.clear();
                        for (int b = 0; b < (int)problem->getDomainInitSize(vars[i]); b++) {
                            if (b == currentval) {
                                costs.push_back(MAX_COST);
                            } else {
                                costs.push_back(MIN_COST);
                            }
                        }
                        for (int b = 0; b < (int)problem->getDomainInitSize(vars[i]); b++) {
                            if (b == currentval) {
                                costs.push_back(MIN_COST);
                            } else {
                                costs.push_back(MAX_COST);
                            }
                        }
                        problem->postBinaryConstraint(Bvar[distance(Diffvalue.begin(), it)], vars[i], costs);
                    }
                }
            }
        }
        switch (cond.op) {
            case OrderType::EQ:
                switch (cond.operandType) {
                    case OperandType::VARIABLE:
                        Bvar.push_back(getMyVar(cond.var));
                        params = to_string(0);
                        for (unsigned int i=0; i < Bvar.size()-1; i++) {
                            params+=" 2 0 0 1 1";
                        }
                        domsize=problem->getDomainInitSize(Bvar.back());
                        params += " " + to_string(domsize);
                        for (int idval=0; idval < domsize; idval++) {
                            int value = problem->toValue(Bvar.back(), idval);
                            params += " " + to_string(value) + " " + to_string(-value);
                        }
                        problem->postKnapsackConstraint(Bvar, params, false, true, false);
                        params = to_string(0);
                        for (unsigned int i=0; i < Bvar.size()-1; i++) {
                            params+=" 2 0 0 1 -1";
                        }
                        domsize=problem->getDomainInitSize(Bvar.back());
                        params += " " + to_string(domsize);
                        for (int idval=0; idval < domsize; idval++) {
                            int value = problem->toValue(Bvar.back(), idval);
                            params += " " + to_string(value) + " " + to_string(value);
                        }
                        problem->postKnapsackConstraint(Bvar, params, false, true, false);
                        break;
                    case OperandType::INTEGER:
                        rightcoef = cond.val;
                        params = to_string(rightcoef);
                        for (unsigned int i=0; i < Bvar.size(); i++) {
                            params+=" 2 0 0 1 1";
                        }
                        problem->postKnapsackConstraint(Bvar, params, false, true, false);
                        params = to_string(-rightcoef);
                        for (unsigned int i=0; i < Bvar.size(); i++) {
                            params+=" 2 0 0 1 -1";
                        }
                        problem->postKnapsackConstraint(Bvar, params, false, true, false);
                        break;
                    default:
                        cerr << "Sorry operandType " << cond.operandType << " not implemented in Nvalues constraint!" << endl;
                        throw WrongFileFormat();
                }
                break;
            case OrderType::GT:
                params = to_string(1);
                for (unsigned int i=0; i < Bvar.size(); i++) {
                    params+=" 2 0 0 1 1";
                }
                problem->postKnapsackConstraint(Bvar, params, false, true, false);
                break;
            default:
                cerr << "Sorry operandType " << cond.operandType << " not implemented in Nvalues constraint!" << endl;
                throw WrongFileFormat();
        }
    }

    void buildConstraintNValues(string id, vector<XVariable *> &list, vector<int> &except, XCondition &cond) override {
        vector<int> vars;
        toMyVariables(list,vars);
        buildConstraintNValues(vars, except, cond);
    }

    void buildConstraintNValues(string id, vector<XVariable *> &list, XCondition &cond) override{
        vector<int> vars;
        toMyVariables(list,vars);
        vector<int> except;
        buildConstraintNValues(vars, except, cond);
    }

    void buildConstraintNValues(string id, vector<Tree *> &trees, XCondition &cond) override {
        vector<int> vars,except;
        for (unsigned int i = 0; i < trees.size(); i++) {
            vars.push_back(buildConstraintIntension(trees[i]));
        }
        assert(vars.size() == trees.size());
        buildConstraintNValues(vars,except, cond);
    }

    void buildConstraintCardinality(string id, vector<XVariable *> &list, vector<int> values, vector<int> &occurs, bool closed) {
        vector<int> vars;
        toMyVariables(list,vars);
        string params="";
        string params2="";
        int domsize;
        int nbval;
        for (int k = 0; k < (int)values.size(); ++k) {
            params=to_string(occurs[k]);
            for (int i = 0; i < (int)vars.size(); ++i) {
                domsize = problem->getDomainInitSize(vars[i]);
                nbval=0;
                params2="";
                //params +=" "+to_string(domsize);
                for (int idval=0; idval < domsize; idval++) {
                    int value = problem->toValue(vars[i], idval);
                    if(value==values[k]){
                        params2 += " " + to_string(value) + " 1";
                        nbval++;
                    }
                    else {
                        if(closed)
                        {
                            auto it =find(values.begin(),values.end(),value);
                            if(it==values.end()) {
                                params2 += " " + to_string(value) + " " + to_string(-vars.size());
                                nbval++;
                            }
                        }
                    }
                }
                params+=" "+to_string(nbval)+params2;
            }
            problem->postKnapsackConstraint(vars, params, false, true, false);
            params=to_string(-occurs[k])+ " ";
            for (int i = 0; i < (int)vars.size(); ++i) {
                domsize = problem->getDomainInitSize(vars[i]);
                nbval=0;
                params2="";
                //params +=" "+to_string(domsize);
                for (int idval=0; idval < domsize; idval++) {
                    int value = problem->toValue(vars[i], idval);
                    if(value==values[k]){
                        params2 += " " + to_string(value) + " -1";
                        nbval++;
                    }
                    else {
                        if(closed)
                        {
                            auto it =find(values.begin(),values.end(),value);
                            if(it==values.end()) {
                                params2 += " " + to_string(value) + " " + to_string(-vars.size());
                                nbval++;
                            }
                        }
                    }
                }
                params+=" "+to_string(nbval)+params2;
            }
            problem->postKnapsackConstraint(vars, params, false, true, false);
        }
    }

    void buildConstraintCardinality(string id, vector<XVariable *> &list, vector<int> values, vector<XVariable *> &occurs, bool closed){
        vector<int> vars;
        toMyVariables(list,vars);
        string params;
        string params2;
        int domsize,nbval;
        for (int k = 0; k < (int)values.size(); ++k) {
            params="0 ";
            vars.push_back(getMyVar(occurs[k]));
            for (int i = 0; i < (int)vars.size()-1; ++i) {
                domsize = problem->getDomainInitSize(vars[i]);
                nbval=0;
                params2="";
                //params +=" "+to_string(domsize);
                for (int idval=0; idval < domsize; idval++) {
                    int value = problem->toValue(vars[i], idval);
                    if(value==values[k]){
                        params2 += " " + to_string(value) + " 1";
                        nbval++;
                    }
                    else {
                        if(closed)
                        {
                            auto it =find(values.begin(),values.end(),value);
                            if(it==values.end()) {
                                params2 += " " + to_string(value) + " " + to_string(-vars.size());
                                nbval++;
                            }
                        }
                    }
                }
                params+=" "+to_string(nbval)+params2;
            }
            domsize= problem->getDomainInitSize(vars.back());
            for (int i = 0; i < domsize; ++i) {
                for (int idval=0; idval < domsize; idval++) {
                    int value = problem->toValue(vars.back(), idval);
                        params += " " + to_string(value) + " " + to_string(-value);
                }
            }
            problem->postKnapsackConstraint(vars, params, false, true, false);
            params="0 ";
            for (int i = 0; i < (int)vars.size()-1; ++i) {
                domsize = problem->getDomainInitSize(vars[i]);
                nbval=0;
                params2="";
                //params +=" "+to_string(domsize);
                for (int idval=0; idval < domsize; idval++) {
                    int value = problem->toValue(vars[i], idval);
                    if(value==values[k]){
                        params2 += " " + to_string(value) + " -1";
                        nbval++;
                    }
                    else {
                        if(closed)
                        {
                            auto it =find(values.begin(),values.end(),value);
                            if(it==values.end()) {
                                params2 += " " + to_string(value) + " " + to_string(-vars.size());
                                nbval++;
                            }
                        }
                    }
                }
                params+=" "+to_string(nbval)+params2;
            }
            domsize= problem->getDomainInitSize(vars.back());
            for (int i = 0; i < domsize; ++i) {
                for (int idval=0; idval < domsize; idval++) {
                    int value = problem->toValue(vars.back(), idval);
                    params += " " + to_string(value) + " " + to_string(value);
                }
            }
            problem->postKnapsackConstraint(vars, params, false, true, false);
        }
    }

    void buildConstraintCardinality(string id, vector<XVariable *> &list, vector<int> values, vector<XInterval> &occurs, bool closed){
        vector<int> vars;
        toMyVariables(list,vars);
        string params="";
        string params2="";
        int domsize,nbval;
        for (int k = 0; k < (int)values.size(); ++k) {
            params=to_string(occurs[k].min)+ " ";
            for (int i = 0; i < (int)vars.size(); ++i) {
                domsize = problem->getDomainInitSize(vars[i]);
                nbval=0;
                params2="";
                //params +=" "+to_string(domsize);
                for (int idval=0; idval < domsize; idval++) {
                    int value = problem->toValue(vars[i], idval);
                    if(value==values[k]){
                        params2 += " " + to_string(value) + " 1";
                        nbval++;
                    }
                    else {
                        if(closed)
                        {
                            auto it =find(values.begin(),values.end(),value);
                            if(it==values.end()) {
                                params2 += " " + to_string(value) + " " + to_string(-vars.size());
                                nbval++;
                            }
                        }
                    }
                }
                params+=" "+to_string(nbval)+params2;
            }
            problem->postKnapsackConstraint(vars, params, false, true, false);
            params=to_string(-occurs[k].max)+ " ";
            for (int i = 0; i < (int)vars.size(); ++i) {
                domsize = problem->getDomainInitSize(vars[i]);
                nbval=0;
                params2="";
                //params +=" "+to_string(domsize);
                for (int idval=0; idval < domsize; idval++) {
                    int value = problem->toValue(vars[i], idval);
                    if(value==values[k]){
                        params2 += " " + to_string(value) + " -1";
                        nbval++;
                    }
                    else {
                        if(closed)
                        {
                            auto it =find(values.begin(),values.end(),value);
                            if(it==values.end()) {
                                params2 += " " + to_string(value) + " " + to_string(-vars.size());
                                nbval++;
                            }
                        }
                    }
                }
                params+=" "+to_string(nbval)+params2;
            }
            problem->postKnapsackConstraint(vars, params, false, true, false);
        }
    }

    void buildConstraintMinMax(bool max, vector<int> &vars, int varargmax, XCondition &cond) {
        assert(vars.size() > 0);
        int maxinf = INT_MAX;
        int maxsup = -INT_MAX;
        for (unsigned int i=0; i<vars.size(); i++) {
            if (problem->getInf(vars[i]) < maxinf) {
                maxinf = problem->getInf(vars[i]);
            }
            if (problem->getSup(vars[i]) > maxsup) {
                maxsup = problem->getSup(vars[i]);
            }
        }
        string varmaxname = IMPLICIT_VAR_TAG + ((max)?"max":"min") + to_string(problem->numberOfVariables());
        int varmax = problem->makeEnumeratedVariable(varmaxname, maxinf, maxsup);
        mapping[varmaxname] = varmax;
        buildConstraintElement(OrderType::EQ, vars, varargmax, varmax);
        for (unsigned int i=0; i<vars.size(); i++) {
            buildConstraintPrimitive((max)?OrderType::LE:OrderType::GE, vars[i], 0, varmax);
        }
        switch (cond.operandType) {
        case OperandType::VARIABLE:
            buildConstraintPrimitive(cond.op, varmax, 0, getMyVar(cond.var));
            break;
        case OperandType::INTEGER:
            buildConstraintPrimitive(cond.op, varmax, cond.val);
            break;
        case OperandType::INTERVAL:
            assert(cond.op == OrderType::IN);
            buildConstraintPrimitive(varmax, true, cond.min, cond.max);
            break;
        default:
            cerr << "Sorry operandType " << cond.operandType << " not implemented in min/max constraint!" << endl;
            throw WrongFileFormat();
        }
    }

    void buildConstraintMaximum(string id, vector<XVariable *> &list, XCondition &cond) override {
        vector<int> vars;
        toMyVariables(list,vars);
        string varargmaxname = IMPLICIT_VAR_TAG + "argmax" + to_string(problem->numberOfVariables());
        int varargmax = problem->makeEnumeratedVariable(varargmaxname, 0, vars.size()-1);
        mapping[varargmaxname] = varargmax;
        buildConstraintMinMax(true, vars, varargmax, cond);
    }

    void buildConstraintMinimum(string id, vector<XVariable *> &list, XCondition &cond) override {
        vector<int> vars;
        toMyVariables(list,vars);
        string varargmaxname = IMPLICIT_VAR_TAG + "argmin" + to_string(problem->numberOfVariables());
        int varargmax = problem->makeEnumeratedVariable(varargmaxname, 0, vars.size()-1);
        mapping[varargmaxname] = varargmax;
        buildConstraintMinMax(false, vars, varargmax, cond);
    }

    void buildConstraintMaximum(string id, vector<Tree*> &trees, XCondition &cond) override {
        vector<int> vars;
        for (unsigned int i = 0; i < trees.size(); i++) {
            vars.push_back(buildConstraintIntension(trees[i]));
        }
        assert(vars.size() == trees.size());
        string varargmaxname = IMPLICIT_VAR_TAG + "argmax" + to_string(problem->numberOfVariables());
        int varargmax = problem->makeEnumeratedVariable(varargmaxname, 0, vars.size()-1);
        mapping[varargmaxname] = varargmax;
        buildConstraintMinMax(true, vars, varargmax, cond);
    }

    void buildConstraintMinimum(string id, vector<Tree*> &trees, XCondition &cond) override {
        vector<int> vars;
        for (unsigned int i = 0; i < trees.size(); i++) {
            vars.push_back(buildConstraintIntension(trees[i]));
        }
        assert(vars.size() == trees.size());
        string varargmaxname = IMPLICIT_VAR_TAG + "argmin" + to_string(problem->numberOfVariables());
        int varargmax = problem->makeEnumeratedVariable(varargmaxname, 0, vars.size()-1);
        mapping[varargmaxname] = varargmax;
        buildConstraintMinMax(false, vars, varargmax, cond);
    }

    void buildConstraintMaximum(string id, vector<XVariable *> &list, XVariable *index, int startIndex, RankType rank, XCondition &cond) override {
        assert(startIndex == 0);
        assert(rank == RankType::ANY);
        vector<int> vars;
        toMyVariables(list,vars);
        buildConstraintMinMax(true, vars, getMyVar(index), cond);
    }

    void buildConstraintMinimum(string id, vector<XVariable *> &list, XVariable *index, int startIndex, RankType rank, XCondition &cond) override {
        assert(startIndex == 0);
        assert(rank == RankType::ANY);
        vector<int> vars;
        toMyVariables(list,vars);
        buildConstraintMinMax(false, vars, getMyVar(index), cond);
    }

    void buildConstraintMinMaxArg(bool max, vector<int> &vars, int varargmax, XCondition &cond) {
        assert(vars.size() > 0);
        int maxinf = INT_MAX;
        int maxsup = -INT_MAX;
        for (unsigned int i=0; i<vars.size(); i++) {
            if (problem->getInf(vars[i]) < maxinf) {
                maxinf = problem->getInf(vars[i]);
            }
            if (problem->getSup(vars[i]) > maxsup) {
                maxsup = problem->getSup(vars[i]);
            }
        }
        string varmaxname = IMPLICIT_VAR_TAG + ((max)?"max":"min") + to_string(problem->numberOfVariables());
        int varmax = problem->makeEnumeratedVariable(varmaxname, maxinf, maxsup);
        mapping[varmaxname] = varmax;
        buildConstraintElement(OrderType::EQ, vars, varargmax, varmax);
        for (unsigned int i=0; i<vars.size(); i++) {
            buildConstraintPrimitive((max)?OrderType::LE:OrderType::GE, vars[i], 0, varmax);
        }
        switch (cond.operandType) {
        case OperandType::VARIABLE:
            buildConstraintPrimitive(cond.op, varargmax, 0, getMyVar(cond.var));
            break;
        case OperandType::INTEGER:
            buildConstraintPrimitive(cond.op, varargmax, cond.val);
            break;
        case OperandType::INTERVAL:
            assert(cond.op == OrderType::IN);
            buildConstraintPrimitive(varargmax, true, cond.min, cond.max);
            break;
        default:
            cerr << "Sorry operandType " << cond.operandType << " not implemented in argmin/argmax constraint!" << endl;
            throw WrongFileFormat();
        }
    }

    void buildConstraintMaximumArg(string id, vector<XVariable*> &list, RankType rank, XCondition &cond) override {
        vector<int> vars;
        toMyVariables(list,vars);
        string varargmaxname = IMPLICIT_VAR_TAG + "argmax" + to_string(problem->numberOfVariables());
        int varargmax = problem->makeEnumeratedVariable(varargmaxname, 0, vars.size()-1);
        mapping[varargmaxname] = varargmax;
        buildConstraintMinMaxArg(true, vars, varargmax, cond);
    }

    void buildConstraintMinimumArg(string id, vector<XVariable*> &list, RankType rank, XCondition &cond) override {
        vector<int> vars;
        toMyVariables(list,vars);
        string varargmaxname = IMPLICIT_VAR_TAG + "argmin" + to_string(problem->numberOfVariables());
        int varargmax = problem->makeEnumeratedVariable(varargmaxname, 0, vars.size()-1);
        mapping[varargmaxname] = varargmax;
        buildConstraintMinMaxArg(false, vars, varargmax, cond);
    }

    void buildConstraintMaximumArg(string id, vector<Tree*> &trees, RankType rank, XCondition &cond) override {
        assert(rank == RankType::ANY);
        vector<int> vars;
        for (unsigned int i = 0; i < trees.size(); i++) {
            vars.push_back(buildConstraintIntension(trees[i]));
        }
        assert(vars.size() == trees.size());
        string varargmaxname = IMPLICIT_VAR_TAG + "argmax" + to_string(problem->numberOfVariables());
        int varargmax = problem->makeEnumeratedVariable(varargmaxname, 0, vars.size()-1);
        mapping[varargmaxname] = varargmax;
        buildConstraintMinMaxArg(true, vars, varargmax, cond);
    }

    void buildConstraintMinimumArg(string id, vector<Tree*> &trees, RankType rank, XCondition &cond) override {
        assert(rank == RankType::ANY);
        vector<int> vars;
        for (unsigned int i = 0; i < trees.size(); i++) {
            vars.push_back(buildConstraintIntension(trees[i]));
        }
        assert(vars.size() == trees.size());
        string varargmaxname = IMPLICIT_VAR_TAG + "argmin" + to_string(problem->numberOfVariables());
        int varargmax = problem->makeEnumeratedVariable(varargmaxname, 0, vars.size()-1);
        mapping[varargmaxname] = varargmax;
        buildConstraintMinMaxArg(false, vars, varargmax, cond);
    }

    void buildConstraintElement(OrderType op, vector<int> &vars, int varindex, int varvalue) {
        assert(problem->getInf(varindex) >= 0);
        assert(problem->getSup(varindex) < (Value)vars.size());
        if (varindex != varvalue) {
            for (unsigned int i=0; i<vars.size(); i++) {
                if (problem->canbe(varindex, (Value)i)) {
                    if (varindex != vars[i] && varvalue != vars[i]) {
                        vector<Cost> costs(problem->getDomainInitSize(varvalue) * problem->getDomainInitSize(varindex) * problem->getDomainInitSize(vars[i]), MIN_COST);
                        for (unsigned int a=0; a < problem->getDomainInitSize(varvalue); a++) {
                            for (unsigned int b=0; b < problem->getDomainInitSize(vars[i]); b++) {
                                switch (op) {
                                case OrderType::IN:
                                case OrderType::EQ:
                                    if (problem->toValue(vars[i], b) != problem->toValue(varvalue, a)) {
                                        costs[a * problem->getDomainInitSize(varindex) * problem->getDomainInitSize(vars[i]) + i * problem->getDomainInitSize(vars[i]) + b] = MAX_COST;
                                    }
                                    break;
                                case OrderType::NE:
                                    if (problem->toValue(vars[i], b) == problem->toValue(varvalue, a)) {
                                        costs[a * problem->getDomainInitSize(varindex) * problem->getDomainInitSize(vars[i]) + i * problem->getDomainInitSize(vars[i]) + b] = MAX_COST;
                                    }
                                    break;
                                case OrderType::LE:
                                    if (problem->toValue(vars[i], b) > problem->toValue(varvalue, a)) {
                                        costs[a * problem->getDomainInitSize(varindex) * problem->getDomainInitSize(vars[i]) + i * problem->getDomainInitSize(vars[i]) + b] = MAX_COST;
                                    }
                                    break;
                                case OrderType::LT:
                                    if (problem->toValue(vars[i], b) >= problem->toValue(varvalue, a)) {
                                        costs[a * problem->getDomainInitSize(varindex) * problem->getDomainInitSize(vars[i]) + i * problem->getDomainInitSize(vars[i]) + b] = MAX_COST;
                                    }
                                    break;
                                case OrderType::GE:
                                    if (problem->toValue(vars[i], b) < problem->toValue(varvalue, a)) {
                                        costs[a * problem->getDomainInitSize(varindex) * problem->getDomainInitSize(vars[i]) + i * problem->getDomainInitSize(vars[i]) + b] = MAX_COST;
                                    }
                                    break;
                                case OrderType::GT:
                                    if (problem->toValue(vars[i], b) <= problem->toValue(varvalue, a)) {
                                        costs[a * problem->getDomainInitSize(varindex) * problem->getDomainInitSize(vars[i]) + i * problem->getDomainInitSize(vars[i]) + b] = MAX_COST;
                                    }
                                    break;
                                default:
                                    cerr << "Sorry operator " << op << " not implemented in element constraint!" << endl;
                                    throw WrongFileFormat();
                                }
                            }
                        }
                        problem->postTernaryConstraint(varvalue, varindex, vars[i], costs);
                    } else if (varvalue == vars[i]) {
                        switch (op) {
                        case OrderType::LE:
                        case OrderType::GE:
                        case OrderType::IN:
                        case OrderType::EQ:
                            // element constraint is always satisfied if varindex=i and vars[i]=varvalue
                            break;
                        case OrderType::LT:
                        case OrderType::GT:
                        case OrderType::NE:
                            buildConstraintPrimitive(OrderType::NE, varindex, (Value)i);
                            break;
                        default:
                            cerr << "Sorry operator " << op << " not implemented in element constraint!" << endl;
                            throw WrongFileFormat();
                        }
                    } else if (varindex == vars[i]) {
                        vector<Cost> costs(problem->getDomainInitSize(varvalue) * problem->getDomainInitSize(varindex), MIN_COST);
                        for (unsigned int a=0; a < problem->getDomainInitSize(varvalue); a++) {
                            switch (op) {
                            case OrderType::IN:
                            case OrderType::EQ:
                                if ((Value)i != problem->toValue(varvalue, a)) {
                                    costs[a * problem->getDomainInitSize(varindex) + i] = MAX_COST;
                                }
                                break;
                            case OrderType::NE:
                                if ((Value)i == problem->toValue(varvalue, a)) {
                                    costs[a * problem->getDomainInitSize(varindex) + i] = MAX_COST;
                                }
                                break;
                            case OrderType::LE:
                                if ((Value)i > problem->toValue(varvalue, a)) {
                                    costs[a * problem->getDomainInitSize(varindex) + i] = MAX_COST;
                                }
                                break;
                            case OrderType::LT:
                                if ((Value)i >= problem->toValue(varvalue, a)) {
                                    costs[a * problem->getDomainInitSize(varindex) + i] = MAX_COST;
                                }
                                break;
                            case OrderType::GE:
                                if ((Value)i < problem->toValue(varvalue, a)) {
                                    costs[a * problem->getDomainInitSize(varindex) + i] = MAX_COST;
                                }
                                break;
                            case OrderType::GT:
                                if ((Value)i <= problem->toValue(varvalue, a)) {
                                    costs[a * problem->getDomainInitSize(varindex) + i] = MAX_COST;
                                }
                                break;
                            default:
                                cerr << "Sorry operator " << op << " not implemented in element constraint!" << endl;
                                throw WrongFileFormat();
                            }
                        }
                        problem->postBinaryConstraint(varvalue, varindex, costs);
                    }
                }
            }
        } else {
            for (unsigned int i=0; i<vars.size(); i++) {
                vector<Cost> costs(problem->getDomainInitSize(varindex) * problem->getDomainInitSize(vars[i]), MIN_COST);
                assert(varindex != vars[i]);
                assert(varvalue != vars[i]);
                if (problem->canbe(varindex, (Value)i)) {
                    for (unsigned int b=0; b < problem->getDomainInitSize(vars[i]); b++) {
                        switch (op) {
                        case OrderType::IN:
                        case OrderType::EQ:
                            if (problem->toValue(vars[i], b) != (Value)i) {
                                costs[i * problem->getDomainInitSize(vars[i]) + b] = MAX_COST;
                            }
                            break;
                        case OrderType::NE:
                            if (problem->toValue(vars[i], b) == (Value)i) {
                                costs[i * problem->getDomainInitSize(vars[i]) + b] = MAX_COST;
                            }
                            break;
                        case OrderType::LE:
                            if (problem->toValue(vars[i], b) > (Value)i) {
                                costs[i * problem->getDomainInitSize(vars[i]) + b] = MAX_COST;
                            }
                            break;
                        case OrderType::LT:
                            if (problem->toValue(vars[i], b) >= (Value)i) {
                                costs[i * problem->getDomainInitSize(vars[i]) + b] = MAX_COST;
                            }
                            break;
                        case OrderType::GE:
                            if (problem->toValue(vars[i], b) < (Value)i) {
                                costs[i * problem->getDomainInitSize(vars[i]) + b] = MAX_COST;
                            }
                            break;
                        case OrderType::GT:
                            if (problem->toValue(vars[i], b) <= (Value)i) {
                                costs[i * problem->getDomainInitSize(vars[i]) + b] = MAX_COST;
                            }
                            break;
                        default:
                            cerr << "Sorry operator " << op << " not implemented in element constraint!" << endl;
                            throw WrongFileFormat();
                        }
                    }
                    problem->postBinaryConstraint(varindex, vars[i], costs);
                }
            }
        }
    }

    void buildConstraintElement(string id, OrderType op, vector<XVariable *> &list, int startIndex, XVariable *index, RankType rank, string value) {
        vector<int> vars;
        toMyVariables(list,vars);
        assert(startIndex == 0);
        assert(rank == RankType::ANY);
        int varindex = getMyVar(index->id);
        int varvalue = getMyVar(value);
        buildConstraintElement(op, vars, varindex, varvalue);
    }

    void buildConstraintElement(string id, vector<XVariable *> &list, int startIndex, XVariable *index, RankType rank, XVariable *value) override {
        buildConstraintElement(id, OrderType::EQ, list, startIndex, index, rank, value->id);
    }

    void buildConstraintElement(string id, vector<XVariable *> &list, int startIndex, XVariable *index, RankType rank, int minvalue, int maxvalue) {
        vector<int> vars;
        toMyVariables(list,vars);
        assert(startIndex == 0);
        assert(rank == RankType::ANY);
        int varindex = getMyVar(index->id);
        assert(problem->getDomainInitSize(varindex) == vars.size());
        for (unsigned int i=0; i<vars.size(); i++) {
            vector<Cost> costs(problem->getDomainInitSize(varindex) * problem->getDomainInitSize(vars[i]), MIN_COST);
            assert(varindex != vars[i]);
            assert(problem->toValue(varindex, i) == (Value)i);
            for (unsigned int b=0; b < problem->getDomainInitSize(vars[i]); b++) {
                if (minvalue > problem->toValue(vars[i], b) || maxvalue < problem->toValue(vars[i], b)) {
                    costs[i * problem->getDomainInitSize(vars[i]) + b] = MAX_COST;
                }
            }
            problem->postBinaryConstraint(varindex, vars[i], costs);
        }
    }

    void buildConstraintElement(string id, OrderType op, vector<XVariable *> &list, int startIndex, XVariable *index, RankType rank, int value) {
        vector<int> vars;
        toMyVariables(list,vars);
        assert(startIndex == 0);
        assert(rank == RankType::ANY);
        int varindex = getMyVar(index->id);
        assert(problem->getDomainInitSize(varindex) == vars.size());
        for (unsigned int i=0; i<vars.size(); i++) {
            vector<Cost> costs(problem->getDomainInitSize(varindex) * problem->getDomainInitSize(vars[i]), MIN_COST);
            assert(varindex != vars[i]);
            assert(problem->toValue(varindex, i) == (Value)i);
            for (unsigned int b=0; b < problem->getDomainInitSize(vars[i]); b++) {
                switch (op) {
                case OrderType::EQ:
                    if (problem->toValue(vars[i], b) != value) {
                        costs[i * problem->getDomainInitSize(vars[i]) + b] = MAX_COST;
                    }
                    break;
                case OrderType::NE:
                    if (problem->toValue(vars[i], b) == value) {
                        costs[i * problem->getDomainInitSize(vars[i]) + b] = MAX_COST;
                    }
                    break;
                case OrderType::LE:
                    if (problem->toValue(vars[i], b) > value) {
                        costs[i * problem->getDomainInitSize(vars[i]) + b] = MAX_COST;
                    }
                    break;
                case OrderType::LT:
                    if (problem->toValue(vars[i], b) >= value) {
                        costs[i * problem->getDomainInitSize(vars[i]) + b] = MAX_COST;
                    }
                    break;
                case OrderType::GE:
                    if (problem->toValue(vars[i], b) < value) {
                        costs[i * problem->getDomainInitSize(vars[i]) + b] = MAX_COST;
                    }
                    break;
                case OrderType::GT:
                    if (problem->toValue(vars[i], b) <= value) {
                        costs[i * problem->getDomainInitSize(vars[i]) + b] = MAX_COST;
                    }
                    break;
                default:
                    cerr << "Sorry operator " << op << " not implemented in element constraint!" << endl;
                    throw WrongFileFormat();
                }
            }
            problem->postBinaryConstraint(varindex, vars[i], costs);
        }
    }

    void buildConstraintElement(string id, vector<XVariable *> &list, int startIndex, XVariable *index, RankType rank, int value) override {
        vector<int> vars;
        toMyVariables(list,vars);
        assert(startIndex == 0);
        assert(rank == RankType::ANY);
        int varindex = getMyVar(index->id);
        assert(problem->getDomainInitSize(varindex) == vars.size());
        for (unsigned int i=0; i<vars.size(); i++) {
            vector<Cost> costs(problem->getDomainInitSize(varindex) * problem->getDomainInitSize(vars[i]), MIN_COST);
            assert(varindex != vars[i]);
            assert(problem->toValue(varindex, i) == (Value)i);
            for (unsigned int b=0; b < problem->getDomainInitSize(vars[i]); b++) {
                if (value != problem->toValue(vars[i], b)) {
                    costs[i * problem->getDomainInitSize(vars[i]) + b] = MAX_COST;
                }
            }
            problem->postBinaryConstraint(varindex, vars[i], costs);
        }
    }

    void buildConstraintElement(string id, vector<int> &list, int startIndex, XVariable *index, RankType rank, XVariable *value) override {
        assert(startIndex == 0);
        assert(rank == RankType::ANY);
        int varindex = getMyVar(index->id);
        int varvalue = getMyVar(value->id);
        assert(varindex != varvalue);
        assert(problem->getDomainInitSize(varindex) == list.size());
        vector<Cost> costs(problem->getDomainInitSize(varvalue) * problem->getDomainInitSize(varindex), MIN_COST);
        for (unsigned int a=0; a < problem->getDomainInitSize(varvalue); a++) {
            for (unsigned int b=0; b < problem->getDomainInitSize(varindex); b++) {
                if (problem->toValue(varvalue, a) != list[b]) {
                    costs[a * problem->getDomainInitSize(varindex) + b] = MAX_COST;
                }
            }
        }
        problem->postBinaryConstraint(varvalue, varindex, costs);
    }

    void buildConstraintElement(string id, vector<XVariable *> &list, XVariable *value) override {
        int varvalue = getMyVar(value->id);
        vector<int> vars;
        toMyVariables(list,vars);
        vector<int> myvars = vars;
        myvars.push_back(varvalue);
        for (unsigned int a=0; a < problem->getDomainInitSize(varvalue); a++) {
            string params = to_string(1);
            for (unsigned int i=0; i < vars.size(); i++) {
                params += " 1 " + to_string(problem->toValue(vars[i], a)) + " 1";
            }
            params += " " + to_string(problem->getDomainInitSize(varvalue) - 1);
            for (unsigned int b=0; b < problem->getDomainInitSize(varvalue); b++) {
                if (a != b) {
                    params += " " + to_string(problem->toValue(varvalue, b))  + " -1";
                }
            }
            problem->postKnapsackConstraint(myvars, params, false, true, false);
        }
    }

    void buildConstraintElement(string id, vector<XVariable *> &list, int value) override {
        vector<int> vars;
        toMyVariables(list,vars);
        string params = to_string(1);
        for (unsigned int i=0; i < vars.size(); i++) {
            params += " 1 " + to_string(value) + " 1";
        }
        problem->postKnapsackConstraint(vars, params, false, true, false);
    }

    void buildConstraintElement(string id, vector<XVariable *> &list, XVariable *index, int startIndex, XCondition &xc) override {
            switch (xc.operandType) {
            case OperandType::INTEGER:
                buildConstraintElement(id, xc.op, list, startIndex, index, RankType::ANY, xc.val);
                break;
            case OperandType::VARIABLE:
                buildConstraintElement(id, xc.op, list, startIndex, index, RankType::ANY, xc.var);
                break;
            case OperandType::INTERVAL:
                assert(xc.op == OrderType::IN);
                buildConstraintElement(id, list, startIndex, index, RankType::ANY, xc.min, xc.max);
                break;
            default:
                cerr << "Sorry operand type " << xc.operandType << " not implemented in element constraint!" << endl;
                throw WrongFileFormat();
            }
    }

    void buildConstraintElement(string id, vector<vector<XVariable*> > &matrix, int startRowIndex, XVariable *rowIndex, int startColIndex, XVariable* colIndex, XVariable* value) override {
        vector<vector<int> > vars;
        for (unsigned int i=0; i<matrix.size(); i++) {
            vars.push_back(vector<int>());
            toMyVariables(matrix[i],vars[i]);
        }
        assert(startRowIndex == 0);
        assert(startColIndex == 0);
        int varrowindex = getMyVar(rowIndex->id);
        int varcolindex = getMyVar(colIndex->id);
        int varvalue = getMyVar(value->id);
        assert(varcolindex != varrowindex);
        assert(varrowindex != varvalue);
        assert(varcolindex != varvalue);
        assert(problem->getDomainInitSize(varrowindex) == vars.size());
        for (unsigned int i=0; i<vars.size(); i++) {
            assert(problem->getDomainInitSize(varcolindex) == vars[i].size());
            assert(problem->toValue(varrowindex, i) == (Value)i);
            for (unsigned int j=0; j<vars[i].size(); j++) {
                assert(varrowindex != vars[i][j]);
                assert(varcolindex != vars[i][j]);
                assert(varvalue != vars[i][j]);
                assert(problem->toValue(varcolindex, j) == (Value)j);
                vector<int> scope;
                scope.push_back(varvalue);
                scope.push_back(varrowindex);
                scope.push_back(varcolindex);
                scope.push_back(vars[i][j]);
                int naryctr = problem->postNaryConstraintBegin(scope, MIN_COST, problem->getDomainInitSize(varvalue) * problem->getDomainInitSize(vars[i][j]));
                for (unsigned int a=0; a < problem->getDomainInitSize(varvalue); a++) {
                    for (unsigned int b=0; b < problem->getDomainInitSize(vars[i][j]); b++) {
                        if (problem->toValue(varvalue, a) != problem->toValue(vars[i][j], b)) {
                            vector<Value> tuple;
                            tuple.push_back(problem->toValue(varvalue, a));
                            tuple.push_back(problem->toValue(varrowindex, i));
                            tuple.push_back(problem->toValue(varcolindex, j));
                            tuple.push_back(problem->toValue(vars[i][j], b));
                            problem->postNaryConstraintTuple(naryctr, tuple, MAX_COST);
                        }
                    }
                }
                problem->postNaryConstraintEnd(naryctr);
            }
        }
    }

    void buildConstraintElement(string id, vector<vector<XVariable*> > &matrix, int startRowIndex, XVariable *rowIndex, int startColIndex, XVariable* colIndex, int value) override {
        vector<vector<int> > vars;
        for (unsigned int i=0; i<matrix.size(); i++) {
            vars.push_back(vector<int>());
            toMyVariables(matrix[i],vars[i]);
        }
        assert(startRowIndex == 0);
        assert(startColIndex == 0);
        int varrowindex = getMyVar(rowIndex->id);
        int varcolindex = getMyVar(colIndex->id);
        assert(varcolindex != varrowindex);
        assert(problem->getDomainInitSize(varrowindex) == vars.size());
        for (unsigned int i=0; i<vars.size(); i++) {
            assert(problem->getDomainInitSize(varcolindex) == vars[i].size());
            assert(problem->toValue(varrowindex, i) == (Value)i);
            for (unsigned int j=0; j<vars[i].size(); j++) {
                assert(varrowindex != vars[i][j]);
                assert(varcolindex != vars[i][j]);
                assert(problem->toValue(varcolindex, j) == (Value)j);
                vector<Cost> costs(problem->getDomainInitSize(varrowindex) * problem->getDomainInitSize(varcolindex) * problem->getDomainInitSize(vars[i][j]), MIN_COST);
                for (unsigned int b=0; b < problem->getDomainInitSize(vars[i][j]); b++) {
                    if (value != problem->toValue(vars[i][j], b)) {
                        costs[i * problem->getDomainInitSize(varcolindex) * problem->getDomainInitSize(vars[i][j]) + j * problem->getDomainInitSize(vars[i][j]) + b] = MAX_COST;
                    }
                }
                problem->postTernaryConstraint(varrowindex, varcolindex, vars[i][j], costs);
            }
        }
    }

    void buildConstraintElement(string id, vector<vector<int> > &matrix, int startRowIndex, XVariable *rowIndex, int startColIndex, XVariable* colIndex, XVariable *value) override {
        assert(startRowIndex == 0);
        assert(startColIndex == 0);
        int varrowindex = getMyVar(rowIndex->id);
        int varcolindex = getMyVar(colIndex->id);
        int varvalue = getMyVar(value->id);
        assert(varcolindex != varrowindex);
        assert(varrowindex != varvalue);
        assert(varcolindex != varvalue);
        assert(problem->getDomainInitSize(varrowindex) == matrix.size());
        for (unsigned int i=0; i<matrix.size(); i++) {
            assert(problem->getDomainInitSize(varcolindex) == matrix[i].size());
            assert(problem->toValue(varrowindex, i) == (Value)i);
            for (unsigned int j=0; j<matrix[i].size(); j++) {
                assert(varrowindex != matrix[i][j]);
                assert(varcolindex != matrix[i][j]);
                assert(problem->toValue(varcolindex, j) == (Value)j);
                vector<Cost> costs(problem->getDomainInitSize(varrowindex) * problem->getDomainInitSize(varcolindex) * problem->getDomainInitSize(varvalue), MIN_COST);
                for (unsigned int a=0; a < problem->getDomainInitSize(varvalue); a++) {
                    if (matrix[i][j] != problem->toValue(varvalue, a)) {
                        costs[i * problem->getDomainInitSize(varcolindex) * problem->getDomainInitSize(varvalue) + j * problem->getDomainInitSize(varvalue) + a] = MAX_COST;
                    }
                }
                problem->postTernaryConstraint(varrowindex, varcolindex, varvalue, costs);
            }
        }
    }

    // Semantic of subcircuit (successor variables may be not inserted in the circuit by selecting them-self, i.e. self-loops are authorized)
    void buildConstraintCircuit(string id, vector<XVariable *> &list, int startIndex) override {
        vector<int> succ;
        toMyVariables(list,succ);
        assert(startIndex == 0);
        unsigned int n = succ.size();
        buildConstraintAlldifferent(succ);
        vector<int> order;
        for (unsigned int i = 0; i < n; i++) {
            string extravarname = IMPLICIT_VAR_TAG + "order" + to_string(problem->numberOfVariables());
            int extravar = problem->makeEnumeratedVariable(extravarname, 0, n-1);
            mapping[extravarname] = extravar;
            order.push_back(extravar);
        }
        buildConstraintAlldifferent(order);
        //buildConstraintPrimitive(OrderType::EQ, order[0], 0); // only if circuit of length n
        //forall(i in range(n))(order[x[i]] = if order[i] = n-1 then 0 else order[i] + 1 endif
        for (unsigned int i = 0; i < n; i++) {
            for (unsigned int j = 0; j < n; j++) {
                if (i != j) {
                    vector<Cost> costs;
                    for (unsigned int a = 0; a < problem->getDomainInitSize(succ[i]); a++) {
                        for (unsigned int b = 0; b < problem->getDomainInitSize(order[i]); b++) {
                            for (unsigned int c = 0; c < problem->getDomainInitSize(order[j]); c++) {
                                if ((a != j) || (a==j && b==n-1 && c==0) || (a==j && b<n-1 && c>=b+1)) {
                                    costs.push_back(MIN_COST);
                                } else {
                                    costs.push_back(MAX_COST);
                                }
                            }
                        }
                    }
                    problem->postTernaryConstraint(succ[i], order[i], order[j], costs);
                }
                if (i < j) {
                    vector<Cost> costs;
                    for (unsigned int a = 0; a < problem->getDomainInitSize(succ[i]); a++) {
                        for (unsigned int b = 0; b < problem->getDomainInitSize(succ[j]); b++) {
                            if ((a!=j && b!=i) || (a==j && b!=j) || (b==i && a!=i)) {
                                costs.push_back(MIN_COST);
                            } else {
                                costs.push_back(MAX_COST);
                            }
                        }
                    }
                    problem->postBinaryConstraint(succ[i], succ[j], costs);
                }
            }
            //buildConstraintPrimitive(OrderType::NE, succ[i], i); // only if circuit of length n
        }
    }

    void buildUnaryCostFunction(Value mult, int var) {
        unsigned int domsize = problem->getDomainInitSize(var);
        vector<Cost> costs;
        Cost negcost = min(MIN_COST, (Cost)min((Cost)problem->toValue(var, 0) * mult, (Cost)problem->toValue(var, domsize-1) * mult));
        for (unsigned int a=0; a < domsize; a++) {
            costs.push_back((Cost)problem->toValue(var, a) * mult - negcost);
        }
        if (negcost < MIN_COST) {
            problem->decreaseLb(-negcost);
        }
        problem->postUnaryConstraint(var, costs);
    }

    void buildUnaryCostFunction(Value mult, XVariable *x) {
        int var = getMyVar(x->id);
        buildUnaryCostFunction(mult, var);
    }

    void buildObjectiveMinimizeVariable(XVariable *x) override {
        buildUnaryCostFunction(1, x);
    }

    void buildObjectiveMaximizeVariable(XVariable *x) override {
        ToulBar2::costMultiplier *= -1.0;
        buildUnaryCostFunction(-1, x);
    }

    void buildCostFunction(int coef, Tree *tree) {
        vector<string> list = tree->listOfVariables;
        vector<int> vars;
        toMyVariables(list,vars);
        vector<vector<int> > tuples;
        vector<Cost> tcosts;
        Cost negcost = MIN_COST;
        map<string, int> tuplemap;
        vector<int> tuplevec;
        unsigned int pos = 0;
        for (unsigned int i=0; i<vars.size(); i++) {
            Value val = problem->getInf(vars[i]);
            tuplemap[list[i]] = val;
            tuplevec.push_back(val);
        }
        do {
            tuples.push_back(tuplevec);
            Cost cost = (Cost)coef * tree->evaluate(tuplemap);
            tcosts.push_back(cost);
            if (cost < negcost) {
                negcost = cost;
            }
            while (pos < vars.size() && tuplevec[pos] == problem->getSup(vars[pos])) {
                pos++;
            }
            if (pos < vars.size()) {
                Value nextval = problem->nextValue(vars[pos], tuplevec[pos]);
                assert(nextval != tuplevec[pos]);
                tuplemap[list[pos]] = nextval;
                tuplevec[pos] = nextval;
                for (unsigned int i=0; i<pos; i++) {
                    Value val = problem->getInf(vars[i]);
                    tuplemap[list[i]] = val;
                    tuplevec[i] = val;
                }
                pos = 0;
            }
        } while (pos < vars.size());
        assert(tuples.size() == tcosts.size());
        if (negcost < MIN_COST) {
            problem->decreaseLb(-negcost);
        }
        if (vars.size()==1) {
            vector<Cost> costs(problem->getDomainInitSize(vars[0]), MAX_COST);
            for(unsigned int i=0; i<tuples.size(); i++) {
                costs[problem->toIndex(vars[0], tuples[i][0])] = tcosts[i] - negcost;
            }
            problem->postUnaryConstraint(vars[0], costs);
        } else if (vars.size()==2) {
            vector<Cost> costs(problem->getDomainInitSize(vars[0]) * problem->getDomainInitSize(vars[1]), MAX_COST);
            for(unsigned int i=0; i<tuples.size(); i++) {
                costs[problem->toIndex(vars[0], tuples[i][0]) * problem->getDomainInitSize(vars[1]) + problem->toIndex(vars[1], tuples[i][1])] = tcosts[i] - negcost;
            }
            problem->postBinaryConstraint(vars[0], vars[1], costs);
        } else if (vars.size()==3) {
            vector<Cost> costs(problem->getDomainInitSize(vars[0]) * problem->getDomainInitSize(vars[1]) * problem->getDomainInitSize(vars[2]), MAX_COST);
            for(unsigned int i=0; i<tuples.size(); i++) {
                costs[problem->toIndex(vars[0], tuples[i][0]) * problem->getDomainInitSize(vars[1]) * problem->getDomainInitSize(vars[2]) + problem->toIndex(vars[1], tuples[i][1]) * problem->getDomainInitSize(vars[2]) + problem->toIndex(vars[2], tuples[i][2])] = tcosts[i] - negcost;
            }
            problem->postTernaryConstraint(vars[0], vars[1], vars[2], costs);
        } else {
            int ctridx = problem->postNaryConstraintBegin(vars, MAX_COST, tuples.size());
            for(unsigned int i=0; i<tuples.size(); i++) {
                problem->postNaryConstraintTuple(ctridx, tuples[i], tcosts[i] - negcost);
            }
            problem->postNaryConstraintEnd(ctridx);
        }
    }

    void buildObjective(Cost sign, ExpressionObjective type, vector<XVariable *> &list, vector<int> &coefs) {
        assert(list.size() == coefs.size());
        vector<int> vars;
        toMyVariables(list,vars);
        Cost negcost = MIN_COST;
        vector<WeightedVarValPair> weightFunction;
        switch (type) {
        case ExpressionObjective::PRODUCT_O:
            if (list.size() == 2) {
                for (unsigned int a=0; a < problem->getDomainInitSize(vars[0]); a++) {
                    for (unsigned int b=0; b < problem->getDomainInitSize(vars[1]); b++) {
                        Cost cost = (Cost)coefs[0] * problem->toValue(vars[0], a) * coefs[1] * problem->toValue(vars[1], b) * sign;
                        if (cost < negcost) {
                            negcost = cost;
                        }
                    }
                }
                vector<Cost> costs;
                for (unsigned int a=0; a < problem->getDomainInitSize(vars[0]); a++) {
                    for (unsigned int b=0; b < problem->getDomainInitSize(vars[1]); b++) {
                        costs.push_back((Cost)coefs[0] * problem->toValue(vars[0], a) * coefs[1] * problem->toValue(vars[1], b) * sign - negcost);
                    }
                }
                problem->postBinaryConstraint(vars[0], vars[1], costs);
                break;
            } else if (list.size() >= 3) {
                cerr << "Sorry multiplicative objective for " << list.size() << " variables not implemented!" << endl;
                throw WrongFileFormat();
            }
        case ExpressionObjective::EXPRESSION_O:
            assert(list.size() == 1);
        case ExpressionObjective::SUM_O:
            for (unsigned int i=0; i<list.size(); i++) {
                buildUnaryCostFunction(sign * coefs[i], list[i]);
            }
            break;
        case ExpressionObjective::MINIMUM_O:
            sign *= -UNIT_COST;
        case ExpressionObjective::MAXIMUM_O:
            for (unsigned int i=0; i<list.size(); i++) {
                for (unsigned int a=0; a < problem->getDomainInitSize(vars[i]); a++) {
                    Cost cost = (Cost)problem->toValue(vars[i], a) * coefs[i] * sign;
                    if (cost < negcost) {
                        negcost = cost;
                    }
                }
            }
            for (unsigned int i=0; i<list.size(); i++) {
                for (unsigned int a=0; a < problem->getDomainInitSize(vars[i]); a++) {
                    WeightedVarValPair elt;
                    elt.varIndex = vars[i];
                    elt.val = problem->toValue(vars[i], a);
                    elt.weight = (Cost)elt.val * coefs[i] * sign - negcost;
                    weightFunction.push_back(elt);
                }
            }
            if (negcost < MIN_COST) {
                problem->decreaseLb(-negcost);
            }
            problem->postMaxWeight(vars.data(), vars.size(), "val", "DAG", MIN_COST, weightFunction);
            break;
        case ExpressionObjective::NVALUES_O:
        case ExpressionObjective::LEX_O:
        default:
            cerr << "Sorry objective type " << type << " not implemented!" << endl;
            throw WrongFileFormat();
        }
    }

    void buildObjectiveMinimize(ExpressionObjective type, vector<XVariable *> &list, vector<int> &coefs) override {
        buildObjective(UNIT_COST, type, list, coefs);
    }

    void buildObjectiveMaximize(ExpressionObjective type, vector<XVariable *> &list, vector<int> &coefs) override {
        ToulBar2::costMultiplier *= -1.0;
        buildObjective(-UNIT_COST, type, list, coefs);
    }

    void buildObjectiveMinimize(ExpressionObjective type, vector<XVariable *> &list) override {
        vector<int> coefs(list.size(), 1);
        buildObjectiveMinimize(type, list, coefs);
    }

    void buildObjectiveMaximize(ExpressionObjective type, vector<XVariable *> &list) override {
        vector<int> coefs(list.size(), 1);
        buildObjectiveMaximize(type, list, coefs);
    }

    void buildObjective(Cost sign, ExpressionObjective type, vector<Tree *> &trees, vector<int> &coefs) {
        assert(trees.size() == coefs.size());
        vector<int> vars;
        vector<WeightedVarValPair> weightFunction;
        string varargmaxname;
        int varargmax;
        string varmaxname;
        int varmax;
        XCondition cond;
        bool themax = true;
        switch (type) {
        case ExpressionObjective::EXPRESSION_O:
            assert(trees.size() == 1);
        case ExpressionObjective::SUM_O:
            for (unsigned int i=0; i<trees.size(); i++) {
                buildCostFunction(sign * coefs[i], trees[i]);
            }
            break;
        case ExpressionObjective::MINIMUM_O:
            themax = false;
        case ExpressionObjective::MAXIMUM_O:
            for (unsigned int i=0; i<trees.size(); i++) {
                int var = buildConstraintIntension(trees[i]);
                if (coefs[i] != 1) {
                    string prodname = IMPLICIT_VAR_TAG + "prod" + to_string(problem->numberOfVariables());
                    int varprod = problem->makeEnumeratedVariable(prodname, min(problem->getInf(var) * coefs[i], problem->getSup(var) * coefs[i]), max(problem->getInf(var) * coefs[i], problem->getSup(var) * coefs[i]));
                    mapping[prodname] = varprod;
                    buildConstraintPrimitiveMult(OrderType::EQ, var, coefs[i], varprod);
                    vars.push_back(varprod);
                } else {
                    vars.push_back(var);
                }
            }
            assert(vars.size() == coefs.size());
            varmaxname = IMPLICIT_VAR_TAG + ((themax)?"max":"min") + to_string(problem->numberOfVariables());
            varmax = problem->makeEnumeratedVariable(varmaxname, 0, vars.size()-1);
            mapping[varmaxname] = varmax;
            if ((sign==UNIT_COST && themax) || (sign==-UNIT_COST && !themax)) {
                for (unsigned int i=0; i<vars.size(); i++) {
                    buildConstraintPrimitive((themax)?OrderType::LE:OrderType::GE, vars[i], 0, varmax);
                }
            } else {
                varargmaxname = IMPLICIT_VAR_TAG + ((themax)?"argmax":"argmin") + to_string(problem->numberOfVariables());
                varargmax = problem->makeEnumeratedVariable(varargmaxname, 0, vars.size()-1);
                mapping[varargmaxname] = varargmax;
                cond.operandType = OperandType::VARIABLE;
                cond.op = OrderType::EQ;
                cond.var = varmaxname;
                buildConstraintMinMax(themax, vars, varargmax, cond);
            }
            buildUnaryCostFunction(sign, varmax);
            break;
        case ExpressionObjective::PRODUCT_O:
        case ExpressionObjective::NVALUES_O:
        case ExpressionObjective::LEX_O:
        default:
            cerr << "Sorry objective type " << type << " on trees not implemented!" << endl;
            throw WrongFileFormat();
        }
    }

    void buildObjectiveMinimize(ExpressionObjective type, vector<Tree *> &trees, vector<int> &coefs) override {
        buildObjective(UNIT_COST, type, trees, coefs);
    }

    void buildObjectiveMaximize(ExpressionObjective type, vector<Tree *> &trees, vector<int> &coefs) override {
        ToulBar2::costMultiplier *= -1.0;
        buildObjective(-UNIT_COST, type, trees, coefs);
    }

    void buildObjectiveMinimize(ExpressionObjective type, vector<Tree *> &trees) override {
        vector<int> coefs(trees.size(), 1);
        buildObjectiveMinimize(type, trees, coefs);
    }

    void buildObjectiveMaximize(ExpressionObjective type, vector<Tree *> &trees) override {
        vector<int> coefs(trees.size(), 1);
        buildObjectiveMaximize(type, trees, coefs);
    }

    void buildObjectiveMinimizeExpression(string expr) override {
        Tree tree(expr);
        vector<Tree *> trees;
        trees.push_back(&tree);
        buildObjectiveMinimize(ExpressionObjective::EXPRESSION_O, trees);
    }

    void buildObjectiveMaximizeExpression(string expr) override {
        Tree tree(expr);
        vector<Tree *> trees;
        trees.push_back(&tree);
        buildObjectiveMaximize(ExpressionObjective::EXPRESSION_O, trees);
    }

    void buildConstraintInstantiation(string id, vector<XVariable *> &list, vector<int> &values) override {
        vector<int> vars;
        toMyVariables(list,vars);
        problem->assignLS(vars, values);
    }

    void buildAnnotationDecision(vector<XVariable *> &list) override {}
};

#endif /* SRC_XCSP3_XMLCSP3_H_ */

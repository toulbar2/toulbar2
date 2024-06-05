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
        vector<int> assignedVars;
        vector<Value> assignedValues;

    virtual ~MySolverCallbacks() {}

    const vector<string> orderString = {"le", "lt", "ge", "gt", "eq", "eq", "ne"};
    const vector<string> lexString = {"le", "le", "ge", "ge", "eq", "eq", "ne"};

    void beginInstance(InstanceType type) override {
        ToulBar2::xmlcop = false;
        XCSP3CoreCallbacks::intensionUsingString = false;
        XCSP3CoreCallbacks::recognizeSpecialIntensionCases = true;
        XCSP3CoreCallbacks::recognizeSpecialCountCases = false;
        XCSP3CoreCallbacks::recognizeNValuesCases = true;
        problem->updateUb(MAX_COST);
    }

    void endInstance() override {
        if (!ToulBar2::xmlcop) problem->updateUb(UNIT_COST);
        problem->sortConstraints();
        if (assignedVars.size() > 0) {
            problem->assignLS(assignedVars, assignedValues);
        }
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
        if (mapping.find(var) == mapping.end()) {
            buildVariableInteger(var, std::stoi(var), std::stoi(var));
        }
        assert(mapping.find(var) != mapping.end());
        assert(mapping[var] >= 0 && mapping[var] < (int)problem->numberOfVariables());
        return mapping[var];
    }

    int getMyVar(XVariable *var) {
        return getMyVar(var->id);
    }

    // transforms a vector of XVariable in vector of toulbar2 variable indices and add it to dest (assuming only one occurrence of each variable)
    void toMyVariables(vector<XVariable*> &src, vector<int> &dest) {
        //set<int> control;
#ifndef NDEBUG
        size_t initsize = dest.size();
#endif
        for(unsigned int i = 0;i<src.size();i++) {
            dest.push_back(getMyVar(src[i]));
            //control.insert(getMyVar(src[i]));
        }
        assert(dest.size() > initsize);
        //assert(dest.size() - initsize == control.size());
    }

    // transforms a vector of string Variable name in vector of toulbar2 variable indices and add it to dest
    void toMyVariables(vector<string> &src, vector<int> &dest) {
        //set<int> control;
        for(unsigned int i = 0;i<src.size();i++) {
            dest.push_back(getMyVar(src[i]));
            //control.insert(getMyVar(src[i]));
        }
        //assert(dest.size() == control.size());
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

    // bug correction: vector<vector<int> > &tuples generates a wrong solution on BeerJugs-table-01_c23.xml in MiniCOP23
    void buildConstraintExtension(vector<int> vars, vector<vector<int> > tuples, bool isSupport, bool hasStar) {
        if(hasStar) {
            vector<vector<int> > newtuples;
            for(auto& tuple:tuples) {
                recursiveExpanse(vars, tuple, 0, tuple[0], newtuples, vector<int>());
            }
            assert(newtuples.size() >= tuples.size());
            tuples = newtuples; // warning it will change the input tuples
        }

        // remove multiple occurrences of the same variable
        for (unsigned int i = 0; i < vars.size(); i++) {
            for (unsigned int j = i + 1; j < vars.size(); j++) {
                if (vars[i] == vars[j]) {
                    vector<vector<int> > newtuples;
                    for(auto& tuple:tuples) {
                           if (tuple[i] == tuple[j]) {
                               tuple.erase(tuple.begin()+j);
                               newtuples.push_back(tuple);
                           }
                    }
                    tuples = newtuples; // warning it will change the input tuples
                    vars.erase(vars.begin()+j);
                    j--;
                }
            }
        }

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
            vector<Cost> costs((size_t)problem->getDomainInitSize(vars[0]) * (size_t)problem->getDomainInitSize(vars[1]) * (size_t)problem->getDomainInitSize(vars[2]), (isSupport)?MAX_COST:MIN_COST);
            for(auto& tuple:tuples) {
                if (isSupport) {
//                    costs[(size_t)problem->toIndex(vars[0], tuple[0]) * problem->getDomainInitSize(vars[1]) * problem->getDomainInitSize(vars[2]) + (size_t)problem->toIndex(vars[1], tuple[1]) * problem->getDomainInitSize(vars[2]) + problem->toIndex(vars[2], tuple[2])] = MIN_COST;
                    costs[(size_t)problem->toIndex(vars[2], tuple[2]) * problem->getDomainInitSize(vars[0]) * problem->getDomainInitSize(vars[1]) + (size_t)problem->toIndex(vars[0], tuple[0]) * problem->getDomainInitSize(vars[1]) + problem->toIndex(vars[1], tuple[1])] = MIN_COST;
                } else {
//                    costs[(size_t)problem->toIndex(vars[0], tuple[0]) * problem->getDomainInitSize(vars[1]) * problem->getDomainInitSize(vars[2]) + (size_t)problem->toIndex(vars[1], tuple[1]) * problem->getDomainInitSize(vars[2]) + problem->toIndex(vars[2], tuple[2])] = MAX_COST;
                    costs[(size_t)problem->toIndex(vars[2], tuple[2]) * problem->getDomainInitSize(vars[0]) * problem->getDomainInitSize(vars[1]) + (size_t)problem->toIndex(vars[0], tuple[0]) * problem->getDomainInitSize(vars[1]) + problem->toIndex(vars[1], tuple[1])] = MAX_COST;
                }
            }
//            problem->postTernaryConstraint(vars[0], vars[1], vars[2], costs);
            problem->postTernaryConstraint(vars[2], vars[0], vars[1], costs);
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
        lastTuples = tuples; // copy tuples for future usage by buildConstraintExtensionAs before tuples expansion
        vector<string> thelist;
        for (XVariable *var : list) {
            thelist.push_back(var->id);
        }
        buildConstraintExtension(id, thelist, tuples, isSupport, hasStar);
    }

    void buildConstraintExtension(int var, vector<int> &tuples, bool isSupport, bool hasStar) {
        // This function is called for unary extensional constraint.
        assert(tuples[0] != STAR || tuples.size()==1);
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
        lastTuples = vector<vector<int> >();
        for (auto value:tuples) {
            vector<int> tuple;
            tuple.push_back(value);
            lastTuples.push_back(tuple);
        }
        buildConstraintExtension(id, variable->id, tuples, isSupport, hasStar);
    }

    // This function is called with group of constraint where the set of tuples is exactly the same
    // than the previous one (then, you can save time/memory using the same set of tuples.
    void buildConstraintExtensionAs(string id, vector<XVariable *> list, bool isSupport, bool hasStar) override {
        vector<int> vars;
        toMyVariables(list,vars);
        buildConstraintExtension(vars, lastTuples, isSupport, hasStar);
    }

    // returns a WCSP variable index corresponding to the evaluation of Tree expression
    int buildConstraintIntensionVar(Tree *tree) {
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
                for (unsigned int c = 0; c < problem->getDomainInitSize(varz); c++) {
                    for (unsigned int a = 0; a < problem->getDomainInitSize(varx); a++) {
                        for (unsigned int b = 0; b < problem->getDomainInitSize(vary); b++) {
                            costs.push_back((problem->toValue(varx, a) * problem->toValue(vary, b) == problem->toValue(varz, c))?MIN_COST:MAX_COST);
                        }
                    }
                }
//                problem->postTernaryConstraint(varx, vary, varz, costs);
                problem->postTernaryConstraint(varz, varx, vary, costs);
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

    // warning! vars may be modified
    void buildConstraintAlldifferent(vector<int> &vars) {
        // remove multiple occurrences of the same variable
        for (unsigned int i = 0; i < vars.size(); i++) {
            for (unsigned int j = i + 1; j < vars.size(); j++) {
                if (vars[i] == vars[j]) {
                    vars.erase(vars.begin()+j);
                    j--;
                }
            }
        }
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
            vars.push_back(buildConstraintIntensionVar(list[i]));
        }
        assert(vars.size() == list.size());
        buildConstraintAlldifferent(vars);
    }

    void buildConstraintAlldifferentList(string id, vector<vector<XVariable *>> &lists) override {
        for (unsigned int i = 0; i < lists.size(); i++) {
            for (unsigned int j = 0; j < i; j++) {
                assert(lists[i].size() == lists[j].size());
                vector<int> vars;
                vector<Long> coefs;
                for (unsigned int k = 0; k < lists[i].size(); k++) {
                    Tree tree("ne(" + lists[i][k]->id + "," + lists[j][k]->id + ")");
                    vars.push_back(buildConstraintIntensionVar(&tree));
                    coefs.push_back(1);
                }
                XCondition cond;
                cond.op = OrderType::GE;
                cond.val = 1;
                cond.operandType = OperandType::INTEGER;
                buildConstraintSum(vars, coefs, cond);
            }
        }
    }

    void buildConstraintAlldifferentMatrix(string id, vector<vector<XVariable *>> &matrix) override {
        for (unsigned int i = 0; i < matrix.size(); i++) {
            buildConstraintAlldifferent(id, matrix[i]);
        }
        for (unsigned int j = 0; j < matrix[0].size(); j++) {
            vector<XVariable *> vars;
            for (unsigned int i = 0; i < matrix[0].size(); i++) {
                vars.push_back(matrix[i][j]);
            }
            buildConstraintAlldifferent(id, vars);
        }
    }

    void buildConstraintAlldifferentExcept(vector<int> &vars, vector<int> &except) {
        // remove multiple occurrences of the same variable
        for (unsigned int i = 0; i < vars.size(); i++) {
            for (unsigned int j = i + 1; j < vars.size(); j++) {
                if (vars[i] == vars[j]) {
                    vars.erase(vars.begin()+j);
                    j--;
                }
            }
        }
        if (vars.size() <= 50) {
            for (unsigned int i = 0; i < vars.size(); i++) {
                for (unsigned int j = i+1; j < vars.size(); j++) {
                    vector<Cost> costs;
                    for (unsigned int a = 0; a < problem->getDomainInitSize(vars[i]); a++) {
                        for (unsigned int b = 0; b < problem->getDomainInitSize(vars[j]); b++) {
                            costs.push_back((problem->toValue(vars[i], a)!=problem->toValue(vars[j], b) || std::find(except.begin(), except.end(), problem->toValue(vars[i], a))!=except.end() || std::find(except.begin(), except.end(), problem->toValue(vars[j], b))!=except.end())?MIN_COST:MAX_COST);
                        }
                    }
                    problem->postBinaryConstraint(vars[i], vars[j], costs);
                }
            }
        } else {
            set<Value> values;
            for (unsigned int i = 0; i < vars.size(); i++) {
                ((EnumeratedVariable*)((WCSP*)problem)->getVar(vars[i]))->getDomain(values);
            }
            for (Value value : values) {
                if (std::find(except.begin(), except.end(), value)==except.end()) {
                    string params = to_string(-1);
                    for (unsigned int i = 0; i < vars.size(); i++) {
                        if (problem->canbe(vars[i], value)) {
                            params += " 1 " + to_string(value) + " -1";
                        } else {
                            params += " 0";
                        }
                    }
                    problem->postKnapsackConstraint(vars, params, false, true, false);
                }
            }
        }
    }

    void buildConstraintAlldifferentExcept(string id, vector<string> &list, vector<int> &except) {
        vector<int> vars;
        toMyVariables(list,vars);
        buildConstraintAlldifferentExcept(vars, except);
    }

    void buildConstraintAlldifferentExcept(string id, vector<XVariable *> &list, vector<int> &except) override {
        vector<string> slist;
        for (auto& x:list) {
            slist.push_back(x->id);
        }
        buildConstraintAlldifferentExcept(id, slist, except);
    }

    // warning! vars and coefs may be modified
    void buildConstraintSum(vector<int> &vars, vector<Long> &coefs, XCondition &cond, Long valLong = LONG_MAX, Long minLong = LONG_MAX, Long maxLong = LONG_MAX) {
        assert(vars.size() == coefs.size());
        // remove multiple occurrences of the same variable
        for (unsigned int i = 0; i < vars.size(); i++) {
            for (unsigned int j = i + 1; j < vars.size(); j++) {
                if (vars[i] == vars[j]) {
                    vars.erase(vars.begin()+j);
                    coefs[i] += coefs[j];
                    coefs.erase(coefs.begin()+j);
                    j--;
                }
            }
        }
        Long rightcoef = 0;
        string extravarname = "";
        int extravar = -1;
        string params = "";
        switch (cond.operandType) {
        case OperandType::VARIABLE:
            vars.push_back(getMyVar(cond.var));
            coefs.push_back(-1);
            for (unsigned int i = 0; i < vars.size()-1; i++) {
                if (vars[i] == vars.back()) {
                    vars.pop_back();
                    coefs[i] += coefs.back();
                    coefs.pop_back();
                    break;
                }
            }
            rightcoef -= ((valLong<LONG_MAX)?valLong:cond.val);
        case OperandType::INTEGER:
            rightcoef += ((valLong<LONG_MAX)?valLong:cond.val);
            switch (cond.op) {
            case OrderType::LE:
                params = to_string(-rightcoef);
                for (unsigned int i=0; i < vars.size(); i++) {
                    int domsize = problem->getDomainInitSize(vars[i]);
                    params += " " + to_string(domsize);
                    for (int idval=0; idval < domsize; idval++) {
                        int value = problem->toValue(vars[i], idval);
                        params += " " + to_string(value) + " " + to_string(-coefs[i] * value);
                    }
                }
                problem->postKnapsackConstraint(vars, params, false, true, false);
                break;
            case OrderType::LT:
                params = to_string(-rightcoef + 1);
                for (unsigned int i=0; i < vars.size(); i++) {
                    int domsize = problem->getDomainInitSize(vars[i]);
                    params += " " + to_string(domsize);
                    for (int idval=0; idval < domsize; idval++) {
                        int value = problem->toValue(vars[i], idval);
                        params += " " + to_string(value) + " " + to_string(-coefs[i] * value);
                    }
                }
                problem->postKnapsackConstraint(vars, params, false, true, false);
                break;
            case OrderType::GE:
                params = to_string(rightcoef);
                for (unsigned int i=0; i < vars.size(); i++) {
                    int domsize = problem->getDomainInitSize(vars[i]);
                    params += " " + to_string(domsize);
                    for (int idval=0; idval < domsize; idval++) {
                        int value = problem->toValue(vars[i], idval);
                        params += " " + to_string(value) + " " + to_string(coefs[i] * value);
                    }
                }
                problem->postKnapsackConstraint(vars, params, false, true, false);
                break;
            case OrderType::GT:
                params = to_string(rightcoef + 1);
                for (unsigned int i=0; i < vars.size(); i++) {
                    int domsize = problem->getDomainInitSize(vars[i]);
                    params += " " + to_string(domsize);
                    for (int idval=0; idval < domsize; idval++) {
                        int value = problem->toValue(vars[i], idval);
                        params += " " + to_string(value) + " " + to_string(coefs[i] * value);
                    }
                }
                problem->postKnapsackConstraint(vars, params, false, true, false);
                break;
            case OrderType::NE:
                extravarname = IMPLICIT_VAR_TAG + "sum" + to_string(problem->numberOfVariables());
                extravar = problem->makeEnumeratedVariable(extravarname, 0, 1);
                mapping[extravarname] = extravar;
                vars.push_back(extravar);
                coefs.push_back(INT_MAX);
                params = to_string(rightcoef + 1);
                for (unsigned int i=0; i < vars.size(); i++) {
                    int domsize = problem->getDomainInitSize(vars[i]);
                    params += " " + to_string(domsize);
                    for (int idval=0; idval < domsize; idval++) {
                        int value = problem->toValue(vars[i], idval);
                        params += " " + to_string(value) + " " + to_string(coefs[i] * value);
                    }
                }
                problem->postKnapsackConstraint(vars, params, false, true, false);
                params = to_string(-rightcoef + 1 - INT_MAX);
                for (unsigned int i=0; i < vars.size(); i++) {
                    int domsize = problem->getDomainInitSize(vars[i]);
                    params += " " + to_string(domsize);
                    for (int idval=0; idval < domsize; idval++) {
                        int value = problem->toValue(vars[i], idval);
                        params += " " + to_string(value) + " " + to_string(-coefs[i] * value);
                    }
                }
                problem->postKnapsackConstraint(vars, params, false, true, false);
                break;
            case OrderType::EQ:
                params = to_string(rightcoef);
                for (unsigned int i=0; i < vars.size(); i++) {
                    int domsize = problem->getDomainInitSize(vars[i]);
                    params += " " + to_string(domsize);
                    for (int idval=0; idval < domsize; idval++) {
                        int value = problem->toValue(vars[i], idval);
                        params += " " + to_string(value) + " " + to_string(coefs[i] * value);
                    }
                }
                problem->postKnapsackConstraint(vars, params, false, true, false);
                params = to_string(-rightcoef);
                for (unsigned int i=0; i < vars.size(); i++) {
                    int domsize = problem->getDomainInitSize(vars[i]);
                    params += " " + to_string(domsize);
                    for (int idval=0; idval < domsize; idval++) {
                        int value = problem->toValue(vars[i], idval);
                        params += " " + to_string(value) + " " + to_string(-coefs[i] * value);
                    }
                }
                problem->postKnapsackConstraint(vars, params, false, true, false);
                break;
            default:
                cerr << "Sorry operator " << cond.op << " not implemented in sum constraint!" << endl;
                throw WrongFileFormat();
            }
            break;
        case OperandType::INTERVAL:
            assert(cond.op == OrderType::IN);
            // sum >= cond.min
            params = to_string((minLong<LONG_MAX)?minLong:cond.min);
            for (unsigned int i=0; i < vars.size(); i++) {
                int domsize = problem->getDomainInitSize(vars[i]);
                params += " " + to_string(domsize);
                for (int idval=0; idval < domsize; idval++) {
                    int value = problem->toValue(vars[i], idval);
                    params += " " + to_string(value) + " " + to_string(coefs[i] * value);
                }
            }
            problem->postKnapsackConstraint(vars, params, false, true, false);
            // sum <= cond.max
            params = to_string(-((maxLong<LONG_MAX)?maxLong:cond.max));
            for (unsigned int i=0; i < vars.size(); i++) {
                int domsize = problem->getDomainInitSize(vars[i]);
                params += " " + to_string(domsize);
                for (int idval=0; idval < domsize; idval++) {
                    int value = problem->toValue(vars[i], idval);
                    params += " " + to_string(value) + " " + to_string(-coefs[i] * value);
                }
            }
            problem->postKnapsackConstraint(vars, params, false, true, false);
            break;
        default:
            cerr << "Sorry operandType " << cond.operandType << " not implemented in sum constraint!" << endl;
            throw WrongFileFormat();
        }
    }

    void buildConstraintSum(string id, vector<XVariable *> &list, vector<int> &coefs, XCondition &cond) override {
        vector<int> vars;
        toMyVariables(list,vars);
        vector<Long> mycoefs(coefs.begin(), coefs.end());
        buildConstraintSum(vars, mycoefs, cond);
    }

    void buildConstraintSum(string id, vector<XVariable *> &list, XCondition &cond) override {
        vector<int> coefs(list.size(), 1);
        buildConstraintSum(id, list, coefs, cond);
    }

    void buildConstraintSum(string id, vector<XVariable *> &list, vector<XVariable *> &coefs, XCondition &cond) override {
        vector<int> vars;
        toMyVariables(list,vars);
        vector<int> varscoefs;
        toMyVariables(coefs,varscoefs);
        assert(vars.size() == varscoefs.size());
        vector<Long> mycoefs(list.size(), 1);
        vector<int> varsprod;
        for (unsigned int i = 0; i < list.size(); i++) {
            string prodname = IMPLICIT_VAR_TAG + "prod" + to_string(problem->numberOfVariables());
            Value prodinf = min(min(problem->getInf(vars[i]) * problem->getInf(varscoefs[i]), problem->getInf(vars[i]) * problem->getSup(varscoefs[i])), min(problem->getSup(vars[i]) * problem->getInf(varscoefs[i]), problem->getSup(vars[i]) * problem->getSup(varscoefs[i])));
            Value prodsup = max(max(problem->getInf(vars[i]) * problem->getInf(varscoefs[i]), problem->getInf(vars[i]) * problem->getSup(varscoefs[i])), max(problem->getSup(vars[i]) * problem->getInf(varscoefs[i]), problem->getSup(vars[i]) * problem->getSup(varscoefs[i])));
            assert(prodinf <= prodsup);
            int varprod = problem->makeEnumeratedVariable(prodname, prodinf, prodsup);
            mapping[prodname] = varprod;
            buildConstraintMult(vars[i], varscoefs[i], varprod);
            varsprod.push_back(varprod);
        }
        assert(varsprod.size() == list.size());
        buildConstraintSum(varsprod, mycoefs, cond);

    }

    void buildConstraintSum(string id, vector<Tree *> &trees, vector<int> &coefs, XCondition &cond) override {
        vector<int> vars;
        for (unsigned int i = 0; i < trees.size(); i++) {
            vars.push_back(buildConstraintIntensionVar(trees[i]));
        }
        assert(vars.size() == trees.size());
        vector<Long> mycoefs(coefs.begin(), coefs.end());
        buildConstraintSum(vars, mycoefs, cond);
    }

    void buildConstraintSum(string id, vector<Tree *> &trees, XCondition &cond) override {
        vector<int> coefs(trees.size(), 1);
        buildConstraintSum(id, trees, coefs, cond);
    }

    void buildConstraintCount(vector<int> &vars, vector<int> &values, XCondition &cond) {
        int rightcoef = 0;
        string params = "";
        string params2 = "";
        int domsize;
        int nbval;
        switch (cond.operandType) {
            case OperandType::VARIABLE:
                vars.push_back(getMyVar(cond.var));
                switch (cond.op) {
                    case OrderType::LE:
                        params = to_string(0);
                        params2="";
                        for (unsigned int i=0; i < vars.size()-1; i++) {
                            domsize = problem->getDomainInitSize(vars[i]);
                            //params += " " + to_string(domsize);
                            nbval=0;
                            params2="";
                            for (int idval=0; idval < domsize; idval++) {
                                int value = problem->toValue(vars[i], idval);
                                auto it = find(values.begin(), values.end(), value);
                                if (it != values.end()) {
                                    nbval+=1;
                                    params2 += " " + to_string(value) + " " + to_string(-1);
                                }
                            }
                            params+=" "+to_string(nbval)+params2;
                        }
                        domsize=problem->getDomainInitSize(vars.back());
                        params += " " + to_string(domsize);
                        for (int idval=0; idval < domsize; idval++) {
                            int value = problem->toValue(vars.back(), idval);
                            params += " " + to_string(value) + " " + to_string(value);
                        }
                        problem->postKnapsackConstraint(vars, params, false, true, false);
                        break;
                    case OrderType::LT:
                        params = to_string(1);
                        params2="";
                        for (unsigned int i=0; i < vars.size()-1; i++) {
                            domsize = problem->getDomainInitSize(vars[i]);
                            //params += " " + to_string(domsize);
                            nbval=0;
                            params2="";
                            for (int idval=0; idval < domsize; idval++) {
                                int value = problem->toValue(vars[i], idval);
                                auto it = find(values.begin(), values.end(), value);
                                if (it != values.end()) {
                                    nbval+=1;
                                    params2 += " " + to_string(value) + " " + to_string(-1);
                                }
                            }
                            params+=" "+to_string(nbval)+params2;
                        }
                        domsize=problem->getDomainInitSize(vars.back());
                        params += " " + to_string(domsize);
                        for (int idval=0; idval < domsize; idval++) {
                            int value = problem->toValue(vars.back(), idval);
                            params += " " + to_string(value) + " " + to_string(value);
                        }
                        problem->postKnapsackConstraint(vars, params, false, true, false);
                        break;
                    case OrderType::GE:
                        params = to_string(0);
                        params2="";
                        for (unsigned int i=0; i < vars.size()-1; i++) {
                            domsize = problem->getDomainInitSize(vars[i]);
                            //params += " " + to_string(domsize);
                            nbval=0;
                            params2="";
                            for (int idval=0; idval < domsize; idval++) {
                                int value = problem->toValue(vars[i], idval);
                                auto it = find(values.begin(), values.end(), value);
                                if (it != values.end()) {
                                    nbval+=1;
                                    params2 += " " + to_string(value) + " " + to_string(1);
                                }
                            }
                            params+=" "+to_string(nbval)+params2;
                        }
                        domsize=problem->getDomainInitSize(vars.back());
                        params += " " + to_string(domsize);
                        for (int idval=0; idval < domsize; idval++) {
                            int value = problem->toValue(vars.back(), idval);
                            params += " " + to_string(value) + " " + to_string(-value);
                        }
                        problem->postKnapsackConstraint(vars, params, false, true, false);
                        break;
                    case OrderType::GT:
                        params = to_string(1);
                        params2="";
                        for (unsigned int i=0; i < vars.size()-1; i++) {
                            domsize = problem->getDomainInitSize(vars[i]);
                            //params += " " + to_string(domsize);
                            nbval=0;
                            params2="";
                            for (int idval=0; idval < domsize; idval++) {
                                int value = problem->toValue(vars[i], idval);
                                auto it = find(values.begin(), values.end(), value);
                                if (it != values.end()) {
                                    nbval+=1;
                                    params2 += " " + to_string(value) + " " + to_string(1);
                                }
                            }
                            params+=" "+to_string(nbval)+params2;
                        }
                        domsize=problem->getDomainInitSize(vars.back());
                        params += " " + to_string(domsize);
                        for (int idval=0; idval < domsize; idval++) {
                            int value = problem->toValue(vars.back(), idval);
                            params += " " + to_string(value) + " " + to_string(-value);
                        }
                        problem->postKnapsackConstraint(vars, params, false, true, false);
                        break;
                    case OrderType::NE:
                        params = to_string(1);
                        params2="";
                        for (unsigned int i=0; i < vars.size()-1; i++) {
                            domsize = problem->getDomainInitSize(vars[i]);
                            //params += " " + to_string(domsize);
                            nbval=0;
                            params2="";
                            for (int idval=0; idval < domsize; idval++) {
                                int value = problem->toValue(vars[i], idval);
                                auto it = find(values.begin(), values.end(), value);
                                if (it != values.end()) {
                                    nbval+=1;
                                    params2 += " " + to_string(value) + " " + to_string(1);
                                }
                            }
                            params+=" "+to_string(nbval)+params2;
                        }
                        domsize=problem->getDomainInitSize(vars.back());
                        params += " " + to_string(domsize);
                        for (int idval=0; idval < domsize; idval++) {
                            int value = problem->toValue(vars.back(), idval);
                            params += " " + to_string(value) + " " + to_string(-value);
                        }
                        problem->postKnapsackConstraint(vars, params, false, true, false);
                        params = to_string(1);
                        params2="";
                        for (unsigned int i=0; i < vars.size()-1; i++) {
                            domsize = problem->getDomainInitSize(vars[i]);
                            //params += " " + to_string(domsize);
                            nbval=0;
                            params2="";
                            for (int idval=0; idval < domsize; idval++) {
                                int value = problem->toValue(vars[i], idval);
                                auto it = find(values.begin(), values.end(), value);
                                if (it != values.end()) {
                                    nbval+=1;
                                    params2 += " " + to_string(value) + " " + to_string(-1);
                                }
                            }
                            params+=" "+to_string(nbval)+params2;
                        }
                        domsize=problem->getDomainInitSize(vars.back());
                        params += " " + to_string(domsize);
                        for (int idval=0; idval < domsize; idval++) {
                            int value = problem->toValue(vars.back(), idval);
                            params += " " + to_string(value) + " " + to_string(value);
                        }
                        problem->postKnapsackConstraint(vars, params, false, true, false);
                        break;
                    case OrderType::EQ:
                        params = to_string(0);
                        params2="";
                        for (unsigned int i=0; i < vars.size()-1; i++) {
                            domsize = problem->getDomainInitSize(vars[i]);
                            //params += " " + to_string(domsize);
                            nbval=0;
                            params2="";
                            params2="";
                            for (int idval=0; idval < domsize; idval++) {
                                int value = problem->toValue(vars[i], idval);
                                auto it = find(values.begin(), values.end(), value);
                                if (it != values.end()) {
                                    nbval+=1;
                                    params2 += " " + to_string(value) + " " + to_string(1);
                                }
                            }
                            params+=" "+to_string(nbval)+params2;
                        }
                        domsize=problem->getDomainInitSize(vars.back());
                        params += " " + to_string(domsize);
                        for (int idval=0; idval < domsize; idval++) {
                            int value = problem->toValue(vars.back(), idval);
                            params += " " + to_string(value) + " " + to_string(-value);
                        }
                        problem->postKnapsackConstraint(vars, params, false, true, false);
                        params = to_string(0);
                        params2="";
                        for (unsigned int i=0; i < vars.size()-1; i++) {
                            domsize = problem->getDomainInitSize(vars[i]);
                            //params += " " + to_string(domsize);
                            nbval=0;
                            params2="";
                            for (int idval=0; idval < domsize; idval++) {
                                int value = problem->toValue(vars[i], idval);
                                auto it = find(values.begin(), values.end(), value);
                                if (it != values.end()) {
                                    nbval+=1;
                                    params2 += " " + to_string(value) + " " + to_string(-1);
                                }
                            }
                            params+=" "+to_string(nbval)+params2;
                        }
                        domsize=problem->getDomainInitSize(vars.back());
                        params += " " + to_string(domsize);
                        for (int idval=0; idval < domsize; idval++) {
                            int value = problem->toValue(vars.back(), idval);
                            params += " " + to_string(value) + " " + to_string(value);
                        }
                        problem->postKnapsackConstraint(vars, params, false, true, false);
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
                        for (unsigned int i=0; i < vars.size(); i++) {
                            domsize = problem->getDomainInitSize(vars[i]);
                            //params += " " + to_string(domsize);
                            nbval=0;
                            params2="";
                            for (int idval=0; idval < domsize; idval++) {
                                int value = problem->toValue(vars[i], idval);
                                auto it = find(values.begin(), values.end(), value);
                                if (it != values.end()) {
                                    nbval+=1;
                                    params2 += " " + to_string(value) + " " + to_string(-1);
                                }
                            }
                            params+=" "+to_string(nbval)+params2;
                        }
                        problem->postKnapsackConstraint(vars, params, false, true, false);
                        break;
                    case OrderType::LT:
                        params = to_string(-rightcoef+1);
                        params2="";
                        for (unsigned int i=0; i < vars.size(); i++) {
                            domsize = problem->getDomainInitSize(vars[i]);
                            //params += " " + to_string(domsize);
                            nbval=0;
                            params2="";
                            for (int idval=0; idval < domsize; idval++) {
                                int value = problem->toValue(vars[i], idval);
                                auto it = find(values.begin(), values.end(), value);
                                if (it != values.end()) {
                                    nbval+=1;
                                    params2 += " " + to_string(value) + " " + to_string(-1);
                                }
                            }
                            params+=" "+to_string(nbval)+params2;
                        }
                        problem->postKnapsackConstraint(vars, params, false, true, false);
                        break;
                    case OrderType::GE:
                        params = to_string(rightcoef);
                        params2="";
                        for (unsigned int i=0; i < vars.size(); i++) {
                            domsize = problem->getDomainInitSize(vars[i]);
                            //params += " " + to_string(domsize);
                            nbval=0;
                            params2="";
                            for (int idval=0; idval < domsize; idval++) {
                                int value = problem->toValue(vars[i], idval);
                                auto it = find(values.begin(), values.end(), value);
                                if (it != values.end()) {
                                    nbval+=1;
                                    params2 += " " + to_string(value) + " " + to_string(1);
                                }
                            }
                            params+=" "+to_string(nbval)+params2;
                        }
                        problem->postKnapsackConstraint(vars, params, false, true, false);
                        break;
                    case OrderType::GT:
                        params = to_string(rightcoef+1);
                        params2="";
                        for (unsigned int i=0; i < vars.size(); i++) {
                            domsize = problem->getDomainInitSize(vars[i]);
                            //params += " " + to_string(domsize);
                            nbval=0;
                            params2="";
                            for (int idval=0; idval < domsize; idval++) {
                                int value = problem->toValue(vars[i], idval);
                                auto it = find(values.begin(), values.end(), value);
                                if (it != values.end()) {
                                    nbval+=1;
                                    params2 += " " + to_string(value) + " " + to_string(1);
                                }
                            }
                            params+=" "+to_string(nbval)+params2;
                        }
                        problem->postKnapsackConstraint(vars, params, false, true, false);
                        break;
                    case OrderType::NE:
                        params = to_string(rightcoef+1);
                        params2="";
                        for (unsigned int i=0; i < vars.size(); i++) {
                            domsize = problem->getDomainInitSize(vars[i]);
                            //params += " " + to_string(domsize);
                            nbval=0;
                            params2="";
                            for (int idval=0; idval < domsize; idval++) {
                                int value = problem->toValue(vars[i], idval);
                                auto it = find(values.begin(), values.end(), value);
                                if (it != values.end()) {
                                    nbval+=1;
                                    params2 += " " + to_string(value) + " " + to_string(1);
                                }
                            }
                            params+=" "+to_string(nbval)+params2;
                        }
                        problem->postKnapsackConstraint(vars, params, false, true, false);
                        params = to_string(-rightcoef+1);
                        params2="";
                        for (unsigned int i=0; i < vars.size(); i++) {
                            domsize = problem->getDomainInitSize(vars[i]);
                            //params += " " + to_string(domsize);
                            nbval=0;
                            params2="";
                            for (int idval=0; idval < domsize; idval++) {
                                int value = problem->toValue(vars[i], idval);
                                auto it = find(values.begin(), values.end(), value);
                                if (it != values.end()) {
                                    nbval+=1;
                                    params2 += " " + to_string(value) + " " + to_string(-1);
                                }
                            }
                            params+=" "+to_string(nbval)+params2;
                        }
                        problem->postKnapsackConstraint(vars, params, false, true, false);
                        break;
                    case OrderType::EQ:
                        params = to_string(rightcoef);
                        params2="";
                        for (unsigned int i=0; i < vars.size(); i++) {
                            domsize = problem->getDomainInitSize(vars[i]);
                            //params += " " + to_string(domsize);
                            nbval=0;
                            params2="";
                            for (int idval=0; idval < domsize; idval++) {
                                int value = problem->toValue(vars[i], idval);
                                auto it = find(values.begin(), values.end(), value);
                                if (it != values.end()) {
                                    nbval+=1;
                                    params2 += " " + to_string(value) + " " + to_string(1);
                                }
                            }
                            params+=" "+to_string(nbval)+params2;
                        }
                        problem->postKnapsackConstraint(vars, params, false, true, false);
                        params = to_string(-rightcoef);
                        params2="";
                        for (unsigned int i=0; i < vars.size(); i++) {
                            domsize = problem->getDomainInitSize(vars[i]);
                            //params += " " + to_string(domsize);
                            nbval=0;
                            params2="";
                            for (int idval=0; idval < domsize; idval++) {
                                int value = problem->toValue(vars[i], idval);
                                auto it = find(values.begin(), values.end(), value);
                                if (it != values.end()) {
                                    nbval+=1;
                                    params2 += " " + to_string(value) + " " + to_string(-1);
                                }
                            }
                            params+=" "+to_string(nbval)+params2;
                        }
                        problem->postKnapsackConstraint(vars, params, false, true, false);
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
                for (unsigned int i=0; i < vars.size(); i++) {
                    domsize = problem->getDomainInitSize(vars[i]);
                    //params += " " + to_string(domsize);
                    nbval=0;
                    params2="";
                    for (int idval=0; idval < domsize; idval++) {
                        int value = problem->toValue(vars[i], idval);
                        auto it = find(values.begin(), values.end(), value);
                        if (it != values.end()) {
                            nbval+=1;
                            params2 += " " + to_string(value) + " " + to_string(1);
                        }
                    }
                    params+=" "+to_string(nbval)+params2;
                }
                problem->postKnapsackConstraint(vars, params, false, true, false);
                // sum <= cond.max
                params = to_string(-cond.max);
                params2="";
                for (unsigned int i=0; i < vars.size(); i++) {
                    domsize = problem->getDomainInitSize(vars[i]);
                    //params += " " + to_string(domsize);
                    nbval=0;
                    params2="";
                    for (int idval=0; idval < domsize; idval++) {
                        int value = problem->toValue(vars[i], idval);
                        auto it = find(values.begin(), values.end(), value);
                        if (it != values.end()) {
                            nbval+=1;
                            params2 += " " + to_string(value) + " " + to_string(-1);
                        }
                    }
                    params+=" "+to_string(nbval)+params2;
                }
                problem->postKnapsackConstraint(vars, params, false, true, false);
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
            vars.push_back(buildConstraintIntensionVar(trees[i]));
        }
        assert(vars.size() == trees.size());
        buildConstraintCount(vars, values, cond);
    }

    void buildConstraintNValues(vector<int> &vars, vector<int> &except, XCondition &cond) {
        vector<int> Bvar;
        int n=vars.size();
        vector<Value> Diffvalue;
        int domsize;
        vector <Cost> costs;
        bool emptyExcept=false;
        string params;
        int rightcoef;
        if(except.empty())
            emptyExcept= true;
        for (int i = 0; i < n; i++) {
            for (Value currentval : problem->getEnumDomain(vars[i])) {
                if (emptyExcept || currentval != except[0]){
                    auto it = find(Diffvalue.begin(), Diffvalue.end(), currentval);
                    if (it == Diffvalue.end()) {
                        string extravarname = IMPLICIT_VAR_TAG + "Nval" + to_string(problem->numberOfVariables());
                        int extravar = problem->makeEnumeratedVariable(extravarname, 0, 1);
                        mapping[extravarname] = extravar;
                        Bvar.push_back(extravar);
                        Diffvalue.push_back(currentval);
                        costs.clear();
                        for (int b = 0; b < (int)problem->getDomainInitSize(vars[i]); b++) {
                            if (problem->toValue(vars[i], b) == currentval) {
                                costs.push_back(MAX_COST);
                            } else {
                                costs.push_back(MIN_COST);
                            }
                        }
                        for (int b = 0; b < (int)problem->getDomainInitSize(vars[i]); b++) {
                            costs.push_back(MIN_COST);
                        }
                        problem->postBinaryConstraint(Bvar.back(), vars[i], costs);
                    }else{
                        costs.clear();
                        for (int b = 0; b < (int)problem->getDomainInitSize(vars[i]); b++) {
                            if (problem->toValue(vars[i], b) == currentval) {
                                costs.push_back(MAX_COST);
                            } else {
                                costs.push_back(MIN_COST);
                            }
                        }
                        for (int b = 0; b < (int)problem->getDomainInitSize(vars[i]); b++) {
                            costs.push_back(MIN_COST);
                        }
                        problem->postBinaryConstraint(Bvar[distance(Diffvalue.begin(), it)], vars[i], costs);
                    }
                }
            }
        }
        for (unsigned int i=0; i < Diffvalue.size(); i++) {
            Value diffvalue = Diffvalue[i];
            vector<int> tempvars;
            int nbvar=0;
            params="0 ";
            for (int j = 0; j < n; ++j) {
                if (problem->canbe(vars[j], diffvalue)) {
                    nbvar++;
                    tempvars.push_back(vars[j]);
                    params+=" 1 "+to_string(diffvalue)+" 1";
                }
            }
            if (nbvar > 0) {
                tempvars.push_back(Bvar[i]);
                params+=" 1 1 -1";
                problem->postKnapsackConstraint(tempvars, params, false, true, false);
            } // assert(nbvar > 0);
        }
        switch (cond.op) {
        case OrderType::EQ:
            switch (cond.operandType) {
            case OperandType::VARIABLE:
                Bvar.push_back(getMyVar(cond.var));
                params = to_string(0);
                for (unsigned int i=0; i < Bvar.size()-1; i++) {
                    params+=" 1 1 1";
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
                    params+=" 1 1 -1";
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
                    params+=" 1 1 1";
                }
                problem->postKnapsackConstraint(Bvar, params, false, true, false);
                params = to_string(-rightcoef);
                for (unsigned int i=0; i < Bvar.size(); i++) {
                    params+=" 1 1 -1";
                }
                problem->postKnapsackConstraint(Bvar, params, false, true, false);
                break;
            default:
                cerr << "Sorry operandType " << cond.operandType << " not implemented in Nvalues constraint!" << endl;
                throw WrongFileFormat();
            }
            break;
//        case OrderType::GT:
//            params = to_string(2);
//            for (unsigned int i=0; i < Bvar.size(); i++) {
//                params+=" 1 1 1";
//            }
//            problem->postKnapsackConstraint(Bvar, params, false, true, false);
//            break;
        default:
            cerr << "Sorry orderType " << cond.op << " not implemented in Nvalues constraint!" << endl;
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
            vars.push_back(buildConstraintIntensionVar(trees[i]));
        }
        assert(vars.size() == trees.size());
        buildConstraintNValues(vars,except, cond);
    }

    void buildConstraintCardinality(string id, vector<XVariable *> &list, vector<int> values, vector<int> &occurs, bool closed) override {
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
                                params2 += " " + to_string(value) + " " + to_string(-(int)vars.size());
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
                                params2 += " " + to_string(value) + " " + to_string(-(int)vars.size());
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

    void buildConstraintCardinality(string id, vector<XVariable *> &list, vector<int> values, vector<XVariable *> &occurs, bool closed) override {
        vector<int> vars;
        toMyVariables(list,vars);
        vector<int> vars_copy(vars);
        string params;
        string params2;
        int domsize,nbval;
        for (int k = 0; k < (int)values.size(); ++k) {
            vars = vars_copy;
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
                                params2 += " " + to_string(value) + " " + to_string(-(int)vars.size());
                                nbval++;
                            }
                        }
                    }
                }
                params+=" "+to_string(nbval)+params2;
            }
            domsize= problem->getDomainInitSize(vars.back());
            params+=" "+to_string(domsize);
            for (int idval=0; idval < domsize; idval++) {
                int value = problem->toValue(vars.back(), idval);
                params += " " + to_string(value) + " " + to_string(-value);
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
                                params2 += " " + to_string(value) + " " + to_string(-(int)vars.size());
                                nbval++;
                            }
                        }
                    }
                }
                params+=" "+to_string(nbval)+params2;
            }
            domsize= problem->getDomainInitSize(vars.back());
            params+=" "+to_string(domsize);
            for (int idval=0; idval < domsize; idval++) {
                int value = problem->toValue(vars.back(), idval);
                params += " " + to_string(value) + " " + to_string(value);
            }
            problem->postKnapsackConstraint(vars, params, false, true, false);
        }
    }

    void buildConstraintCardinality(string id, vector<XVariable *> &list, vector<int> values, vector<XInterval> &occurs, bool closed) override {
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
                                params2 += " " + to_string(value) + " " + to_string(-(int)vars.size());
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
                                params2 += " " + to_string(value) + " " + to_string(-(int)vars.size());
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
            vars.push_back(buildConstraintIntensionVar(trees[i]));
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
            vars.push_back(buildConstraintIntensionVar(trees[i]));
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
            vars.push_back(buildConstraintIntensionVar(trees[i]));
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
            vars.push_back(buildConstraintIntensionVar(trees[i]));
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
                        vector<Cost> costs((size_t)problem->getDomainInitSize(varvalue) * (size_t)problem->getDomainInitSize(varindex) * (size_t)problem->getDomainInitSize(vars[i]), MIN_COST);
                        for (unsigned int a=0; a < problem->getDomainInitSize(varvalue); a++) {
                            for (unsigned int b=0; b < problem->getDomainInitSize(vars[i]); b++) {
                                switch (op) {
                                case OrderType::IN:
                                case OrderType::EQ:
                                    if (problem->toValue(vars[i], b) != problem->toValue(varvalue, a)) {
                                        costs[(size_t)a * problem->getDomainInitSize(varindex) * problem->getDomainInitSize(vars[i]) + (size_t)problem->toIndex(varindex, (Value)i) * problem->getDomainInitSize(vars[i]) + b] = MAX_COST;
                                    }
                                    break;
                                case OrderType::NE:
                                    if (problem->toValue(vars[i], b) == problem->toValue(varvalue, a)) {
                                        costs[(size_t)a * problem->getDomainInitSize(varindex) * problem->getDomainInitSize(vars[i]) + (size_t)problem->toIndex(varindex, (Value)i) * problem->getDomainInitSize(vars[i]) + b] = MAX_COST;
                                    }
                                    break;
                                case OrderType::LE:
                                    if (problem->toValue(vars[i], b) > problem->toValue(varvalue, a)) {
                                        costs[(size_t)a * problem->getDomainInitSize(varindex) * problem->getDomainInitSize(vars[i]) + (size_t)problem->toIndex(varindex, (Value)i) * problem->getDomainInitSize(vars[i]) + b] = MAX_COST;
                                    }
                                    break;
                                case OrderType::LT:
                                    if (problem->toValue(vars[i], b) >= problem->toValue(varvalue, a)) {
                                        costs[(size_t)a * problem->getDomainInitSize(varindex) * problem->getDomainInitSize(vars[i]) + (size_t)problem->toIndex(varindex, (Value)i) * problem->getDomainInitSize(vars[i]) + b] = MAX_COST;
                                    }
                                    break;
                                case OrderType::GE:
                                    if (problem->toValue(vars[i], b) < problem->toValue(varvalue, a)) {
                                        costs[(size_t)a * problem->getDomainInitSize(varindex) * problem->getDomainInitSize(vars[i]) + (size_t)problem->toIndex(varindex, (Value)i) * problem->getDomainInitSize(vars[i]) + b] = MAX_COST;
                                    }
                                    break;
                                case OrderType::GT:
                                    if (problem->toValue(vars[i], b) <= problem->toValue(varvalue, a)) {
                                        costs[(size_t)a * problem->getDomainInitSize(varindex) * problem->getDomainInitSize(vars[i]) + (size_t)problem->toIndex(varindex, (Value)i) * problem->getDomainInitSize(vars[i]) + b] = MAX_COST;
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
                                    costs[a * problem->getDomainInitSize(varindex) + problem->toIndex(varindex, (Value)i)] = MAX_COST;
                                }
                                break;
                            case OrderType::NE:
                                if ((Value)i == problem->toValue(varvalue, a)) {
                                    costs[a * problem->getDomainInitSize(varindex) + problem->toIndex(varindex, (Value)i)] = MAX_COST;
                                }
                                break;
                            case OrderType::LE:
                                if ((Value)i > problem->toValue(varvalue, a)) {
                                    costs[a * problem->getDomainInitSize(varindex) + problem->toIndex(varindex, (Value)i)] = MAX_COST;
                                }
                                break;
                            case OrderType::LT:
                                if ((Value)i >= problem->toValue(varvalue, a)) {
                                    costs[a * problem->getDomainInitSize(varindex) + problem->toIndex(varindex, (Value)i)] = MAX_COST;
                                }
                                break;
                            case OrderType::GE:
                                if ((Value)i < problem->toValue(varvalue, a)) {
                                    costs[a * problem->getDomainInitSize(varindex) + problem->toIndex(varindex, (Value)i)] = MAX_COST;
                                }
                                break;
                            case OrderType::GT:
                                if ((Value)i <= problem->toValue(varvalue, a)) {
                                    costs[a * problem->getDomainInitSize(varindex) + problem->toIndex(varindex, (Value)i)] = MAX_COST;
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
        vars.push_back(varvalue);
        for (unsigned int a=0; a < problem->getDomainInitSize(varvalue); a++) {
            string params = to_string(1);
            for (unsigned int i=0; i < vars.size()-1; i++) {
                params += " 1 " + to_string(problem->toValue(vars[i], a)) + " 1";
            }
            params += " " + to_string(problem->getDomainInitSize(varvalue) - 1);
            for (unsigned int b=0; b < problem->getDomainInitSize(varvalue); b++) {
                if (a != b) {
                    params += " " + to_string(problem->toValue(varvalue, b))  + " -1";
                }
            }
            problem->postKnapsackConstraint(vars, params, false, true, false);
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

    void buildConstraintElement(string id, vector<int> &list, XVariable *index, int startIndex, XCondition &xc) {
        assert(startIndex == 0);
        int varindex = getMyVar(index->id);
        if (xc.operandType == OperandType::VARIABLE) {
            int varvalue = getMyVar(xc.var);
            assert(varindex != varvalue);
            assert(problem->getDomainInitSize(varindex) == list.size());
            vector<Cost> costs(problem->getDomainInitSize(varvalue) * problem->getDomainInitSize(varindex), MIN_COST);
            for (unsigned int a=0; a < problem->getDomainInitSize(varvalue); a++) {
                for (unsigned int b=0; b < problem->getDomainInitSize(varindex); b++) {
                    switch (xc.op) {
                    case OrderType::LE:
                        if (list[b] > problem->toValue(varvalue, a)) {
                            costs[a * problem->getDomainInitSize(varindex) + b] = MAX_COST;
                        }
                        break;
                    case OrderType::LT:
                        if (list[b] >= problem->toValue(varvalue, a)) {
                            costs[a * problem->getDomainInitSize(varindex) + b] = MAX_COST;
                        }
                        break;
                    case OrderType::GE:
                        if (list[b] < problem->toValue(varvalue, a)) {
                            costs[a * problem->getDomainInitSize(varindex) + b] = MAX_COST;
                        }
                        break;
                    case OrderType::GT:
                        if (list[b] <= problem->toValue(varvalue, a)) {
                            costs[a * problem->getDomainInitSize(varindex) + b] = MAX_COST;
                        }
                        break;
                    case OrderType::IN:
                    case OrderType::EQ:
                        if (list[b] != problem->toValue(varvalue, a)) {
                            costs[a * problem->getDomainInitSize(varindex) + b] = MAX_COST;
                        }
                        break;
                    case OrderType::NE:
                        if (list[b] == problem->toValue(varvalue, a)) {
                            costs[a * problem->getDomainInitSize(varindex) + b] = MAX_COST;
                        }
                        break;
                    default:
                        cerr << "Sorry operator " << xc.op << " not implemented in element constraint!" << endl;
                        throw WrongFileFormat();
                    }
                }
            }
            problem->postBinaryConstraint(varvalue, varindex, costs);
        } else {
            cerr << "Sorry operand type " << xc.operandType << " not implemented in element constraint!" << endl;
            throw WrongFileFormat();
        }
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
                vector<vector<int> > tuples;
                for (unsigned int a=0; a < problem->getDomainInitSize(varvalue); a++) {
                    for (unsigned int b=0; b < problem->getDomainInitSize(vars[i][j]); b++) {
                        if (problem->toValue(varvalue, a) != problem->toValue(vars[i][j], b)) {
                            vector<Value> tuple;
                            tuple.push_back(problem->toValue(varvalue, a));
                            tuple.push_back(problem->toValue(varrowindex, i));
                            tuple.push_back(problem->toValue(varcolindex, j));
                            tuple.push_back(problem->toValue(vars[i][j], b));
                            tuples.push_back(tuple);
                        }
                    }
                }
                buildConstraintExtension(scope, tuples, false, false);
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
                vector<Cost> costs((size_t)problem->getDomainInitSize(varrowindex) * (size_t)problem->getDomainInitSize(varcolindex) * (size_t)problem->getDomainInitSize(vars[i][j]), MIN_COST);
                for (unsigned int b=0; b < problem->getDomainInitSize(vars[i][j]); b++) {
                    if (value != problem->toValue(vars[i][j], b)) {
                        costs[(size_t)i * problem->getDomainInitSize(varcolindex) * problem->getDomainInitSize(vars[i][j]) + (size_t)j * problem->getDomainInitSize(vars[i][j]) + b] = MAX_COST;
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
                assert(problem->toValue(varcolindex, j) == (Value)j);
                vector<Cost> costs((size_t)problem->getDomainInitSize(varrowindex) * (size_t)problem->getDomainInitSize(varcolindex) * (size_t)problem->getDomainInitSize(varvalue), MIN_COST);
                for (unsigned int a=0; a < problem->getDomainInitSize(varvalue); a++) {
                    if (matrix[i][j] != problem->toValue(varvalue, a)) {
                        costs[(size_t)i * problem->getDomainInitSize(varcolindex) * problem->getDomainInitSize(varvalue) + (size_t)j * problem->getDomainInitSize(varvalue) + a] = MAX_COST;
                    }
                }
                problem->postTernaryConstraint(varrowindex, varcolindex, varvalue, costs);
            }
        }
    }

    void buildConstraintChannel(string id, vector<XVariable *> &list1, int startIndex1, vector<XVariable *> &list2, int startIndex2) override {
        assert(startIndex1 == 0);
        assert(startIndex2 == 0);
        vector<int> vars1;
        toMyVariables(list1,vars1);
        vector<int> vars2;
        toMyVariables(list2,vars2);
        if (vars1.size() <= vars2.size()) {
           for (unsigned int j = 0; j < vars2.size(); j++) {
               buildConstraintPrimitive(vars2[j], true, 0, vars1.size()-1);
           }
           for (unsigned int i = 0; i < vars1.size(); i++) {
                buildConstraintPrimitive(vars1[i], true, 0, vars2.size()-1);
                for (unsigned int j = 0; j < vars2.size(); j++) {
                    vector<Cost> costs;
                    for (unsigned int a = 0; a < problem->getDomainInitSize(vars1[i]); a++) {
                        for (unsigned int b = 0; b < problem->getDomainInitSize(vars2[j]); b++) {
                            if ((problem->toValue(vars1[i], a) != (Value)j || problem->toValue(vars2[j], b) == (Value)i) && (vars1.size() < vars2.size() || problem->toValue(vars2[j], b) != (Value)i || problem->toValue(vars1[i], a) == (Value)j)) {
                                costs.push_back(MIN_COST);
                            } else {
                                costs.push_back(MAX_COST);
                            }
                        }
                    }
                    problem->postBinaryConstraint(vars1[i], vars2[j], costs);
                }
            }
        } else {
            cerr << "Sorry wrong list sizes in channel constraint!" << endl;
            throw WrongFileFormat();

        }
    }

    void buildConstraintChannel(string id, vector<XVariable *> &list, int startIndex) override {
        buildConstraintChannel(id, list, startIndex, list, startIndex);
    }

    void buildConstraintChannel(string id, vector<XVariable *> &list, int startIndex, XVariable *value) override {
        assert(startIndex == 0);
        vector<int> vars;
        toMyVariables(list,vars);
        int varvalue =  getMyVar(value);
        buildConstraintPrimitive(varvalue, true, 0, list.size()-1);
        for (unsigned int i = 0; i < vars.size(); i++) {
            vector<Cost> costs(problem->getDomainInitSize(vars[i]) * problem->getDomainInitSize(varvalue), MIN_COST);
            for (unsigned int a = 0; a < problem->getDomainInitSize(vars[i]); a++) {
                for (unsigned int b = 0; b < problem->getDomainInitSize(varvalue); b++) {
                    if ((problem->toValue(vars[i],a) == 0 && problem->toValue(varvalue,b) == (Value)i) || (problem->toValue(vars[i],a) == 1 && problem->toValue(varvalue,b) != (Value)i)) {
                        costs[a * problem->getDomainInitSize(varvalue) + b] = MAX_COST;
                    }
                }
            }
            problem->postBinaryConstraint(vars[i], varvalue, costs);
        }
    }

    void buildConstraintOrdered(string id, vector<XVariable *> &list, OrderType order) override {
        vector<int> lengths(list.size()-1, 0);
        buildConstraintOrdered(id, list, lengths, order);
    }

    void buildConstraintOrdered(string id, vector<XVariable *> &list, vector<int> &lengths, OrderType order) override {
        vector<int> vars;
        toMyVariables(list,vars);
        assert(vars.size() <= lengths.size() + 1);
        for (unsigned int i = 0; i < vars.size()-1; i++) {
            buildConstraintPrimitive(order, vars[i], lengths[i], vars[i+1]);
        }
    }

    void buildConstraintLex(string id, vector<vector<XVariable *>> &lists, OrderType order) override {
        for (unsigned int l = 0; l < lists.size() - 1; l++) {
            assert(lists[l].size() == lists[l+1].size());
            vector<int> vars1;
            toMyVariables(lists[l],vars1);
            vector<int> vars2;
            toMyVariables(lists[l+1],vars2);
            unsigned int n = vars1.size();
            if (order == OrderType::EQ ||order == OrderType::IN) {
                for (unsigned int i = 0; i < n; i++) {
                    buildConstraintPrimitive(OrderType::EQ, vars1[i], 0, vars2[i]);
                }
            } else {
                Value maxdom = 0;
                for (unsigned int i = 0; i < n; i++) {
                    maxdom = max(maxdom, max(problem->getSup(vars1[i]), problem->getSup(vars2[i])) - min(problem->getInf(vars1[i]), problem->getInf(vars2[i])) + 1);
                }
                if (Pow(maxdom, lists[l].size()) < (MAX_COST / maxdom / n)) {
                    vector<int> vars;
                    vector<Long> coefs;
                    Long rightcoef = 0;
                    for (unsigned int i = 0; i < n; i++) {
                        vars.push_back(vars1[i]);
                        coefs.push_back((Long)Pow(maxdom, n-1-i));
                        rightcoef += (Long)Pow(maxdom, n-1-i) * min(problem->getInf(vars1[i]),problem->getInf(vars2[i]));
                    }
                    for (unsigned int i = 0; i < n; i++) {
                        vars.push_back(vars2[i]);
                        coefs.push_back((Long)-Pow(maxdom, n-1-i));
                        rightcoef -= (Long)Pow(maxdom, n-1-i) * min(problem->getInf(vars1[i]),problem->getInf(vars2[i]));
                    }
                    XCondition cond;
                    cond.op = order;
                    cond.val = INT_MAX; //possible integer overflow, do not use it.
                    cond.operandType = OperandType::INTEGER;
                    buildConstraintSum(vars, coefs, cond, rightcoef);
                } else {
                    if (order == OrderType::NE) {
                        vector<int> varseq;
                        vector<Long> coefseq;
                        for (unsigned int i = 0; i < n; i++) {
                            Tree treeeq("eq(" + lists[l][i]->id + "," + lists[l+1][i]->id + ")");
                            int vareq = buildConstraintIntensionVar(&treeeq);
                            varseq.push_back(vareq);
                            coefseq.push_back(1);
                        }
                        XCondition cond;
                        cond.op = OrderType::LE;
                        cond.val = n-1;
                        cond.operandType = OperandType::INTEGER;
                        buildConstraintSum(varseq, coefseq, cond);
                    } else {
                        vector<int> varseq;
                        vector<int> varscond;
                        for (unsigned int i = 0; i < n; i++) {
                            Tree treeeq("eq(" + lists[l][i]->id + "," + lists[l+1][i]->id + ")");
                            int vareq = buildConstraintIntensionVar(&treeeq);
                            varseq.push_back(vareq);
                            Tree treecond(((i<n-1)?lexString[order]:orderString[order]) + "(" + lists[l][i]->id + "," + lists[l+1][i]->id + ")");
                            int varcond = buildConstraintIntensionVar(&treecond);
                            varscond.push_back(varcond);
                        }
                        buildConstraintPrimitive(OrderType::EQ, varscond[0],1);
                        for (unsigned int i = 1; i < n; i++) {
                            vector<int> vars;
                            vector<Long> coefs;
                            int rightcoef = 1;
                            for (unsigned int j = 0; j < i; j++) {
                                vars.push_back(varseq[j]);
                                coefs.push_back(-1);
                                rightcoef -= 1;
                            }
                            vars.push_back(varscond[i]);
                            coefs.push_back(1);
                            XCondition cond;
                            cond.op = OrderType::GE;
                            cond.val = rightcoef;
                            cond.operandType = OperandType::INTEGER;
                            buildConstraintSum(vars, coefs, cond);
                        }
                    }
                }
            }
        }
    }

    void buildConstraintLexMatrix(string id, vector<vector<XVariable *>> &matrix, OrderType order) override {
        buildConstraintLex(id, matrix, order); // row lex
        vector<vector<XVariable *>> transpose;
        for (unsigned int c = 0; c < matrix[0].size(); c++) {
            transpose.push_back(vector<XVariable *>());
            for (unsigned int l = 0; l < matrix.size(); l++) {
                transpose.back().push_back(matrix[l][c]);
            }
        }
        buildConstraintLex(id, transpose, order); // column lex
    }

    void buildConstraintAllEqual(vector<int> &vars) {
        // remove multiple occurrences of the same variable
        for (unsigned int i = 0; i < vars.size(); i++) {
            for (unsigned int j = i + 1; j < vars.size(); j++) {
                if (vars[i] == vars[j]) {
                    vars.erase(vars.begin()+j);
                    j--;
                }
            }
        }
        for (unsigned int i = 0; i < vars.size()-1; i++) {
            buildConstraintPrimitive(OrderType::EQ, vars[i], 0, vars[i+1]);
        }
    }

    void buildConstraintAllEqual(string id, vector<XVariable *> &list) override {
        vector<int> vars;
        toMyVariables(list,vars);
        buildConstraintAllEqual(vars);
    }

    void buildConstraintAllEqual(string id, vector<Tree *> &trees) override {
        vector<int> vars;
        for (unsigned int i = 0; i < trees.size(); i++) {
            vars.push_back(buildConstraintIntensionVar(trees[i]));
        }
        assert(vars.size() == trees.size());
        buildConstraintAllEqual(vars);
    }

    void buildConstraintNotAllEqual(string id, vector<XVariable *> &list) override {
        vector<int> vars;
        toMyVariables(list,vars);
        // remove multiple occurrences of the same variable
        for (unsigned int i = 0; i < vars.size(); i++) {
            for (unsigned int j = i + 1; j < vars.size(); j++) {
                if (vars[i] == vars[j]) {
                    vars.erase(vars.begin()+j);
                    j--;
                }
            }
        }
        assert(vars.size() >= 2);
        set<Value> values;
        for (unsigned int i = 0; i < vars.size(); i++) {
            for (unsigned int a = 0; a < problem->getDomainInitSize(vars[i]); a++) {
                values.insert(problem->toValue(vars[i], a));
            }
        }
        for (Value val: values) {
            string params = to_string(-(int)vars.size() + 1);
            for (unsigned int i = 0; i < vars.size(); i++) {
                params += " 1 " + to_string(val) + " -1";
            }
            problem->postKnapsackConstraint(vars, params, false, true, false);
        }
    }

    // Semantic of subcircuit (successor variables may be not inserted in the circuit by selecting them-self, i.e. self-loops are authorized)
    void buildConstraintCircuit(string id, vector<XVariable *> &list, int startIndex) override {
        vector<int> succ;
        toMyVariables(list,succ);
        assert(startIndex == 0);
        buildConstraintAlldifferent(succ); // keep only one occurrence of each variable
        unsigned int n = succ.size();
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

    void buildConstraintMDD(string id, vector<XVariable*>& list, vector<XTransition>& transitions) override {
        vector<int> vars;
        toMyVariables(list, vars);
        vector<Cost> trCosts(transitions.size(), 0);
        set<string> states;
        map<string, int> stateids;
        map<int, string> statenames;
        map<int, vector<int>> incoming, outgoing;
        map<int, int> statelvl; // from state to level, 0 is the root
        map<int, vector<int>> lvlstates; // from level to vector of states
        map<int, int> stateval; // from state to the value in the
        // corresponding state variable
        int root{-1}, sink{-1};
        int maxid{0};
        for (auto& tr : transitions) {
            states.insert(tr.from);
            states.insert(tr.to);
            if (stateids.count(tr.from) == 0) {
                stateids[tr.from] = maxid;
                statenames[maxid] = tr.from;
                ++maxid;
            }
            if (stateids.count(tr.to) == 0) {
                stateids[tr.to] = maxid;
                statenames[maxid] = tr.to;
                ++maxid;
            }
            outgoing[stateids[tr.from]].push_back(stateids[tr.to]);
            incoming[stateids[tr.to]].push_back(stateids[tr.from]);
        }
        // one more pass over the states to identify the root and sink
        for (auto&& sname : states) {
            int s = stateids[sname];
            if (incoming[s].empty()) {
                if (root >= 0)
                    throw runtime_error("Too many root nodes in MDD constraint " + id);
                root = s;
            }
            if (outgoing[s].empty()) {
                if (sink >= 0)
                    throw runtime_error("Too many sink nodes in MDD constraint " + id);
                sink = s;
            }
        }
        if (root < 0)
            throw runtime_error("No root node in MDD constraint " + id);
        if (sink < 0)
            throw runtime_error("No sink node in MDD constraint " + id);

        // a Breadth-first search to compute levels for each state
        statelvl[root] = 0;
        lvlstates[0].push_back(0);
        stateval[root] = 0;
        vector<int> Q;
        set<int> inq;
        Q.push_back(root);
        inq.insert(root);
        size_t qhead = 0;
        while(qhead != Q.size()) {
            int s = Q[qhead];
            ++qhead;

            int slvl = statelvl[s];
            for(int n : outgoing[s]) {
                if (inq.count(n) == 0) {
                    Q.push_back(n);
                    inq.insert(n);
                    int& nlvl = statelvl[n];
                    nlvl = slvl + 1;
                    stateval[n] = lvlstates[nlvl].size();
                    lvlstates[nlvl].push_back(n);
                }
            }

            // some sanity checking
            for(int p : incoming[s]) {
                if (statelvl[p] != slvl - 1)
                    throw runtime_error("MDD constraint " + id + " state "
                            + statenames[p] + " has long incoming edges");
            }
        }

        // construct the state variables
        if(lvlstates.size() != list.size() + 1)
            throw runtime_error("MDD constraint " + id
                    + " has different number of variables and levels");
        vector<int> statevars(lvlstates.size());
        for (size_t i = 0; i != lvlstates.size(); ++i) {
            string varname = id + "_Q" + to_string(i);
            statevars[i] = problem->makeEnumeratedVariable(varname, 0, lvlstates[i].size());
        }

        // and now the ternary transition constraints
        for (size_t i = 0; i != lvlstates.size()-1; ++i) {
            int tvars[3] = {statevars[i], vars[i], statevars[i + 1]};
            vector<Cost> costs((size_t)problem->getDomainInitSize(tvars[0])
                    * (size_t)problem->getDomainInitSize(tvars[1])
                    * (size_t)problem->getDomainInitSize(tvars[2]),
                    MAX_COST);
            auto idx = [&](int i, int j, int k) {
                return problem->toIndex(tvars[0], i) * problem->getDomainInitSize(tvars[1]) * problem->getDomainInitSize(tvars[2])
                        + problem->toIndex(tvars[1], j) * problem->getDomainInitSize(tvars[2])
                        + problem->toIndex(tvars[2], k);
            };
            for (auto&& tr : transitions) {
                if (statelvl[stateids[tr.from]] != static_cast<int>(i))
                    continue;
                int from = stateids[tr.from];
                int to = stateids[tr.to];
                int val = tr.val;
                costs[idx(stateval[from], val, stateval[to])] = MIN_COST;
            }
            problem->postTernaryConstraint(statevars[i], vars[i], statevars[i + 1], costs);
        }
    }

    void buildConstraintRegular(string id, vector<XVariable *> &list, string start, vector<string> &final, vector<XTransition> &transitions) override {
        vector<int> vars;
        toMyVariables(list, vars);
        vector<Cost> trCosts(transitions.size(), 0);
        set<string> states;
        map<string, int> stateids;
        vector<pair<int, Cost>> initial, accepting;
        vector<int> flat_transitions;
        vector<int*> transition_indices;
        vector<Cost> trcosts;
        int maxid{0};
        for (auto& tr : transitions) {
            states.insert(tr.from);
            states.insert(tr.to);
            if (stateids.count(tr.from) == 0)
                stateids[tr.from] = maxid++;
            if (stateids.count(tr.to) == 0)
                stateids[tr.to] = maxid++;
            flat_transitions.push_back(stateids[tr.from]);
            flat_transitions.push_back(stateids[tr.to]);
            flat_transitions.push_back(tr.val);
            trcosts.push_back(0);
        }

        initial.push_back({stateids[start], 0});
        for (auto &&acc : final)
            accepting.push_back({stateids[acc], 0});
        // populate transition_indices (transitions is now stable, so
        // we can get pointers to it)
        for (size_t i = 0; i != transitions.size(); ++i)
            transition_indices.push_back(&(flat_transitions[i*3]));

        problem->postWRegular(&vars[0], vars.size(), states.size(), initial,
                accepting, &(transition_indices[0]), trCosts);
    }

    void buildConstraintNoOverlap(string id, vector<XVariable *> &origins, vector<int> &lengths, bool zeroIgnored) override {
        assert(origins.size() == lengths.size());
        unsigned int n = origins.size();
        for (unsigned int i = 0; i < n; i++) {
            for (unsigned int j = i +1 ; j < n; j++) {
                if (!zeroIgnored || (lengths[i]>0 && lengths[j]>0)) {
                    Tree tree("or(le(add(" + origins[i]->id + "," + to_string(lengths[i]) + ")," + origins[j]->id + "),le(add(" + origins[j]->id + "," + to_string(lengths[j]) + ")," + origins[i]->id + "))");
                    buildConstraintIntension(id, &tree);
                }
            }
        }
    }

    void buildConstraintNoOverlap(string id, vector<vector<XVariable *>> &origins, vector<vector<int>> &lengths, bool zeroIgnored) override {
        assert(origins.size() == lengths.size());
        if (origins[0].size()==1) {
            unsigned int n = origins.size();
            for (unsigned int i = 0; i < n; i++) {
                for (unsigned int j = i +1 ; j < n; j++) {
                    if (!zeroIgnored || (lengths[i][0]>0 && lengths[j][0]>0)) {
                        Tree tree("or(le(add(" + origins[i][0]->id + "," + to_string(lengths[i][0]) + ")," + origins[j][0]->id + "),le(add(" + origins[j][0]->id + "," + to_string(lengths[j][0]) + ")," + origins[i][0]->id + "))");
                        buildConstraintIntension(id, &tree);
                    }
                }
            }
        } else if (origins[0].size()==2) {
            unsigned int n = origins.size();
            for (unsigned int i = 0; i < n; i++) {
                for (unsigned int j = i +1 ; j < n; j++) {
                    if (!zeroIgnored || (lengths[i][0]>0 && lengths[j][0]>0 && lengths[i][1]>0 && lengths[j][1]>0)) {
                        Tree tree("or(or(le(add(" + origins[i][0]->id + "," + to_string(lengths[i][0]) + ")," + origins[j][0]->id + "),le(add(" + origins[j][0]->id + "," + to_string(lengths[j][0]) + ")," + origins[i][0]->id + ")),or(le(add(" + origins[i][1]->id + "," + to_string(lengths[i][1]) + ")," + origins[j][1]->id + "),le(add(" + origins[j][1]->id + "," + to_string(lengths[j][1]) + ")," + origins[i][1]->id + ")))");
                        buildConstraintIntension(id, &tree);
                    }
                }
            }
        } else if (origins[0].size()==3) {
            unsigned int n = origins.size();
            for (unsigned int i = 0; i < n; i++) {
                for (unsigned int j = i +1 ; j < n; j++) {
                    if (!zeroIgnored || (lengths[i][0]>0 && lengths[j][0]>0 && lengths[i][1]>0 && lengths[j][1]>0 && lengths[i][2]>0 && lengths[j][2]>0)) {
                        Tree tree("or(or(or(le(add(" + origins[i][0]->id + "," + to_string(lengths[i][0]) + ")," + origins[j][0]->id + "),le(add(" + origins[j][0]->id + "," + to_string(lengths[j][0]) + ")," + origins[i][0]->id + ")),or(le(add(" + origins[i][1]->id + "," + to_string(lengths[i][1]) + ")," + origins[j][1]->id + "),le(add(" + origins[j][1]->id + "," + to_string(lengths[j][1]) + ")," + origins[i][1]->id + "))),or(le(add(" + origins[i][2]->id + "," + to_string(lengths[i][2]) + ")," + origins[j][2]->id + "),le(add(" + origins[j][2]->id + "," + to_string(lengths[j][2]) + ")," + origins[i][2]->id + "))))");
                        buildConstraintIntension(id, &tree);
                    }
                }
            }
        } else {
            cerr << "Sorry " << origins[0].size() << " dimension not implemented in NoOverlap constraint!" << endl;
            throw WrongFileFormat();
        }
    }

    void buildConstraintNoOverlap(string id, vector<XVariable *> &origins, vector<XVariable *> &lengths, bool zeroIgnored) override {
        assert(origins.size() == lengths.size());
        unsigned int n = origins.size();
        for (unsigned int i = 0; i < n; i++) {
            for (unsigned int j = i +1 ; j < n; j++) {
                if (!zeroIgnored) {
                    Tree tree("or(le(add(" + origins[i]->id + "," + lengths[i]->id + ")," + origins[j]->id + "),le(add(" + origins[j]->id + "," + lengths[j]->id + ")," + origins[i]->id + "))");
                    buildConstraintIntension(id, &tree);
                } else {
                    Tree tree("or(eq(" + lengths[i]->id + ",0),eq(" + lengths[j]->id + ",0),le(add(" + origins[i]->id + "," + lengths[i]->id + ")," + origins[j]->id + "),le(add(" + origins[j]->id + "," + lengths[j]->id + ")," + origins[i]->id + "))");
                    buildConstraintIntension(id, &tree);
                }
            }
        }
    }

    void buildConstraintNoOverlap(string id, vector<vector<XVariable *>> &origins, vector<vector<XVariable *>> &lengths, bool zeroIgnored) override {
        assert(origins.size() == lengths.size());
        if (origins[0].size()==1) {
            unsigned int n = origins.size();
            for (unsigned int i = 0; i < n; i++) {
                for (unsigned int j = i +1 ; j < n; j++) {
                    if (!zeroIgnored) {
                        Tree tree("or(le(add(" + origins[i][0]->id + "," + lengths[i][0]->id + ")," + origins[j][0]->id + "),le(add(" + origins[j][0]->id + "," + lengths[j][0]->id + ")," + origins[i][0]->id + "))");
                        buildConstraintIntension(id, &tree);
                    } else {
                        Tree tree("or(eq(" + lengths[i][0]->id + ",0),eq(" + lengths[j][0]->id + ",0),le(add(" + origins[i][0]->id + "," + lengths[i][0]->id + ")," + origins[j][0]->id + "),le(add(" + origins[j][0]->id + "," + lengths[j][0]->id + ")," + origins[i][0]->id + "))");
                        buildConstraintIntension(id, &tree);
                    }
                }
            }
        } else if (origins[0].size()==2) {
            unsigned int n = origins.size();
            for (unsigned int i = 0; i < n; i++) {
                for (unsigned int j = i +1 ; j < n; j++) {
                    if (!zeroIgnored) {
                        Tree tree("or(or(le(add(" + origins[i][0]->id + "," + lengths[i][0]->id + ")," + origins[j][0]->id + "),le(add(" + origins[j][0]->id + "," + lengths[j][0]->id + ")," + origins[i][0]->id + ")),or(le(add(" + origins[i][1]->id + "," + lengths[i][1]->id + ")," + origins[j][1]->id + "),le(add(" + origins[j][1]->id + "," + lengths[j][1]->id + ")," + origins[i][1]->id + ")))");
                        buildConstraintIntension(id, &tree);
                    } else {
                        Tree tree("or(eq(" + lengths[i][0]->id + ",0),eq(" + lengths[j][0]->id + ",0),eq(" + lengths[i][1]->id + ",0),eq(" + lengths[j][1]->id + ",0),or(le(add(" + origins[i][0]->id + "," + lengths[i][0]->id + ")," + origins[j][0]->id + "),le(add(" + origins[j][0]->id + "," + lengths[j][0]->id + ")," + origins[i][0]->id + ")),or(le(add(" + origins[i][1]->id + "," + lengths[i][1]->id + ")," + origins[j][1]->id + "),le(add(" + origins[j][1]->id + "," + lengths[j][1]->id + ")," + origins[i][1]->id + ")))");
                        buildConstraintIntension(id, &tree);
                    }
                }
            }
        } else if (origins[0].size()==3) {
            unsigned int n = origins.size();
            for (unsigned int i = 0; i < n; i++) {
                for (unsigned int j = i +1 ; j < n; j++) {
                    if (!zeroIgnored) {
                        Tree tree("or(or(or(le(add(" + origins[i][0]->id + "," + lengths[i][0]->id + ")," + origins[j][0]->id + "),le(add(" + origins[j][0]->id + "," + lengths[j][0]->id + ")," + origins[i][0]->id + ")),or(le(add(" + origins[i][1]->id + "," + lengths[i][1]->id + ")," + origins[j][1]->id + "),le(add(" + origins[j][1]->id + "," + lengths[j][1]->id + ")," + origins[i][1]->id + "))),or(le(add(" + origins[i][2]->id + "," + lengths[i][2]->id + ")," + origins[j][2]->id + "),le(add(" + origins[j][2]->id + "," + lengths[j][2]->id + ")," + origins[i][2]->id + "))))");
                        buildConstraintIntension(id, &tree);
                    } else {
                        Tree tree("or(eq(" + lengths[i][0]->id + ",0),eq(" + lengths[j][0]->id + ",0),eq(" + lengths[i][1]->id + ",0),eq(" + lengths[j][1]->id + ",0),eq(" + lengths[i][2]->id + ",0),eq(" + lengths[j][2]->id + ",0),or(or(le(add(" + origins[i][0]->id + "," + lengths[i][0]->id + ")," + origins[j][0]->id + "),le(add(" + origins[j][0]->id + "," + lengths[j][0]->id + ")," + origins[i][0]->id + ")),or(le(add(" + origins[i][1]->id + "," + lengths[i][1]->id + ")," + origins[j][1]->id + "),le(add(" + origins[j][1]->id + "," + lengths[j][1]->id + ")," + origins[i][1]->id + "))),or(le(add(" + origins[i][2]->id + "," + lengths[i][2]->id + ")," + origins[j][2]->id + "),le(add(" + origins[j][2]->id + "," + lengths[j][2]->id + ")," + origins[i][2]->id + "))))");
                        buildConstraintIntension(id, &tree);
                    }
                }
            }
        } else {
            cerr << "Sorry " << origins[0].size() << " dimension not implemented in NoOverlap constraint!" << endl;
            throw WrongFileFormat();
        }
    }

    void buildConstraintCumulative(string id, vector<XVariable *> &origins, vector<int> &lengths, vector<int> &heights, XCondition &cond) override {
        vector<int> vars;
        toMyVariables(origins, vars);
        int mininf = INT_MAX;
        int maxsup = -INT_MAX;
        for (unsigned int i=0; i<vars.size(); i++) {
            if (problem->getInf(vars[i]) < mininf) {
                mininf = problem->getInf(vars[i]);
            }
            if (problem->getSup(vars[i]) > maxsup) {
                maxsup = problem->getSup(vars[i]);
            }
        }
        for (Value time = mininf; time <= maxsup; time++) {
            string paramspos;
            string paramsneg;
            vector<int> scope;
            int nbvar = 0;
            for (unsigned int i = 0; i < vars.size(); i++) {
                string valparamspos;
                string valparamsneg;
                int nbval = 0;
                for (Value value : problem->getEnumDomain(vars[i])) {
                    if (time >= value && time < value + lengths[i]) {
                        valparamspos += " " + to_string(value) + " " + to_string(heights[i]);
                        valparamsneg += " " + to_string(value) + " " + to_string(-heights[i]);
                        nbval++;
                    }
                }
                if (nbval > 0) {
                    paramspos += " " + to_string(nbval) + valparamspos;
                    paramsneg += " " + to_string(nbval) + valparamsneg;
                    scope.push_back(vars[i]);
                    nbvar++;
                }
            }
            if (nbvar > 0) {
                switch (cond.operandType) {
                case OperandType::INTEGER:
                    switch (cond.op) {
                    case OrderType::LE:
                        problem->postKnapsackConstraint(scope, to_string(-cond.val) + paramsneg, false, true, false);
                        break;
                    case OrderType::LT:
                        problem->postKnapsackConstraint(scope, to_string(-cond.val + 1) + paramsneg, false, true, false);
                        break;
                    case OrderType::GE:
                        problem->postKnapsackConstraint(scope, to_string(cond.val) + paramspos, false, true, false);
                        break;
                    case OrderType::GT:
                        problem->postKnapsackConstraint(scope, to_string(cond.val + 1) + paramspos, false, true, false);
                        break;
                    case OrderType::IN:
                    case OrderType::EQ:
                        problem->postKnapsackConstraint(scope, to_string(cond.val) + paramspos, false, true, false);
                        problem->postKnapsackConstraint(scope, to_string(-cond.val) + paramsneg, false, true, false);
                        break;
                    case OrderType::NE:
                    default:
                        cerr << "Sorry operator " << cond.op << " not implemented in cumulative constraint!" << endl;
                        throw WrongFileFormat();
                    }
                    break;
                case OperandType::INTERVAL:
                    switch (cond.op) {
                    case OrderType::IN:
                        problem->postKnapsackConstraint(scope, to_string(cond.min) + paramspos, false, true, false);
                        problem->postKnapsackConstraint(scope, to_string(-cond.max) + paramsneg, false, true, false);
                        break;
                    default:
                        cerr << "Sorry operator " << cond.op << " not implemented in cumulative constraint with interval!" << endl;
                        throw WrongFileFormat();
                    }
                    break;
                case OperandType::VARIABLE:
                default:
                    cerr << "Sorry operand type VARIABLE not implemented in cumulative constraint!" << endl;
                    throw WrongFileFormat();
                }
            }
        }
    }

    void buildConstraintBinPacking(string id, vector<XVariable *> &list, vector<int> &sizes, XCondition &cond) override {
        vector<int> vars;
        toMyVariables(list,vars);
        set<Value> svalues;
        for (unsigned int i = 0; i < vars.size(); i++) {
            for (Value currentval : problem->getEnumDomain(vars[i])) {
                auto it = svalues.find(currentval);
                if (it == svalues.end()) {
                    svalues.insert(currentval);
                }
            }
        }
        for (Value value : svalues) { // for each bin
            vector<Tree> mytrees; // keep in memory each Tree object
            for (unsigned int i=0; i<list.size(); i++) {
                Tree tree("eq(" + list[i]->id + "," + to_string(value) + ")" );
                mytrees.push_back(tree);
            }
            vector<Tree *> trees(list.size(), NULL);
            for (unsigned int i=0; i<list.size(); i++) {
                trees[i] = &mytrees[i]; // give access to each Tree object using a pointer // weird bug occurring with the following code trees.push_back(&mytrees.back())
            }
            assert(list.size() == trees.size());
            buildConstraintSum(id, trees, sizes, cond);
        }
    }

    void buildConstraintBinPacking(string id, vector<XVariable *> &list, vector<int> &sizes, vector<int> &capacities, bool load) override {
        vector<int> vars;
        toMyVariables(list,vars);
        set<Value> svalues;
        for (unsigned int i = 0; i < vars.size(); i++) {
            for (Value currentval : problem->getEnumDomain(vars[i])) {
                auto it = svalues.find(currentval);
                if (it == svalues.end()) {
                    svalues.insert(currentval);
                }
            }
        }
        vector<Value> values;
        for (Value currentval : svalues) {
            values.push_back(currentval);
        }
        for (unsigned int v = 0; v < values.size() ; v++) { // for each bin
            Value value = values[v];
            vector<Tree> mytrees; // keep in memory each Tree object
            for (unsigned int i=0; i<list.size(); i++) {
                Tree tree("eq(" + list[i]->id + "," + to_string(value) + ")" );
                mytrees.push_back(tree);
            }
            vector<Tree *> trees(list.size(), NULL);
            for (unsigned int i=0; i<list.size(); i++) {
                trees[i] = &mytrees[i]; // give access to each Tree object using a pointer // weird bug occurring with the following code trees.push_back(&mytrees.back())
            }
            assert(list.size() == trees.size());
            XCondition cond;
            if (load) {
                cond.op = OrderType::EQ;
            } else {
                cond.op = OrderType::LE;
            }
            cond.val = capacities[v];
            cond.operandType = OperandType::INTEGER;
            buildConstraintSum(id, trees, sizes, cond);
        }
    }

    void buildConstraintBinPacking(string id, vector<XVariable *> &list, vector<int> &sizes, vector<XVariable*> &capacities, bool load) override {
        vector<int> vars;
        toMyVariables(list,vars);
        set<Value> svalues;
        for (unsigned int i = 0; i < vars.size(); i++) {
            for (Value currentval : problem->getEnumDomain(vars[i])) {
                auto it = svalues.find(currentval);
                if (it == svalues.end()) {
                    svalues.insert(currentval);
                }
            }
        }
        vector<Value> values;
        for (Value currentval : svalues) {
            values.push_back(currentval);
        }
        for (unsigned int v = 0; v < values.size() ; v++) { // for each bin
            Value value = values[v];
            vector<Tree> mytrees; // keep in memory each Tree object
            for (unsigned int i=0; i<list.size(); i++) {
                Tree tree("eq(" + list[i]->id + "," + to_string(value) + ")" );
                mytrees.push_back(tree);
            }
            vector<Tree *> trees(list.size(), NULL);
            for (unsigned int i=0; i<list.size(); i++) {
                trees[i] = &mytrees[i]; // give access to each Tree object using a pointer // weird bug occurring with the following code trees.push_back(&mytrees.back())
            }
            assert(list.size() == trees.size());
            XCondition cond;
            if (load) {
                cond.op = OrderType::EQ;
            } else {
                cond.op = OrderType::LE;
            }
            cond.var = capacities[v]->id;
            cond.operandType = OperandType::VARIABLE;
            buildConstraintSum(id, trees, sizes, cond);
        }
    }

    void buildConstraintBinPacking(string id, vector<XVariable *> &list, vector<int> &sizes, vector<XCondition> &conditions, int startindex) override {
        assert(startindex == 0);
        vector<int> vars;
        toMyVariables(list,vars);
        set<Value> svalues;
        for (unsigned int i = 0; i < vars.size(); i++) {
            for (Value currentval : problem->getEnumDomain(vars[i])) {
                auto it = svalues.find(currentval);
                if (it == svalues.end()) {
                    svalues.insert(currentval);
                }
            }
        }
        vector<Value> values;
        for (Value currentval : svalues) {
            values.push_back(currentval);
        }
        for (unsigned int v = 0; v < values.size() ; v++) { // for each bin
            Value value = values[v];
            vector<Tree> mytrees; // keep in memory each Tree object
            for (unsigned int i=0; i<list.size(); i++) {
                Tree tree("eq(" + list[i]->id + "," + to_string(value) + ")" );
                mytrees.push_back(tree);
            }
            vector<Tree *> trees(list.size(), NULL);
            for (unsigned int i=0; i<list.size(); i++) {
                trees[i] = &mytrees[i]; // give access to each Tree object using a pointer // weird bug occurring with the following code trees.push_back(&mytrees.back())
            }
            assert(list.size() == trees.size());
            buildConstraintSum(id, trees, sizes, conditions[v]);
        }
    }

    void buildConstraintKnapsack(string id, vector<XVariable *> &list, vector<int> &weights, vector<int> &profits, XCondition weightsCondition, XCondition &profitCondition) override {
        buildConstraintSum(id, list, weights, weightsCondition);
        buildConstraintSum(id, list, profits, profitCondition);
    }

    void buildConstraintPrecedence(string id, vector<XVariable *> &list, vector<int> values, bool covered) override {
        vector<int> vars;
        toMyVariables(list,vars);
        vector<Cost> costs(problem->getDomainInitSize(vars[0]), MIN_COST);
        for (unsigned int v = 1; v < values.size(); v++) {
            costs[problem->toIndex(vars[0], values[v])] = MAX_COST;
        }
        problem->postUnaryConstraint(vars[0], costs);  // forbid all values for the first variables except values[0]
        for (int v = 0; v < (int)values.size() - 1; v++) {
            for (unsigned int i = 1; i < vars.size(); i++) {
                vector<int> scope;
                string params = "0";
                for (unsigned int j = 0; j < i; j++) {
                    scope.push_back(vars[j]);
                    params += " 1 " + to_string(values[v]) + " 1";
                }
                scope.push_back(vars[i]);
                params += " 1 " + to_string(values[v+1]) + " -1";
                problem->postKnapsackConstraint(scope, params, false, true, false); // sum( x[j] = values[v] , j < i ) >= (x[i] == values[v+1])
            }
        }
        if (covered) {
            for (unsigned int v = 0; v < values.size(); v++) {
                string params = "1 ";
                for (unsigned int i = 0; i < vars.size(); i++) {
                    params += " 1 " + to_string(values[v]) + " 1";
                }
                problem->postKnapsackConstraint(vars, params, false, true, false); // each value must be selected at least once
            }
        }
    }

    void buildConstraintPrecedence(string id, vector<XVariable *> &list, bool covered) override {
        vector<int> vars;
        toMyVariables(list,vars);
        set<Value> svalues;
        for (unsigned int i = 0; i < vars.size(); i++) {
            for (Value currentval : problem->getEnumDomain(vars[i])) {
                auto it = svalues.find(currentval);
                if (it == svalues.end()) {
                    svalues.insert(currentval);
                }
            }
        }
        vector<Value> values;
//        values.assign(svalues.begin(), svalues.end());
        for (Value currentval : svalues) {
            values.push_back(currentval);
        }

        buildConstraintPrecedence(id, list, values, covered);
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
        ToulBar2::xmlcop = true;
        buildUnaryCostFunction(1, x);
    }

    void buildObjectiveMaximizeVariable(XVariable *x) override {
        ToulBar2::xmlcop = true;
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
            vector<Cost> costs((size_t)problem->getDomainInitSize(vars[0]) * (size_t)problem->getDomainInitSize(vars[1]) * (size_t)problem->getDomainInitSize(vars[2]), MAX_COST);
            for(unsigned int i=0; i<tuples.size(); i++) {
//                costs[(size_t)problem->toIndex(vars[0], tuples[i][0]) * problem->getDomainInitSize(vars[1]) * problem->getDomainInitSize(vars[2]) + (size_t)problem->toIndex(vars[1], tuples[i][1]) * problem->getDomainInitSize(vars[2]) + problem->toIndex(vars[2], tuples[i][2])] = tcosts[i] - negcost;
                costs[(size_t)problem->toIndex(vars[2], tuples[i][2]) * problem->getDomainInitSize(vars[0]) * problem->getDomainInitSize(vars[1]) + (size_t)problem->toIndex(vars[0], tuples[i][0]) * problem->getDomainInitSize(vars[1]) + problem->toIndex(vars[1], tuples[i][1])] = tcosts[i] - negcost;
            }
//            problem->postTernaryConstraint(vars[0], vars[1], vars[2], costs);
            problem->postTernaryConstraint(vars[2], vars[0], vars[1], costs);
        } else {
            int ctridx = problem->postNaryConstraintBegin(vars, MAX_COST, tuples.size());
            for(unsigned int i=0; i<tuples.size(); i++) {
                problem->postNaryConstraintTuple(ctridx, tuples[i], tcosts[i] - negcost);
            }
            problem->postNaryConstraintEnd(ctridx);
        }
    }

    void buildObjective(Cost sign, ExpressionObjective type, vector<int> &vars, vector<int> &coefs) {
        assert(vars.size() == coefs.size());
//        Cost negcost = MIN_COST;
//        vector<WeightedVarValPair> weightFunction;
        vector<int> objvars;
        string varargmaxname;
        int varargmax;
        string varmaxname;
        int varmax;
        string varnvaluename;
        int varnvalue;
        vector<int> except;
        XCondition cond;
        bool themax = true;
        int maxinf = INT_MAX;
        int maxsup = -INT_MAX;
        switch (type) {
//        case ExpressionObjective::PRODUCT_O:
//            if (list.size() == 2) {
//                for (unsigned int a=0; a < problem->getDomainInitSize(vars[0]); a++) {
//                    for (unsigned int b=0; b < problem->getDomainInitSize(vars[1]); b++) {
//                        Cost cost = (Cost)coefs[0] * problem->toValue(vars[0], a) * coefs[1] * problem->toValue(vars[1], b) * sign;
//                        if (cost < negcost) {
//                            negcost = cost;
//                        }
//                    }
//                }
//                vector<Cost> costs;
//                for (unsigned int a=0; a < problem->getDomainInitSize(vars[0]); a++) {
//                    for (unsigned int b=0; b < problem->getDomainInitSize(vars[1]); b++) {
//                        costs.push_back((Cost)coefs[0] * problem->toValue(vars[0], a) * coefs[1] * problem->toValue(vars[1], b) * sign - negcost);
//                    }
//                }
//                problem->postBinaryConstraint(vars[0], vars[1], costs);
//                break;
//            } else if (list.size() >= 3) {
//                cerr << "Sorry multiplicative objective for " << list.size() << " variables not implemented!" << endl;
//                throw WrongFileFormat();
//            }
        case ExpressionObjective::EXPRESSION_O:
            assert(vars.size() == 1);
        case ExpressionObjective::SUM_O:
            for (unsigned int i=0; i<vars.size(); i++) {
                buildUnaryCostFunction(sign * coefs[i], vars[i]);
            }
            break;
        case ExpressionObjective::MINIMUM_O:
            themax = false;
        case ExpressionObjective::MAXIMUM_O:
            for (unsigned int i=0; i<vars.size(); i++) {
                int var = vars[i];
                if (coefs[i] != 1) {
                    string prodname = IMPLICIT_VAR_TAG + "prod" + to_string(problem->numberOfVariables());
                    int varprod = problem->makeEnumeratedVariable(prodname, min(problem->getInf(var) * coefs[i], problem->getSup(var) * coefs[i]), max(problem->getInf(var) * coefs[i], problem->getSup(var) * coefs[i]));
                    mapping[prodname] = varprod;
                    buildConstraintPrimitiveMult(OrderType::EQ, var, coefs[i], varprod);
                    objvars.push_back(varprod);
                } else {
                    objvars.push_back(var);
                }
            }
            assert(objvars.size() == coefs.size());
            for (unsigned int i=0; i<objvars.size(); i++) {
                if (problem->getInf(objvars[i]) < maxinf) {
                    maxinf = problem->getInf(objvars[i]);
                }
                if (problem->getSup(objvars[i]) > maxsup) {
                    maxsup = problem->getSup(objvars[i]);
                }
            }
            varmaxname = IMPLICIT_VAR_TAG + ((themax)?"max":"min") + to_string(problem->numberOfVariables());
            varmax = problem->makeEnumeratedVariable(varmaxname, maxinf, maxsup);
            mapping[varmaxname] = varmax;
            if ((sign==UNIT_COST && themax) || (sign==-UNIT_COST && !themax)) {
                for (unsigned int i=0; i<objvars.size(); i++) {
                    buildConstraintPrimitive((themax)?OrderType::LE:OrderType::GE, objvars[i], 0, varmax);
                }
            } else {
                varargmaxname = IMPLICIT_VAR_TAG + ((themax)?"argmax":"argmin") + to_string(problem->numberOfVariables());
                varargmax = problem->makeEnumeratedVariable(varargmaxname, 0, objvars.size()-1);
                mapping[varargmaxname] = varargmax;
                cond.operandType = OperandType::VARIABLE;
                cond.op = OrderType::EQ;
                cond.var = varmaxname;
                buildConstraintMinMax(themax, objvars, varargmax, cond);
            }
            buildUnaryCostFunction(sign, varmax);
            break;
//        case ExpressionObjective::MINIMUM_O:
//            sign *= -UNIT_COST;
//        case ExpressionObjective::MAXIMUM_O:
//            for (unsigned int i=0; i<list.size(); i++) {
//                for (unsigned int a=0; a < problem->getDomainInitSize(vars[i]); a++) {
//                    Cost cost = (Cost)problem->toValue(vars[i], a) * coefs[i] * sign;
//                    if (cost < negcost) {
//                        negcost = cost;
//                    }
//                }
//            }
//            for (unsigned int i=0; i<list.size(); i++) {
//                for (unsigned int a=0; a < problem->getDomainInitSize(vars[i]); a++) {
//                    WeightedVarValPair elt;
//                    elt.varIndex = vars[i];
//                    elt.val = problem->toValue(vars[i], a);
//                    elt.weight = (Cost)elt.val * coefs[i] * sign - negcost;
//                    weightFunction.push_back(elt);
//                }
//            }
//            if (negcost < MIN_COST) {
//                problem->decreaseLb(-negcost);
//            }
//            problem->postMaxWeight(vars.data(), vars.size(), "val", "DAG", MIN_COST, weightFunction);
//            break;
        case ExpressionObjective::NVALUES_O:
            varnvaluename = IMPLICIT_VAR_TAG + "nvalue" + to_string(problem->numberOfVariables());
            varnvalue = problem->makeEnumeratedVariable(varnvaluename, 0, vars.size());
            mapping[varnvaluename] = varnvalue;
            cond.operandType = OperandType::VARIABLE;
            cond.op = OrderType::EQ;
            cond.var = varnvaluename;
            buildConstraintNValues(vars, except, cond);
            buildUnaryCostFunction(sign, varnvalue);
            break;
        case ExpressionObjective::LEX_O:
        default:
            cerr << "Sorry objective type " << type << " not implemented!" << endl;
            throw WrongFileFormat();
        }
    }

    void buildObjective(Cost sign, ExpressionObjective type, vector<XVariable *> &list, vector<int> &coefs) {
        assert(list.size() == coefs.size());
        vector<int> vars;
        toMyVariables(list,vars);
        buildObjective(sign, type, vars, coefs);
    }

    void buildObjectiveMinimize(ExpressionObjective type, vector<XVariable *> &list, vector<int> &coefs) override {
        ToulBar2::xmlcop = true;
        buildObjective(UNIT_COST, type, list, coefs);
    }

    void buildObjectiveMaximize(ExpressionObjective type, vector<XVariable *> &list, vector<int> &coefs) override {
        ToulBar2::xmlcop = true;
        ToulBar2::costMultiplier *= -1.0;
        buildObjective(-UNIT_COST, type, list, coefs);
    }

//    void buildObjectiveMinimize(ExpressionObjective type, vector<XVariable *> &list, vector<XVariable*> &coefs) override {
//        ToulBar2::xmlcop = true;
//        vector<int> vars;
//        toMyVariables(list,vars);
//        vector<int> varscoefs;
//        toMyVariables(coefs,varscoefs);
//        assert(vars.size() == varscoefs.size());
//        vector<int> mycoefs(list.size(), 1);
//        vector<int> varsprod;
//        for (unsigned int i = 0; i < list.size(); i++) {
//            string prodname = IMPLICIT_VAR_TAG + "prod" + to_string(problem->numberOfVariables());
//            Value prodinf = min(min(problem->getInf(vars[i]) * problem->getInf(varscoefs[i]), problem->getInf(vars[i]) * problem->getSup(varscoefs[i])), min(problem->getSup(vars[i]) * problem->getInf(varscoefs[i]), problem->getSup(vars[i]) * problem->getSup(varscoefs[i])));
//            Value prodsup = max(max(problem->getInf(vars[i]) * problem->getInf(varscoefs[i]), problem->getInf(vars[i]) * problem->getSup(varscoefs[i])), max(problem->getSup(vars[i]) * problem->getInf(varscoefs[i]), problem->getSup(vars[i]) * problem->getSup(varscoefs[i])));
//            assert(prodinf <= prodsup);
//            int varprod = problem->makeEnumeratedVariable(prodname, prodinf, prodsup);
//            mapping[prodname] = varprod;
//            buildConstraintMult(vars[i], varscoefs[i], varprod);
//            varsprod.push_back(varprod);
//        }
//        assert(varsprod.size() == list.size());
//        buildObjective(UNIT_COST, type, varsprod, mycoefs);
//    }
    void buildObjectiveMinimize(ExpressionObjective type, vector<XVariable *> &list, vector<XVariable*> &coefs) override {
        assert(list.size() == coefs.size());
        vector<Tree> mytrees; // keep in memory each Tree object
        for (unsigned int i=0; i<list.size(); i++) {
            Tree tree("mul(" + coefs[i]->id + "," + list[i]->id + ")" );
            mytrees.push_back(tree);
        }
        vector<Tree *> trees(list.size(), NULL);
        for (unsigned int i=0; i<list.size(); i++) {
            trees[i] = &mytrees[i]; // give access to each Tree object using a pointer // weird bug occuring with the following code trees.push_back(&mytrees.back())
        }
        assert(list.size() == trees.size());
        buildObjectiveMinimize(type, trees);
    }

    void buildObjectiveMaximize(ExpressionObjective type, vector<XVariable *> &list, vector<XVariable*> &coefs) override {
        assert(list.size() == coefs.size());
        vector<Tree> mytrees;
        for (unsigned int i=0; i<list.size(); i++) {
            Tree tree("mul(" + coefs[i]->id + "," + list[i]->id + ")" );
            mytrees.push_back(tree);
        }
        vector<Tree *> trees(list.size(), NULL);
        for (unsigned int i=0; i<list.size(); i++) {
            trees[i] = &mytrees[i];
        }
        assert(list.size() == trees.size());
        buildObjectiveMaximize(type, trees);
    }

    void buildObjectiveMinimize(ExpressionObjective type, vector<XVariable *> &list) override {
        ToulBar2::xmlcop = true;
        vector<int> coefs(list.size(), 1);
        buildObjectiveMinimize(type, list, coefs);
    }

    void buildObjectiveMaximize(ExpressionObjective type, vector<XVariable *> &list) override {
        ToulBar2::xmlcop = true;
        vector<int> coefs(list.size(), 1);
        buildObjectiveMaximize(type, list, coefs);
    }

    void buildObjective(Cost sign, ExpressionObjective type, vector<Tree *> &trees, vector<int> &coefs) {
        assert(trees.size() == coefs.size());
        vector<int> vars;
        string varargmaxname;
        int varargmax;
        string varmaxname;
        int varmax;
        string varnvaluename;
        int varnvalue;
        vector<int> except;
        XCondition cond;
        bool themax = true;
        int maxinf = INT_MAX;
        int maxsup = -INT_MAX;
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
                int var = buildConstraintIntensionVar(trees[i]);
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
            for (unsigned int i=0; i<vars.size(); i++) {
                if (problem->getInf(vars[i]) < maxinf) {
                    maxinf = problem->getInf(vars[i]);
                }
                if (problem->getSup(vars[i]) > maxsup) {
                    maxsup = problem->getSup(vars[i]);
                }
            }
            varmaxname = IMPLICIT_VAR_TAG + ((themax)?"max":"min") + to_string(problem->numberOfVariables());
            varmax = problem->makeEnumeratedVariable(varmaxname, maxinf, maxsup);
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
        case ExpressionObjective::NVALUES_O:
            varnvaluename = IMPLICIT_VAR_TAG + "nvalue" + to_string(problem->numberOfVariables());
            varnvalue = problem->makeEnumeratedVariable(varnvaluename, 0, trees.size());
            mapping[varnvaluename] = varnvalue;
            cond.operandType = OperandType::VARIABLE;
            cond.op = OrderType::EQ;
            cond.var = varnvaluename;
            buildConstraintNValues("NValue", trees, cond);
            buildUnaryCostFunction(sign, varnvalue);
            break;
        case ExpressionObjective::PRODUCT_O:
        case ExpressionObjective::LEX_O:
        default:
            cerr << "Sorry objective type " << type << " on trees not implemented!" << endl;
            throw WrongFileFormat();
        }
    }

    void buildObjectiveMinimize(ExpressionObjective type, vector<Tree *> &trees, vector<int> &coefs) override {
        ToulBar2::xmlcop = true;
        buildObjective(UNIT_COST, type, trees, coefs);
    }

    void buildObjectiveMaximize(ExpressionObjective type, vector<Tree *> &trees, vector<int> &coefs) override {
        ToulBar2::xmlcop = true;
        ToulBar2::costMultiplier *= -1.0;
        buildObjective(-UNIT_COST, type, trees, coefs);
    }

    void buildObjectiveMinimize(ExpressionObjective type, vector<Tree *> &trees, vector<XVariable*> &coefs) override {
        assert(trees.size() == coefs.size());
        vector<Tree> mytrees; // keep in memory each Tree object
        for (unsigned int i=0; i<trees.size(); i++) {
            Tree tree("mul(" + coefs[i]->id + "," + trees[i]->toString() + ")" );
            mytrees.push_back(tree);
        }
        vector<Tree *> ptrtrees(trees.size(), NULL);
        for (unsigned int i=0; i<trees.size(); i++) {
            ptrtrees[i] = &mytrees[i]; // give access to each Tree object using a pointer // weird bug occuring with the following code trees.push_back(&mytrees.back())
        }
        assert(ptrtrees.size() == trees.size());
        buildObjectiveMinimize(type, ptrtrees);
    }

    void buildObjectiveMaximize(ExpressionObjective type, vector<Tree *> &trees, vector<XVariable *> &coefs) override {
        assert(trees.size() == coefs.size());
        vector<Tree> mytrees; // keep in memory each Tree object
        for (unsigned int i=0; i<trees.size(); i++) {
            Tree tree("mul(" + coefs[i]->id + "," + trees[i]->toString() + ")" );
            mytrees.push_back(tree);
        }
        vector<Tree *> ptrtrees(trees.size(), NULL);
        for (unsigned int i=0; i<trees.size(); i++) {
            ptrtrees[i] = &mytrees[i]; // give access to each Tree object using a pointer // weird bug occuring with the following code trees.push_back(&mytrees.back())
        }
        assert(ptrtrees.size() == trees.size());
        buildObjectiveMaximize(type, ptrtrees);
    }

    void buildObjectiveMinimize(ExpressionObjective type, vector<Tree *> &trees) override {
        ToulBar2::xmlcop = true;
        vector<int> coefs(trees.size(), 1);
        buildObjectiveMinimize(type, trees, coefs);
    }

    void buildObjectiveMaximize(ExpressionObjective type, vector<Tree *> &trees) override {
        ToulBar2::xmlcop = true;
        vector<int> coefs(trees.size(), 1);
        buildObjectiveMaximize(type, trees, coefs);
    }

    void buildObjectiveMinimizeExpression(string expr) override {
        ToulBar2::xmlcop = true;
        Tree tree(expr);
        vector<Tree *> trees;
        trees.push_back(&tree);
        buildObjectiveMinimize(ExpressionObjective::EXPRESSION_O, trees);
    }

    void buildObjectiveMaximizeExpression(string expr) override {
        ToulBar2::xmlcop = true;
        Tree tree(expr);
        vector<Tree *> trees;
        trees.push_back(&tree);
        buildObjectiveMaximize(ExpressionObjective::EXPRESSION_O, trees);
    }

    void buildConstraintInstantiation(string id, vector<XVariable *> &list, vector<int> &values) override {
        toMyVariables(list,assignedVars);
        assignedValues.insert(assignedValues.end(), values.begin(), values.end());
    }

    void buildAnnotationDecision(vector<XVariable *> &list) override {}
};

#endif /* SRC_XCSP3_XMLCSP3_H_ */

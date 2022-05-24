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

class MySolverCallbacks : public XCSP3CoreCallbacks {
    public:
        WeightedCSP *problem;
        std::map<string,int> mapping;
        vector<vector<int> > lastTuples;

    void beginInstance(InstanceType type) {
//            problem = makeWeightedCSP(MAX_COST);
        XCSP3CoreCallbacks::intensionUsingString = false;
        XCSP3CoreCallbacks::recognizeSpecialIntensionCases = true;
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

    // transforms a vector of XVariable in vector of toulbar2 variable indices
    void toMyVariables(vector<XVariable*> &src, vector<int> &dest) {
            for(unsigned int i = 0;i<src.size();i++)
                dest.push_back(mapping[src[i]->id]);
    }

    // transforms a vector of string Variable name in vector of toulbar2 variable indices
    void toMyVariables(vector<string> &src, vector<int> &dest) {
            for(unsigned int i = 0;i<src.size();i++)
                dest.push_back(mapping[src[i]]);
    }

    void buildConstraintExtension(string id, vector<string> list, vector<vector<int> > &tuples, bool isSupport, bool hasStar) {
        if(hasStar) {
            cerr << "Sorry tuples with stars not implemented!" << endl;
            throw WrongFileFormat();
        }
        lastTuples = tuples;
        vector<int> vars;
        toMyVariables(list,vars);
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

    void buildConstraintExtension(string id, vector<XVariable *> list, vector<vector<int> > &tuples, bool isSupport, bool hasStar) override {
        vector<string> thelist;
        for (XVariable *var : list) {
            thelist.push_back(var->id);
        }
        buildConstraintExtension(id, thelist, tuples, isSupport, hasStar);
    }

    void buildConstraintExtension(string id, XVariable *variable, vector<int> &tuples, bool isSupport, bool hasStar) override {
        // This function is called for unary extensional constraint.
        if(hasStar) {
            cerr << "Sorry unary tuples with stars not implemented!" << endl;
            throw WrongFileFormat();
        }
        int var = mapping[variable->id];
        vector<Cost> costs(problem->getDomainInitSize(var), (isSupport)?MAX_COST:MIN_COST);
        for(auto value:tuples) {
            if (isSupport) {
                costs[problem->toIndex(var, value)] = MIN_COST;
            } else {
                costs[problem->toIndex(var, value)] = MAX_COST;
            }
        }
        problem->postUnaryConstraint(var, costs);
    }

    // This function is called with group of constraint where the set of tuples is exactly the same
    // than the previous one (then, you can save time/memory using the same set of tuples.
    void buildConstraintExtensionAs(string id, vector<XVariable *> list, bool isSupport, bool hasStar) override {
        buildConstraintExtension(id, list, lastTuples, isSupport, hasStar);
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

    void buildConstraintPrimitive(string id, OrderType op, XVariable *x, int k, XVariable *y) override {
        int varx = mapping[x->id];
        int vary = mapping[y->id];
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
                case OrderType::EQ:
                    costs.push_back((problem->toValue(varx, a) + k == problem->toValue(vary, b))?MIN_COST:MAX_COST);
                    break;
                case OrderType::IN:
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

    void buildConstraintPrimitive(string id, OrderType op, XVariable *x, int k) override {
        int varx = mapping[x->id];
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
            case OrderType::EQ:
                costs.push_back((problem->toValue(varx, a) == k)?MIN_COST:MAX_COST);
                break;
            case OrderType::IN:
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

    void buildConstraintPrimitive(string id, XVariable *x, bool in, int min, int max) override {
        int varx = mapping[x->id];
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

    void buildConstraintMult(string id, XVariable *x, XVariable *y, XVariable *z) override {
        int varx = mapping[x->id];
        int vary = mapping[y->id];
        int varz = mapping[z->id];
        vector<Cost> costs;
        if (varx != vary) {
            assert(varx != varz);
            assert(vary != varz);
            for (unsigned int a = 0; a < problem->getDomainInitSize(varx); a++) {
                for (unsigned int b = 0; b < problem->getDomainInitSize(vary); b++) {
                    for (unsigned int c = 0; c < problem->getDomainInitSize(varz); c++) {
                        costs.push_back((problem->toValue(varx, a) * problem->toValue(vary, b) == problem->toValue(varz, c))?MIN_COST:MAX_COST);
                    }
                }
            }
            problem->postTernaryConstraint(varx, vary, varz, costs);
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

//    Remark: Note that the argument list of each function does not always provides the scope of the constraint: a variable can appear twice
//    in the list for example or a variable can appear in the list and also elsewhere (for example in the values array).
    void buildConstraintAlldifferent(string id, vector<XVariable *> &list) override {
        vector<int> vars;
        toMyVariables(list,vars);
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

    void buildConstraintSum(string id, vector<XVariable *> &list, vector<int> &coeffs, XCondition &cond) override {
        vector<int> vars;
        toMyVariables(list,vars);
        vector<int> myvars = vars;
        vector<int> mycoeffs = coeffs;
        int rightcoef = 0;
        string extravarname = "";
        int extravar = -1;
        string params = "";
        switch (cond.operandType) {
        case OperandType::VARIABLE:
            myvars.push_back(mapping[cond.var]);
            mycoeffs.push_back(-1);
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
                        params += " " + to_string(value) + " " + to_string(-mycoeffs[i] * value);
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
                        params += " " + to_string(value) + " " + to_string(-mycoeffs[i] * value);
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
                        params += " " + to_string(value) + " " + to_string(mycoeffs[i] * value);
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
                        params += " " + to_string(value) + " " + to_string(mycoeffs[i] * value);
                    }
                }
                problem->postKnapsackConstraint(myvars, params, false, true, false);
                break;
            case OrderType::NE:
                extravarname = IMPLICIT_VAR_TAG + "sum" + to_string(problem->numberOfVariables());
                extravar = problem->makeEnumeratedVariable(extravarname, 0, 1);
                mapping[extravarname] = extravar;
                myvars.push_back(extravar);
                mycoeffs.push_back(INT_MAX);
                params = to_string(rightcoef + 1);
                for (unsigned int i=0; i < myvars.size(); i++) {
                    int domsize = problem->getDomainInitSize(myvars[i]);
                    params += " " + to_string(domsize);
                    for (int idval=0; idval < domsize; idval++) {
                        int value = problem->toValue(myvars[i], idval);
                        params += " " + to_string(value) + " " + to_string(mycoeffs[i] * value);
                    }
                }
                problem->postKnapsackConstraint(myvars, params, false, true, false);
                params = to_string((Long)-rightcoef + 1 - INT_MAX);
                for (unsigned int i=0; i < myvars.size(); i++) {
                    int domsize = problem->getDomainInitSize(myvars[i]);
                    params += " " + to_string(domsize);
                    for (int idval=0; idval < domsize; idval++) {
                        int value = problem->toValue(myvars[i], idval);
                        params += " " + to_string(value) + " " + to_string(-mycoeffs[i] * value);
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
                        params += " " + to_string(value) + " " + to_string(mycoeffs[i] * value);
                    }
                }
                problem->postKnapsackConstraint(myvars, params, false, true, false);
                params = to_string(-rightcoef);
                for (unsigned int i=0; i < myvars.size(); i++) {
                    int domsize = problem->getDomainInitSize(myvars[i]);
                    params += " " + to_string(domsize);
                    for (int idval=0; idval < domsize; idval++) {
                        int value = problem->toValue(myvars[i], idval);
                        params += " " + to_string(value) + " " + to_string(-mycoeffs[i] * value);
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
                    params += " " + to_string(value) + " " + to_string(mycoeffs[i] * value);
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
                    params += " " + to_string(value) + " " + to_string(-mycoeffs[i] * value);
                }
            }
            problem->postKnapsackConstraint(myvars, params, false, true, false);
            break;
        default:
            cerr << "Sorry operandType " << cond.operandType << " not implemented in sum constraint!" << endl;
            throw WrongFileFormat();
        }
    }

    void buildConstraintSum(string id, vector<XVariable *> &list, XCondition &cond) override {
        vector<int> coeffs(list.size(), 1);
        buildConstraintSum(id, list, coeffs, cond);
    }

    void buildConstraintElement(string id, OrderType op, vector<XVariable *> &list, int startIndex, XVariable *index, RankType rank, string value) {
        vector<int> vars;
        toMyVariables(list,vars);
        assert(startIndex == 0);
        assert(rank == RankType::ANY);
        int varindex = mapping[index->id];
        assert(problem->getDomainInitSize(varindex) == vars.size());
        int varvalue = mapping[value];
        if (varindex != varvalue) {
            for (unsigned int i=0; i<vars.size(); i++) {
                vector<Cost> costs(problem->getDomainInitSize(varvalue) * problem->getDomainInitSize(varindex) * problem->getDomainInitSize(vars[i]), MIN_COST);
                assert(varindex != vars[i]);
                assert(varvalue != vars[i]);
                assert(problem->toValue(varindex, i) == (Value)i);
                for (unsigned int a=0; a < problem->getDomainInitSize(varvalue); a++) {
                    for (unsigned int b=0; b < problem->getDomainInitSize(vars[i]); b++) {
                        switch (op) {
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
            }
        } else {
            for (unsigned int i=0; i<vars.size(); i++) {
                vector<Cost> costs(problem->getDomainInitSize(varindex) * problem->getDomainInitSize(vars[i]), MIN_COST);
                assert(varindex != vars[i]);
                assert(varvalue != vars[i]);
                assert(problem->toValue(varindex, i) == (Value)i);
                for (unsigned int b=0; b < problem->getDomainInitSize(vars[i]); b++) {
                    switch (op) {
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
    void buildConstraintElement(string id, vector<XVariable *> &list, int startIndex, XVariable *index, RankType rank, XVariable *value) override {
        buildConstraintElement(id, OrderType::EQ, list, startIndex, index, rank, value->id);
    }

    void buildConstraintElement(string id, vector<XVariable *> &list, int startIndex, XVariable *index, RankType rank, int minvalue, int maxvalue) {
        vector<int> vars;
        toMyVariables(list,vars);
        assert(startIndex == 0);
        assert(rank == RankType::ANY);
        int varindex = mapping[index->id];
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
        int varindex = mapping[index->id];
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
        int varindex = mapping[index->id];
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
        int varindex = mapping[index->id];
        int varvalue = mapping[value->id];
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
        int varvalue = mapping[value->id];
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
        int varrowindex = mapping[rowIndex->id];
        int varcolindex = mapping[colIndex->id];
        int varvalue = mapping[value->id];
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
        int varrowindex = mapping[rowIndex->id];
        int varcolindex = mapping[colIndex->id];
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
        int varrowindex = mapping[rowIndex->id];
        int varcolindex = mapping[colIndex->id];
        int varvalue = mapping[value->id];
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

    void buildUnaryCostFunction(Value mult, XVariable *x) {
        int var = mapping[x->id];
        unsigned int domsize = problem->getDomainInitSize(var);
        vector<Cost> costs(domsize, MIN_COST);
        Cost negcost = min(MIN_COST, (Cost)min((Cost)problem->toValue(var, 0) * mult, (Cost)problem->toValue(var, domsize-1) * mult));
        for (unsigned int a=0; a < domsize; a++) {
            costs[a] = (Cost)problem->toValue(var, a) * mult - negcost;
        }
        if (negcost < MIN_COST) {
            problem->decreaseLb(-negcost);
        }
        problem->postUnaryConstraint(var, costs);
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
        vector<Cost> costs;
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
            costs.push_back((Cost)coef * tree->evaluate(tuplemap));
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
//        problem->postCostFunction(vars, tuples); // TODO
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
                vector<Cost> costs(problem->getDomainInitSize(vars[0]) * problem->getDomainInitSize(vars[1]), MIN_COST);
                for (unsigned int a=0; a < problem->getDomainInitSize(vars[0]); a++) {
                    for (unsigned int b=0; b < problem->getDomainInitSize(vars[1]); b++) {
                        Cost cost = (Cost)coefs[0] * problem->toValue(vars[0], a) * coefs[1] * problem->toValue(vars[1], b) * sign;
                        if (cost < negcost) {
                            negcost = cost;
                        }
                    }
                }
                for (unsigned int a=0; a < problem->getDomainInitSize(vars[0]); a++) {
                    for (unsigned int b=0; b < problem->getDomainInitSize(vars[1]); b++) {
                        costs[a * problem->getDomainInitSize(vars[1]) + b] = (Cost)coefs[0] * problem->toValue(vars[0], a) * coefs[1] * problem->toValue(vars[1], b) * sign - negcost;
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
//        vector<int> vars;
//        toMyVariables(list,vars);
//        Cost negcost = MIN_COST;
        vector<WeightedVarValPair> weightFunction;
        switch (type) {
        case ExpressionObjective::EXPRESSION_O:
            assert(trees.size() == 1);
        case ExpressionObjective::SUM_O:
            for (unsigned int i=0; i<trees.size(); i++) {
                buildCostFunction(sign * coefs[i], trees[i]);
            }
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

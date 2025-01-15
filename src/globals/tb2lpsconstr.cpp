#include "tb2lpsconstr.hpp"
#include "core/tb2wcsp.hpp"

LPSConstraint::LPSConstraint(WCSP* wcsp, EnumeratedVariable** scope_in,
    int arity_in, int* constrcounter)
    : LinearConstraint(wcsp, scope_in, arity_in)
    , wcspconstrcounter(constrcounter)
{
    buildIndex();
}

void LPSConstraint::buildIndex()
{
    vector<Value> D;
    count = 0; // total number of domains (number of vars * number of domains)
    count2 = 0; // number of possible domains (union all possible domains)
    mapvar = new map<Value, int>[arity_];
    for (int i = 0; i < arity_; i++) {
        EnumeratedVariable* x = (EnumeratedVariable*)getVar(i);
        for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
            D.push_back(*iterx);
            mapvar[i].insert(pair<Value, int>(*iterx, count++));
        }
    }
    sort(D.begin(), D.end());
    D.erase(unique(D.begin(), D.end()), D.end());

    domainSize = D.size();

    for (vector<Value>::iterator i = D.begin(); i != D.end(); i++) {
        count2++;
    }
}

void LPSConstraint::read(istream& file, bool mult)
{

    string str;
    int nvalues, low, high, windowsi, d;
    file >> nwindows;
#ifndef ILOGCPLEX
    if (wcspconstrcounter)
        (*wcspconstrcounter) = (*wcspconstrcounter) + nwindows;
    cout << "Warning! slinear global cost function skipped... (recompile ToulBar2 with ILOGCPLEX flag)" << endl;
    deconnect();
    if (arity_ != (int)wcsp->numberOfVariables())
        throw WrongFileFormat();
    return;
#endif
    nrows = 0;
    nslacks = 0;
    windowVars = (int**)malloc(sizeof(int*) * nwindows);
    group = (int**)malloc(sizeof(int*) * nwindows);
    kpcoefs = (int***)malloc(sizeof(int**) * nwindows);

    string typeID, defstr;
    for (int i = 0; i < nwindows; i++) {
        file >> windowsi;
        windowSize.push_back(windowsi);
        windowVars[i] = (int*)malloc(sizeof(int) * windowsi);
        for (int j = 0; j < windowsi; j++) {
            file >> d;
            windowVars[i][j] = d;
        }
        file >> d;
        if (d != -1) {
            cerr << "Error occurred in reading slinear" << endl;
            throw WrongFileFormat();
        } else {
            file >> typeID;
            windowType.push_back(typeID);
            if (strcmp(windowType[i].c_str(), "knapsack") != 0 && strcmp(windowType[i].c_str(), "knapsackp") != 0) {
                file >> defstr;
                if (strcmp(defstr.c_str(), "def") == 0 || strcmp(defstr.c_str(), "var") == 0) {
                    file >> d;
                    subdef.push_back(d);
                } else {
                    cerr << "Error occurred in reading slinear def|var" << endl;
                    throw WrongFileFormat();
                }
            }
            if (strcmp(windowType[i].c_str(), "salldiff") == 0) {
                nrows += domainSize;
                nslacks += domainSize;
                sumlow.push_back(domainSize);
                sumhigh.push_back(domainSize);
            } else if (strcmp(windowType[i].c_str(), "samong") == 0) {
                nrows += 1;
                nslacks += 2;
                file >> low >> high >> nvalues;
                if (high < low) {
                    cerr << "Error occurred in reading samong" << endl;
                    throw WrongFileFormat();
                }
                sumlow.push_back(low);
                sumhigh.push_back(high);

                group[i] = (int*)malloc(sizeof(int) * count2);
                for (int j = 0; j < count2; j++) {
                    group[i][j] = 0;
                }
                for (int j = 0; j < nvalues; j++) {
                    file >> d;
                    group[i][d] = 1;
                }
            } else if (strcmp(windowType[i].c_str(), "sgcc") == 0) {
                file >> d;
                for (int j = 0; j < d; j++) {
                    nrows += 1;
                    nslacks += 2;
                    windowType[i] = "samong";
                    file >> nvalues >> low >> high;
                    if (high < low) {
                        cerr << "Error occurred in reading sgcc" << endl;
                        throw WrongFileFormat();
                    }
                    group[i] = (int*)malloc(sizeof(int) * count2);
                    for (int k = 0; k < count2; k++) {
                        group[i][k] = 0;
                    }
                    group[i][nvalues] = 1;
                    sumlow.push_back(low);
                    sumhigh.push_back(high);
                    if (j != d - 1) {
                        i++;
                        nwindows++;
                        windowSize.push_back(windowsi);
                        windowVars[i] = (int*)malloc(sizeof(int) * windowsi);
                        for (int k = 0; k < windowsi; k++) {
                            windowVars[i][k] = windowVars[i - 1][k];
                        }
                        subdef[i] = subdef[i - 1];
                    }
                }
            } else if (strcmp(windowType[i].c_str(), "ssame") == 0) {
                nrows += domainSize;
                nslacks += domainSize * 2;
                file >> low >> high;
                sumlow.push_back(low);
                sumhigh.push_back(high);

                group[i] = (int*)malloc(sizeof(int) * (low + high));
                for (int k = 0; k < low + high; k++) {
                    group[i][k] = -1;
                }
                for (int j = 0; j < low; j++) {
                    file >> d;
                    for (int k = 0; k < windowsi; k++) {
                        if (windowVars[i][k] == d) {
                            if (group[i][k] == -1) {
                                group[i][k] = 0;
                            } else {
                                cerr << "Error occurred in reading ssame" << endl;
                                throw WrongFileFormat();
                            }
                            break;
                        }
                    }
                }
                for (int j = 0; j < high; j++) {
                    file >> d;
                    for (int k = 0; k < windowsi; k++) {
                        if (windowVars[i][k] == d) {
                            if (group[i][k] == -1) {
                                group[i][k] = 1;
                            } else {
                                cerr << "Error occurred in reading ssame" << endl;
                                throw WrongFileFormat();
                            }
                            break;
                        }
                    }
                }
            } else if (strcmp(windowType[i].c_str(), "ssum") == 0) {
                nrows += 1;
                nslacks += 2;
                file >> low >> high >> nvalues;
                if (high < low) {
                    cerr << "Error occurred in reading sum" << endl;
                    throw WrongFileFormat();
                }
                sumlow.push_back(low);
                sumhigh.push_back(high);
                group[i] = (int*)malloc(sizeof(int) * count2);
                for (int j = 0; j < count2; j++) {
                    group[i][j] = 0;
                }
                for (int j = 0; j < nvalues; j++) {
                    file >> d;
                    group[i][j] = d;
                }
            } else if (strcmp(windowType[i].c_str(), "knapsack") == 0) {
                nrows += 1;
                file >> low;
                sumlow.push_back(low);
                sumhigh.push_back(std::numeric_limits<int>::max());
                group[i] = NULL;
                kpcoefs[i] = (int**)malloc(sizeof(int*) * windowSize[i]);
                for (int j = 0; j < windowSize[i]; j++) {
                    kpcoefs[i][j] = (int*)malloc(sizeof(int) * count2);
                    for (int d = 0; d < count2; d++) {
                        kpcoefs[i][j][d] = 0;
                    }
                }
                for (int j = 0; j < windowSize[i]; j++) {
                    file >> d;
                    kpcoefs[i][j][1] = d;
                }
            } else if (strcmp(windowType[i].c_str(), "segcc") == 0) {
                nrows += 1;
                nslacks += 2;
                file >> d >> nvalues;
                sumlow.push_back(d); // abused the sumlow array
                sumhigh.push_back(0);
                group[i] = (int*)malloc(sizeof(int) * count2);
                for (int j = 0; j < count2; j++) {
                    group[i][j] = 0;
                }
                for (int j = 0; j < nvalues; j++) {
                    file >> d;
                    group[i][d] = 1;
                }
            } else if (strcmp(windowType[i].c_str(), "sdisjunctive") == 0) {
                nrows += count2;
                nslacks += count2;
                file >> d;
                sumlow.push_back(d); // abused the sumlow array
                sumhigh.push_back(0);
                group[i] = (int*)malloc(sizeof(int) * d * 2);
                for (int j = 0; j < d; j++) {
                    file >> low >> high;
                    group[i][j] = low;
                    group[i][j + d] = high;
                }
            } else {
                cerr << "Error occurred in reading slinear: no linearization method for: " << typeID << endl;
                throw WrongFileFormat();
            }
        }
    }
}

Cost LPSConstraint::evalOriginal(const Tuple& s)
{

    Cost cost = 0;
    for (int i = 0; i < nwindows; i++) {
        if (strcmp(windowType[i].c_str(), "salldiff") == 0) {
            set<char> count;
            for (int j = 0; j < windowSize[i]; j++) {
                count.insert(s[windowVars[i][j]]);
            }
            cost += (windowSize[i] - count.size()) * subdef[i];
        } else if (strcmp(windowType[i].c_str(), "samong") == 0) {
            int appear = 0;
            for (int j = 0; j < windowSize[i]; j++) {
                if (group[i][s[windowVars[i][j]]]) {
                    appear++;
                }
            }
            if (appear < sumlow[i]) {
                cost += subdef[i] * (sumlow[i] - appear);
            } else if (appear > sumhigh[i]) {
                cost += subdef[i] * (appear - sumhigh[i]);
            }
        } else if (strcmp(windowType[i].c_str(), "ssame") == 0) {
            map<char, int> appear;
            for (int j = 0; j < windowSize[i]; j++) {
                if (group[i][j] == 0) {
                    appear[s[windowVars[i][j]]] += subdef[i];
                } else if (group[i][j] == 1) {
                    appear[s[windowVars[i][j]]] -= subdef[i];
                } else {
                    cerr << "Error occurred in reading ssame()" << endl;
                    throw WrongFileFormat();
                }
            }
            int sum = 0;
            for (map<char, int>::iterator it = appear.begin(); it != appear.end(); it++) {
                sum += (it->second < 0) ? (-(it->second)) : it->second;
            }
            cost += sum / 2;
        } else if (strcmp(windowType[i].c_str(), "ssum") == 0) {
            int tmpSum = 0;
            for (int j = 0; j < windowSize[i]; j++) {
                tmpSum += group[i][s[windowVars[i][j]]];
            }
            if (tmpSum < sumlow[i]) {
                cost += subdef[i] * (sumlow[i] - tmpSum);
            } else if (tmpSum > sumhigh[i]) {
                cost += subdef[i] * (tmpSum - sumhigh[i]);
            }
        } else if (strcmp(windowType[i].c_str(), "knapsack") == 0) {
            int tmpSum = 0;
            for (int j = 0; j < windowSize[i]; j++) {
                tmpSum += kpcoefs[i][j][s[windowVars[i][j]]];
            }
            if (tmpSum < sumlow[i]) {
                cost = wcsp->getUb();
            }
        } else if (strcmp(windowType[i].c_str(), "segcc") == 0) {
            int appear = 0;
            int tmpSum = s[windowVars[i][sumlow[i]]]; // abused the sumlow array
            for (int j = 0; j < windowSize[i]; j++) {
                // cout << s[j];
                if (group[i][s[windowVars[i][j]]]) {
                    appear++;
                }
            }
            if (appear < tmpSum) {
                cost += subdef[i] * (tmpSum - appear);
            } else if (appear > tmpSum) {
                cost += subdef[i] * (appear - tmpSum);
            }
        } else if (strcmp(windowType[i].c_str(), "sdisjunctive") == 0) {
            int tmpSum;
            for (int k = 0; k < count2; k++) {
                tmpSum = 0;

                for (int j = 0; j < sumlow[i]; j++) {
                    if (k >= s[windowVars[i][group[i][j]]] && k < s[windowVars[i][group[i][j]]] + group[i][j + sumlow[i]]) {
                        tmpSum++;
                    }
                }

                if (tmpSum > 1) {
                    cost += (tmpSum - 1) * subdef[i];
                }
            }
        } else {
            cerr << "Error occurred in evaloriginal: Unknown ID" << endl;
            throw WrongFileFormat();
        }
    }
    return cost;
}

Cost LPSConstraint::buildMIP(MIP& mip)
{

    mip.clear();
    mip.addRows(nrows);
    mip.addBool(count);
    mip.addInt(nslacks);

    double* vars;
    double* var2;
    int* idxs;
    vars = new double[count + nslacks];
    var2 = new double[count + nslacks];
    idxs = new int[count + nslacks];

    for (int i = 0; i < count + nslacks; i++) {
        idxs[i] = i;
    }

    int rowCount = 0;
    int slackCount = 0;
    for (int i = 0; i < nwindows; i++) {
        if (strcmp(windowType[i].c_str(), "salldiff") == 0) {
            for (int k = 0; k < domainSize; k++) {
                for (int j = 0; j < count + nslacks; j++) {
                    vars[j] = 0;
                }
                for (int j = 0; j < windowSize[i]; j++) {
                    EnumeratedVariable* x = (EnumeratedVariable*)getVar(windowVars[i][j]);
                    for (EnumeratedVariable::iterator v = x->begin(); v != x->end(); ++v) {
                        if (k == *v) {
                            vars[mapvar[windowVars[i][j]][*v]] = 1;
                        }
                    }
                }
                vars[count + slackCount] = -1;
                mip.objCoeff(count + slackCount, subdef[i]);
                slackCount++;
                mip.rowCoeff(rowCount, count + nslacks, idxs, vars);
                mip.rowBound(rowCount, 0, 1);
                rowCount++;
            }
        } else if (strcmp(windowType[i].c_str(), "samong") == 0) {
            for (int j = 0; j < count + nslacks; j++) {
                vars[j] = 0;
            }
            for (int j = 0; j < windowSize[i]; j++) {
                EnumeratedVariable* x = (EnumeratedVariable*)getVar(windowVars[i][j]);
                for (EnumeratedVariable::iterator v = x->begin(); v != x->end(); ++v) {

                    vars[mapvar[windowVars[i][j]][*v]] = group[i][*v];
                }
            }
            vars[count + slackCount] = 1;
            mip.objCoeff(count + slackCount, subdef[i]);
            slackCount++;
            vars[count + slackCount] = -1;
            mip.objCoeff(count + slackCount, subdef[i]);
            slackCount++;
            mip.rowCoeff(rowCount, count + nslacks, idxs, vars);
            mip.rowBound(rowCount, sumlow[i], sumhigh[i]);
            rowCount++;
        } else if (strcmp(windowType[i].c_str(), "ssame") == 0) {
            for (int k = 0; k < domainSize; k++) {
                for (int j = 0; j < count + nslacks; j++) {
                    vars[j] = 0;
                }
                for (int j = 0; j < windowSize[i]; j++) {
                    EnumeratedVariable* x = (EnumeratedVariable*)getVar(windowVars[i][j]);
                    for (EnumeratedVariable::iterator v = x->begin(); v != x->end(); ++v) {
                        if (k == *v) {
                            if (group[i][j] == 0) {
                                vars[mapvar[windowVars[i][j]][*v]] = 1;
                            } else {
                                vars[mapvar[windowVars[i][j]][*v]] = -1;
                            }
                        }
                    }
                }
                vars[count + slackCount] = -1;
                mip.objCoeff(count + slackCount, subdef[i]);
                slackCount++;
                vars[count + slackCount] = 1;
                mip.objCoeff(count + slackCount, subdef[i]);
                slackCount++;
                mip.rowCoeff(rowCount, count + nslacks, idxs, vars);
                mip.rowBound(rowCount, 0, 0);
                rowCount++;
            }
        } else if (strcmp(windowType[i].c_str(), "ssum") == 0) {
            for (int j = 0; j < count + nslacks; j++) {
                vars[j] = 0;
            }
            for (int j = 0; j < windowSize[i]; j++) {
                EnumeratedVariable* x = (EnumeratedVariable*)getVar(windowVars[i][j]);
                for (EnumeratedVariable::iterator v = x->begin(); v != x->end(); ++v) {
                    vars[mapvar[windowVars[i][j]][*v]] = group[i][*v];
                }
            }
            vars[count + slackCount] = 1;
            mip.objCoeff(count + slackCount, subdef[i]);
            slackCount++;
            vars[count + slackCount] = -1;
            mip.objCoeff(count + slackCount, subdef[i]);
            slackCount++;
            mip.rowCoeff(rowCount, count + nslacks, idxs, vars);
            mip.rowBound(rowCount, sumlow[i], sumhigh[i]);
            rowCount++;
        } else if (strcmp(windowType[i].c_str(), "knapsack") == 0) {
            for (int j = 0; j < count + nslacks; j++) {
                vars[j] = 0;
            }
            for (int j = 0; j < windowSize[i]; j++) {
                EnumeratedVariable* x = (EnumeratedVariable*)getVar(windowVars[i][j]);
                for (EnumeratedVariable::iterator v = x->begin(); v != x->end(); ++v) {
                    vars[mapvar[windowVars[i][j]][*v]] = kpcoefs[i][j][x->toIndex(*v)];
                }
            }
            mip.rowCoeff(rowCount, count + nslacks, idxs, vars);
            mip.rowBound(rowCount, sumlow[i], sumhigh[i]);
            rowCount++;
        } else if (strcmp(windowType[i].c_str(), "segcc") == 0) {
            for (int j = 0; j < count + nslacks; j++) {
                vars[j] = 0;
            }
            for (int j = 0; j < windowSize[i]; j++) {
                EnumeratedVariable* x = (EnumeratedVariable*)getVar(windowVars[i][j]);
                for (EnumeratedVariable::iterator v = x->begin(); v != x->end(); ++v) {

                    vars[mapvar[windowVars[i][j]][*v]] = group[i][*v];
                }
            }
            EnumeratedVariable* x = (EnumeratedVariable*)getVar(sumlow[i]);
            for (EnumeratedVariable::iterator v = x->begin(); v != x->end(); ++v) {
                vars[mapvar[sumlow[i]][*v]] = -(*v);
            }
            vars[count + slackCount] = 1;
            mip.objCoeff(count + slackCount, subdef[i]);
            slackCount++;
            vars[count + slackCount] = -1;
            mip.objCoeff(count + slackCount, subdef[i]);
            slackCount++;
            mip.rowCoeff(rowCount, count + nslacks, idxs, vars);
            mip.rowBound(rowCount, 0, 0);
            rowCount++;
        } else if (strcmp(windowType[i].c_str(), "sdisjunctive") == 0) {
            for (int l = 0; l < count2; l++) {
                for (int j = 0; j < count + nslacks; j++) {
                    vars[j] = 0;
                }
                for (int j = 0; j < sumlow[i]; j++) {
                    int k = l - group[i][j + sumlow[i]];
                    while (k < l) {
                        EnumeratedVariable* x = (EnumeratedVariable*)getVar(windowVars[i][group[i][j]]);
                        for (EnumeratedVariable::iterator v = x->begin(); v != x->end(); ++v) {
                            if (k == *v) {
                                vars[mapvar[windowVars[i][group[i][j]]][*v]] = 1;
                                k++;
                            } else {
                                while (k < l && k != *v) {
                                    k++;
                                }
                                if (k == *v) {
                                    vars[mapvar[windowVars[i][group[i][j]]][*v]] = 1;
                                    k++;
                                }
                            }
                        }
                    }
                }
                vars[count + slackCount] = -1;
                mip.objCoeff(count + slackCount, subdef[i]);
                slackCount++;
                mip.rowCoeff(rowCount, count + nslacks, idxs, vars);
                mip.rowBound(rowCount, 0, 1);
                rowCount++;
            }
        } else {
            cerr << "Error occurred in building mip: Unknown ID" << windowType[i] << endl;
            throw WrongFileFormat();
        }
    }

    for (int i = 0; i < arity_; i++) {
        mip.addRows(1);
        for (int j = 0; j < count + nslacks; j++) {
            var2[j] = 0;
        }
        EnumeratedVariable* x = (EnumeratedVariable*)getVar(i);
        for (EnumeratedVariable::iterator v = x->begin(); v != x->end(); ++v) {
            mip.colUpperBound(mapvar[i][*v], 1);
            var2[mapvar[i][*v]] = 1;
            mip.objCoeff(mapvar[i][*v], -deltaCost[i][*v]);
        }
        mip.rowCoeff(i + nrows, count, idxs, var2);
        mip.rowBound(i + nrows, 1, 1);
    }

    mip.end();
    mip.solve();
    return mip.solValue();
}

Cost LPSConstraint::solveMIP(MIP& mip)
{

    return mip.solValue();
}

void LPSConstraint::dump(ostream& os, bool original)
{
    if (original) {
        os << arity_;
        for (int i = 0; i < arity_; i++)
            os << " " << scope[i]->wcspIndex;
    } else {
        os << getNonAssigned();
        for (int i = 0; i < arity_; i++)
            if (scope[i]->unassigned())
                os << " " << scope[i]->getCurrentVarId();
    }
    os << " -1 slinear " << nwindows << endl;
    for (int i = 0; i < nwindows; i++) {
        os << windowSize[i];
        for (int j = 0; j < windowSize[i]; j++) {
            os << " " << windowVars[i][j];
        }
        os << " -1 " << windowType[i];
        if (windowType[i] == "knapsack") {
            os << " " << sumlow[i];
            for (int j = 0; j < windowSize[i]; j++) {
                os << " " << kpcoefs[i][j][1];
            }
        } else {
            os << " " << subdef[i] << " " << sumlow[i] << " " << sumhigh[i];
            cerr << endl
                 << "sorry, dump function for slinear not fully implemented!!!" << endl;
            throw WrongFileFormat();
        }
        os << endl;
    }
}

void LPSConstraint::print(ostream& os)
{
    int nvalues = 1;

    os << "slinear(";
    for (int i = 0; i < arity_; i++) {
        os << scope[i]->wcspIndex;
        if (i < arity_ - 1)
            os << ",";
    }
    os << ")[" << ((mode == VAR) ? "var" : "dec") << "," << def << "," << nvalues;

    os << ","
       << "bounds"
       << "," << sumlow[0] << "," << sumhigh[0];

    os << "," << nwindows << "," << nrows << "," << nslacks;

    os << "]";
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

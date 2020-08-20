/*
 * **************** Generic solver *******************
 *
 */

#include "tb2solver.hpp"
#include "core/tb2vac.hpp"
#include "core/tb2domain.hpp"
#include "applis/tb2pedigree.hpp"
#include "applis/tb2haplotype.hpp"
#include "applis/tb2bep.hpp"
#include "tb2clusters.hpp"
#include "vns/tb2vnsutils.hpp"
#include "vns/tb2dgvns.hpp"
#ifdef OPENMPI
#include "vns/tb2cpdgvns.hpp"
#include "vns/tb2rpdgvns.hpp"
#endif
#include <unistd.h>

extern void setvalue(int wcspId, int varIndex, Value value, void* solver);

const string Solver::CPOperation[CP_MAX] = { "ASSIGN", "REMOVE", "INCREASE", "DECREASE", "RANGEREMOVAL" };

/*
 * Solver constructors
 *
 */

WeightedCSPSolver* WeightedCSPSolver::makeWeightedCSPSolver(Cost ub)
{
#ifdef OPENMPI
    MPIEnv env0;
    MPI_Comm_size(MPI_COMM_WORLD, &env0.ntasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &env0.myrank);
#endif
    WeightedCSPSolver* solver = NULL;
    switch (ToulBar2::searchMethod) {
    case VNS:
    case DGVNS:
#ifdef BOOST
        solver = new VNSSolver(ub);
#else
        cerr << "Error: compiling with Boost graph library is needed to allow VNS-like search methods." << endl;
        exit(EXIT_FAILURE);
#endif
        break;
#ifdef OPENMPI
    case CPDGVNS:
#ifdef BOOST
        solver = new CooperativeParallelDGVNS(ub, env0);
#else
        cerr << "Error: compiling with Boost graph library is needed to allow VNS-like search methods." << endl;
        exit(EXIT_FAILURE);
#endif
        break;
    case RPDGVNS:
#ifdef BOOST
        solver = new ReplicatedParallelDGVNS(ub, env0);
#else
        cerr << "Error: compiling with Boost graph library is needed to allow VNS-like search methods." << endl;
        exit(EXIT_FAILURE);
#endif
        break;
#endif
    case TREEDEC:
#ifdef BOOST
        solver = new TreeDecRefinement(ub);
#else
        cerr << "Error: compiling with Boost graph library is needed to allow VNS-like search methods." << endl;
        exit(EXIT_FAILURE);
#endif
        break;
    default:
        solver = new Solver(ub);
        break;
    };
    return solver;
}

Solver::Solver(Cost initUpperBound)
    : nbNodes(0)
    , nbBacktracks(0)
    , nbBacktracksLimit(LONGLONG_MAX)
    , wcsp(NULL)
    , allVars(NULL)
    , unassignedVars(NULL)
    , lastConflictVar(-1)
    , nbSol(0.)
    , nbSGoods(0)
    , nbSGoodsUse(0)
    , tailleSep(0)
    , cp(NULL)
    , open(NULL)
    , hbfsLimit(LONGLONG_MAX)
    , nbHybrid(0)
    , nbHybridContinue(0)
    , nbHybridNew(0)
    , nbRecomputationNodes(0)
    , initialLowerBound(MIN_COST)
    , globalLowerBound(MIN_COST)
    , globalUpperBound(MAX_COST)
    , initialDepth(0)
    , prevDivSolutionCost(MIN_COST)
{
    searchSize = new StoreCost(MIN_COST);
    wcsp = WeightedCSP::makeWeightedCSP(initUpperBound, (void*)this);
}

Solver::~Solver()
{
    delete cp;
    delete open;
    delete unassignedVars;
    delete[] allVars;
    delete wcsp;
    delete ((StoreCost*)searchSize);
}

void Solver::initVarHeuristic()
{
    unassignedVars = new BTList<Value>(&Store::storeDomain);
    allVars = new DLink<Value>[wcsp->numberOfVariables()];
    for (unsigned int j = 0; j < wcsp->numberOfVariables(); j++) {
        unsigned int i = wcsp->getDACOrder(j);
        allVars[i].content = j;
    }
    for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
        unassignedVars->push_back(&allVars[i], false);
        if (wcsp->assigned(allVars[i].content) || (ToulBar2::nbDecisionVars > 0 && allVars[i].content >= ToulBar2::nbDecisionVars))
            unassignedVars->erase(&allVars[i], false);
        else
            wcsp->resetWeightedDegree(allVars[i].content);
    }
    // Now function setvalue can be called safely!
    ToulBar2::setvalue = setvalue;
}

Cost Solver::read_wcsp(const char* fileName)
{
    ToulBar2::setvalue = NULL;
    return wcsp->read_wcsp(fileName);
}

void Solver::read_random(int n, int m, vector<int>& p, int seed, bool forceSubModular, string globalname)
{
    ToulBar2::setvalue = NULL;
    wcsp->read_random(n, m, p, seed, forceSubModular, globalname);
}

void Solver::read_solution(const char* filename, bool updateValueHeuristic)
{
    // open the file
    ifstream file(filename);
    if (!file) {
        cout << "Solution file " << filename << " not found.." << endl;
        return;
    }

    wcsp->propagate();

    int depth = Store::getDepth();
    Store::store();

    vector<int> variables;
    vector<Value> values;
    int i = 0;
    while (!file.eof()) {
        if ((unsigned int)i >= wcsp->numberOfVariables())
            break;
        int i_copy = i;
        Value value = 0;
        string token;
        file >> token;
        if (token.length() == 0)
            break;
        if (not isdigit(token[0])) {
            size_t operation = token.find("=");
            if (operation != string::npos) {
                i = wcsp->getVarIndex(token.substr(0, operation));
                if ((unsigned int)i >= wcsp->numberOfVariables()) {
                    cerr << "Solution file incorrect! " << i_copy << " " << token << " " << i << endl;
                    exit(EXIT_FAILURE);
                }
                operation++;
            } else {
                operation = 0;
            }
            if (not isdigit(token[operation])) {
                unsigned int idx = wcsp->toIndex(i, token.substr(operation, string::npos));
                if (idx >= wcsp->getDomainInitSize(i)) {
                    cerr << "Solution file incorrect! " << i_copy << " " << token << " " << i << " " << idx << endl;
                    exit(EXIT_FAILURE);
                }
                value = wcsp->toValue(i, idx);
            } else {
                value = atoi(token.substr(operation, string::npos).c_str());
            }
        } else {
            value = atoi(token.c_str());
        }

        if (ToulBar2::sortDomains && ToulBar2::sortedDomains.find(i) != ToulBar2::sortedDomains.end()) {
            int j = wcsp->getDomainInitSize(i) - 1;
            while (j >= 0) {
                if (ToulBar2::sortedDomains[i][j].value == value)
                    break;
                j--;
            }
            assert(j >= 0);
            value = j;
        }
        if (!file)
            break;
        variables.push_back(i);
        values.push_back(value);
        // side-effect: remember last solution
        if (updateValueHeuristic)
            wcsp->setBestValue(i, value);
        //        if (wcsp->unassigned(i)) {
        //		  assign(i, value);
        //		  // side-effect: remember last solution
        //		  wcsp->setBestValue(i, value);
        //        } else {
        //		    if (wcsp->getValue(i) != value) {
        //			  THROWCONTRADICTION;
        //			} else {
        //			  wcsp->setBestValue(i, value); // side-effect: remember last solution
        //			}
        //        }
        i = i_copy;
        i++;
    }
    bool contradiction = false;
    try {
        wcsp->assignLS(variables, values);
    } catch (const Contradiction&) {
        contradiction = true;
    }
    assert(wcsp->numberOfUnassignedVariables() == 0);
    if (contradiction) {
        if (ToulBar2::verbose >= 0)
            cout << " Input complete assignment " << filename << " is not a valid solution!" << endl;
    } else {
        if (ToulBar2::verbose >= 0)
            cout << " Input solution cost: " << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->getDDualBound() << std::setprecision(DECIMAL_POINT) << " (nb. of unassigned variables: " << wcsp->numberOfUnassignedVariables() << ")" << endl;
        if (ToulBar2::verifyOpt) {
            ToulBar2::verifiedOptimum = wcsp->getLb();
        } else {
            wcsp->updateUb(wcsp->getLb() + ((updateValueHeuristic) ? UNIT_COST : MIN_COST));
        }
    }
    Store::restore(depth);
    if (ToulBar2::verifyOpt) {
        wcsp->setIsPartOfOptimalSolution(true); // must be done after restoring the original problem
    }
}

void Solver::parse_solution(const char* certificate, bool updateValueHeuristic)
{
    wcsp->propagate();

    //  int depth = Store::getDepth();
    //    Store::store();

    //certif2 = index(certif2,',');
    char* certif2;
    char sep[] = ",";
    certif2 = strdup(certificate);
    certif2 = strstr(certif2, sep);

    if (certif2)
        certif2++;

    vector<int> variables;
    vector<Value> values;
    char svar[1024];
    char svalue[1024];
    int var;
    Value value;
    char operation = '\0';
    while ((certif2 != NULL) && (certif2[0] != '\0')) {
        int items = 0;
        char *ope = strpbrk(certif2+1, "=#<>"); // avoid first character of a variable name (can be an operation char)
        if (ope) {
            items++;
            operation = *ope;
            items++;
            strncpy(svar, certif2, ope - certif2);
            svar[ope - certif2] = '\0';
            char *nextsep = strpbrk(ope+2, sep);
            if (nextsep) {
                items++;
                strncpy(svalue, ope+1, nextsep-ope-1);
                svalue[nextsep-ope-1] = '\0';
            } else {
                if (strlen(ope+1) > 0) {
                    items++;
                    strcpy(svalue, ope+1);
                }
            }

        }
        if (items != 3) {
            cerr << "Certificate " << certif2 << " incorrect! " << items << endl;
            exit(EXIT_FAILURE);
        }
        certif2 = strstr(certif2, sep);
        if (certif2)
            certif2++;

        if (not isdigit(svar[0])) {
            var = wcsp->getVarIndex(to_string(svar));
            if ((unsigned int)var >= wcsp->numberOfVariables()) {
                cerr << "Certificate " << certif2 << " incorrect!" << endl;
                exit(EXIT_FAILURE);
            }
        } else {
            var = atoi(svar);
        }

        if (not isdigit(svalue[0])) {
            unsigned int idx = wcsp->toIndex(var, to_string(svalue));
            if (idx >= wcsp->getDomainInitSize(var)) {
                cerr << "Certificate " << certif2 << " incorrect!" << endl;
                exit(EXIT_FAILURE);
            }
            value = wcsp->toValue(var, idx);
        } else {
            value = atoi(svalue);
        }

        if (ToulBar2::sortDomains && ToulBar2::sortedDomains.find(var) != ToulBar2::sortedDomains.end()) {
            int j = wcsp->getDomainInitSize(var) - 1;
            while (j >= 0) {
                if (ToulBar2::sortedDomains[var][j].value == value)
                    break;
                j--;
            }
            assert(j >= 0);
            value = j;
        }

        switch (operation) {
        case '=': {
            variables.push_back(var);
            values.push_back(value);
            // side-effect: remember last solution
            if (updateValueHeuristic)
                wcsp->setBestValue(var, value);
            break;
        }
        case '#': {
            wcsp->remove(var, value);
            break;
        }
        case '>': {
            wcsp->increase(var, value + 1);
            break;
        }
        case '<': {
            wcsp->decrease(var, value - 1);
            break;
        }
        default: {
            cerr << "unknown choice point '" << operation << "' for partial assignment!!!" << endl;
            exit(EXIT_FAILURE);
        }
        }
        //        if (wcsp->unassigned(var)) {
        //          assign(var, value);
        //          // side-effect: remember last solution
        //          wcsp->setBestValue(var, value);
        //        } else {
        //		  if (wcsp->getValue(var) != value) {
        //			THROWCONTRADICTION;
        //		  } else {
        //			wcsp->setBestValue(var, value); // side-effect: remember last solution
        //		  }
        //        }
    }
    wcsp->assignLS(variables, values);
    wcsp->propagate();
    if (ToulBar2::verbose >= 0)
        cout << " Input (partial) assignment bounds: [" << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->getDLb() << "," << wcsp->getDUb() << std::setprecision(DECIMAL_POINT) << "] (nb. of unassigned variables: " << wcsp->numberOfUnassignedVariables() << ")" << endl;

    //    if (ToulBar2::btdMode>=2) wcsp->updateUb(wcsp->getLb()+UNIT_COST);
    //    Store::restore(depth);
}

void Solver::dump_wcsp(const char* fileName, bool original, ProblemFormat format)
{
    ofstream pb(fileName);
    switch (format) {
    case WCSP_FORMAT:
        wcsp->dump(pb, original);
        break;
    case CFN_FORMAT:
        wcsp->dump_CFN(pb, original);
        break;
    default:
        cerr << "Cannot save in this problem format! " << format << endl;
        exit(EXIT_FAILURE);
    }
}

set<int> Solver::getUnassignedVars() const
{
    assert(unassignedVars);
    set<int> res;
    for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
        res.insert(*iter);
    }
    return res;
}
unsigned int Solver::numberOfUnassignedVariables() const
{
    assert(unassignedVars);
    return unassignedVars->getSize();
}

/*
 * Link between solver and wcsp: maintain a backtrackable list of unassigned variable indexes
 *
 */

void setvalue(int wcspId, int varIndex, Value value, void* _solver_)
{
    //    assert(wcspId == 0); // WARNING! assert not compatible with sequential execution of solve() method
    Solver* solver = (Solver*)_solver_;
    unsigned int i = solver->getWCSP()->getDACOrder(varIndex);
    if (!solver->allVars[i].removed) {
        solver->unassignedVars->erase(&solver->allVars[i], true);
    }
}

/*
 * Variable ordering heuristics
 *
 */

/// \defgroup heuristics Variable and value search ordering heuristics
/// \see <em> Boosting Systematic Search by Weighting Constraints </em>. Frederic Boussemart, Fred Hemery, Christophe Lecoutre, Lakhdar Sais. Proc. of ECAI 2004, pages 146-150. Valencia, Spain, 2004.
/// \see <em> Last Conflict Based Reasoning </em>. Christophe Lecoutre, Lakhdar Sais, Sebastien Tabary, Vincent Vidal. Proc. of ECAI 2006, pages 133-137. Trentino, Italy, 2006.

int Solver::getNextUnassignedVar()
{
    //    if (lastConflictVar != -1 && wcsp->unassigned(lastConflictVar)) return lastConflictVar;
    return (unassignedVars->empty()) ? -1 : (*unassignedVars->begin());
}

int Solver::getVarMinDomainDivMaxDegree()
{
    int varIndex = -1;
    Cost worstUnaryCost = MIN_COST;
    double best = MAX_VAL - MIN_VAL;

    for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
        double heuristic = (double)wcsp->getDomainSize(*iter) / (double)(wcsp->getDegree(*iter) + 1);
        if (varIndex < 0 || heuristic < best - epsilon * best
            || (heuristic < best + epsilon * best && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
            best = heuristic;
            varIndex = *iter;
            worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
        }
    }
    return varIndex;
}

int Solver::getVarMinDomainDivMaxDegreeRandomized()
{
    int varIndex = -1;
    Cost worstUnaryCost = MIN_COST;
    double best = MAX_VAL - MIN_VAL;
    int ties[unassignedVars->getSize()];
    int nbties = 0;

    for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
        double heuristic = (double)wcsp->getDomainSize(*iter) / (double)(wcsp->getDegree(*iter) + 1);
        if (varIndex < 0 || heuristic < best - epsilon * best
            || (heuristic < best + epsilon * best && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
            best = heuristic;
            varIndex = *iter;
            nbties = 1;
            ties[0] = varIndex;
            worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
        } else if (heuristic < best + epsilon * best && wcsp->getMaxUnaryCost(*iter) == worstUnaryCost) {
            ties[nbties] = *iter;
            nbties++;
        }
    }
    if (nbties > 1)
        return ties[myrand() % nbties];
    else
        return varIndex;
}

int Solver::getVarMinDomainDivMaxDegreeLastConflict()
{
    if (lastConflictVar != -1 && wcsp->unassigned(lastConflictVar))
        return lastConflictVar;
    int varIndex = -1;
    Cost worstUnaryCost = MIN_COST;
    double best = MAX_VAL - MIN_VAL;
    for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
        // remove following "+1" when isolated variables are automatically assigned
        double heuristic = (double)wcsp->getDomainSize(*iter) / (double)(wcsp->getDegree(*iter) + 1);
        if (varIndex < 0 || heuristic < best - epsilon * best
            || (heuristic < best + epsilon * best && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
            best = heuristic;
            varIndex = *iter;
            worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
        }
    }
    return varIndex;
}

int Solver::getVarMinDomainDivMaxDegreeLastConflictRandomized()
{
    if (lastConflictVar != -1 && wcsp->unassigned(lastConflictVar))
        return lastConflictVar;
    int varIndex = -1;
    Cost worstUnaryCost = MIN_COST;
    double best = MAX_VAL - MIN_VAL;
    int ties[unassignedVars->getSize()];
    int nbties = 0;

    for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
        // remove following "+1" when isolated variables are automatically assigned
        double heuristic = (double)wcsp->getDomainSize(*iter) / (double)(wcsp->getDegree(*iter) + 1);
        if (varIndex < 0 || heuristic < epsilon * best
            || (heuristic < best + epsilon * best && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
            best = heuristic;
            varIndex = *iter;
            nbties = 1;
            ties[0] = varIndex;
            worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
            //        } else if ((heuristic < best + epsilon * best && wcsp->getMaxUnaryCost(*iter) == worstUnaryCost) || ((myrand()%100)==0)) {
        } else if (heuristic < best + epsilon * best && wcsp->getMaxUnaryCost(*iter) == worstUnaryCost) {
            ties[nbties] = *iter;
            nbties++;
        }
    }
    if (nbties > 1) {
        return ties[myrand() % nbties];
    } else
        return varIndex;
}

int Solver::getVarMinDomainDivMaxWeightedDegree()
{
    int varIndex = -1;
    Cost worstUnaryCost = MIN_COST;
    double best = MAX_VAL - MIN_VAL;

    for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
        Cost unarymediancost = MIN_COST;
        int domsize = wcsp->getDomainSize(*iter);
        if (ToulBar2::weightedTightness) {
            ValueCost array[domsize];
            wcsp->getEnumDomainAndCost(*iter, array);
            unarymediancost = stochastic_selection<ValueCost>(array, 0, domsize - 1, domsize / 2).cost;
        }
        double heuristic = (double)domsize / (double)(wcsp->getWeightedDegree(*iter) + 1 + unarymediancost);
        if (varIndex < 0 || heuristic < best - epsilon * best
            || (heuristic < best + epsilon * best && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
            best = heuristic;
            varIndex = *iter;
            worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
        }
    }
    return varIndex;
}

int Solver::getVarMinDomainDivMaxWeightedDegreeRandomized()
{
    int varIndex = -1;
    Cost worstUnaryCost = MIN_COST;
    double best = MAX_VAL - MIN_VAL;
    int ties[unassignedVars->getSize()];
    int nbties = 0;

    for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
        Cost unarymediancost = MIN_COST;
        int domsize = wcsp->getDomainSize(*iter);
        if (ToulBar2::weightedTightness) {
            ValueCost array[domsize];
            wcsp->getEnumDomainAndCost(*iter, array);
            unarymediancost = stochastic_selection<ValueCost>(array, 0, domsize - 1, domsize / 2).cost;
        }
        double heuristic = (double)domsize / (double)(wcsp->getWeightedDegree(*iter) + 1 + unarymediancost);
        if (varIndex < 0 || heuristic < best - epsilon * best
            || (heuristic < best + epsilon * best && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
            best = heuristic;
            varIndex = *iter;
            nbties = 1;
            ties[0] = varIndex;
            worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
        } else if (heuristic < best + epsilon * best && wcsp->getMaxUnaryCost(*iter) == worstUnaryCost) {
            ties[nbties] = *iter;
            nbties++;
        }
    }
    if (nbties > 1)
        return ties[myrand() % nbties];
    else
        return varIndex;
}

int Solver::getVarMinDomainDivMaxWeightedDegreeLastConflict()
{
    if (lastConflictVar != -1 && wcsp->unassigned(lastConflictVar))
        return lastConflictVar;
    int varIndex = -1;
    Cost worstUnaryCost = MIN_COST;
    double best = MAX_VAL - MIN_VAL;
    for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
        Cost unarymediancost = MIN_COST;
        int domsize = wcsp->getDomainSize(*iter);
        if (ToulBar2::weightedTightness) {
            ValueCost array[domsize];
            wcsp->getEnumDomainAndCost(*iter, array);
            unarymediancost = stochastic_selection<ValueCost>(array, 0, domsize - 1, domsize / 2).cost;
        }
        //remove following "+1" when isolated variables are automatically assigned
        Long wdeg = wcsp->getWeightedDegree(*iter);
        double heuristic = (double)domsize / (double)(wdeg + 1 + unarymediancost);
        //double heuristic = 1. / (double) (wcsp->getMaxUnaryCost(*iter) + 1);
        if (ToulBar2::FullEAC) {
            EnumeratedVariable* var = (EnumeratedVariable*)((WCSP*)wcsp)->getVar(*iter);
            if (!var->isFullEAC()
                && ((varIndex < 0)
                       || (heuristic < best - epsilon * best)
                       || (heuristic < best + epsilon * best && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost))) {
                best = heuristic;
                varIndex = *iter;
                worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
            }
        } else {
            if ((varIndex < 0)
                || (heuristic < best - epsilon * best)
                || (heuristic < best + epsilon * best && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
                best = heuristic;
                varIndex = *iter;
                worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
            }
        }
    }

    if (varIndex == -1) {
        if (ToulBar2::FullEAC) {
            if (!unassignedVars->empty()) {
                if (ToulBar2::verbose >= 2)
                    cout << "Fast greedy assignment for " << unassignedVars->getSize() << " variables!" << endl;
                Cost currentUb = wcsp->getUb();
                Cost newUb = currentUb;
                int depth = Store::getDepth();
                try {
                    Store::store();
                    vector<int> variables;
                    vector<Value> values;
                    for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
                        variables.push_back(*iter);
                        values.push_back(((EnumeratedVariable*)((WCSP*)wcsp)->getVar(*iter))->getSupport());
                    }
                    // Fast Greedy Assignment
                    wcsp->assignLS(variables, values);
                    nbNodes++;
                    newSolution(); /* it will update ub */
                    newUb = wcsp->getUb();
                } catch (const Contradiction&) {
                    wcsp->whenContradiction();
                }
                Store::restore(depth);
                if (newUb < currentUb) { /* a better solution has been found */
                    wcsp->enforceUb(); /* it will generate a contradiction if lb >= ub */
                    wcsp->propagate(); /* it will generate a contradiction if lb >= ub */
                }
                if (unassignedVars->empty()) // a new solution was found and all vars assigned by propagation
                    varIndex = -1;
                else {
                    // Wrong heuristic guess
                    varIndex = getVarMinDomainDivMaxWeightedDegreeLastConflictRandomized();
                    assert(varIndex != -1);
                }
            }
        }
    }
    return varIndex;
}

int Solver::getVarMinDomainDivMaxWeightedDegreeLastConflictRandomized()
{
    if (lastConflictVar != -1 && wcsp->unassigned(lastConflictVar))
        return lastConflictVar;
    int varIndex = -1;
    Cost worstUnaryCost = MIN_COST;
    double best = MAX_VAL - MIN_VAL;
    int ties[unassignedVars->getSize()];
    int nbties = 0;

    for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
        Cost unarymediancost = MIN_COST;
        int domsize = wcsp->getDomainSize(*iter);
        if (ToulBar2::weightedTightness) {
            ValueCost array[domsize];
            wcsp->getEnumDomainAndCost(*iter, array);
            unarymediancost = stochastic_selection<ValueCost>(array, 0, domsize - 1, domsize / 2).cost;
        }
        // remove following "+1" when isolated variables are automatically assigned
        double heuristic = (double)domsize / (double)(wcsp->getWeightedDegree(*iter) + 1 + unarymediancost);
        if (varIndex < 0 || heuristic < best - epsilon * best
            || (heuristic < best + epsilon * best && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
            best = heuristic;
            varIndex = *iter;
            nbties = 1;
            ties[0] = varIndex;
            worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
            //       } else if ((heuristic < best + epsilon * best && wcsp->getMaxUnaryCost(*iter) == worstUnaryCost) || ((myrand()%100)==0)) {
        } else if (heuristic < best + epsilon * best && wcsp->getMaxUnaryCost(*iter) == worstUnaryCost) {
            ties[nbties] = *iter;
            nbties++;
        }
    }
    if (nbties > 1) {
        return ties[myrand() % nbties];
    } else
        return varIndex;
}

int Solver::getMostUrgent()
{
    int varIndex = -1;
    Value best = MAX_VAL;
    Cost worstUnaryCost = MIN_COST;

    for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
        if (varIndex < 0 || wcsp->getInf(*iter) < best || (wcsp->getInf(*iter) == best && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
            best = wcsp->getInf(*iter);
            worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
            varIndex = *iter;
        }
    }
    return varIndex;
}

/*
 * Choice points
 *
 */

/// \brief Enforce WCSP upper-bound and backtrack if ub <= lb or in the case of probabilistic inference if the contribution is too small
void Solver::enforceUb()
{
    wcsp->enforceUb();
    if (ToulBar2::isZ) {
        Cost newCost = wcsp->getLb() + wcsp->getNegativeLb();
        for (BTList<Value>::iterator iter_variable = unassignedVars->begin(); iter_variable != unassignedVars->end(); ++iter_variable) {
            if (wcsp->enumerated(*iter_variable)) {
                EnumeratedVariable* var = (EnumeratedVariable*)((WCSP*)wcsp)->getVar(*iter_variable);
                Cost sumUnaryCosts = MAX_COST;
                for (EnumeratedVariable::iterator iter_value = var->begin(); iter_value != var->end(); ++iter_value) {
                    sumUnaryCosts = wcsp->LogSumExp(sumUnaryCosts, var->getCost(*iter_value));
                }
                newCost += sumUnaryCosts;
            } else {
                newCost += wcsp->LogProb2Cost(Log(wcsp->getDomainSize(*iter_variable)));
            }
        }
        TLogProb newlogU = wcsp->LogSumExp(ToulBar2::logU, newCost);
        if (newlogU < ToulBar2::logepsilon + ToulBar2::logZ) {
            if (ToulBar2::verbose >= 1)
                cout << "ZCUT " << newlogU << " " << ToulBar2::logZ << " " << Store::getDepth() << endl;
            ToulBar2::logU = newlogU;
            THROWCONTRADICTION;
        }
    }
}

void Solver::increase(int varIndex, Value value, bool reverse)
{
    enforceUb();
    nbNodes++;
    if (ToulBar2::verbose >= 1) {
        if (ToulBar2::verbose >= 2)
            cout << *wcsp;
        if (ToulBar2::debug >= 3) {
            string pbname = "problem" + to_string(nbNodes) + ".wcsp";
            ofstream pb(pbname.c_str());
            wcsp->dump(pb);
            cout << " #" << nbNodes;
        }
        cout << "[" << Store::getDepth() << "," << wcsp->getLb() << "," << wcsp->getUb() << "," << wcsp->getDomainSizeSum();
        if (wcsp->getTreeDec())
            cout << ",C" << wcsp->getTreeDec()->getCurrentCluster()->getId();
        cout << "] Try " << wcsp->getName(varIndex) << " >= " << value << " (s:" << wcsp->getSupport(varIndex) << ")" << endl;
    }
    wcsp->increase(varIndex, value);
    wcsp->propagate();
    if (ToulBar2::hbfs)
        addChoicePoint(CP_INCREASE, varIndex, value, reverse);
}

void Solver::decrease(int varIndex, Value value, bool reverse)
{
    enforceUb();
    nbNodes++;
    if (ToulBar2::verbose >= 1) {
        if (ToulBar2::verbose >= 2)
            cout << *wcsp;
        if (ToulBar2::debug >= 3) {
            string pbname = "problem" + to_string(nbNodes) + ".wcsp";
            ofstream pb(pbname.c_str());
            wcsp->dump(pb);
            cout << " #" << nbNodes;
        }
        cout << "[" << Store::getDepth() << "," << wcsp->getLb() << "," << wcsp->getUb() << "," << wcsp->getDomainSizeSum();
        if (wcsp->getTreeDec())
            cout << ",C" << wcsp->getTreeDec()->getCurrentCluster()->getId();
        cout << "] Try " << wcsp->getName(varIndex) << " <= " << value << " (s:" << wcsp->getSupport(varIndex) << ")" << endl;
    }
    wcsp->decrease(varIndex, value);
    wcsp->propagate();
    if (ToulBar2::hbfs)
        addChoicePoint(CP_DECREASE, varIndex, value, reverse);
}

void Solver::assign(int varIndex, Value value, bool reverse)
{
    enforceUb();
    nbNodes++;
    if (ToulBar2::debug && ((nbNodes % 128) == 0)) {
        if (isatty(fileno(stdout)))
            cout << "\r";
        cout << Store::getDepth();
        if (ToulBar2::hbfs) {
            if (wcsp->getTreeDec()) {
                Cost delta = wcsp->getTreeDec()->getCurrentCluster()->getCurrentDelta();
                if (wcsp->getTreeDec()->getCurrentCluster()->open->size() > 0)
                    cout << " [" << wcsp->getTreeDec()->getCurrentCluster()->open->getLb(delta) << "," << wcsp->getUb() << "]/" << wcsp->getTreeDec()->getCurrentCluster()->open->size() << "/" << wcsp->getTreeDec()->getCurrentCluster()->cp->size() << " " << (100. * (wcsp->getUb() - wcsp->getTreeDec()->getCurrentCluster()->open->getLb(delta)) / wcsp->getUb()) << "%";
            } else {
                if (open->size() > 0)
                    cout << " [" << open->getLb() << "," << wcsp->getUb() << "]/" << open->size() << "/" << cp->size() << "/" << nbNodes << " " << (100. * (wcsp->getUb() - open->getLb()) / wcsp->getUb()) << "%";
            }
        } else if (ToulBar2::vnsKmax > 0) {
            cout << " " << ToulBar2::vnsKcur << " " << ToulBar2::vnsLDScur;
        }
        cout << " " << Exp(((Cost)(*((StoreCost*)searchSize))) / 10e6);
        if (wcsp->getTreeDec())
            cout << " C" << wcsp->getTreeDec()->getCurrentCluster()->getId();
        if (isatty(fileno(stdout)))
            cout << "             ";
        else
            cout << endl;
        cout.flush();
    }
    if (ToulBar2::verbose >= 1) {
        if (ToulBar2::verbose >= 2)
            cout << *wcsp;
        if (ToulBar2::debug >= 3) {
            string pbname = "problem" + to_string(nbNodes) + ".wcsp";
            ofstream pb(pbname.c_str());
            wcsp->dump(pb);
            cout << " #" << nbNodes;
        }
        cout << "[" << Store::getDepth() << "," << wcsp->getLb() << "," << wcsp->getUb() << "," << wcsp->getDomainSizeSum();
        if (wcsp->getTreeDec())
            cout << ",C" << wcsp->getTreeDec()->getCurrentCluster()->getId();
        cout << "] Try " << wcsp->getName(varIndex) << " == " << value << endl;
    }
    wcsp->assign(varIndex, value);
    wcsp->propagate();
    if (ToulBar2::hbfs)
        addChoicePoint(CP_ASSIGN, varIndex, value, reverse);
}

void Solver::remove(int varIndex, Value value, bool reverse)
{
    enforceUb();
    nbNodes++;
    if (ToulBar2::verbose >= 1) {
        if (ToulBar2::verbose >= 2)
            cout << *wcsp;
        if (ToulBar2::debug >= 3) {
            string pbname = "problem" + to_string(nbNodes) + ".wcsp";
            ofstream pb(pbname.c_str());
            wcsp->dump(pb);
            cout << " #" << nbNodes;
        }
        cout << "[" << Store::getDepth() << "," << wcsp->getLb() << "," << wcsp->getUb() << "," << wcsp->getDomainSizeSum();
        if (wcsp->getTreeDec())
            cout << ",C" << wcsp->getTreeDec()->getCurrentCluster()->getId();
        cout << "] Try " << wcsp->getName(varIndex) << " != " << value << endl;
    }
    wcsp->remove(varIndex, value);
    wcsp->propagate();
    if (ToulBar2::hbfs)
        addChoicePoint(CP_REMOVE, varIndex, value, reverse);
}

void Solver::remove(int varIndex, ValueCost* array, int first, int last, bool reverse)
{
    enforceUb();
    nbNodes++;
    if (ToulBar2::verbose >= 1) {
        if (ToulBar2::verbose >= 2)
            cout << *wcsp;
        if (ToulBar2::debug >= 3) {
            string pbname = "problem" + to_string(nbNodes) + ".wcsp";
            ofstream pb(pbname.c_str());
            wcsp->dump(pb);
            cout << " #" << nbNodes;
        }
        cout << "[" << Store::getDepth() << "," << wcsp->getLb() << "," << wcsp->getUb() << "," << wcsp->getDomainSizeSum();
        if (wcsp->getTreeDec())
            cout << ",C" << wcsp->getTreeDec()->getCurrentCluster()->getId();
        cout << "] Try " << wcsp->getName(varIndex) << " !=";
        for (int i = first; i <= last; i++)
            cout << " " << array[i].value;
        cout << endl;
    }
    for (int i = first; i <= last; i++)
        wcsp->remove(varIndex, array[i].value);
    wcsp->propagate();
    if (ToulBar2::hbfs)
        addChoicePoint(CP_REMOVE_RANGE, varIndex, array[first].value, reverse); // Warning! only first value memorized!
}

int cmpValueCost(const void* p1, const void* p2)
{
    Cost c1 = ((ValueCost*)p1)->cost;
    Cost c2 = ((ValueCost*)p2)->cost;
    Value v1 = ((ValueCost*)p1)->value;
    Value v2 = ((ValueCost*)p2)->value;
    if (c1 < c2)
        return -1;
    else if (c1 > c2)
        return 1;
    else if (v1 < v2)
        return -1;
    else if (v1 > v2)
        return 1;
    else
        return 0;
}

void Solver::initGap(Cost newLb, Cost newUb)
{
    initialLowerBound = newLb;
    globalLowerBound = newLb;
    globalUpperBound = newUb;
    initialDepth = Store::getDepth();
}

void Solver::showGap(Cost newLb, Cost newUb)
{
    if (newLb > newUb)
        newLb = newUb;
    if (newUb > initialLowerBound && Store::getDepth() == initialDepth) {
        int oldgap = (int)(100. - 100. * (globalLowerBound - initialLowerBound) / (globalUpperBound - initialLowerBound));
        globalLowerBound = MAX(globalLowerBound, newLb);
        globalUpperBound = MIN(globalUpperBound, newUb);
        int newgap = (int)(100. - 100. * (globalLowerBound - initialLowerBound) / (globalUpperBound - initialLowerBound));
        if (ToulBar2::verbose >= 0 && newgap < oldgap) {
            Double Dglb = (ToulBar2::costMultiplier >= 0 ? wcsp->Cost2ADCost(globalLowerBound) : wcsp->Cost2ADCost(globalUpperBound));
            Double Dgub = (ToulBar2::costMultiplier >= 0 ? wcsp->Cost2ADCost(globalUpperBound) : wcsp->Cost2ADCost(globalLowerBound));
            std::ios_base::fmtflags f(cout.flags());
            cout << "Optimality gap: [" << std::fixed << std::setprecision(ToulBar2::decimalPoint) << Dglb << ", " << Dgub << "] " << std::setprecision(DECIMAL_POINT) << (100. * (Dgub - Dglb)) / max(fabsl(Dglb), fabsl(Dgub)) << " % (" << nbBacktracks << " backtracks, " << nbNodes << " nodes)" << endl;
            cout.flags(f);
        }
    }
}

void Solver::binaryChoicePoint(int varIndex, Value value, Cost lb)
{
    assert(wcsp->unassigned(varIndex));
    assert(wcsp->canbe(varIndex, value));
    if (ToulBar2::interrupted)
        throw TimeOut();
    unsigned int domsize = wcsp->getDomainSize(varIndex);
    bool dichotomic = (ToulBar2::dichotomicBranching && ToulBar2::dichotomicBranchingSize < domsize);
    Value middle = domsize / 2;
    bool increasing = true;
    ValueCost sorted[domsize];
    //	bool reverse = true; // (ToulBar2::restart>0);
    if (dichotomic) {
        if (ToulBar2::dichotomicBranching == 1) {
            middle = (wcsp->getInf(varIndex) + wcsp->getSup(varIndex)) / 2;
            //          if (value <= middle || reverse) increasing = true;
            if (value <= middle)
                increasing = true;
            else
                increasing = false;
        } else if (ToulBar2::dichotomicBranching == 2) {
            wcsp->getEnumDomainAndCost(varIndex, sorted);
            qsort(sorted, domsize, sizeof(ValueCost), cmpValueCost);
        }
        //    } else if (reverse) {
        //    	value = wcsp->getMaxUnaryCostValue(varIndex);
        //		assert(wcsp->canbe(varIndex,value));
    }
    int storedepth = Store::getDepth();
    try {
        Store::store();
        lastConflictVar = varIndex;
        if (dichotomic) {
            if (ToulBar2::dichotomicBranching == 1) {
                if (increasing)
                    decrease(varIndex, middle);
                else
                    increase(varIndex, middle + 1);
            } else if (ToulBar2::dichotomicBranching == 2) {
                if (increasing)
                    remove(varIndex, sorted, middle, domsize - 1);
                else
                    remove(varIndex, sorted, 0, middle - 1);
            }
            //    	} else if (reverse) {
            //    		remove(varIndex, value);
        } else
            assign(varIndex, value);
        lastConflictVar = -1;
        recursiveSolve(lb);
    } catch (const Contradiction&) {
        wcsp->whenContradiction();
    }
    Store::restore(storedepth);
    enforceUb();
    nbBacktracks++;
    if (nbBacktracks > nbBacktracksLimit)
        throw NbBacktracksOut();
#ifdef OPENMPI
    if (ToulBar2::vnsParallel && ((nbBacktracks % 128) == 0) && MPI_interrupted())
        throw TimeOut();
#endif
    if (dichotomic) {
        if (ToulBar2::dichotomicBranching == 1) {
            if (increasing)
                increase(varIndex, middle + 1, nbBacktracks >= hbfsLimit);
            else
                decrease(varIndex, middle, nbBacktracks >= hbfsLimit);
        } else if (ToulBar2::dichotomicBranching == 2) {
            if (increasing)
                remove(varIndex, sorted, 0, middle - 1, nbBacktracks >= hbfsLimit);
            else
                remove(varIndex, sorted, middle, domsize - 1, nbBacktracks >= hbfsLimit);
        }
        //    } else if (reverse) {
        //    	assign(varIndex, value, nbBacktracks >= hybridBFSLimit);
    } else
        remove(varIndex, value, nbBacktracks >= hbfsLimit);
    if (!ToulBar2::hbfs)
        showGap(wcsp->getLb(), wcsp->getUb());
    if (nbBacktracks >= hbfsLimit)
        addOpenNode(*cp, *open, MAX(lb, wcsp->getLb()));
    else
        recursiveSolve(lb);
}

void Solver::binaryChoicePointLDS(int varIndex, Value value, int discrepancy)
{
    assert(wcsp->unassigned(varIndex));
    assert(wcsp->canbe(varIndex, value));
    if (ToulBar2::interrupted)
        throw TimeOut();
    unsigned int domsize = wcsp->getDomainSize(varIndex);
    bool dichotomic = (ToulBar2::dichotomicBranching && ToulBar2::dichotomicBranchingSize < domsize);
    Value middle = domsize / 2;
    bool increasing = true;
    ValueCost sorted[domsize];
    if (dichotomic) {
        if (ToulBar2::dichotomicBranching == 1) {
            middle = (wcsp->getInf(varIndex) + wcsp->getSup(varIndex)) / 2;
            if (value <= middle)
                increasing = true;
            else
                increasing = false;
        } else if (ToulBar2::dichotomicBranching == 2) {
            wcsp->getEnumDomainAndCost(varIndex, sorted);
            qsort(sorted, domsize, sizeof(ValueCost), cmpValueCost);
        }
    }
    if (discrepancy > 0) {
        int storedepth = Store::getDepth();
        try {
            Store::store();
            lastConflictVar = varIndex;
            if (dichotomic) {
                if (ToulBar2::dichotomicBranching == 1) {
                    if (increasing)
                        increase(varIndex, middle + 1);
                    else
                        decrease(varIndex, middle);
                } else if (ToulBar2::dichotomicBranching == 2) {
                    if (increasing)
                        remove(varIndex, sorted, 0, middle - 1);
                    else
                        remove(varIndex, sorted, middle, domsize - 1);
                }
            } else
                remove(varIndex, value);
            lastConflictVar = -1;
            recursiveSolveLDS(discrepancy - 1);
        } catch (const Contradiction&) {
            wcsp->whenContradiction();
        }
        Store::restore(storedepth);
        enforceUb();
        nbBacktracks++;
        if (nbBacktracks > nbBacktracksLimit)
            throw NbBacktracksOut();
#ifdef OPENMPI
        if (ToulBar2::vnsParallel && ((nbBacktracks % 128) == 0) && MPI_interrupted())
            throw TimeOut();
#endif
        if (dichotomic) {
            if (ToulBar2::dichotomicBranching == 1) {
                if (increasing)
                    decrease(varIndex, middle);
                else
                    increase(varIndex, middle + 1);
            } else if (ToulBar2::dichotomicBranching == 2) {
                if (increasing)
                    remove(varIndex, sorted, middle, domsize - 1);
                else
                    remove(varIndex, sorted, 0, middle - 1);
            }
        } else
            assign(varIndex, value);
        if (!ToulBar2::limited)
            showGap(wcsp->getLb(), wcsp->getUb());
        recursiveSolveLDS(discrepancy);
    } else {
        ToulBar2::limited = true;
        lastConflictVar = varIndex;
        if (dichotomic) {
            if (ToulBar2::dichotomicBranching == 1) {
                if (increasing)
                    decrease(varIndex, middle);
                else
                    increase(varIndex, middle + 1);
            } else if (ToulBar2::dichotomicBranching == 2) {
                if (increasing)
                    remove(varIndex, sorted, middle, domsize - 1);
                else
                    remove(varIndex, sorted, 0, middle - 1);
            }
        } else
            assign(varIndex, value);
        lastConflictVar = -1;
        recursiveSolveLDS(0);
    }
}

Value Solver::postponeRule(int varIndex)
{
    assert(ToulBar2::bep);
    Value best = ToulBar2::bep->latest[varIndex] + 1;

    for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
        if (*iter != varIndex) {
            Value time = wcsp->getInf(*iter) + ToulBar2::bep->duration[*iter] + ToulBar2::bep->delay[*iter * ToulBar2::bep->size + varIndex];
            if (time < best) {
                best = time;
            }
        }
    }
    return best;
}

void Solver::scheduleOrPostpone(int varIndex)
{
    assert(wcsp->unassigned(varIndex));
    if (ToulBar2::interrupted)
        throw TimeOut();
    Value xinf = wcsp->getInf(varIndex);
    Value postponeValue = postponeRule(varIndex);
    postponeValue = max(postponeValue, xinf + 1);
    assert(postponeValue <= ToulBar2::bep->latest[varIndex] + 1);
    bool reverse = (wcsp->getUnaryCost(varIndex, xinf) > MIN_COST) ? true : false;
    int storedepth = Store::getDepth();
    try {
        Store::store();
        if (reverse)
            increase(varIndex, postponeValue);
        else
            assign(varIndex, xinf);
        recursiveSolve();
    } catch (const Contradiction&) {
        wcsp->whenContradiction();
    }
    Store::restore(storedepth);
    enforceUb();
    nbBacktracks++;
    if (nbBacktracks > nbBacktracksLimit)
        throw NbBacktracksOut();
#ifdef OPENMPI
    if (ToulBar2::vnsParallel && ((nbBacktracks % 128) == 0) && MPI_interrupted())
        throw TimeOut();
#endif
    if (reverse)
        assign(varIndex, xinf);
    else
        increase(varIndex, postponeValue);
    recursiveSolve();
}

void Solver::narySortedChoicePoint(int varIndex, Cost lb)
{
    assert(wcsp->enumerated(varIndex));
    int size = wcsp->getDomainSize(varIndex);
    ValueCost sorted[size];
    //ValueCost* sorted = new ValueCost [size];
    wcsp->getEnumDomainAndCost(varIndex, sorted);
    qsort(sorted, size, sizeof(ValueCost), cmpValueCost);
    for (int v = 0; wcsp->getLb() < wcsp->getUb() && v < size; v++) {
        if (ToulBar2::interrupted)
            throw TimeOut();
        int storedepth = Store::getDepth();
        try {
            Store::store();
            assign(varIndex, sorted[v].value);
            recursiveSolve(lb);
        } catch (const Contradiction&) {
            wcsp->whenContradiction();
        }
        Store::restore(storedepth);
    }
    //delete [] sorted;
    enforceUb();
    nbBacktracks++;
    if (nbBacktracks > nbBacktracksLimit)
        throw NbBacktracksOut();
#ifdef OPENMPI
    if (ToulBar2::vnsParallel && ((nbBacktracks % 128) == 0) && MPI_interrupted())
        throw TimeOut();
#endif
}

void Solver::narySortedChoicePointLDS(int varIndex, int discrepancy)
{
    assert(wcsp->enumerated(varIndex));
    int size = wcsp->getDomainSize(varIndex);
    ValueCost sorted[size];
    //ValueCost* sorted = new ValueCost [size];
    wcsp->getEnumDomainAndCost(varIndex, sorted);
    qsort(sorted, size, sizeof(ValueCost), cmpValueCost);
    if (discrepancy < size - 1)
        ToulBar2::limited = true;
    for (int v = min(size - 1, discrepancy); wcsp->getLb() < wcsp->getUb() && v >= 0; v--) {
        if (ToulBar2::interrupted)
            throw TimeOut();
        int storedepth = Store::getDepth();
        try {
            Store::store();
            assign(varIndex, sorted[v].value);
            recursiveSolveLDS(discrepancy - v);
        } catch (const Contradiction&) {
            wcsp->whenContradiction();
        }
        Store::restore(storedepth);
    }
    //delete [] sorted;
    enforceUb();
    nbBacktracks++;
    if (nbBacktracks > nbBacktracksLimit)
        throw NbBacktracksOut();
#ifdef OPENMPI
    if (ToulBar2::vnsParallel && ((nbBacktracks % 128) == 0) && MPI_interrupted())
        throw TimeOut();
#endif
}

void Solver::singletonConsistency()
{
    bool deadend;
    bool done = false;
    while (!done) {
        done = true;
        for (unsigned int varIndex = 0; varIndex < ((ToulBar2::nbDecisionVars > 0) ? ToulBar2::nbDecisionVars : wcsp->numberOfVariables()); varIndex++) {
            int size = wcsp->getDomainSize(varIndex);
            ValueCost sorted[size];
            //ValueCost* sorted = new ValueCost [size];
            wcsp->iniSingleton();
            wcsp->getEnumDomainAndCost(varIndex, sorted);
            qsort(sorted, size, sizeof(ValueCost), cmpValueCost);
            for (int a = 0; a < size; a++) {
                deadend = false;
                int storedepth = Store::getDepth();
                try {
                    Store::store();
                    assign(varIndex, sorted[a].value);
                } catch (const Contradiction&) {
                    wcsp->whenContradiction();
                    deadend = true;
                    done = false;
                }
                Store::restore(storedepth);
                wcsp->updateSingleton();
                //cout << "(" << varIndex << "," << a <<  ")" << endl;
                if (deadend) {
                    remove(varIndex, sorted[a].value);
                    if (ToulBar2::verbose >= 0) {
                        cout << ".";
                        flush(cout);
                    }
                    // WARNING!!! can we stop if the variable is assigned, what about removeSingleton after???
                }
            }
            wcsp->removeSingleton();
            //delete [] sorted;
        }
    }
    if (ToulBar2::verbose >= 0)
        cout << "Done Singleton Consistency" << endl;
}

/*
 * Hybrid Depth-First and Best-First Branch and Bound
 *
 */

void Solver::newSolution()
{
#ifndef NDEBUG
    bool allVarsAssigned = true;
    for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
        if (wcsp->unassigned(i)) {
            allVarsAssigned = false;
            break;
        }
    }
    assert(allVarsAssigned);
#endif
    if (!ToulBar2::allSolutions && !ToulBar2::isZ) {
        ToulBar2::deltaUb = max(ToulBar2::deltaUbAbsolute, (Cost)(ToulBar2::deltaUbRelativeGap * (Double)wcsp->getLb()));
        wcsp->updateUb(wcsp->getLb());
    } else if (!ToulBar2::btdMode)
        nbSol += 1.;
    if (ToulBar2::isZ) {
        ToulBar2::logZ = wcsp->LogSumExp(ToulBar2::logZ, wcsp->getLb() + wcsp->getNegativeLb());
        if (ToulBar2::debug && (nbBacktracks % 10000LL) == 0 && ToulBar2::logepsilon > -numeric_limits<TLogProb>::infinity())
            cout << (ToulBar2::logZ + ToulBar2::markov_log) << " , " << (wcsp->LogSumExp(ToulBar2::logZ, ToulBar2::logU) + ToulBar2::markov_log) << endl;
    }
    if ((!ToulBar2::allSolutions && !ToulBar2::isZ) || ToulBar2::debug >= 2) {
        if (ToulBar2::verbose >= 0 || ToulBar2::showSolutions) {
            if (ToulBar2::haplotype)
                cout << "***New solution: " << wcsp->getLb() << " log10like: " << ToulBar2::haplotype->Cost2LogProb(wcsp->getLb()) / Log(10.) << " logProb: " << ToulBar2::haplotype->Cost2LogProb(wcsp->getLb()) << " (" << nbBacktracks << " backtracks, " << nbNodes << " nodes, depth " << Store::getDepth() << ")" << endl;
            else if (!ToulBar2::bayesian)
                cout << "New solution: " << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->getDDualBound() << std::setprecision(DECIMAL_POINT) << " (" << nbBacktracks << " backtracks, " << nbNodes << " nodes, depth " << Store::getDepth() << ")" << endl;
            else
                cout << "New solution: " << wcsp->getLb() << " energy: " << -(wcsp->Cost2LogProb(wcsp->getLb() + wcsp->getNegativeLb()) + ToulBar2::markov_log) << " prob: " << std::scientific << wcsp->Cost2Prob(wcsp->getLb() + wcsp->getNegativeLb()) * Exp(ToulBar2::markov_log) << std::fixed << " (" << nbBacktracks << " backtracks, " << nbNodes << " nodes, depth " << Store::getDepth() << ")" << endl;
        }
    }

    wcsp->restoreSolution();
    if (!ToulBar2::isZ)
        wcsp->setSolution(wcsp->getLb());

    if (ToulBar2::showSolutions) {

        if (ToulBar2::verbose >= 2)
            cout << *wcsp << endl;

        if (ToulBar2::allSolutions) {
            cout << std::setprecision(0) << nbSol << " solution(" << std::setprecision(ToulBar2::decimalPoint) << wcsp->getDDualBound() << "): ";
        }

        for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
            cout << " ";
            if (ToulBar2::pedigree) {
                cout << wcsp->getName(i) << ":";
                ToulBar2::pedigree->printGenotype(cout, wcsp->getValue(i));
            } else if (ToulBar2::haplotype) {
                ToulBar2::haplotype->printHaplotype(cout, wcsp->getValue(i), i);
            } else if (wcsp->enumerated(i) && ((EnumeratedVariable*)((WCSP*)wcsp)->getVar(i))->isValueNames()) {
                EnumeratedVariable* myvar = (EnumeratedVariable*)((WCSP*)wcsp)->getVar(i);
                Value myvalue = ((ToulBar2::sortDomains && ToulBar2::sortedDomains.find(i) != ToulBar2::sortedDomains.end()) ? ToulBar2::sortedDomains[i][myvar->toIndex(myvar->getValue())].value : myvar->getValue());
                string valuelabel = myvar->getValueName(myvar->toIndex(myvalue));
                string varlabel = myvar->getName();

                switch (ToulBar2::showSolutions) {
                case 1:
                    cout << myvalue;
                    break;
                case 2:
                    cout << valuelabel;
                    break;
                case 3:
                    cout << varlabel << "=" << valuelabel;
                    break;
                default:
                    break;
                }
            } else {
                cout << ((ToulBar2::sortDomains && ToulBar2::sortedDomains.find(i) != ToulBar2::sortedDomains.end()) ? ToulBar2::sortedDomains[i][wcsp->toIndex(i, wcsp->getValue(i))].value : wcsp->getValue(i));
            }
        }
        cout << endl;
        if (ToulBar2::bep)
            ToulBar2::bep->printSolution((WCSP*)wcsp);
    }
    if (ToulBar2::pedigree) {
        ToulBar2::pedigree->printCorrection((WCSP*)wcsp);
    }
    if (ToulBar2::writeSolution) {
        if (ToulBar2::pedigree) {
            string problemname = ToulBar2::problemsaved_filename;
            if (problemname.rfind(".wcsp") != string::npos)
                problemname.replace(problemname.rfind(".wcsp"), 5, ".pre");
            else if (problemname.rfind(".cfn") != string::npos)
                problemname.replace(problemname.rfind(".wcsp"), 4, ".pre");
            ToulBar2::pedigree->save((problemname.find("problem.pre") == 0) ? "problem_corrected.pre" : problemname.c_str(), (WCSP*)wcsp, true, false);
            ToulBar2::pedigree->printSol((WCSP*)wcsp);
            ToulBar2::pedigree->printCorrectSol((WCSP*)wcsp);
        } else if (ToulBar2::haplotype) {
            ToulBar2::haplotype->printSol((WCSP*)wcsp);
        }
        if (ToulBar2::solutionFile != NULL) {
            if (!ToulBar2::allSolutions)
                fseek(ToulBar2::solutionFile, ToulBar2::solutionFileRewindPos, SEEK_SET);
            wcsp->printSolution(ToulBar2::solutionFile);
            fprintf(ToulBar2::solutionFile, "\n");
        }
    }

    if (ToulBar2::xmlflag) {
        cout << "o " << wcsp->getLb() << endl; //" ";
    }
    if (ToulBar2::maxsateval) {
        cout << "o " << wcsp->getLb() << endl;
    }
    if (ToulBar2::uaieval && !ToulBar2::isZ) {
        ((WCSP*)wcsp)->solution_UAI(wcsp->getLb());
    }

    if (ToulBar2::newsolution)
        (*ToulBar2::newsolution)(wcsp->getIndex(), wcsp->getSolver());

    if (ToulBar2::restart == 0 && !ToulBar2::lds && !ToulBar2::isZ)
        throw NbBacktracksOut();
    if (ToulBar2::allSolutions && nbSol >= ToulBar2::allSolutions)
        throw NbSolutionsOut();
    if (ToulBar2::divNbSol > 1 && wcsp->getLb() <= prevDivSolutionCost)
        throw DivSolutionOut();
}

void Solver::recursiveSolve(Cost lb)
{
    int varIndex = -1;
    if (ToulBar2::bep)
        varIndex = getMostUrgent();
    else if (ToulBar2::Static_variable_ordering)
        varIndex = getNextUnassignedVar();
    else if (ToulBar2::weightedDegree && ToulBar2::lastConflict)
        varIndex = ((ToulBar2::restart > 0) ? getVarMinDomainDivMaxWeightedDegreeLastConflictRandomized() : getVarMinDomainDivMaxWeightedDegreeLastConflict());
    else if (ToulBar2::lastConflict)
        varIndex = ((ToulBar2::restart > 0) ? getVarMinDomainDivMaxDegreeLastConflictRandomized() : getVarMinDomainDivMaxDegreeLastConflict());
    else if (ToulBar2::weightedDegree)
        varIndex = ((ToulBar2::restart > 0) ? getVarMinDomainDivMaxWeightedDegreeRandomized() : getVarMinDomainDivMaxWeightedDegree());
    else
        varIndex = ((ToulBar2::restart > 0) ? getVarMinDomainDivMaxDegreeRandomized() : getVarMinDomainDivMaxDegree());
    if (varIndex >= 0) {
        *((StoreCost*)searchSize) += ((Cost)(10e6 * Log(wcsp->getDomainSize(varIndex))));
        if (ToulBar2::bep)
            scheduleOrPostpone(varIndex);
        else if (wcsp->enumerated(varIndex)) {
            if (ToulBar2::binaryBranching) {
                assert(wcsp->canbe(varIndex, wcsp->getSupport(varIndex)));
                // Reuse last solution found if available
                Value bestval = ((ToulBar2::verifyOpt) ? (wcsp->getSup(varIndex) + 1) : wcsp->getBestValue(varIndex));
                binaryChoicePoint(varIndex, (wcsp->canbe(varIndex, bestval)) ? bestval : wcsp->getSupport(varIndex), lb);
            } else
                narySortedChoicePoint(varIndex, lb);
        } else {
            return binaryChoicePoint(varIndex, wcsp->getInf(varIndex), lb);
        }
    } else {
        assert(ToulBar2::isZ || lb <= wcsp->getLb());
        newSolution();
    }
}

void Solver::recursiveSolveLDS(int discrepancy)
{
    int varIndex = -1;
    if (ToulBar2::bep)
        varIndex = getMostUrgent();
    else if (ToulBar2::weightedDegree && ToulBar2::lastConflict)
        varIndex = ((ToulBar2::restart > 0) ? getVarMinDomainDivMaxWeightedDegreeLastConflictRandomized() : getVarMinDomainDivMaxWeightedDegreeLastConflict());
    else if (ToulBar2::lastConflict)
        varIndex = ((ToulBar2::restart > 0) ? getVarMinDomainDivMaxDegreeLastConflictRandomized() : getVarMinDomainDivMaxDegreeLastConflict());
    else if (ToulBar2::weightedDegree)
        varIndex = ((ToulBar2::restart > 0) ? getVarMinDomainDivMaxWeightedDegreeRandomized() : getVarMinDomainDivMaxWeightedDegree());
    else
        varIndex = ((ToulBar2::restart > 0) ? getVarMinDomainDivMaxDegreeRandomized() : getVarMinDomainDivMaxDegree());
    if (varIndex >= 0) {
        if (ToulBar2::bep)
            scheduleOrPostpone(varIndex);
        else if (wcsp->enumerated(varIndex)) {
            if (ToulBar2::binaryBranching) {
                assert(wcsp->canbe(varIndex, wcsp->getSupport(varIndex)));
                // Reuse last solution found if available
                Value bestval = wcsp->getBestValue(varIndex);
                binaryChoicePointLDS(varIndex, (wcsp->canbe(varIndex, bestval)) ? bestval : wcsp->getSupport(varIndex), discrepancy);
            } else {
                narySortedChoicePointLDS(varIndex, discrepancy);
            }
        } else {
            binaryChoicePointLDS(varIndex, wcsp->getInf(varIndex), discrepancy);
        }
    } else {
        newSolution();
    }
}

pair<Cost, Cost> Solver::hybridSolve(Cluster* cluster, Cost clb, Cost cub)
{
    if (ToulBar2::verbose >= 1 && cluster)
        cout << "hybridSolve C" << cluster->getId() << " " << clb << " " << cub << endl;
    assert(clb < cub);
    assert(wcsp->getUb() == cub);
    assert(wcsp->getLb() <= clb);
    if (ToulBar2::hbfs) {
        CPStore* cp_ = NULL;
        OpenList* open_ = NULL;
        Cost delta = MIN_COST;
        if (cluster) {
            // BFS with BTD on current cluster (can be root or not)
            assert(cluster->cp);
            cp_ = cluster->cp;
            if (cluster == wcsp->getTreeDec()->getRoot()) {
                if (!cluster->open)
                    cluster->open = new OpenList();
                cluster->setUb(cub); // global problem upper bound
            } else {
                delta = cluster->getCurrentDelta();
                if (!cluster->open) {
                    cluster->nogoodRec(clb, MAX_COST, &cluster->open); // create an initial empty open list
                    cluster->setUb(MAX_COST); // no initial solution found for this cluster
                }
            }
            assert(cluster->open);
            open_ = cluster->open;
#ifndef NDEBUG
            OpenList* prevopen = cluster->open;
            Cost tmplb = MIN_COST;
            Cost tmpub = MAX_COST;
            assert(cluster == wcsp->getTreeDec()->getRoot() || cluster->nogoodGet(tmplb, tmpub, &cluster->open));
            assert(prevopen == cluster->open);
            assert(cluster == wcsp->getTreeDec()->getRoot() || tmpub == cluster->getUb());
            assert(cluster != wcsp->getTreeDec()->getRoot() || cub == cluster->getUb());
#endif
        } else {
            // normal BFS without BTD, i.e., hybridSolve is not reentrant
            if (cp != NULL)
                delete cp;
            cp = new CPStore();
            cp_ = cp;
            if (open != NULL)
                delete open;
            open = new OpenList();
            open_ = open;
        }
        cp_->store();
        if (open_->size() == 0 || (cluster && (clb >= open_->getClosedNodesLb(delta) || cub > open_->getUb(delta)))) { // start a new list of open nodes if needed
            if (open_->size() == 0 && (!cluster || cluster->getNbVars() > 0))
                nbHybridNew++;
            // reinitialize current open list and insert empty node
            *open_ = OpenList(MAX(MIN_COST, cub + delta), MAX(MIN_COST, cub + delta));
            addOpenNode(*cp_, *open_, clb, delta);
        } else if (!cluster || cluster->getNbVars() > 0)
            nbHybridContinue++;
        if (!cluster || cluster->getNbVars() > 0)
            nbHybrid++; // do not count empty root cluster
        if (cluster)
            cluster->hbfsGlobalLimit = ((ToulBar2::hbfsGlobalLimit > 0) ? (nbBacktracks + ToulBar2::hbfsGlobalLimit) : LONGLONG_MAX);
        Cost initiallb = clb;
        Cost initialub = cub;
        open_->updateUb(cub, delta);
        clb = MAX(clb, open_->getLb(delta));
        if (ToulBar2::verbose >= 1 && cluster)
            cout << "hybridSolve-2 C" << cluster->getId() << " " << clb << " " << cub << " " << delta << " " << open_->size() << " " << open_->top().getCost(delta) << " " << open_->getClosedNodesLb(delta) << " " << open_->getUb(delta) << endl;
        while (clb < cub && !open_->finished() && (!cluster || (clb == initiallb && cub == initialub && nbBacktracks <= cluster->hbfsGlobalLimit))) {
            if (cluster) {
                cluster->hbfsLimit = ((ToulBar2::hbfs > 0) ? (cluster->nbBacktracks + ToulBar2::hbfs) : LONGLONG_MAX);
                assert(wcsp->getTreeDec()->getCurrentCluster() == cluster);
                wcsp->setUb(cub);
                assert(cluster->isActive());
                assert(cluster->getLbRec() == wcsp->getLb());
            } else
                hbfsLimit = ((ToulBar2::hbfs > 0) ? (nbBacktracks + ToulBar2::hbfs) : LONGLONG_MAX);
            int storedepthBFS = Store::getDepth();
            int storeVAC = ToulBar2::vac;
            try {
                Store::store();
                OpenNode nd = open_->top();
                open_->pop();
                if (ToulBar2::verbose >= 3) {
                    if (wcsp->getTreeDec())
                        cout << "[C" << wcsp->getTreeDec()->getCurrentCluster()->getId() << "] ";
                    cout << "[ " << nd.getCost(delta) << ", " << cub << "] ( " << open_->size() << "+1 still open)" << endl;
                }
                if (ToulBar2::vac < 0 && Store::getDepth() + (nd.last - nd.first) >= abs(ToulBar2::vac))
                    ToulBar2::vac = 0;
                restore(*cp_, nd);
                Cost bestlb = MAX(nd.getCost(delta), wcsp->getLb());
                bestlb = MAX(bestlb, clb);
                if (cluster) {
                    pair<Cost, Cost> res = recursiveSolve(cluster, bestlb, cub);
                    assert(res.first <= res.second);
                    assert(res.first >= bestlb);
                    assert(res.second <= cub);
                    assert(res.second == cub || cluster->getUb() == res.second);
                    assert(open_->empty() || open_->top().getCost(delta) >= nd.getCost(delta));
                    open_->updateClosedNodesLb(res.first, delta);
                    open_->updateUb(res.second, delta);
                    cub = MIN(cub, res.second);
                } else {
                    if (ToulBar2::vac < 0)
                        ToulBar2::vac = 0;
                    recursiveSolve(bestlb);
                }
            } catch (const Contradiction&) {
                wcsp->whenContradiction();
            }
            if (!cluster) { // synchronize current upper bound with DFS (without tree decomposition)
                cub = wcsp->getUb();
                open_->updateUb(cub);
            }
            Store::restore(storedepthBFS);
            ToulBar2::vac = storeVAC;
            cp_->store();
            if (cp_->size() >= static_cast<std::size_t>(ToulBar2::hbfsCPLimit) || open_->size() >= static_cast<std::size_t>(ToulBar2::hbfsOpenNodeLimit)) {
                ToulBar2::hbfs = 0;
                ToulBar2::hbfsGlobalLimit = 0;
                if (cluster) {
                    cluster->hbfsGlobalLimit = LONGLONG_MAX;
                    cluster->hbfsLimit = LONGLONG_MAX;
                } else
                    hbfsLimit = LONGLONG_MAX;
            }
            clb = MAX(clb, open_->getLb(delta));
            showGap(clb, cub);
            if (ToulBar2::hbfs && nbRecomputationNodes > 0) { // wait until a nonempty open node is restored (at least after first global solution is found)
                assert(nbNodes > 0);
                if (nbRecomputationNodes > nbNodes / ToulBar2::hbfsBeta && ToulBar2::hbfs <= ToulBar2::hbfsGlobalLimit)
                    ToulBar2::hbfs *= 2;
                else if (nbRecomputationNodes < nbNodes / ToulBar2::hbfsAlpha && ToulBar2::hbfs >= 2)
                    ToulBar2::hbfs /= 2;
                if (ToulBar2::debug >= 2)
                    cout << "HBFS backtrack limit: " << ToulBar2::hbfs << endl;
            }
        }
        assert(clb >= initiallb && cub <= initialub);
    } else {
        if (cluster) {
            cluster->hbfsGlobalLimit = LONGLONG_MAX;
            cluster->hbfsLimit = LONGLONG_MAX;
            pair<Cost, Cost> res = recursiveSolve(cluster, clb, cub);
            clb = MAX(clb, res.first);
            cub = MIN(cub, res.second);
        } else {
            hbfsLimit = LONGLONG_MAX;
            recursiveSolve();
            cub = wcsp->getUb();
            clb = cub;
        }
    }
    assert(clb <= cub);
    return make_pair(clb, cub);
}

void Solver::beginSolve(Cost ub)
{
    // Last-minute compatibility checks for ToulBar2 selected options
    if (ub <= MIN_COST) {
        cerr << "Error: wrong initial primal bound (negative or zero)." << endl;
        exit(1);
    }
    if (ToulBar2::allSolutions && ToulBar2::btdMode == 1 && ub > 1) {
        cerr << "Error: Solution enumeration by BTD-like search methods is only possible for feasability (use -ub=1 and integer costs only)." << endl;
        exit(1);
    }
    if (ToulBar2::allSolutions && ToulBar2::btdMode == 1 && ub == 1 && ToulBar2::hbfs) {
        cerr << "Error: Hybrid best-first search cannot currently look for all solutions when BTD mode is activated. Shift to DFS (use -hbfs:)." << endl;
        exit(1);
    }
    if (ToulBar2::FullEAC && ToulBar2::vac > 1 && wcsp->numberOfConnectedConstraints() > wcsp->numberOfConnectedBinaryConstraints()) {
        cerr << "Warning: VAC during search and Full EAC variable ordering heuristic not implemented with non binary cost functions (remove -vacint option)." << endl;
        exit(1);
    }

    if (ToulBar2::searchMethod != DFBB) {
        if (!ToulBar2::lds || ToulBar2::vnsLDSmax < 0)
            ToulBar2::vnsLDSmax = wcsp->getDomainSizeSum() - wcsp->numberOfUnassignedVariables();
        if (!ToulBar2::lds)
            ToulBar2::vnsLDSmin = wcsp->getDomainSizeSum() - wcsp->numberOfUnassignedVariables();
        if (ToulBar2::vnsKmax <= 0)
            ToulBar2::vnsKmax = wcsp->numberOfUnassignedVariables();
    }
    if (wcsp->isGlobal() && ToulBar2::btdMode >= 1) {
        cout << "Error: cannot use BTD-like search methods with monolithic global cost functions (remove -B option)." << endl;
        exit(1);
    }
    if (wcsp->isGlobal() && (ToulBar2::elimDegree_preprocessing >= 1 || ToulBar2::elimDegree_preprocessing < -1)) {
        cout << "Warning! Cannot use generic variable elimination with global cost functions." << endl;
        ToulBar2::elimDegree_preprocessing = -1;
    }
    if (ToulBar2::incop_cmd.size() > 0) {
        for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
            if (wcsp->unassigned(i) && !wcsp->enumerated(i)) {
                cout << "Warning! Cannot use INCOP local search with bounds arc propagation (non enumerated variable domains)." << endl;
                ToulBar2::incop_cmd = "";
                break;
            }
        }
    }
    if (((WCSP*)wcsp)->isAlreadyTreeDec(ToulBar2::varOrder)) {
        if (ToulBar2::btdMode >= 3) {
            cout << "Warning! Cannot apply path decomposition with a given tree decomposition file." << endl;
            ToulBar2::btdMode = 2;
        }
        if (ToulBar2::btdMode >= 1) {
            if (ToulBar2::elimDegree_preprocessing > 0) {
                cout << "Warning! Cannot apply variable elimination in preprocessing with a given tree decomposition file." << endl;
                ToulBar2::elimDegree_preprocessing = 0;
            }
            if (ToulBar2::elimDegree > 0) {
                cout << "Warning! Cannot apply variable elimination during search with a given tree decomposition file." << endl;
                ToulBar2::elimDegree = 0;
            }
            if (ToulBar2::preprocessFunctional > 0) {
                cout << "Warning! Cannot apply functional variable elimination with a given tree decomposition file." << endl;
                ToulBar2::preprocessFunctional = 0;
            }
        }
    }

    nbBacktracks = 0;
    nbBacktracksLimit = ToulBar2::backtrackLimit;
    nbNodes = 0;
    nbRecomputationNodes = 0;
    lastConflictVar = -1;
    tailleSep = 0;
    ToulBar2::limited = false;

    if (ToulBar2::DEE)
        ToulBar2::DEE_ = ToulBar2::DEE; // enforces PSNS after closing the model

    if (CSP(wcsp->getLb(), wcsp->getUb())) {
        ToulBar2::LcLevel = LC_AC;
    }

    if (ToulBar2::isZ) {
        ToulBar2::logZ = -numeric_limits<TLogProb>::infinity();
        ToulBar2::logU = -numeric_limits<TLogProb>::infinity();
    }

    // reactivate on-the-fly variable elimination and dead-end elimination if needed
    for (int i = wcsp->numberOfVariables() - 1; i >= 0; i--) {
        if (wcsp->unassigned(i)) {
            ((WCSP *)wcsp)->getVar(i)->queueEliminate();
            ((WCSP *)wcsp)->getVar(i)->queueDEE();
        }
    }
}

Cost Solver::preprocessing(Cost initialUpperBound)
{
    Long hbfs_ = ToulBar2::hbfs;
    ToulBar2::hbfs = 0; // do not perform hbfs operations in preprocessing except for building tree decomposition
    if (!ToulBar2::isZ) {
        Cost finiteUb = wcsp->finiteUb(); // find worst-case assignment finite cost plus one as new upper bound
        if (finiteUb < initialUpperBound) {
            initialUpperBound = finiteUb;
            ToulBar2::deltaUb = max(ToulBar2::deltaUbAbsolute, (Cost)(ToulBar2::deltaUbRelativeGap * (Double)min(finiteUb, wcsp->getUb())));
            wcsp->updateUb(finiteUb + ToulBar2::deltaUb);
        }
        wcsp->setInfiniteCost(); // shrink forbidden costs based on problem lower and upper bounds to avoid integer overflow errors when summing costs
    }
    Cost initialLowerBound = wcsp->getLb();
    wcsp->enforceUb();
    wcsp->propagate(); // initial propagation
    if (!ToulBar2::isZ) {
        Cost finiteUb = wcsp->finiteUb(); // find worst-case assignment finite cost plus one as new upper bound
        if (finiteUb < initialUpperBound || wcsp->getLb() > initialLowerBound) {
            if (finiteUb < initialUpperBound) {
                ToulBar2::deltaUb = max(ToulBar2::deltaUbAbsolute, (Cost)(ToulBar2::deltaUbRelativeGap * (Double)min(finiteUb, wcsp->getUb())));
                wcsp->updateUb(finiteUb + ToulBar2::deltaUb);
            }
            wcsp->setInfiniteCost();
            if (finiteUb < initialUpperBound) {
                wcsp->enforceUb();
                wcsp->propagate();
                initialUpperBound = finiteUb;
            }
        }
    }
    wcsp->preprocessing(); // preprocessing after initial propagation
    if (!ToulBar2::isZ) {
        Cost finiteUb = wcsp->finiteUb(); // find worst-case assignment finite cost plus one as new upper bound
        if (finiteUb < initialUpperBound) {
            ToulBar2::deltaUb = max(ToulBar2::deltaUbAbsolute, (Cost)(ToulBar2::deltaUbRelativeGap * (Double)min(finiteUb, wcsp->getUb())));
            wcsp->updateUb(finiteUb + ToulBar2::deltaUb);
        }
        wcsp->setInfiniteCost();
        if (finiteUb < initialUpperBound) {
            wcsp->enforceUb();
            wcsp->propagate();
            initialUpperBound = finiteUb;
        }
    }
    if (ToulBar2::verbose >= 0)
        cout << "Preprocessing time: " << cpuTime() - ToulBar2::startCpuTime << " seconds." << endl;

    // special data structure to be initialized for variable ordering heuristics
    initVarHeuristic();

    int lds = ToulBar2::lds;
    ToulBar2::lds = 0; // avoid TimeOut exception when new solutions found
    if (ToulBar2::incop_cmd.size() > 0) {
        double incopStartTime = cpuTime();
        vector<int> bestsol(getWCSP()->numberOfVariables(), 0);
        for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++)
            bestsol[i] = (wcsp->canbe(i, wcsp->getBestValue(i)) ? wcsp->getBestValue(i) : wcsp->getSupport(i));
        narycsp(ToulBar2::incop_cmd, bestsol);
        if (ToulBar2::verbose >= 0)
            cout << "INCOP solving time: " << cpuTime() - incopStartTime << " seconds." << endl;
    }
    ToulBar2::lds = lds;

    if (ToulBar2::singletonConsistency) {
        singletonConsistency();
        wcsp->propagate();
    }

    ToulBar2::hbfs = hbfs_; // do not perform hbfs operations in preprocessing except for building tree decomposition

    if (ToulBar2::verbose >= 0)
        cout << wcsp->numberOfUnassignedVariables() << " unassigned variables, " << wcsp->getDomainSizeSum() << " values in all current domains (med. size:" << wcsp->medianDomainSize() << ", max size:" << wcsp->getMaxCurrentDomainSize() << ") and " << wcsp->numberOfConnectedConstraints() << " non-unary cost functions (med. arity:" << wcsp->medianArity() << ", med. degree:" << wcsp->medianDegree() << ")" << endl;
    if (ToulBar2::verbose >= 0) {
        Double Dlb = wcsp->getDLb();
        Double Dub = wcsp->getDUb();
        cout << "Initial lower and upper bounds: [" << std::fixed << std::setprecision(ToulBar2::decimalPoint) << Dlb << ", " << Dub << "] " << std::setprecision(DECIMAL_POINT) << (100.0 * (Dub - Dlb)) / max(fabsl(Dlb), fabsl(Dub)) << "%" << endl;
    }
    initGap(wcsp->getLb(), wcsp->getUb());

    if (ToulBar2::DEE == 4) {
        ToulBar2::DEE_ = 0; // only PSNS in preprocessing
        ToulBar2::DEE = 0; // avoid doing DEE later in the case of incremental solving
    }

    if (ToulBar2::isZ && ToulBar2::verbose >= 1)
        cout << "NegativeShiftingCost= " << wcsp->getNegativeLb() << endl;

    if (ToulBar2::btdMode) {
        if (wcsp->numberOfUnassignedVariables() == 0 || wcsp->numberOfConnectedConstraints() == 0)
            ToulBar2::approximateCountingBTD = 0;
        ToulBar2::vac = 0; // VAC is not compatible with restricted tree decomposition propagation
        wcsp->buildTreeDecomposition();
    } else if (ToulBar2::weightedDegree && (((Long)wcsp->numberOfConnectedConstraints()) >= ((Long)ToulBar2::weightedDegree))) {
        if (ToulBar2::verbose >= 0)
            cout << "Weighted degree heuristic disabled (#costfunctions=" << wcsp->numberOfConnectedConstraints() << " >= " << ToulBar2::weightedDegree << ")" << endl;
        ToulBar2::weightedDegree = 0;
    }

    if (ToulBar2::dumpWCSP) {
        dump_wcsp(ToulBar2::problemsaved_filename.c_str(), false, static_cast<ProblemFormat>((ToulBar2::dumpWCSP >> 1)+(ToulBar2::dumpWCSP & 1)));
        cout << "end." << endl;
        exit(0);
    }

    return initialUpperBound;
}

bool Solver::solve(bool first)
{
    beginSolve(wcsp->getUb());

    Cost initialUpperBound = wcsp->getUb();

    //        Store::store();       // if uncomment then solve() does not change the problem but all preprocessing operations will allocate in backtrackable memory
    int initdepth = Store::getDepth();
    try {
        try {
            if (first) {
                initialUpperBound = preprocessing(initialUpperBound);
            } else {
                if (ToulBar2::elimDegree >= 0)
                    ToulBar2::elimDegree_ = ToulBar2::elimDegree;
                initGap(wcsp->getLb(), wcsp->getUb());
            }

            Cost upperbound = MAX_COST;
            if (ToulBar2::restart >= 0) {
                if (ToulBar2::restart > 0)
                    nbBacktracksLimit = 1;
                upperbound = wcsp->getUb();
            }
            bool nbbacktracksout = false;
            int nbrestart = 0;
            Long currentNbBacktracksLimit = 1;
            Long nbBacktracksLimitTop = 1;
            int storedepth = Store::getDepth();
            do {
                Store::store();
                if (ToulBar2::restart >= 0) {
                    nbbacktracksout = false;
                    nbrestart++;
                    // currentNbBacktracksLimit = max(currentNbBacktracksLimit + 1, (Long) (1.2 * (Double) currentNbBacktracksLimit + 0.5));
                    // if (ToulBar2::lds) currentNbBacktracksLimit *= 4;
                    currentNbBacktracksLimit = luby(nbrestart);
                    if (currentNbBacktracksLimit > nbBacktracksLimitTop || (wcsp->getUb() < upperbound)) {
                        nbBacktracksLimitTop = currentNbBacktracksLimit;
                        currentNbBacktracksLimit = 1;
                    }
                    //			if (!(wcsp->getUb() < upperbound) && nbNodes >= ToulBar2::restart) {
                    if (nbNodes >= ToulBar2::restart) {
                        nbBacktracksLimit = ToulBar2::backtrackLimit;
                        ToulBar2::restart = 0;
                        if (ToulBar2::verbose >= 0)
                            cout << "****** Restart " << nbrestart << " with no backtrack limit and UB=" << wcsp->getUb() << " ****** (" << nbNodes << " nodes)" << endl;
                        if (ToulBar2::debug >= 1 && ToulBar2::weightedDegree > 0) {
                            int size = unassignedVars->getSize();
                            ValueCost sorted[size];
                            int i = 0;
                            for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
                                sorted[i].value = *iter;
                                sorted[i].cost = wcsp->getWeightedDegree(*iter);
                                i++;
                            }
                            qsort(sorted, size, sizeof(ValueCost), cmpValueCost);
                            for (i = 0; i < size; i++) {
                                cout << wcsp->getName(sorted[i].value) << " " << wcsp->getDomainSize(sorted[i].value) << " " << sorted[i].cost << endl;
                            }
                        }
                    } else {
                        nbBacktracksLimit = min(ToulBar2::backtrackLimit, nbBacktracks + currentNbBacktracksLimit * 100);
                        if (ToulBar2::verbose >= 0)
                            cout << "****** Restart " << nbrestart << " with " << currentNbBacktracksLimit * 100 << " backtracks max and UB=" << wcsp->getUb() << " ****** (" << nbNodes << " nodes)" << endl;
                    }
                    upperbound = wcsp->getUb();
                    enforceUb();
                    wcsp->propagate();
                    Store::store();
                    if (ToulBar2::isZ) {
                        ToulBar2::logZ = -numeric_limits<TLogProb>::infinity();
                        ToulBar2::logU = -numeric_limits<TLogProb>::infinity();
                    }
                }
                try {
                    if (ToulBar2::restart <= 0 && ToulBar2::lds) {
                        int discrepancy = 0;
                        do {
                            if (discrepancy > abs(ToulBar2::lds)) {
                                if (ToulBar2::verbose >= 0)
                                    cout << "--- [" << Store::getDepth() << "] Search with no discrepancy limit --- (" << nbNodes << " nodes)" << endl;
                            } else {
                                if (ToulBar2::verbose >= 0)
                                    cout << "--- [" << Store::getDepth() << "] LDS " << discrepancy << " --- (" << nbNodes << " nodes)" << endl;
                            }
                            ToulBar2::limited = false;
                            enforceUb();
                            wcsp->propagate();
                            if (ToulBar2::isZ) {
                                ToulBar2::logZ = -numeric_limits<TLogProb>::infinity();
                                ToulBar2::logU = -numeric_limits<TLogProb>::infinity();
                            }
                            if (discrepancy > abs(ToulBar2::lds)) {
                                if (ToulBar2::lds < 0) {
                                    ToulBar2::limited = true;
                                    THROWCONTRADICTION;
                                }
                                ToulBar2::lds = 0;
                                initialDepth = Store::getDepth();
                                hybridSolve();
                            } else {
                                int storedepth = Store::getDepth();
                                try {
                                    Store::store();
                                    initialDepth = Store::getDepth();
                                    recursiveSolveLDS(discrepancy);
                                } catch (const Contradiction&) {
                                    wcsp->whenContradiction();
                                }
                                Store::restore(storedepth);
                            }
                            if (discrepancy > 0)
                                discrepancy *= 2;
                            else
                                discrepancy++;
                        } while (ToulBar2::limited);
                    } else {
                        TreeDecomposition* td = wcsp->getTreeDec();
                        if (td) {
                            Cost ub = wcsp->getUb();
                            Cluster* start = td->getRoot();
                            assert(start->getLbRec() == MIN_COST); // local lower bounds (and delta costs) must be zero!
                            if (ToulBar2::btdSubTree >= 0)
                                start = td->getCluster(ToulBar2::btdSubTree);
                            td->setCurrentCluster(start);
                            if (start == td->getRoot())
                                start->setLb(wcsp->getLb()); // initial lower bound found by propagation is associated to tree decompostion root cluster
                            switch (ToulBar2::btdMode) {
                            case 0:
                            case 1: {
                                if (ToulBar2::allSolutions) {
                                    timeDeconnect = 0.;
                                    BigInteger cartesianProduct = 1;
                                    nbSol = (wcsp->numberOfConnectedConstraints() == 0) ? (wcsp->cartProd(cartesianProduct), cartesianProduct) : sharpBTD(start);
                                    if (ToulBar2::approximateCountingBTD && nbSol > 0. && td->getRoot()->getNbVars() == 0) { //if there are several parts
                                        approximate(nbSol, td);
                                    }
                                    // computation of maximal separator size
                                    for (int i = 0; i < td->getNbOfClusters(); i++) {
                                        if (td->getCluster(i)->sepSize() > tailleSep)
                                            tailleSep = td->getCluster(i)->sepSize();
                                    }
                                } else {
                                    pair<Cost, Cost> res = make_pair(wcsp->getLb(), ub);
                                    do {
                                        int storedepth = Store::getDepth();
                                        try {
                                            Store::store();
                                            td->setCurrentCluster(start);
                                            enforceUb();
                                            wcsp->propagate();
                                            initialDepth = Store::getDepth();
                                            res = hybridSolve(start, MAX(wcsp->getLb(), res.first), res.second);
                                            //				                if (res.first < res.second) cout << "Optimality gap: [ " <<  res.first << " , " << res.second << " ] " << (100. * (res.second-res.first)) / res.second << " % (" << nbBacktracks << " backtracks, " << nbNodes << " nodes)" << endl;
                                        } catch (const Contradiction&) {
                                            wcsp->whenContradiction();
                                            res.first = res.second;
                                        }
                                        Store::restore(storedepth);
                                        ub = res.second;
                                        wcsp->setUb(ub);
                                    } while (res.first < res.second);
                                    assert(res.first == res.second);
                                }
                                break;
                            }
                            case 2:
                            case 3: {
                                pair<Cost, Cost> res = make_pair(wcsp->getLb(), ub);
                                do { //TODO: set up for optimality gap pretty print
                                    res = russianDollSearch(start, res.second);
                                    //				        if (res.first < res.second) cout << "Optimality gap: [ " <<  res.first << " , " << res.second << " ] " << (100. * (res.second-res.first)) / res.second << " % (" << nbBacktracks << " backtracks, " << nbNodes << " nodes)" << endl;
                                } while (res.first < res.second);
                                assert(res.first == res.second);
                                ub = start->getLbRDS();
                                assert(ub == res.second);
                                wcsp->setUb(ub);
                                break;
                            }
                            default: {
                                cerr << "Unknown search method B" << ToulBar2::btdMode << endl;
                                exit(EXIT_FAILURE);
                            }
                            }
                            if (ToulBar2::debug)
                                start->printStatsRec();
                            if (ToulBar2::verbose >= 0 && nbHybrid >= 1)
                                cout << "HBFS open list restarts: " << (100. * (nbHybrid - nbHybridNew - nbHybridContinue) / nbHybrid) << " % and reuse: " << (100. * nbHybridContinue / nbHybrid) << " % of " << nbHybrid << endl;
                        } else {
                            if (ToulBar2::useRASPS) {
                                enforceUb();
                                wcsp->propagate();
                                ToulBar2::RASPS = true;
                                ((WCSP*)wcsp)->vac->iniThreshold(ToulBar2::RASPSlastitThreshold);
                                ((WCSP*)wcsp)->vac->propagate(); // VAC done again
                                ToulBar2::RASPS = false;
                                if (ToulBar2::verbose >= 0)
                                    cout << "RASPS done in preprocessing (backtrack: " << nbBacktracks << " nodes: " << nbNodes << ")" << endl;
                                enforceUb();
                                wcsp->propagate();
                                if (ToulBar2::RASPSreset) {
                                    for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
                                        wcsp->resetWeightedDegree(i);
                                        wcsp->setBestValue(i, wcsp->getSup(i) + 1);
                                    }
                                }
                            }
                            if (ToulBar2::divNbSol > 1) {
                                int initDepth = Store::getDepth();
                                Cost initUb = initialUpperBound;
                                prevDivSolutionCost = MIN_COST;
                                bool incrementalSearch = true;
                                solTrie.init(wcsp->getDivVariables());
                                vector<Cost> energies;
                                bool extrapolatedBound = false;

                                try {
                                    do {
                                        if (ToulBar2::verbose >= 0)
                                            cout << "+++++++++ Search for solution " << energies.size() + 1 << " +++++++++" << endl;
                                        wcsp->setUb(initUb); // (re-)start search with an initial upper bound
                                        wcsp->initSolutionCost(); // set solution cost to infinity but do not forget solution vector

                                        //get solution from previous solve and add pairwise Hamming distance constraint
                                        if (energies.size() > 0 && !extrapolatedBound) {
                                            switch (ToulBar2::divMethod) {
                                            case 0:
                                                wcsp->addDivConstraint(wcsp->getSolution(), energies.size() - 1, initUb);
                                                break;
                                            case 1:
                                                wcsp->addHDivConstraint(wcsp->getSolution(), energies.size() - 1, initUb);
                                                break;
                                            case 2:
                                                wcsp->addTDivConstraint(wcsp->getSolution(), energies.size() - 1, initUb);
                                                break;
                                            default:
                                                cerr << "Error: no such diversity encoding method: " << ToulBar2::divMethod << endl;
                                                exit(EXIT_FAILURE);
                                            }
                                            wcsp->propagate();
                                        }

                                        Cost eUpperBound = initUb;
                                        if (!extrapolatedBound) { // previous extrapolated bound was fine
                                            if (energies.size() >= 2) {
                                                auto back = energies.end() - 1;
                                                Cost maxDelta = *(back) - *(--back);
                                                int maxCount = 5;
                                                while (back != energies.begin() && --maxCount >= 0) {
                                                    maxDelta = max(maxDelta, *(back) - *(--back));
                                                }
                                                Cost newUb = energies.back() + max(UNIT_COST, 2 * maxDelta);
                                                if (initUb > newUb) {
                                                    extrapolatedBound = true;
                                                    eUpperBound = newUb;
                                                    wcsp->setUb(eUpperBound); // start search with an extrapolated upper bound
                                                    if (ToulBar2::verbose >= 0)
                                                        cout << "+++++++++ predictive bounding: " << wcsp->Cost2ADCost(eUpperBound) << endl;
                                                }
                                            }
                                        } else {
                                            extrapolatedBound = false;
                                        }

                                        Store::store(); // protect the current CFN from changes by search or new cost functions
                                        try {
                                            try {
                                                if (ToulBar2::divWidth > 0 && energies.size() > 1) {
                                                    if (ToulBar2::verbose >= 1)
                                                        cout << "computing MDD.." << endl;
                                                    Mdd mdd = computeMDD(&solTrie, initUb);
                                                    if (ToulBar2::verbose >= 1)
                                                        cout << "MDD computed." << endl;
                                                    //ofstream os(to_string(this) + "-wregular.dot");
                                                    //printLayers(os, mdd);
                                                    //os.close();
                                                    switch (ToulBar2::divMethod) {
                                                    case 0:
                                                        wcsp->addMDDConstraint(mdd, ToulBar2::divNbSol - 1); //ToulBar2::divNbSol = index of the relaxed constraint
                                                        break;
                                                    case 1:
                                                        wcsp->addHMDDConstraint(mdd, ToulBar2::divNbSol - 1);
                                                        break;
                                                    case 2:
                                                        wcsp->addTMDDConstraint(mdd, ToulBar2::divNbSol - 1);
                                                        break;
                                                    default:
                                                        cerr << "Error: no such diversity encoding method: " << ToulBar2::divMethod << endl;
                                                        exit(EXIT_FAILURE);
                                                    }
                                                }
                                                // reactivate on-the-fly variable elimination and dead-end elimination
                                                for (int i = wcsp->numberOfVariables() - 1; i >= 0; i--) {
                                                    if (wcsp->unassigned(i)) {
                                                        ((WCSP *)wcsp)->getVar(i)->queueEliminate();
                                                        ((WCSP *)wcsp)->getVar(i)->queueDEE();
                                                    }
                                                }
                                                wcsp->propagate();
                                                initialDepth = Store::getDepth();
                                                hybridSolve(); // do not give prevDivSolutionCost as initial lower bound because it will generate too many open nodes with the same lower bound
                                            } catch (const DivSolutionOut&) {
                                                ToulBar2::limited = false;
                                            }
                                        } catch (const Contradiction&) {
                                            wcsp->whenContradiction();
                                        }
                                        Store::restore(initDepth); // undo search
                                        if (wcsp->getSolutionCost() < eUpperBound) {
                                            assert(wcsp->getSolutionCost() >= prevDivSolutionCost);
                                            prevDivSolutionCost = wcsp->getSolutionCost();
                                            energies.push_back(prevDivSolutionCost);
                                            vector<Value> divSol;
                                            for (auto const& var : wcsp->getDivVariables()) {
                                                divSol.push_back(wcsp->getSolution()[var->wcspIndex]);
                                            }
                                            solTrie.insertSolution(divSol);
                                            if (ToulBar2::solutionFile) {
                                                ToulBar2::solutionFileRewindPos = ftell(ToulBar2::solutionFile);
                                            }
                                            incrementalSearch = (energies.size() < ToulBar2::divNbSol && wcsp->getDivVariables().size() > 0 && ToulBar2::divBound <= wcsp->getDivVariables().size());
                                            extrapolatedBound = false;
                                        } else { // No solution found
                                            if (!extrapolatedBound) {
                                                incrementalSearch = false;
                                            }
                                        }
                                        endSolve(wcsp->getSolutionCost() < eUpperBound, wcsp->getSolutionCost(), !ToulBar2::limited);
                                    } while (incrementalSearch); // this or an exception (no solution)
                                } catch (const Contradiction&) {
                                    wcsp->whenContradiction();
                                    endSolve(wcsp->getSolutionCost() < initUb, wcsp->getSolutionCost(), !ToulBar2::limited);
                                }
                            } else {
                                initialDepth = Store::getDepth();
                                hybridSolve();
                            }
                        }
                    }
                } catch (const NbBacktracksOut&) {
                    if (nbBacktracks > ToulBar2::backtrackLimit)
                        throw NbBacktracksOut();
                    nbbacktracksout = true;
                    ToulBar2::limited = false;
                }
                Store::restore(storedepth);
            } while (nbbacktracksout);
        } catch (const Contradiction&) {
            wcsp->whenContradiction();
        }
    } catch (const SolverOut&) {
    }
    Store::restore(initdepth);
    //  Store::restore();         // see above for Store::store()

    if (ToulBar2::divNbSol <= 1)
        endSolve(wcsp->getSolutionCost() < initialUpperBound, wcsp->getSolutionCost(), !ToulBar2::limited);

    return (ToulBar2::isZ || ToulBar2::allSolutions || wcsp->getSolutionCost() < initialUpperBound);
}

void Solver::endSolve(bool isSolution, Cost cost, bool isComplete)
{
    ToulBar2::DEE_ = 0;
    ToulBar2::elimDegree_ = -1;

    static string solType[4] = { "Optimum: ", "Primal bound: ", "guaranteed primal bound: ", "Primal bound: " };

    int isLimited = (!isComplete) | ((ToulBar2::deltaUb != MIN_COST) << 1);

    if (ToulBar2::isZ) {
        if (ToulBar2::verbose >= 1)
            cout << "NegativeShiftingCost= " << wcsp->getNegativeLb() << endl;
        if (ToulBar2::uaieval) {
            rewind(ToulBar2::solution_uai_file);
            fprintf(ToulBar2::solution_uai_file, "PR\n");
            fprintf(ToulBar2::solution_uai_file, PrintFormatProb, (wcsp->LogSumExp(ToulBar2::logZ, ToulBar2::logU) + ToulBar2::markov_log) / Log(10.));
            fprintf(ToulBar2::solution_uai_file, "\n");
        }
        cout << (ToulBar2::logZ + ToulBar2::markov_log) << " <= Log(Z) <= ";
        cout << (wcsp->LogSumExp(ToulBar2::logZ, ToulBar2::logU) + ToulBar2::markov_log) << " in " << nbBacktracks << " backtracks and " << nbNodes << " nodes and " << cpuTime() - ToulBar2::startCpuTime << " seconds" << endl;
        cout << (ToulBar2::logZ + ToulBar2::markov_log) / Log(10.) << " <= Log10(Z) <= ";
        cout << (wcsp->LogSumExp(ToulBar2::logZ, ToulBar2::logU) + ToulBar2::markov_log) / Log(10.) << " in " << nbBacktracks << " backtracks and " << nbNodes << " nodes and " << cpuTime() - ToulBar2::startCpuTime << " seconds" << endl;
        return;
    }
    if (ToulBar2::allSolutions) {
        if (ToulBar2::approximateCountingBTD)
            cout << "Number of solutions    : ~= " << std::fixed << std::setprecision(0) << nbSol << std::setprecision(DECIMAL_POINT) << endl;
        else {
            if (!isComplete)
                cout << "Number of solutions    : >=  " << std::fixed << std::setprecision(0) << nbSol << std::setprecision(DECIMAL_POINT) << endl;
            else
                cout << "Number of solutions    : =  " << std::fixed << std::setprecision(0) << nbSol << std::setprecision(DECIMAL_POINT) << endl;
        }
        if (ToulBar2::btdMode >= 1) {
            cout << "Number of #goods       :    " << nbSGoods << endl;
            cout << "Number of used #goods  :    " << nbSGoodsUse << endl;
            cout << "Size of sep            :    " << tailleSep << endl;
        }
        cout << "Time                   :    " << cpuTime() - ToulBar2::startCpuTime << " seconds" << endl;
        cout << "... in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << ((ToulBar2::DEE) ? (" ( " + to_string(wcsp->getNbDEE()) + " removals by DEE)") : "") << endl;
        return;
    }

    if (ToulBar2::vac)
        wcsp->printVACStat();

    if (ToulBar2::verbose >= 0 && nbHybrid >= 1 && nbNodes > 0)
        cout << "Node redundancy during HBFS: " << 100. * nbRecomputationNodes / nbNodes << " %" << endl;

    if (isSolution) {
        if (ToulBar2::verbose >= 0 && !ToulBar2::uai && !ToulBar2::xmlflag && !ToulBar2::maxsateval) {
            if (ToulBar2::haplotype)
                cout << endl;

            if (isLimited == 2)
                cout << "(" << ToulBar2::deltaUbS << "," << std::scientific << ToulBar2::deltaUbRelativeGap << std::fixed << ")-";
            if (ToulBar2::haplotype)
                cout << solType[isLimited] << cost << " log10like: " << ToulBar2::haplotype->Cost2LogProb(cost) / Log(10.) << " loglike: " << ToulBar2::haplotype->Cost2LogProb(cost) << " in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << ((ToulBar2::DEE) ? (" ( " + to_string(wcsp->getNbDEE()) + " removals by DEE)") : "") << " and " << cpuTime() - ToulBar2::startCpuTime << " seconds." << endl;
            else if (!ToulBar2::bayesian)
                cout << solType[isLimited] << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->Cost2ADCost(cost) << std::setprecision(DECIMAL_POINT) << " in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << ((ToulBar2::DEE) ? (" ( " + to_string(wcsp->getNbDEE()) + " removals by DEE)") : "") << " and " << cpuTime() - ToulBar2::startCpuTime << " seconds." << endl;
            else
                cout << solType[isLimited] << cost << " energy: " << -(wcsp->Cost2LogProb(cost) + ToulBar2::markov_log) << std::scientific << " prob: " << wcsp->Cost2Prob(cost) * Exp(ToulBar2::markov_log) << std::fixed << " in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << ((ToulBar2::DEE) ? (" ( " + to_string(wcsp->getNbDEE()) + " removals by DEE)") : "") << " and " << cpuTime() - ToulBar2::startCpuTime << " seconds." << endl;
        } else {
            if (ToulBar2::xmlflag) {
                ((WCSP*)wcsp)->solution_XML(true);
            } else if (ToulBar2::uai && !ToulBar2::isZ) {
                if (isLimited == 2)
                    cout << "(" << ToulBar2::deltaUbS << "," << std::scientific << ToulBar2::deltaUbRelativeGap << std::fixed << ")-";
                cout << solType[isLimited] << cost << " energy: " << -(wcsp->Cost2LogProb(cost) + ToulBar2::markov_log) << std::scientific << " prob: " << wcsp->Cost2Prob(cost) * Exp(ToulBar2::markov_log) << std::fixed << " in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << ((ToulBar2::DEE) ? (" ( " + to_string(wcsp->getNbDEE()) + " removals by DEE)") : "") << " and " << cpuTime() - ToulBar2::startCpuTime << " seconds." << endl;
            } else if (ToulBar2::maxsateval && !isLimited) {
                cout << "o " << cost << endl;
                cout << "s OPTIMUM FOUND" << endl;
                ((WCSP*)wcsp)->printSolutionMaxSAT(cout);
            }
        }
    } else {
        if (ToulBar2::verbose >= 0)
            cout << "No solution" << ((!isLimited) ? "" : " found") << " in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << ((ToulBar2::DEE) ? (" ( " + to_string(wcsp->getNbDEE()) + " removals by DEE)") : "") << " and " << cpuTime() - ToulBar2::startCpuTime << " seconds." << endl;
        if (ToulBar2::maxsateval && !isLimited) {
            cout << "o " << cost << endl;
            cout << "s UNSATISFIABLE" << endl;
        }
    }
}

void Solver::approximate(BigInteger& nbsol, TreeDecomposition* td)
{
    BigInteger cartesianProduct = 1;
    wcsp->cartProd(cartesianProduct);
    for (map<int, BigInteger>::iterator it = ubSol.begin(); it != ubSol.end(); ++it) {
        (it->second) *= cartesianProduct;
    }
    BigInteger nbSolInter = nbsol * cartesianProduct;
    BigInteger subCartesianProduct = 1.;
    for (int i = 0; i < td->getNbOfClusters(); i++) {
        BigInteger ssCartProd = 1.;
        if ((td->getCluster(i)->getParent() != NULL) && (td->getCluster(i)->getParent()->getParent() == NULL)) {
            /* on considere seulement les clusters fils de la racine */
            Cluster* c = td->getCluster(i);
            c->cartProduct(ssCartProd);
            subCartesianProduct *= ssCartProd;
            (ubSol.find(c->getPart())->second) /= ssCartProd;
        }
    }
    nbsol = (nbSolInter / subCartesianProduct);
    if (nbsol < 1)
        nbsol = 1;
    // the minimum upper bound of solutions number
    cout << "\nCartesian product \t\t   :    " << std::fixed << std::setprecision(0) << cartesianProduct << std::setprecision(DECIMAL_POINT) << endl;
    BigInteger minUBsol = cartesianProduct;
    for (map<int, BigInteger>::iterator it = ubSol.begin(); it != ubSol.end(); ++it) {
        if (it->second < minUBsol)
            minUBsol = it->second;
    }
    cout << "Upper bound of number of solutions : <= " << std::fixed << std::setprecision(0) << minUBsol << std::setprecision(DECIMAL_POINT) << endl;
}

// Maximize h' W h where W is expressed by all its
// non-zero half squared matrix costs (can be positive or negative, with posx <= posy)
// notice that costs for posx <> posy are multiplied by 2 by this method

// convention: h = 1 <=> x = 0 and h = -1 <=> x = 1

// warning! does not allow infinite costs (no forbidden assignments)

// returns true if at least one solution has been found (array sol being filled with the best solution)
bool Solver::solve_symmax2sat(int n, int m, int* posx, int* posy, double* cost, int* sol)
{
    if (n == 0 || m == 0)
        return true;
    ToulBar2::setvalue = NULL;

    // create Boolean variables
    for (int i = 0; i < n; i++) {
        wcsp->makeEnumeratedVariable(to_string(i), 0, 1);
    }

    vector<Cost> unaryCosts0(n, 0);
    vector<Cost> unaryCosts1(n, 0);

    // find total cost
    Double sumcost = 0.;
    for (int e = 0; e < m; e++) {
        sumcost += 2. * abs(cost[e]);
    }
    Double multiplier = ((Double)MAX_COST) / sumcost;
    multiplier /= MEDIUM_COST;

    // create weighted binary clauses
    for (int e = 0; e < m; e++) {
        if (posx[e] != posy[e]) {
            vector<Cost> costs(4, 0);
            if (cost[e] > 0) {
                costs[1] = (Cost)(multiplier * 2. * cost[e]);
                costs[2] = costs[1];
            } else {
                costs[0] = (Cost)(multiplier * -2. * cost[e]);
                costs[3] = costs[0];
            }
            wcsp->postBinaryConstraint(posx[e] - 1, posy[e] - 1, costs);
        } else {
            if (cost[e] > 0) {
                unaryCosts1[posx[e] - 1] += (Cost)(multiplier * cost[e]);
            } else {
                unaryCosts0[posx[e] - 1] += (Cost)(multiplier * -cost[e]);
            }
        }
    }

    // create weighted unary clauses
    for (int i = 0; i < n; i++) {
        if (unaryCosts0[i] > 0 || unaryCosts1[i] > 0) {
            vector<Cost> costs(2, 0);
            costs[0] = unaryCosts0[i];
            costs[1] = unaryCosts1[i];
            wcsp->postUnary(i, costs);
        }
    }
    wcsp->sortConstraints();

    if (ToulBar2::verbose >= 0)
        cout << "Read " << n << " variables, with " << 2 << " values at most, and " << m << " cost functions." << endl;
    // dump_wcsp("mydebug.wcsp", true);

    // solve using BTD exploiting a lexicographic elimination order with a path decomposition

    ToulBar2::btdMode = 3;
    ToulBar2::minProperVarSize = 4;
    ToulBar2::elimDegree_preprocessing = 12; // Prefer variable elimination than search (do not impose a limit on maximum separator size)

    bool res = solve();
    if (res) {
        vector<Value> solution = getSolution();
        assert(solution.size() == getWCSP()->numberOfVariables());
        for (unsigned int i = 0; i < getWCSP()->numberOfVariables(); i++) {
            if (solution[i] == 0) {
                sol[i] = 1;
            } else {
                sol[i] = -1;
            }
        }
    }
    return res;
}

/// \brief interface for Fortran call
/// \code
/// integer     ,dimension(sW),intent(out)        :: H
/// integer     ,dimension(:)    ,allocatable       :: posx,posy ! On exit dimension value is m
/// real(kind=dp),dimension(:)   ,allocatable   :: cost
/// logical :: ok
/// allocate (posx(sW*sW),posy(sW*sW),cost(sW*sW))
/// ret = solvesymmax2sat_(n,m,posx,posy,cost,H)
/// ok = ( ret /= 0 )
/// deallocate(posx,posy,cost)
/// \endcode
int solvesymmax2sat_(int* n, int* m, int* posx, int* posy, double* cost, int* sol)
{
    return solveSymMax2SAT(*n, *m, posx, posy, cost, sol);
}

int solveSymMax2SAT(int n, int m, int* posx, int* posy, double* cost, int* sol)
{
    // select verbosity during search
    ToulBar2::verbose = -1;

    initCosts();
    Solver solver(MAX_COST);

    ToulBar2::startCpuTime = cpuTime();
    return solver.solve_symmax2sat(n, m, posx, posy, cost, sol);
}

/* Hybrid Best-First/Depth-First Search */

void Solver::CPStore::addChoicePoint(ChoicePointOp op, int varIndex, Value value, bool reverse)
{
    if (ToulBar2::verbose >= 1)
        cout << "add choice point " << CPOperation[op] << ((reverse) ? "*" : "") << " (" << varIndex << ", " << value << ") at position " << index << endl;
    if ((size_t)index >= size()) {
        assert((size_t)index == size());
        push_back(ChoicePoint(op, varIndex, value, reverse));
    } else {
        operator[](index) = ChoicePoint(op, varIndex, value, reverse);
    }
    index = index + 1;
}

void Solver::addChoicePoint(ChoicePointOp op, int varIndex, Value value, bool reverse)
{
    TreeDecomposition* td = wcsp->getTreeDec();
    if (td) {
        if (ToulBar2::verbose >= 1)
            cout << "[C" << td->getCurrentCluster()->getId() << "] ";
        CPStore* cp_ = td->getCurrentCluster()->cp;
        CPStore::size_type before = cp_->capacity();
        cp_->addChoicePoint(op, varIndex, value, reverse);
        CPStore::size_type after = cp_->capacity();
        if (ToulBar2::verbose >= 0 && after > before && after > (1 << STORE_SIZE))
            cout << "c " << after * sizeof(ChoicePointOp) + td->getCurrentCluster()->open->capacity() * sizeof(OpenNode) << " Bytes allocated for hybrid best-first search open nodes at cluster " << td->getCurrentCluster()->getId() << "." << endl;
    } else {
        CPStore::size_type before = cp->capacity();
        cp->addChoicePoint(op, varIndex, value, reverse);
        CPStore::size_type after = cp->capacity();
        if (ToulBar2::verbose >= 0 && after > before && after > (1 << STORE_SIZE))
            cout << "c " << after * sizeof(ChoicePointOp) + open->capacity() * sizeof(OpenNode) << " Bytes allocated for hybrid best-first search open nodes." << endl;
    }
}

void Solver::addOpenNode(CPStore& cp, OpenList& open, Cost lb, Cost delta)
{
    ptrdiff_t idx = cp.index;
    if (ToulBar2::verbose >= 1) {
        if (wcsp->getTreeDec())
            cout << "[C" << wcsp->getTreeDec()->getCurrentCluster()->getId() << "] ";
        cout << "add open node " << lb << " + " << delta << " (" << cp.start << ", " << idx << ")" << endl;
    }
    assert(cp.start <= idx);
    open.push(OpenNode(MAX(MIN_COST, lb + delta), cp.start, idx));

    cp.stop = max(cp.stop, idx);
}

//// BUG: not compatible with boosting search by variable elimination (default dummy assignment may be incompatible with restored choice point)
//void Solver::restore(CPStore &cp, OpenNode nd)
//{
//    if (ToulBar2::verbose >= 1) {
//        if (wcsp->getTreeDec()) cout << "[C" << wcsp->getTreeDec()->getCurrentCluster()->getId() << "] ";
//        cout << "restore open node " << nd.getCost(((wcsp->getTreeDec())?wcsp->getTreeDec()->getCurrentCluster()->getCurrentDelta():MIN_COST)) << " (" << nd.first << ", " << nd.last << ")" << endl;
//    }
//    for (ptrdiff_t idx = nd.first; idx < nd.last; ++idx) {
//        assert(idx < cp.size());
//        if (ToulBar2::verbose >= 1) cout << "retrieve choice point " << CPOperation[cp[idx].op] << ((cp[idx].reverse)?"*":"") << " (" << wcsp->getName(cp[idx].varIndex) << ", " << cp[idx].value << ") at position " << idx  << endl;
//        if (ToulBar2::verbose >= 1) cout << *((WCSP *) wcsp)->getVar(cp[idx].varIndex) << endl;
//        assert(!wcsp->getTreeDec() || wcsp->getTreeDec()->getCurrentCluster()->isVar(cp[idx].varIndex));
//        switch (cp[idx].op) {
//        case CP_ASSIGN:
//            if (cp[idx].reverse && idx < nd.last-1) remove(cp[idx].varIndex, cp[idx].value);
//            else assign(cp[idx].varIndex, cp[idx].value);
//            break;
//        case CP_REMOVE:
//            if (cp[idx].reverse && idx < nd.last-1) assign(cp[idx].varIndex, cp[idx].value);
//            else remove(cp[idx].varIndex, cp[idx].value);
//            break;
//        case CP_INCREASE:
//            if (cp[idx].reverse && idx < nd.last-1) decrease(cp[idx].varIndex, cp[idx].value - 1);
//            else increase(cp[idx].varIndex, cp[idx].value);
//            break;
//        case CP_DECREASE:
//            if (cp[idx].reverse && idx < nd.last-1) increase(cp[idx].varIndex, cp[idx].value + 1);
//            else decrease(cp[idx].varIndex, cp[idx].value);
//            break;
//        default:
//            cerr << "unknown choice point for hybrid best first search!!!" << endl;
//            exit(EXIT_FAILURE);
//        }
//    }
//}

void Solver::restore(CPStore& cp, OpenNode nd)
{
    if (ToulBar2::verbose >= 1) {
        if (wcsp->getTreeDec()) {
            cout << "[C" << wcsp->getTreeDec()->getCurrentCluster()->getId() << "] restore open node " << nd.getCost(wcsp->getTreeDec()->getCurrentCluster()->getCurrentDelta()) << " (" << nd.first << ", " << nd.last << ")" << endl;
        } else {
            cout << "restore open node " << nd.getCost(MIN_COST) << " (" << nd.first << ", " << nd.last << ")" << endl;
        }
    }
    assert(nd.last >= nd.first);
    nbRecomputationNodes += nd.last - nd.first;

    ptrdiff_t maxsize = nd.last - nd.first;
    int assignLS[maxsize];
    Value valueLS[maxsize];
    unsigned int size = 0;
    for (ptrdiff_t idx = nd.first; idx < nd.last; ++idx) {
        assert((size_t)idx < cp.size());
        assert(!wcsp->getTreeDec() || wcsp->getTreeDec()->getCurrentCluster()->isVar(cp[idx].varIndex));
        if ((cp[idx].op == CP_ASSIGN && !(cp[idx].reverse && idx < nd.last - 1)) || (cp[idx].op == CP_REMOVE && cp[idx].reverse && idx < nd.last - 1)) {
            assignLS[size] = cp[idx].varIndex;
            valueLS[size] = cp[idx].value;
            size++;
        }
    }
    wcsp->enforceUb();
    wcsp->assignLS(assignLS, valueLS, size, false); // fast multiple assignments
    for (ptrdiff_t idx = nd.first; idx < nd.last; ++idx) {
        assert((size_t)idx < cp.size());
        if (ToulBar2::verbose >= 1)
            cout << "retrieve choice point " << CPOperation[cp[idx].op] << ((cp[idx].reverse) ? "*" : "") << " (" << wcsp->getName(cp[idx].varIndex) << ", " << cp[idx].value << ") at position " << idx << endl;
        if (ToulBar2::verbose >= 1)
            cout << *((WCSP*)wcsp)->getVar(cp[idx].varIndex) << endl;
        nbNodes++;
        switch (cp[idx].op) { //TODO: some operations (remove,increase,decrease) are useless because of all assigns previously done
        case CP_ASSIGN: {
            if (cp[idx].reverse && idx < nd.last - 1) {
                wcsp->remove(cp[idx].varIndex, cp[idx].value);
                addChoicePoint(CP_REMOVE, cp[idx].varIndex, cp[idx].value, false);
            } else
                addChoicePoint(CP_ASSIGN, cp[idx].varIndex, cp[idx].value, false);
            break;
        }
        case CP_REMOVE: {
            if (cp[idx].reverse && idx < nd.last - 1) {
                addChoicePoint(CP_ASSIGN, cp[idx].varIndex, cp[idx].value, false);
            } else {
                wcsp->remove(cp[idx].varIndex, cp[idx].value);
                addChoicePoint(CP_REMOVE, cp[idx].varIndex, cp[idx].value, false);
            }
            break;
        }
        case CP_INCREASE: {
            if (cp[idx].reverse && idx < nd.last - 1) {
                wcsp->decrease(cp[idx].varIndex, cp[idx].value - 1);
                addChoicePoint(CP_DECREASE, cp[idx].varIndex, cp[idx].value - 1, false);
            } else {
                wcsp->increase(cp[idx].varIndex, cp[idx].value);
                addChoicePoint(CP_INCREASE, cp[idx].varIndex, cp[idx].value, false);
            }
            break;
        }
        case CP_DECREASE: {
            if (cp[idx].reverse && idx < nd.last - 1) {
                wcsp->increase(cp[idx].varIndex, cp[idx].value + 1);
                addChoicePoint(CP_INCREASE, cp[idx].varIndex, cp[idx].value + 1, false);
            } else {
                wcsp->decrease(cp[idx].varIndex, cp[idx].value);
                addChoicePoint(CP_DECREASE, cp[idx].varIndex, cp[idx].value, false);
            }
            break;
        }
        default: {
            cerr << "unknown choice point for hybrid best first search!!!" << endl;
            exit(EXIT_FAILURE);
        }
        }
    }
    wcsp->propagate();
    //if (wcsp->getLb() != nd.getCost(((wcsp->getTreeDec())?wcsp->getTreeDec()->getCurrentCluster()->getCurrentDelta():MIN_COST))) cout << "***** node cost: " << nd.getCost(((wcsp->getTreeDec())?wcsp->getTreeDec()->getCurrentCluster()->getCurrentDelta():MIN_COST)) << " but lb: " << wcsp->getLb() << endl;
}

Solver::SolutionTrie::TrieNode::TrieNode(size_t w)
{
    sons.resize(w, NULL);
}

Solver::SolutionTrie::TrieNode::~TrieNode()
{
    for (size_t i = 0; i < sons.size(); i++)
        delete sons[i];
}

vector<size_t> Solver::SolutionTrie::TrieNode::widths;

bool Solver::SolutionTrie::TrieNode::present(Value v)
{
    return (sons[v] != NULL);
}

vector<vector<Solver::SolutionTrie::TrieNode*>> Solver::SolutionTrie::TrieNode::insertNode(Value v, unsigned int pos, vector<vector<TrieNode*>> nodesAtPos)
{
    sons[v] = new TrieNode(widths[pos + 1]);
    nodesAtPos[pos + 1].push_back(sons[v]);
    return nodesAtPos;
}

vector<vector<Solver::SolutionTrie::TrieNode*>> Solver::SolutionTrie::TrieNode::insertSolution(const vector<Value>& sol, unsigned int pos, vector<vector<TrieNode*>> nodesAtPos)
{
    if (pos < sol.size()) {
        if (!present(sol[pos])) {
            nodesAtPos = insertNode(sol[pos], pos, nodesAtPos);
        }
        assert((unsigned)sol[pos] < sons.size());
        return sons[sol[pos]]->insertSolution(sol, pos + 1, nodesAtPos);
    } else {
        return nodesAtPos;
    }
}

void Solver::SolutionTrie::init(const vector<Variable*>& vv)
{
    for (auto var : vv) {
        Solver::SolutionTrie::TrieNode::widths.push_back(((EnumeratedVariable*)var)->getDomainInitSize());
    }
    Solver::SolutionTrie::TrieNode::widths.push_back(0); // for leaf nodes
    if (!vv.empty()) {
        root.sons.resize(Solver::SolutionTrie::TrieNode::widths[0], NULL);
        nodesAtPos.resize(vv.size() + 1);
        nodesAtPos[0].push_back(&root);
    }
}
void Solver::SolutionTrie::TrieNode::printTrie(vector<Value>& sol)
{
    if (sons.size() == 0) {
        cout << sol << endl;
    } else {
        for (size_t i = 0; i < sons.size(); i++)
            if (sons[i] != NULL) {
                sol.push_back(i);
                sons[i]->printTrie(sol);
                sol.pop_back();
            }
    }
}

size_t Solver::SolutionTrie::TrieNode::nbSolutions = 0;

void Solver::SolutionTrie::insertSolution(const vector<Value>& sol)
{
    nodesAtPos = root.insertSolution(sol, 0, nodesAtPos);
}

void Solver::SolutionTrie::printTrie()
{
    vector<Value> sol;
    root.printTrie(sol);
}

Mdd Solver::computeMDD(Solver::SolutionTrie* solTrie, Cost cost)
{
    //The SolutionTrie computes solution Prefix Tree, in divVariables order.
    // To merge equivalent nodes, we need the suffix tree, so we reverse the variables order.
    // variables
    vector<Variable*> varReverse;
    for (unsigned v = 0; v < wcsp->getDivVariables().size(); v++) {
        varReverse.push_back(wcsp->getDivVariables()[wcsp->getDivVariables().size() - v - 1]);
    }
    int nLayers = varReverse.size();
    Mdd mdd(nLayers);
    vector<int> layerWidth;
    layerWidth.push_back(1);
    //solTrie
    vector<vector<Solver::SolutionTrie::TrieNode*>> nodesAtLayer;

    for (unsigned pos = 0; pos < solTrie->getNodesAtPos().size(); pos++) {
        nodesAtLayer.push_back(solTrie->getNodesAtPos()[solTrie->getNodesAtPos().size() - pos - 1]);
    }
    map<vector<int>, int> DistCountsA;
    vector<int> initCount(nodesAtLayer[0].size(), 0);
    DistCountsA[initCount] = 0;
    map<vector<int>, int> DistCountsB;
    map<vector<int>, int>& prevDistCounts = DistCountsA;
    map<vector<int>, int>& nextDistCounts = DistCountsB;

    //for relaxation
    vector<vector<Cost>> oldArcs; // for arcs redirection during relaxation
    vector<vector<Cost>> alphap(nLayers + 1); // if divRelax=3, alphap[layer][node] = smallest path weight from root to node, with unary costs.
    alphap[0].push_back(0); // alphap at root = 0

    for (int layer = 0; layer < nLayers; layer++) { // layer = arcs
        EnumeratedVariable* x = (EnumeratedVariable*)varReverse[layer];

        mdd[layer].resize(prevDistCounts.size());
        for (unsigned source = 0; source < prevDistCounts.size(); source++) {
            mdd[layer][source].resize(ToulBar2::divWidth);
            for (unsigned target = 0; target < ToulBar2::divWidth; target++) {
                mdd[layer][source][target].resize(x->getDomainInitSize(), wcsp->getUb());
            }
        }

        unsigned source;
        unsigned target;
        Cost toPay;
        for (auto const& nodeState : prevDistCounts) {
            source = nodeState.second;
            for (unsigned val = 0; val < x->getDomainInitSize(); val++) {
                vector<int> nextCount(nodesAtLayer[layer + 1].size(), -1);
                for (unsigned node_index = 0; node_index < nodesAtLayer[layer + 1].size(); node_index++) {
                    auto node = nodesAtLayer[layer + 1][node_index];
                    for (unsigned sol = 0; sol < node->sons.size(); sol++) { // node <---sol----nodep
                        auto nodep = node->sons[sol];
                        if (nodep != NULL) {
                            auto nodep_it = find(nodesAtLayer[layer].begin(), nodesAtLayer[layer].end(), nodep);
                            unsigned nodep_index = distance(nodesAtLayer[layer].begin(), nodep_it);
                            assert(nodep_index < nodesAtLayer[layer].size());

                            if (nextCount[node_index] == -1) {
                                nextCount[node_index] = nodeState.first[nodep_index] + (val != sol);
                            } else {
                                nextCount[node_index] = min(nextCount[node_index], nodeState.first[nodep_index] + (val != sol));
                            }
                        }
                    }
                    nextCount[node_index] = min((unsigned int)nextCount[node_index], ToulBar2::divBound);
                }
                //merge nodes that won't lead to a satisfying solution:
                bool sat = true;
                for (int count : nextCount) {
                    if (count < (int)ToulBar2::divBound + layer - nLayers) {
                        sat = false;
                        break;
                    }
                }
                if (!sat)
                    nextCount.resize(nextCount.size(), 0);
                if (layer != nLayers - 1) {
                    toPay = MIN_COST;
                    map<vector<int>, int>::iterator it;
                    std::tie(it, std::ignore) = nextDistCounts.insert(pair<vector<int>, int>(nextCount, nextDistCounts.size()));
                    target = (*it).second;
                } else {
                    toPay = MIN_COST;
                    target = 0;
                    assert(layer + 1 == nLayers);
                    for (int count : nextCount) {
                        if (count < (int)ToulBar2::divBound) {
                            toPay = cost;
                            break;
                        }
                    }
                }
                if (target < mdd[layer][source].size()) {
                    mdd[layer][source][target][val] = toPay;
                } else {
                    vector<Cost> newTarget(x->getDomainInitSize(), wcsp->getUb());
                    // when a new target appears in the mdd, we need to add all arcs from ALL sources!!
                    for (unsigned s = 0; s < mdd[layer].size(); s++) {
                        mdd[layer][s].push_back(newTarget);
                    }
                    mdd[layer][source][target][val] = toPay;
                }
            }
        }
        unsigned nTargets = nextDistCounts.size();
        if (nTargets > ToulBar2::divWidth) {
            // select nodes for merging
            int n_merge = nTargets - ToulBar2::divWidth + 1;
            vector<int> to_merge(nextDistCounts.size(), -1);
            iota(to_merge.begin(), to_merge.end(), 0);
            if (ToulBar2::divRelax == 0) {
                for (int i = 0; i < n_merge; ++i) {
                    int j = myrand() % (nextDistCounts.size() - i);
                    std::swap(to_merge[i], to_merge[i + j]);
                }
            } else if (ToulBar2::divRelax == 1) {
                vector<int> stateDiv(nextDistCounts.size());
                for (const auto& node : nextDistCounts) {
                    stateDiv[node.second] = accumulate(node.first.begin(), node.first.end(), 0);
                }
                auto comparator = [stateDiv](int a, int b) { return stateDiv[a] > stateDiv[b]; };
                std::sort(to_merge.begin(), to_merge.end(), comparator);
            } else if (ToulBar2::divRelax == 2) {
                vector<int> stateDiv(nextDistCounts.size());
                for (const auto& node : nextDistCounts) {
                    stateDiv[node.second] = accumulate(node.first.begin(), node.first.end(), 0);
                }
                auto comparator = [stateDiv](int a, int b) { return stateDiv[a] < stateDiv[b]; };
                std::sort(to_merge.begin(), to_merge.end(), comparator);
            } else if (ToulBar2::divRelax == 3) {
                vector<Cost> alphaptmp(nTargets, wcsp->getUb());
                for (unsigned source = 0; source < mdd[layer].size(); source++) {
                    for (unsigned target = 0; target < mdd[layer][source].size(); target++) {
                        for (unsigned val = 0; val < x->getDomainInitSize(); val++) {
                            alphaptmp[target] = min(alphaptmp[target], alphap[layer][source] + mdd[layer][source][target][val] + x->getCost(x->toValue(val)));
                        }
                    }
                }
                auto comparator = [alphaptmp](int a, int b) { return alphaptmp[a] > alphaptmp[b]; };
                std::sort(to_merge.begin(), to_merge.end(), comparator);
            } else {
                cerr << "Error: no such relaxing method: " << ToulBar2::divRelax;
            }

            to_merge.resize(n_merge);

            //Merging nodes:TODO
            //Computing new state for merged nodes
            vector<int> newCount(nodesAtLayer[layer + 1].size(), -1);

            vector<int> newTarget(nextDistCounts.size(), -1); //vector with new state nodes ids
            for (int state_index : to_merge) {
                auto state_it = std::find_if(nextDistCounts.begin(), nextDistCounts.end(), [state_index](const pair<vector<int>, int>& mo) { return mo.second == state_index; });
                assert(state_it != nextDistCounts.end());
                for (unsigned nodeid = 0; nodeid < nodesAtLayer[layer + 1].size(); nodeid++) {
                    if (newCount[nodeid] == -1) {
                        newCount[nodeid] = state_it->first[nodeid];
                    } else {
                        newCount[nodeid] = max(newCount[nodeid], state_it->first[nodeid]);
                        // we are allowing more solutions - we don't want to remove any solution (!! relaxation)
                        // exact mdds for each single solution are required
                    }
                }
                nextDistCounts.erase(state_it->first);
            }
            //The nodes need to be renumbered - we want nodeids = 0, 1 , ... , divWidth
            map<vector<int>, int>::iterator it;
            std::tie(it, std::ignore) = nextDistCounts.insert(pair<vector<int>, int>(newCount, to_merge[0]));
            int newNode = (*it).second;
            unsigned nodeid = 0;
            for (auto state : nextDistCounts) {
                if (state.second == newNode) {
                    newTarget[newNode] = nodeid;
                    for (auto node : to_merge) {
                        newTarget[node] = nodeid;
                    }
                } else {
                    newTarget[state.second] = nodeid;
                }
                nextDistCounts[state.first] = nodeid;
                nodeid++;
            }
            //redirecting arcs in mdd[layer] from each source to new targets
            for (unsigned source = 0; source < mdd[layer].size(); source++) {
                oldArcs.clear();
                oldArcs = mdd[layer][source];
                mdd[layer][source].resize(nodeid);
                for (unsigned target = 0; target < nodeid; target++) {
                    mdd[layer][source][target].resize(x->getDomainInitSize());
                    for (unsigned val = 0; val < x->getDomainInitSize(); val++) {
                        mdd[layer][source][target][val] = wcsp->getUb(); // erase all arcs from source
                    }
                }
                for (unsigned oldTarget = 0; oldTarget < oldArcs.size(); oldTarget++) {
                    for (unsigned val = 0; val < x->getDomainInitSize(); val++) {
                        mdd[layer][source][newTarget[oldTarget]][val] = min(mdd[layer][source][newTarget[oldTarget]][val], oldArcs[oldTarget][val]);
                    }
                }
            }
        }
        if (ToulBar2::divWidth > 0 && ToulBar2::divRelax == 3 && layer != nLayers - 1) {
            //Computing alphap[layer+1]
            alphap[layer + 1].resize(ToulBar2::divWidth, wcsp->getUb());
            for (unsigned source = 0; source < mdd[layer].size(); source++) {
                for (unsigned target = 0; target < mdd[layer][source].size(); target++) {
                    for (unsigned val = 0; val < x->getDomainInitSize(); val++) {
                        alphap[layer + 1][target] = min(alphap[layer + 1][target], alphap[layer][source] + mdd[layer][source][target][val] + x->getCost(x->toValue(val)));
                    }
                }
            }
        }

        layerWidth.push_back((layer != nLayers - 1) ? nextDistCounts.size() : 1);
        prevDistCounts.clear();
        swap(prevDistCounts, nextDistCounts);
    }
    return mdd;
}

std::ostream& Solver::printLayers(std::ostream& os, Mdd mdd)
{

    os << "digraph \"wregular\" {" << endl;
    os << "\tgraph [hierarchic=1];" << endl;
    // Draw vertices
    int nodeShift = 0;
    for (unsigned layer = 0; layer < mdd.size(); layer++) {
        for (unsigned node = 0; node < mdd[layer].size(); node++) {
            os << "\t" << nodeShift + node << " [name=\"" << layer << "," << node << "\"];" << endl;
        }
        nodeShift += mdd[layer].size();
    }
    // and Arcs
    nodeShift = 0;
    for (unsigned layer = 0; layer < mdd.size(); layer++) {
        for (unsigned source = 0; source < mdd[layer].size(); source++) {
            for (unsigned target = 0; target < mdd[layer][source].size(); target++) {
                for (unsigned val = 0; val < mdd[layer][source][target].size(); val++) {
                    if (mdd[layer][source][target][val] < wcsp->getUb()) {
                        os << "\t" << nodeShift + source << " -> " << nodeShift + mdd[layer].size() + target << " [label=\"";
                        os << val << "," << mdd[layer][source][target][val] << "\"];" << endl;
                    }
                }
            }
        }
        nodeShift += mdd[layer].size();
    }
    os << "}";
    return os;
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

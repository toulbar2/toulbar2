/*
 * **************** Generic solver *******************
 *
 * Parallel HBFS version thanks to Abdelkader Beldjilali
 *
 */

#include "tb2solver.hpp"
#include "core/tb2vac.hpp"
#include "core/tb2domain.hpp"
#include "core/tb2globalwcsp.hpp"
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
#include <thread>

extern void setvalue(int wcspId, int varIndex, Value value, void* solver);
extern void tb2setvalue(int wcspId, int varIndex, Value value, void* solver);
extern void newsolution(int wcspId, void* solver);

const string Solver::CPOperation[CP_MAX] = { "ASSIGN", "REMOVE", "INCREASE", "DECREASE", "RANGEREMOVAL" };

Solver* Solver::CurrentSolver;

/*
 * Solver constructors
 *
 */

WeightedCSPSolver* WeightedCSPSolver::makeWeightedCSPSolver(Cost ub, WeightedCSP* wcsp)
{
    WeightedCSPSolver* solver = NULL;

    if (wcsp != NULL && ToulBar2::searchMethod != DFBB) {
        cerr << "Error: provided WeightedCSP object not taken into account by this solver method " << ToulBar2::searchMethod << endl;
        throw BadConfiguration();
    }

    switch (ToulBar2::searchMethod) {
    case VNS:
    case DGVNS:
#ifdef BOOST
        solver = new VNSSolver(ub);
#else
        cerr << "Error: compiling with Boost graph library is needed to allow VNS-like search methods." << endl;
        throw BadConfiguration();
#endif
        break;
#ifdef OPENMPI
    case CPDGVNS:
#ifdef BOOST
        solver = new CooperativeParallelDGVNS(ub);
#else
        cerr << "Error: compiling with Boost graph library is needed to allow VNS-like search methods." << endl;
        throw BadConfiguration();
#endif
        break;
    case RPDGVNS:
#ifdef BOOST
        solver = new ReplicatedParallelDGVNS(ub);
#else
        cerr << "Error: compiling with Boost graph library is needed to allow VNS-like search methods." << endl;
        throw BadConfiguration();
#endif
        break;
#endif
    case TREEDEC:
#ifdef BOOST
        solver = new TreeDecRefinement(ub);
#else
        cerr << "Error: compiling with Boost graph library is needed to allow VNS-like search methods." << endl;
        throw BadConfiguration();
#endif
        break;
    default:
        solver = new Solver(ub, wcsp);
        break;
    };
    return solver;
}

Solver::Solver(Cost initUpperBound, WeightedCSP* wcsp)
    : self_wcsp(true)
    , nbNodes(0)
    , nbBacktracks(0)
    , nbBacktracksLimit(LONGLONG_MAX)
    , wcsp(NULL)
    , unassignedVars(NULL)
    , lastConflictVar(-1)
    , nbSol(0.)
    , nbSGoods(0)
    , nbSGoodsUse(0)
    , timeDeconnect(0.)
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
    , nbChoices(0)
    , nbForcedChoices(0)
    , nbForcedChoiceChange(0)
    , nbChoiceChange(0)
    , nbReadOnly(0)
    , solveDepth(0)
#ifdef OPENMPI
    , hbfsWaitingTime(0.)
    , initWorkerNbNodes(0)
    , initWorkerNbBacktracks(0)
    , initWorkerNbDEE(0)
    , initWorkerNbRecomputationNodes(0)
#endif
{
    searchSize = new StoreInt(0);

    if (wcsp == NULL) {
        this->wcsp = WeightedCSP::makeWeightedCSP(initUpperBound, (void*)this);
    } else { /* the provided wcsp will be solved */
        self_wcsp = false;
        this->wcsp = wcsp;
        dynamic_cast<WCSP*>(wcsp)->setSolver((void*)this);
    }
    CurrentSolver = this;
}

Solver::~Solver()
{
    if (cp) {
        delete cp;
    }
    if (open) {
        delete open;
    }
    if (wcsp->getTreeDec()) {
        Cluster *cluster = wcsp->getTreeDec()->getRoot();
        if (cluster && cluster->open) {
            delete cluster->open;
        }
    }
    delete unassignedVars;
    for (unsigned int i = 0; i < allVars.size(); i++) {
        delete allVars[i];
    }

    // if the wcsp has been created internally, then it is deleted else must be done by the user
    if (self_wcsp) {
        delete wcsp;
    }

    for (auto elt : WCSP::CollectionOfWCSP) {
        delete elt.second;
    }
    WCSP::CollectionOfWCSP.clear();

    delete ((StoreInt*)searchSize);
}

void Solver::initVarHeuristic()
{
    for (unsigned int j = 0; j < wcsp->numberOfVariables(); j++) {
        allVars.push_back(new DLink<Value>());
    }
    unassignedVars = new BTList<Value>(&Store::storeDomain);
    for (unsigned int j = 0; j < wcsp->numberOfVariables(); j++) {
        unsigned int i = wcsp->getDACOrder(j);
        allVars[i]->content = j;
    }
    heuristics.clear();
    for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
        heuristics.push_back(wcsp->getDegree(i));
        unassignedVars->push_back(allVars[i], false);
        if (wcsp->assigned(allVars[i]->content) || (ToulBar2::nbDecisionVars > 0 && allVars[i]->content >= ToulBar2::nbDecisionVars))
            unassignedVars->erase(allVars[i], false);
    }
    wcsp->resetTightnessAndWeightedDegree();
    // Now function setvalue can be called safely!
    if (ToulBar2::setvalue == NULL) {
        ToulBar2::setvalue = setvalue;
    }
}

// keep consistent order of allVars with getDACOrder, needed by setvalue function
void Solver::updateVarHeuristic()
{
    sort(allVars.begin(), allVars.end(),
        [this](const DLink<Value>* v1, const DLink<Value>* v2) -> bool {
            return (wcsp->getDACOrder(v1->content) < wcsp->getDACOrder(v2->content));
        });
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
                    throw WrongFileFormat();
                }
                operation++;
            } else {
                operation = 0;
            }
            if (not isdigit(token[operation])) {
                unsigned int idx = wcsp->toIndex(i, token.substr(operation, string::npos));
                if (idx >= wcsp->getDomainInitSize(i)) {
                    cerr << "Solution file incorrect! " << i_copy << " " << token << " " << i << " " << idx << endl;
                    throw WrongFileFormat();
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
        if (ToulBar2::verbose >= 0) {
            if (ToulBar2::bayesian) {
                cout << " Input solution cost: " << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->getDDualBound() << std::setprecision(DECIMAL_POINT) << " (nb. of unassigned variables: " << wcsp->numberOfUnassignedVariables() << ")"
                     << " energy: " << -(wcsp->Cost2LogProb(wcsp->getLb() + wcsp->getNegativeLb()) + ToulBar2::markov_log) << " prob: " << std::scientific << wcsp->Cost2Prob(wcsp->getLb() + wcsp->getNegativeLb()) * Exp(ToulBar2::markov_log) << std::fixed << endl;
                if (ToulBar2::uaieval)
                    cout << "Energy_log10: " << (wcsp->Cost2LogProb(wcsp->getLb() + wcsp->getNegativeLb()) + ToulBar2::markov_log) / log(10) << endl;
            } else {
                cout << " Input solution cost: " << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->getDDualBound() << std::setprecision(DECIMAL_POINT) << " (nb. of unassigned variables: " << wcsp->numberOfUnassignedVariables() << ")" << endl;
            }
        }
        if (ToulBar2::verifyOpt) {
            ToulBar2::verifiedOptimum = wcsp->getLb();
        } else {
            if (wcsp->numberOfUnassignedVariables() == 0) {
                wcsp->updateUb(wcsp->getLb() + ((updateValueHeuristic) ? UNIT_COST : MIN_COST));
            }
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

    // certif2 = index(certif2,',');
    char *certif2, *certif_copy;
    char sep[] = ",";
    certif_copy = strdup(certificate);
    certif2 = strstr(certif_copy, sep);

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
        char* ope = strpbrk(certif2 + 1, "=#<>"); // avoid first character of a variable name (can be an operation char)
        if (ope) {
            items++;
            operation = *ope;
            items++;
            strncpy(svar, certif2, ope - certif2);
            svar[ope - certif2] = '\0';
            char* nextsep = strpbrk(ope + 2, sep);
            if (nextsep) {
                items++;
                strncpy(svalue, ope + 1, nextsep - ope - 1);
                svalue[nextsep - ope - 1] = '\0';
            } else {
                if (strlen(ope + 1) > 0) {
                    items++;
                    strcpy(svalue, ope + 1);
                }
            }
        }
        if (items != 3) {
            cerr << "Certificate " << certif2 << " incorrect! " << items << endl;
            throw WrongFileFormat();
        }
        certif2 = strstr(certif2, sep);
        if (certif2)
            certif2++;

        if (not isdigit(svar[0])) {
            var = wcsp->getVarIndex(to_string(svar));
            if ((unsigned int)var >= wcsp->numberOfVariables()) {
                cerr << "Certificate " << svar << " incorrect!" << endl;
                throw WrongFileFormat();
            }
        } else {
            var = atoi(svar);
        }

        if (not isdigit(svalue[0])) {
            unsigned int idx = wcsp->toIndex(var, to_string(svalue));
            if (idx >= wcsp->getDomainInitSize(var)) {
                cerr << "Certificate " << svalue << " incorrect!" << endl;
                throw WrongFileFormat();
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
            throw InternalError();
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

    free(certif_copy);
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
        throw WrongFileFormat();
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

int Solver::numberOfUnassignedVariables() const
{
    return ((unassignedVars) ? unassignedVars->getSize() : -1);
}

/*
 * Link between solver and wcsp: maintain a backtrackable list of unassigned variable indexes
 *
 */

void newsolution(int wcspId, void* _solver_)
{
    assert(_solver_);
    Solver* solver = (Solver*)_solver_;
    if (ToulBar2::searchMethod == DFBB && ToulBar2::vnsOptimum > 0 && solver->getWCSP()->getLb() > 0 && solver->getWCSP()->getLb() <= ToulBar2::vnsOptimum) {
        ToulBar2::limited = true;
        throw BestSolFound();
    }
}

void setvalue(int wcspId, int varIndex, Value value, void* _solver_)
{
    //    assert(wcspId == 0); // WARNING! assert not compatible with sequential execution of solve() method
    assert(_solver_);
    Solver* solver = (Solver*)_solver_;
    unsigned int i = solver->getWCSP()->getDACOrder(varIndex);
    if (i < solver->allVars.size() && !solver->allVars[i]->removed) {
        solver->unassignedVars->erase(solver->allVars[i], true);
    }
}

void tb2setvalue(int wcspId, int varIndex, Value value, void* solver)
{
    assert(WeightedCSPConstraint::MasterWeightedCSP);
    assert(wcspId == WeightedCSPConstraint::MasterWeightedCSP->getIndex() || WeightedCSPConstraint::WeightedCSPConstraints.find(wcspId) != WeightedCSPConstraint::WeightedCSPConstraints.end());
    WCSP* wcsp = NULL;
    Variable* masterVar = NULL;
    bool activeState = true;
    if (wcspId != WeightedCSPConstraint::MasterWeightedCSP->getIndex()) { // we came from a slave, wake up the master
        WeightedCSPConstraint* gc = WeightedCSPConstraint::WeightedCSPConstraints[wcspId];
        if (gc->problem && gc->problem->getIndex() == wcspId) {
            wcsp = gc->problem;
        } else {
            assert(gc->negproblem && gc->negproblem->getIndex() == wcspId);
            wcsp = gc->negproblem;
        }
        activeState = wcsp->isactivatePropagate();
        wcsp->deactivatePropagate();
        masterVar = WeightedCSPConstraint::WeightedCSPConstraints[wcspId]->getVar(varIndex);
        WeightedCSPConstraint::unprotect();
        masterVar->assign(value);
        WeightedCSPConstraint::protect(false);
    } else {
        // we came from the master
        wcsp = WeightedCSPConstraint::MasterWeightedCSP;
        activeState = wcsp->isactivatePropagate();
        wcsp->deactivatePropagate();
        masterVar = WeightedCSPConstraint::MasterWeightedCSP->getVar(varIndex);
        setvalue(WeightedCSPConstraint::MasterWeightedCSP->getIndex(), masterVar->wcspIndex, value, WeightedCSPConstraint::MasterWeightedCSP->getSolver());
        WeightedCSPConstraint::protect(true);
    }
    if (ToulBar2::verbose >= 2)
        cout << "EVENT: x" << varIndex << "_" << wcspId << " = " << value << endl;
    for (auto gc : WeightedCSPConstraint::WeightedCSPConstraints)
        if (gc.second->connected()) {
            int varCtrIndex = gc.second->getIndex(masterVar);
            if (varCtrIndex != -1) { // only for slave problems which are concerned by this variable
                if (gc.second->problem && wcspId != gc.second->problem->getIndex()) { // do not reenter inside the same problem as the one we came
                    assert(WeightedCSPConstraint::_protected_);
                    try {
                        gc.second->problem->getVar(varCtrIndex)->assign(value);
                    } catch (const Contradiction&) {
                        gc.second->problem->whenContradiction();
                        WeightedCSPConstraint::unprotect();
                        throw Contradiction();
                    }
                }
                if (gc.second->connected() && gc.second->negproblem && wcspId != gc.second->negproblem->getIndex()) { // do not reenter inside the same problem as the one we came
                    assert(WeightedCSPConstraint::_protected_);
                    try {
                        gc.second->negproblem->getVar(varCtrIndex)->assign(value);
                    } catch (const Contradiction&) {
                        gc.second->negproblem->whenContradiction();
                        WeightedCSPConstraint::unprotect();
                        throw Contradiction();
                    }
                }
            }
        }
    assert(!wcsp->isactivatePropagate());
    if (activeState)
        wcsp->reactivatePropagate();
    if (wcspId == WeightedCSPConstraint::MasterWeightedCSP->getIndex()) {
        WeightedCSPConstraint::unprotect();
    }
}

void tb2removevalue(int wcspId, int varIndex, Value value, void* solver)
{
    assert(WeightedCSPConstraint::MasterWeightedCSP);
    assert(wcspId == WeightedCSPConstraint::MasterWeightedCSP->getIndex() || WeightedCSPConstraint::WeightedCSPConstraints.find(wcspId) != WeightedCSPConstraint::WeightedCSPConstraints.end());
    WCSP* wcsp = NULL;
    Variable* masterVar = NULL;
    bool activeState = true;
    if (wcspId != WeightedCSPConstraint::MasterWeightedCSP->getIndex()) { // we came from a slave, wake up the master
        WeightedCSPConstraint* gc = WeightedCSPConstraint::WeightedCSPConstraints[wcspId];
        if (gc->problem && gc->problem->getIndex() == wcspId) {
            wcsp = gc->problem;
        } else {
            assert(gc->negproblem && gc->negproblem->getIndex() == wcspId);
            wcsp = gc->negproblem;
        }
        activeState = wcsp->isactivatePropagate();
        wcsp->deactivatePropagate();
        WeightedCSPConstraint::unprotect();
        masterVar = WeightedCSPConstraint::WeightedCSPConstraints[wcspId]->getVar(varIndex);
        masterVar->remove(value);
        WeightedCSPConstraint::protect(false);
    } else {
        // we came from the master
        wcsp = WeightedCSPConstraint::MasterWeightedCSP;
        activeState = wcsp->isactivatePropagate();
        wcsp->deactivatePropagate();
        masterVar = WeightedCSPConstraint::MasterWeightedCSP->getVar(varIndex);
        WeightedCSPConstraint::protect(true);
    }
    if (ToulBar2::verbose >= 2)
        cout << "EVENT: x" << varIndex << "_" << wcspId << " != " << value << endl;
    for (auto gc : WeightedCSPConstraint::WeightedCSPConstraints)
        if (gc.second->connected()) {
            int varCtrIndex = gc.second->getIndex(masterVar);
            if (varCtrIndex != -1) {
                if (gc.second->problem && wcspId != gc.second->problem->getIndex()) { // do not reenter inside the same problem as the one we came
                    assert(WeightedCSPConstraint::_protected_);
                    try {
                        gc.second->problem->getVar(varCtrIndex)->remove(value);
                    } catch (const Contradiction&) {
                        gc.second->problem->whenContradiction();
                        WeightedCSPConstraint::unprotect();
                        throw Contradiction();
                    }
                }
                if (gc.second->connected() && gc.second->negproblem && wcspId != gc.second->negproblem->getIndex()) { // do not reenter inside the same problem as the one we came
                    assert(WeightedCSPConstraint::_protected_);
                    try {
                        gc.second->negproblem->getVar(varCtrIndex)->remove(value);
                    } catch (const Contradiction&) {
                        gc.second->negproblem->whenContradiction();
                        WeightedCSPConstraint::unprotect();
                        throw Contradiction();
                    }
                }
            }
        }
    assert(!wcsp->isactivatePropagate());
    if (activeState)
        wcsp->reactivatePropagate();
    if (wcspId == WeightedCSPConstraint::MasterWeightedCSP->getIndex()) {
        WeightedCSPConstraint::unprotect();
    }
}

void tb2setmin(int wcspId, int varIndex, Value value, void* solver)
{
    assert(WeightedCSPConstraint::MasterWeightedCSP);
    assert(wcspId == WeightedCSPConstraint::MasterWeightedCSP->getIndex() || WeightedCSPConstraint::WeightedCSPConstraints.find(wcspId) != WeightedCSPConstraint::WeightedCSPConstraints.end());
    WCSP* wcsp = NULL;
    Variable* masterVar = NULL;
    bool activeState = true;
    if (wcspId != WeightedCSPConstraint::MasterWeightedCSP->getIndex()) { // we came from a slave, wake up the master
        WeightedCSPConstraint* gc = WeightedCSPConstraint::WeightedCSPConstraints[wcspId];
        if (gc->problem && gc->problem->getIndex() == wcspId) {
            wcsp = gc->problem;
        } else {
            assert(gc->negproblem && gc->negproblem->getIndex() == wcspId);
            wcsp = gc->negproblem;
        }
        activeState = wcsp->isactivatePropagate();
        wcsp->deactivatePropagate();
        WeightedCSPConstraint::unprotect();
        masterVar = WeightedCSPConstraint::WeightedCSPConstraints[wcspId]->getVar(varIndex);
        masterVar->increase(value);
        WeightedCSPConstraint::protect(false);
    } else {
        // we came from the master
        wcsp = WeightedCSPConstraint::MasterWeightedCSP;
        activeState = wcsp->isactivatePropagate();
        wcsp->deactivatePropagate();
        masterVar = WeightedCSPConstraint::MasterWeightedCSP->getVar(varIndex);
        WeightedCSPConstraint::protect(true);
    }
    if (ToulBar2::verbose >= 2)
        cout << "EVENT: x" << varIndex << "_" << wcspId << " >= " << value << endl;
    for (auto gc : WeightedCSPConstraint::WeightedCSPConstraints)
        if (gc.second->connected()) {
            int varCtrIndex = gc.second->getIndex(masterVar);
            if (varCtrIndex != -1) {
                if (gc.second->problem && wcspId != gc.second->problem->getIndex()) { // do not reenter inside the same problem as the one we came
                    assert(WeightedCSPConstraint::_protected_);
                    try {
                        gc.second->problem->getVar(varCtrIndex)->increase(value);
                    } catch (const Contradiction&) {
                        gc.second->problem->whenContradiction();
                        WeightedCSPConstraint::unprotect();
                        throw Contradiction();
                    }
                }
                if (gc.second->connected() && gc.second->negproblem && wcspId != gc.second->negproblem->getIndex()) { // do not reenter inside the same problem as the one we came
                    assert(WeightedCSPConstraint::_protected_);
                    try {
                        gc.second->negproblem->getVar(varCtrIndex)->increase(value);
                    } catch (const Contradiction&) {
                        gc.second->negproblem->whenContradiction();
                        WeightedCSPConstraint::unprotect();
                        throw Contradiction();
                    }
                }
            }
        }
    assert(!wcsp->isactivatePropagate());
    if (activeState)
        wcsp->reactivatePropagate();
    if (wcspId == WeightedCSPConstraint::MasterWeightedCSP->getIndex()) {
        WeightedCSPConstraint::unprotect();
    }
}

void tb2setmax(int wcspId, int varIndex, Value value, void* solver)
{
    assert(WeightedCSPConstraint::MasterWeightedCSP);
    WCSP* wcsp = NULL;
    assert(wcspId == WeightedCSPConstraint::MasterWeightedCSP->getIndex() || WeightedCSPConstraint::WeightedCSPConstraints.find(wcspId) != WeightedCSPConstraint::WeightedCSPConstraints.end());
    Variable* masterVar = NULL;
    bool activeState = true;
    if (wcspId != WeightedCSPConstraint::MasterWeightedCSP->getIndex()) { // we came from a slave, wake up the master
        WeightedCSPConstraint* gc = WeightedCSPConstraint::WeightedCSPConstraints[wcspId];
        if (gc->problem && gc->problem->getIndex() == wcspId) {
            wcsp = gc->problem;
        } else {
            assert(gc->negproblem && gc->negproblem->getIndex() == wcspId);
            wcsp = gc->negproblem;
        }
        activeState = wcsp->isactivatePropagate();
        wcsp->deactivatePropagate();
        WeightedCSPConstraint::unprotect();
        masterVar = WeightedCSPConstraint::WeightedCSPConstraints[wcspId]->getVar(varIndex);
        masterVar->decrease(value);
        WeightedCSPConstraint::protect(false);
    } else {
        // we came from the master
        wcsp = WeightedCSPConstraint::MasterWeightedCSP;
        activeState = wcsp->isactivatePropagate();
        wcsp->deactivatePropagate();
        masterVar = WeightedCSPConstraint::MasterWeightedCSP->getVar(varIndex);
        WeightedCSPConstraint::protect(true);
    }
    if (ToulBar2::verbose >= 2)
        cout << "EVENT: x" << varIndex << "_" << wcspId << " <= " << value << endl;
    for (auto gc : WeightedCSPConstraint::WeightedCSPConstraints)
        if (gc.second->connected()) {
            int varCtrIndex = gc.second->getIndex(masterVar);
            if (varCtrIndex != -1) {
                if (gc.second->problem && wcspId != gc.second->problem->getIndex()) { // do not reenter inside the same problem as the one we came
                    assert(WeightedCSPConstraint::_protected_);
                    try {
                        gc.second->problem->getVar(varCtrIndex)->decrease(value);
                    } catch (const Contradiction&) {
                        gc.second->problem->whenContradiction();
                        WeightedCSPConstraint::unprotect();
                        throw Contradiction();
                    }
                }
                if (gc.second->connected() && gc.second->negproblem && wcspId != gc.second->negproblem->getIndex()) { // do not reenter inside the same problem as the one we came
                    assert(WeightedCSPConstraint::_protected_);
                    try {
                        gc.second->negproblem->getVar(varCtrIndex)->decrease(value);
                    } catch (const Contradiction&) {
                        gc.second->negproblem->whenContradiction();
                        WeightedCSPConstraint::unprotect();
                        throw Contradiction();
                    }
                }
            }
        }
    assert(!wcsp->isactivatePropagate());
    if (activeState)
        wcsp->reactivatePropagate();
    if (wcspId == WeightedCSPConstraint::MasterWeightedCSP->getIndex()) {
        WeightedCSPConstraint::unprotect();
    }
}

/*
 * Variable ordering heuristics
 *
 */

/// \defgroup heuristics Variable and value search ordering heuristics
///
/// See : <em> Boosting Systematic Search by Weighting Constraints </em>. Frederic Boussemart, Fred Hemery, Christophe Lecoutre, Lakhdar Sais. Proc. of ECAI 2004, pages 146-150. Valencia, Spain, 2004.
///
/// See : <em> Last Conflict Based Reasoning </em>. Christophe Lecoutre, Lakhdar Sais, Sebastien Tabary, Vincent Vidal. Proc. of ECAI 2006, pages 133-137. Trentino, Italy, 2006.
///
/// See : <em> Solution-based phase saving for CP: A value-selection heuristic to simulate local search behavior in complete solvers </em>. Emir Demirovic, Geoffrey Chu, and Peter Stuckey. Proc. of CP-18, pages 99–108. Lille, France, 2018.

int Solver::getNextUnassignedVar()
{
    //    if (lastConflictVar != -1 && wcsp->unassigned(lastConflictVar)) return lastConflictVar;
    return (unassignedVars->empty()) ? -1 : (*unassignedVars->begin());
}

int Solver::greedy(intFunctionCall_t varHeurFunc)
{
    int varIndex = -1;
    if (!unassignedVars->empty()) {
        if (ToulBar2::verbose >= 2)
            cout << "Fast greedy assignment for " << unassignedVars->getSize() << " variables!" << endl;
        Cost currentUb = wcsp->getUb();
        Cost newUb = currentUb;
        int weightedDegree = ToulBar2::weightedDegree;
        ToulBar2::weightedDegree = 0; // do not update weighted degrees inside greedy method
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
        ToulBar2::weightedDegree = weightedDegree;
        if (newUb < currentUb) { /* a better solution has been found */
            wcsp->enforceUb(); /* it will generate a contradiction if lb >= ub */
            wcsp->propagate(); /* it will generate a contradiction if lb >= ub */
        }
        if (unassignedVars->empty()) // a new solution was found and all vars assigned by propagation
            varIndex = -1;
        else {
            // Wrong heuristic guess
            ToulBar2::FullEAC = false;
            varIndex = (this->*varHeurFunc)();
            ToulBar2::FullEAC = true;
            assert(varIndex != -1);
        }
    }
    return varIndex;
}

int Solver::getVarMinDomainDivMaxDegree()
{
    int varIndex = -1;
    Cost worstUnaryCost = MIN_COST;
    double best = MAX_VAL - MIN_VAL;

    for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
        if (ToulBar2::FullEAC) {
            EnumeratedVariable* var = (EnumeratedVariable*)((WCSP*)wcsp)->getVar(*iter);
            if (var->isFullEAC()) {
                continue;
            }
        }
        double heuristic = (double)wcsp->getDomainSize(*iter) / (double)(wcsp->getDegree(*iter) + 1);
        if (varIndex < 0 || heuristic < best - (double)ToulBar2::epsilon * best
            || (heuristic < best + (double)ToulBar2::epsilon * best && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
            best = heuristic;
            varIndex = *iter;
            worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
        }
    }

    if (varIndex == -1) {
        if (ToulBar2::FullEAC) {
            varIndex = greedy(&Solver::getVarMinDomainDivMaxDegree);
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
        if (ToulBar2::FullEAC) {
            EnumeratedVariable* var = (EnumeratedVariable*)((WCSP*)wcsp)->getVar(*iter);
            if (var->isFullEAC()) {
                continue;
            }
        }
        double heuristic = (double)wcsp->getDomainSize(*iter) / (double)(wcsp->getDegree(*iter) + 1);
        if (varIndex < 0 || heuristic < best - (double)ToulBar2::epsilon * best
            || (heuristic < best + (double)ToulBar2::epsilon * best && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
            best = heuristic;
            varIndex = *iter;
            nbties = 1;
            ties[0] = varIndex;
            worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
        } else if (heuristic < best + (double)ToulBar2::epsilon * best && wcsp->getMaxUnaryCost(*iter) == worstUnaryCost) {
            ties[nbties] = *iter;
            nbties++;
        }
    }
    if (nbties > 1)
        return ties[myrand() % nbties];
    else {
        if (varIndex == -1) {
            if (ToulBar2::FullEAC) {
                varIndex = greedy(&Solver::getVarMinDomainDivMaxDegreeRandomized);
            }
        }
        return varIndex;
    }
}

int Solver::getVarMinDomainDivMaxDegreeLastConflict()
{
    if (lastConflictVar != -1 && wcsp->unassigned(lastConflictVar))
        return lastConflictVar;
    int varIndex = -1;
    Cost worstUnaryCost = MIN_COST;
    double best = MAX_VAL - MIN_VAL;
    for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
        if (ToulBar2::FullEAC) {
            EnumeratedVariable* var = (EnumeratedVariable*)((WCSP*)wcsp)->getVar(*iter);
            if (var->isFullEAC()) {
                continue;
            }
        }
        // remove following "+1" when isolated variables are automatically assigned
        double heuristic = (double)wcsp->getDomainSize(*iter) / (double)(wcsp->getDegree(*iter) + 1);
        if (varIndex < 0 || heuristic < best - (double)ToulBar2::epsilon * best
            || (heuristic < best + (double)ToulBar2::epsilon * best && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
            best = heuristic;
            varIndex = *iter;
            worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
        }
    }

    if (varIndex == -1) {
        if (ToulBar2::FullEAC) {
            varIndex = greedy(&Solver::getVarMinDomainDivMaxDegreeLastConflict);
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
        if (ToulBar2::FullEAC) {
            EnumeratedVariable* var = (EnumeratedVariable*)((WCSP*)wcsp)->getVar(*iter);
            if (var->isFullEAC()) {
                continue;
            }
        }
        // remove following "+1" when isolated variables are automatically assigned
        double heuristic = (double)wcsp->getDomainSize(*iter) / (double)(wcsp->getDegree(*iter) + 1);
        if (varIndex < 0 || heuristic < (double)ToulBar2::epsilon * best
            || (heuristic < best + (double)ToulBar2::epsilon * best && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
            best = heuristic;
            varIndex = *iter;
            nbties = 1;
            ties[0] = varIndex;
            worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
            //        } else if ((heuristic < best + (double)ToulBar2::epsilon * best && wcsp->getMaxUnaryCost(*iter) == worstUnaryCost) || ((myrand()%100)==0)) {
        } else if (heuristic < best + (double)ToulBar2::epsilon * best && wcsp->getMaxUnaryCost(*iter) == worstUnaryCost) {
            ties[nbties] = *iter;
            nbties++;
        }
    }
    if (nbties > 1) {
        return ties[myrand() % nbties];
    } else {
        if (varIndex == -1) {
            if (ToulBar2::FullEAC) {
                varIndex = greedy(&Solver::getVarMinDomainDivMaxDegreeLastConflictRandomized);
            }
        }
        return varIndex;
    }
}

int Solver::getVarMinDomainDivMaxWeightedDegree()
{
    int varIndex = -1;
    Cost worstUnaryCost = MIN_COST;
    double best = MAX_VAL - MIN_VAL;

    for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
        if (ToulBar2::FullEAC) {
            EnumeratedVariable* var = (EnumeratedVariable*)((WCSP*)wcsp)->getVar(*iter);
            if (var->isFullEAC()) {
                continue;
            }
        }
        Cost unarymediancost = MIN_COST;
        int domsize = wcsp->getDomainSize(*iter);
        if (ToulBar2::weightedTightness) {
            ValueCost array[domsize];
            wcsp->getEnumDomainAndCost(*iter, array);
            unarymediancost = stochastic_selection<ValueCost>(array, 0, domsize - 1, domsize / 2).cost;
        }
        Long wdeg = wcsp->getWeightedDegree(*iter);
        heuristics[*iter] = max(wdeg, heuristics[*iter]); //cout << "write var " << *iter << " " << wdeg << " " << heuristic[*iter] << endl;
        double heuristic = (double)domsize / (double)(wdeg + 1 + unarymediancost);
        if (varIndex < 0 || heuristic < best - (double)ToulBar2::epsilon * best
            || (heuristic < best + (double)ToulBar2::epsilon * best && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
            best = heuristic;
            varIndex = *iter;
            worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
        }
    }
    if (varIndex == -1) {
        if (ToulBar2::FullEAC) {
            varIndex = greedy(&Solver::getVarMinDomainDivMaxWeightedDegree);
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
        if (ToulBar2::FullEAC) {
            EnumeratedVariable* var = (EnumeratedVariable*)((WCSP*)wcsp)->getVar(*iter);
            if (var->isFullEAC()) {
                continue;
            }
        }
        Cost unarymediancost = MIN_COST;
        int domsize = wcsp->getDomainSize(*iter);
        if (ToulBar2::weightedTightness) {
            ValueCost array[domsize];
            wcsp->getEnumDomainAndCost(*iter, array);
            unarymediancost = stochastic_selection<ValueCost>(array, 0, domsize - 1, domsize / 2).cost;
        }
        Long wdeg = wcsp->getWeightedDegree(*iter);
        heuristics[*iter] = max(wdeg, heuristics[*iter]); //cout << "write var " << *iter << " " << wdeg << " " << heuristic[*iter] << endl;
        double heuristic = (double)domsize / (double)(wdeg + 1 + unarymediancost);
        if (varIndex < 0 || heuristic < best - (double)ToulBar2::epsilon * best
            || (heuristic < best + (double)ToulBar2::epsilon * best && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
            best = heuristic;
            varIndex = *iter;
            nbties = 1;
            ties[0] = varIndex;
            worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
        } else if (heuristic < best + (double)ToulBar2::epsilon * best && wcsp->getMaxUnaryCost(*iter) == worstUnaryCost) {
            ties[nbties] = *iter;
            nbties++;
        }
    }
    if (nbties > 1)
        return ties[myrand() % nbties];
    else {
        if (varIndex == -1) {
            if (ToulBar2::FullEAC) {
                varIndex = greedy(&Solver::getVarMinDomainDivMaxWeightedDegreeRandomized);
            }
        }
        return varIndex;
    }
}

int Solver::getVarMinDomainDivMaxWeightedDegreeLastConflict()
{
    if (lastConflictVar != -1 && wcsp->unassigned(lastConflictVar))
        return lastConflictVar;
    int varIndex = -1;
    Cost worstUnaryCost = MIN_COST;
    double best = MAX_VAL - MIN_VAL;
    for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
        if (ToulBar2::FullEAC) {
            EnumeratedVariable* var = (EnumeratedVariable*)((WCSP*)wcsp)->getVar(*iter);
            if (var->isFullEAC()) {
                continue;
            }
        }
        Cost unarymediancost = MIN_COST;
        int domsize = wcsp->getDomainSize(*iter);
        if (ToulBar2::weightedTightness) {
            ValueCost array[domsize];
            wcsp->getEnumDomainAndCost(*iter, array);
            unarymediancost = stochastic_selection<ValueCost>(array, 0, domsize - 1, domsize / 2).cost;
        }
        // remove following "+1" when isolated variables are automatically assigned
        Long wdeg = wcsp->getWeightedDegree(*iter);
        heuristics[*iter] = max(wdeg, heuristics[*iter]); //cout << "write var " << *iter << " " << wdeg << " " << heuristic[*iter] << endl;
        double heuristic = (double)domsize / (double)(wdeg + 1 + unarymediancost);
        // double heuristic = 1. / (double) (wcsp->getMaxUnaryCost(*iter) + 1);
        if ((varIndex < 0)
            || (heuristic < best - (double)ToulBar2::epsilon * best)
            || (heuristic < best + (double)ToulBar2::epsilon * best && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
            best = heuristic;
            varIndex = *iter;
            worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
        }
    }

    if (varIndex == -1) {
        if (ToulBar2::FullEAC) {
            varIndex = greedy(&Solver::getVarMinDomainDivMaxWeightedDegreeLastConflict);
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
        if (ToulBar2::FullEAC) {
            EnumeratedVariable* var = (EnumeratedVariable*)((WCSP*)wcsp)->getVar(*iter);
            if (var->isFullEAC()) {
                continue;
            }
        }
        Cost unarymediancost = MIN_COST;
        int domsize = wcsp->getDomainSize(*iter);
        if (ToulBar2::weightedTightness) {
            ValueCost array[domsize];
            wcsp->getEnumDomainAndCost(*iter, array);
            unarymediancost = stochastic_selection<ValueCost>(array, 0, domsize - 1, domsize / 2).cost;
        }
        // remove following "+1" when isolated variables are automatically assigned
        Long wdeg = wcsp->getWeightedDegree(*iter);
        heuristics[*iter] = max(wdeg, heuristics[*iter]); //cout << "write var " << *iter << " " << wdeg << " " << heuristic[*iter] << endl;
        double heuristic = (double)domsize / (double)(wdeg + 1 + unarymediancost);
        if ((varIndex < 0)
            || (heuristic < best - (double)ToulBar2::epsilon * best)
            || (heuristic < best + (double)ToulBar2::epsilon * best && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
            best = heuristic;
            varIndex = *iter;
            nbties = 1;
            ties[0] = varIndex;
            worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
            //       } else if ((heuristic < best + (double)ToulBar2::epsilon * best && wcsp->getMaxUnaryCost(*iter) == worstUnaryCost) || ((myrand()%100)==0)) {
        } else if (heuristic < best + (double)ToulBar2::epsilon * best && wcsp->getMaxUnaryCost(*iter) == worstUnaryCost) {
            ties[nbties] = *iter;
            nbties++;
        }
    }
    if (nbties > 1) {
        return ties[myrand() % nbties];
    } else {
        if (varIndex == -1) {
            if (ToulBar2::FullEAC) {
                varIndex = greedy(&Solver::getVarMinDomainDivMaxWeightedDegreeLastConflictRandomized);
            }
        }
        return varIndex;
    }
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

Cost Solver::logZCurrentEstimate()
{
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
    return newCost;
}

/// \brief Enforce WCSP upper-bound and backtrack if ub <= lb or in the case of probabilistic inference if the contribution is too small
void Solver::enforceUb()
{
    wcsp->enforceUb();
    if (ToulBar2::isZ && ToulBar2::logepsilon != -numeric_limits<TLogProb>::infinity()) {
        TLogProb newlogU = wcsp->LogSumExp(ToulBar2::logU, logZCurrentEstimate());
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
            string pbname = to_string("problem") + to_string(nbNodes) + to_string(".wcsp");
            ofstream pb(pbname.c_str());
            wcsp->dump(pb);
            cout << " #" << nbNodes;
        }
        cout << "[" << Store::getDepth() << "," << wcsp->getLb() << "," << wcsp->getUb() << "," << wcsp->getDomainSizeSum();
        if (wcsp->getTreeDec()) {
            cout << ",C" << wcsp->getTreeDec()->getCurrentCluster()->getId();
            if (ToulBar2::heuristicFreedom) {
                cout << "," << wcsp->getTreeDec()->getCurrentCluster()->getFreedom();
            }
        }
        string valname = wcsp->getValueName(varIndex, value);
        if (valname.size() > 0) {
            cout << "] Try " << wcsp->getName(varIndex) << " >= " << value << " (" << valname << ") (s:" << wcsp->getSupport(varIndex) << ")" << endl;
        } else {
            cout << "] Try " << wcsp->getName(varIndex) << " >= " << value << " (s:" << wcsp->getSupport(varIndex) << ")" << endl;
        }
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
            string pbname = to_string("problem") + to_string(nbNodes) + to_string(".wcsp");
            ofstream pb(pbname.c_str());
            wcsp->dump(pb);
            cout << " #" << nbNodes;
        }
        cout << "[" << Store::getDepth() << "," << wcsp->getLb() << "," << wcsp->getUb() << "," << wcsp->getDomainSizeSum();
        if (wcsp->getTreeDec()) {
            cout << ",C" << wcsp->getTreeDec()->getCurrentCluster()->getId();
            if (ToulBar2::heuristicFreedom) {
                cout << "," << wcsp->getTreeDec()->getCurrentCluster()->getFreedom();
            }
        }
        string valname = wcsp->getValueName(varIndex, value);
        if (valname.size() > 0) {
            cout << "] Try " << wcsp->getName(varIndex) << " <= " << value << " (" << valname << ") (s:" << wcsp->getSupport(varIndex) << ")" << endl;
        } else {
            cout << "] Try " << wcsp->getName(varIndex) << " <= " << value << " (s:" << wcsp->getSupport(varIndex) << ")" << endl;
        }
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
                Cost delta = wcsp->getTreeDec()->getCurrentCluster()->getCurrentDeltaUb();
                if (wcsp->getTreeDec()->getCurrentCluster()->open->size() > 0)
                    cout << " [" << wcsp->getTreeDec()->getCurrentCluster()->open->getLb(delta) << "," << wcsp->getUb() << "]/" << wcsp->getTreeDec()->getCurrentCluster()->open->size() << "/" << wcsp->getTreeDec()->getCurrentCluster()->cp->size() << " " << (100. * (wcsp->getUb() - wcsp->getTreeDec()->getCurrentCluster()->open->getLb(delta)) / wcsp->getUb()) << "%";
            } else {
                if (open->size() > 0)
                    cout << " [" << open->getLb() << "," << wcsp->getUb() << "]/" << open->size() << "/" << cp->size() << "/" << nbNodes << " " << (100. * (wcsp->getUb() - open->getLb()) / wcsp->getUb()) << "%";
            }
        } else if (ToulBar2::vnsKmax > 0) {
            cout << " " << ToulBar2::vnsKcur << " " << ToulBar2::vnsLDScur;
        }
        cout << " " << Exp(((int)(*((StoreInt*)searchSize))) / 10e3);
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
            string pbname = to_string("problem") + to_string(nbNodes) + to_string(".wcsp");
            ofstream pb(pbname.c_str());
            wcsp->dump(pb);
            cout << " #" << nbNodes;
        }
        cout << "[" << Store::getDepth() << "," << wcsp->getLb() << "," << wcsp->getUb() << "," << wcsp->getDomainSizeSum();
        if (wcsp->getTreeDec()) {
            cout << ",C" << wcsp->getTreeDec()->getCurrentCluster()->getId();
            if (ToulBar2::heuristicFreedom) {
                cout << "," << wcsp->getTreeDec()->getCurrentCluster()->getFreedom();
            }
        }
        string valname = wcsp->getValueName(varIndex, value);
        if (valname.size() > 0) {
            cout << "] Try " << wcsp->getName(varIndex) << " == " << value << " (" << valname << ")" << endl;
        } else {
            cout << "] Try " << wcsp->getName(varIndex) << " == " << value << endl;
        }
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
            string pbname = to_string("problem") + to_string(nbNodes) + to_string(".wcsp");
            ofstream pb(pbname.c_str());
            wcsp->dump(pb);
            cout << " #" << nbNodes;
        }
        cout << "[" << Store::getDepth() << "," << wcsp->getLb() << "," << wcsp->getUb() << "," << wcsp->getDomainSizeSum();
        if (wcsp->getTreeDec()) {
            cout << ",C" << wcsp->getTreeDec()->getCurrentCluster()->getId();
            if (ToulBar2::heuristicFreedom) {
                cout << "," << wcsp->getTreeDec()->getCurrentCluster()->getFreedom();
            }
        }
        string valname = wcsp->getValueName(varIndex, value);
        if (valname.size() > 0) {
            cout << "] Try " << wcsp->getName(varIndex) << " != " << value << " (" << valname << ")" << endl;
        } else {
            cout << "] Try " << wcsp->getName(varIndex) << " != " << value << endl;
        }
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
            string pbname = to_string("problem") + to_string(nbNodes) + to_string(".wcsp");
            ofstream pb(pbname.c_str());
            wcsp->dump(pb);
            cout << " #" << nbNodes;
        }
        cout << "[" << Store::getDepth() << "," << wcsp->getLb() << "," << wcsp->getUb() << "," << wcsp->getDomainSizeSum();
        if (wcsp->getTreeDec()) {
            cout << ",C" << wcsp->getTreeDec()->getCurrentCluster()->getId();
            if (ToulBar2::heuristicFreedom) {
                cout << "," << wcsp->getTreeDec()->getCurrentCluster()->getFreedom();
            }
        }
        cout << "] Try " << wcsp->getName(varIndex) << " !=";
        for (int i = first; i <= last; i++) {
            string valname = wcsp->getValueName(varIndex, array[i].value);
            if (valname.size() > 0) {
                cout << " " << array[i].value << " (" << valname << ")";
            } else {
                cout << " " << array[i].value;
            }
        }
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
            cout << "Optimality gap: [" << std::fixed << std::setprecision(ToulBar2::decimalPoint) << Dglb << ", " << Dgub << "] " << std::setprecision(DECIMAL_POINT) << (100. * (Dgub - Dglb)) / max(fabsl(Dglb), fabsl(Dgub)) << " % (" << nbBacktracks << " backtracks, " << nbNodes << " nodes, " << ((ToulBar2::parallel) ? (realTime() - ToulBar2::startRealTime) : (cpuTime() - ToulBar2::startCpuTime)) << " seconds)" << endl;
            cout.flags(f);
        }
    }
}

void Solver::binaryChoicePoint(int varIndex, Value value, Cost lb)
{
    assert(wcsp->unassigned(varIndex));
    assert(wcsp->canbe(varIndex, value));
    if (ToulBar2::interrupted && !ToulBar2::isZ)
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
    if (ToulBar2::parallel && ((nbBacktracks % 128) == 0) && MPI_interrupted())
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
    if (ToulBar2::interrupted && ToulBar2::isZ) {
        ToulBar2::logU = wcsp->LogSumExp(ToulBar2::logU, logZCurrentEstimate());
    } else {
        if (nbBacktracks >= hbfsLimit) {
            assert(ToulBar2::hbfs);
            addOpenNode(*cp, *open, MAX(lb, wcsp->getLb()));
#ifdef OPENMPI
            if (ToulBar2::parallel && ToulBar2::burst && world.rank() != MASTER) {
                vector<Value> emptySol;
                Work work2(*cp, *open, nbNodes - initWorkerNbNodes, nbBacktracks - initWorkerNbBacktracks, wcsp->getNbDEE() - initWorkerNbDEE, nbRecomputationNodes - initWorkerNbRecomputationNodes, MIN_COST, MAX_COST, emptySol);
                if (ToulBar2::verbose >= 1)
                    cout << ">>> worker " << world.rank() << " send open-node message to master " << work2 << endl;
                double beginWaiting = realTime();
                mpi::request req = world.isend(MASTER, WORKTAG, work2); // non-blocking send to master
                while (!req.test().is_initialized() && !MPI_interrupted())
                    ;
                hbfsWaitingTime += realTime() - beginWaiting;
                assert(open->empty());
                initWorkerNbNodes = nbNodes;
                initWorkerNbBacktracks = nbBacktracks;
                initWorkerNbDEE = wcsp->getNbDEE();
                initWorkerNbRecomputationNodes = nbRecomputationNodes;
            }
#endif
        } else {
            recursiveSolve(lb);
        }
    }
}

void Solver::binaryChoicePointLDS(int varIndex, Value value, int discrepancy)
{
    assert(wcsp->unassigned(varIndex));
    assert(wcsp->canbe(varIndex, value));
    if (ToulBar2::interrupted && !ToulBar2::isZ)
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
        if (ToulBar2::parallel && ((nbBacktracks % 128) == 0) && MPI_interrupted())
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
        if (ToulBar2::interrupted && ToulBar2::isZ) {
            ToulBar2::logU = wcsp->LogSumExp(ToulBar2::logU, logZCurrentEstimate());
        } else {
            recursiveSolveLDS(discrepancy);
        }
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
    if (ToulBar2::parallel && ((nbBacktracks % 128) == 0) && MPI_interrupted())
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
    // ValueCost* sorted = new ValueCost [size];
    wcsp->getEnumDomainAndCost(varIndex, sorted);
    qsort(sorted, size, sizeof(ValueCost), cmpValueCost);
    for (int v = 0; wcsp->getLb() < wcsp->getUb() && v < size; v++) {
        if (ToulBar2::interrupted && !ToulBar2::isZ)
            throw TimeOut();
        int storedepth = Store::getDepth();
        try {
            Store::store();
            assign(varIndex, sorted[v].value);
            if (ToulBar2::interrupted && ToulBar2::isZ) {
                ToulBar2::logU = wcsp->LogSumExp(ToulBar2::logU, logZCurrentEstimate());
            } else {
                recursiveSolve(lb);
            }
        } catch (const Contradiction&) {
            wcsp->whenContradiction();
        }
        Store::restore(storedepth);
    }
    // delete [] sorted;
    enforceUb();
    nbBacktracks++;
    if (nbBacktracks > nbBacktracksLimit)
        throw NbBacktracksOut();
#ifdef OPENMPI
    if (ToulBar2::parallel && ((nbBacktracks % 128) == 0) && MPI_interrupted())
        throw TimeOut();
#endif
}

void Solver::narySortedChoicePointLDS(int varIndex, int discrepancy)
{
    assert(wcsp->enumerated(varIndex));
    int size = wcsp->getDomainSize(varIndex);
    ValueCost sorted[size];
    // ValueCost* sorted = new ValueCost [size];
    wcsp->getEnumDomainAndCost(varIndex, sorted);
    qsort(sorted, size, sizeof(ValueCost), cmpValueCost);
    if (discrepancy < size - 1)
        ToulBar2::limited = true;
    for (int v = min(size - 1, discrepancy); wcsp->getLb() < wcsp->getUb() && v >= 0; v--) {
        if (ToulBar2::interrupted && !ToulBar2::isZ)
            throw TimeOut();
        int storedepth = Store::getDepth();
        try {
            Store::store();
            assign(varIndex, sorted[v].value);
            if (ToulBar2::interrupted && ToulBar2::isZ) {
                ToulBar2::logU = wcsp->LogSumExp(ToulBar2::logU, logZCurrentEstimate());
            } else {
                recursiveSolveLDS(discrepancy - v);
            }
        } catch (const Contradiction&) {
            wcsp->whenContradiction();
        }
        Store::restore(storedepth);
    }
    // delete [] sorted;
    enforceUb();
    nbBacktracks++;
    if (nbBacktracks > nbBacktracksLimit)
        throw NbBacktracksOut();
#ifdef OPENMPI
    if (ToulBar2::parallel && ((nbBacktracks % 128) == 0) && MPI_interrupted())
        throw TimeOut();
#endif
}

void Solver::singletonConsistency(int restricted)
{
    bool deadend;
    bool done = false;
    int vac = ToulBar2::vac;
    if (ToulBar2::vac) {
        ToulBar2::vac = Store::getDepth() + 2; // make sure VAC is performed if requested
    }
    vector<int> revelimorder(wcsp->numberOfVariables(), -1);
    for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
        revelimorder[wcsp->getDACOrder(i)] = i;
    }
    while (!done && restricted) {
        done = true;
        unsigned int revelimpos = 0;
        while (revelimpos < wcsp->numberOfVariables() && restricted) {
            assert(revelimorder[revelimpos] >= 0 && revelimorder[revelimpos] < (int)wcsp->numberOfVariables());
            unsigned int varIndex = revelimorder[revelimpos++];
            unsigned int size = wcsp->getDomainSize(varIndex);
            if (size > 1 && (ToulBar2::nbDecisionVars <= 0 || varIndex < (unsigned int)ToulBar2::nbDecisionVars)) {
                restricted--;
                ValueCost sorted[size];
                // ValueCost* sorted = new ValueCost [size];
                // wcsp->iniSingleton(); //Warning! constructive disjunction is not compatible with variable elimination
                wcsp->getEnumDomainAndCost(varIndex, sorted);
                qsort(sorted, size, sizeof(ValueCost), cmpValueCost);
                Cost initlb = wcsp->getLb();
                Cost minlambda = MAX_COST;
                for (int a = size - 1; a >= 0; --a) {
                    if (wcsp->canbe(varIndex, sorted[a].value)) {
                        deadend = false;
                        int storedepth = Store::getDepth();
                        try {
                            Store::store();
                            wcsp->assign(varIndex, sorted[a].value);
                            wcsp->propagate();
                            if (wcsp->getLb() - initlb < minlambda) {
                                minlambda = wcsp->getLb() - initlb;
                            }
                            if (ToulBar2::verbose >= 1 && wcsp->getLb() - initlb > sorted[a].cost) {
                                cout << "singleton consistency may increase unary cost of variable " << wcsp->getName(varIndex) << " value " << sorted[a].value << " from " << sorted[a].cost  << " to " << wcsp->getLb() - initlb << endl;
                            }
                        } catch (const Contradiction&) {
                            wcsp->whenContradiction();
                            deadend = true;
                            done = false;
                        }
                        Store::restore(storedepth);
                        // wcsp->updateSingleton();
                        // cout << "(" << varIndex << "," << a <<  ")" << endl;
                        if (deadend) {
                            wcsp->remove(varIndex, sorted[a].value);
                            wcsp->propagate();
                            if (ToulBar2::verbose >= 0) {
                                cout << ".";
                                flush(cout);
                            }
                            // WARNING!!! can we stop if the variable is assigned, what about removeSingleton after???
                        }
                    }
                }
//                if (minlambda > MIN_COST) {
//                    if (ToulBar2::verbose >= 0) {
//                        cout << "singleton consistency may increase lower bound by " << minlambda << " from variable " << wcsp->getName(varIndex) << endl;
//                    }
//                    initlb = wcsp->getLb();
//                    deadend = false;
//                    unsigned int initsize = wcsp->getDomainInitSize(varIndex);
//                    vector<Cost> costs(initsize, minlambda);
//                    Value support = wcsp->getSupport(varIndex);
//                    unsigned int supportIndex = wcsp->toIndex(varIndex, support);
//                    int storedepth = Store::getDepth();
//                    bool stealing = false;
//                    try {
//                        Store::store();
//                        costs[supportIndex] = MIN_COST;
//                        wcsp->postUnaryConstraint(varIndex, costs);
//                        wcsp->propagate();
//                        assert(wcsp->getUnaryCost(varIndex, support) == MIN_COST);
//                        for (unsigned int a = 0; a < size; a++) {
//                            Value value = sorted[a].value;
//                            Cost ucost = wcsp->getUnaryCost(varIndex, value);
//                            if (value != support && wcsp->canbe(varIndex, value) && ucost < minlambda && minlambda - ucost > sorted[a].cost) {
//                                if (ToulBar2::verbose >= 0) {
//                                    cout << "singleton consistency has stealing cost " << (minlambda - ucost - sorted[a].cost) << " from value " << value << endl;
//                                }
//                                stealing = true;
//                            }
//                        }
//                        if (ToulBar2::verbose >= 0 && wcsp->getLb() > initlb && !stealing) {
//                            cout << "singleton consistency should increase unary cost of variable " << wcsp->getName(varIndex) << " value " << support << " from 0 to " << (wcsp->getLb() - initlb) << endl;
//                        }
//                    } catch (const Contradiction&) {
//                        wcsp->whenContradiction();
//                        deadend = true;
//                    }
//                    Store::restore(storedepth);
//                    // wcsp->updateSingleton();
//                    // cout << "(" << varIndex << "," << a <<  ")" << endl;
//                    if (!deadend) {
//                        if (ToulBar2::verbose >= 0) {
//                            cout << "singleton consistency will increase unary cost of variable " << wcsp->getName(varIndex) << endl;
//                        }
//                        // WARNING!!! can we stop if the variable is assigned, what about removeSingleton after???
//                    }
//                }
                // wcsp->removeSingleton();
                // delete [] sorted;
            }
        }
    }
    ToulBar2::vac = vac;
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
        ToulBar2::logZ = wcsp->LogSumExp(ToulBar2::logZ, (Cost)(wcsp->getLb() + wcsp->getNegativeLb()));
        if (ToulBar2::debug && (nbBacktracks % 10000LL) == 0 && ToulBar2::logepsilon != -numeric_limits<TLogProb>::infinity())
            cout << (ToulBar2::logZ + ToulBar2::markov_log - (TLogProb)wcsp->getMaxDomainSize() * Exp10(-(TLogProb)ToulBar2::resolution)) << " , " << (wcsp->LogSumExp(ToulBar2::logZ, ToulBar2::logU) + ToulBar2::markov_log + (TLogProb)wcsp->getMaxDomainSize() * Exp10(-(TLogProb)ToulBar2::resolution)) << " in " << ((ToulBar2::parallel) ? (realTime() - ToulBar2::startRealTime) : (cpuTime() - ToulBar2::startCpuTime)) << " seconds" << endl;
    }
    if ((!ToulBar2::allSolutions && !ToulBar2::isZ) || ToulBar2::debug >= 2) {
        if (ToulBar2::verbose >= 0 || (!ToulBar2::parallel && ToulBar2::showSolutions)) {
            if (ToulBar2::haplotype)
                cout << "***New solution: " << wcsp->getLb() << " log10like: " << ToulBar2::haplotype->Cost2LogProb(wcsp->getLb()) / Log(10.) << " logProb: " << ToulBar2::haplotype->Cost2LogProb(wcsp->getLb()) << " (" << nbBacktracks << " backtracks, " << nbNodes << " nodes, depth " << Store::getDepth() << ", " << ((ToulBar2::parallel) ? (realTime() - ToulBar2::startRealTime) : (cpuTime() - ToulBar2::startCpuTime)) << " seconds)" << endl;
            else if (!ToulBar2::bayesian)
                cout << "New solution: " << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->getDDualBound() << std::setprecision(DECIMAL_POINT) << " (" << nbBacktracks << " backtracks, " << nbNodes << " nodes, depth " << Store::getDepth() << ", " << ((ToulBar2::parallel) ? (realTime() - ToulBar2::startRealTime) : (cpuTime() - ToulBar2::startCpuTime)) << " seconds)" << endl;
            else
                cout << "New solution: " << wcsp->getLb() << " energy: " << -(wcsp->Cost2LogProb(wcsp->getLb() + wcsp->getNegativeLb()) + ToulBar2::markov_log) << " prob: " << std::scientific << wcsp->Cost2Prob(wcsp->getLb() + wcsp->getNegativeLb()) * Exp(ToulBar2::markov_log) << std::fixed << " (" << nbBacktracks << " backtracks, " << nbNodes << " nodes, depth " << Store::getDepth() << ", " << ((ToulBar2::parallel) ? (realTime() - ToulBar2::startRealTime) : (cpuTime() - ToulBar2::startCpuTime)) << " seconds)" << endl;
        }
    }

    wcsp->restoreSolution(); // update all variables to be in the current assignment (necessary when some variables are eliminated)
    if (!ToulBar2::isZ)
        wcsp->setSolution(wcsp->getLb()); // take current assignment and put it in solution (stl c++ map)

#ifdef OPENMPI
    if (!ToulBar2::parallel || world.rank() == MASTER) {
#endif
        if (ToulBar2::showSolutions) {
            if (ToulBar2::verbose >= 2)
                cout << *wcsp << endl;

            if (ToulBar2::allSolutions) {
                cout << std::setprecision(0) << nbSol << " solution(" << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->getDDualBound() << std::setprecision(DECIMAL_POINT) << "): ";
            }

            for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
                if (ToulBar2::pedigree) {
                    cout << " ";
                    cout << wcsp->getName(i) << ":";
                    ToulBar2::pedigree->printGenotype(cout, wcsp->getValue(i));
                } else if (ToulBar2::haplotype) {
                    cout << " ";
                    ToulBar2::haplotype->printHaplotype(cout, wcsp->getValue(i), i);
                } else if (wcsp->enumerated(i) && ((EnumeratedVariable*)((WCSP*)wcsp)->getVar(i))->isValueNames()) {
                    EnumeratedVariable* myvar = (EnumeratedVariable*)((WCSP*)wcsp)->getVar(i);
                    Value myvalue = ((ToulBar2::sortDomains && ToulBar2::sortedDomains.find(i) != ToulBar2::sortedDomains.end()) ? ToulBar2::sortedDomains[i][myvar->toIndex(myvar->getValue())].value : myvar->getValue());
                    string valuelabel = myvar->getValueName(myvar->toIndex(myvalue));
                    string varlabel = myvar->getName();

                    if (ToulBar2::showHidden || (varlabel.rfind(HIDDEN_VAR_TAG, 0) != 0)) {
                        switch (ToulBar2::showSolutions) {
                        case 1:
                            cout << " ";
                            cout << myvalue;
                            break;
                        case 2:
                            cout << " ";
                            cout << valuelabel;
                            break;
                        case 3:
                            cout << " ";
                            cout << varlabel << "=" << valuelabel;
                            break;
                        default:
                            break;
                        }
                    }
                } else {
                    cout << " ";
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
        if (!ToulBar2::uaieval && ToulBar2::writeSolution) {
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
            cout << "o " << std::fixed << std::setprecision(0) << wcsp->getDDualBound() << std::setprecision(DECIMAL_POINT) << endl; //" ";
            ((WCSP*)wcsp)->solution_XML(false);
        }
        if (ToulBar2::maxsateval) {
            cout << "o " << wcsp->getLb() << endl;
        }
        if (ToulBar2::uaieval && !ToulBar2::isZ) {
            ((WCSP*)wcsp)->solution_UAI(wcsp->getLb());
        }
#ifdef OPENMPI
    }
#endif

    if (ToulBar2::newsolution)
        (*ToulBar2::newsolution)(wcsp->getIndex(), wcsp->getSolver());

    if (ToulBar2::restart == 0 && !ToulBar2::lds && !ToulBar2::isZ)
        throw NbBacktracksOut();
    if (ToulBar2::allSolutions && nbSol >= ToulBar2::allSolutions)
        throw NbSolutionsOut();
    if (ToulBar2::divNbSol > 1 && wcsp->getLb() <= prevDivSolutionCost)
        throw DivSolutionOut();
#ifdef OPENMPI
    if (ToulBar2::parallel && ToulBar2::searchMethod == DFBB && ToulBar2::burst && world.rank() != MASTER) { // HBFS may be turn-off due to open list memory-out and switch to DFS
        Cost newWorkerUb = wcsp->getSolutionCost();
        vector<Value> workerSol = wcsp->getSolution();
        assert(open && open->empty());
        Work work2(*cp, *open, nbNodes - initWorkerNbNodes, nbBacktracks - initWorkerNbBacktracks, wcsp->getNbDEE() - initWorkerNbDEE, nbRecomputationNodes - initWorkerNbRecomputationNodes, MIN_COST, newWorkerUb, workerSol);
        if (ToulBar2::verbose >= 1)
            cout << ">>> worker " << world.rank() << " send solution message to master " << work2 << endl;
        double beginWaiting = realTime();
        mpi::request req = world.isend(MASTER, WORKTAG, work2); // non-blocking send to master
        while (!req.test().is_initialized() && !MPI_interrupted())
            ;
        hbfsWaitingTime += realTime() - beginWaiting;
        initWorkerNbNodes = nbNodes;
        initWorkerNbBacktracks = nbBacktracks;
        initWorkerNbDEE = wcsp->getNbDEE();
        initWorkerNbRecomputationNodes = nbRecomputationNodes;
    }
#endif
}

void Solver::recursiveSolve(Cost lb)
{
    assert(ToulBar2::nbDecisionVars > 0 || numberOfUnassignedVariables() == (int)getWCSP()->numberOfUnassignedVariables());
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
        *((StoreInt*)searchSize) += ((int)(10e3 * Log(wcsp->getDomainSize(varIndex))));
        if (ToulBar2::bep)
            scheduleOrPostpone(varIndex);
        else if (wcsp->enumerated(varIndex)) {
            if (ToulBar2::binaryBranching) {
                assert(wcsp->canbe(varIndex, wcsp->getSupport(varIndex)));
                // Reuse last solution found if available and if not valid look at a (bi-objective) bounding constraint support if any
                Value bestval = ((ToulBar2::verifyOpt) ? (wcsp->getSup(varIndex) + 1) : wcsp->getBestValue(varIndex));
                if (wcsp->cannotbe(varIndex, bestval) && ToulBar2::bisupport != 0. && WeightedCSPConstraint::WeightedCSPConstraints.size() > 0) {
                    WeightedCSPConstraint* objective2 = WeightedCSPConstraint::WeightedCSPConstraints.begin()->second; // TODO: choose the WeightedCSPConstraints with the highest weighted degree
                    int sign = 1;
                    bool mingap = true;
                    if (ToulBar2::bisupport < 0. || (wcsp->getUnaryCost(varIndex, objective2->getSupport(varIndex, sign, mingap)) < (Double)ToulBar2::bisupport * objective2->getUnaryCost(varIndex, wcsp->getSupport(varIndex), sign))) {
                        if (ToulBar2::bisupport < 0.) {
                            switch ((int)(-ToulBar2::bisupport)) {
                            case BISUPPORT_HEUR_LB:
                                sign = 1;
                                break;
                            case BISUPPORT_HEUR_UB:
                                sign = -1;
                                break;
                            case BISUPPORT_HEUR_MINGAP:
                                sign = 0;
                                mingap = true;
                                break;
                            case BISUPPORT_HEUR_MAXGAP:
                                sign = 0;
                                mingap = false;
                                break;
                            default:
                                cerr << "Unknown bisupport heuristic! " << ToulBar2::bisupport << endl;
                                throw BadConfiguration();
                            }
                        }
                        bestval = objective2->getSupport(varIndex, sign, mingap);
                    }
                }
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
#ifdef OPENMPI
        if (ToulBar2::parallel && (!cluster || cluster == wcsp->getTreeDec()->getRoot())) {
            world.barrier(); /* IMPORTANT */
            ToulBar2::startRealTimeAfterPreProcessing = realTime();
            hbfsWaitingTime = 0.;
            if (world.rank() == MASTER) {
                return hybridSolveMaster(cluster, clb, cub);
            } else {
                return hybridSolveWorker(cluster, clb, cub);
            }
        }
#endif
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
                delta = cluster->getCurrentDeltaUb();
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
        Cost lastUb = MAX_COST;
        Long iterBFS = 0;
        Long sortBFS = ToulBar2::sortBFS;
        while (clb < cub && !open_->finished() && (!cluster || (clb == initiallb && cub == initialub && nbBacktracks <= cluster->hbfsGlobalLimit))) {
            iterBFS++;
            if (cluster) {
                cluster->hbfsLimit = ((ToulBar2::hbfs > 0) ? (cluster->nbBacktracks + ToulBar2::hbfs) : LONGLONG_MAX);
                assert(wcsp->getTreeDec()->getCurrentCluster() == cluster);
                wcsp->setUb(cub);
                assert(cluster->isActive());
                assert(cluster->getLbRec() == wcsp->getLb());
            } else {
                hbfsLimit = ((ToulBar2::hbfs > 0) ? (nbBacktracks + ToulBar2::hbfs) : LONGLONG_MAX);
            }
            if (cub < lastUb && (!cluster || cluster == wcsp->getTreeDec()->getRoot())) {
                // enforces best upper bound at current depth to avoid redoing its propagation work when restoring each stored open node
                int saveElimDegree = ToulBar2::elimDegree_;
                int saveDEE = ToulBar2::DEE_;
                ToulBar2::elimDegree_ = -1; // variable elimination or dead-end elimination may be incompatible with stored open nodes!!!
                ToulBar2::DEE_ = 0;
                wcsp->enforceUb();
                wcsp->propagate();
                lastUb = cub;
                // reactivate on-the-fly variable elimination and dead-end elimination (undone at this level but at level+1)
                ToulBar2::elimDegree_ = saveElimDegree;
                ToulBar2::DEE_ = saveDEE;
                for (int i = wcsp->numberOfVariables() - 1; i >= 0; i--) {
                    if (wcsp->unassigned(i) && (!cluster || wcsp->getTreeDec()->getCluster(((WCSP*)wcsp)->getVar(i)->getCluster())->isActive())) {
                        ((WCSP*)wcsp)->getVar(i)->queueEliminate();
                        ((WCSP*)wcsp)->getVar(i)->queueDEE();
                    }
                }
                if (ToulBar2::sortBFS) {
                    //re-sort open nodes at every new solution
                    if (ToulBar2::verbose >= 1) {
                        cout << "Sort open nodes using heuristics.." << endl;
                    }
                    OpenList *resort = new OpenList;
                    for (auto iter = open_->begin(); iter != open_->end(); ++iter) {
                        resort->push(*iter);
                    }
                    delete open_;
                    open_ = resort;
                    open = resort;
                }
            } else if (ToulBar2::sortBFS && iterBFS >= sortBFS) {
                //re-sort open nodes
                if (ToulBar2::verbose >= 1) {
                    cout << "Sort open nodes using heuristics.." << endl;
                }
                OpenList *resort = new OpenList;
                for (auto iter = open_->begin(); iter != open_->end(); ++iter) {
                   resort->push(*iter);
                }
                delete open_;
                open_ = resort;
                open = resort;
                while (sortBFS <= iterBFS) {
                    sortBFS *= 2;
                }
                if (ToulBar2::verbose >= 1) {
                    cout << "Limit before next sorting: " << sortBFS << " (" << nbNodes << ")" << endl;
                }
            }
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
                    if (ToulBar2::interrupted && ToulBar2::isZ) {
                        ToulBar2::logU = wcsp->LogSumExp(ToulBar2::logU, logZCurrentEstimate());
                    } else {
                        recursiveSolve(bestlb);
                    }
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
            if (ToulBar2::eps && open_->size() >= static_cast<std::size_t>(ToulBar2::eps)) {
                epsDumpSubProblems(*cp_, *open_);
                ToulBar2::interrupted = true;
                throw TimeOut();
            }
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
    assert(ToulBar2::bilevel || clb <= cub);
    return make_pair(clb, cub);
}

#ifdef OPENMPI
pair<Cost, Cost> Solver::hybridSolveMaster(Cluster* cluster, Cost clb, Cost cub)
{
    if (ToulBar2::verbose >= 1) {
        if (cluster)
            cout << "hybridSolveMaster C" << cluster->getId() << " " << clb << " " << cub << endl;
        else
            cout << "hybridSolveMaster " << clb << " " << cub << endl;
    }
    assert(clb < cub);
    assert(wcsp->getUb() == cub);
    assert(wcsp->getLb() <= clb);
    assert(ToulBar2::hbfs);
    CPStore* cp_ = NULL;
    OpenList* open_ = NULL;
    if (cluster) {
        // BFS with BTD on current cluster (can be root or not)
        assert(cluster->cp);
        cp_ = cluster->cp;
        assert(cluster == wcsp->getTreeDec()->getRoot());
        if (!cluster->open) {
            cluster->open = new OpenList();
            assert(idleQ.size() == (size_t)(world.size() - 1));
        }
        cluster->setUb(cub); // global problem upper bound
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
        assert(idleQ.size() == (size_t)(world.size() - 1));
    }
    cp_->store();
    if (open_->size() == 0 || (cluster && (clb >= open_->getClosedNodesLb() || cub > open_->getUb()))) { // start a new list of open nodes if needed
        if (open_->size() == 0 && (!cluster || cluster->getNbVars() > 0))
            nbHybridNew++;
        // reinitialize current open list and insert empty node
        *open_ = OpenList(MAX(MIN_COST, cub), MAX(MIN_COST, cub));
        addOpenNode(*cp_, *open_, clb);
        // wait for messages from workers in order to reinitialize the idleQ //TODO: send a message to ask to stop current search but not die
        while (idleQ.size() < (size_t)(world.size() - 1)) {
            Work work; // dummy work
            mpi::status status = world.recv(mpi::any_source, mpi::any_tag, work); // blocking recv to wait for matching messages from any worker
            if (status.tag() == IDLETAG) {
                activeWork.erase(status.source());
                idleQ.push(status.source());
            }
        }
    } else if (!cluster || cluster->getNbVars() > 0)
        nbHybridContinue++;
    if (!cluster || cluster->getNbVars() > 0)
        nbHybrid++; // do not count empty root cluster
    if (cluster)
        cluster->hbfsGlobalLimit = ((ToulBar2::hbfsGlobalLimit > 0) ? (nbBacktracks + ToulBar2::hbfsGlobalLimit) : LONGLONG_MAX);
    Cost initiallb = clb;
    Cost initialub = cub;
    open_->updateUb(cub);
    clb = MAX(clb, open_->getLb());
    if (ToulBar2::verbose >= 1 && cluster)
        cout << "hybridSolve-2 C" << cluster->getId() << " " << clb << " " << cub << " " << open_->size() << " " << open_->top().getCost() << " " << open_->getClosedNodesLb() << " " << open_->getUb() << endl;
    while (clb < cub && (!open_->finished() || !activeWork.empty()) && (!cluster || (clb == initiallb && cub == initialub && nbBacktracks <= cluster->hbfsGlobalLimit))) {
        if (cluster) {
            assert(wcsp->getTreeDec()->getCurrentCluster() == cluster);
            wcsp->setUb(cub);
            assert(cluster->isActive());
            assert(cluster->getLbRec() == wcsp->getLb());
        }
        Cost initub = wcsp->getUb();
        // loop to distribute jobs to workers
        vector<mpi::request> reqs;
        while (!open_->finished() && !idleQ.empty()) { // while there is work to do and workers to do it
            int worker = idleQ.front(); // get the first worker in the queue
            vector<Value> masterSol;
            Cost masterUb = wcsp->getSolutionCost();
            if (masterUb < MAX_COST && (bestsolWork.find(worker) == bestsolWork.end() || masterUb < bestsolWork[worker])) {
                masterSol = wcsp->getSolution(); // the master sends the best solution or nothing
                bestsolWork[worker] = masterUb;
            }
            assert(masterUb >= MAX_COST || wcsp->getUb() == masterUb);
            Work work(*cp_, *open_, clb, wcsp->getUb(), masterSol);
            if (!ToulBar2::hbfs) {
                work.hbfs = false;
            }
            activeWork[worker] = work.open[0];
            idleQ.pop(); // pop it, hence the worker is considered active

            if (ToulBar2::verbose >= 1)
                cout << ">>> master send to worker " << worker << " a message " << work << endl;
            reqs.push_back(world.isend(worker, WORKTAG, work)); // non-blocking send: the master send work to an idle worker
        }
        double beginWaiting = realTime();
        mpi::wait_all(reqs.begin(), reqs.end());

        Work work2; // object work2 will be populated with workers' best solution ub and other information from this worker after it has performed a (can be partial in burst mode) DFS
        mpi::status status2 = world.recv(mpi::any_source, mpi::any_tag, work2); // blocking recv to wait for matching messages from any worker
        hbfsWaitingTime += realTime() - beginWaiting;

        wcsp->updateUb(work2.ub);
        open_->updateUb(work2.ub);
        nbNodes += work2.nbNodes;
        nbBacktracks += work2.nbBacktracks;
        ((WCSP*)wcsp)->incNbDEE(work2.nbDEE);
        nbRecomputationNodes += work2.nbRecomputationNodes;

        if (!work2.open.empty()) { // the master updates its CPStore with the decisions associated with the nodes sent by the worker
            for (ptrdiff_t i = 0; i < (ptrdiff_t)work2.cp.size(); i++) {
                cp_->push_back(work2.cp[i]);
            }
            for (size_t i = 0; i < work2.open.size(); i++) { // push the nodes sent by the worker
                OpenNode node = work2.open[i];
                node.first += cp_->start;
                node.last += cp_->start;
                open_->push(node);
            }
            cp_->stop = cp_->start + work2.cp.size();
            cp_->store();
        }

        if (cluster) {
            assert(work2.lb <= work2.ub);
            open_->updateClosedNodesLb(work2.lb);
        }
        cub = MIN(cub, work2.ub);

        if (status2.tag() == IDLETAG) {
            activeWork.erase(status2.source());
            idleQ.push(status2.source());
        }

        Cost minLbWorkers = MAX_COST;
        for (std::unordered_map<int, OpenNode>::const_iterator it = activeWork.begin(); it != activeWork.end(); ++it) { // compute the min of lb among those of active workers
            if (it->second.getCost() < minLbWorkers)
                minLbWorkers = it->second.getCost();
        }

        if (ToulBar2::eps && open_->size() >= static_cast<std::size_t>(ToulBar2::eps)) {
            for (std::unordered_map<int, OpenNode>::const_iterator it = activeWork.begin(); it != activeWork.end(); ++it) { // push back unfinished work
                open_->push(it->second);
            }
            epsDumpSubProblems(*cp_, *open_);
            vector<mpi::request> reqs;
            for (int i = 0; i < world.size(); i++)
                if (i != MASTER) {
                    reqs.push_back(world.isend(i, DIETAG, Work()));
                }
            mpi::wait_all(reqs.begin(), reqs.end());
            ToulBar2::interrupted = true;
            throw TimeOut();
        }

        if (cp_->size() >= static_cast<std::size_t>(ToulBar2::hbfsCPLimit) || open_->size() >= static_cast<std::size_t>(ToulBar2::hbfsOpenNodeLimit)) {
            ToulBar2::hbfs = 0;
            ToulBar2::hbfsGlobalLimit = 0;
            if (cluster) {
                cluster->hbfsGlobalLimit = LONGLONG_MAX;
                cluster->hbfsLimit = LONGLONG_MAX;
            } else
                hbfsLimit = LONGLONG_MAX;
        }

        clb = MAX(clb, MIN(minLbWorkers, open_->getLb()));

        showGap(clb, cub);

        if (work2.ub < initub) { // if the master receives an improving solution from the worker
            assert(work2.sol.size() == wcsp->numberOfVariables());

            // transformation of vector to map necessary because wcsp->getSolution returns a vector
            map<int, Value> workerSolMap;
            for (int i = 0; i < int(work2.sol.size()); i++) { // convert vector to map
                workerSolMap[i] = work2.sol[i];
            }
            wcsp->setSolution(work2.ub, &workerSolMap); // take current assignment and stock it in solution

            assert(work2.ub == wcsp->getUb());
            bestsolWork[status2.source()] = work2.ub; // the worker knows this solution already

            if (ToulBar2::verbose >= 0 || ToulBar2::showSolutions) {
                if (ToulBar2::haplotype)
                    cout << "***New solution: " << work2.ub << " log10like: " << ToulBar2::haplotype->Cost2LogProb(work2.ub) / Log(10.) << " logProb: " << ToulBar2::haplotype->Cost2LogProb(work2.ub) << " (" << nbBacktracks << " backtracks, " << nbNodes << " nodes, depth " << Store::getDepth() << ", " << ((ToulBar2::parallel) ? (realTime() - ToulBar2::startRealTime) : (cpuTime() - ToulBar2::startCpuTime)) << " seconds)" << endl;
                else if (!ToulBar2::bayesian)
                    cout << "New solution: " << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->getDPrimalBound() << std::setprecision(DECIMAL_POINT) << " (" << nbBacktracks << " backtracks, " << nbNodes << " nodes, depth " << Store::getDepth() << ", " << ((ToulBar2::parallel) ? (realTime() - ToulBar2::startRealTime) : (cpuTime() - ToulBar2::startCpuTime)) << " seconds)" << endl;
                else
                    cout << "New solution: " << std::setprecision(ToulBar2::decimalPoint) << wcsp->getDPrimalBound() << std::setprecision(DECIMAL_POINT) << " energy: " << -(wcsp->Cost2LogProb(work2.ub) + ToulBar2::markov_log) << " prob: " << std::scientific << wcsp->Cost2Prob(work2.ub) * Exp(ToulBar2::markov_log) << std::fixed << " (" << nbBacktracks << " backtracks, " << nbNodes << " nodes, depth " << Store::getDepth() << ", " << ((ToulBar2::parallel) ? (realTime() - ToulBar2::startRealTime) : (cpuTime() - ToulBar2::startCpuTime)) << " seconds)" << endl;
            }

            if (ToulBar2::showSolutions) {
                wcsp->printSolution(cout);
                cout << endl;
            }
            if (!ToulBar2::uaieval && ToulBar2::writeSolution && ToulBar2::solutionFile != NULL) {
                rewind(ToulBar2::solutionFile);
                wcsp->printSolution(ToulBar2::solutionFile);
                fprintf(ToulBar2::solutionFile, "\n");
            }
            if (ToulBar2::xmlflag) {
                cout << "o " << std::fixed << std::setprecision(0) << wcsp->getDPrimalBound() << std::setprecision(DECIMAL_POINT) << endl; //" ";
                ((WCSP*)wcsp)->solution_XML(false);
            }
            if (ToulBar2::maxsateval) {
                cout << "o " << work2.ub << endl;
            }
            if (ToulBar2::uaieval && !ToulBar2::isZ) {
                ((WCSP*)wcsp)->solution_UAI(work2.ub);
            }
        }
    }
    assert(clb >= initiallb && cub <= initialub);
    assert(clb <= cub);

    if (clb == cub) {
        vector<mpi::request> reqs;
        for (int i = 0; i < world.size(); i++)
            if (i != MASTER) {
                reqs.push_back(world.isend(i, DIETAG, Work()));
            }
        mpi::wait_all(reqs.begin(), reqs.end());
    }
    return make_pair(clb, cub);
}

pair<Cost, Cost> Solver::hybridSolveWorker(Cluster* cluster, Cost clb, Cost cub)
{
    if (ToulBar2::verbose >= 1) {
        if (cluster)
            cout << "hybridSolveWorker#" << world.rank() << " C" << cluster->getId() << " " << clb << " " << cub << endl;
        else
            cout << "hybridSolveWorker#" << world.rank() << " " << clb << " " << cub << endl;
    }
    assert(clb < cub);
    assert(wcsp->getUb() == cub);
    assert(wcsp->getLb() <= clb);
    assert(ToulBar2::hbfs);
    CPStore* cp_ = NULL;
    OpenList* open_ = NULL;

    if (cluster) {
        // BFS with BTD on current cluster (can be root or not)
        assert(cluster->cp);
        cp_ = cluster->cp;
        assert(cluster == wcsp->getTreeDec()->getRoot());
        if (!cluster->open)
            cluster->open = new OpenList();
        cluster->setUb(cub); // global problem upper bound
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

    nbHybrid++;

    while (true) {
        assert(open_->empty());
        cp_->stop = 0; // to have stop=start=index=0
        cp_->store();

        initWorkerNbNodes = nbNodes;
        initWorkerNbBacktracks = nbBacktracks;
        initWorkerNbDEE = wcsp->getNbDEE();
        initWorkerNbRecomputationNodes = nbRecomputationNodes;

        Work work;
        double beginWaiting = realTime();
        mpi::status status = world.recv(MASTER, mpi::any_tag, work); // blocking recv from the master
        hbfsWaitingTime += realTime() - beginWaiting;

        if (status.tag() == DIETAG) {
            ToulBar2::limited = true;
            throw TimeOut();
        }

        if (!work.sol.empty()) { // the master is supposed to always sends an optimal solution or nothing
            map<int, Value> masterSolMap; // convert vector to map
            for (int i = 0; i < int(work.sol.size()); i++) {
                masterSolMap[i] = work.sol[i];
            }
            assert(work.ub <= cub);
            wcsp->setSolution(work.ub, &masterSolMap); // take current assignment and stock it in solution
        }

        cub = MIN(cub, work.ub);

        wcsp->updateUb(work.ub); // update global UB in worker's wcsp object
        open_->updateUb(work.ub); // update cub and clb that are attributes of worker's open queue
        // TODO: enforceUb() and propagate() at each better solution received (see sequential HBFS)

        assert(work.open.size() == 1); // only one open node from the master

        for (size_t i = 0; i < work.cp.size(); i++) {
            addChoicePoint(work.cp[i].op, work.cp[i].varIndex, work.cp[i].value, work.cp[i].reverse); // update work cp->index
        }

        addOpenNode(*cp_, *open_, work.open[0].getCost()); // update of cp->stop and push node with first= cp-> start and last= cp->index
        assert(open_->size() == 1);
        assert(open_->top().first == 0);
        assert((size_t)open_->top().last == work.cp.size());

        cp_->store();

        if (!work.hbfs) { // master is full, must switch to DFS only
            ToulBar2::hbfs = 0;
            ToulBar2::hbfsGlobalLimit = 0;
            if (cluster) {
                cluster->hbfsGlobalLimit = LONGLONG_MAX;
                cluster->hbfsLimit = LONGLONG_MAX;
            } else {
                hbfsLimit = LONGLONG_MAX;
            }
        }
        if (cluster) {
            cluster->hbfsGlobalLimit = ((ToulBar2::hbfsGlobalLimit > 0) ? (nbBacktracks + ToulBar2::hbfsGlobalLimit) : LONGLONG_MAX);
            cluster->hbfsLimit = ((ToulBar2::hbfs > 0) ? (cluster->nbBacktracks + ToulBar2::hbfs) : LONGLONG_MAX);
            assert(wcsp->getTreeDec()->getCurrentCluster() == cluster);
            wcsp->setUb(cub);
            assert(cluster->isActive());
            assert(cluster->getLbRec() == wcsp->getLb());
        } else {
            hbfsLimit = ((ToulBar2::hbfs > 0) ? (nbBacktracks + ToulBar2::hbfs) : LONGLONG_MAX);
        }

        int storedepthBFS = Store::getDepth();
        int storeVAC = ToulBar2::vac;
        try {
            Store::store(); // store the depth of the DFS search
            OpenNode nd = open_->top(); // get a reference on the best node (min lower bound or, in case of equality, max depth)
            open_->pop(); // best node is taken from priority queue open
            if (ToulBar2::verbose >= 3) {
                if (wcsp->getTreeDec())
                    cout << "[C" << wcsp->getTreeDec()->getCurrentCluster()->getId() << "] ";
                cout << "[ " << nd.getCost() << ", " << cub << "] ( " << open_->size() << "+1 still open)" << endl;
            }
            if (ToulBar2::vac < 0 && Store::getDepth() + (nd.last - nd.first) >= abs(ToulBar2::vac))
                ToulBar2::vac = 0;
            restore(*cp_, nd); // replay the sequence of decisions and recompute soft arc consistency
            Cost bestlb = MAX(nd.getCost(), wcsp->getLb());
            bestlb = MAX(bestlb, work.lb);
            if (cluster) {
                pair<Cost, Cost> res = recursiveSolve(cluster, bestlb, cub);
                assert(res.first <= res.second);
                assert(res.first >= bestlb);
                assert(res.second <= cub);
                assert(res.second == cub || cluster->getUb() == res.second);
                assert(open_->empty() || open_->top().getCost() >= nd.getCost());
                work.lb = MAX(work.lb, res.first);
                cub = MIN(cub, res.second);
            } else {
                if (ToulBar2::vac < 0)
                    ToulBar2::vac = 0;
                recursiveSolve(bestlb); // call DFS (can generate a Contradiction, even after finding a better solution)
                work.lb = MAX(work.lb, bestlb);
                cub = MIN(cub, wcsp->getUb());
            }
        } catch (const Contradiction&) {
            wcsp->whenContradiction();
            work.lb = MAX(work.lb, wcsp->getUb());
        }
        if (!cluster) { // synchronize current upper bound with DFS (without tree decomposition)
            cub = MIN(cub, wcsp->getUb());
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
            } else {
                hbfsLimit = LONGLONG_MAX;
            }
        }

        if (ToulBar2::hbfs && nbRecomputationNodes > 0) { // wait until a nonempty open node is restored (at least after first global solution is found)
            if (nbRecomputationNodes > nbNodes / ToulBar2::hbfsBeta && ToulBar2::hbfs * 2 <= ToulBar2::hbfsGlobalLimit)
                ToulBar2::hbfs *= 2;
            else if (nbRecomputationNodes < nbNodes / ToulBar2::hbfsAlpha && ToulBar2::hbfs >= 2)
                ToulBar2::hbfs /= 2;
            if (ToulBar2::debug >= 2)
                cout << "HBFS backtrack limit for worker#" << world.rank() << ": " << ToulBar2::hbfs << endl;
        }

        vector<Value> workerSol;
        if (!ToulBar2::burst) {
            Cost newWorkerUb = wcsp->getSolutionCost();
            if (newWorkerUb < work.ub) {
                workerSol = wcsp->getSolution();
                assert(cub == newWorkerUb);
            }
        } else {
            open_->init(); // clb=cub=MAX_COST  method added to init openList attributes
            cp_->clear(); // size = 0  added to put new cp out of the while(1)
            assert(open_->empty());
        }
        Work work2(*cp_, *open_, nbNodes - initWorkerNbNodes, nbBacktracks - initWorkerNbBacktracks, wcsp->getNbDEE() - initWorkerNbDEE, nbRecomputationNodes - initWorkerNbRecomputationNodes, work.lb, cub, workerSol);
        if (ToulBar2::verbose >= 1)
            cout << ">>> worker " << world.rank() << " send closing-node message to master " << work2 << endl;
        beginWaiting = realTime();
        mpi::request req = world.isend(MASTER, IDLETAG, work2); // non-blocking send to master saying we have finished exploring its open node
        while (!req.test().is_initialized() && !MPI_interrupted())
            ;
        hbfsWaitingTime += realTime() - beginWaiting;
    }

    assert(clb <= cub);
    assert(false); // a worker should never stop, except by throwing an exception (TimeOut) when receiving DIETAG
    return make_pair(clb, cub); // used only to avoid warnings: a worker should not return anything!
}
#endif

void Solver::beginSolve(Cost ub)
{
    // Last-minute compatibility checks for ToulBar2 selected options
    assert(ub > MIN_COST);
    if (ToulBar2::allSolutions && ToulBar2::btdMode == 1 && ub > 1) {
        cerr << "Error: Solution enumeration by BTD-like search methods is only possible for feasability (use -ub=1 and integer costs only)." << endl;
        throw BadConfiguration();
    }
    if (ToulBar2::allSolutions && ToulBar2::btdMode == 1 && ub == 1 && ToulBar2::hbfs) {
        cerr << "Error: Hybrid best-first search cannot currently look for all solutions when BTD mode is activated. Shift to DFS (use -hbfs:)." << endl;
        throw BadConfiguration();
    }
    if ((ToulBar2::hve <= 0 || ToulBar2::pwc < 0) && ToulBar2::FullEAC && ToulBar2::vac > 1 && (wcsp->numberOfConnectedConstraints() > ((ToulBar2::VAClin)?wcsp->numberOfConnectedKnapsackConstraints():0) + wcsp->numberOfConnectedBinaryConstraints() || ToulBar2::elimDegree_preprocessing >= 3 || ToulBar2::preprocessTernaryRPC != 0)) {
        cerr << "Error: VAC during search and Full EAC variable ordering heuristic not implemented with non binary cost functions in extension (remove -vacint option)." << endl;
        throw BadConfiguration();
    }
    if (ToulBar2::searchMethod != DFBB) {
        if (ToulBar2::vnsLDSmax < 0)
            ToulBar2::vnsLDSmax = wcsp->getDomainSizeSum() - wcsp->numberOfUnassignedVariables();
        if (!ToulBar2::lds) {
            ToulBar2::vnsLDSmin = wcsp->getDomainSizeSum() - wcsp->numberOfUnassignedVariables();
            ToulBar2::vnsLDSmax = wcsp->getDomainSizeSum() - wcsp->numberOfUnassignedVariables();
        }
        if (ToulBar2::vnsKmax <= 0)
            ToulBar2::vnsKmax = wcsp->numberOfUnassignedVariables();
    }
    if (wcsp->isGlobal() && ToulBar2::btdMode >= 1) {
        cout << "Error: cannot use BTD-like search methods with monolithic global cost functions (remove -B or -F options)." << endl;
        throw BadConfiguration();
    }
#ifdef ILOGCPLEX
    if (wcsp->isPLPS() && ToulBar2::hbfs) {
        cout << "Error: cannot use Hybrid best-first search with PLPS slinear global cost functions (add -hbfs: option)." << endl;
        throw BadConfiguration();
    }
#endif
    if (ToulBar2::incop_cmd.size() > 0) {
        for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
            if (wcsp->unassigned(i) && !wcsp->enumerated(i)) {
                cout << "Warning! Cannot use INCOP local search with bounds arc propagation (non enumerated variable domains)." << endl;
                ToulBar2::incop_cmd = "";
                break;
            }
        }
    }
    if (ToulBar2::pils_cmd.size() > 0) {
        for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
            if (wcsp->unassigned(i) && !wcsp->enumerated(i)) {
                cout << "Warning! Cannot use PILS local search with bounds arc propagation (non enumerated variable domains)." << endl;
                ToulBar2::pils_cmd = "";
                break;
            }
        }
    }
    if (ToulBar2::lrBCD_cmd.size() > 0) {
        for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
            if (wcsp->unassigned(i) && !wcsp->enumerated(i)) {
                cout << "Warning! Cannot use LR-BCD local search with bounds arc propagation (non enumerated variable domains)." << endl;
                ToulBar2::pils_cmd = "";
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
#ifdef OPENMPI
    hbfsWaitingTime = 0.;
    initWorkerNbNodes = 0;
    initWorkerNbBacktracks = 0;
    initWorkerNbDEE = 0;
    initWorkerNbRecomputationNodes = 0;
    if (ToulBar2::parallel && ToulBar2::hbfs) {
        activeWork.clear();
        bestsolWork.clear();
        assert(idleQ.empty());
        for (int i = 0; i < world.size(); i++)
            if (i != MASTER) {
                idleQ.push(i);
            }
    }
#endif

    if (ToulBar2::DEE)
        ToulBar2::DEE_ = ToulBar2::DEE; // enforces PSNS after closing the model

    //    if (ToulBar2::addAMOConstraints != -1)
    //        ToulBar2::addAMOConstraints_ = true; // only bound propagation for knapsack constraints before adding them AMO constraints

    if (!ToulBar2::isZ && CSP(wcsp->getLb(), wcsp->getUb()) && ToulBar2::setvalue != tb2setvalue) { // do not modify (weakening) local consistency if there are global weighted CSP constraints
        ToulBar2::LcLevel = LC_AC;
        ToulBar2::vac = 0;
        ToulBar2::useRASPS = 0;
    }

    // reactivate on-the-fly variable elimination and dead-end elimination if needed
    for (int i = wcsp->numberOfVariables() - 1; i >= 0; i--) {
        if (wcsp->unassigned(i)) {
            ((WCSP*)wcsp)->getVar(i)->queueEliminate();
            ((WCSP*)wcsp)->getVar(i)->queueDEE();
        }
    }
}

Cost Solver::preprocessing(Cost initialUpperBound)
{
    Long hbfs_ = ToulBar2::hbfs;
    ToulBar2::hbfs = 0; // do not perform hbfs operations in preprocessing except for building tree decomposition
    if (!ToulBar2::isZ) {
        wcsp->enforceUb();
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
    initGap(wcsp->getLb(), wcsp->getUb());
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
                initGap(wcsp->getLb(), wcsp->getUb());
            }
        }

        vector<int> elimorder(wcsp->numberOfVariables(), -1);
        vector<int> revelimorder(wcsp->numberOfVariables(), -1);
        for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
            revelimorder[i] = i;
            elimorder[wcsp->numberOfVariables() - i - 1] = i;
        }
        Cost previouslb = wcsp->getLb();
        do {
            previouslb = wcsp->getLb();
            wcsp->setDACOrder(revelimorder);
            wcsp->setDACOrder(elimorder);
            initGap(wcsp->getLb(), wcsp->getUb());
            if (ToulBar2::verbose >= 0 && wcsp->getLb() > previouslb) {
                if (ToulBar2::uai)
                    cout << "Reverse original DAC dual bound: " << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->getDDualBound() << std::setprecision(DECIMAL_POINT) << " energy: " << -(wcsp->Cost2LogProb(wcsp->getLb()) + ToulBar2::markov_log) << " (+" << 100. * (wcsp->getLb() - previouslb) / wcsp->getLb() << "%)" << endl;
                else
                    cout << "Reverse original DAC dual bound: " << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->getDDualBound() << std::setprecision(DECIMAL_POINT) << " (+" << 100. * (wcsp->getLb() - previouslb) / wcsp->getLb() << "%)" << endl;
            }
        } while (wcsp->getLb() > previouslb && (Double)100. * (wcsp->getLb() - previouslb) / wcsp->getLb() > (Double)0.5);
    }
    wcsp->preprocessing(); // preprocessing after initial propagation
    initGap(wcsp->getLb(), wcsp->getUb());
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
            initGap(wcsp->getLb(), wcsp->getUb());
        }
    }

    // special data structure to be initialized for variable ordering heuristics including weighted degrees and tightness
    initVarHeuristic();

    int lds = ToulBar2::lds;
    ToulBar2::lds = 0; // avoid TimeOut exception when new solutions found
    if (ToulBar2::incop_cmd.size() > 0 && getWCSP()->numberOfUnassignedVariables() > 0) {
        double incopStartTime = cpuTime();
        vector<Value> bestsol(getWCSP()->numberOfVariables(), 0);
        for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++)
            bestsol[i] = (wcsp->canbe(i, wcsp->getBestValue(i)) ? wcsp->getBestValue(i) : wcsp->getSupport(i));
        narycsp(ToulBar2::incop_cmd, bestsol);
        if (ToulBar2::verbose >= 0)
            cout << "INCOP solving time: " << cpuTime() - incopStartTime << " seconds." << endl;
    }
    if (ToulBar2::pils_cmd.size() > 0 && getWCSP()->numberOfUnassignedVariables() > 2 && getWCSP()->numberOfConnectedBinaryConstraints() > 1) {
        double pilsStartTime = cpuTime();
        vector<Value> bestsol(getWCSP()->numberOfVariables(), 0);
        for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++)
            bestsol[i] = (wcsp->canbe(i, wcsp->getBestValue(i)) ? wcsp->getBestValue(i) : wcsp->getSupport(i));
        pils(ToulBar2::pils_cmd, bestsol);
        if (ToulBar2::verbose >= 0)
            cout << "PILS solving time: " << cpuTime() - pilsStartTime << " seconds." << endl;
    }
    if (ToulBar2::lrBCD_cmd.size() > 0 && getWCSP()->numberOfUnassignedVariables() > 2 && getWCSP()->numberOfConnectedBinaryConstraints() > 1) {
        double lrBCDStartTime = cpuTime();
        vector<Value> bestsol(getWCSP()->numberOfVariables(), 0);
        for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++)
            bestsol[i] = (wcsp->canbe(i, wcsp->getBestValue(i)) ? wcsp->getBestValue(i) : wcsp->getSupport(i));
        lrBCD(ToulBar2::lrBCD_cmd, bestsol);
        if (ToulBar2::verbose >= 0)
            cout << "LR-BCD solving time: " << cpuTime() - lrBCDStartTime << " seconds." << endl;
    }
    ToulBar2::lds = lds;

#ifdef BOOST
    if (ToulBar2::addAMOConstraints != -1) {
        ToulBar2::addAMOConstraints_ = true;
        wcsp->addAMOConstraints();
        ToulBar2::addAMOConstraints_ = false;
    }
#endif

    if (ToulBar2::singletonConsistency) {
        singletonConsistency(ToulBar2::singletonConsistency);
        wcsp->propagate();
        wcsp->resetTightnessAndWeightedDegree();
    }

    ToulBar2::hbfs = hbfs_; // do not perform hbfs operations in preprocessing except for building tree decomposition

    if (ToulBar2::verbose >= 0)
        cout << "Preprocessing time: " << ((ToulBar2::parallel) ? (realTime() - ToulBar2::startRealTime) : (cpuTime() - ToulBar2::startCpuTime)) << " seconds." << endl;
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
        int nbdelayedblp = 0;
        if (ToulBar2::bilevel) {
            for (auto s : ((WCSP*)wcsp)->delayedCtrBLP) {
                nbdelayedblp += s.size();
            }
        }
        if (wcsp->numberOfUnassignedVariables() == 0 || (wcsp->numberOfConnectedConstraints() == 0 && nbdelayedblp == 0)) {
            ToulBar2::approximateCountingBTD = false;
            ToulBar2::btdMode = 0;
        } else {
            // ToulBar2::vac = 0; // VAC is not compatible with restricted tree decomposition propagation
            wcsp->buildTreeDecomposition();
            if (ToulBar2::bilevel) {
                Cluster* problem0 = wcsp->getTreeDec()->getRoot();
                auto iter = problem0->beginEdges();
#ifndef NDEBUG
                Cluster* problem1 = *iter;
#endif
                ++iter;
                Cluster* problem2 = *iter;
                ++iter;
#ifndef NDEBUG
                Cluster* negproblem2 = *iter;
#endif
                // problem2.isused = false //FIXME???
                problem2->deactivate(); // avoid future propagation (NC*) in left child Problem2 when branching on Problem0
                assert(problem1->getLb() == MIN_COST);
                assert(problem1->getCurrentDeltaUb() == MIN_COST);
                assert(problem2->getLb() == MIN_COST);
                assert(problem2->getCurrentDeltaUb() == MIN_COST);
                assert(negproblem2->getLb() == MIN_COST);
                assert(negproblem2->getCurrentDeltaUb() == MIN_COST);
            }
        }
    } else if (ToulBar2::weightedDegree && (((Long)wcsp->numberOfConnectedConstraints()) >= ((Long)abs(ToulBar2::weightedDegree)))) {
        if (ToulBar2::verbose >= 0)
            cout << "Weighted degree heuristic disabled (#costfunctions=" << wcsp->numberOfConnectedConstraints() << " >= " << abs(ToulBar2::weightedDegree) << ")" << endl;
        ToulBar2::weightedDegree = 0;
    }

    if (ToulBar2::dumpWCSP) {
        dump_wcsp(ToulBar2::problemsaved_filename.c_str(), ToulBar2::dumpOriginalAfterPreprocessing, static_cast<ProblemFormat>((ToulBar2::dumpWCSP >> 1) + (ToulBar2::dumpWCSP & 1)));
        cout << "Problem saved after preprocessing, then program stopped." << endl;
        throw TimeOut();
    }

    if (ToulBar2::knapsackDP == -1)
        ToulBar2::knapsackDP = -2;
    return initialUpperBound;
}

bool Solver::solve(bool first)
{
    beginSolve(wcsp->getUb());
    initGap(wcsp->getLb(), wcsp->getUb());
    int _DEE_ = ToulBar2::DEE;
    int _elimDegree_ = ToulBar2::elimDegree;

    Cost initialUpperBound = wcsp->getUb();

    //        Store::store();   // if uncommented, solve() does not change the problem but all preprocessing operations will allocate in backtrackable memory
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
                        if (ToulBar2::debug >= 1 && ToulBar2::weightedDegree) {
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
                                pair<Cost, Cost> res = hybridSolve();
                                globalLowerBound = res.first;
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
                                start->setLb(wcsp->getLb()); // initial lower bound found by propagation is associated to tree decomposition root cluster
                            if (ToulBar2::bilevel) {
                                // propagate channeling constraints between leader (Problem0) and negative follower (NegProblem2) problems only
                                for (int ctrIndex : ((WCSP*)wcsp)->delayedCtrBLP[2]) {
                                    Constraint* ctr = ((WCSP*)wcsp)->getCtr(ctrIndex);
                                    assert(ctr->deconnected());
                                    static vector<Cost> costs;
                                    Constraint* incCtr = NULL;
                                    if (ctr->isBinary()) {
                                        costs.resize(ctr->getDomainInitSizeProduct(), MIN_COST);
                                        incCtr = ((WCSP*)wcsp)->getCtr(wcsp->postIncrementalBinaryConstraint(ctr->getVar(0)->wcspIndex, ctr->getVar(1)->wcspIndex, costs));
                                        ((BinaryConstraint*)incCtr)->addCosts((BinaryConstraint*)ctr);
                                    } else if (ctr->isTernary()) {
                                        costs.resize(ctr->getDomainInitSizeProduct(), MIN_COST);
                                        incCtr = ((WCSP*)wcsp)->getCtr(wcsp->postIncrementalTernaryConstraint(ctr->getVar(0)->wcspIndex, ctr->getVar(1)->wcspIndex, ctr->getVar(2)->wcspIndex, costs));
                                        ((TernaryConstraint*)incCtr)->addCosts((TernaryConstraint*)ctr);
                                    } else {
                                        cerr << "Sorry, bilevel optimization not implemented for this type of channeling cost function:" << *ctr << endl;
                                        throw WrongFileFormat();
                                    }
                                    // incCtr->sumScopeIncluded(ctr);
                                    // incCtr->assignCluster(); //SdG: done before inside postIncrementalXXXConstraint
                                    incCtr->propagate();
                                }
                            }
                            switch (ToulBar2::btdMode) {
                            case 0:
                            case 1: {
                                if (ToulBar2::allSolutions) {
                                    timeDeconnect = 0.;
                                    BigInteger cartesianProduct = 1;
                                    nbSol = (wcsp->numberOfConnectedConstraints() == 0) ? (wcsp->cartProd(cartesianProduct), cartesianProduct) : sharpBTD(start);
                                    if (ToulBar2::approximateCountingBTD && nbSol > 0. && td->getRoot()->getNbVars() == 0) { // if there are several parts
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
                                            if (ToulBar2::bilevel && ToulBar2::verbose >= 0) {
                                                Cluster* problem0 = td->getRoot();
                                                auto iter = problem0->beginEdges();
                                                Cluster* problem1 = *iter;
                                                ++iter;
#ifndef NDEBUG
                                                Cluster* problem2 = *iter;
#endif
                                                ++iter;
                                                Cluster* negproblem2 = *iter;
                                                assert(problem2->getLb() == MIN_COST);
                                                assert(problem2->getCurrentDeltaUb() == MIN_COST);
                                                // cout << "C0.lb: " << problem0->getLb() << " C1.lb: " << problem1->getLb() << " C2.lb: " << problem2->getLb() << " NegC2.lb: " << negproblem2->getLb() << " NegC2.delta >= " << negproblem2->getCurrentDeltaLb() << " NegC2.delta <= " << negproblem2->getCurrentDeltaUb() << endl;
                                                Cost lbP1 = problem0->getLb() + problem1->getLb() - ToulBar2::initialLbBLP[2] - negproblem2->getCurrentDeltaUb();
                                                Cost lbP2 = ToulBar2::initialLbBLP[1];
                                                Cost lbNegP2 = negproblem2->getLb() + ToulBar2::initialLbBLP[2] + negproblem2->getCurrentDeltaLb();
                                                cout << "Initial lower bound for the restricted leader problem (without subtracting the follower objective): " << wcsp->Cost2RDCost(lbP1 - ToulBar2::negCostBLP[0]) << endl;
                                                cout << "Initial lower bound for the follower problem: " << wcsp->Cost2RDCost(lbP2 - ToulBar2::negCostBLP[1]) << endl;
                                                cout << "Initial strict upper bound for the follower problem: " << wcsp->Cost2RDCost(-(lbNegP2 - ToulBar2::negCostBLP[2]) + UNIT_COST) << endl;
                                                cout << "Initial lower bound for the leader problem: " << wcsp->Cost2RDCost(wcsp->getLb() - wcsp->getNegativeLb()) << endl;
                                            }
                                            res = hybridSolve(start, MAX(wcsp->getLb(), res.first), res.second);
                                            //				                if (res.first < res.second) cout << "Optimality gap: [ " <<  res.first << " , " << res.second << " ] " << (100. * (res.second-res.first)) / res.second << " % (" << nbBacktracks << " backtracks, " << nbNodes << " nodes)" << endl;
                                            globalLowerBound = res.first;
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
                                do { // TODO: set up for optimality gap pretty print
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
                                throw BadConfiguration();
                            }
                            }
                            if (ToulBar2::verbose >= 0 && nbHybrid >= 1)
                                cout << "HBFS open list restarts: " << (100. * (nbHybrid - nbHybridNew - nbHybridContinue) / nbHybrid) << " % and reuse: " << (100. * nbHybridContinue / nbHybrid) << " % of " << nbHybrid << endl;
                        } else {
                            if (ToulBar2::useRASPS) {
                                int idthres = ToulBar2::RASPSitThresholds.size() - 1;
                                Cost lastUb = wcsp->getUb();
                                do {
                                    enforceUb();
                                    wcsp->propagate();
                                    ToulBar2::RASPS = true;
                                    ((WCSP*)wcsp)->vac->iniThreshold(ToulBar2::RASPSlastitThreshold);
                                    ((WCSP*)wcsp)->vac->propagate(); // VAC done again
                                    ToulBar2::RASPS = false;
                                    if (ToulBar2::verbose >= 0)
                                        cout << "RASPS done in preprocessing at threshold " << ToulBar2::RASPSlastitThreshold << " (backtrack: " << nbBacktracks << " nodes: " << nbNodes << ")" << endl;
                                    while (idthres >= 0) {
                                        if (ToulBar2::RASPSitThresholds[idthres].first > ToulBar2::RASPSlastitThreshold) {
                                            ToulBar2::RASPSlastitThreshold = ToulBar2::RASPSitThresholds[idthres].first;
                                            break;
                                        } else {
                                            idthres--;
                                        }
                                    }
                                } while (nbBacktracks < ToulBar2::RASPSnbBacktracks && idthres >= 0 && wcsp->getUb() == lastUb);
                                enforceUb();
                                wcsp->propagate();
                                if (ToulBar2::RASPSreset) {
                                    wcsp->resetWeightedDegree();
                                    for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
                                        wcsp->setBestValue(i, wcsp->getSup(i) + 1);
                                        heuristics[i] = wcsp->getDegree(i);
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

                                        // get solution from previous solve and add pairwise Hamming distance constraint
                                        if (energies.size() > 0 && !extrapolatedBound) {
                                            vector<Value> solutionvec = wcsp->getSolution();
                                            map<int, Value> solutionmap;
                                            for (Variable* x : wcsp->getDivVariables()) {
                                                solutionmap[x->wcspIndex] = solutionvec[x->wcspIndex];
                                            }
                                            int arity = wcsp->getDivVariables().size();
                                            vector<int> scope;
                                            string parameters;
                                            switch (ToulBar2::divMethod) {
                                            case 0:
                                                ((WCSP*)wcsp)->addDivConstraint(wcsp->getDivVariables(), ToulBar2::divBound, solutionmap, ((WCSP*)wcsp)->divVarsId[energies.size() - 1], true);
                                                break;
                                            case 1:
                                                ((WCSP*)wcsp)->addHDivConstraint(wcsp->getDivVariables(), ToulBar2::divBound, solutionmap, ((WCSP*)wcsp)->divVarsId[energies.size() - 1], ((WCSP*)wcsp)->divHVarsId[energies.size() - 1], true);
                                                break;
                                            case 2:
                                                ((WCSP*)wcsp)->addTDivConstraint(wcsp->getDivVariables(), ToulBar2::divBound, solutionmap, ((WCSP*)wcsp)->divHVarsId[energies.size() - 1], true);
                                                break;
                                            case 3:
                                                parameters.append(to_string(-(arity - (int)ToulBar2::divBound)));
                                                for (Variable* x : wcsp->getDivVariables()) {
                                                    scope.push_back(x->wcspIndex);
                                                    parameters.append(" 1 ");
                                                    parameters.append(to_string(solutionmap[x->wcspIndex]));
                                                    parameters.append(" -1");
                                                }
                                                ((WCSP*)wcsp)->postKnapsackConstraint(scope, parameters, false, true);
                                                break;
                                            default:
                                                cerr << "Error: no such diversity encoding method: " << ToulBar2::divMethod << endl;
                                                throw BadConfiguration();
                                            }
                                            // wcsp->propagate(); it will propagate with a possibly modified DAC order at previous iterations resulting in undefined behavior
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
                                                    // ofstream os(to_string(this) + to_string("-wregular.dot"));
                                                    // printLayers(os, mdd);
                                                    // os.close();
                                                    switch (ToulBar2::divMethod) {
                                                    case 0:
                                                        ((WCSP*)wcsp)->addMDDConstraint(mdd, ToulBar2::divNbSol - 1); // ToulBar2::divNbSol = index of the relaxed constraint
                                                        break;
                                                    case 1:
                                                        ((WCSP*)wcsp)->addHMDDConstraint(mdd, ToulBar2::divNbSol - 1);
                                                        break;
                                                    case 2:
                                                        ((WCSP*)wcsp)->addTMDDConstraint(mdd, ToulBar2::divNbSol - 1);
                                                        break;
                                                    default:
                                                        cerr << "Error: no such MDD diversity encoding method: " << ToulBar2::divMethod << endl;
                                                        throw BadConfiguration();
                                                    }
                                                }
                                                // reactivate on-the-fly variable elimination and dead-end elimination (undone in preprocessing on diversity variables and later deactivated by endSolve call)
                                                ToulBar2::DEE_ = _DEE_;
                                                ToulBar2::elimDegree_ = _elimDegree_;
                                                for (int i = wcsp->numberOfVariables() - 1; i >= 0; i--) {
                                                    if (wcsp->unassigned(i)) {
                                                        ((WCSP*)wcsp)->getVar(i)->queueEliminate();
                                                        ((WCSP*)wcsp)->getVar(i)->queueDEE();
                                                    }
                                                }
                                                int vac = ToulBar2::vac;
                                                if (ToulBar2::vac) {
                                                    ToulBar2::vac = Store::getDepth() + 1; // enforces VAC at each new diverse solution found if VAC required
                                                }
                                                vector<int> revdac = wcsp->getBergeDecElimOrder();
                                                wcsp->enforceUb();
                                                wcsp->setDACOrder(revdac);
                                                ToulBar2::vac = vac;
                                                initGap(wcsp->getLb(), wcsp->getUb());
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
                                pair<Cost, Cost> res = hybridSolve();
                                globalLowerBound = res.first;
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
        wcsp->whenContradiction();
        ToulBar2::interrupted = false;
    }
    Store::restore(initdepth);

    if (ToulBar2::divNbSol <= 1)
        endSolve(wcsp->getSolutionCost() < initialUpperBound, wcsp->getSolutionCost(), !ToulBar2::limited);

    return (ToulBar2::isZ || ToulBar2::allSolutions || wcsp->getSolutionCost() < initialUpperBound);
}

void Solver::endSolve(bool isSolution, Cost cost, bool isComplete)
{
    ToulBar2::DEE_ = 0;
    ToulBar2::elimDegree_ = -1;

    static string solType[4] = { "Optimum: ", "Primal bound: ", "Guaranteed primal bound: ", "Primal bound: " };

    int isLimited = (!isComplete) | ((ToulBar2::deltaUb != MIN_COST) << 1);

    if (ToulBar2::debug && wcsp->getTreeDec()) {
        wcsp->getTreeDec()->getRoot()->printStatsRec();
    }

    if (ToulBar2::isZ) {
        if (ToulBar2::verbose >= 1)
            cout << "NegativeShiftingCost= " << wcsp->getNegativeLb() << endl;
        if (ToulBar2::uaieval) {
            rewind((ToulBar2::writeSolution) ? ToulBar2::solutionFile : ToulBar2::solution_uai_file);
            fprintf((ToulBar2::writeSolution) ? ToulBar2::solutionFile : ToulBar2::solution_uai_file, "PR\n");
            fprintf((ToulBar2::writeSolution) ? ToulBar2::solutionFile : ToulBar2::solution_uai_file, PrintFormatProb, (wcsp->LogSumExp(ToulBar2::logZ, ToulBar2::logU) + ToulBar2::markov_log) / Log(10.));
            fprintf((ToulBar2::writeSolution) ? ToulBar2::solutionFile : ToulBar2::solution_uai_file, "\n");
        }
        cout << std::setprecision(ToulBar2::resolution);
        cout << (ToulBar2::logZ + ToulBar2::markov_log - (TLogProb)wcsp->getMaxDomainSize() * Exp10(-(TLogProb)ToulBar2::resolution)) << " <= Log(Z) <= ";
        cout << (wcsp->LogSumExp(ToulBar2::logZ, ToulBar2::logU) + ToulBar2::markov_log + (TLogProb)wcsp->getMaxDomainSize() * Exp10(-(TLogProb)ToulBar2::resolution)) << " in " << nbBacktracks << " backtracks and " << nbNodes << " nodes and " << ((ToulBar2::parallel) ? (realTime() - ToulBar2::startRealTime) : (cpuTime() - ToulBar2::startCpuTime)) << " seconds" << endl;
        cout << (ToulBar2::logZ + ToulBar2::markov_log - (TLogProb)wcsp->getMaxDomainSize() * Exp10(-(TLogProb)ToulBar2::resolution)) / Log(10.) << " <= Log10(Z) <= ";
        cout << (wcsp->LogSumExp(ToulBar2::logZ, ToulBar2::logU) + ToulBar2::markov_log + (TLogProb)wcsp->getMaxDomainSize() * Exp10(-(TLogProb)ToulBar2::resolution)) / Log(10.) << " in " << nbBacktracks << " backtracks and " << nbNodes << " nodes and " << ((ToulBar2::parallel) ? (realTime() - ToulBar2::startRealTime) : (cpuTime() - ToulBar2::startCpuTime)) << " seconds" << endl;
        cout << std::setprecision(DECIMAL_POINT);
        return;
    }
    if (ToulBar2::allSolutions) {
        if (ToulBar2::verbose >= 0) {
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
            cout << "Time                   :    " << ((ToulBar2::parallel) ? (realTime() - ToulBar2::startRealTime) : (cpuTime() - ToulBar2::startCpuTime)) << " seconds" << endl;
            cout << "... in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << ((ToulBar2::DEE) ? (to_string(" ( ") + to_string(wcsp->getNbDEE()) + to_string(" removals by DEE)")) : to_string("")) << endl;
        }
        return;
    }

    if (ToulBar2::vac)
        wcsp->printVACStat();

#ifdef OPENMPI
    if (((ToulBar2::verbose >= 0 && !ToulBar2::parallel) || (ToulBar2::parallel && ToulBar2::verbose >= -1 && !ToulBar2::uai && !ToulBar2::xmlflag && !ToulBar2::maxsateval)) && nbHybrid >= 1 && nbNodes > 0) {
        cout << "Node redundancy during HBFS: " << 100. * nbRecomputationNodes / nbNodes;
        if (ToulBar2::parallel) {
            cout << " % (#pid: " << world.rank() << " wait: " << hbfsWaitingTime << " seconds)";
        }
        cout << endl;
    }
#else
    if (ToulBar2::verbose >= 0 && nbHybrid >= 1 && nbNodes > 0)
        cout << "Node redundancy during HBFS: " << 100. * nbRecomputationNodes / nbNodes << endl;
#endif

    if (ToulBar2::verbose >= 0 && ToulBar2::heuristicFreedom && wcsp->getTreeDec()) {
        cout << "Summary of adaptive BTD: " << endl;
        double sum = nbChoices + nbForcedChoices + nbForcedChoiceChange + nbReadOnly;
        if (sum != 0)
            cout << "% of new positive choices (i.e. cluster subtree is merged): " << (nbChoices / sum) * 100.0 << endl;
        else
            cout << "% of new positive choices: NA" << endl;

        if (nbChoices != 0)
            cout << "% of transitions (from positive/merged to negative/unmerged) due to an unproductive exploration (w.r.t positive choices): " << (nbChoiceChange / nbChoices) * 100.0 << endl;
        else
            cout << "% of transitions (from positive/merged to negative/unmerged) due to an unproductive exploration (w.r.t positive choices): NA" << endl;

        if (sum != 0)
            cout << "% of transitions due to nogood propagation: " << (nbForcedChoiceChange / sum) * 100.0 << endl;
        else
            cout << "% of transitions due to nogood propagation: NA" << endl;

        if (sum != 0)
            cout << "% of forced negative choices due to nogood propagation: " << (nbForcedChoices / sum) * 100.0 << endl;
        else
            cout << "% of forced negative choices due to nogood propagation: NA" << endl;

        if (sum != 0)
            cout << "% of choices without change: " << (nbReadOnly / sum) * 100.0 << endl;
        else
            cout << "% of choices without change: NA" << endl;

        cout << "Maximum cluster depth visited during search / maximum cluster depth of the original tree decomposition (except the root): " << solveDepth << " / " << wcsp->getTreeDec()->getMaxDepth() << endl;
        ;
    }

    if (isSolution) {
        if (ToulBar2::verbose >= 0 && !ToulBar2::uai && !ToulBar2::xmlflag && !ToulBar2::maxsateval) {
            if (ToulBar2::haplotype)
                cout << endl;

            if (isLimited == 2) {
                cout << "(" << ToulBar2::deltaUbS << "," << std::scientific << ToulBar2::deltaUbRelativeGap << std::fixed << ")-";
            }
            if (ToulBar2::haplotype) {
                cout << solType[isLimited] << cost << " log10like: " << ToulBar2::haplotype->Cost2LogProb(cost) / Log(10.) << " loglike: " << ToulBar2::haplotype->Cost2LogProb(cost) << " in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << ((ToulBar2::DEE) ? (to_string(" ( ") + to_string(wcsp->getNbDEE()) + to_string(" removals by DEE)")) : to_string("")) << " and " << ((ToulBar2::parallel) ? (realTime() - ToulBar2::startRealTime) : (cpuTime() - ToulBar2::startCpuTime)) << " seconds." << endl;
            } else if (!ToulBar2::bayesian) {
                if (!isComplete) {
                    cout << "Dual bound: " << std::fixed << std::setprecision(ToulBar2::decimalPoint) << getDDualBound() << std::setprecision(DECIMAL_POINT) << endl;
                }
                cout << solType[isLimited] << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->Cost2ADCost(cost) << std::setprecision(DECIMAL_POINT) << " in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << ((ToulBar2::DEE) ? (to_string(" ( ") + to_string(wcsp->getNbDEE()) + to_string(" removals by DEE)")) : to_string("")) << " and " << ((ToulBar2::parallel) ? (realTime() - ToulBar2::startRealTime) : (cpuTime() - ToulBar2::startCpuTime)) << " seconds." << endl;
            } else {
                cout << solType[isLimited] << cost << " energy: " << -(wcsp->Cost2LogProb(cost) + ToulBar2::markov_log) << std::scientific << " prob: " << wcsp->Cost2Prob(cost) * Exp(ToulBar2::markov_log) << std::fixed << " in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << ((ToulBar2::DEE) ? (to_string(" ( ") + to_string(wcsp->getNbDEE()) + to_string(" removals by DEE)")) : to_string("")) << " and " << ((ToulBar2::parallel) ? (realTime() - ToulBar2::startRealTime) : (cpuTime() - ToulBar2::startCpuTime)) << " seconds." << endl;
            }
        } else {
#ifdef OPENMPI
            if (ToulBar2::xmlflag && (!ToulBar2::parallel || world.rank() == MASTER)) {
#else
            if (ToulBar2::xmlflag) {
#endif
                ((WCSP*)wcsp)->solution_XML(!isLimited);
            } else if (ToulBar2::verbose >= 0 && ToulBar2::uai && !ToulBar2::isZ) {
                if (isLimited == 2)
                    cout << "(" << ToulBar2::deltaUbS << "," << std::scientific << ToulBar2::deltaUbRelativeGap << std::fixed << ")-";
                cout << solType[isLimited] << cost << " energy: " << -(wcsp->Cost2LogProb(cost) + ToulBar2::markov_log) << std::scientific << " prob: " << wcsp->Cost2Prob(cost) * Exp(ToulBar2::markov_log) << std::fixed << " in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << ((ToulBar2::DEE) ? (to_string(" ( ") + to_string(wcsp->getNbDEE()) + to_string(" removals by DEE)")) : to_string("")) << " and " << ((ToulBar2::parallel) ? (realTime() - ToulBar2::startRealTime) : (cpuTime() - ToulBar2::startCpuTime)) << " seconds." << endl;
#ifdef OPENMPI
            } else if (ToulBar2::maxsateval && !isLimited && (!ToulBar2::parallel || world.rank() == MASTER)) {
#else
            } else if (ToulBar2::maxsateval && !isLimited) {
#endif
                cout << "o " << cost << endl;
                cout << "s OPTIMUM FOUND" << endl;
                ((WCSP*)wcsp)->printSolutionMaxSAT(cout);
            }
        }
    } else {
        if (ToulBar2::verbose >= 0) {
            cout << "No solution" << ((!isLimited) ? "" : " found") << " in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << ((ToulBar2::DEE) ? (to_string(" ( ") + to_string(wcsp->getNbDEE()) + to_string(" removals by DEE)")) : to_string("")) << " and " << ((ToulBar2::parallel) ? (realTime() - ToulBar2::startRealTime) : (cpuTime() - ToulBar2::startCpuTime)) << " seconds." << endl;
        }
#ifdef OPENMPI
        if ((ToulBar2::maxsateval || ToulBar2::xmlflag) && !isLimited && (!ToulBar2::parallel || world.rank() == MASTER)) {
#else
        if ((ToulBar2::maxsateval || ToulBar2::xmlflag) && !isLimited) {
#endif
            //            cout << "o " << cost << endl;
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
        sumcost += 2. * std::abs(cost[e]);
    }
    Double multiplier = ((Double)MAX_COST) / sumcost;
    multiplier /= MEDIUM_COST;

    // create weighted binary clauses
    for (int e = 0; e < m; e++) {
        if (posx[e] != posy[e]) {
            vector<Cost> costs(4, 0);
            if (cost[e] > 0) {
                costs[1] = (Cost)roundl(multiplier * 2. * cost[e]);
                costs[2] = costs[1];
            } else {
                costs[0] = (Cost)roundl(multiplier * -2. * cost[e]);
                costs[3] = costs[0];
            }
            wcsp->postBinaryConstraint(posx[e] - 1, posy[e] - 1, costs);
        } else {
            if (cost[e] > 0) {
                unaryCosts1[posx[e] - 1] += (Cost)roundl(multiplier * cost[e]);
            } else {
                unaryCosts0[posx[e] - 1] += (Cost)roundl(multiplier * -cost[e]);
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
    ToulBar2::startRealTime = realTime();
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
            cout << "c " << after * sizeof(ChoicePoint) + td->getCurrentCluster()->open->capacity() * sizeof(OpenNode) << " Bytes allocated for hybrid best-first search open nodes at cluster " << td->getCurrentCluster()->getId() << "." << endl;
    } else {
        CPStore::size_type before = cp->capacity();
        cp->addChoicePoint(op, varIndex, value, reverse);
        CPStore::size_type after = cp->capacity();
        if (ToulBar2::verbose >= 0 && after > before && after > (1 << STORE_SIZE))
            cout << "c " << after * sizeof(ChoicePoint) + open->capacity() * sizeof(OpenNode) << " Bytes allocated for hybrid best-first search open nodes." << endl;
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
// void Solver::restore(CPStore &cp, OpenNode nd)
//{
//     if (ToulBar2::verbose >= 1) {
//         if (wcsp->getTreeDec()) cout << "[C" << wcsp->getTreeDec()->getCurrentCluster()->getId() << "] ";
//         cout << "restore open node " << nd.getCost(((wcsp->getTreeDec())?wcsp->getTreeDec()->getCurrentCluster()->getCurrentDeltaUb():MIN_COST)) << " (" << nd.first << ", " << nd.last << ")" << endl;
//     }
//     for (ptrdiff_t idx = nd.first; idx < nd.last; ++idx) {
//         assert(idx < cp.size());
//         if (ToulBar2::verbose >= 1) cout << "retrieve choice point " << CPOperation[cp[idx].op] << ((cp[idx].reverse)?"*":"") << " (" << wcsp->getName(cp[idx].varIndex) << ", " << cp[idx].value << ") at position " << idx  << endl;
//         if (ToulBar2::verbose >= 1) cout << *((WCSP *) wcsp)->getVar(cp[idx].varIndex) << endl;
//         assert(!wcsp->getTreeDec() || wcsp->getTreeDec()->getCurrentCluster()->isVar(cp[idx].varIndex));
//         switch (cp[idx].op) {
//         case CP_ASSIGN:
//             if (cp[idx].reverse && idx < nd.last-1) remove(cp[idx].varIndex, cp[idx].value);
//             else assign(cp[idx].varIndex, cp[idx].value);
//             break;
//         case CP_REMOVE:
//             if (cp[idx].reverse && idx < nd.last-1) assign(cp[idx].varIndex, cp[idx].value);
//             else remove(cp[idx].varIndex, cp[idx].value);
//             break;
//         case CP_INCREASE:
//             if (cp[idx].reverse && idx < nd.last-1) decrease(cp[idx].varIndex, cp[idx].value - 1);
//             else increase(cp[idx].varIndex, cp[idx].value);
//             break;
//         case CP_DECREASE:
//             if (cp[idx].reverse && idx < nd.last-1) increase(cp[idx].varIndex, cp[idx].value + 1);
//             else decrease(cp[idx].varIndex, cp[idx].value);
//             break;
//         default:
//             cerr << "unknown choice point for hybrid best first search!!!" << endl;
//             throw InternalError();
//         }
//     }
// }

void Solver::restore(CPStore& cp, OpenNode nd)
{
    if (ToulBar2::verbose >= 1) {
        if (wcsp->getTreeDec()) {
            cout << "[C" << wcsp->getTreeDec()->getCurrentCluster()->getId() << "] restore open node " << nd.getCost(wcsp->getTreeDec()->getCurrentCluster()->getCurrentDeltaUb()) << " (" << nd.first << ", " << nd.last << ")" << endl;
        } else {
            cout << "restore open node " << nd.getCost(MIN_COST) << " (" << nd.first << ", " << nd.last << ")" << endl;
        }
    }
    assert(nd.last >= nd.first);

    ptrdiff_t maxsize = nd.last - nd.first;
    if (maxsize == 0) {
        wcsp->enforceUb();
        assert(wcsp->isactivatePropagate());
        wcsp->propagate();
        return;
    }

    wcsp->deactivatePropagate();
    nbRecomputationNodes += maxsize;
    ChoicePoint* permute[maxsize];
    int assignLS[maxsize];
    Value valueLS[maxsize];
    unsigned int size = 0;
    bool randomOrder = (abs(ToulBar2::constrOrdering) == CONSTR_ORDER_RANDOM);
    if (randomOrder) {
        unsigned int pos = 0;
        for (ptrdiff_t idx = nd.first; idx < nd.last; ++idx) {
            permute[pos] = &cp[idx];
            ++pos;
        }
        assert(pos == maxsize);
        shuffle(&permute[0], &permute[maxsize], myrandom_generator);
    }
    for (ptrdiff_t idx = nd.first; idx < nd.last; ++idx) {
        assert((size_t)idx < cp.size());
        ChoicePoint* cp_ptr = ((randomOrder) ? permute[idx - nd.first] : &cp[idx]);
        assert(!wcsp->getTreeDec() || wcsp->getTreeDec()->getCurrentCluster()->isVar(cp_ptr->varIndex) || (wcsp->getTreeDec()->getCurrentCluster()->getFreedom() && wcsp->getTreeDec()->getCurrentCluster()->isVarTree(cp_ptr->varIndex)));
        if ((cp_ptr->op == CP_ASSIGN && !(cp_ptr->reverse && cp_ptr != &cp[nd.last - 1])) || (cp_ptr->op == CP_REMOVE && cp_ptr->reverse && cp_ptr != &cp[nd.last - 1])) {
            assignLS[size] = cp_ptr->varIndex;
            valueLS[size] = cp_ptr->value;
            size++;
        }
    }
    wcsp->enforceUb();
    wcsp->assignLS(assignLS, valueLS, size, false); // fast multiple assignments
    for (ptrdiff_t idx = nd.first; idx < nd.last; ++idx) {
        assert((size_t)idx < cp.size());
        ChoicePoint* cp_ptr = ((randomOrder) ? permute[idx - nd.first] : &cp[idx]);
        if (ToulBar2::verbose >= 1)
            cout << "retrieve choice point " << CPOperation[cp_ptr->op] << ((cp_ptr->reverse) ? "*" : "") << " (" << wcsp->getName(cp_ptr->varIndex) << ", " << cp_ptr->value << ") at position " << (cp_ptr - &cp[nd.first]) << endl;
        if (ToulBar2::verbose >= 1)
            cout << *((WCSP*)wcsp)->getVar(cp_ptr->varIndex) << endl;
        nbNodes++;
        switch (cp_ptr->op) { // TODO: some operations (remove,increase,decrease) are useless because of all assigns previously done
        case CP_ASSIGN: {
            if (cp_ptr->reverse && cp_ptr != &cp[nd.last - 1]) {
                wcsp->remove(cp_ptr->varIndex, cp_ptr->value);
                addChoicePoint(CP_REMOVE, cp_ptr->varIndex, cp_ptr->value, false);
            } else
                addChoicePoint(CP_ASSIGN, cp_ptr->varIndex, cp_ptr->value, false);
            break;
        }
        case CP_REMOVE: {
            if (cp_ptr->reverse && cp_ptr != &cp[nd.last - 1]) {
                addChoicePoint(CP_ASSIGN, cp_ptr->varIndex, cp_ptr->value, false);
            } else {
                wcsp->remove(cp_ptr->varIndex, cp_ptr->value);
                addChoicePoint(CP_REMOVE, cp_ptr->varIndex, cp_ptr->value, false);
            }
            break;
        }
        case CP_INCREASE: {
            if (cp_ptr->reverse && cp_ptr != &cp[nd.last - 1]) {
                wcsp->decrease(cp_ptr->varIndex, cp_ptr->value - 1);
                addChoicePoint(CP_DECREASE, cp_ptr->varIndex, cp_ptr->value - 1, false);
            } else {
                wcsp->increase(cp_ptr->varIndex, cp_ptr->value);
                addChoicePoint(CP_INCREASE, cp_ptr->varIndex, cp_ptr->value, false);
            }
            break;
        }
        case CP_DECREASE: {
            if (cp_ptr->reverse && cp_ptr != &cp[nd.last - 1]) {
                wcsp->increase(cp_ptr->varIndex, cp_ptr->value + 1);
                addChoicePoint(CP_INCREASE, cp_ptr->varIndex, cp_ptr->value + 1, false);
            } else {
                wcsp->decrease(cp_ptr->varIndex, cp_ptr->value);
                addChoicePoint(CP_DECREASE, cp_ptr->varIndex, cp_ptr->value, false);
            }
            break;
        }
        default: {
            cerr << "unknown choice point for hybrid best first search!!!" << endl;
            throw InternalError();
        }
        }
    }
    wcsp->reactivatePropagate();
    wcsp->propagate();
    // if (wcsp->getLb() != nd.getCost(((wcsp->getTreeDec())?wcsp->getTreeDec()->getCurrentCluster()->getCurrentDeltaUb():MIN_COST))) cout << "***** node cost: " << nd.getCost(((wcsp->getTreeDec())?wcsp->getTreeDec()->getCurrentCluster()->getCurrentDeltaUb():MIN_COST)) << " but lb: " << wcsp->getLb() << endl;
}

/**
 *
 * \brief Convert a choice point operation into a symbol taking into account boolean reverse stuff
 *        with the following  convention:
 * ASSIGN:     =        ( equivalent to assign(variable, value) )
 * REMOVE:     #        ( equivalent to remove(variable, value) )
 * INCREASE:   >       (warning  : equivalent to increase(variable, value+1) )
 * DECREASE:   <       (warning  :  equivalent to decrease(variable, value-1) )
 */
char Solver::opSymbol(const CPStore& cp, const ptrdiff_t idx, OpenNode nd)
{
    switch (cp[idx].op) {
    case CP_ASSIGN: {
        if (cp[idx].reverse && idx < nd.last - 1) {
            return '#';
        } else
            return '=';
        break;
    }
    case CP_REMOVE: {
        if (cp[idx].reverse && idx < nd.last - 1) {
            return '=';
        } else {
            return '#';
        }
        break;
    }
    case CP_INCREASE: {
        if (cp[idx].reverse && idx < nd.last - 1) {
            return '<';
        } else {
            return '>';
        }
        break;
    }
    case CP_DECREASE: {
        if (cp[idx].reverse && idx < nd.last - 1) {
            return '>';
        } else {
            return '<';
        }
        break;
    }
    default: {
        cerr << "unknown choice point for hybrid best first search!!!" << endl;
        throw InternalError();
    }
    }
}

void Solver::epsDumpSubProblems(CPStore& cp, OpenList& open)
{
    ofstream epsfile(ToulBar2::epsFilename);
    if (!epsfile) {
        cerr << "Cannot open file for EPS! " << ToulBar2::epsFilename << endl;
        throw WrongFileFormat();
    }
    Long nbsp = 0;
    while (!open.finished()) {
        OpenNode nd = open.top();
        open.pop();
        if (nd.getCost() < wcsp->getUb()) {
            string epsSubProblem = "-x=\"";
            for (ptrdiff_t idx = nd.first; idx < nd.last; ++idx) {
                epsSubProblem += to_string(",") + to_string(cp[idx].varIndex) + opSymbol(cp, idx, nd) + to_string(cp[idx].value + ((cp[idx].op == CP_INCREASE) ? -1 : 0) + ((cp[idx].op == CP_DECREASE) ? 1 : 0));
            }
            epsfile << epsSubProblem << "\""
                    << " -best=" << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->Cost2ADCost(nd.getCost()) << std::setprecision(DECIMAL_POINT) << endl;
            nbsp++;
        }
    }
    cout << endl
         << "Generating " << nbsp << " subproblems for EPS done in file " << ToulBar2::epsFilename << "." << endl
         << "Run them in parallel using " << std::thread::hardware_concurrency() << " cores, e.g. (*check pathname for original problem file*):" << endl
         << "./misc/script/eps.sh " << std::thread::hardware_concurrency() << " ./" << ToulBar2::epsFilename << " ./toulbar2 ./" << wcsp->getName() << " -ub=" << wcsp->getUb() << endl
         << endl;
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
    // The SolutionTrie computes solution Prefix Tree, in divVariables order.
    //  To merge equivalent nodes, we need the suffix tree, so we reverse the variables order.
    //  variables
    vector<Variable*> varReverse;
    for (unsigned v = 0; v < wcsp->getDivVariables().size(); v++) {
        varReverse.push_back(wcsp->getDivVariables()[wcsp->getDivVariables().size() - v - 1]);
    }
    int nLayers = varReverse.size();
    Mdd mdd(nLayers);
    vector<int> layerWidth;
    layerWidth.push_back(1);
    // solTrie
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

    // for relaxation
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
                            unsigned nodep_index = std::distance(nodesAtLayer[layer].begin(), nodep_it);
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
                // merge nodes that won't lead to a satisfying solution:
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

            // Merging nodes:TODO
            // Computing new state for merged nodes
            vector<int> newCount(nodesAtLayer[layer + 1].size(), -1);

            vector<int> newTarget(nextDistCounts.size(), -1); // vector with new state nodes ids
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
            // The nodes need to be renumbered - we want nodeids = 0, 1 , ... , divWidth
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
            // redirecting arcs in mdd[layer] from each source to new targets
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
            // Computing alphap[layer+1]
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

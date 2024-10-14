#include <cerrno>
#include <stdio.h>
#include <list>
#include <vector>
#include <string>
#include <set>
#include <algorithm>

using namespace std;
#include <iostream>
#include <fstream>

#include "core/tb2naryconstr.hpp"
#include "core/tb2wcsp.hpp"
#include "search/tb2solver.hpp"

#include "incop.h"
#include "incoputil.h"
#include "csproblem.h"
#include "narycsproblem.h"
#include "autotuning2.h"

extern ofstream* ofile; // le fichier de sortie

extern Stat_GWW* Statistiques;

INCOP::NaryCSProblem::NaryCSProblem(int nbvar, int nbconst)
    : CSProblem(nbvar, nbconst)
{
    ;
}

INCOP::NaryConstraint::NaryConstraint(int arit)
{
    arity = arit;
    constr = NULL;
}
INCOP::NaryConstraint::NaryConstraint(int arit, Constraint* ctr)
{
    arity = arit;
    constr = ctr;
}

INCOP::NaryVariable::NaryVariable() { ; }

/** code optimisé pour configuration semi-incrementale IncrCSPConfiguration*/
/*
int NaryCSProblem::move_evaluation
                    (Configuration* configuration,Move* move)
{int var_changee = ((CSPMove*)move)-> variable;
 int val_changee = ((CSPMove*)move)-> value;
 return(configuration->valuation
       +compute_conflict(configuration,var_changee,val_changee)
               -((IncrCSPConfiguration*)configuration)->tabconflicts[var_changee]);
}
*/

// optimisation pour IncrCSPConfiguration
/*
int NaryCSProblem::config_evaluation(Configuration* configuration)
{
  configuration->init_conflicts();
  int value=0;
  for (int i=0; i< (int) naryconstraints->size();i++)
    {int nbconf = (*naryconstraints)[i]->constraint_value (configuration);
      value+= nbconf;
     for (int j=0 ; j< (*naryconstraints)[i]->arity ; j++)
       ((IncrCSPConfiguration*)configuration)->tabconflicts[(*naryconstraints)[i]->constrainedvariables[j]] +=nbconf;
    }
  return value;
}
*/

Long INCOP::NaryCSProblem::config_evaluation(Configuration* configuration)
{
    configuration->init_conflicts();
    Long value = 0;
    for (int i = 0; i < (int)naryconstraints->size(); i++) {
        Long nbconf = (*naryconstraints)[i]->constraint_value(configuration);
        value += nbconf;
    }
    for (int i = 0; i < nbvar; i++)
        for (int j = 0; j < variable_domainsize(i); j++)
            configuration->incr_conflicts(i, j, j, compute_conflict(configuration, i, j));
    return value;
}

Long INCOP::NaryConstraint::constraint_value(Configuration* configuration)
{
    if (!constr) {
        int index = 0;
        for (int i = 0; i < arity; i++)
            index += configuration->config[constrainedvariables[i]] * multiplyers[i];
        return tuplevalues[index];
    } else {
        for (int i = 0; i < arity; i++) {
            tuple[i] = index2index[i][configuration->config[constrainedvariables[i]]];
        }
        assert(constr->isNary() || !constr->extension());
        return ((AbstractNaryConstraint*)constr)->eval(tuple);
    }
}

void INCOP::NaryCSProblem::incr_update_conflicts(IncrCSPConfiguration* configuration, Move* move)
{
    int var = ((CSPMove*)move)->variable;
    int value = ((CSPMove*)move)->value;
    int aval = configuration->config[var];
    Long actvalue, nctvalue;
    NaryVariable* varobjct = (*naryvariables)[var];
    INCOP::NaryConstraint* ct;
    for (int i = 0; i < (int)(varobjct->constraints).size(); i++) {
        ct = (varobjct->constraints)[i];
        actvalue = ct->constraint_value(configuration);
        configuration->config[var] = value;
        nctvalue = ct->constraint_value(configuration);
        configuration->config[var] = aval;
        for (int j = 0; j < ct->arity; j++) {
            configuration->tabconflicts[ct->constrainedvariables[j]] += nctvalue - actvalue;
        }
    }
}

void INCOP::NaryCSProblem::fullincr_update_conflicts(FullincrCSPConfiguration* configuration, Move* move)
{
    int var = ((CSPMove*)move)->variable;
    int value = ((CSPMove*)move)->value;
    int aval = configuration->config[var];
    Long actvalue, nctvalue;
    int var1, aval1;
    INCOP::NaryVariable* varobjct = (*naryvariables)[var];
    INCOP::NaryConstraint* ct;
    for (int i = 0; i < (int)(varobjct->constraints).size(); i++) {
        ct = (varobjct->constraints)[i];
        for (int j = 0; j < ct->arity; j++) {
            var1 = ct->constrainedvariables[j];
            if (var1 != var) {
                aval1 = configuration->config[var1];
                for (int k = 0; k < variable_domainsize(var1); k++) {
                    configuration->config[var1] = k;
                    actvalue = ct->constraint_value(configuration);
                    configuration->config[var] = value;
                    nctvalue = ct->constraint_value(configuration);
                    configuration->config[var] = aval;
                    //	       configuration->incr_conflicts(var1,k,k,nctvalue - actvalue);
                    configuration->tabconflicts[var1][k] += nctvalue - actvalue;
                }
                configuration->config[var1] = aval1;
            }
        }
    }
}

int INCOP::NaryConstraint::compute_index(int* values, vector<Value>* tabdomaines)
{
    int index = 0;
    for (int i = 0; i < arity; i++)
        index += compute_indexpart(i, values[i], tabdomaines);
    return index;
}

int INCOP::NaryConstraint::compute_indexpart(int i, int vali, vector<Value>* tabdomaines)
{
    int factor = 1;
    for (int j = i + 1; j < arity; j++) {
        factor = factor * tabdomaines[constrainedvariables[j]].size();
    }
    return vali * factor;
}

/** nombre de n-uplets d'une contrainte */
/* number of tuples of a constraint */
int INCOP::NaryConstraint::nbtuples(vector<Value>* tabdomaines)
{
    int nbtuples = 1;
    for (int j = 0; j < arity; j++) {
        nbtuples = nbtuples * tabdomaines[constrainedvariables[j]].size();
    }
    return nbtuples;
}

void INCOP::NaryConstraint::compute_indexmultiplyers(vector<Value>* tabdomaines)
{
    for (int i = 0; i < arity; i++)
        multiplyers.push_back(compute_indexmultiplyer(i, tabdomaines));
}

int INCOP::NaryConstraint::compute_indexmultiplyer(int i, vector<Value>* tabdomaines)
{
    int factor = 1;
    for (int j = i + 1; j < arity; j++) {
        factor = factor * tabdomaines[constrainedvariables[j]].size();
    }
    return factor;
}

/** calcul du nombre de conflits d'une affectation - appele par l'évaluation d'un mouvement (cas incr)*/

Long INCOP::NaryCSProblem::compute_conflict(Configuration* configuration, int var, int val)
{
    Long value = 0;
    int aval = configuration->config[var];
    configuration->config[var] = val;
    for (int i = 0; i < (int)(*naryvariables)[var]->constraints.size(); i++)
        value += (*naryvariables)[var]->constraints[i]->constraint_value(configuration);
    configuration->config[var] = aval;
    return value;
}

/** utilisation des configurations "semi-incrementales"IncrCSPConfiguration - les conflits des valeurs courantes des variables
    sont stockés dans le tableau tabconflicts
    ou tout-incrémentales  FullincrCSPConfiguration  : les conflits de toutes les valeurs avec la configuration courante
sont maintenus dans tabconflicts */
Configuration* INCOP::NaryCSProblem::create_configuration()
{
    return (new FullincrCSPConfiguration(nbvar, domainsize));
    // return (new IncrCSPConfiguration(nbvar,domainsize));
}

INCOP::NaryCSProblem* weighted_narycsp_creation(int nbvar, int nbconst, int maxdomsize,
    vector<INCOP::NaryVariable*>* vv, vector<INCOP::NaryConstraint*>* vct)
{
    INCOP::NaryCSProblem* p1 = new INCOP::NaryCSProblem(nbvar, nbconst);
    p1->domainsize = maxdomsize;
    p1->naryconstraints = vct;
    p1->naryvariables = vv;
    return p1;
}

/** lecture du debut du fichier : le probleme et les variables */
void wcspdomaines_file_read(WCSP* wcsp, int nbvar, vector<Value>* tabdomaines, vector<Value>& initsolution, vector<int>& initconfig)
{
    assert(initsolution.size() == wcsp->numberOfVariables());
    assert(initconfig.size() == (unsigned int)nbvar);
    int size = 0;
    for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++)
        if (wcsp->unassigned(i)) {
            for (EnumeratedVariable::iterator it = ((EnumeratedVariable*)wcsp->getVar(i))->begin(); it != ((EnumeratedVariable*)wcsp->getVar(i))->end(); ++it) {
                if (initsolution[i] == *it)
                    initconfig[size] = tabdomaines[size].size();
                tabdomaines[size].push_back(*it);
            }
            size++;
        }
    assert(size == nbvar);
}

/** lecture des contraintes */
int wcspdata_constraint_read(WCSP* wcsp, int nbconst, vector<INCOP::NaryVariable*>* vv, vector<INCOP::NaryConstraint*>* vct,
    vector<int>* connexions, vector<Value>* tabdomaines)
{
    static Tuple tuple;
    assert(wcsp->getUb() > wcsp->getLb());
    Cost gap = wcsp->getUb() - wcsp->getLb();
    int nbconst_ = 0;
    for (unsigned int i = 0; i < wcsp->numberOfConstraints(); i++) {
        if (wcsp->getCtr(i)->connected() && !wcsp->getCtr(i)->isSep() && wcsp->getCtr(i)->arity() <= ToulBar2::preprocessNary && wcsp->getCtr(i)->getDomainSizeProduct() <= MAX_NB_TUPLES) {
            int arity = 0;
            for (int j = 0; j < wcsp->getCtr(i)->arity(); j++)
                if (wcsp->getCtr(i)->getVar(j)->unassigned())
                    arity++;
            INCOP::NaryConstraint* ct = new INCOP::NaryConstraint(arity);
            vct->push_back(ct);
            int numvar = 0;
            for (int j = 0; j < wcsp->getCtr(i)->arity(); j++)
                if (wcsp->getCtr(i)->getVar(j)->unassigned()) {
                    numvar = wcsp->getCtr(i)->getVar(j)->getCurrentVarId();
                    ct->constrainedvariables.push_back(numvar);
                    (*vv)[numvar]->constraints.push_back(ct);
                }
            assert(ct->constrainedvariables.size() == (unsigned int)arity);
            ct->compute_indexmultiplyers(tabdomaines);
            Cost cost;
            wcsp->getCtr(i)->firstlex();
            while (wcsp->getCtr(i)->nextlex(tuple, cost)) {
                ct->tuplevalues.push_back(min(gap, cost));
            }
            nbconst_++;
        } else if (wcsp->getCtr(i)->connected() && !wcsp->getCtr(i)->isSep() && (wcsp->getCtr(i)->isNary() || !wcsp->getCtr(i)->extension())) {
            int arity = 0;
            for (int j = 0; j < wcsp->getCtr(i)->arity(); j++)
                if (wcsp->getCtr(i)->getVar(j)->unassigned())
                    arity++;
            INCOP::NaryConstraint* ct = new INCOP::NaryConstraint(arity, wcsp->getCtr(i));
            vct->push_back(ct);
            int numvar = 0;
            for (int j = 0; j < wcsp->getCtr(i)->arity(); j++) {
                EnumeratedVariable* myvar = (EnumeratedVariable*)wcsp->getCtr(i)->getVar(j);
                if (myvar->unassigned()) {
                    numvar = myvar->getCurrentVarId();
                    ct->constrainedvariables.push_back(numvar);
                    (*vv)[numvar]->constraints.push_back(ct);
                    ct->index2index.push_back(vector<tValue>());
                    for (auto iter = myvar->begin(); iter != myvar->end(); ++iter) {
                        ct->index2index.back().push_back(myvar->toIndex(*iter));
                    }
                }
            }
            assert(ct->constrainedvariables.size() == (unsigned int)arity);
            assert(ct->index2index.size() == (unsigned int)arity);
            ct->tuple.resize(arity, 0);
            nbconst_++;
        }
    }
    for (int i = 0; i < wcsp->getElimBinOrder(); i++) {
        Constraint* ctr = wcsp->getElimBinCtr(i);
        if (ctr->connected() && !ctr->isSep()) {
            int arity = 2;
            int numvar = 0;
            INCOP::NaryConstraint* ct = new INCOP::NaryConstraint(arity);
            vct->push_back(ct);
            for (int j = 0; j < arity; j++)
                if (ctr->getVar(j)->unassigned()) {
                    numvar = ctr->getVar(j)->getCurrentVarId();
                    ct->constrainedvariables.push_back(numvar);
                    (*vv)[numvar]->constraints.push_back(ct);
                }
            ct->compute_indexmultiplyers(tabdomaines);
            Cost cost;
            ctr->firstlex();
            while (ctr->nextlex(tuple, cost)) {
                ct->tuplevalues.push_back(min(gap, cost));
            }
            nbconst_++;
        }
    }
    for (int i = 0; i < wcsp->getElimTernOrder(); i++) {
        Constraint* ctr = wcsp->getElimTernCtr(i);
        if (ctr->connected() && !ctr->isSep()) {
            int arity = 3;
            int numvar = 0;
            INCOP::NaryConstraint* ct = new INCOP::NaryConstraint(arity);
            vct->push_back(ct);
            for (int j = 0; j < arity; j++)
                if (ctr->getVar(j)->unassigned()) {
                    numvar = ctr->getVar(j)->getCurrentVarId();
                    ct->constrainedvariables.push_back(numvar);
                    (*vv)[numvar]->constraints.push_back(ct);
                }
            ct->compute_indexmultiplyers(tabdomaines);
            Cost cost;
            ctr->firstlex();
            while (ctr->nextlex(tuple, cost)) {
                ct->tuplevalues.push_back(min(gap, cost));
            }
            nbconst_++;
        }
    }
    for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
        if (wcsp->unassigned(i) && wcsp->getMaxUnaryCost(i) > MIN_COST) {
            int arity = 1;
            int numvar = 0;
            INCOP::NaryConstraint* ct = new INCOP::NaryConstraint(arity);
            vct->push_back(ct);
            numvar = wcsp->getVar(i)->getCurrentVarId();
            ct->constrainedvariables.push_back(numvar);
            (*vv)[numvar]->constraints.push_back(ct);
            ct->compute_indexmultiplyers(tabdomaines);
            for (EnumeratedVariable::iterator it = ((EnumeratedVariable*)wcsp->getVar(i))->begin(); it != ((EnumeratedVariable*)wcsp->getVar(i))->end(); ++it) {
                ct->tuplevalues.push_back(min(gap, wcsp->getUnaryCost(i, *it)));
            }
            nbconst_++;
        }
    }
    assert(nbconst_ <= nbconst);
    return nbconst_;
}

int split(char* str, char c, char*** arr)
{
    int count = 1;
    int token_len = 1;
    int i = 0;
    char* p;
    char* t;

    p = str;
    while (*p != '\0') {
        if (*p == c)
            count++;
        p++;
    }

    *arr = (char**)malloc(sizeof(char*) * count);
    if (*arr == NULL)
        throw InternalError();

    p = str;
    while (*p != '\0') {
        if (*p == c) {
            (*arr)[i] = (char*)malloc(sizeof(char) * token_len);
            if ((*arr)[i] == NULL)
                throw InternalError();

            token_len = 0;
            i++;
        }
        p++;
        token_len++;
    }
    (*arr)[i] = (char*)malloc(sizeof(char) * token_len);
    if ((*arr)[i] == NULL)
        throw InternalError();

    i = 0;
    p = str;
    t = ((*arr)[i]);
    while (*p != '\0') {
        if (*p != c) {
            *t = *p;
            t++;
        } else {
            *t = '\0';
            i++;
            t = ((*arr)[i]);
        }
        p++;
    }
    *t = '\0';
    return count;
}

void removeSpaces(string& str)
{
    /* remove multiple spaces */
    int k = 0;
    for (unsigned int j = 0; j < str.size(); ++j) {
        if ((str[j] != ' ') || (str[j] == ' ' && str[j + 1] != ' ')) {
            str[k] = str[j];
            ++k;
        }
    }
    str.resize(k);

    /* remove space at the end */
    if (str[k - 1] == ' ')
        str.erase(str.end() - 1);
    /* remove space at the begin */
    if (str[0] == ' ')
        str.erase(str.begin());
}

/// \brief solves the current problem using INCOP local search solver by Bertrand Neveu
/// \return best solution cost found
/// \param cmd command line argument for narycsp INCOP local search solver (cmd format: lowerbound randomseed nbiterations method nbmoves neighborhoodchoice neighborhoodchoice2 minnbneighbors maxnbneighbors  neighborhoodchoice3 autotuning tracemode)
/// \param solution best solution assignment found (MUST BE INITIALIZED WITH A DEFAULT ASSIGNMENT)
Cost Solver::narycsp(string cmd, vector<Value>& bestsolution)
{
    Long result = MAX_COST;
    Cost initialUpperBound = wcsp->getUb();
    if (wcsp->isKnapsack() || wcsp->isGlobal()) {
        wcsp->setUb(MAX_COST / MEDIUM_COST / MEDIUM_COST / MEDIUM_COST);
    }

    string filename = "/dev/stdin";
    string outputfile = "/dev/stdout";
    int verbose = ToulBar2::verbose;
    char line[1024];
    char** argv = NULL;
    int tuningmode = 0; // no automatic tuning
    int argc = 0;

    // remove leading space  from incop command line
    cmd.erase(0, cmd.find_first_not_of(" \t\n\r\f\v"));

    // remove multiples space in cmd
    removeSpaces(cmd);

    sprintf(line, "bin/Linux/narycsp %s %s %s", outputfile.c_str(), filename.c_str(), cmd.c_str());

    argc = split(line, ' ', &argv);

    if (verbose > 0) {
        cout << "---------------------------" << endl;
        cout << "number of arguments for narycsp: " << argc << endl;
        cout << "---------------------------" << endl;
        for (int i = 0; i < argc; i++)
            cout << "arg #" << i << " --> " << argv[i] << endl;
        if (ToulBar2::verbose >= 3)
            cout << *wcsp;
    }

    INCOP::NaryCSProblem* problem; // pointeur sur le probleme

    // les divers arguments lus dans la ligne de commande
    int nbvar, nbconst, domsize;
    Long lbound = 0;
    int taille, nbessais;
    int graine1;
    int narg = 2; // compteur des arguments

    arguments_borneinf(argv, narg, lbound);
    // lecture des paramètres de l'algo et création de l'objet algo
    IncompleteAlgorithm* algo = algo_creation(argv, narg, taille, graine1, nbessais);

    // allocation de l'objet pour les stats
    Statistiques = new Stat_GWW(1, nbessais);

    // argument pour la trace
    arguments_tracemode(argv, narg);
    // pour la recuperation du signal 10
    //  sigaction();

    // argument de temps maximum
    double maxtime;
    if (tuningmode)
        arguments_tempscpu(argv, narg, maxtime);

    // Declaration des variables contenant les structures de données des problemes
    string pbname;
    Long upperbound;
    pbname = wcsp->getName();
    nbvar = wcsp->numberOfUnassignedVariables();
    domsize = 0;
    int nbunarycosts = 0;
    vector<int> tabvars;
    for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
        if (wcsp->unassigned(i)) {
            assert(wcsp->enumerated(i));
            tabvars.push_back(i);
            if ((int)wcsp->getDomainSize(i) > domsize)
                domsize = wcsp->getDomainSize(i);
            if (wcsp->getMaxUnaryCost(i) > MIN_COST)
                nbunarycosts++;
        }
    }
    nbconst = wcsp->numberOfConnectedConstraints() + nbunarycosts;
    upperbound = wcsp->getUb();

    vector<int> initconfig(nbvar, 0);
    vector<Value>* tabdomaines; // les différents types de domaines
    tabdomaines = new vector<Value>[nbvar];
    wcspdomaines_file_read((WCSP*)wcsp, nbvar, tabdomaines, bestsolution, initconfig);

    int domaines[nbvar]; // 1 domaine par variable
    for (int i = 0; i < nbvar; i++) {
        domaines[i] = i;
    }

    // Initialisation des structures de données des problémes
    vector<INCOP::NaryConstraint*> constraints;
    vector<INCOP::NaryVariable*> variables;
    vector<int>* connexions;
    connexions = new vector<int>[nbvar];

    for (int i = 0; i < nbvar; i++) {
        INCOP::NaryVariable* nv = new INCOP::NaryVariable();
        variables.push_back(nv);
    }

    nbconst = wcspdata_constraint_read((WCSP*)wcsp, nbconst, &variables, &constraints, connexions, tabdomaines);
    int pbnumber = 0;
    Statistiques->init_pb(pbnumber);
    problem = weighted_narycsp_creation(nbvar, nbconst, domsize, &variables, &constraints);

    problem->lower_bound = lbound;
    // mise en place des domaines
    problem->set_domains_connections(domaines, tabdomaines, connexions);

    // creation de la population et initialisation
    // La population : tableau de configurations
    Configuration* population[taille];

    problem->init_population(population, taille);

    problem->allocate_moves();

    if (tuningmode)
        autosolving((LSAlgorithm*)algo, population, problem, 0, graine1, nbessais, maxtime, 1000000);
    else {
        // boucle sur les essais
        for (int nessai = 0; nessai < nbessais; nessai++) {
            if (ToulBar2::interrupted) {
                throw TimeOut();
            }
            executer_essai(problem, algo, population, taille, graine1, nessai, &initconfig);
            if (wcsp->getLb() + problem->best_config->valuation < min(upperbound, initialUpperBound)) {
                Cost previousUb = MIN_COST;
                int depth = Store::getDepth();
                try {
                    Store::store();
                    vector<Value> solution(problem->best_config->nbvar);
                    for (int i = 0; i < problem->best_config->nbvar; i++) {
                        solution[i] = tabdomaines[i][problem->best_config->config[i]];
                    }
                    if (wcsp->getUb() > initialUpperBound) {
                        previousUb = wcsp->getUb();
                        wcsp->setUb(initialUpperBound);
                        wcsp->enforceUb();
                    }
                    wcsp->assignLS(tabvars, solution);
                    newSolution();
                    result = wcsp->getUb();
                    upperbound = result;
                    for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
                        bestsolution[i] = wcsp->getValue(i);
                        wcsp->setBestValue(i, bestsolution[i]);
                    }
                } catch (const Contradiction&) {
                    wcsp->whenContradiction();
                }
                Store::restore(depth);
                if (previousUb > MIN_COST) {
                    wcsp->setUb(previousUb);
                }
            }
        }
        // ecriture statistiques
        Statistiques->current_try++;
        //      ecriture_stat_probleme();
    }
   
    for(size_t pop_id = 0; pop_id < static_cast<size_t>(taille); pop_id ++) {
        delete population[pop_id];
    }

    delete problem;
    delete[] tabdomaines;
    delete[] connexions;

    for(auto& cst: constraints) {
        delete cst;
    }
    for(auto& var: variables) {
        delete var;
    }

    delete Statistiques;

    delete algo;

    if (wcsp->getUb() < initialUpperBound) {
        wcsp->enforceUb();
        wcsp->propagate();
    } else {
        assert(wcsp->getUb() == initialUpperBound || wcsp->isKnapsack() || wcsp->isGlobal());
        wcsp->setUb(initialUpperBound);
    }

    // free the options string
    for(int i = 0; i < argc; i ++) {
        free(argv[i]);
    }
    free(argv);


    return result;
}

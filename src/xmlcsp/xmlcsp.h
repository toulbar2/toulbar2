

#include "core/tb2wcsp.hpp"
// We suppose that the XMLCSP library is placed in the directory
// xmlcsp from toulbar2, so tb2wcsp.hpp is in the parent directory

#include <iostream>
#include <fstream>
#include <list>
#include <map>
#include <vector>

#include "XMLParser_libxml2.hh"

using namespace CSPXMLParser;

/*
 * A simple demo to illustrate how to use the XML CSP parser 
 */

#define MAXDOMS 100000
#define MAXDOMSIZE MAX_DOMAIN_SIZE
#define MAXDOMSIZEZERO 1000
#define MAX_COST_XML MAX_COST

/**
 * The methods of this class will be called by the XML parser to
 * report the elements of the definition of the CSP instance to your
 * own solver.
 *
 * This sample callback merely prints the data it receives. You must
 * modify it to transmit the informations to your solver.
 *
 * The description of each method can be found in
 * C++/include/CSPParserCallback.hh
 */
class MyCallback : public CSPParserCallback {
public:
    WCSP* wcsp;

    bool convertWCSP;
    string fname;
    ofstream f;
    bool intension;

private:
    int nvars;
    int wcsp_nvars;
    map<int, string> varName;
    map<int, int> xml2wcspIndex;

    typedef struct {
        int arity;
        int def;
        list<int*> tuples;
        RelType type;
        int id;
        int defaultCost;
    } relation;

    typedef struct {
        int id;
        AST* exp;
        string expstr;
    } predicate;

    typedef struct {
        list<int> scope;
        int arity;
        int type;
        string ref;
    } ctr;

    map<string, relation*> rels;
    map<string, predicate*> preds;
    list<ctr*> ctrs;

    int currentDom;
    relation* currentR;
    predicate* currentP;
    ctr* currentC;
    int currentIndexArity;

    vector<int> nDoms;
    vector<vector<int> > DomsToIndex;

    int maxDomSize;

    bool firstDomainValue;
    bool firstTuple;
    string constraintReference;

    Cost initialLowerBound;

public:
    virtual void beginInstance(const string& name)
    {
        wcsp->Doms.resize(MAXDOMS);
        nDoms.resize(MAXDOMS);
        DomsToIndex.resize(MAXDOMS);

        for (int i = 0; i < MAXDOMS; i++)
            nDoms[i] = 0;
        maxDomSize = 0;

        wcsp->updateUb(MAX_COST_XML);
        initialLowerBound = MIN_COST;
        nvars = 0;
        wcsp_nvars = wcsp->numberOfVariables();
    }

    virtual void beginDomainsSection(int nbDomains)
    {
    }

    virtual void beginDomain(const string& name, int idDomain, int nbValue)
    {
        firstDomainValue = true;
        currentDom = idDomain;
        wcsp->Doms[currentDom].resize(nbValue);
        DomsToIndex[currentDom].resize(MAXDOMSIZE);
    }

    virtual void addDomainValue(int v)
    {
        wcsp->Doms[currentDom][nDoms[currentDom]] = v;
        DomsToIndex[currentDom][MAXDOMSIZEZERO + v] = nDoms[currentDom];
        nDoms[currentDom]++;
        if (!firstDomainValue)
            firstDomainValue = false;
    }

    virtual void addDomainValue(int first, int last)
    {
        for (int i = first; i <= last; i++)
            addDomainValue(i);
    }

    virtual void endDomain() {}
    virtual void endDomainsSection() {}

    virtual void beginVariablesSection(int nbVariables) {}

    virtual void addVariable(const string& name, int idVar,
        const string& domain, int idDomain)
    {
        int varindex = wcsp->getVarIndex(name);
        if (varindex == wcsp->numberOfVariables()) {
            varindex = wcsp_nvars;
            wcsp_nvars++;
        }
        wcsp->varsDom[varindex] = idDomain;
        varName[nvars] = name;
        xml2wcspIndex[nvars] = varindex;
        if (maxDomSize < nDoms[idDomain])
            maxDomSize = nDoms[idDomain];
        nvars++;
    }

    virtual void endVariablesSection() {}
    virtual void beginRelationsSection(int nbRelations) {}

    virtual void beginRelation(const string& name, int idRel,
        int arity, int nbTuples, RelType relType)
    {
        relation* r = new relation;
        r->id = idRel;
        r->arity = arity;
        r->type = relType;
        r->defaultCost = 0;
        currentR = r;
        rels[name] = r;

        switch (relType) {
        case REL_SUPPORT:
            break;
        case REL_CONFLICT:
            break;
        case REL_SOFT:
            break;
        default:
            throw runtime_error("unknown relation type");
        }

        firstTuple = true;
    }

    virtual void relationDefaultCost(const string& name, int idRel,
        int defaultCost)
    {
        currentR->defaultCost = defaultCost;
    }

    virtual void addRelationTuple(int arity, int tuple[])
    {
        firstTuple = false;
        int* t = new int[arity];
        for (int i = 0; i < arity; ++i)
            t[i] = tuple[i];
        currentR->tuples.push_back(t);
    }

    virtual void addRelationTuple(int arity, int tuple[], int cost)
    {
        firstTuple = false;
        int* t = new int[arity + 1];
        for (int i = 0; i < arity; ++i)
            t[i] = tuple[i];
        t[arity] = cost;
        currentR->tuples.push_back(t);
    }

    virtual void endRelation() {}
    virtual void endRelationsSection() {}

    virtual void beginPredicatesSection(int nbPredicates)
    {
        //cout << "<predicates nbPredicates='" << nbPredicates << "'>" << endl;
    }

    virtual void beginPredicate(const string& name, int idPred)
    {
        predicate* p = new predicate;
        p->id = idPred;
        p->exp = NULL;
        preds[name] = p;
        currentP = p;
        intension = true;
        cerr << "Warning!!! Predicate " << name << " not implemented, just skip it..." << endl;
    }

    virtual void addFormalParameter(int pos, const string& name,
        const string& type)
    {
        //cout << "   formal parameter " << pos << ": " << type << " " << name << endl;
    }

    virtual void predicateExpression(AST* tree)
    {
        //cout << "   predicate definition (AST) = ";
        //tree->prefixExpression(cout);
        //cout << endl;
    }

    virtual void predicateExpression(const string& expr)
    {
        currentP->expstr = expr;
    }

    virtual void endPredicate() {}
    virtual void endPredicatesSection() {}
    virtual void beginConstraintsSection(int nbConstraints) {}

    virtual void beginConstraint(const string& name, int idConstr,
        int arity,
        const string& reference,
        CSPDefinitionType type, int id,
        const ASTList& scope)
    {
        currentC = new ctr;
        currentC->arity = arity;
        currentC->type = type;
        currentC->ref = reference;

        for (int i = 0; i < scope.size(); i++) {
            currentC->scope.push_back(scope[i].getVarId());
        }

        constraintReference = reference;
    }

    virtual void constraintsMaximalCost(int maximalCost)
    {
        wcsp->updateUb(maximalCost);
    }

    virtual void constraintsInitialCost(int initialCost)
    {
        initialLowerBound = initialCost;
    }

    virtual void constraintParameters(const ASTList& args)
    {
        if (constraintReference == "global:cumulative") {
            const AST& tasks = args[0];
            const AST& limit = args[1];
            for (int i = 0; i < tasks.size(); ++i) {
                const AST& desc = tasks[i];
                if (desc.hasKey("origin")) {
                    if (desc["origin"].isVar()) {
                    } else if (desc["origin"].isInteger()) {
                    }
                }
                if (desc.hasKey("duration")) {
                    if (desc["duration"].isVar()) {
                    } else if (desc["duration"].isInteger()) {
                    }
                }
                if (desc.hasKey("end")) {
                    if (desc["end"].isVar()) {
                    } else if (desc["end"].isInteger()) {
                    }
                }
                if (desc.hasKey("height")) {
                    if (desc["height"].isVar()) {
                    } else if (desc["height"].isInteger()) {
                    }
                }
            }
            if (limit.isVar()) {
            } else if (limit.isInteger()) {
            }
        } else if (constraintReference == "global:element") {
            const AST& index = args[0];
            const AST& table = args[1];
            const AST& value = args[2];

            if (index.isVar()) {
            } else {
            }
            for (int i = 0; i < table.size(); ++i) {
                if (table[i].isVar()) {
                } else {
                }
            }
            if (value.isVar()) {
            } else {
            }
        } else if (constraintReference == "global:weightedsum") {
            /*const AST &sum=args[0];
		const AST &op=args[1];
		const AST &rhs=args[2];

		cout << sum[0]["coef"].getInteger() 
			 << "*" << sum[0]["var"].getVarName()
			 << showpos;
		for(int i=1;i<sum.size();++i)
		  cout << sum[i]["coef"].getInteger() 
			   << "*" << sum[i]["var"].getVarName();

		cout << noshowpos;
		op.infixExpression(cout);
		if (rhs.isVar()) {}
		else {}
		*/
        } else {
            //cout << "constraint parameters=";
            //args.postfixExpression(cout);
            //cout << endl;
        }
    }

    virtual void endConstraint()
    {
        ctrs.push_back(currentC);
    }

    /**
   * end the definition of all constraints
   */
    virtual void endConstraintsSection()
    {
    }

    /********************************************************************/

    void createWCSP()
    {
        int i = 0;
        map<int, string>::iterator it = varName.begin();
        while (it != varName.end()) {
            int domsize = nDoms[wcsp->varsDom[xml2wcspIndex[it->first]]];
            string varname = it->second;
            int varindex = wcsp->getVarIndex(varname);
            if (ToulBar2::verbose >= 3)
                cout << "read " << ((varindex < wcsp->numberOfVariables())?"known":"new") << " variable " << i << " of size " << domsize << endl;
            if (varindex == wcsp->numberOfVariables()) {
                if (domsize >= 0) {
                    varindex = wcsp->makeEnumeratedVariable(varname, 0, domsize - 1);
                    for (int idx = 0; idx < domsize; idx++) {
                        int v = wcsp->Doms[wcsp->varsDom[varindex]][idx];
                        ((EnumeratedVariable *) wcsp->getVar(varindex))->addValueName("v" + to_string(v));
                    }
                } else {
                    varindex = wcsp->makeIntervalVariable(varname, 0, -domsize - 1);
                }
            } else {
                if (wcsp->enumerated(varindex)) {
                    assert(domsize >= 0);
                    if (domsize != wcsp->getDomainInitSize(varindex)) {
                        cerr << "wrong domain size " << domsize << " compared to previous one " << wcsp->getDomainInitSize(varindex) << endl;
                        exit(EXIT_FAILURE);
                    }
                    for (int idx = 0; idx < domsize; idx++) {
                        int v = wcsp->Doms[wcsp->varsDom[varindex]][idx];
                        string vname = "v" + to_string(v);
                        if (((EnumeratedVariable *) wcsp->getVar(varindex))->getValueName(idx) != vname) {
                            cerr << "wrong domain value name " << vname << " compared to previous one " << ((EnumeratedVariable *) wcsp->getVar(varindex))->getValueName(idx) << endl;
                            exit(EXIT_FAILURE);
                        }
                    }
                } else {
                    wcsp->increase(varindex, 0);
                    wcsp->decrease(varindex, abs(domsize) - 1);
                }
            }
            assert(varindex == xml2wcspIndex[it->first]);
            ++it;
            i++;
        }

        unsigned int a;
        unsigned int b;
        unsigned int c;
        Cost inclowerbound = initialLowerBound;
        vector<TemporaryUnaryConstraint> unaryconstrs;
        int indexUnary = -1;

        list<ctr*>::iterator itc = ctrs.begin();
        while (itc != ctrs.end()) {
            ctr* cxml = *itc;
            int arity = cxml->arity;
            relation* r = NULL;

            if (cxml->type == RelationType) {
                EnumeratedVariable** scopeVar = new EnumeratedVariable*[arity];
                int* scopeIndex = new int[arity];
                int* values = new int[arity];
                Tuple strvalues(arity, 0);
                map<int, int> scopeOrder;

                int index = 0;
                list<int>::iterator its = cxml->scope.begin();
                while (its != cxml->scope.end()) {
                    int id = xml2wcspIndex[*its];
                    scopeOrder[id] = index;
                    ++its;
                    index++;
                }

                index = 0;
                map<int, int>::iterator ito = scopeOrder.begin();
                while (ito != scopeOrder.end()) {
                    int id = ito->first;
                    scopeIndex[index] = id;
                    scopeVar[index] = (EnumeratedVariable*)wcsp->getVar(id);
                    index++;
                    ++ito;
                }

                map<string, relation*>::iterator it = rels.find(cxml->ref);
                if (it != rels.end()) {
                    r = it->second;
                    int ntuples = r->tuples.size();
                    Cost defval;

                    if (r->type == REL_SUPPORT) {
#ifdef MAXCSP
                        defval = UNIT_COST;
#else
                        defval = MAX_COST_XML;
#endif
                    } else if (r->type == REL_SOFT) {
                        defval = r->defaultCost;
                    } else {
                        defval = MIN_COST;
                    }

                    int ctrIndex = -1;
                    Constraint* ctr = NULL;
                    if (arity > 3) {
                        ctrIndex = wcsp->postNaryConstraintBegin(scopeIndex, arity, defval, r->tuples.size());
                        ctr = wcsp->getCtr(ctrIndex);
                        list<int*>::iterator itl = r->tuples.begin();
                        while (itl != r->tuples.end()) {
                            int* t = *itl;
                            for (int i = 0; i < r->arity; i++) {
                                int pos = scopeOrder[scopeIndex[i]];
                                strvalues[i] = DomsToIndex[wcsp->varsDom[scopeIndex[i]]][MAXDOMSIZEZERO + t[pos]];
                            }
                            if (r->type == REL_SUPPORT) {
                                wcsp->postNaryConstraintTuple(ctrIndex, strvalues, MIN_COST);
                            } else if (r->type == REL_CONFLICT) {
#ifdef MAXCSP
                                wcsp->postNaryConstraintTuple(ctrIndex, strvalues, UNIT_COST);
#else
                                wcsp->postNaryConstraintTuple(ctrIndex, strvalues, MAX_COST_XML);
#endif
                            } else if (r->type == REL_SOFT) {
                                wcsp->postNaryConstraintTuple(ctrIndex, strvalues, t[arity]);
                            } else {
                                wcsp->postNaryConstraintTuple(ctrIndex, strvalues, MAX_COST_XML);
                            }

                            ++itl;
                        }
                        wcsp->postNaryConstraintEnd(ctrIndex);
                    } else if (arity == 3) {
                        int i = scopeIndex[0];
                        int j = scopeIndex[1];
                        int k = scopeIndex[2];
                        EnumeratedVariable* x = scopeVar[0];
                        EnumeratedVariable* y = scopeVar[1];
                        EnumeratedVariable* z = scopeVar[2];
                        if (ToulBar2::verbose >= 3)
                            cout << "read ternary constraint "
                                 << " on " << i << "," << j << "," << k << endl;
                        vector<Cost> costs;
                        for (a = 0; a < x->getDomainInitSize(); a++) {
                            for (b = 0; b < y->getDomainInitSize(); b++) {
                                for (c = 0; c < z->getDomainInitSize(); c++) {
                                    costs.push_back(defval);
                                }
                            }
                        }
                        list<int*>::iterator itl = r->tuples.begin();
                        while (itl != r->tuples.end()) {
                            int* t = *itl;
                            for (int i = 0; i < r->arity; i++) {
                                int pos = scopeOrder[scopeIndex[i]];
                                values[i] = DomsToIndex[wcsp->varsDom[scopeIndex[i]]][MAXDOMSIZEZERO + t[pos]];
                            }
                            int pos = values[0] * y->getDomainInitSize() * z->getDomainInitSize() + values[1] * z->getDomainInitSize() + values[2];
                            if (r->type == REL_SUPPORT) {
                                costs[pos] = MIN_COST;
                            } else if (r->type == REL_CONFLICT) {
#ifdef MAXCSP
                                costs[pos] = UNIT_COST;
#else
                                costs[pos] = MAX_COST_XML;
#endif
                            } else if (r->type == REL_SOFT) {
                                costs[pos] = t[arity];
                            } else {
                                costs[pos] = MAX_COST_XML;
                            }
                            ++itl;
                        }
                        if ((defval != MIN_COST) || (ntuples > 0))
                            ctrIndex = wcsp->postTernaryConstraint(i, j, k, costs);
                    } else if (arity == 2) {
                        int i = scopeIndex[0];
                        int j = scopeIndex[1];
                        EnumeratedVariable* x = scopeVar[0];
                        EnumeratedVariable* y = scopeVar[1];
                        if (ToulBar2::verbose >= 3)
                            cout << "read binary constraint "
                                 << " on " << i << "," << j << endl;
                        vector<Cost> costs;
                        for (a = 0; a < x->getDomainInitSize(); a++) {
                            for (b = 0; b < y->getDomainInitSize(); b++) {
                                costs.push_back(defval);
                            }
                        }
                        list<int*>::iterator itl = r->tuples.begin();
                        while (itl != r->tuples.end()) {
                            int* t = *itl;
                            for (int i = 0; i < r->arity; i++) {
                                int pos = scopeOrder[scopeIndex[i]];
                                values[i] = DomsToIndex[wcsp->varsDom[scopeIndex[i]]][MAXDOMSIZEZERO + t[pos]];
                            }
                            int pos = values[0] * y->getDomainInitSize() + values[1];
                            if (r->type == REL_SUPPORT) {
                                costs[pos] = MIN_COST;
                            } else if (r->type == REL_CONFLICT) {
#ifdef MAXCSP
                                costs[pos] = UNIT_COST;
#else
                                costs[pos] = MAX_COST_XML;
#endif
                            } else if (r->type == REL_SOFT) {
                                costs[pos] = t[arity];
                            } else {
                                costs[pos] = MAX_COST_XML;
                            }
                            ++itl;
                        }
                        if ((defval != MIN_COST) || (ntuples > 0))
                            ctrIndex = wcsp->postBinaryConstraint(i, j, costs);
                    } else if (arity == 1) {
                        int i = scopeIndex[0];
                        EnumeratedVariable* x = (EnumeratedVariable*)wcsp->getVar(i);
                        if (ToulBar2::verbose >= 3)
                            cout << "read unary constraint on " << i << endl;
                        if (x->enumerated()) {
                            TemporaryUnaryConstraint unaryconstr;
                            unaryconstr.var = x;
                            for (a = 0; a < x->getDomainInitSize(); a++) {
                                unaryconstr.costs.push_back(defval);
                            }
                            indexUnary = unaryconstrs.size();
                            unaryconstrs.push_back(unaryconstr);
                            x->queueNC();

                            list<int*>::iterator itl = r->tuples.begin();
                            while (itl != r->tuples.end()) {
                                int* t = *itl;
                                for (int i = 0; i < r->arity; i++) {
                                    int pos = scopeOrder[scopeIndex[i]];
                                    values[i] = DomsToIndex[wcsp->varsDom[scopeIndex[i]]][MAXDOMSIZEZERO + t[pos]];
                                }
                                if (r->type == REL_SUPPORT) {
                                    unaryconstrs[indexUnary].costs[values[0]] = MIN_COST;
                                } else if (r->type == REL_CONFLICT) {
#ifdef MAXCSP
                                    unaryconstrs[indexUnary].costs[values[0]] = UNIT_COST;
#else
                                    unaryconstrs[indexUnary].costs[values[0]] = MAX_COST_XML;
#endif
                                } else if (r->type == REL_SOFT) {
                                    unaryconstrs[indexUnary].costs[values[0]] = t[arity];
                                } else {
                                    unaryconstrs[indexUnary].costs[values[0]] = MAX_COST_XML;
                                }
                                ++itl;
                            }
                        }
                    } else if (arity == 0) {
                    }

                    if (ctrIndex >= 0)
                        ctr = wcsp->getCtr(ctrIndex);
                    if (ctr && ToulBar2::verbose >= 4)
                        cout << *ctr << endl;
                }
                delete[] scopeVar;
                delete[] scopeIndex;
                delete[] values;
            }
            if (cxml->type == PredicateType)
                intension = true;
            itc++;
        }

        map<string, relation*>::iterator itr = rels.begin();
        while (itr != rels.end()) {
            relation* r = itr->second;
            list<int*>::iterator itl = r->tuples.begin();
            while (itl != r->tuples.end()) {
                delete[] * itl;
                ++itl;
            }
            ++itr;
        }

        // apply basic initial propagation AFTER complete network loading
        wcsp->increaseLb(inclowerbound);

        for (unsigned int u = 0; u < unaryconstrs.size(); u++) {
            wcsp->postUnaryConstraint(unaryconstrs[u].var->wcspIndex, unaryconstrs[u].costs);
        }
        wcsp->sortConstraints();
        f.close();

        nDoms.clear();
        DomsToIndex.clear();
        rels.clear();
        preds.clear();
        ctrs.clear();
    }

    virtual void endInstance()
    {
        createWCSP();
    }
};

/*
 * ****** Global soft constraint representing a weighted CSP ********
 * 
 * Contains also ToulBar2 global variable definitions
 */

#include "tb2system.hpp"
#include "tb2wcsp.hpp"
#include "tb2enumvar.hpp"
#include "tb2intervar.hpp"
#include "tb2binconstr.hpp"
#include "tb2ternaryconstr.hpp"
#include "tb2naryconstr.hpp"
#include "tb2arithmetic.hpp"
#include "tb2pedigree.hpp"
#include "tb2vac.hpp"


/*
 * Global variables
 * 
 */
 
double ToulBar2::version  = 0.5;
int  ToulBar2::verbose  = 0;
bool ToulBar2::showSolutions  = false;
bool ToulBar2::writeSolution  = false;
int  ToulBar2::elimDegree = -1;
int  ToulBar2::elimDegree_preprocessing  = -1;
bool ToulBar2::preprocessTernaryHeuristic  = false;
bool ToulBar2::FDAComplexity = false;
bool ToulBar2::FDAC = false;
bool ToulBar2::dichotomicBranching = false;
bool ToulBar2::lds = false;
bool ToulBar2::limited = false;
unsigned int ToulBar2::dichotomicBranchingSize = 10;
bool ToulBar2::generation = false;
#ifdef MENDELSOFT
bool ToulBar2::binaryBranching = true;
bool ToulBar2::elimVarWithSmallDegree  = true;
bool ToulBar2::preprocessTernary  = true;
bool ToulBar2::lastConflict = true;
#else
bool ToulBar2::binaryBranching = false;
bool ToulBar2::preprocessTernary  = false;
bool ToulBar2::lastConflict = false;
#endif

externalevent ToulBar2::setvalue = NULL;
externalevent ToulBar2::setmin = NULL;
externalevent ToulBar2::setmax = NULL;
externalevent ToulBar2::removevalue = NULL;
externalcostevent ToulBar2::setminobj = NULL;
Pedigree *ToulBar2::pedigree = NULL;

bool ToulBar2::bayesian = false;
int ToulBar2::resolution = 7;
TProb ToulBar2::errorg = 0.05;
TProb ToulBar2::NormFactor = 1;
int ToulBar2::foundersprob_class = 0;    // 0: 			equal frequencies
										 // 1: 			probs depending on the frequencies
										 // otherwise:  read from file
bool ToulBar2::consecutiveAllele = false;

bool ToulBar2::vac = false;
bool ToulBar2::vacAlternative = false;
bool ToulBar2::vacDecomposition = false;

int WCSP::wcspCounter = 0;



/*
 * WCSP constructors
 * 
 */
 
 

WCSP::WCSP(Store *s, Cost upperBound) : 
 		storeData(s), 
        lb(0, &s->storeCost), 
        ub(upperBound),
        NCBucketSize(cost2log2(getUb()) + 1),
        NCBuckets(NCBucketSize, VariableList(&s->storeVariable)),
        lastConflictConstr(NULL),
		elimOrder(0, &s->storeValue),
		binaryOrder(0, &s->storeValue),
		ternaryOrder(0, &s->storeValue),
		naryOrder(0, &s->storeValue)
{ 
    objectiveChanged = false;
    nbNodes = 0;
    maxdomainsize = 0;
    isternary = false;
    instance = wcspCounter++;

    if (ToulBar2::vac) {
	  vac = new VACExtension(this, ToulBar2::vacAlternative, ToulBar2::vacDecomposition);
	} else {
	  vac = NULL;
	}
}


WCSP::~WCSP()
{
    for (unsigned int i=0; i<vars.size(); i++) delete vars[i];
    for (unsigned int i=0; i<constrs.size(); i++) delete constrs[i];
}

WeightedCSP *WeightedCSP::makeWeightedCSP(Store *s, Cost upperBound)
{
    WeightedCSP * W = new WCSP(s, upperBound);
    return W;
}

int WCSP::makeEnumeratedVariable(string n, Value iinf, Value isup)
{
    EnumeratedVariable *x;
    if (!ToulBar2::vac) {
      x = new EnumeratedVariable(this, n, iinf, isup);
    }
    else {
      x = new VACVariable(this, n, iinf, isup);
    }
    if(maxdomainsize < isup - iinf + 1) maxdomainsize = isup - iinf + 1;
    return x->wcspIndex;
}

int WCSP::makeEnumeratedVariable(string n, Value *d, int dsize)
{
    EnumeratedVariable *x;
    if (!ToulBar2::vac) {
      x = new EnumeratedVariable(this, n, d, dsize);
    }
    else {
      x = new VACVariable(this, n, d, dsize);
    }
    if(maxdomainsize < dsize) maxdomainsize = dsize;
    return x->wcspIndex;
}

int WCSP::makeIntervalVariable(string n, Value iinf, Value isup)
{
    IntervalVariable *x = new IntervalVariable(this, n, iinf, isup);
    return x->wcspIndex;
}


// postBinaryConstraint looks for an existing constraint 
// (even not connected constraints). It also alocates
// memory for a new constraints. For this two reasons it should 
// ONLY be called when creating the initial wcsp
int WCSP::postBinaryConstraint(int xIndex, int yIndex, vector<Cost> &costs)
{
	EnumeratedVariable* x =  (EnumeratedVariable *) vars[xIndex];
	EnumeratedVariable* y =  (EnumeratedVariable *) vars[yIndex];		

	BinaryConstraint* ctr = x->getConstr(y);   		
	if(ctr)	{
        ctr->reconnect();
		ctr->addCosts(x,y,costs);
        ctr->propagate();
	}
    else {
      if (!ToulBar2::vac) {
        ctr = new BinaryConstraint(this, (EnumeratedVariable *) vars[xIndex], (EnumeratedVariable *) vars[yIndex], costs, &storeData->storeCost);
      }
      else {
        ctr = new VACConstraint(this, (EnumeratedVariable *) vars[xIndex], (EnumeratedVariable *) vars[yIndex], costs, &storeData->storeCost);
      }
    }

    return ctr->wcspIndex;
}


int WCSP::postTernaryConstraint(int xIndex, int yIndex, int zIndex, vector<Cost> &costs)
{
    isternary = true;

	EnumeratedVariable* x =  (EnumeratedVariable *) vars[xIndex];
	EnumeratedVariable* y =  (EnumeratedVariable *) vars[yIndex];
	EnumeratedVariable* z =  (EnumeratedVariable *) vars[zIndex];

	TernaryConstraint* ctr = x->getConstr(y,z);    		

	if(!ctr)
	{
		unsigned int a,b;
		vector<Cost> zerocostsxy;
		vector<Cost> zerocostsxz;
		vector<Cost> zerocostsyz;
		
	    for (a = 0; a < x->getDomainInitSize(); a++) { for (b = 0; b < y->getDomainInitSize(); b++) { zerocostsxy.push_back(0); } }
	    for (a = 0; a < x->getDomainInitSize(); a++) { for (b = 0; b < z->getDomainInitSize(); b++) { zerocostsxz.push_back(0); } }
	    for (a = 0; a < y->getDomainInitSize(); a++) { for (b = 0; b < z->getDomainInitSize(); b++) { zerocostsyz.push_back(0); } }

		BinaryConstraint* xy = x->getConstr(y);
		BinaryConstraint* xz = x->getConstr(z);
		BinaryConstraint* yz = y->getConstr(z);
	       
		if(!xy) { xy = new BinaryConstraint(this, x, y, zerocostsxy, &storeData->storeCost); xy->deconnect(); }
		if(!xz) { xz = new BinaryConstraint(this, x, z, zerocostsxz, &storeData->storeCost); xz->deconnect(); }
		if(!yz) { yz = new BinaryConstraint(this, y, z, zerocostsyz, &storeData->storeCost); yz->deconnect(); }

	    ctr = new TernaryConstraint(this,x,y,z, xy, xz, yz, costs, &storeData->storeCost);  
	}
	else ctr->addCosts(x,y,z,costs);
	
	return ctr->wcspIndex;
}


int WCSP::postNaryConstraint(EnumeratedVariable** scopeVars, int arity, Cost defval)
{
     NaryConstraint *ctr = new NaryConstraint(this,scopeVars,arity,defval);  
     return ctr->wcspIndex;	
}

int WCSP::postNaryConstraint(int* scopeIndex, int arity, Cost defval)
{
     EnumeratedVariable** scopeVars = new EnumeratedVariable* [arity];
     for(int i = 0; i < arity; i++) scopeVars[i] = (EnumeratedVariable *) vars[scopeIndex[i]];
     NaryConstraint *ctr = new NaryConstraint(this,scopeVars,arity,defval);  
     delete [] scopeVars;
     return ctr->wcspIndex;
}



void WCSP::postSupxyc(int xIndex, int yIndex, Value cste)
{
    new Supxyc(this, vars[xIndex], vars[yIndex], cste, &storeData->storeCost, &storeData->storeValue);
}

void WCSP::sortConstraints()
{
    for (unsigned int i=0; i<vars.size(); i++) {
        vars[i]->sortConstraints();
    }
}




void WCSP::processTernary()
{
	double maxtight = 0;
	Variable* var;
	TernaryConstraint *tctr1max = NULL, *tctr2max = NULL;

    if (ToulBar2::preprocessTernaryHeuristic) {

        for (unsigned int i=0; i<vars.size(); i++) {
    		TernaryConstraint *tctr1, *tctr2;
            double tight = vars[i]->strongLinkedby(var,tctr1,tctr2);
            if(tight > maxtight) { maxtight = tight; tctr1max = tctr1; tctr2max = tctr2; }
    		if(tctr1 && tctr2 && (tctr1 != tctr2)) {
    			tctr1->extendTernary();
    			tctr2->extendTernary();
    			BinaryConstraint* b = tctr1->commonBinary(tctr2);
    			if(b->connected())
    			{
	    			tctr1->projectTernaryBinary(b);	
	    			tctr2->projectTernaryBinary(b);
	    			b->propagate();	
    			}
    		}
        }
        
        if(ToulBar2::verbose > 0) {
            cout << "Strongest part has mean cost: " << maxtight;
            if(var) cout << "  Variable: " << var->wcspIndex;  
            if(tctr1max) cout << ", 1. ternary with tight: " << tctr1max->getTightness();
            if(tctr2max) cout << ", 2. ternary with tight: " << tctr2max->getTightness();
            cout << endl;
        }
    }
    
    for (unsigned int i=0; i<constrs.size(); i++) 
    	if(constrs[i]->arity() == 3)  
    	{
    		TernaryConstraint* t = (TernaryConstraint*) constrs[i];
    		t->extendTernary();
    		t->projectTernary();
    	}
}

void WCSP::preprocessing()
{
	
	if (ToulBar2::preprocessTernary) {
        cout << "Preproject ternary constraints to binary constraints" << endl;
        processTernary();
    }

	ToulBar2::elimDegree_preprocessing = -ToulBar2::elimDegree_preprocessing;
    for (unsigned int i=0; i<vars.size(); i++) vars[i]->queueEliminate();
    cout << "Variable elimination in preprocessing of real degree <= " << ToulBar2::elimDegree_preprocessing << endl; 

    propagate();		
    
    ToulBar2::elimDegree = -ToulBar2::elimDegree;
    if(ToulBar2::elimDegree >= 0) {
    	initElimConstrs();
        cout << "Variable elimination during search of degree <= " << ToulBar2::elimDegree << endl; 		
    }
     
    if (getenv("TB2DEGREE")) degreeDistribution();
}

Value WCSP::getDomainSizeSum()
{
    Value sum = 0;
    for (unsigned int i=0; i<vars.size(); i++) {
        if (vars[i]->unassigned()) sum += vars[i]->getDomainSize();
    }
    return sum;
}

bool WCSP::getEnumDomain(int varIndex, Value *array)
{
	if (EnumeratedVariable *var = dynamic_cast<EnumeratedVariable*>(vars[varIndex])) {
		var->getDomain(array);
		return true;
	} else return false;
}

bool WCSP::getEnumDomainAndCost(int varIndex, ValueCost *array)
{
	if (EnumeratedVariable *var = dynamic_cast<EnumeratedVariable*>(vars[varIndex])) {
		var->getDomainAndCost(array);
		return true;
	} else return false;
}

unsigned int WCSP::numberOfConnectedConstraints() const
{
    int res = 0; 
    for (unsigned int i=0; i<constrs.size(); i++) if (constrs[i]->connected()) res++;
    for (int i=0; i<elimOrder; i++) if (elimConstrs[i]->connected()) res++;
    return res;
}

void WCSP::degreeDistribution()
{
  int* degDistrib = new int [vars.size()];
 
  for (unsigned int i=0; i<vars.size(); i++) degDistrib[i] = 0;
  for (unsigned int i=0; i<vars.size(); i++) if (unassigned(i)) degDistrib[getRealDegree(i)]++;
 
  unsigned int lastnonzero = 0;
  for (unsigned int i=0; i<vars.size(); i++)
      if(degDistrib[i]) lastnonzero = i;
 
  ofstream file("problem.deg");
  for (unsigned int i=0; i<=lastnonzero; i++) if (degDistrib[i]) file << i << " " << degDistrib[i] << endl;
  delete [] degDistrib;
  
  int limit = atoi(getenv("TB2DEGREE"));
  if (limit > 0) {
      cout << "remove variables with degree > " << limit << endl;
      for (unsigned int i=0; i<vars.size(); i++) if (unassigned(i) && getRealDegree(i) > limit) {
            cout << "remove variable " << getName(i) << endl;
            for (ConstraintList::iterator iter=vars[i]->getConstrs()->begin(); iter != vars[i]->getConstrs()->end(); ++iter) {
                (*iter).constr->deconnect();
            }
            assign(i, getSupport(i));
      }
  }
  
  ofstream pb("problem.wcsp");
  dump(pb);
  exit(0);
}

void WCSP::printNCBuckets()
{
    for (int bucket = 0; bucket < NCBucketSize; bucket++) {
        cout << "NC " << bucket << ":";
        for (VariableList::iterator iter = NCBuckets[bucket].begin (); iter != NCBuckets[bucket].end(); ++iter) {
           cout << " " << (*iter)->getName() << "," << (*iter)->getMaxCostValue() << "," << (*iter)->getMaxCost();
           assert((*iter)->canbe((*iter)->getMaxCostValue()));
           assert((*iter)->getCost((*iter)->getMaxCostValue()) == (*iter)->getMaxCost());
           assert((bucket) ? ((*iter)->getMaxCost() >= (Long) pow(2.,bucket)) : ((*iter)->getMaxCost() > 0));
           //assert((*iter)->getMaxCost() >= (Long) pow(2.,bucket));
           assert((*iter)->getMaxCost() < (Long) pow(2.,bucket+1));
        }
        cout << endl;
    }
}

void WCSP::print(ostream& os)
{
    os << "Objective: [" << getLb() << "," << getUb() << "]" << endl;
    os << "Variables:" << endl;
    for (unsigned int i=0; i<vars.size(); i++) os << *vars[i] << endl;
    if (ToulBar2::verbose >= 4) {
        os << "Constraints:" << endl;
        for (unsigned int i=0; i<constrs.size(); i++) if (constrs[i]->connected()) os << *constrs[i];
        for (int i=0; i<elimOrder; i++) if (elimConstrs[i]->connected()) os << *elimConstrs[i];
    }
}

void printClique(ostream& os, int arity, Constraint *ctr) {
    assert(arity >= 2);
    for (int i=0;i<arity-1;i++) {
        for (int j=i+1;j<arity;j++) {
            os << ctr->getVar(i)->wcspIndex + 1 << " " << ctr->getVar(j)->wcspIndex + 1 << endl;
        }
    }
}

// Warning! make the assumption that all initial domains start at zero!!!
void WCSP::dump(ostream& os)
{
    Value maxdomsize = 0;
    Value xcosts = 0;
    if (getLb() > 0) xcosts++;
    for (unsigned int i=0; i<vars.size(); i++) {
        if (vars[i]->getInf() < 0) {
            cerr << "Cannot save domain of variable " << vars[i]->getName() << " with negative values!!!" << endl;
            exit(EXIT_FAILURE);
        }
        if (vars[i]->getSup()+1 > maxdomsize) maxdomsize = vars[i]->getSup()+1;
        if (vars[i]->enumerated()) xcosts++;
        else if (vars[i]->getInfCost() > 0 || vars[i]->getSupCost() > 0) {
            cerr << "Cannot save interval variable " << vars[i]->getName() << " with bound unary costs!!!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    os << "wcsp " << numberOfVariables() << " " << maxdomsize << " " << numberOfConnectedConstraints()+xcosts << " " << getUb() << endl;
    for (unsigned int i=0; i<vars.size(); i++) {
        if (!vars[i]->enumerated()) os << "-";
        os << vars[i]->getSup()+1;
        if (i < vars.size()-1) os << " ";
    }
    os << endl;
    for (unsigned int i=0; i<constrs.size(); i++) if (constrs[i]->connected()) constrs[i]->dump(os);
    for (int i=0; i<elimOrder; i++) if (elimConstrs[i]->connected()) elimConstrs[i]->dump(os);
    for (unsigned int i=0; i<vars.size(); i++) {
        if (vars[i]->enumerated()) {
            int size = vars[i]->getDomainSize();
            ValueCost domcost[size]; // replace size by MAX_DOMAIN_SIZE in case of compilation problem
            getEnumDomainAndCost(i, domcost);
            os << "1 " << i << " " << getUb() << " " << size << endl;
            for (int v=0; v<size; v++) {
                os << domcost[v].value << " " << domcost[v].cost << endl;
            }
        }
    }
    if (getLb() > 0) os << "0 " << getLb() << " 0" << endl;
    
    if (getenv("TB2GRAPH")) {
        ofstream pb("problem.graph");
        int res = 0; 
        for (unsigned int i=0; i<constrs.size(); i++) if (constrs[i]->connected()) res+=(constrs[i]->arity()*(constrs[i]->arity()-1)/2);
        for (int i=0; i<elimOrder; i++) if (elimConstrs[i]->connected()) res+=(elimConstrs[i]->arity()*(elimConstrs[i]->arity()-1)/2);
        pb << res << " " << numberOfVariables() << endl;
        for (unsigned int i=0; i<constrs.size(); i++) if (constrs[i]->connected()) {
//            pb << constrs[i]->getVar(0)->wcspIndex + 1;
//            for (int j=1; j<constrs[i]->arity(); j++) {
//                pb << " " << constrs[i]->getVar(j)->wcspIndex + 1;
//            }
//            pb << endl;
            printClique(pb, constrs[i]->arity(), constrs[i]);
        }
        for (int i=0; i<elimOrder; i++) if (elimConstrs[i]->connected()) {
//            pb << elimConstrs[i]->getVar(0)->wcspIndex + 1;
//            for (int j=1; j<elimConstrs[i]->arity(); j++) {
//                pb << " " << elimConstrs[i]->getVar(j)->wcspIndex + 1;
//            }
//            pb << endl;
            printClique(pb, elimConstrs[i]->arity(), elimConstrs[i]);
        }
    }
    
    if (ToulBar2::pedigree) ToulBar2::pedigree->save("problem.pre", this, false, true);
}

ostream& operator<<(ostream& os, WCSP &wcsp)
{
    wcsp.print(os);
    return os;
}

ostream& operator<<(ostream& os, WeightedCSP &wcsp)
{
    wcsp.print(os);
    return os;
}


/*
 * WCSP propagation methods
 * 
 */

bool WCSP::verify()
{
    for (unsigned int i=0; i<vars.size(); i++) {
        if (vars[i]->unassigned() && !vars[i]->verifyNC()) return false;
        // warning! in the CSP case, EDAC is no equivalent to GAC on ternary constraints due to the combination with binary constraints
        if (!ToulBar2::FDAC && vars[i]->unassigned() && (!isternary || getUb()-getLb() > 1) && !vars[i]->isEAC()) {
            cout << "variable " << vars[i]->getName() << " not EAC!" << endl;
            return false;
        }
    }
    for (unsigned int i=0; i<constrs.size(); i++) {
        if (constrs[i]->connected() &&  !constrs[i]->verify()) return false;
    }
    return true;
}

void WCSP::whenContradiction()
{
    NC.clear();
    IncDec.clear();
    AC.clear();
    DAC.clear();
    EAC1.clear();
    EAC2.clear();
    Eliminate.clear();
    if (ToulBar2::vac) {
      vac->clear();
    }
    objectiveChanged = false;
    nbNodes++;
}

void WCSP::propagateNC()
{
    if (ToulBar2::verbose >= 2) cout << "NCQueue size: " << NC.getSize() << endl;
    while (!NC.empty()) {
        Variable *x = NC.pop();
        if (x->unassigned()) x->propagateNC();
    }
    if (ToulBar2::verbose >= 3) {
        for (unsigned int i=0; i<vars.size(); i++) cout << *vars[i] << endl;
    }
    if (ToulBar2::verbose >= 2) printNCBuckets();

    if (objectiveChanged) {
        objectiveChanged = false;
        int bucket = cost2log2(getUb() - getLb());
        if (bucket < 0) bucket = 0;
        for (; bucket < NCBucketSize; bucket++) {
            for (VariableList::iterator iter = NCBuckets[bucket].begin(); iter != NCBuckets[bucket].end();) {
                Variable *x = *iter;
                ++iter; // Warning! the iterator could be moved to another place by propagateNC
                if (x->unassigned() && x->getMaxCost() + getLb() >= getUb()) x->propagateNC();
            }
        }
    }
}

void WCSP::propagateIncDec()
{
    if (ToulBar2::verbose >= 2) cout << "IncDecQueue size: " << IncDec.getSize() << endl;
    while (!IncDec.empty()) {
        int incdec;
        Variable *x = IncDec.pop(&incdec);
        if (x->unassigned()) x->propagateIncDec(incdec);
    }
}

void WCSP::propagateAC()
{
    if (ToulBar2::verbose >= 2) cout << "ACQueue size: " << AC.getSize() << endl;
    while (!AC.empty()) {
        EnumeratedVariable *x = (EnumeratedVariable *) ((ToulBar2::FDAComplexity)?AC.pop_min():AC.pop());
        if (x->unassigned()) x->propagateAC();
        // Warning! propagateIncDec() necessary to transform inc/dec event into remove event
        propagateIncDec();          // always examine inc/dec events before remove events
    }
}

void WCSP::propagateDAC()
{
    if (ToulBar2::verbose >= 2) cout << "DACQueue size: " << DAC.getSize() << endl;
    while (!DAC.empty()) {
        EnumeratedVariable *x = (EnumeratedVariable *) ((ToulBar2::FDAComplexity)?DAC.pop_max():DAC.pop());
        if (x->unassigned()) x->propagateDAC();
        propagateIncDec();          // always examine inc/dec events before projectFromZero events
    }
}

void WCSP::fillEAC2()
{
  assert(EAC2.empty());
  if (ToulBar2::verbose >= 2) cout << "EAC1Queue size: " << EAC1.getSize() << endl;
  while(!EAC1.empty()) {
    EnumeratedVariable *x = (EnumeratedVariable *) ((ToulBar2::FDAComplexity)?EAC1.pop_min():EAC1.pop());
    if (x->unassigned()) x->fillEAC2(true);
  }
}

void WCSP::propagateEAC()
{
    fillEAC2();
    if (ToulBar2::verbose >= 2) cout << "EAC2Queue size: " << EAC2.getSize() << endl;
    while (!EAC2.empty()) {
        EnumeratedVariable *x = (EnumeratedVariable *) ((ToulBar2::FDAComplexity)?EAC2.pop_min():EAC2.pop());
        if (x->unassigned()) x->propagateEAC();
        propagateIncDec();          // always examine inc/dec events before projectFromZero events
    }
}



void WCSP::eliminate()
{
	//if(!Eliminate.empty())
	while(!Eliminate.empty())
	{
	    EnumeratedVariable *x = (EnumeratedVariable *) Eliminate.pop_first();
	    if (x->unassigned()) x->eliminate();
	}
}



void WCSP::propagate()
{    
//    revise(NULL);
	
  do {
    do {
      eliminate();
      while (objectiveChanged || !NC.empty() || !IncDec.empty() || !AC.empty() || !DAC.empty()
          || (((ToulBar2::vac) || (!ToulBar2::FDAC && (getUb()-getLb() > 1))) && !EAC1.empty())) {
        propagateIncDec();
        if ((ToulBar2::vac) || (!ToulBar2::FDAC && getUb()-getLb() > 1)) propagateEAC();
        assert(IncDec.empty());
        propagateDAC();
        assert(IncDec.empty());
        propagateAC();
        assert(IncDec.empty());
        propagateNC();
      }
    } while (!Eliminate.empty());
    if (ToulBar2::FDAC || getUb()-getLb() <= 1) EAC1.clear();
    if (ToulBar2::vac) vac->propagate();
  } while ((ToulBar2::vac) && (!vac->remainedIdle()));

  //    revise(NULL);

  assert(verify());
  assert(!objectiveChanged);
  assert(NC.empty());
  assert(IncDec.empty());
  assert(AC.empty());
  assert(DAC.empty());
  assert(EAC1.empty());
  assert(EAC2.empty());
  assert(Eliminate.empty());
  nbNodes++;
}


void WCSP::restoreSolution()
{
	int elimo = getElimOrder();
	for(int i = elimo-1; i>=0; i--)
	{
		elimInfo ei = elimInfos[i];
		EnumeratedVariable* x = (EnumeratedVariable*) ei.x;
		EnumeratedVariable* y = (EnumeratedVariable*) ei.y;
		EnumeratedVariable* z = (EnumeratedVariable*) ei.z;
		BinaryConstraint* xy = ei.xy;
		BinaryConstraint* xz = ei.xz;
		TernaryConstraint* xyz = ei.xyz;
		//   if (xy) cout << "xy=" << *xy << endl;
		//   if (xz) cout << "xz=" << *xz << endl;
		//   if (xyz) cout << "xyz=" << *xyz << endl;
    
		Value vy = -1;
		Value vz = -1;
			
		if(y) vy = getValue(y->wcspIndex);
		if(z) vz = getValue(z->wcspIndex);		
	
		int minv = -1;		    
	    Cost mincost = MAX_COST;
	
	    for (unsigned int vxi = 0; vxi < x->getDomainInitSize(); vxi++)
	    {
	    	Value vx = x->toValue(vxi);
	    	if(!x->canbeAfterElim(vx)) continue; 
	    	Cost cxy = 0;
	    	Cost cxz = 0;
	    	Cost cxyz = 0;
	    	
	    	if(xy) {
	    		if(xy->getIndex(y) == 0) cxy = xy->getCost(vy, vx);
	    		else cxy = xy->getCost(vx, vy);
	    	}

	    	if(xz) {
	    		if(xz->getIndex(z) == 0) cxz = xz->getCost(vz, vx);
	    		else cxz = xz->getCost(vx, vz);
	    	}

			if(xyz) cxyz = xyz->getCost(x,y,z, vx, vy, vz); 

	    	Cost c = x->getCost(vx) + cxy + cxz + cxyz;
			//cout << "test " << vx << "," << x->getCost(vx) << "," << cxy << "," << cxz << "," << cxyz << endl;
	    	    	
	    	if(c < mincost)	{
	    		mincost = c;
	    		minv = vx; 
	    	}
	    }
		//cout << i << ": elim " << x->getName() << "_" << minv << ", y= " << ((y)?y->getName():"-") << "_" << vy << ", z= " << ((z)?z->getName():"-") << "_" << vz << endl;
	    x->assignWhenEliminated( minv );
	}
}



// Creates n fake empty constraints and puts them in the pool 'elimConstrs'
void WCSP::initElimConstrs()
{
	unsigned int i;
    BinaryConstraint* xy;
    for (i=0; i<vars.size(); i++) {	
	   xy = new BinaryConstraint(this, &storeData->storeCost); 
	   elimConstrs.push_back(xy);
	   elimInfo ei = { NULL, NULL, NULL, NULL, NULL, NULL };
	   elimInfos.push_back(ei);
       vars[i]->queueEliminate();
    }
    /*TernaryConstraint* xyz;
    for (i=0; i<vars.size(); i++) {	
	   xyz = new TernaryConstraint(this, &storeData->storeCost); 
	   elimTernaryConstrs.push_back(xyz);
	   elimInfo ei = { NULL, NULL, NULL, NULL, NULL, NULL };
	   elimInfos.push_back(ei);
    }
    NaryConstraint* ctr;
    for (i=0; i<vars.size(); i++) {	
	   ctr = new NaryConstraint(this); 
	   elimNaryConstrs.push_back(ctr);
	   elimInfo ei = { NULL, NULL, NULL, NULL, NULL, NULL };
	   elimInfos.push_back(ei);
    }*/

}


// Function that adds a new binary constraint from the pool of fake constraints
BinaryConstraint* WCSP::newBinaryConstr( EnumeratedVariable* x, EnumeratedVariable* y )
{
	int newIndex = (int) elimOrder;
	BinaryConstraint* ctr = (BinaryConstraint*) elimConstrs[newIndex]; 
	ctr->fillElimConstr(x,y);
	return ctr;
}

// for managing variable elimination of any degree during search. Not used yet
TernaryConstraint* WCSP::newTernaryConstr( EnumeratedVariable* x, EnumeratedVariable* y, EnumeratedVariable* z )
{
	//int newIndex = (int) ternaryOrder;
	TernaryConstraint* ctr = NULL; 
	return ctr;
}

// for managing variable elimination of any degree during search. Not used yet
NaryConstraint* WCSP::newNaryConstr( TSCOPE& scope )
{
	//int newIndex = (int) naryOrder;
	NaryConstraint* ctr = NULL; 
	return ctr;
}

void WCSP::deleteTmpConstraint( Constraint* ctr )
{
}


Constraint* WCSP::sum( Constraint* ctr1, Constraint* ctr2  )
{
	if (ToulBar2::verbose >= 1) cout << endl << "Sum of constraints: " << *ctr1 << " " << *ctr2 << endl;
		
	if(ctr1->order(ctr2) < 0) { Constraint* ctraux = ctr1; ctr1 = ctr2; ctr2 = ctraux; }
		
	ctr1->deconnect();
	ctr2->deconnect();
		
	TSCOPE scopeUinv;
	TSCOPE scopeIinv;
	ctr1->scopeUnion(scopeUinv, ctr2);
	ctr1->scopeCommon(scopeIinv, ctr2);
	int arityU = scopeUinv.size();
	int arityI = scopeIinv.size();
	EnumeratedVariable** scopeU = new EnumeratedVariable* [arityU];
	EnumeratedVariable** scopeI = new EnumeratedVariable* [arityI];


	int i = 0;
	TSCOPE::iterator it = scopeIinv.begin();
	while(it != scopeIinv.end()) {
		int xi = it->first;
		scopeU[i] = (EnumeratedVariable*) vars[xi];
		scopeI[i] = scopeU[i];
		it++; i++;
		scopeUinv.erase(xi);
	}	
    it = scopeUinv.begin();
	while(it != scopeUinv.end()) {
		int xi = it->first;
		scopeU[i] = (EnumeratedVariable*) vars[xi];
		it++; i++;
	}	
			
	EnumeratedVariable* x = scopeU[0]; 
	EnumeratedVariable* y = scopeU[1]; 

	Cost Top = getUb(); 
	unsigned int vx,vy,vz;
	string tuple,tuple1,tuple2;
	Cost cost,cost1,cost2;
	int ctrIndex = -1; 
	Constraint* ctr = NULL;
	vector<Cost> costs;
		
	if(arityU > 3) {
		ctrIndex= postNaryConstraint(scopeU, arityU, Top);
		ctr =  constrs[ctrIndex];
		NaryConstraint* nary = (NaryConstraint*) ctr;
		
		bool tupleXtuple = (ctr1->getDefCost() >= Top) && (ctr2->getDefCost() >= Top);

		if(tupleXtuple) {
			ctr1->first();
		    while(ctr1->next(tuple1,cost1)) {
				ctr2->first();
			    while(ctr2->next(tuple2,cost2)) {
			    	if(cost1 + cost2 < Top) nary->insertSum(tuple1,cost1,ctr1,tuple2,cost2,ctr2);
			    }
			}
		} else {
			bool notfin = true;
			nary->firstlex(tuple); 
		    while(notfin) {
		    	cost1 = ctr1->evalsubstr(tuple, nary );	
		    	cost2 = ctr2->evalsubstr(tuple, nary );	
			    if(cost1 + cost2 < Top) nary->setTuple(tuple, cost1 + cost2);
			    notfin = nary->nextlex(tuple,cost);
			}
		}    
	}
	else if(arityU == 3) {
		EnumeratedVariable* z = scopeU[2]; 
		for (vx = 0; vx < x->getDomainInitSize(); vx++)
		  for (vy = 0; vy < y->getDomainInitSize(); vy++)
		    for (vz = 0; vz < z->getDomainInitSize(); vz++) {
		    	if(arityI == 1)       costs.push_back( ((BinaryConstraint*)ctr1)->getCost(x,y,vx,vy) + ((BinaryConstraint*)ctr2)->getCost(x,z,vx,vz));		
		    	else if(arityI == 2)  costs.push_back( ((BinaryConstraint*)ctr1)->getCost(x,y,vx,vy) + ((TernaryConstraint*)ctr2)->getCost(x,y,z,vx,vy,vz));
		    	else if(arityI == 3)  costs.push_back( ((TernaryConstraint*)ctr1)->getCost(x,y,z,vx,vy,vz) + ((TernaryConstraint*)ctr2)->getCost(x,y,z,vx,vy,vz));
		    	else assert(false);
		    }		   
		ctrIndex = postTernaryConstraint( x->wcspIndex, y->wcspIndex, z->wcspIndex, costs );  
	}
	else if(arityU == 2) {
		BinaryConstraint* bctr1 = (BinaryConstraint*) ctr1;
		BinaryConstraint* bctr2 = (BinaryConstraint*) ctr2;
		for (vx = 0; vx < x->getDomainInitSize(); vx++)
		  for (vy = 0; vy < y->getDomainInitSize(); vy++)
		   	 costs.push_back( bctr1->getCost(x,y,vx,vy) + bctr2->getCost(x,y,vx,vy) );		
		   	 
		ctrIndex= postBinaryConstraint( x->wcspIndex, y->wcspIndex, costs );  
	}
	assert(ctrIndex >= 0);
	delete [] scopeU;
	delete [] scopeI;
	ctr =  constrs[ctrIndex];
	ctr->propagate();	
	
	if (ToulBar2::verbose >= 1) cout << endl << "Has result: " << *ctr << endl;
	
	return ctr;
}
  
  
  
       
void WCSP::project( Constraint* &ctr_inout, EnumeratedVariable* var  )
{
	unsigned int vx,vy,vz,vvar;

	int arity = ctr_inout->arity();
	if(ctr_inout->getIndex(var) < 0) return;	

	if(arity-1 > 3) { 
		((NaryConstraint*)ctr_inout)->project(var); 
		ctr_inout->propagate();
		return; 
	}

	ctr_inout->deconnect();


	int i,j;	
	int ivars[3];
	EnumeratedVariable* evars[3];
	
	j = 0;
	for(i=0; i<arity; i++) {
		EnumeratedVariable* v = (EnumeratedVariable*) ctr_inout->getVar(i);
		if(v != var) { ivars[j] = v->wcspIndex; evars[j] = v; j++; } 
	}
	
	Constraint* ctr = NULL;
	TernaryConstraint* tctr = NULL;
	Cost Top = getUb(); 
	int ctrIndex;
	char t[5];
	vector<Cost> costs;


	
	switch(arity-1)
	{
		case 3:  for (vx = 0; vx < evars[0]->getDomainInitSize(); vx++)
				  for (vy = 0; vy < evars[1]->getDomainInitSize(); vy++)
				   for (vz = 0; vz < evars[2]->getDomainInitSize(); vz++) {  	
				   	 t[ ctr_inout->getIndex(evars[0]) ] = vx + CHAR_FIRST;
				   	 t[ ctr_inout->getIndex(evars[1]) ] = vy + CHAR_FIRST;
				   	 t[ ctr_inout->getIndex(evars[2]) ] = vz + CHAR_FIRST;

				   	 Cost mincost = Top;
				  	 for (vvar = 0; vvar < var->getDomainInitSize(); vvar++) {
					   	t[ctr_inout->getIndex(var)] = vvar + CHAR_FIRST;
			   	   		t[4] =  '\0';
				   	 	Cost c = ((NaryConstraint*)ctr_inout)->eval(string(t)) + var->getCost(vvar);
				   	 	if(c < mincost) mincost = c;
				   	 }  	
				   	 costs.push_back(mincost);	 
				   }
			  	  ctrIndex = postTernaryConstraint(ivars[0],ivars[1],ivars[2],costs);  
				  ctr = constrs[ctrIndex];
				  break;

		case 2:   tctr = (TernaryConstraint*)ctr_inout;
				  for (vx = 0; vx < evars[0]->getDomainInitSize(); vx++)
					for (vy = 0; vy < evars[1]->getDomainInitSize(); vy++) {
					  Cost mincost = Top;
					  for (vvar = 0; vvar < var->getDomainInitSize(); vvar++) {
				   	   	Cost c = ((TernaryConstraint*)ctr_inout)->getCost(evars[0], evars[1], var, vx, vy, vvar)  + var->getCost(vvar);
	 				    if(c < mincost) mincost = c;
					  }  		 
				      costs.push_back(mincost);
				   }				 	  
				   ctrIndex = postBinaryConstraint(ivars[0],ivars[1],costs);  
				   ctr = constrs[ctrIndex];
				break;						
			
		case 1: for (vx = 0; vx < evars[0]->getDomainInitSize(); vx++) {
				   	 Cost mincost = Top;
				  	 for (vvar = 0; vvar < var->getDomainInitSize(); vvar++) {
			   	   		Cost c = ((BinaryConstraint*)ctr_inout)->getCost(evars[0], var, vx, vvar)  + var->getCost(vvar);
 				   	 	if(c < mincost) mincost = c;
				   	 }  		 
			     	 costs.push_back(mincost);
				  }
				for(unsigned int a = 0; a < evars[0]->getDomainInitSize(); a++)
					if(costs[a] > 0) evars[0]->project(a, costs[a]);	

				evars[0]->findSupport();
				break;						
							
		default:;
	}
	ctr_inout = ctr;	
	if(ctr) ctr->propagate();
}


void WCSP::variableElimination( EnumeratedVariable* var )
{
	if(var->wcspIndex == 6) return; 
	if (ToulBar2::verbose >= 1) cout << endl << "Variable General Elimination de " << var->getName() << "    degree: " << var->getDegree() << " real: " << var->getRealDegree() << endl;
	list<Constraint*> ctrs;

	if(var->getDegree() > 0) {
		for (ConstraintList::iterator it=var->getConstrs()->begin(); it != var->getConstrs()->end(); ++it)  {
			Constraint* ctr = (*it).constr;
			ctrs.push_back( ctr );
		}
		list<Constraint*>::iterator it = ctrs.begin();
		Constraint* csum1 = *it;
		if(ctrs.size() == 1) csum1->reconnect();
		Constraint* csum2;
		++it;
		
		for (; it != ctrs.end(); ++it) {
	    	Constraint* ctr = *it;
	    	csum2 = sum(ctr,csum1);
			csum1 = csum2;    	
	   	}
	   	project(csum1, var);
	}
   	assert(var->getDegree() == 0);
   	
   	//elimInfo ei = {this,y,z,(BinaryConstraint*) links[(flag_rev)?1:0].constr, (BinaryConstraint*) links[(flag_rev)?0:1].constr, xyz};
	//elimInfos[wcsp->getElimOrder()] = ei;
	//elimination();

	var->assign(var->getSupport());          // warning! dummy assigned value
}


BinaryConstraint* WCSP::getArbitraryBinaryCtr()
{
	vector<BinaryConstraint*> l;
    for (unsigned int i=0; i<constrs.size(); i++) if (constrs[i]->arity() == 2) l.push_back( (BinaryConstraint*) constrs[i] );
	return l[ rand()%(l.size()) ];
}

TernaryConstraint* WCSP::getArbitraryTernaryCtr()
{
	vector<TernaryConstraint*> l;
    for (unsigned int i=0; i<constrs.size(); i++) if (constrs[i]->arity() == 3) l.push_back( (TernaryConstraint*) constrs[i] );
	return l[ rand()%(l.size()) ];
}

NaryConstraint* WCSP::getArbitraryNaryCtr()
{
	vector<NaryConstraint*> l;
    for (unsigned int i=0; i<constrs.size(); i++) if (constrs[i]->arity() > 3) l.push_back( (NaryConstraint*) constrs[i] );
	return l[ rand()%(l.size()) ];
}




Cost WCSP::Prob2Cost(TProb p) const {
	if (p == 0.0) return getUb();
    TProb res = -Log10(p)*ToulBar2::NormFactor;
    if (res > to_double(MAX_COST)) {
      cerr << "Overflow when converting probability to cost." << endl;
      abort();
    }
	Cost c = (Cost) ((Long) res);
	if(c > getUb()) return getUb();
	return c;
}

TProb WCSP::Cost2LogLike(Cost c) const { return -to_double(c)/ToulBar2::NormFactor; }
TProb WCSP::Cost2Prob(Cost c) const { return Pow((TProb)10., -to_double(c)/ToulBar2::NormFactor); }

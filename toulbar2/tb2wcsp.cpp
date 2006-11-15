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


/*
 * Global variables
 * 
 */
 
double ToulBar2::version  = 0.4;
int ToulBar2::verbose  = 0;
bool ToulBar2::showSolutions  = false;
bool ToulBar2::writeSolution  = false;
bool ToulBar2::elimVarWithSmallDegree_  = false;
bool ToulBar2::only_preprocessing  = false;
bool ToulBar2::preprocessTernaryHeuristic  = false;
bool ToulBar2::FDAComplexity = false;
bool ToulBar2::dichotomicBranching = false;
bool ToulBar2::lds = false;
bool ToulBar2::limited = false;
unsigned int ToulBar2::dichotomicBranchingSize = 10;
#ifdef MENDELSOFT
bool ToulBar2::binaryBranching = true;
bool ToulBar2::elimVarWithSmallDegree  = true;
bool ToulBar2::preprocessTernary  = true;
bool ToulBar2::lastConflict = true;
#else
bool ToulBar2::binaryBranching = false;
bool ToulBar2::elimVarWithSmallDegree  = false;
bool ToulBar2::preprocessTernary  = false;
bool ToulBar2::lastConflict = false;
#endif

externalevent ToulBar2::setvalue = NULL;
externalevent ToulBar2::setmin = NULL;
externalevent ToulBar2::setmax = NULL;
externalevent ToulBar2::removevalue = NULL;
externalcostevent ToulBar2::setminobj = NULL;
Pedigree *ToulBar2::pedigree = NULL;

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
		elimOrder(0, &s->storeValue)
{ 
    objectiveChanged = false;
    nbNodes = 0;
    maxdomainsize = 0;
    isternary = false;
    instance = wcspCounter++;
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
    EnumeratedVariable *x = new EnumeratedVariable(this, n, iinf, isup);
    if(maxdomainsize < isup - iinf + 1) maxdomainsize = isup - iinf + 1;
    return x->wcspIndex;
}

int WCSP::makeEnumeratedVariable(string n, Value *d, int dsize)
{
    EnumeratedVariable *x = new EnumeratedVariable(this, n, d, dsize);
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
void WCSP::postBinaryConstraint(int xIndex, int yIndex, vector<Cost> &costs)
{
	EnumeratedVariable* x =  (EnumeratedVariable *) vars[xIndex];
	EnumeratedVariable* y =  (EnumeratedVariable *) vars[yIndex];		

	BinaryConstraint* ctr = x->getConstr(y);   		
	if(ctr)	{
		ctr->addCosts(x,y,costs);
		if(!ctr->connected()) ctr->reconnect();
	}
    else ctr = new BinaryConstraint(this, (EnumeratedVariable *) vars[xIndex], (EnumeratedVariable *) vars[yIndex], costs, &storeData->storeCost);
}


void WCSP::postTernaryConstraint(int xIndex, int yIndex, int zIndex, vector<Cost> &costs)
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

	    ctr = new TernaryConstraint(this,x,y,z, costs, &storeData->storeCost);  
	    ctr->setBinaries(xy,xz,yz);
	}
	else ctr->addCosts(x,y,z,costs);
}



NaryConstraint* WCSP::postNaryConstraint(EnumeratedVariable** scope, int arity, Cost defval)
{
	 return new NaryConstraint(this,scope,arity,defval);  
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

// Creates n fake empty constraints and puts them in the pool 'elimConstrs'
void WCSP::initElimConstrs()
{
    BinaryConstraint* xy;
    for (unsigned int i=0; i<vars.size(); i++) {	
	   xy = new BinaryConstraint(this, &storeData->storeCost); 
	   elimConstrs.push_back(xy);

	   elimInfo ei = { NULL, NULL, NULL, NULL, NULL, NULL };
	   elimInfos.push_back(ei);

       vars[i]->queueEliminate();
    }
}


// Function that adds a new binary constraint from the pool of fake constraints
BinaryConstraint* WCSP::newBinaryConstr( EnumeratedVariable* x, EnumeratedVariable* y )
{
	int newIndex = (int) elimOrder;
	BinaryConstraint* ctr = elimConstrs[newIndex]; 
	ctr->fillElimConstr(x,y);
	return ctr;
}


void WCSP::processTernary()
{
	double maxtight = 0;
	Variable* var;
	TernaryConstraint *tctr1max = NULL, *tctr2max = NULL;

    if (ToulBar2::preprocessTernaryHeuristic) {
        for (unsigned int i=0; i<vars.size(); i++) {
    //        vars[i]->queueEliminate();
    		TernaryConstraint *tctr1, *tctr2;
            double tight = vars[i]->strongLinkedby(var,tctr1,tctr2);
            if(tight > maxtight) { maxtight = tight; tctr1max = tctr1; tctr2max = tctr2; }
    
    		if(tctr1 && tctr2)
    		{
    			tctr1->extendTernary();
    			tctr2->extendTernary();
    
    			BinaryConstraint* b = tctr1->commonBinary(tctr2);
    			if(!b->connected()) b->reconnect();
    	
    			tctr1->projectTernaryBinary(b);	
    			tctr2->projectTernaryBinary(b);
    			b->propagate();	
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
    		//t->extendTernary();
    		t->projectTernary();
    	}
}

void WCSP::preprocessing()
{
	if (ToulBar2::preprocessTernary) {
        cout << "Preproject ternary constraints to binary constraints" << endl;
        processTernary();
    }

    if (ToulBar2::elimVarWithSmallDegree) {

		initElimConstrs();
        
        cout << "Eliminates variables with small degree"; flush(cout);		
        ToulBar2::elimVarWithSmallDegree_ = true;
        propagate();
		
		if(ToulBar2::only_preprocessing) { ToulBar2::elimVarWithSmallDegree_ = false; cout << "  only in preprocessing"; }
		else cout << "  during search";
		cout << endl;
		flush(cout);		
    }
    else propagate();
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

void WCSP::printNCBuckets()
{
    for (int bucket = 0; bucket < NCBucketSize; bucket++) {
        cout << "NC " << bucket << ":";
        for (VariableList::iterator iter = NCBuckets[bucket].begin (); iter != NCBuckets[bucket].end(); ++iter) {
           cout << " " << (*iter)->getName() << "," << (*iter)->getMaxCostValue() << "," << (*iter)->getMaxCost();
           assert((*iter)->canbe((*iter)->getMaxCostValue()));
           assert((*iter)->getCost((*iter)->getMaxCostValue()) == (*iter)->getMaxCost());
           assert((*iter)->getMaxCost() >= (Long) pow(2.,bucket));
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
    }
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
        if (vars[i]->unassigned() && (!isternary || getUb()-getLb() > 1) && !vars[i]->isEAC()) {
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
    if (x->unassigned()) x->fillEAC2(true); // test unassigned ???????????????????????????????????????????????????
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


void WCSP::propagate()
{    
//    revise(NULL);
	
  do {
    eliminate();
    while (objectiveChanged || !NC.empty() || !IncDec.empty() || !AC.empty() || !DAC.empty()
	       || ((getUb()-getLb() > 1) && !EAC1.empty())) {
      propagateIncDec();
      if (getUb()-getLb() > 1) propagateEAC();
      assert(IncDec.empty());
      propagateDAC();
      assert(IncDec.empty());
      propagateAC();
      assert(IncDec.empty());
      propagateNC();
    }
  } while (!Eliminate.empty());
  if (getUb()-getLb() <= 1) EAC1.clear();
  
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

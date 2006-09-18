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
#include "tb2arithmetic.hpp"


/*
 * Global variables
 * 
 */
 
int ToulBar2::verbose  = 0;
bool ToulBar2::showSolutions  = false;
bool ToulBar2::binaryBranching = false;
int ToulBar2::elimLevel  = -1;
bool ToulBar2::only_preprocessing  = false;
externalevent ToulBar2::setvalue = NULL;
externalevent ToulBar2::setmin = NULL;
externalevent ToulBar2::setmax = NULL;
externalevent ToulBar2::removevalue = NULL;
externalevent ToulBar2::setminobj = NULL;

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


// Inefficient function that looks in the list of all constraints
BinaryConstraint* WCSP::existBinaryConstraint( int xIndex, int yIndex )
{
	Constraint* ctr = NULL;    		
    for (unsigned int i=0; i<constrs.size(); i++) {
    	if(constrs[i]->arity() == 2) {
    		ctr = constrs[i];    		
			
    		Variable* x = ctr->getVar(0);
    		Variable* y = ctr->getVar(1);
			bool exist = ((x->wcspIndex == xIndex) && (y->wcspIndex == yIndex)) ||
						 ((x->wcspIndex == yIndex) && (y->wcspIndex == xIndex));

			if(exist) return (BinaryConstraint*) ctr;			 
    	}
    }
    return NULL;
}

TernaryConstraint* WCSP::existTernaryConstraint( int xIndex, int yIndex, int zIndex )
{
	TernaryConstraint* ctr = NULL;    		
    for (unsigned int i=0; i<constrs.size(); i++) {
    	if(constrs[i]->arity() == 3) {
    		ctr = (TernaryConstraint*) constrs[i];    		
			
    		Variable* x = (EnumeratedVariable *) vars[xIndex];
    		Variable* y = (EnumeratedVariable *) vars[yIndex];
    		Variable* z = (EnumeratedVariable *) vars[zIndex];
	
			int neq = 0;
			if(ctr->getIndex(x) >= 0) neq++;
			if(ctr->getIndex(y) >= 0) neq++;
			if(ctr->getIndex(z) >= 0) neq++;
			
			if(neq == 3) return ctr;			 
    	}
    }
    return NULL;
}

// postBinaryConstraint looks for an existing constraint in the whole
// list of constraints (even not connected constraints). It also alocates
// memory for a new constraints. For this two reasons it should 
// ONLY be called when creating the initial wcsp
void WCSP::postBinaryConstraint(int xIndex, int yIndex, vector<Cost> &costs)
{
	BinaryConstraint* ctr = existBinaryConstraint( xIndex, yIndex );    		
	if(ctr)	{
		EnumeratedVariable* x =  (EnumeratedVariable *) vars[xIndex];
		EnumeratedVariable* y =  (EnumeratedVariable *) vars[yIndex];		
		ctr->addCosts(x,y,costs);
		if(!ctr->connected()) ctr->reconnect();
	}
    else ctr = new BinaryConstraint(this, (EnumeratedVariable *) vars[xIndex], (EnumeratedVariable *) vars[yIndex], costs, &storeData->storeCost);
}


void WCSP::postTernaryConstraint(int xIndex, int yIndex, int zIndex, vector<Cost> &costs)
{
    isternary = true;

	TernaryConstraint* ctr = existTernaryConstraint( xIndex, yIndex, zIndex );    		

	EnumeratedVariable* x =  (EnumeratedVariable *) vars[xIndex];
	EnumeratedVariable* y =  (EnumeratedVariable *) vars[yIndex];
	EnumeratedVariable* z =  (EnumeratedVariable *) vars[zIndex];

	if(!ctr)
	{
		unsigned int a,b;
		vector<Cost> zerocostsxy;
		vector<Cost> zerocostsxz;
		vector<Cost> zerocostsyz;
		
	    for (a = 0; a < x->getDomainInitSize(); a++) { for (b = 0; b < y->getDomainInitSize(); b++) { zerocostsxy.push_back(0); } }
	    for (a = 0; a < x->getDomainInitSize(); a++) { for (b = 0; b < z->getDomainInitSize(); b++) { zerocostsxz.push_back(0); } }
	    for (a = 0; a < y->getDomainInitSize(); a++) { for (b = 0; b < z->getDomainInitSize(); b++) { zerocostsyz.push_back(0); } }

		BinaryConstraint* xy = existBinaryConstraint( xIndex, yIndex );
		BinaryConstraint* xz = existBinaryConstraint( xIndex, zIndex );
		BinaryConstraint* yz = existBinaryConstraint( yIndex, zIndex );
	       
		if(!xy) { xy = new BinaryConstraint(this, x, y, zerocostsxy, &storeData->storeCost); xy->deconnect(); }
		if(!xz) { xz = new BinaryConstraint(this, x, z, zerocostsxz, &storeData->storeCost); xz->deconnect(); }
		if(!yz) { yz = new BinaryConstraint(this, y, z, zerocostsyz, &storeData->storeCost); yz->deconnect(); }

	    ctr = new TernaryConstraint(this,x,y,z, costs, &storeData->storeCost);  
	    ctr->setBinaries(xy,xz,yz);
	}
	else ctr->addCosts(x,y,z,costs);
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



void WCSP::preprocessing()
{
    if (ToulBar2::elimLevel >= 0) {

		initElimConstrs();
		
	    for (unsigned int i=0; i<vars.size(); i++) {
	        vars[i]->queueEliminate();
	    }

        cout << "Eliminates variables with small degree <= " << ToulBar2::elimLevel;  flush(cout);		
        propagate();
		
		if(ToulBar2::only_preprocessing) { ToulBar2::elimLevel = -1; cout << "  only in preprocessing"; }
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
        EnumeratedVariable *x = (EnumeratedVariable *) AC.pop_min();
        if (x->unassigned()) x->propagateAC();
        // Warning! propagateIncDec() necessary to transform inc/dec event into remove event
        propagateIncDec();          // always examine inc/dec events before remove events
    }
}

void WCSP::propagateDAC()
{
    if (ToulBar2::verbose >= 2) cout << "DACQueue size: " << DAC.getSize() << endl;
    while (!DAC.empty()) {
        EnumeratedVariable *x = (EnumeratedVariable *) DAC.pop_max();
        if (x->unassigned()) x->propagateDAC();
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
	
    while (!Eliminate.empty() || !IncDec.empty() || !AC.empty() || !DAC.empty() || !NC.empty() || objectiveChanged) 
    {
	    eliminate();  

	    while (!IncDec.empty() || !AC.empty() || !DAC.empty() || !NC.empty() || objectiveChanged) {
	        propagateIncDec();
	        propagateAC();
	        assert(IncDec.empty());
	        propagateDAC();
	        assert(IncDec.empty());
	        propagateNC();
	    }
	}

	
    assert(verify());
    assert(!objectiveChanged);
    assert(NC.empty());
    assert(IncDec.empty());
    assert(AC.empty());
    assert(DAC.empty());
    assert(Eliminate.empty());
    nbNodes++;
}

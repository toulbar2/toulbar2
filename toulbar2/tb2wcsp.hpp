/** \file tb2wcsp.hpp
 *  \brief Global soft constraint representing a weighted CSP
 * 
 */
 
#ifndef TB2WCSP_HPP_
#define TB2WCSP_HPP_

#include "toulbar2.hpp"
#include "tb2variable.hpp"
#include "tb2constraint.hpp"

class BinaryConstraint;
class TernaryConstraint;
class NaryConstraint;
class VACExtension;
class TreeDecomposition;

class WCSP : public WeightedCSP {
    static int wcspCounter; // count the number of instantiations of WCSP
    int instance; // instantiation occurence number
    Store *storeData;
    StoreCost lb;
    Cost ub;
    vector<Variable *> vars;
    vector<Constraint *> constrs;
    int NCBucketSize;
    vector< VariableList > NCBuckets;         // vector of backtrackable lists
    Queue NC;                                 // non backtrackable list
    Queue IncDec;                             // non backtrackable list
    Queue AC;                                 // non backtrackable list
    Queue DAC;                                // non backtrackable list
    Queue EAC1;                               // non backtrackable list
    Queue EAC2;                               // non backtrackable list
    Queue Eliminate;                          // non backtrackable list
    bool objectiveChanged;
    Long nbNodes;                             // used as a time-stamp by Queue methods
    Constraint *lastConflictConstr;

	// make it private because we don't want copy nor assignment
    WCSP(const WCSP &wcsp);
    WCSP& operator=(const WCSP &wcsp);
    
    
public:

    WCSP(Store *s, Cost upperBound);
    
    virtual ~WCSP();
    

    // -----------------------------------------------------------
    // General API for weighted CSP global constraint
    
    int getIndex() const {return instance;} // instantiation occurence number of current WCSP object
    
    Cost getLb() const {return lb;}
    Cost getUb() const {return ub;}

    void updateUb(Cost newUb) {
        if (newUb < ub) {
            ub = newUb;
    		if (vars.size()==0) NCBucketSize = cost2log2gub(ub) + 1;
        }
    }
    
	void enforceUb() {
       if (CUT(((lb % ToulBar2::costMultiplier) != MIN_COST)?(lb + ToulBar2::costMultiplier):lb, ub)) THROWCONTRADICTION;
       objectiveChanged=true;
    }
    void increaseLb(Cost newLb) {
       if (newLb > lb) {
           if (CUT(((newLb % ToulBar2::costMultiplier) != MIN_COST)?(newLb + ToulBar2::costMultiplier):newLb, ub)) THROWCONTRADICTION;
           lb = newLb;
           objectiveChanged=true;
           if (ToulBar2::setminobj) (*ToulBar2::setminobj)(getIndex(), -1, newLb);
       }
    }
    void decreaseUb(Cost newUb) {
       if (newUb < ub) {
           if (CUT(((lb % ToulBar2::costMultiplier) != MIN_COST)?(lb + ToulBar2::costMultiplier):lb, newUb)) THROWCONTRADICTION;
           ub = newUb;
           objectiveChanged=true;
       }
    }     
    
    void setUb(Cost newUb) { ub = newUb; }
    void setLb(Cost newLb) { lb = newLb; }
 
    bool enumerated(int varIndex) const {return vars[varIndex]->enumerated();}
    
    string getName(int varIndex) const {return vars[varIndex]->getName();}
    Value getInf(int varIndex) const {return vars[varIndex]->getInf();}
    Value getSup(int varIndex) const {return vars[varIndex]->getSup();}
    Value getValue(int varIndex) const {return vars[varIndex]->getValue();}
    unsigned int getDomainSize(int varIndex) const {return vars[varIndex]->getDomainSize();}
    bool getEnumDomain(int varIndex, Value *array);
    bool getEnumDomainAndCost(int varIndex, ValueCost *array);

    bool assigned(int varIndex) const {return vars[varIndex]->assigned();}
    bool unassigned(int varIndex) const {return vars[varIndex]->unassigned();}
    bool canbe(int varIndex, Value v) const {return vars[varIndex]->canbe(v);}
    bool cannotbe(int varIndex, Value v) const {return vars[varIndex]->cannotbe(v);}
       
    void increase(int varIndex, Value newInf) {vars[varIndex]->increase(newInf);}
    void decrease(int varIndex, Value newSup) {vars[varIndex]->decrease(newSup);}
    void assign(int varIndex, Value newValue) {vars[varIndex]->assign(newValue);}
    void remove(int varIndex, Value remValue) {vars[varIndex]->remove(remValue);}
        
    Cost getUnaryCost(int varIndex, Value v) const {return vars[varIndex]->getCost(v);}
    Cost getMaxUnaryCost(int varIndex) const {return vars[varIndex]->getMaxCost();}
    Value getMaxUnaryCostValue(int varIndex) const {return vars[varIndex]->getMaxCostValue();}
    Value getSupport(int varIndex) const {return vars[varIndex]->getSupport();}
    
    int getDegree(int varIndex) const {return vars[varIndex]->getDegree();}
    int getTrueDegree(int varIndex) const {return vars[varIndex]->getTrueDegree();}
    double getWeightedDegree(int varIndex) const {return vars[varIndex]->getWeightedDegree();}
    void revise(Constraint *c) {lastConflictConstr = c;}
    void conflict() {
        if (lastConflictConstr) {
            if (ToulBar2::verbose>=2) cout << "Last conflict on " << *lastConflictConstr << endl;
            lastConflictConstr->incConflictWeight();
            lastConflictConstr=NULL;
        }
    }
    
    void whenContradiction();       // after a contradiction, reset propagation queues and increase nbNodes
    void propagate();               // propagate until a fix point and increase nbNodes
    bool verify();

    unsigned int numberOfVariables() const {return vars.size();}
    unsigned int numberOfUnassignedVariables() const {
        int res = 0; 
        for (unsigned int i=0; i<vars.size(); i++) if (unassigned(i)) res++;
        return res;}
    unsigned int numberOfConstraints() const {return constrs.size();}
    unsigned int numberOfConnectedConstraints() const;
    Value getDomainSizeSum();       // total current number of values
    void degreeDistribution();
#ifdef BOOST
    int diameter();
    int connectedComponents();
    int biConnectedComponents();
    void minimumDegreeOrdering();
#endif
    int makeEnumeratedVariable(string n, Value iinf, Value isup);
    int makeEnumeratedVariable(string n, Value *d, int dsize);
    int makeIntervalVariable(string n, Value iinf, Value isup);

    void processTernary();
    
    void postUnary(int xIndex, Value *d, int dsize, Cost penalty);
    void postSupxyc(int xIndex, int yIndex, Value cst, Value deltamax = MAX_VAL-MIN_VAL);
    void postDisjunction(int xIndex, int yIndex, Value cstx, Value csty, Cost penalty);
    void postSpecialDisjunction(int xIndex, int yIndex, Value cstx, Value csty, Value xinfty, Value yinfty, Cost costx, Cost costy);
	int postBinaryConstraint(int xIndex, int yIndex, vector<Cost> &costs);
    int postTernaryConstraint(int xIndex, int yIndex, int zIndex, vector<Cost> &costs);
    int postNaryConstraint(int* scopeIndex, int arity, Cost defval);
    
    void read_wcsp(const char *fileName);
    void read_random(int n, int m, vector<int>& p, int seed, bool forceSubModular = false );

    void print(ostream& os);
    void dump(ostream& os);
    friend ostream& operator<<(ostream& os, WCSP &wcsp);


    // -----------------------------------------------------------
    // Specific API for Variable and Constraint classes

    Store *getStore() {return storeData;}
    Variable   *getVar(int varIndex) const {return vars[varIndex];}
    Constraint *getCtr(int ctrIndex) const {return constrs[ctrIndex];}

    void link(Variable *x) {vars.push_back(x);}
    void link(Constraint *c) {constrs.push_back(c);}

	VariableList* getNCBucket( int ibucket ) { return &NCBuckets[ibucket]; }
    int getNCBucketSize() const {return NCBucketSize;}
    void changeNCBucket(int oldBucket, int newBucket, DLink<Variable *> *elt) {
        if (oldBucket >= 0) NCBuckets[oldBucket].erase(elt, true);
        if (newBucket >= 0) NCBuckets[newBucket].push_back(elt, true);
    }
    void printNCBuckets();

    Long getNbNodes() {return nbNodes;}
    void queueNC(DLink<VariableWithTimeStamp> *link) {NC.push(link, nbNodes);}
    void queueInc(DLink<VariableWithTimeStamp> *link) {IncDec.push(link, INCREASE_EVENT, nbNodes);}
    void queueDec(DLink<VariableWithTimeStamp> *link) {IncDec.push(link, DECREASE_EVENT, nbNodes);}
    void queueAC(DLink<VariableWithTimeStamp> *link) {AC.push(link, nbNodes);}
    void queueDAC(DLink<VariableWithTimeStamp> *link) {DAC.push(link, nbNodes);}
    void queueEAC1(DLink<VariableWithTimeStamp> *link) {EAC1.push(link, nbNodes);}
    void queueEAC2(DLink<VariableWithTimeStamp> *link) {EAC2.push(link, nbNodes);}
    void queueEliminate(DLink<VariableWithTimeStamp> *link) { Eliminate.push(link, nbNodes);  }  
    
    void propagateNC();
    void propagateIncDec();
    void propagateAC();
    void propagateDAC();
    void fillEAC2();
    Queue *getQueueEAC1() {return &EAC1;}
    void propagateEAC();

    void sortVariables();
    void sortConstraints();
    void preprocessing();


    // -----------------------------------------------------------
    // Data and methods for Variable Elimination

	typedef struct {
		EnumeratedVariable* x;
		EnumeratedVariable* y;
		EnumeratedVariable* z;
		BinaryConstraint*   xy;
		BinaryConstraint*   xz;
		TernaryConstraint*  xyz;
	} elimInfo;

	int maxdomainsize;	                              						   

	void initElimConstrs();

	StoreInt elimOrder;    	 				        // used to count the order in which variables are eliminated
	vector<elimInfo> elimInfos; 
	vector<Constraint*>  elimConstrs;

	StoreInt elimBinOrder;    	 				        
	StoreInt elimTernOrder;    	 				        
	vector<Constraint*>  elimBinConstrs;
	vector<Constraint*>  elimTernConstrs;
	
	int getElimOrder() { return (int) elimOrder; } 
	int getElimBinOrder() { return (int) elimBinOrder; } 
	int getElimTernOrder() { return (int) elimTernOrder; } 
	void elimOrderInc() { elimOrder = elimOrder + 1; }   
	void elimBinOrderInc() { elimBinOrder = elimBinOrder + 1; }   
	void elimTernOrderInc() { elimTernOrder = elimTernOrder + 1; }   
    Constraint *getElimCtr(int elimIndex) const {return elimConstrs[elimIndex];}
    Constraint *getElimBinCtr(int elimBinIndex) const {return elimBinConstrs[elimBinIndex];}
    Constraint *getElimTernCtr(int elimTernIndex) const {return elimTernConstrs[elimTernIndex];}
	
	BinaryConstraint*  newBinaryConstr( EnumeratedVariable* x, EnumeratedVariable* y );
	TernaryConstraint*  newTernaryConstr( EnumeratedVariable* x, EnumeratedVariable* y, EnumeratedVariable* z );

	void eliminate();
	void restoreSolution(); 
    
	Constraint* sum( Constraint* ctr1, Constraint* ctr2  );
	void project( Constraint* &ctr_inout, EnumeratedVariable* var  );
	void variableElimination( EnumeratedVariable* var );

	void ternaryCompletion();
	bool  kconsistency(int xIndex, int yIndex, int zIndex, BinaryConstraint* xy, BinaryConstraint* yz, BinaryConstraint* xz );

    // -----------------------------------------------------------
    // Data and methods for Virtual Arc Consistency

    VACExtension*  vac;

    void histogram( Cost c );
    void histogram();
    void iniSingleton();
    void updateSingleton();
    void removeSingleton();
    void printVACStat();
	int  getVACHeuristic();

    // -----------------------------------------------------------
    // Data and methods for Tree Decomposition

    TreeDecomposition* td;
	void buildTreeDecomposition();

    // -----------------------------------------------------------
    // Functions for dealing with probabilities
    // Warning: ToulBar2::NormFactor has to be initialized

    Cost Prob2Cost(TProb p) const;
    TProb Cost2Prob(Cost c) const;
    TProb Cost2LogLike(Cost c) const;    
};

#endif /*TB2WCSP_HPP_*/

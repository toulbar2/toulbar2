/** \file toulbar2.hpp
 *  \brief Main protocol class of a global soft constraint representing a weighted CSP
 * 
 */
 
#ifndef TOULBAR2_HPP_
#define TOULBAR2_HPP_

#include "tb2types.hpp"

class WeightedCSP {
public:
    static WeightedCSP *makeWeightedCSP(Store *s, Cost upperBound);      // WCSP factory
    
    virtual ~WeightedCSP() {}
    
    virtual int getIndex() const = 0;       // instantiation occurence number of current WCSP object
    
    virtual Cost getLb() const = 0;
    virtual Cost getUb() const =0;

    virtual Cost Prob2Cost(TProb p)const =0;
    virtual TProb Cost2LogLike(Cost c) const =0;
    virtual TProb Cost2Prob(Cost c) const =0;

    virtual void updateUb(Cost newUb) =0;
    virtual void enforceUb() =0;
    virtual void increaseLb(Cost newLb) =0;
    virtual void decreaseUb(Cost newUb) =0;
    virtual void setLb(Cost newLb) =0;
    virtual void setUb(Cost newUb) =0;

    virtual bool enumerated(int varIndex) const =0;
    
    virtual string getName(int varIndex) const =0;
    virtual Value getInf(int varIndex) const =0;
    virtual Value getSup(int varIndex) const =0;
    virtual Value getValue(int varIndex) const =0;
    virtual unsigned int getDomainSize(int varIndex) const =0;
    virtual bool getEnumDomain(int varIndex, Value *array) =0;
    virtual bool getEnumDomainAndCost(int varIndex, ValueCost *array) =0;

    virtual bool assigned(int varIndex) const =0;
    virtual bool unassigned(int varIndex) const =0;
    virtual bool canbe(int varIndex, Value v) const =0;
    virtual bool cannotbe(int varIndex, Value v) const =0;
       
    virtual void increase(int varIndex, Value newInf) =0;
    virtual void decrease(int varIndex, Value newSup) =0;
    virtual void assign(int varIndex, Value newValue) =0;
    virtual void remove(int varIndex, Value remValue) =0;
        
    virtual Cost getUnaryCost(int varIndex, Value v) const =0;
    virtual Cost getMaxUnaryCost(int varIndex) const =0;
    virtual Value getMaxUnaryCostValue(int varIndex) const =0;
    virtual Value getSupport(int varIndex) const =0;
    
    virtual int getDegree(int varIndex) const =0;
    virtual double getWeightedDegree(int varIndex) const =0;
    virtual void preprocessing() =0;

    virtual void whenContradiction() =0;       // after a contradiction, reset propagation queues
    virtual void propagate() =0;               // propagate until a fix point
    virtual bool verify() =0;

    virtual unsigned int numberOfVariables() const =0;
    virtual unsigned int numberOfUnassignedVariables() const =0;
    virtual unsigned int numberOfConstraints() const =0;
    virtual unsigned int numberOfConnectedConstraints() const =0;
    virtual Value getDomainSizeSum() =0;       // total current number of values

    virtual int makeEnumeratedVariable(string n, Value iinf, Value isup) =0;
    virtual int makeEnumeratedVariable(string n, Value *d, int dsize) =0;
    virtual int makeIntervalVariable(string n, Value iinf, Value isup) =0;
    
    virtual int postBinaryConstraint(int xIndex, int yIndex, vector<Cost> &costs) =0;
    virtual int postTernaryConstraint(int xIndex, int yIndex, int zIndex, vector<Cost> &costs) =0;
    virtual int postNaryConstraint(int* scope, int arity, Cost defval) =0;
    virtual void postUnary(int xIndex, Value *d, int dsize, Cost penalty) =0;
    virtual void postSupxyc(int xIndex, int yIndex, Value cst, Value deltamax = MAX_VAL-MIN_VAL) =0;
    virtual void postDisjunction(int xIndex, int yIndex, Value cstx, Value csty, Cost penalty) =0;
    virtual void postSpecialDisjunction(int xIndex, int yIndex, Value cstx, Value csty, Value xinfty, Value yinfty, Cost costx, Cost costy) =0;
    
    virtual void read_wcsp(const char *fileName) =0;
    virtual void read_random(int n, int m, vector<int>& p, int seed, bool forceSubModular = false) =0;
    
    virtual int getElimOrder() =0;
    virtual void restoreSolution( Cluster* c = NULL ) =0;

	virtual void buildTreeDecomposition() = 0;
    virtual TreeDecomposition* getTreeDec() = 0;

	virtual void iniSingleton() = 0;
	virtual	void updateSingleton() = 0;
	virtual void removeSingleton() = 0;
	virtual int  getVACHeuristic() = 0;
	virtual void printVACStat() = 0;    
    virtual void print(ostream& os) =0;
    virtual void dump(ostream& os) =0;
};

ostream& operator<<(ostream& os, WeightedCSP &wcsp);

#endif /*TOULBAR2_HPP_*/

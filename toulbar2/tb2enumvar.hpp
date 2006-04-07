/** \file tb2enumvar.hpp
 *  \brief Variable with domain represented by an enumerated domain.
 * 
 */

#ifndef TB2ENUMVAR_HPP_
#define TB2ENUMVAR_HPP_

#include "tb2variable.hpp"
#include "tb2domain.hpp"

class EnumeratedVariable : public Variable
{
    Domain domain;
    vector<StoreCost> costs;
    StoreValue support;     // Warning! the unary support has to be backtrackable 

    DLink<VariableWithTimeStamp> linkACQueue;
    DLink<VariableWithTimeStamp> linkDACQueue;

    void init();
        
    void increaseFast(Value newInf);        // Do not check for a support nor insert in NC and DAC queue
    void decreaseFast(Value newSup);        // Do not check for a support nor insert in NC and DAC queue
    void removeFast(Value val);             // Do not check for a support nor insert in NC and DAC queue
    
public:    
    EnumeratedVariable(WCSP *wcsp, string n, Value iinf, Value isup);
    EnumeratedVariable(WCSP *wcsp, string n, Value *d, int dsize);
    
    bool enumerated() const {return true;}

    unsigned int getDomainInitSize() const {return domain.getInitSize();}
#ifdef WCSPFORMATONLY
    int toIndex(int v) const {return v;}
    int toValue(int idx) const {return idx;}
#else
    unsigned int toIndex(Value v) const {return domain.toIndex(v);}
    Value toValue(int idx) const {return domain.toValue(idx);}
#endif
    unsigned int getDomainSize() const {
        if (assigned()) return 1; 
        else return domain.getSize();
    }

    bool canbe(Value v) const {return v >= inf && v <= sup && domain.canbe(v);}
    bool cannotbe(Value v) const {return v < inf || v > sup || domain.cannotbe(v);}
    
    void increase(Value newInf);
    void decrease(Value newSup);
    void remove(Value value);
    void assign(Value newValue);

    void project(Value value, Cost cost);
    void extend(Value value, Cost cost);
    Value getSupport() const {return support;}
    void setSupport(Value val) {support = val;}    
    Cost getCost(const Value value) const {
        return costs[toIndex(value)] - deltaCost;
    }

    Cost getInfCost() const {return costs[toIndex(getInf())] - deltaCost;}
    Cost getSupCost() const {return costs[toIndex(getSup())] - deltaCost;}
    void projectInfCost(Cost cost);
    void projectSupCost(Cost cost);

    void propagateNC();    
    bool verifyNC();
    void queueAC();                     // public method used also by tb2binconstr.hpp
    void queueDAC();
    void propagateAC();
    void propagateDAC();
    void findSupport();

    class iterator;
    friend class iterator;
    class iterator { //: public Variable::iterator {
        EnumeratedVariable &var;
        Domain::iterator diter;
    public:
        iterator(EnumeratedVariable &v, Domain::iterator iter) : var(v), diter(iter) {}

        Value operator*() const {return *diter;}
        
        iterator &operator++() {    // Prefix form
            if (var.unassigned()) ++diter;
            else {
                if (*diter < var.getValue()) diter = var.domain.lower_bound(var.getValue());
                else diter = var.domain.end();
            }
            return *this;
        }
        
        iterator &operator--() {    // Prefix form
            if (var.unassigned()) --diter;
            else {
                if (*diter > var.getValue()) diter = var.domain.lower_bound(var.getValue());
                else diter = var.domain.end();
            }
            return *this;
        }

        // To see if you're at the end:
        bool operator==(const iterator &iter) const {return diter == iter.diter;}
        bool operator!=(const iterator &iter) const {return diter != iter.diter;}
    };
    iterator begin() {
        if (assigned()) return iterator(*this, domain.lower_bound(getValue()));
        else return iterator(*this, domain.begin());
    }
    iterator end() {return iterator(*this, domain.end());}
    iterator rbegin() {
        if (assigned()) return iterator(*this, domain.upper_bound(getValue()));
        else return iterator(*this, domain.rbegin());
    }
    iterator rend() {return end();}

    //Finds the first available element whose value is greater or equal to v
    iterator lower_bound(Value v) {
        if (assigned()) {
            if (v <= getValue()) return iterator(*this, domain.lower_bound(getValue()));
            else return end();
        } else return iterator(*this, domain.lower_bound(v));
    }

    //Finds the first available element whose value is lower or equal to v
    iterator upper_bound(Value v) {
        if (assigned()) {
            if (v >= getValue()) return iterator(*this, domain.upper_bound(getValue()));
            else return end();
        } else return iterator(*this, domain.upper_bound(v));
    }
    
    void print(ostream& os);
};

#endif /*TB2ENUMVAR_HPP_*/

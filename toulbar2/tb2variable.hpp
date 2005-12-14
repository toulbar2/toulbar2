/** \file tb2variable.hpp
 *  \brief Variable with domain represented by an interval or an enumerated domain.
 * 
 */
 
#ifndef TB2VARIABLE_HPP_
#define TB2VARIABLE_HPP_

#include "tb2domain.hpp"

/*
 * Backtrack exception
 * 
 */

class Contradiction
{
public:
    Contradiction() {if (ToulBar2::verbose >= 2) cout << "... contradiction!" << endl;}
};

/*
 * Main class
 * 
 */

typedef enum {OBJ_VAR=0, AUX_VAR=1, DECISION_VAR=2} VariableType;

class Variable
{
    bool enumerated;            // should be a constant
    VariableType type;
    string name;
    StoreValue inf;
    StoreValue sup;
    StoreValue value;
    Domain domain;
    vector< WCSPLink > wcsps;   // in which soft global constraint the variable occurs
    Solver *solver;
    DLink< Variable * > linkUnassignedVars;
    
    void addVarToSolver();   
    
public:    
    Variable(string n, Value iinf, Value isup, Solver *s, VariableType t);
    Variable(string n, Value iinf, Value isup, Solver *s, VariableType t, bool enumerate);
    Variable(string n, Value *d, int dsize, Solver *s, VariableType t, bool enumerate);

    bool getEnumerated() const {return enumerated;}
    // Warning! Only valid if the variable is represented by an enumerated domain
    unsigned int getDomainInitSize() const {assert(enumerated); return domain.getInitSize();}
    unsigned int toIndex(Value v) const {assert(enumerated); return domain.toIndex(v);}
    Value toValue(int idx) const {assert(enumerated); return domain.toValue(idx);}
   
    string getName() const {return name;}
    Value getInf() const {return inf;}
    Value getSup() const {return sup;}
    Value getValue() const {return value;}
    unsigned int getDomainSize() const {
        if (enumerated) {if (assigned()) return 1; else return domain.getSize();}
        else return sup - inf + 1;
    }

    bool assigned() const {return inf == sup;}
    bool unassigned() const {return inf != sup;}
    bool canbe(Value v) const {return v >= inf && v <= sup && (!enumerated || domain.canbe(v));}
    bool cannotbe(Value v) const {return v < inf || v > sup || (enumerated && domain.cannotbe(v));}
    
    void increase(WCSP *wcsp, Value newInf);
    void decrease(WCSP *wcsp, Value newSup);
    void remove(WCSP *wcsp, Value val);
    void increase(Value newInf) {increase(NULL, newInf);}
    void decrease(Value newSup) {decrease(NULL, newSup);}
    void remove(Value val) {remove(NULL, val);}
    void assign(Value newValue);

    class iterator;
    friend class iterator;
    class iterator {
        bool enumerated;            // should be a constant
        Variable &var;
        Value value;
        Domain::iterator diter;
    public:
        iterator(Variable &v, Value vv) : enumerated(false), var(v), value(vv), diter(NULL) {}
        iterator(Variable &v, Domain::iterator iter) : enumerated(true), var(v), value(v.sup+1), diter(iter) {}

        Value operator*() const {if (enumerated) return *diter; else return value;}
        
        inline iterator &operator++() {    // Prefix form
            if (enumerated) {
                if (var.inf != var.sup) ++diter;
                else {
                    if (*diter < var.value) diter = var.domain.lower_bound(var.value);
                    else diter = var.domain.end();
                }
            } else {
                if (value < var.sup) ++value;
                else value = var.sup + 1;
            }
            return *this;
        }
        
        iterator &operator--() {    // Prefix form
            if (enumerated) {
                if (var.unassigned()) --diter;
                else {
                    if (*diter > var.value) diter = var.domain.lower_bound(var.value);
                    else diter = var.domain.end();
                }
            } else {
                if (value > var.inf) --value;
                else value = var.sup + 1;
            }
            return *this;
        }

        // To see if you're at the end:
        bool operator==(const iterator &iter) const {
            if (enumerated) return diter == iter.diter;
            else return value == iter.value;
        }
        bool operator!=(const iterator &iter) const {
            if (enumerated) return diter != iter.diter;
            else return value != iter.value;
        }
    };
    iterator begin() {
        if (enumerated) {
            if (assigned()) return iterator(*this, domain.lower_bound(value));
            else return iterator(*this, domain.begin());
        } else return iterator(*this, inf);
    }
    iterator end() {
        if (enumerated) return iterator(*this, domain.end()); 
        else return iterator(*this, sup + 1);
    }
    iterator rbegin() {
        if (enumerated) {
            if (assigned()) return iterator(*this, domain.upper_bound(value));
            else return iterator(*this, domain.rbegin());
        } else return iterator(*this, sup);
    }
    iterator rend() {return end();}

    //Finds the first available element whose value is greater or equal to v
    iterator lower_bound(Value v) {
        if (enumerated) {
            if (assigned()) {
                if (v <= value) return iterator(*this, domain.lower_bound(value));
                else return end();
            } else return iterator(*this, domain.lower_bound(v));
        } else {
            if (v <= sup) return iterator(*this, v);
            else return end();
        }
    }

    //Finds the first available element whose value is lower or equal to v
    iterator upper_bound(Value v) {
        if (enumerated) {
            if (assigned()) {
                if (v >= value) return iterator(*this, domain.upper_bound(value));
                else return end();
            } else return iterator(*this, domain.upper_bound(v));
        } else {
            if (v >= inf) return iterator(*this, v);
            else return end();
        }
    }

    int findWCSPIndex(WCSP *wcsp) {
        for (unsigned int i=0; i<wcsps.size(); i++) {
            if (wcsps[i].wcsp == wcsp) return wcsps[i].wcspIndex;
        }
        return -1;
    }
    void addWCSP(WCSP *wcsp, int index) {WCSPLink w; w.wcsp = wcsp; w.wcspIndex = index; wcsps.push_back(w);}
    
    int getDegree();

    friend ostream& operator<<(ostream& os, Variable &var) {
        os << var.name;
        if (var.enumerated && var.unassigned()) {
            os << " " << var.domain;
        } else {
            os << " [" << var.inf << "," << var.sup << "]";
        }
        os << "/" << var.getDegree();
        return os;
    }
};

/*
 * For internal use only! Interaction between tb2store and tb2btlist
 * 
 */

template <class T, class V> template <class Q> void StoreStack<T,V>::restore(BTList<Q> **l, DLink<Q> **elt, int &x)
{
    if (elt[x] == NULL) {
        l[x]->undoPushBack();
    } else {
        assert(l[x] == l[x-1]);
        l[x]->undoErase(elt[x],elt[x-1]);
        x--;
    }
}

#endif /*TB2VARIABLE_HPP_*/

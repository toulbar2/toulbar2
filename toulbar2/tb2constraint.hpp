/** \file tb2constraint.hpp
 *  \brief Generic virtual constraint.
 *
 */

#ifndef TB2CONSTRAINT_HPP_
#define TB2CONSTRAINT_HPP_

#include "tb2types.hpp"

class Constraint : public WCSPLink
{
    // make it private because we don't want copy nor assignment
    Constraint(const Constraint &c);
    Constraint& operator=(const Constraint &c);

public:
    Constraint() {}
    virtual ~Constraint() {}

    // remove a constraint from the set of active constraints
    virtual void deconnect() {cout << "dummy deconnect on (" << this << ")!" << endl;}
    virtual void reconnect() {cout << "dummy reconnect on (" << this << ")!" << endl;}
    virtual bool connected() {cout << "dummy connected on (" << this << ")!" << endl;return true;}
    virtual bool deconnected() {cout << "dummy deconnected on (" << this << ")!" << endl;return false;}

    virtual int arity() const = 0;
    virtual CostVariable *getCostVar(int index) const = 0;
    // return the smallest wcsp index in the constraint scope except for one variable having a forbidden scope index
    virtual int getSmallestVarIndexInScope(int forbiddenScopeIndex) = 0;

    virtual void propagate() = 0;
    virtual void increase(int index) {propagate();}
    virtual void decrease(int index) {propagate();}
    virtual void remove(int index) {propagate();}
    virtual void projectFromZero(int index) {}
    virtual void assign(int index) {propagate();}

    virtual bool verify() {return true;};
    
    virtual void print(ostream& os) {os << this << " Unknown constraint!";}

    friend ostream& operator<<(ostream& os, Constraint &c) {
        c.print(os);
        return os;
    }
};

#endif /*TB2CONSTRAINT_HPP_*/

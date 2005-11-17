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

    virtual void propagate() = 0;
    virtual void propagate(int varIndex) = 0;
    virtual bool verify() {return true;};

    virtual void increase(int index) {cout << "dummy increase on (" << this << "," << index << ")!" << endl;}
    virtual void decrease(int index) {cout << "dummy decrease on (" << this << "," << index << ")!" << endl;}
    virtual void assign(int index) {cout << "dummy assign on (" << this << "," << index << ")!" << endl;}
    virtual void remove(int index) {cout << "dummy remove on (" << this << "," << index << ")!" << endl;}
    
    virtual void print(ostream& os) {os << this << " Unknown constraint!";}

    friend ostream& operator<<(ostream& os, Constraint &c) {
        c.print(os);
        return os;
    }
};

#endif /*TB2CONSTRAINT_HPP_*/

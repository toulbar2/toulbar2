/** \file tb2abstractconstr.hpp
 *  \brief Abstract constraints of predefined arities
 *
 */

#ifndef TB2ABSTRACTCONSTR_HPP_
#define TB2ABSTRACTCONSTR_HPP_

#include "tb2constraint.hpp"
#include "tb2variable.hpp"

template<class T1, class T2>
class AbstractBinaryConstraint : public Constraint
{
protected:
    T1 *x;
    T2 *y;
    DLink<ConstraintLink> *linkX;
    DLink<ConstraintLink> *linkY;
    
public:
    AbstractBinaryConstraint(WCSP *wcsp, T1 *xx, T2 *yy) : Constraint(wcsp), x(xx), y(yy), linkX(NULL), linkY(NULL) {
        assert(xx != yy);
        linkX = xx->link(this,0);
        linkY = yy->link(this,1);
    }

    virtual ~AbstractBinaryConstraint() {delete linkX; delete linkY;}

    bool connected() {return !linkX->removed && !linkY->removed;}
    bool deconnected() {return linkX->removed || linkY->removed;}
    void deconnect() {
        if (connected()) {
            if (ToulBar2::verbose >= 3) cout << "deconnect " << this << endl; 
            x->getConstrs()->erase(linkX, true); 
            y->getConstrs()->erase(linkY, true);
        }
    }
    void reconnect() {
        if (deconnected()) {
            x->getConstrs()->push_back(linkX, true); 
            y->getConstrs()->push_back(linkY, true);
        }
    }

    int arity() const {return 2;}
    
    Variable *getVar(int varIndex) const {return (varIndex == 0)?x:y;}

    int getSmallestVarIndexInScope(int forbiddenScopeIndex) {return (forbiddenScopeIndex)?x->wcspIndex:y->wcspIndex;}
};

#endif /*TB2ABSTRACTCONSTR_HPP_*/

/** \file tb2abstractconstr.hpp
 *  \brief Abstract constraints of predefined arities
 *
 */

#ifndef TB2ABSTRACTCONSTR_HPP_
#define TB2ABSTRACTCONSTR_HPP_

#include "tb2constraint.hpp"
#include "tb2wcsp.hpp"

class AbstractBinaryConstraint : public Constraint
{
protected:
    CostVariable *x;
    CostVariable *y;
    DLink<ConstraintLink> *linkX;
    DLink<ConstraintLink> *linkY;
    
    CostVariable *getX() const {return x;}
    CostVariable *getY() const {return y;}
    
public:
    AbstractBinaryConstraint(CostVariable *xx, CostVariable *yy);

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
    
    CostVariable *getCostVar(int varIndex) const {return (varIndex == 0)?getX():getY();}

    int getSmallestVarIndexInScope(int forbiddenScopeIndex);    
};

#endif /*TB2ABSTRACTCONSTR_HPP_*/

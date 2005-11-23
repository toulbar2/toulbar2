/** \file tb2queue.hpp
 *  \brief Propagation queue with time stamping.
 * 
 */

#ifndef TB2QUEUE_HPP_
#define TB2QUEUE_HPP_

#include "tb2btlist.hpp"

struct CostVariableWithTimeStamp 
{
    CostVariable *var;
    long long timeStamp;
};

class Queue : private BTList<CostVariableWithTimeStamp>
{  
    // make it private because we don't want copy nor assignment
    Queue(const Queue &s);
    Queue& operator=(const Queue &s);
    
public:
    Queue() : BTList<CostVariableWithTimeStamp>(NULL) {}
    
    int getSize() const {return BTList<CostVariableWithTimeStamp>::getSize();}
    bool empty() const {return BTList<CostVariableWithTimeStamp>::empty();}

    void clear() {BTList<CostVariableWithTimeStamp>::clear();}
    
    void push(DLink<CostVariableWithTimeStamp> *elt, long long curTimeStamp) {
        if (elt->content.timeStamp < curTimeStamp) {
            elt->content.timeStamp = curTimeStamp;
            push_back(elt, false);
        }
    }
    
    CostVariable *pop() {
        assert(!empty());
        DLink<CostVariableWithTimeStamp> *elt = pop_back(false);
        elt->content.timeStamp = -1;
        return elt->content.var;
    }
    
    CostVariable *pop_min();
};

#endif /*TB2QUEUE_HPP_*/

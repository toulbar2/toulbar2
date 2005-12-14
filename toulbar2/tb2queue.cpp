/*
 * ****** Propagation queue with time stamping *******
 */
 
#include "tb2queue.hpp"
#include "tb2costvar.hpp"

CostVariable *Queue::pop_min()
{
    assert(!empty());
    iterator iter=begin();
    DLink<CostVariableWithTimeStamp> *elt = iter.getElt();
    int pos = (*iter).var->wcspIndex;
    for (++iter; iter != end(); ++iter) {
        if ((*iter).var->wcspIndex < pos) {
            elt = iter.getElt();
            pos = (*iter).var->wcspIndex;
        }
    }
    erase(elt, false);
    elt->content.timeStamp = -1;
    elt->content.incdec = NOTHING_EVENT;
    return elt->content.var;
}

CostVariable *Queue::pop_min(int *incdec)
{
    assert(!empty());
    iterator iter=begin();
    DLink<CostVariableWithTimeStamp> *elt = iter.getElt();
    int pos = (*iter).var->wcspIndex;
    for (++iter; iter != end(); ++iter) {
        if ((*iter).var->wcspIndex < pos) {
            elt = iter.getElt();
            pos = (*iter).var->wcspIndex;
        }
    }
    erase(elt, false);
    elt->content.timeStamp = -1;
    *incdec = elt->content.incdec;
    elt->content.incdec = NOTHING_EVENT;
    return elt->content.var;
}

CostVariable *Queue::pop_max()
{
    assert(!empty());
    iterator iter=begin();
    DLink<CostVariableWithTimeStamp> *elt = iter.getElt();
    int pos = (*iter).var->wcspIndex;
    for (++iter; iter != end(); ++iter) {
        if ((*iter).var->wcspIndex > pos) {
            elt = iter.getElt();
            pos = (*iter).var->wcspIndex;
        }
    }
    erase(elt, false);
    elt->content.timeStamp = -1;
    elt->content.incdec = NOTHING_EVENT;
    return elt->content.var;
}

CostVariable *Queue::pop_max(int *incdec)
{
    assert(!empty());
    iterator iter=begin();
    DLink<CostVariableWithTimeStamp> *elt = iter.getElt();
    int pos = (*iter).var->wcspIndex;
    for (++iter; iter != end(); ++iter) {
        if ((*iter).var->wcspIndex > pos) {
            elt = iter.getElt();
            pos = (*iter).var->wcspIndex;
        }
    }
    erase(elt, false);
    elt->content.timeStamp = -1;
    *incdec = elt->content.incdec;
    elt->content.incdec = NOTHING_EVENT;
    return elt->content.var;
}

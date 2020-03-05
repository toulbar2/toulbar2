/*
 * ****** Propagation queue with time stamping *******
 */

#include "tb2queue.hpp"
#include "core/tb2variable.hpp"

void Queue::push(DLink<VariableWithTimeStamp>* elt, Long curTimeStamp)
{
    if (elt->content.timeStamp < curTimeStamp) {
        elt->content.timeStamp = curTimeStamp;
        push_back(elt, false);
    }
}

void Queue::push(DLink<VariableWithTimeStamp>* elt, EventType incdec, Long curTimeStamp)
{
    elt->content.incdec |= incdec;
    push(elt, curTimeStamp);
}

void Queue::remove(DLink<VariableWithTimeStamp>* elt)
{
    elt->content.timeStamp = -1;
    elt->content.incdec = NOTHING_EVENT;
    erase(elt, false);
}

Variable* Queue::pop()
{
    assert(!empty());
    DLink<VariableWithTimeStamp>* elt = pop_back(false);
    elt->content.timeStamp = -1;
    elt->content.incdec = NOTHING_EVENT;
    return elt->content.var;
}

Variable* Queue::pop(int* incdec)
{
    assert(!empty());
    *incdec = (*rbegin()).incdec;
    return pop();
}

Variable* Queue::pop_min()
{
    assert(!empty());
    iterator iter = begin();
    DLink<VariableWithTimeStamp>* elt = iter.getElt();
    int pos = (*iter).var->getDACOrder();
    for (++iter; iter != end(); ++iter) {
        if ((*iter).var->getDACOrder() < pos) {
            elt = iter.getElt();
            pos = (*iter).var->getDACOrder();
        }
    }
    erase(elt, false);
    elt->content.timeStamp = -1;
    elt->content.incdec = NOTHING_EVENT;
    return elt->content.var;
}

Variable* Queue::pop_min(int* incdec)
{
    assert(!empty());
    iterator iter = begin();
    DLink<VariableWithTimeStamp>* elt = iter.getElt();
    int pos = (*iter).var->getDACOrder();
    for (++iter; iter != end(); ++iter) {
        if ((*iter).var->getDACOrder() < pos) {
            elt = iter.getElt();
            pos = (*iter).var->getDACOrder();
        }
    }
    erase(elt, false);
    elt->content.timeStamp = -1;
    *incdec = elt->content.incdec;
    elt->content.incdec = NOTHING_EVENT;
    return elt->content.var;
}

Variable* Queue::pop_max()
{
    assert(!empty());
    iterator iter = begin();
    DLink<VariableWithTimeStamp>* elt = iter.getElt();
    int pos = (*iter).var->getDACOrder();
    for (++iter; iter != end(); ++iter) {
        if ((*iter).var->getDACOrder() > pos) {
            elt = iter.getElt();
            pos = (*iter).var->getDACOrder();
        }
    }
    erase(elt, false);
    elt->content.timeStamp = -1;
    elt->content.incdec = NOTHING_EVENT;
    return elt->content.var;
}

Variable* Queue::pop_max(int* incdec)
{
    assert(!empty());
    iterator iter = begin();
    DLink<VariableWithTimeStamp>* elt = iter.getElt();
    int pos = (*iter).var->getDACOrder();
    for (++iter; iter != end(); ++iter) {
        if ((*iter).var->getDACOrder() > pos) {
            elt = iter.getElt();
            pos = (*iter).var->getDACOrder();
        }
    }
    erase(elt, false);
    elt->content.timeStamp = -1;
    *incdec = elt->content.incdec;
    elt->content.incdec = NOTHING_EVENT;
    return elt->content.var;
}

Variable* Queue::pop_first()
{
    assert(!empty());
    iterator iter = begin();
    DLink<VariableWithTimeStamp>* elt = iter.getElt();
    erase(elt, false);
    elt->content.timeStamp = -1;
    elt->content.incdec = NOTHING_EVENT;
    return elt->content.var;
}

int cmpVariableDAC(const void* p1, const void* p2)
{
    DLink<VariableWithTimeStamp>* c1 = *((DLink<VariableWithTimeStamp>**)p1);
    DLink<VariableWithTimeStamp>* c2 = *((DLink<VariableWithTimeStamp>**)p2);
    int v1 = c1->content.var->getDACOrder();
    int v2 = c2->content.var->getDACOrder();
    if (v1 > v2)
        return 1;
    else if (v1 < v2)
        return -1;
    else
        return 0;
}

int cmpVariableRevDAC(const void* p1, const void* p2)
{
    DLink<VariableWithTimeStamp>* c1 = *((DLink<VariableWithTimeStamp>**)p1);
    DLink<VariableWithTimeStamp>* c2 = *((DLink<VariableWithTimeStamp>**)p2);
    int v1 = c1->content.var->getDACOrder();
    int v2 = c2->content.var->getDACOrder();
    if (v1 < v2)
        return 1;
    else if (v1 > v2)
        return -1;
    else
        return 0;
}

void Queue::sort(bool increase)
{
    int size = getSize();
    DLink<VariableWithTimeStamp>** sorted = new DLink<VariableWithTimeStamp>*[size]; // replace size by MAX_BRANCH_SIZE in case of compilation problem
    int i = 0;
    for (iterator iter = begin(); iter != end(); ++iter) {
        sorted[i++] = iter.getElt();
    }
    qsort(sorted, size, sizeof(DLink<VariableWithTimeStamp>*), (increase) ? cmpVariableDAC : cmpVariableRevDAC);
    for (int i = 0; i < size; i++) {
        erase(sorted[i], false);
        push_back(sorted[i], false);
    }
    delete[] sorted;
}

void Queue::print(ostream& os)
{
    os << "Queue: ";
    iterator iter = begin();
    if (iter != end()) {
        VariableWithTimeStamp vts = iter.getElt()->content;
        os << "<var:" << vts.var->getName() << ",node:" << vts.timeStamp << "> ";
        for (++iter; iter != end(); ++iter) {
            vts = iter.getElt()->content;
            os << "<var:" << vts.var->getName() << ",node:" << vts.timeStamp << "> ";
        }
    }
    os << endl;
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

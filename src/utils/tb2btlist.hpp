/** \file tb2btlist.hpp
 *  \brief Backtrackable double-linked list.
 *
 * Convention:
 *
 * elements can be inserted at the end of the list only
 * these insertions can be undone in the reverse order of their insertion
 *
 * an exception to the previous rule is when inserting an element permanently during search (using freeze mode),
 * this element will be inserted at the beginning of the list.
 *
 * elements can be removed in any order
 * these removals can be undone in the reverse order of their removal.
 *
 */

#ifndef TB2BTLIST_HPP_
#define TB2BTLIST_HPP_

#include "tb2store.hpp"

template <class T>
class DLink {
public:
    bool removed; // true if the corresponding element has been removed
    DLink* next;
    DLink* prev;
    T content;

public:
    DLink()
        : removed(true)
        , next(NULL)
        , prev(NULL)
    {
    }
};

template <class T>
class BTList {
    StoreStack<BTList, DLink<T>*>* storeUndo;
    int size;
    DLink<T>* head;
    DLink<T>* last;

public:
    BTList(StoreStack<BTList, DLink<T>*>* s)
        : storeUndo(s)
        , size(0)
        , head(NULL)
        , last(NULL)
    {
    }

    int getSize() const { return size; }
    bool empty() const { return size == 0; }

    // Warning! clear() is not a backtrackable operation
    void clear()
    {
        size = 0;
        head = NULL;
        last = NULL;
    }

    bool inBTList(DLink<T>* elt)
    {
        for (iterator iter = begin(); iter != end(); ++iter) {
            if (elt == iter.getElt())
                return !elt->removed;
        }
        return false;
    }

    void push_front(DLink<T>* elt)
    {
        assert(!inBTList(elt));
        size++;
        elt->removed = false;
        if (head != NULL) {
            head->prev = elt;
            elt->next = head;
        } else {
            last = elt;
            elt->next = NULL;
        }
        head = elt;
        head->prev = NULL;
        assert(head != NULL || last == NULL);
        assert(last != NULL || head == NULL);
    }

    void push_back(DLink<T>* elt, bool backtrack)
    {
        assert(!backtrack || storeUndo);
        if (storeUndo && storeUndo->isFrozen() && (Store::getDepth() > 0)) { // learned constraints during search are pushed in front to not interfere with backtrackable constraints (coming from variable elimination or projectNary) which are pushed back.
            push_front(elt);
            return;
        }
        assert(!inBTList(elt));
        size++;
        elt->removed = false;
        if (last != NULL) {
            last->next = elt;
            elt->prev = last;
        } else {
            head = elt;
            elt->prev = NULL;
        }
        last = elt;
        last->next = NULL;
        if (backtrack)
            storeUndo->store(this, NULL);
        assert(head != NULL || last == NULL);
        assert(last != NULL || head == NULL);
    }

    void undoPushBack()
    {
        assert(last != NULL);
        size--;
        last->removed = true;
        if (last->prev != NULL) {
            last = last->prev;
            last->next->prev = NULL;
            last->next = NULL;
        } else {
            head = NULL;
            last = NULL;
        }
        assert(head != NULL || last == NULL);
        assert(last != NULL || head == NULL);
    }

    void erase(DLink<T>* elt, bool backtrack)
    {
        assert(!elt->removed);
        size--;
        elt->removed = true;
        if (elt->prev != NULL) {
            assert(!elt->prev->removed);
            assert(elt->prev->next == elt);
            elt->prev->next = elt->next;
        } else
            head = elt->next;
        if (elt->next != NULL) {
            assert(!elt->next->removed);
            assert(elt->next->prev == elt);
            elt->next->prev = elt->prev;
        } else
            last = elt->prev;
        if (backtrack) {
            storeUndo->store(this, elt->next);
            storeUndo->store(this, elt);
        }
        assert(head != NULL || last == NULL);
        assert(last != NULL || head == NULL);
    }

    // older version with no learned constraints and only push_back operations
    //    void undoErase(DLink<T>* elt, DLink<T>* prev)
    //    {
    //        assert(elt->removed);
    //        size++;
    //        elt->removed = false;
    //        if (prev != NULL) {
    //            assert(!prev->removed);
    //            elt->prev = prev;
    //            elt->next = prev->next;
    //            if (prev->next != NULL)
    //                prev->next->prev = elt;
    //            else
    //                last = elt;
    //            prev->next = elt;
    //        } else {
    //            if (head != NULL)
    //                head->prev = elt;
    //            else
    //                last = elt;
    //            elt->prev = NULL;
    //            elt->next = head;
    //            head = elt;
    //        }
    //        assert(head != NULL || last == NULL);
    //        assert(last != NULL || head == NULL);
    //    }

    void undoErase(DLink<T>* elt, DLink<T>* next)
    {
        assert(elt->removed);
        size++;
        elt->removed = false;
        if (next != NULL) {
            assert(!next->removed);
            elt->next = next;
            elt->prev = next->prev;
            if (next->prev != NULL) {
                assert(!next->prev->removed);
                next->prev->next = elt;
            } else {
                assert(head == next);
                head = elt;
            }
            next->prev = elt;
        } else {
            if (last != NULL) {
                assert(!last->removed);
                last->next = elt;
            } else {
                head = elt;
            }
            elt->prev = last;
            elt->next = NULL;
            last = elt;
        }
        assert(head != NULL || last == NULL);
        assert(last != NULL || head == NULL);
    }

    // does not work because elt->prev and elt->next may be null???
    //    void undoErase(DLink<T> *elt) {
    //        assert(elt->removed);
    //        size++;
    //        elt->removed = false;
    //        if (elt->prev != NULL) {
    //            assert(!elt->prev->removed);
    //            assert(elt->prev->next == elt->next);
    //            elt->prev->next = elt;
    //        } else {
    //            assert(head == elt->next);
    //            head = elt;
    //        }
    //        if (elt->next != NULL) {
    //            assert(!elt->next->removed);
    //            assert(elt->next->prev == elt->prev);
    //            elt->next->prev = elt;
    //        } else {
    //            assert(last == elt->prev);
    //            last = elt;
    //        }
    //        assert(head != NULL || last == NULL);
    //        assert(last != NULL || head == NULL);
    //    }

    DLink<T>* pop_back(bool backtrack)
    {
        assert(last != NULL);
        DLink<T>* oldlast = last;
        erase(last, backtrack);
        return oldlast;
    }

    class iterator {
        DLink<T>* elt;

    public:
        iterator() { elt = NULL; }
        iterator(DLink<T>* e)
            : elt(e)
        {
        }

        T operator*() const
        {
            assert(elt != NULL);
            return elt->content;
        }

        DLink<T>* getElt() const { return elt; }

        iterator& operator++()
        { // Prefix form
            if (elt != NULL) {
                while (elt->next != NULL && elt->next->removed) {
                    elt = elt->next;
                }
                elt = elt->next;
            }
            assert(elt == NULL || !elt->removed);
            return *this;
        }

        iterator& operator--()
        { // Prefix form
            if (elt != NULL) {
                while (elt->prev != NULL && elt->prev->removed) {
                    elt = elt->prev;
                }
                elt = elt->prev;
            }
            assert(elt == NULL || !elt->removed);
            return *this;
        }

        // To see if you're at the end:
        bool operator==(const iterator& iter) const { return elt == iter.elt; }
        bool operator!=(const iterator& iter) const { return elt != iter.elt; }
    };

    iterator begin() { return iterator(head); }
    iterator end() { return iterator(NULL); }
    iterator rbegin() { return iterator(last); }
    iterator rend() { return end(); }
};

typedef BTList<ConstraintLink> ConstraintList;
typedef BTList<Variable*> VariableList;
typedef BTList<Separator*> SeparatorList;
typedef BTList<KnapsackConstraint*> KnapsackList;

/*
 * For internal use only! Interaction between tb2store and tb2btlist
 *
 */

template <class T, class V>
template <class Q>
void StoreStack<T, V>::restore(BTList<Q>** l, DLink<Q>** elt, ptrdiff_t& x)
{
    if (elt[x] == NULL) {
        l[x]->undoPushBack();
    } else {
        assert(l[x] == l[x - 1]);
        l[x]->undoErase(elt[x], elt[x - 1]);
        x--;
    }
}

#endif /*TB2BTLIST_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

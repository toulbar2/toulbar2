/** \file tb2globaldecomposable.hpp
 *  \brief Decomposable global cost functions : WeightedRegular, WeightedAmong
 */

#ifndef TB2WFA_HPP_
#define TB2WFA_HPP_

#include "core/tb2wcsp.hpp"
#include "core/tb2types.hpp"
#include "core/tb2enumvar.hpp"

struct WTransition {
    unsigned int start;
    unsigned int end;
    unsigned int symbol;
    Cost weight;

    WTransition(unsigned int _start, unsigned int _end, unsigned int _symbol, Cost _weight)
    {
        start = _start;
        end = _end;
        symbol = _symbol;
        weight = _weight;
    }

    void display()
    {
        cout << start << " x " << symbol << " --(" << weight << ")--> " << end << endl;
    }
};

class WFA {
private:
    unsigned int nbStates;
    list<pair<int, Cost>> initialStates;
    list<pair<int, Cost>> acceptingStates;
    list<WTransition*> transitions;

public:
    WFA();
    WFA(int _nbStates);
    WFA(istream& file, bool mult = true);
    WFA(int nbSymbol, string forbiddenPattern, Cost cost);
    
    ~WFA();

    inline unsigned int getNbStates() { return nbStates; }
    inline list<pair<int, Cost>>& getInitialStates() { return initialStates; }
    inline list<pair<int, Cost>>& getAcceptingStates() { return acceptingStates; }
    inline list<WTransition*>& getTransitions() { return transitions; }

    void display();
};

#endif

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

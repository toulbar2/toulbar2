/** \file tb2vac.hpp
 *  \brief Enforce VAC in a WCSP.
 */

#ifndef TB2VAC_HPP_
#define TB2VAC_HPP_

#include "utils/tb2queue.hpp"

class tVACStat;
class VACVariable;
class VACBinaryConstraint;
class VACTernaryConstraint;

typedef map<Cost, int> tScale;

/**
 * The class that enforces VAC
 */
class VACExtension {

private:
    WCSP* wcsp; /**< Reference to the WCSP that will be processed */
    Queue VAC; /**< Non backtrackable list; AC2001 queue used for Pass1 inside VAC */
#ifdef INCREMENTALVAC
    Queue VAC2; /**< Non backtrackable list; AC2001 queue used for incremental VAC */
#endif
    Queue SeekSupport; /**< Non backtrackable list; collect all variables with a deletion caused by binary constraints during Pass1 */
    Long nbIterations; /**< Incremented at each pass, used as a TimeStamp */
    int inconsistentVariable; /**< WipeOut variable, Used also to check after enforcePass1() if the network is VAC */
    int PBconflict;
    vector<tuple<VACVariable*, Value, bool>> acSupport;
    vector<pair<int, Value>> PBkillersctr;
    vector<pair<int, Value>> killers;
    vector<pair<int, Value>> killed;
    vector<pair<pair<int, Value>, Cost>> EPT;

    Cost prevItThreshold; /**< The previous cost threshold (theta) for the iterative threshold descent */
    Cost itThreshold; /**< The cost threshold (theta) for the iterative threshold descent */
    int breakCycles; /**< Number of iterations with no c0 increase */
    tScale scaleCost; /**w The list of all costs used in the WCSP ? */
    list<Cost> scaleVAC; /**< The scale of costs used for the thresold descent */

    Cost minlambda; /**< The amount of cost that will go to c0 (lambda) */

    stack<pair<int, Value>>* queueP; /**< Values removed by hard AC (created in Pass1, used in Pass2) */
    stack<pair<int, Value>>* queueR; /**< Minimal set of deletions needed to increase c0 (created in Pass2, used in Pass3) */

    void enforcePass1(); /**< Enforces instrumented hard AC (Phase 1) */
    bool enforcePass1(VACVariable* xj, VACBinaryConstraint* cij); /**< Revises /a xj wrt /a cij and updates /a k */
    bool checkPass1() const; /**< Checks if Bool(P) is AC */
    void enforcePass2(); /**< Finds a minimal set of deletions needed for wipeout and computes k and lambda */
    bool enforcePass3(); /**< Project and extends costs to increase c0 according to the plan */
    void enforcePass3VACDecomposition(); /**< Enforces VAC decomposition pass 3 (substract cost and decrease top) */

    void clear(); /**< empty all VAC queues */
    void resetSupports(); /**< reset binary supports for AC2001 optimal O(ed^2) complexity */
    bool enqueueVAC(Cost threshold, Cost previousThreshold); /**< selects more variables for AC2001 queue having unary costs between threshold and previousThreshold ; returns false if nothing to be done */

    Cost sumlb;
    Long nlb;
    Long sumvars;
    int sumk;
    int theMaxK;

    int bneckVar;
    VACBinaryConstraint* bneckCF;
    Cost bneckCost;

    set<int> tempvars;

public:
    VACExtension(WCSP* w);
    ~VACExtension();

    bool firstTime() { return nbIterations == 0; } /**< Is it the first iteration ? */

    bool isVAC() const { return (inconsistentVariable == -1); } /**< Is the WCSP VAC-epsilon ? Pass1 must be enforced */
    bool propagate(); /**< Starts one VAC iteration */
    bool isNull(Cost c) const { return (c < itThreshold); } /**< is the Cost significant (above itThreshold aka theta) */

    void queueVAC(DLink<VariableWithTimeStamp>* link);
#ifdef INCREMENTALVAC
    void queueVAC2(DLink<VariableWithTimeStamp>* link);
#endif
    void queueSeekSupport(DLink<VariableWithTimeStamp>* link);

    void init();
    void iniThreshold(); /**< Initialize itThreshold to the strongest cost in the cost scale */
    void iniThreshold(Cost c);
    Cost getThreshold() { return itThreshold; }
    void nextScaleCost(); /**< Sets ItThreshold to the next scale */
    void histogram(Cost c);
    void histogram(); /**< Computes the ScaleVAC splitting the cost scale in 20 buckets or less */

    set<Long> singletonI;
    set<Long> singleton;

    void iniSingleton();
    void updateSingleton();
    void removeSingleton();

    void printStat(bool ini = false);
    void printTightMatrix();

    void minsumDiffusion(); /**< MinSumDiffusion implementation */

    Cost RASPSFindItThreshold();
};

#endif /*TB2VAC_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

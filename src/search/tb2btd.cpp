/*
 * **************** Backtrack and Russian Doll Search with Tree Decomposition *******************
 *
 * \warning EDAC is not applied on duplicated binary or ternary constraints
 */

#include "tb2solver.hpp"
#include "core/tb2domain.hpp"
#include "tb2clusters.hpp"
#ifdef OPENMPI
#include "vns/tb2cpdgvns.hpp"
#endif

/*
 * Variable ordering heuristics
 *
 */

int Solver::getNextUnassignedVar(Cluster* cluster)
{
    if (unassignedVars->empty())
        return -1;
    TVarsSorted::iterator iter, iter_begin, iter_end;
    if (!cluster->getFreedom()) {
        iter_begin = cluster->beginSortedVars();
        iter_end = cluster->endSortedVars();
    } else {
        iter_begin = cluster->beginSortedVarsTree();
        iter_end = cluster->endSortedVarsTree();
    }

    for (iter = iter_begin; iter != iter_end; ++iter) {
        if (wcsp->unassigned(*iter))
            return *iter;
    }
    return -1;
}

int Solver::getVarMinDomainDivMaxDegree(Cluster* cluster)
{
    if (unassignedVars->empty())
        return -1;
    TVarsSorted::iterator iter, iter_begin, iter_end;
    if (!cluster->getFreedom()) {
        iter_begin = cluster->beginSortedVars();
        iter_end = cluster->endSortedVars();
    } else {
        iter_begin = cluster->beginSortedVarsTree();
        iter_end = cluster->endSortedVarsTree();
    }
    int varIndex = -1;
    Cost worstUnaryCost = MIN_COST;
    double best = MAX_VAL - MIN_VAL;

    for (iter = iter_begin; iter != iter_end; ++iter) {
        if (wcsp->unassigned(*iter)) {
            int deg = wcsp->getDegree(*iter) + 1; // - ((WCSP *)wcsp)->getVar(*iter)->nbSeparators();
            double heuristic = (double)wcsp->getDomainSize(*iter) / (double)max(deg, 1);
            if (varIndex < 0 || heuristic < best - (double)ToulBar2::epsilon * best
                || (heuristic < best + (double)ToulBar2::epsilon * best && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
                best = heuristic;
                varIndex = *iter;
                worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
            }
        }
    }
    return varIndex;
}

int Solver::getVarMinDomainDivMaxDegreeRandomized(Cluster* cluster)
{
    if (unassignedVars->empty())
        return -1;
    TVarsSorted::iterator iter, iter_begin, iter_end;
    if (!cluster->getFreedom()) {
        iter_begin = cluster->beginSortedVars();
        iter_end = cluster->endSortedVars();
    } else {
        iter_begin = cluster->beginSortedVarsTree();
        iter_end = cluster->endSortedVarsTree();
    }
    int varIndex = -1;
    Cost worstUnaryCost = MIN_COST;
    double best = MAX_VAL - MIN_VAL;
    int ties[cluster->getNbVars()];
    int nbties = 0;

    for (iter = iter_begin; iter != iter_end; ++iter) {
        if (wcsp->unassigned(*iter)) {
            int deg = wcsp->getDegree(*iter) + 1; // - ((WCSP *)wcsp)->getVar(*iter)->nbSeparators();
            double heuristic = (double)wcsp->getDomainSize(*iter) / (double)max(deg, 1);
            if (varIndex < 0 || heuristic < best - (double)ToulBar2::epsilon * best
                || (heuristic < best + (double)ToulBar2::epsilon * best && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
                best = heuristic;
                varIndex = *iter;
                nbties = 1;
                ties[0] = varIndex;
                worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
            } else if (heuristic < best + (double)ToulBar2::epsilon * best && wcsp->getMaxUnaryCost(*iter) == worstUnaryCost) {
                ties[nbties] = *iter;
                nbties++;
            }
        }
    }
    if (nbties > 1)
        return ties[myrand() % nbties];
    else
        return varIndex;
}

int Solver::getVarMinDomainDivMaxDegreeLastConflict(Cluster* cluster)
{
    if (unassignedVars->empty())
        return -1;

    if (lastConflictVar != -1 && wcsp->unassigned(lastConflictVar) && ((!cluster->getFreedom() && cluster->isVar(lastConflictVar)) || (cluster->getFreedom() && cluster->isVarTree(lastConflictVar))))
        return lastConflictVar;

    TVarsSorted::iterator iter, iter_begin, iter_end;
    if (!cluster->getFreedom()) {
        iter_begin = cluster->beginSortedVars();
        iter_end = cluster->endSortedVars();
    } else {
        iter_begin = cluster->beginSortedVarsTree();
        iter_end = cluster->endSortedVarsTree();
    }
    int varIndex = -1;
    Cost worstUnaryCost = MIN_COST;
    double best = MAX_VAL - MIN_VAL;

    for (iter = iter_begin; iter != iter_end; ++iter) {
        if (wcsp->unassigned(*iter)) {
            int deg = wcsp->getDegree(*iter) + 1; // - ((WCSP *)wcsp)->getVar(*iter)->nbSeparators();
            double heuristic = (double)wcsp->getDomainSize(*iter) / (double)max(deg, 1);
            if (varIndex < 0 || heuristic < best - (double)ToulBar2::epsilon * best
                || (heuristic < best + (double)ToulBar2::epsilon * best && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
                best = heuristic;
                varIndex = *iter;
                worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
            }
        }
    }
    return varIndex;
}

int Solver::getVarMinDomainDivMaxDegreeLastConflictRandomized(Cluster* cluster)
{
    if (unassignedVars->empty())
        return -1;

    if (lastConflictVar != -1 && wcsp->unassigned(lastConflictVar) && ((!cluster->getFreedom() && cluster->isVar(lastConflictVar)) || (cluster->getFreedom() && cluster->isVarTree(lastConflictVar))))
        return lastConflictVar;

    TVarsSorted::iterator iter, iter_begin, iter_end;
    if (!cluster->getFreedom()) {
        iter_begin = cluster->beginSortedVars();
        iter_end = cluster->endSortedVars();
    } else {
        iter_begin = cluster->beginSortedVarsTree();
        iter_end = cluster->endSortedVarsTree();
    }

    int varIndex = -1;
    Cost worstUnaryCost = MIN_COST;
    double best = MAX_VAL - MIN_VAL;
    int ties[cluster->getNbVars()];
    int nbties = 0;

    for (iter = iter_begin; iter != iter_end; ++iter) {
        if (wcsp->unassigned(*iter)) {
            int deg = wcsp->getDegree(*iter) + 1; // - ((WCSP *)wcsp)->getVar(*iter)->nbSeparators();
            double heuristic = (double)wcsp->getDomainSize(*iter) / (double)max(deg, 1);
            if (varIndex < 0 || heuristic < best - (double)ToulBar2::epsilon * best
                || (heuristic < best + (double)ToulBar2::epsilon * best && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
                best = heuristic;
                varIndex = *iter;
                nbties = 1;
                ties[0] = varIndex;
                worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
            } else if (heuristic < best + (double)ToulBar2::epsilon * best && wcsp->getMaxUnaryCost(*iter) == worstUnaryCost) {
                ties[nbties] = *iter;
                nbties++;
            }
        }
    }
    if (nbties > 1)
        return ties[myrand() % nbties];
    else
        return varIndex;
}

int Solver::getVarMinDomainDivMaxWeightedDegree(Cluster* cluster)
{
    if (unassignedVars->empty())
        return -1;
    TVarsSorted::iterator iter, iter_begin, iter_end;
    if (!cluster->getFreedom()) {
        iter_begin = cluster->beginSortedVars();
        iter_end = cluster->endSortedVars();
    } else {
        iter_begin = cluster->beginSortedVarsTree();
        iter_end = cluster->endSortedVarsTree();
    }
    int varIndex = -1;
    Cost worstUnaryCost = MIN_COST;
    double best = MAX_VAL - MIN_VAL;

    for (iter = iter_begin; iter != iter_end; ++iter) {
        if (wcsp->unassigned(*iter)) {
            Cost unarymediancost = MIN_COST;
            int domsize = wcsp->getDomainSize(*iter);
            if (ToulBar2::weightedTightness) {
                ValueCost array[domsize];
                wcsp->getEnumDomainAndCost(*iter, array);
                unarymediancost = stochastic_selection<ValueCost>(array, 0, domsize - 1, domsize / 2).cost;
            }
            Long deg = wcsp->getWeightedDegree(*iter) + 1 + unarymediancost; // - ((WCSP *)wcsp)->getVar(*iter)->nbSeparators();
            double heuristic = (double)domsize / (double)max(deg, (Long)1);
            if (varIndex < 0 || heuristic < best - (double)ToulBar2::epsilon * best
                || (heuristic < best + (double)ToulBar2::epsilon * best && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
                best = heuristic;
                varIndex = *iter;
                worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
            }
        }
    }
    return varIndex;
}

int Solver::getVarMinDomainDivMaxWeightedDegreeRandomized(Cluster* cluster)
{
    if (unassignedVars->empty())
        return -1;
    TVarsSorted::iterator iter, iter_begin, iter_end;
    if (!cluster->getFreedom()) {
        iter_begin = cluster->beginSortedVars();
        iter_end = cluster->endSortedVars();
    } else {
        iter_begin = cluster->beginSortedVarsTree();
        iter_end = cluster->endSortedVarsTree();
    }
    int varIndex = -1;
    Cost worstUnaryCost = MIN_COST;
    double best = MAX_VAL - MIN_VAL;
    int ties[cluster->getNbVars()];
    int nbties = 0;

    for (iter = iter_begin; iter != iter_end; ++iter) {
        if (wcsp->unassigned(*iter)) {
            Cost unarymediancost = MIN_COST;
            int domsize = wcsp->getDomainSize(*iter);
            if (ToulBar2::weightedTightness) {
                ValueCost array[domsize];
                wcsp->getEnumDomainAndCost(*iter, array);
                unarymediancost = stochastic_selection<ValueCost>(array, 0, domsize - 1, domsize / 2).cost;
            }
            Long deg = wcsp->getWeightedDegree(*iter) + 1 + unarymediancost; // - ((WCSP *)wcsp)->getVar(*iter)->nbSeparators();
            double heuristic = (double)domsize / (double)max(deg, (Long)1);
            if (varIndex < 0 || heuristic < best - (double)ToulBar2::epsilon * best
                || (heuristic < best + (double)ToulBar2::epsilon * best && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
                best = heuristic;
                varIndex = *iter;
                nbties = 1;
                ties[0] = varIndex;
                worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
            } else if (heuristic < best + (double)ToulBar2::epsilon * best && wcsp->getMaxUnaryCost(*iter) == worstUnaryCost) {
                ties[nbties] = *iter;
                nbties++;
            }
        }
    }
    if (nbties > 1)
        return ties[myrand() % nbties];
    else
        return varIndex;
}

int Solver::getVarMinDomainDivMaxWeightedDegreeLastConflict(Cluster* cluster)
{
    if (unassignedVars->empty())
        return -1;

    if (lastConflictVar != -1 && wcsp->unassigned(lastConflictVar) && ((!cluster->getFreedom() && cluster->isVar(lastConflictVar)) || (cluster->getFreedom() && cluster->isVarTree(lastConflictVar))))
        return lastConflictVar;

    TVarsSorted::iterator iter, iter_begin, iter_end;
    if (!cluster->getFreedom()) {
        iter_begin = cluster->beginSortedVars();
        iter_end = cluster->endSortedVars();
    } else {
        iter_begin = cluster->beginSortedVarsTree();
        iter_end = cluster->endSortedVarsTree();
    }
    int varIndex = -1;
    Cost worstUnaryCost = MIN_COST;
    double best = MAX_VAL - MIN_VAL;

    for (iter = iter_begin; iter != iter_end; ++iter) {
        if (wcsp->unassigned(*iter)) {
            Cost unarymediancost = MIN_COST;
            int domsize = wcsp->getDomainSize(*iter);
            if (ToulBar2::weightedTightness) {
                ValueCost array[domsize];
                wcsp->getEnumDomainAndCost(*iter, array);
                unarymediancost = stochastic_selection<ValueCost>(array, 0, domsize - 1, domsize / 2).cost;
            }
            double heuristic = (double)domsize / (double)(wcsp->getWeightedDegree(*iter) + 1 + unarymediancost);
            //	      Long deg = wcsp->getWeightedDegree(*iter) + 1; // - ((WCSP *)wcsp)->getVar(*iter)->nbSeparators();
            //        double heuristic = (double) wcsp->getDomainSize(*iter) / (double) max(deg,(Long)1);
            if (varIndex < 0 || heuristic < best - (double)ToulBar2::epsilon * best
                || (heuristic < best + (double)ToulBar2::epsilon * best && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
                best = heuristic;
                varIndex = *iter;
                worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
            }
        }
    }
    return varIndex;
}

int Solver::getVarMinDomainDivMaxWeightedDegreeLastConflictRandomized(Cluster* cluster)
{
    if (unassignedVars->empty())
        return -1;

    if (lastConflictVar != -1 && wcsp->unassigned(lastConflictVar) && ((!cluster->getFreedom() && cluster->isVar(lastConflictVar)) || (cluster->getFreedom() && cluster->isVarTree(lastConflictVar))))
        return lastConflictVar;

    TVarsSorted::iterator iter, iter_begin, iter_end;
    if (!cluster->getFreedom()) {
        iter_begin = cluster->beginSortedVars();
        iter_end = cluster->endSortedVars();
    } else {
        iter_begin = cluster->beginSortedVarsTree();
        iter_end = cluster->endSortedVarsTree();
    }
    int varIndex = -1;
    Cost worstUnaryCost = MIN_COST;
    double best = MAX_VAL - MIN_VAL;
    int ties[cluster->getNbVars()];
    int nbties = 0;

    for (iter = iter_begin; iter != iter_end; ++iter) {
        if (wcsp->unassigned(*iter)) {
            Cost unarymediancost = MIN_COST;
            int domsize = wcsp->getDomainSize(*iter);
            if (ToulBar2::weightedTightness) {
                ValueCost array[domsize];
                wcsp->getEnumDomainAndCost(*iter, array);
                unarymediancost = stochastic_selection<ValueCost>(array, 0, domsize - 1, domsize / 2).cost;
            }
            double heuristic = (double)domsize / (double)(wcsp->getWeightedDegree(*iter) + 1 + unarymediancost);
            //	      Long deg = wcsp->getWeightedDegree(*iter) + 1; // - ((WCSP *)wcsp)->getVar(*iter)->nbSeparators();
            //        double heuristic = (double) wcsp->getDomainSize(*iter) / (double) max(deg,(Long)1);
            if (varIndex < 0 || heuristic < best - (double)ToulBar2::epsilon * best
                || (heuristic < best + (double)ToulBar2::epsilon * best && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
                best = heuristic;
                varIndex = *iter;
                nbties = 1;
                ties[0] = varIndex;
                worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
            } else if (heuristic < best + (double)ToulBar2::epsilon * best && wcsp->getMaxUnaryCost(*iter) == worstUnaryCost) {
                ties[nbties] = *iter;
                nbties++;
            }
        }
    }
    if (nbties > 1)
        return ties[myrand() % nbties];
    else
        return varIndex;
}

// Managing freedom of clusters
void Solver::Manage_Freedom(Cluster* c)
{
    if (c->getNbVars() == 0) {
        c->setFreedom(false);
        return;
    }

    // via separator
    bool recorded = false;
    recorded = c->freeGet();

    if (c->isLeaf()) {
        c->freeRec(false);
    } else {
        if (!recorded || c->getFreedom()) {
            bool found = false;
            for (TClusters::iterator iter = c->beginDescendants(); iter != c->endDescendants() && !found; ++iter) {
                if ((*iter)->getId() != c->getId()) {
                    if ((*iter)->getSep()->used()) {
                        found = true;
                    }
                }
            }
            if (found) { // some separator inside the current subtree has been used to compute the current lower bound
                if (ToulBar2::verbose >= 1)
                    cout << " propagation has already used some separator nogood inside the current subtree of cluster " << c->getId() << endl;
                if (recorded && c->getFreedom()) {
                    if (c->open) {
                        *(c->open) = OpenList();
                    }
                    nbForcedChoiceChange++;
                } else {
                    nbForcedChoices++;
                }
                c->setFreedom(false);
            } else {
                if (!recorded) {
                    nbChoices++;
                }

                bool before = c->getFreedom();

                c->freeRecInc();

                if (recorded && before) {
                    if (!c->getFreedom()) {
                        if (c->open) {
                            *(c->open) = OpenList();
                        }
                        nbChoiceChange++;
                    } else {
                        nbReadOnly++;
                    }
                }
            }
        } else {
            nbReadOnly++;
        }
    }
}

/*
 * Choice points
 *
 */

pair<Cost, Cost> Solver::binaryChoicePoint(Cluster* cluster, Cost lbgood, Cost cub, int varIndex, Value value)
{
    assert(lbgood < cub);
    if (ToulBar2::interrupted)
        throw TimeOut();
    Cost clb = cub;
    TreeDecomposition* td = wcsp->getTreeDec();
    assert(wcsp->unassigned(varIndex));
    assert(wcsp->canbe(varIndex, value));
    bool dichotomic = (ToulBar2::dichotomicBranching && ToulBar2::dichotomicBranchingSize < wcsp->getDomainSize(varIndex));
    Value middle = value;
    bool increasing = true;
    if (dichotomic) {
        middle = (wcsp->getInf(varIndex) + wcsp->getSup(varIndex)) / 2;
        if (value <= middle)
            increasing = true;
        else
            increasing = false;
    }
    try {
        Store::store();
        assert(td->getCurrentCluster() == cluster);
        assert(wcsp->getLb() == cluster->getLbRec());
        wcsp->setUb(cub);
        Cost bestlb = lbgood;
        if (CUT(bestlb, cub))
            THROWCONTRADICTION;
        if (ToulBar2::btdMode >= 2) {
            Cost rds = td->getLbRecRDS();
            bestlb = MAX(bestlb, rds);
            if (CUT(bestlb, cub))
                THROWCONTRADICTION;
        }
        lastConflictVar = varIndex;
        if (dichotomic) {
            if (increasing)
                decrease(varIndex, middle);
            else
                increase(varIndex, middle + 1);
        } else
            assign(varIndex, value);
        lastConflictVar = -1;
        bestlb = MAX(bestlb, wcsp->getLb());
        pair<Cost, Cost> res = recursiveSolve(cluster, bestlb, cub);
        clb = MIN(res.first, clb);
        cub = MIN(res.second, cub);
    } catch (const Contradiction&) {
        wcsp->whenContradiction();
    }
    Store::restore();
    nbBacktracks++;
    if (nbBacktracks > nbBacktracksLimit)
        throw NbBacktracksOut();
#ifdef OPENMPI
    if (ToulBar2::parallel && ((nbBacktracks % 128) == 0) && MPI_interrupted())
        throw TimeOut();
#endif
    cluster->nbBacktracks++;
    try {
        Store::store();
        assert(wcsp->getTreeDec()->getCurrentCluster() == cluster);
        assert(wcsp->getLb() == cluster->getLbRec());
        wcsp->setUb(cub);
        Cost bestlb = lbgood;
        if (CUT(bestlb, cub))
            THROWCONTRADICTION;
        if (ToulBar2::btdMode >= 2) {
            Cost rds = td->getLbRecRDS();
            bestlb = MAX(bestlb, rds);
            if (CUT(bestlb, cub))
                THROWCONTRADICTION;
        }
        if (dichotomic) {
            if (increasing)
                increase(varIndex, middle + 1, cluster->nbBacktracks >= cluster->hbfsLimit || nbBacktracks >= cluster->hbfsGlobalLimit);
            else
                decrease(varIndex, middle, cluster->nbBacktracks >= cluster->hbfsLimit || nbBacktracks >= cluster->hbfsGlobalLimit);
        } else
            remove(varIndex, value, cluster->nbBacktracks >= cluster->hbfsLimit || nbBacktracks >= cluster->hbfsGlobalLimit);
        bestlb = MAX(bestlb, wcsp->getLb());
        if (!ToulBar2::hbfs && cluster == td->getRoot() && initialDepth + 1 == Store::getDepth()) {
            initialDepth++;
            showGap(bestlb, cub);
        };
        if (cluster->nbBacktracks >= cluster->hbfsLimit || nbBacktracks >= cluster->hbfsGlobalLimit) {
            addOpenNode(*(cluster->cp), *(cluster->open), bestlb, cluster->getCurrentDeltaUb());
        } else {
            pair<Cost, Cost> res = recursiveSolve(cluster, bestlb, cub);
            clb = MIN(res.first, clb);
            cub = MIN(res.second, cub);
        }
    } catch (const Contradiction&) {
        wcsp->whenContradiction();
    }
    Store::restore();
    assert(lbgood <= clb);
    assert(ToulBar2::bilevel || clb <= cub);
    return make_pair(clb, cub);
}

/*
 * Choice points for solution counting
 *
 */

BigInteger Solver::binaryChoicePointSBTD(Cluster* cluster, int varIndex, Value value)
{
    if (ToulBar2::interrupted)
        throw TimeOut();
    Cost cub = 1;
    Cost lbgood = 0;
    BigInteger nbSol = 0, nb = 0;
    assert(wcsp->unassigned(varIndex));
    assert(wcsp->canbe(varIndex, value));
    bool dichotomic = (ToulBar2::dichotomicBranching && ToulBar2::dichotomicBranchingSize < wcsp->getDomainSize(varIndex));
    Value middle = value;
    bool increasing = true;
    if (dichotomic) {
        middle = (wcsp->getInf(varIndex) + wcsp->getSup(varIndex)) / 2;
        if (value <= middle)
            increasing = true;
        else
            increasing = false;
    }
    try {
        Store::store();
        assert(wcsp->getTreeDec()->getCurrentCluster() == cluster);

        wcsp->setUb(cub);
        if (CUT(lbgood, cub))
            THROWCONTRADICTION;
        lastConflictVar = varIndex;
        if (dichotomic) {
            if (increasing)
                decrease(varIndex, middle);
            else
                increase(varIndex, middle + 1);
        } else
            assign(varIndex, value);
        lastConflictVar = -1;
        nb = sharpBTD(cluster);
        nbSol += nb;
    } catch (const Contradiction&) {
        wcsp->whenContradiction();
    }
    Store::restore();
    nbBacktracks++;
    if (nbBacktracks > nbBacktracksLimit)
        throw NbBacktracksOut();
#ifdef OPENMPI
    if (ToulBar2::parallel && ((nbBacktracks % 128) == 0) && MPI_interrupted())
        throw TimeOut();
#endif
    try {
        Store::store();
        assert(wcsp->getTreeDec()->getCurrentCluster() == cluster);
        //		assert(wcsp->getLb() == cluster->getLbRec());
        wcsp->setUb(cub);
        if (CUT(lbgood, cub))
            THROWCONTRADICTION;
        if (dichotomic) {
            if (increasing)
                increase(varIndex, middle + 1);
            else
                decrease(varIndex, middle);
        } else
            remove(varIndex, value);

        nb = sharpBTD(cluster);
        nbSol += nb;
    } catch (const Contradiction&) {
        wcsp->whenContradiction();
    }
    Store::restore();
    return nbSol;
}

/*
 * Backtrack with Tree Decomposition
 *
 */

/// \defgroup bilevel Bilevel optimization
/// We assume a tree decomposition with four clusters: root cluster 0 is Problem0, left child cluster 1 is Problem1, middle child cluster 2 is Problem2, and right child cluster 3 is NegProblem2 (i.e., -Problem2)
/// Our goal is to minimize the leader problem (Problem0 + Problem1 - Problem2) such that the follower problem (Problem2) takes its minimum cost w.r.t. the leader decisions (locally minimizing Problem2 after Problem0 is fixed)
/// For each child cluster, we have:
/// - an initial lower bound computed in preprocessing (see \ref Toulbar2::initialLbBLP)
/// - a shifting cost value in order to deal with non-negative costs (see \ref ToulBar2::negCostBLP)
///
/// We maintain the following properties during search:
/// - bilevelShiftP0 = wcsp->getNegativeLb() - ToulBar2::initialLbP1 - Toulbar2::initialLbP2 - \ref Toulbar2::initialLbNegP2
/// - deltaP1 = upper bound on costs that have been moved from Problem1 to Problem0 by soft local consistencies
/// - deltaNegP2 = lower bound on costs that has been moved from NegProblem2 to Problem0 by soft local consistencies
/// - cluster0.lb <= Problem0 + initialLbNegP2 + deltaNegP2 + bilevelShiftP0 + bilevelShiftP1 + bilevelShiftP2
/// - cluster1.lb <= Problem1 - deltaP1 + bilevelShiftP1 -initialLbP1
/// - cluster2.lb <= Problem2 + bilevelShiftP2 - initialLbP2
/// - cluster3.lb <= -Problem2 + bilevelShiftNegP2 - initialLbNegP2 - deltaNegP2
/// Thus, during search in root cluster, we always have:
/// - wcsp->getLb() = cluster0.lb + cluster1.lb + cluster3.lb <= Problem0 + Problem1 -Problem2 + bilevelShiftP0 + bilevelShiftP1 + bilevelShiftP2 + bilevelShiftNegP2 <= wcsp->getUb()
/// - wcsp->getLb() = cluster0.lb + cluster2.lb <= Problem1 - Problem2  + wcsp->getNegativeLb() <= wcsp->getUb()
///
/// \note Maintains the best (monotonically increasing) lower bound of the cluster in parameter lbgood
pair<Cost, Cost> Solver::recursiveSolve(Cluster* cluster, Cost lbgood, Cost cub)
{
    if (ToulBar2::verbose >= 1)
        cout << "[" << Store::getDepth() << "] recursive solve     cluster: " << cluster->getId() << "     clb: " << lbgood << "     cub: " << cub << "     clb0: " << cluster->getLb() << "     wcsp->lb: " << wcsp->getLb() << "     wcsp->ub: " << wcsp->getUb() << endl;
    assert(lbgood <= cub);
    TreeDecomposition* td = wcsp->getTreeDec();
    int varIndex = -1;

    // dummy virtual root has no freedom
    assert(!ToulBar2::heuristicFreedom || cluster != td->getRoot() || cluster->getNbVars() == 0);
    assert(!ToulBar2::heuristicFreedom || cluster != td->getRoot() || cluster->getFreedom() == false);

    if (ToulBar2::heuristicFreedom && cluster->getDepth() > solveDepth) {
        solveDepth = cluster->getDepth();
    }

    if (ToulBar2::Static_variable_ordering)
        varIndex = getNextUnassignedVar(cluster);
    else if (ToulBar2::weightedDegree && ToulBar2::lastConflict)
        varIndex = ((ToulBar2::restart > 0) ? getVarMinDomainDivMaxWeightedDegreeLastConflictRandomized(cluster) : getVarMinDomainDivMaxWeightedDegreeLastConflict(cluster));
    else if (ToulBar2::lastConflict)
        varIndex = ((ToulBar2::restart > 0) ? getVarMinDomainDivMaxDegreeLastConflictRandomized(cluster) : getVarMinDomainDivMaxDegreeLastConflict(cluster));
    else if (ToulBar2::weightedDegree)
        varIndex = ((ToulBar2::restart > 0) ? getVarMinDomainDivMaxWeightedDegreeRandomized(cluster) : getVarMinDomainDivMaxWeightedDegree(cluster));
    else
        varIndex = ((ToulBar2::restart > 0) ? getVarMinDomainDivMaxDegreeRandomized(cluster) : getVarMinDomainDivMaxDegree(cluster));

    if (varIndex < 0) {
        // Current cluster is completely assigned (if no freedom), current cluster descendant is completely assigned (if freedom)
        Cost clb = wcsp->getLb();
        assert(clb <= cub);
        Cost csol = clb;

        if (!cluster->getFreedom()) {
            for (TClusters::iterator iter = cluster->beginSortedEdges(); clb < cub && iter != cluster->endSortedEdges();) {
                // Solves each cluster son with local lower and upper bounds
                Cluster* c = *iter;
                ++iter;
                if (ToulBar2::bilevel && td->getRoot() == cluster && c == *cluster->rbeginEdges())
                    continue; // skip the last cluster son representing the negative follower problem
                Cost lbSon = MIN_COST;
                Cost ubSon = MAX_COST;
                bool good = false;
                if (!c->isActive()) {
                    c->reactivate();
                    c->nogoodGet(lbSon, ubSon, &c->open);
                    good = true;
                } else {
                    lbSon = c->getLbRec();
                    if (!ToulBar2::bilevel) { // TODO: why not reusing nogoods on the leader problem?
                        ubSon = c->getUb();
#ifndef NDEBUG
                        Cost dummylb = -MAX_COST;
                        Cost tmpub = -MAX_COST;
                        c->nogoodGet(dummylb, tmpub, &c->open);
                        assert(tmpub == ubSon);
#endif
                    }
                }
                if (ToulBar2::verbose >= 2)
                    cout << "lbson: " << lbSon << " ubson: " << ubSon << " lbgood:" << lbgood << " clb: " << clb << " csol: " << csol << " cub: " << cub << " cluster->lb: " << c->getLbRec() << endl;
                if (lbSon < ubSon) { // we do not have an optimality proof
                    if (ToulBar2::bilevel || clb <= lbgood || (csol < MAX_COST && ubSon >= cub - csol + lbSon)) { // we do not know a good enough son's solution or the currently reconstructed father's solution is not working or the currently reconstructed father's lower bound is not increasing
                        if (ToulBar2::heuristicFreedom) {
                            Manage_Freedom(c);
                            td->updateInTD(c);
                        }
                        bool csolution = (csol < MAX_COST && ubSon < cub - csol + lbSon);
                        assert(!csolution || ubSon < cub - clb + lbSon);
                        ubSon = MIN(ubSon, cub - clb + lbSon); // this rule is not valid for bilevel optimization
                        td->setCurrentCluster(c);
                        wcsp->setUb(ubSon);
                        wcsp->setLb((good) ? c->getLbRec() : lbSon); // FIXME: good=true in bilevel but c->getLbRec is always zero???
                        // Compute an initial bound for the follower problem
                        assert(!ToulBar2::bilevel || td->getRoot() != cluster || c != *cluster->rbeginEdges()); // child cluster c cannot be NegP2
                        Cost bestLbP2 = MIN_COST;
                        if (ToulBar2::bilevel && td->getRoot() == cluster && c != *cluster->beginEdges()) { // initialize current bounds for the follower problem
                            assert(c != *cluster->rbeginEdges());
                            Cost deltaNegP2Lb = (*cluster->rbeginEdges())->getCurrentDeltaLb();
                            Cost deltaNegP2Ub = (*cluster->rbeginEdges())->getCurrentDeltaUb();
                            assert(c->getLbRec() == MIN_COST);
                            Cost lbP1 = cluster->getLb() + (*cluster->beginEdges())->getLb() - ToulBar2::initialLbBLP[2] - deltaNegP2Ub; // FIXME: P1 is completely solved (use cluster son recorded lb instead of propagation lb)
                            Cost lbNegP2 = (*cluster->rbeginEdges())->getLb() + ToulBar2::initialLbBLP[2] + deltaNegP2Lb;
                            ubSon = -ToulBar2::initialLbBLP[1] - lbNegP2 + ToulBar2::negCostBLP[1] + ToulBar2::negCostBLP[2] + UNIT_COST;
                            assert(ubSon >= UNIT_COST);
                            bestLbP2 = max(MIN_COST, lbP1 - cub + ToulBar2::negCostBLP[1] + ToulBar2::negCostBLP[2]);
                            assert(bestLbP2 < ubSon);
                            // if (ToulBar2::verbose>=1 && bestLbP2>MIN_COST) cout << bestLbP2 << " <= P2 < "<< ubSon << endl;
                            lbSon = bestLbP2;
                            wcsp->setUb(ubSon);
                            wcsp->setLb(MIN_COST);
                        }
                        int depth = Store::getDepth();
                        try {
                            Store::store();
                            wcsp->enforceUb();
                            if (ToulBar2::bilevel && td->getRoot() == cluster && c != *cluster->beginEdges()) {
                                // add channeling constraints between leader (Problem0) and follower (Problem2) problems
                                for (int ctrIndex : ((WCSP*)wcsp)->delayedCtrBLP[1]) {
                                    Constraint* ctr = ((WCSP*)wcsp)->getCtr(ctrIndex);
                                    assert(ctr->deconnected());
                                    static vector<Cost> costs;
                                    Constraint* incCtr = NULL;
                                    if (ctr->isBinary()) {
                                        costs.resize(ctr->getDomainInitSizeProduct(), MIN_COST);
                                        incCtr = ((WCSP*)wcsp)->getCtr(wcsp->postIncrementalBinaryConstraint(ctr->getVar(0)->wcspIndex, ctr->getVar(1)->wcspIndex, costs));
                                        ((BinaryConstraint*)incCtr)->addCosts((BinaryConstraint*)ctr);
                                    } else if (ctr->isTernary()) {
                                        costs.resize(ctr->getDomainInitSizeProduct(), MIN_COST);
                                        incCtr = ((WCSP*)wcsp)->getCtr(wcsp->postIncrementalTernaryConstraint(ctr->getVar(0)->wcspIndex, ctr->getVar(1)->wcspIndex, ctr->getVar(2)->wcspIndex, costs));
                                        ((TernaryConstraint*)incCtr)->addCosts((TernaryConstraint*)ctr);
                                    } else {
                                        cerr << "Sorry, bilevel optimization not implemented for this type of channeling cost function:" << *ctr << endl;
                                        throw WrongFileFormat();
                                    }
                                    // incCtr->sumScopeIncluded(ctr);
                                    // incCtr->assignCluster(); //SdG: done before inside postIncrementalXXXConstraint
                                    incCtr->propagate();
                                }
                            }
                            wcsp->propagate();
                            Cost bestlb = MAX(wcsp->getLb(), lbSon);
                            if (!ToulBar2::bilevel && csol < MAX_COST && iter == cluster->endSortedEdges())
                                bestlb = MAX(bestlb, lbgood - csol + lbSon); // simple trick to provide a better initial lower bound for the last son
                            if (ToulBar2::btdMode >= 2) {
                                Cost rds = td->getLbRecRDS();
                                bestlb = MAX(bestlb, rds);
                                if (CUT(bestlb, ubSon))
                                    THROWCONTRADICTION;
                            }
                            pair<Cost, Cost> res = hybridSolve(c, bestlb, ubSon);
                            assert(res.first >= bestlb && res.second <= ubSon);
                            if (ToulBar2::heuristicFreedom && c->getFreedom() && res.first == bestlb && res.second == ubSon && res.first != res.second) {
                                c->freeInc();
                            }
                            c->nogoodRec(res.first, ((res.second < ubSon) ? res.second : MAX_COST), &c->open);
                            if (ToulBar2::bilevel && td->getRoot() == cluster && c != *cluster->beginEdges()) {
                                assert(res.first >= res.second); // leader and follower problems are solved
                                assert(c->getCurrentDeltaUb() == MIN_COST); // we assume no cost moves from the follower to leader problem
                                assert((*cluster->rbeginEdges())->getCurrentDeltaLb() == (*cluster->rbeginEdges())->getCurrentDeltaUb());
                                // cout << "clb: " << clb << " - C3.lb: " << (*cluster->rbeginEdges())->getLb() << " - C3.initlb: " << ToulBar2::initialLbBLP[2] << " - C3.deltalb: " << (*cluster->rbeginEdges())->getCurrentDeltaLb() << " - C2.opt: " << res.first << " - C2.initlb: " << ToulBar2::initialLbBLP[1] << " + C2.negcost: " << ToulBar2::negCostBLP[1] << " + C3.negcost: " << ToulBar2::negCostBLP[2] << endl;
                                csol = clb - (*cluster->rbeginEdges())->getLb() - ToulBar2::initialLbBLP[2] - (*cluster->rbeginEdges())->getCurrentDeltaLb() - res.first - ToulBar2::initialLbBLP[1] + ToulBar2::negCostBLP[1] + ToulBar2::negCostBLP[2];
                                // cout << "csol: " << csol << endl;
                                clb = csol;
                                Store::restore(depth);
                                break;
                            }
                            clb += res.first - lbSon;
                            if (csol < MAX_COST) {
                                if (res.second < ubSon || csolution)
                                    csol += res.second - lbSon;
                                else
                                    csol = MAX_COST;
                            }
                        } catch (const Contradiction&) {
                            wcsp->whenContradiction();
                            c->nogoodRec(ubSon, MAX_COST, &c->open);
                            clb += ubSon - lbSon;
                            if (csolution)
                                csol += ubSon - lbSon;
                            else
                                csol = MAX_COST;
                        }
                        Store::restore(depth);
                    } else {
                        assert(!ToulBar2::bilevel);
                        if (csol < MAX_COST) {
                            assert(ubSon < MAX_COST);
                            csol += ubSon - lbSon;
                        }
                    }
                } else {
                    if (ToulBar2::bilevel && c != *cluster->beginEdges() && c != *cluster->rbeginEdges()) {
                        // the follower problem has been solved already
                        csol = clb - (*cluster->rbeginEdges())->getLb() - ToulBar2::initialLbBLP[2] - (*cluster->rbeginEdges())->getCurrentDeltaLb() - lbSon - ToulBar2::initialLbBLP[1] + ToulBar2::negCostBLP[1] + ToulBar2::negCostBLP[2];
                        clb = csol;
                    }
                }
            }
        }
        assert(csol >= clb);
        if (csol < cub) {
            // A new solution has been found for the current cluster
            cub = csol;
            cluster->solutionRec(csol);
            if (cluster == td->getRoot() || cluster == td->getRootRDS()) {
                if (ToulBar2::verbose >= 0 || ToulBar2::showSolutions) {
                    if (!ToulBar2::bayesian)
                        cout << "New solution: " << std::setprecision(ToulBar2::decimalPoint) << wcsp->Cost2ADCost(csol) << std::setprecision(DECIMAL_POINT) << " (" << nbBacktracks << " backtracks, " << nbNodes << " nodes, depth " << Store::getDepth() << ", " << ((ToulBar2::parallel) ? (realTime() - ToulBar2::startRealTime) : (cpuTime() - ToulBar2::startCpuTime)) << " seconds)" << endl;
                    else
                        cout << "New solution: " << csol << " energy: " << -(wcsp->Cost2LogProb(csol) + ToulBar2::markov_log) << " prob: " << std::scientific << wcsp->Cost2Prob(csol) * Exp(ToulBar2::markov_log) << std::fixed << " (" << nbBacktracks << " backtracks, " << nbNodes << " nodes, depth " << Store::getDepth() << ", " << ((ToulBar2::parallel) ? (realTime() - ToulBar2::startRealTime) : (cpuTime() - ToulBar2::startCpuTime)) << " seconds)" << endl;
                }
                if (cluster == td->getRoot())
                    td->newSolution(csol);
                else {
                    assert(cluster == td->getRootRDS());
                    // Remember current solution for value ordering heuristic
                    wcsp->restoreSolution(cluster);
                    TAssign a;
                    cluster->getSolution(a);
                    if (ToulBar2::showSolutions) {
                        TAssign::iterator it = a.begin();
                        while (it != a.end()) {
                            Value v = it->second;
                            cout << it->first << ":" << v << " ";
                            ++it;
                        }
                        cout << endl;
                    }
                }
            }
        }
        Cost bestlb = MAX(lbgood, clb);
        if (ToulBar2::verbose >= 1)
            cout << "[" << Store::getDepth() << "] C" << cluster->getId() << " return " << bestlb << " " << cub << endl;
        assert(ToulBar2::bilevel || bestlb <= cub);
        if (ToulBar2::hbfs && bestlb < cub) { // keep current node in open list instead of closing it!
            if (cluster->getNbVars() > 0) {
                int varid = -1;
                if (!cluster->getFreedom())
                    varid = *cluster->getVars().begin();
                else
                    varid = *cluster->getVarsTree().begin();
                assert(wcsp->assigned(varid));
                cluster->cp->addChoicePoint(CP_ASSIGN, varid, wcsp->getValue(varid), true); // dummy additional choice point to avoid the reversal of the last effective choice point for this open node
            }
#ifndef NDEBUG
            OpenList* prevopen = cluster->open;
            Cost tmplb = MIN_COST;
            Cost tmpub = MAX_COST;
            Cost tmpclusterub = cluster->getUb();
            assert(cluster == wcsp->getTreeDec()->getRoot() || cluster->nogoodGet(tmplb, tmpub, &cluster->open)); // warning! it can destroy cluster->ub
            cluster->setUb(tmpclusterub);
            assert(prevopen == cluster->open);
#endif
            addOpenNode(*(cluster->cp), *(cluster->open), bestlb, cluster->getCurrentDeltaUb()); // reinsert as a new open node
            cluster->hbfsLimit = cluster->nbBacktracks; // and stop current visited node
            bestlb = cub;
        }
        return make_pair(bestlb, cub);
    } else {
        // Enumerates cluster proper variables (if no freedom), variables in cluster descendants (if freedom)
        *((StoreInt*)searchSize) += ((int)(10e3 * Log(wcsp->getDomainSize(varIndex))));
        pair<Cost, Cost> res = make_pair(MIN_COST, MAX_COST);
        if (wcsp->enumerated(varIndex)) {
            assert(wcsp->canbe(varIndex, wcsp->getSupport(varIndex)));
            // Reuse last solution found if available
            Value bestval = ((ToulBar2::verifyOpt) ? (wcsp->getSup(varIndex) + 1) : wcsp->getBestValue(varIndex));
            res = binaryChoicePoint(cluster, lbgood, cub, varIndex, (wcsp->canbe(varIndex, bestval)) ? bestval : wcsp->getSupport(varIndex));
        } else {
            res = binaryChoicePoint(cluster, lbgood, cub, varIndex, wcsp->getInf(varIndex));
        }
        if (ToulBar2::verbose >= 1)
            cout << "[" << Store::getDepth() << "] C" << cluster->getId() << " return " << res.first << " " << res.second << endl;
        assert(res.first >= lbgood);
        assert(res.second <= cub);
        assert(ToulBar2::bilevel || res.first <= res.second);
        return res;
    }
}

/*
 * Russian Doll Search with Tree Decomposition
 *
 */

pair<Cost, Cost> Solver::russianDollSearch(Cluster* c, Cost cub)
{
    TreeDecomposition* td = wcsp->getTreeDec();
    pair<Cost, Cost> res = make_pair(MIN_COST, cub);

    TClusters::iterator it = c->beginSortedEdges();
    while (it != c->endSortedEdges()) {
        russianDollSearch(*it, cub);
        ++it;
    }

    try {
        Store::store();

        Cost nogoodlb = MIN_COST;
        Cost nogoodub = MAX_COST;
        if (c != td->getRoot()) {
            c->deconnectSep();
            c->nogoodGet(nogoodlb, nogoodub, &c->open); // update c->open and c->ub
            //	      if (nogoodlb == bestub) {
            //	          assert(bestub <= cub);
            //	          Store::restore();
            //	          return make_pair(bestub,bestub);
            //	      }
            assert(c->getLbRec() == MIN_COST);
            c->setLb(MIN_COST);
            wcsp->setLb(MIN_COST);
            td->setCurrentCluster(td->getRoot());
            Cost lbroot = td->getLbRecRDS();
            td->setCurrentCluster(c);
            Cost lbc = td->getLbRecRDS();
            cub = cub - lbroot + lbc;
            cub = MIN(cub, nogoodub);
        }
        wcsp->setUb(cub);
        td->setCurrentCluster(c);
        td->setRootRDS(c);
        lastConflictVar = -1;

        if (ToulBar2::verbose >= 0)
            cout << "--- Solving cluster subtree " << c->getId() << " ..." << endl;

        //	  if(c == td->getRoot()) wcsp->propagate(); // needed if there are connected components
        enforceUb();
        wcsp->propagate();
        Cost bestlb = td->getLbRecRDS();
        bestlb = MAX(bestlb, nogoodlb);
        if (bestlb >= cub)
            THROWCONTRADICTION;
        res = hybridSolve(c, bestlb, cub);
        assert(res.first >= bestlb);
        c->setLbRDS(res.first);
        //	  if (c->sepSize() == 0)  {
        c->nogoodRec(res.first, ((res.second < cub) ? res.second : MAX_COST), &c->open);
        //	  }

        if (ToulBar2::debug || ToulBar2::verbose >= 1)
            c->printStatsRec();
        if (ToulBar2::verbose >= 0)
            cout << "---  done  cost = [" << res.first << "," << res.second << "] (" << nbBacktracks << " backtracks, " << nbNodes << " nodes, depth " << Store::getDepth() << ")" << endl
                 << endl;

    } catch (const Contradiction&) {
        wcsp->whenContradiction();
        res.first = res.second;
        c->setLbRDS(cub);
        //	  if (c->sepSize() == 0) {
        c->nogoodRec(cub, MAX_COST, &c->open);
        //	  }
    }
    Store::restore();
    if (c == td->getRoot()) {
        c->resetLbRec();
    } else {
        if (c->open)
            *(c->open) = OpenList(); // clear current open list
        c->resetUbRec(c);
    }
    return res;
}

/*
 * Backtrack with Tree Decomposition for counting in CSP
 *
 */

BigInteger Solver::sharpBTD(Cluster* cluster)
{

    TreeDecomposition* td = wcsp->getTreeDec();
    BigInteger NbSol = 0, nb = 0;
    TCtrs totalList;
    if (ToulBar2::verbose >= 1)
        cout << "[" << Store::getDepth() << "] recursive solve     cluster: " << cluster->getId() << " **************************************************************" << endl;

    int varIndex = -1;
    if (ToulBar2::Static_variable_ordering)
        varIndex = getNextUnassignedVar(cluster);
    else if (ToulBar2::weightedDegree && ToulBar2::lastConflict)
        varIndex = ((ToulBar2::restart > 0) ? getVarMinDomainDivMaxWeightedDegreeLastConflictRandomized(cluster) : getVarMinDomainDivMaxWeightedDegreeLastConflict(cluster));
    else if (ToulBar2::lastConflict)
        varIndex = ((ToulBar2::restart > 0) ? getVarMinDomainDivMaxDegreeLastConflictRandomized(cluster) : getVarMinDomainDivMaxDegreeLastConflict(cluster));
    else if (ToulBar2::weightedDegree)
        varIndex = ((ToulBar2::restart > 0) ? getVarMinDomainDivMaxWeightedDegreeRandomized(cluster) : getVarMinDomainDivMaxWeightedDegree(cluster));
    else
        varIndex = ((ToulBar2::restart > 0) ? getVarMinDomainDivMaxDegreeRandomized(cluster) : getVarMinDomainDivMaxDegree(cluster));

    if (varIndex < 0) {
        // Current cluster is completely assigned
        Cost lb = wcsp->getLb();
        if (ToulBar2::verbose >= 1)
            cout << "[" << Store::getDepth() << "] C" << cluster->getId() << " lb= " << lb << endl;
        NbSol = 1 * cluster->getCount();

        if (ToulBar2::approximateCountingBTD && cluster->getParent() == NULL)
            totalList = cluster->getCtrsTree();

        for (TClusters::iterator iter = cluster->beginSortedEdges(); NbSol > 0 && iter != cluster->endSortedEdges(); ++iter) {
            // Solves each cluster son
            nb = 0;
            Cluster* c = *iter;
            if ((nb = c->sgoodGet()) != -1) {
                nbSGoodsUse++;
            } else {
                nb = 0;
                td->setCurrentCluster(c);
                try {
                    Store::store();
                    if (ToulBar2::approximateCountingBTD) {
                        if (c->getParent() != NULL && c->getParent()->getParent() == NULL && c->getNbVars() > 1) {
                            // for this son of root, we disconnect the constraints which isn't in intersection
                            TCtrs usefulCtrsList = c->getCtrsTree();
                            c->deconnectDiff(totalList, usefulCtrsList);
                        }
                    }
                    wcsp->propagate();
                    nb = sharpBTD(c);
                    c->sgoodRec(0, nb);
                    nbSGoods++;
                } catch (const Contradiction&) {
                    wcsp->whenContradiction();
                    c->sgoodRec(0, 0); // no solution
                    nbSGoods++;
                }
                Store::restore();
            }
            if (cluster->getParent() == NULL && ToulBar2::approximateCountingBTD) {
                // computation of upper bound of solutions number for each part
                if (ubSol.find(c->getPart()) == ubSol.end()) {
                    ubSol[c->getPart()] = 1;
                }
                ubSol[c->getPart()] *= nb;
            }
            NbSol *= nb;
        }
        return NbSol;
    } else {
        // Enumerates cluster proper variables
        if (wcsp->enumerated(varIndex)) {
            assert(wcsp->canbe(varIndex, wcsp->getSupport(varIndex)));
            // Reuse last solution found if available
            Value bestval = wcsp->getBestValue(varIndex);

            NbSol = binaryChoicePointSBTD(cluster, varIndex, (wcsp->canbe(varIndex, bestval)) ? bestval : wcsp->getSupport(varIndex));
        } else {
            NbSol = binaryChoicePointSBTD(cluster, varIndex, wcsp->getInf(varIndex));
        }
        if (ToulBar2::verbose >= 1)
            cout << "[" << Store::getDepth() << "] C" << cluster->getId() << " return " << NbSol << endl;
        return NbSol;
    }
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

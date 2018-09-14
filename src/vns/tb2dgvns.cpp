/*
 * tb2dgvns.cpp
 *
 *  Created on: 15 December 2016
 *      Author: Abdelkader Ouali
 *      PhD Student: LITIO, University of Oran ; GREYC, University of Caen.
 */

#include "tb2dgvns.hpp"
#include "tb2wcsp.hpp"

bool VNSSolver::solveLS()
{
    ToulBar2::vnsLDSmax = min(ToulBar2::vnsLDSmax, (int) wcsp->numberOfVariables() * (wcsp->getMaxDomainSize()-1));
    ToulBar2::vnsLDSmin = min(ToulBar2::vnsLDSmin, ToulBar2::vnsLDSmax);
    ToulBar2::vnsKmax = min(ToulBar2::vnsKmax, (int) wcsp->numberOfVariables());
    ToulBar2::vnsKmin = min(ToulBar2::vnsKmin, ToulBar2::vnsKmax);
    assert(ToulBar2::vnsLDSmin >= 0);
    assert(ToulBar2::vnsLDSmax >= 0);
    assert(ToulBar2::vnsLDSmin <= ToulBar2::vnsLDSmax);
    assert(ToulBar2::vnsKmin >= 0);
    assert(ToulBar2::vnsKmax >= 0);
    assert(ToulBar2::vnsKmin <= ToulBar2::vnsKmax);

    // Initialization
    if (ToulBar2::DEE)
        ToulBar2::DEE_ = ToulBar2::DEE; // enforces PSNS after closing the model
    wcsp->propagate();
    initVarHeuristic();

    // Compute the Initial solution
    bool complete = false;
    bestUb = generateInitSolution(ToulBar2::vnsInitSol, bestSolution, complete);
    if (ToulBar2::verbose >= 1) cout << "VNS: initial solution with" << ((complete)?" optimal":"") << " cost " << bestUb << endl;

    NeighborhoodStructure* h = NULL;
    switch (ToulBar2::vnsNeighborVarHeur) {
    case RANDOMVAR:
        if (ToulBar2::verbose >= 1) cout << "Random Variables Neighborhood Structure selection" << endl;
        h = new RandomNeighborhoodChoice();
        break;
    case CLUSTERRAND:
        if (ToulBar2::verbose >= 1) cout << "Random Clusters Neighborhood Structure selection" << endl;
        h = new RandomClusterChoice();
        break;
    default:
        cerr << "Unkown Neighborhood Structure" << endl;
        exit(EXIT_FAILURE);
    }
    h->init(wcsp, this);

    //vns/lds+cp
    bool stop = false;
    Long nbRestart = 1;
    Long restart = 1;
    int lds = ToulBar2::vnsLDSmin;
    while (!stop && !complete && bestUb > ToulBar2::vnsOptimum) {
        if (ToulBar2::verbose >= 0 && ToulBar2::restart>1 && ToulBar2::lds) cout << "****** Restart " << nbRestart << " with " << lds << " discrepancies and UB=" << bestUb << " ****** (" << nbNodes << " nodes)" << endl;
        Long rank = 1;
        int k = ToulBar2::vnsKmin;
        while (!complete && k <= ToulBar2::vnsKmax && bestUb > ToulBar2::vnsOptimum) {
            //neighborhood and partial instantiation
            set<int> neighborhood = h->getNeighborhood(k);
            if (ToulBar2::verbose >= 1) {
                cout << "Neighborhood " << k << ": ";
                for (set<int>::iterator it = neighborhood.begin(); it != neighborhood.end(); it++) cout << " " << *it;
                cout << endl;
            }
            vector<int> variables;
            variables.reserve(wcsp->numberOfVariables());
            vector<int> values;
            values.reserve(wcsp->numberOfVariables());
            for (int v = 0; v < (int) wcsp->numberOfVariables(); v++) {
                if (neighborhood.find(v) == neighborhood.end()) {
                    variables.push_back(v);
                    values.push_back(bestSolution[v]);
                }
            }

            //repair
            ToulBar2::vnsKcur = k;
            ToulBar2::vnsLDScur = (ToulBar2::lds)?lds:-1;
            if (ToulBar2::lds) complete = repair_recursiveSolve(lds, variables, values, bestUb);
            else complete = repair_recursiveSolve(variables, values, bestUb);

            //updating
            if (lastUb >= bestUb) {
                if (h->incrementK()) {
                    rank++;
                    if (ToulBar2::restart > 1 && k < ToulBar2::vnsKmax) {
                        switch (ToulBar2::vnsKinc) {
                        case VNS_ADD1:
                            k++;
                            break;
                        case VNS_MULT2:
                            k *= 2;
                            break;
                        case VNS_LUBY:
                            k = ToulBar2::vnsKmin * luby(rank);
                            break;
                        default:
                            cerr << "Unknown neighborhood size increment strategy inside VNS (see option -kinc)!" << endl;
                            exit(EXIT_FAILURE);
                        }
                        k = min(k, ToulBar2::vnsKmax);
                    } else k++;
//                    cout << "rank: " << rank << " luby: " << luby(rank) << " k: " << k << " lds: " << ToulBar2::lds << endl;
                }
            } else {
                rank = 1;
                k = ToulBar2::vnsKmin;
                restart = 1;
                lds = ToulBar2::vnsLDSmin;
                bestUb = lastUb;
                for (int v = 0; v < (int) wcsp->numberOfVariables(); v++) {
                    assert(lastSolution.find(v) != lastSolution.end());
                    bestSolution[v] = lastSolution[v];
                }
                if (ToulBar2::verbose >= 1) cout << "VNS: new solution with cost " << bestUb << endl;
            }
        }
        if (!complete && bestUb > ToulBar2::vnsOptimum) {
            nbRestart++;
            restart++;
            if (nbRestart <= ToulBar2::restart) {
                if (ToulBar2::lds) {
                    switch (ToulBar2::vnsLDSinc) {
                    case VNS_ADD1:
                        lds++;
                        break;
                    case VNS_MULT2:
                        lds *= 2;
                        break;
                    case VNS_LUBY:
                        lds = ToulBar2::vnsLDSmin * luby(restart);
                        break;
                    default:
                        cerr << "Unknown LDS increment strategy inside VNS (see option -linc)!" << endl;
                        exit(EXIT_FAILURE);
                    }
                    lds = min(lds, ToulBar2::vnsLDSmax);
                }
            } else stop = true;
        }
    }

    if (ToulBar2::vnsOutput) ToulBar2::vnsOutput.close();
    if (ToulBar2::verbose >= 0) {
        if (bestUb < MAX_COST) {
            if (ToulBar2::bayesian) cout << (complete?"Optimum: ":"Best upper-bound: ") << bestUb << " log10like: " << (wcsp->Cost2LogProb(bestUb) + ToulBar2::markov_log)/Log(10.) << " prob: " << wcsp->Cost2Prob(bestUb) * Exp(ToulBar2::markov_log) << " in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << ((ToulBar2::DEE)?(" ( "+to_string(wcsp->getNbDEE())+" removals by DEE)"):"") << " and " << cpuTime() - ToulBar2::startCpuTime << " seconds." << endl;
            else cout << (complete?"Optimum: ":"Best upper-bound: ") << bestUb << " in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << ((ToulBar2::DEE)?(" ( "+to_string(wcsp->getNbDEE())+" removals by DEE)"):"") << " and " << cpuTime() - ToulBar2::startCpuTime << " seconds." << endl;
        } else cout << "No solution found by " << ((ToulBar2::vnsNeighborVarHeur==RANDOMVAR)?"":"DG") << "VNS in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << ((ToulBar2::DEE)?(" ( "+to_string(wcsp->getNbDEE())+" removals by DEE)"):"") << " and " << cpuTime() - ToulBar2::startCpuTime << " seconds." << endl;
    }

    return complete && (bestUb < MAX_COST);
}


/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

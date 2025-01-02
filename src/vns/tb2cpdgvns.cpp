/*
 * tb2cpdgvns.cpp
 *
 *  Created on: 3 mars 2015
 *      Author: Abdelkader Ouali
 *      Phd. Student : LITIO, University of Oran. GREYC, University of Caen.
 */

#include "tb2cpdgvns.hpp"
#include "core/tb2wcsp.hpp"
#ifdef OPENMPI

//---------------- Class Definition --------------------------//

bool CooperativeParallelDGVNS::solve(bool first)
{
    // Initialization
    beginSolve(MAX_COST);
    try {
        lastUb = MAX_COST;
        lastSolution.clear();
        if (first) {
            preprocessing(MAX_COST);
        } else {
            if (ToulBar2::elimDegree >= 0)
                ToulBar2::elimDegree_ = ToulBar2::elimDegree;
        }
    } catch (const Contradiction&) {
        wcsp->whenContradiction();
        if (world.rank() == MASTER) {
            if (lastUb < MAX_COST)
                wcsp->setSolution(lastUb, &lastSolution);
            endSolve(lastUb < MAX_COST, lastUb, true);
        }
        return (lastUb < MAX_COST);
    }

    mysrand(abs(ToulBar2::seed) + world.rank());

    world.barrier(); /* IMPORTANT */

    ToulBar2::startRealTimeAfterPreProcessing = realTime();

    if (world.rank() == MASTER) {
        if (ToulBar2::verbose)
            cout << "Run CPDGVNS method." << endl;
        master();
        if (bestUb < MAX_COST)
            wcsp->setSolution(bestUb, &bestSolution);
        endSolve(bestUb < MAX_COST, bestUb, false); // TODO: return complete=true if known by the master
    } else {
        ToulBar2::showSolutions = false;
        slave();
    }

    if (ToulBar2::vnsOutput)
        ToulBar2::vnsOutput << "Search end"
                            << " " << cpuTime() << endl;

    return (bestUb < MAX_COST);
}

//----------------------- Model Master/Slave -----------------------//
void CooperativeParallelDGVNS::master()
{
    // Structure de voisinage basee sur la notion des clusters -  ParallelRandomClusterChoice - et l'initialisation du file - file contenant les clusters pour chaque processus -
    ParallelRandomClusterChoice* h = new ParallelRandomClusterChoice();
    h->init(wcsp, this);

    // verify the number of processes and number of clusters
    int processes = min(world.size(), h->getSize() + 1);

    if (ToulBar2::vnsOutput)
        ToulBar2::vnsOutput << "#param " << ToulBar2::vnsNeighborVarHeur << " " << ToulBar2::vnsKmin << " " << ToulBar2::vnsKmax << " " << ToulBar2::lds << " " << ToulBar2::vnsInitSol << " " << processes << " " << world.size() << " " << h->getSize() << endl;
    // Generation of initial Solution
    bool complete = false;
    bestSolution.clear();
    bestUb = generateInitSolution(ToulBar2::vnsInitSol, bestSolution, complete);
    BestTimeS = (int)cpuTime();
    BestTimeMS = (int)(((Long)cpuTime() * 1000LL) % 1000LL);
    // Get all clusters from the tree decomposition of constraint graph
    file = h->getClustersIndex();

    /* Seed the slaves; send one unit of work (initial solution) to each slave. */
    for (int rank = 0; rank < processes; ++rank)
        if (rank != MASTER) {
            /* choice one free cluster to send it to the processe */
            uint cluster = getCluster();
            int adjcluster = 0;
            int kmax = h->getClustersSize(cluster, adjcluster);
            while (kmax < ToulBar2::vnsKmin) {
                adjcluster++;
                kmax = h->getClustersSize(cluster, adjcluster);
            }
            // cout << kmax << " " << adjcluster << endl ;
            /* Convert initial solution with cluster and kinit parameters in buffer, for each slave process */
            SolutionMessage solmsg(cluster, adjcluster, ToulBar2::vnsKmin, ToulBar2::vnsKmax, BestTimeS, BestTimeMS, bestUb, bestSolution);

            /* Send Initial Solution to each process */
            world.send(rank, WORKTAG, solmsg);
        }

    /* Loop over getting new Best Solutions */
    uint finished = 0;
    uint worker = processes - 1;
    map<int, Value> slastSolution;
    while (finished < worker) {

        /* Receive result (new better solution) from a slave */
        SolutionMessage solmsg;
        mpi::status status = world.recv(mpi::any_source, mpi::any_tag, solmsg);
        /* setting up the best solution in memory, and checkout if we continue or not */
        // Parsing the received buffer
        uint scluster = 0;
        int sk = 0;
        uint snumberclu;
        int skmax;
        int sBestTimeS = 0;
        int sBestTimeMS = 0;
        Cost sbestUb = 0;

        solmsg.get(scluster, snumberclu, sk, skmax, sBestTimeS, sBestTimeMS, sbestUb, slastSolution);

        // getting the cluster back
        file.push_back(scluster);
        scluster = getCluster();
        // Updating
        if (sbestUb <= bestUb
            && !(BestTimeS == sBestTimeS && BestTimeMS == sBestTimeMS)) {
            bestUb = sbestUb;
            BestTimeS = sBestTimeS;
            BestTimeMS = sBestTimeMS;
            snumberclu = 0;
            skmax = h->getClustersSize(scluster, snumberclu);
            while (skmax < ToulBar2::vnsKmin) {
                snumberclu++;
                skmax = h->getClustersSize(scluster, snumberclu);
            }
            for (uint v = 0; v < wcsp->numberOfVariables(); v++) {
                bestSolution[v] = lastSolution[v];
            }
            if (ToulBar2::vnsOutput) {
                ToulBar2::vnsOutput
                    << "# ------------------------------------------"
                    << endl;
                ToulBar2::vnsOutput
                    << "InstanceVnsBestTime "
                    << (BestTimeS + float(float(BestTimeMS) / float(100)))
                    << endl;
                ToulBar2::vnsOutput << "ReceiveTime " << cpuTime()
                                    << "# Confirmation" << endl;
                ToulBar2::vnsOutput << "Cost " << wcsp->Cost2ADCost(bestUb) << endl;
                ToulBar2::vnsOutput << "Solution "
                                    << " (" << sk << " "
                                    << scluster << ") ";
                for (map<int, Value>::iterator it = bestSolution.begin();
                     it != bestSolution.end(); ++it)
                    ToulBar2::vnsOutput << (*it).first << "=" << (*it).second
                                        << " ";
                ToulBar2::vnsOutput << endl;
                ToulBar2::vnsOutput
                    << "# ------------------------------------------"
                    << endl;
            }
        } else {
            if (snumberclu < file.size())
                snumberclu++;
            skmax = h->getClustersSize(scluster, snumberclu);
            while (skmax < ToulBar2::vnsKmin) {
                snumberclu++;
                skmax = h->getClustersSize(scluster, snumberclu);
            }
        }
        if (bestUb > ToulBar2::vnsOptimum) {
            SolutionMessage solmsg(scluster, snumberclu, ToulBar2::vnsKmin, skmax, BestTimeS, BestTimeMS, bestUb, bestSolution);
            /* Send the best solution to the received slave, for next the search */
            world.send(status.source(), WORKTAG, solmsg);
        } else {
            finished++;
            world.send(status.source(), DIETAG, SolutionMessage());
        }
    }
    vector<mpi::request> reqs;
    for (int rank = 0; rank < world.size(); ++rank)
        if (rank != MASTER) {
            // printf("Send finish empty msg to finish with %d\n",rank);
            reqs.push_back(world.isend(rank, DIETAG, SolutionMessage()));
        }
    mpi::wait_all(reqs.begin(), reqs.end());
}

void CooperativeParallelDGVNS::slave()
{
    // Structure de voisinage basée sur la notion des clusters
    // ParallelRandomClusterChoice* h = NeighborhoodStructure::NeighborhoodStructureFactory(VariableHeuristic(hname), static_cast<WCSP*>(wcsp), this);
    ParallelRandomClusterChoice* h = new ParallelRandomClusterChoice();
    h->init(wcsp, this);

    // initializing the timer begin for the hole search of the slave process, in seconds
    double btime = cpuTime();

    wcsp->setUb(MAX_COST);
    while (true) {
        /* Receive a message from the master */
        SolutionMessage solmsg;
        mpi::status status = world.recv(0, mpi::any_tag, solmsg); /* receive from master */

        /* Check the tag of the received message. */
        if (status.tag() == DIETAG) {
            return;
        }

        /* launch Vns/Lds+cp */
        VnsLdsCP(solmsg, btime, h);

        /* send the result back */
        mpi::request req = world.isend(0, 0, solmsg);
        // cout << env0.myrank <<" slave end" << endl ;
        while (!req.test().is_initialized() && !MPI_interrupted())
            ;
    }
}

//------------------- Base functions -----------------------//

void CooperativeParallelDGVNS::VnsLdsCP(SolutionMessage& solmsg, double btime, ParallelRandomClusterChoice* h)
{
    uint currentcluster = 0;
    uint numberclu = 0;
    int kinit = ToulBar2::vnsKmin;
    int kmax = ToulBar2::vnsKmax;

    solmsg.get(currentcluster, numberclu, kinit, kmax, BestTimeS, BestTimeMS, bestUb, bestSolution);

    for (map<int, Value>::iterator it = bestSolution.begin();
         it != bestSolution.end(); ++it)
        lastSolution[(*it).first] = (*it).second;
    if (ToulBar2::verbose >= 1)
        cout << "VNS :: Initial Solution cost " << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->Cost2ADCost(bestUb) << std::setprecision(DECIMAL_POINT) << endl;
    lastUb = bestUb;
    if (ToulBar2::verbose >= 5)
        cout << "VNS :: Initial Solution" << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->Cost2ADCost(bestUb) << std::setprecision(DECIMAL_POINT) << endl;

    // vns/lds+cp
    int k = kinit;
    // cout << "taille maximal du cluster "<< currentcluster<< " "<< numberclu<< " "<< kmax << endl ;
    // cout << env0.myrank <<" slave 1" << endl ;
    // cout << k <<"<="<< kmax <<"&&"<< k <<"<="<< wcsp->numberOfVariables() << "&&"<< ToulBar2::vns_optimum <<"<"  << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->Cost2ADCost(bestUb) << std::setprecision(DECIMAL_POINT) <<"&&"<< (cpuTime()-lbtime)<<endl ;
    for (; k <= kmax && k <= unassignedVars->getSize() && ToulBar2::vnsOptimum < bestUb;) {
        // neighborhood and partial instantiation
        // cout <<"neighborhood"<< " "<<currentcluster<< " "<< numberclu << " " << k << endl ;
        set<int> neighborhood = h->SlaveGetNeighborhood(currentcluster, numberclu, k);

        if (ToulBar2::verbose >= 1) {
            cout << world.rank() << ": LDS " << ToulBar2::lds << " Neighborhood " << k;
            for (set<int>::iterator it = neighborhood.begin(); it != neighborhood.end(); it++)
                cout << " " << *it;
            cout << endl;
        }
        vector<int> variables;
        variables.reserve(unassignedVars->getSize());
        vector<Value> values;
        values.reserve(unassignedVars->getSize());
        for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
            int v = *iter;
            if (wcsp->canbe(v, lastSolution[v]) && neighborhood.find(v) == neighborhood.end()) {
                variables.push_back(v);
                values.push_back(lastSolution[v]);
            }
        }
        // repair
        if (ToulBar2::lds)
            repair_recursiveSolve(ToulBar2::lds, variables, values, bestUb);
        else
            repair_recursiveSolve(variables, values, bestUb);
        // updating
        // cout <<"updating "<< k<< endl ;
        if (lastUb >= bestUb) {
            k++;
        } else {
            bestUb = lastUb;
            BestTimeS = (int)cpuTime();
            BestTimeMS = (int)(((Long)cpuTime() * 1000LL) % 1000LL);
            k = kinit;
            for (uint v = 0; v < wcsp->numberOfVariables(); v++) {
                bestSolution[v] = lastSolution[v];
            }
            try {
                wcsp->updateUb(bestUb);
                wcsp->enforceUb();
                wcsp->propagate();
                if (unassignedVars->getSize() == 0) {
                    lastUb = MAX_COST;
                    lastSolution.clear();
                    ToulBar2::lds = 0;
                    newSolution();
                    if (lastUb < MAX_COST)
                        wcsp->setSolution(lastUb, &lastSolution);
                    bestUb = lastUb;
                    for (uint v = 0; v < wcsp->numberOfVariables(); v++) {
                        bestSolution[v] = lastSolution[v];
                    }
                    break;
                }
            } catch (const Contradiction&) {
                wcsp->whenContradiction();
                break;
            }

        }
    }
    // cout << env0.myrank <<" slave 2" << endl ;
    //  return the best solution found
    solmsg = SolutionMessage(currentcluster, numberclu, k, kmax, BestTimeS, BestTimeMS, bestUb, bestSolution);
}

/*----------------------------- clustring functions ------------------------------- */
uint CooperativeParallelDGVNS::getCluster()
{
    assert(file.size() != 0);
    uint c = file[0];
    file.erase(file.begin());
    return c;
}
#endif

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

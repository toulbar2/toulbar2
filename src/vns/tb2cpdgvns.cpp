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

// Conversion Tools
// Solution to Message
void CooperativeParallelDGVNS::SolToMsg(
    MPIEnv& env0, uint cluster, uint numberclu, int kinit, int kmax,
    Cost bestUb, int sec, int msec, map<int, Value>& bestSolution)
{
    env0.sendbuff[0] = cluster;
    env0.sendbuff[1] = kinit;
    env0.sendbuff[2] = kmax;
    env0.sendbuff[3] = sec;
    env0.sendbuff[4] = msec;
    env0.sendbuff[5] = numberclu;
    stringstream nfile;
    nfile << bestUb;
    string temp = nfile.str();
    env0.sendbuff[6] = temp.size();
    int i = 7;
    for (string::iterator it = temp.begin(); it != temp.end(); ++it) {
        env0.sendbuff[i] = *it - '0';
        i++;
    }
    nfile.clear();
    for (map<int, Value>::iterator it = bestSolution.begin();
         it != bestSolution.end(); ++it) {
        env0.sendbuff[i] = it->second;
        i++;
    }
}

//Message to solution
void CooperativeParallelDGVNS::MsgToSol(
    MPIEnv& env0, int nov, uint& cluster, uint& numberclu, int& k,
    int& kmax, Cost& bestUb, int& sec, int& msec, map<int, Value>& bestSolution)
{
    cluster = env0.recvbuff[0];
    k = env0.recvbuff[1];
    kmax = env0.recvbuff[2];
    sec = env0.recvbuff[3];
    msec = env0.recvbuff[4];
    numberclu = env0.recvbuff[5];
    uint size = env0.recvbuff[6];
    bestUb = env0.recvbuff[7];
    int j = 8;
    for (uint it = 1; it < size; it++) {
        bestUb = (bestUb * 10) + env0.recvbuff[j];
        j++;
    }
    for (int i = 0; i < nov; i++) {
        bestSolution[i] = env0.recvbuff[j];
        j++;
    }
}

//---------------- Class Definition --------------------------//

bool CooperativeParallelDGVNS::solve(bool first)
{
    mysrand(abs(ToulBar2::seed) + env0.myrank);

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
        if (env0.myrank == 0) {
            if (lastUb < MAX_COST)
                wcsp->setSolution(lastUb, &lastSolution);
            endSolve(lastUb < MAX_COST, lastUb, true);
        }
        /* Shut down MPI */
        MPI_Finalize();
        return (lastUb < MAX_COST);
    }

    MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */

    env0.buffsize = (int)wcsp->numberOfVariables() + 107; // 3 : cluster + k + cost, second time, msecond time, localtime,the rest is the size of solution
    env0.sendbuff = new int[env0.buffsize];
    env0.recvbuff = new int[env0.buffsize];
    if (env0.myrank == 0) {
        if (ToulBar2::verbose)
            cout << "Run CPDGVNS method." << endl;
        master();
        if (bestUb < MAX_COST)
            wcsp->setSolution(bestUb, &bestSolution);
        endSolve(bestUb < MAX_COST, bestUb, false); //TODO: return complete=true if known by the master
    } else {
        ToulBar2::showSolutions = false;
        slave();
    }

    if (ToulBar2::vnsOutput)
        ToulBar2::vnsOutput << "Search end"
                            << " " << cpuTime() << endl;

    /* Shut down MPI */
    //    MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
    MPI_Finalize();

    return (bestUb < MAX_COST);
}

//----------------------- Model Master/Slave -----------------------//
void CooperativeParallelDGVNS::master()
{
    // Structure de voisinage bas�e sur la notion des clusters -  ParallelRandomClusterChoice - et l'initialisation du file - file contenant les clusters pour chaque processus -
    ParallelRandomClusterChoice* h = new ParallelRandomClusterChoice();
    h->init(wcsp, this);

    // MPI data
    int rank;
    MPI_Status status;

    // verify the number of processes and number of clusters
    if (env0.ntasks > h->getSize()) {
        env0.processes = h->getSize();
    } else {
        env0.processes = env0.ntasks;
    }
    if (ToulBar2::vnsOutput)
        ToulBar2::vnsOutput << "#param " << ToulBar2::vnsNeighborVarHeur << " " << ToulBar2::vnsKmin << " " << ToulBar2::vnsKmax << " " << ToulBar2::lds << " " << ToulBar2::vnsInitSol << " " << env0.processes << " " << env0.ntasks << " " << h->getSize() << endl;
    // Generation of initial Solution
    bool complete = false;
    bestSolution.clear();
    bestUb = generateInitSolution(ToulBar2::vnsInitSol, bestSolution, complete);
    BestTimeS = (int)cpuTime();
    BestTimeMS = (int)(((Long)cpuTime() * 1000LL) % 1000LL);
    // Get all clusters from the tree decomposition of constraint graph
    file = h->getClustersIndex();

    /* Seed the slaves; send one unit of work (initial solution) to each slave. */
    for (rank = 1; rank < env0.processes; ++rank) {
        /* choice one free cluster to send it to the processe */
        uint cluster = getCluster();
        int adjcluster = 0;
        int kmax = h->getClustersSize(cluster, adjcluster);
        while (kmax < ToulBar2::vnsKmin) {
            adjcluster++;
            kmax = h->getClustersSize(cluster, adjcluster);
        }
        //cout << kmax << " " << adjcluster << endl ;
        /* Convert initial solution with cluster and kinit parameters in buffer, for each slave process */
        SolToMsg(env0, cluster, adjcluster, ToulBar2::vnsKmin, ToulBar2::vnsKmax, bestUb, BestTimeS,
            BestTimeMS, bestSolution);

        /* Send Initial Solution to each process */
        MPI_Send(&env0.sendbuff[0], /* message buffer */
            env0.buffsize, /* buffer size */
            MPI_INT, /* data item is an integer */
            rank, /* destination process rank */
            WORKTAG, /* user chosen message tag */
            MPI_COMM_WORLD); /* default communicator */
    }

    /* Loop over getting new Best Solutions */
    uint finished = 0;
    uint worker = env0.processes - 1;
    map<int, Value> slastSolution;
    while (finished < worker) {

        /* Receive result (new better solution) from a slave */
        MPI_Recv(&env0.recvbuff[0], /* message buffer */
            env0.buffsize, /* size of data */
            MPI_INT, /* data item is an integer */
            MPI_ANY_SOURCE, /* receive from any sender */
            MPI_ANY_TAG, /* any type of message */
            MPI_COMM_WORLD, /* default communicator */
            &status); /* info about the received message */
        /* setting up the best solution in memory, and checkout if we continue or not */
        // Parsing the received buffer
        uint scluster = 0;
        int sk = 0;
        int skmax;
        Cost sbestUb = 0;
        int sBestTimeS = 0;
        int sBestTimeMS = 0;
        uint snumberclu;

        MsgToSol(env0, wcsp->numberOfVariables(), scluster, snumberclu, sk,
            skmax, sbestUb, sBestTimeS, sBestTimeMS, slastSolution);

        // getting the cluster back
        file.push_back(scluster);
        scluster = getCluster();
        //Updating
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
            SolToMsg(env0, scluster, snumberclu, ToulBar2::vnsKmin, skmax, bestUb,
                BestTimeS, BestTimeMS, bestSolution);
            /* Send the best solution to the received slave, for next the search */
            MPI_Send(&env0.sendbuff[0], /* message buffer */
                env0.buffsize, /* size of data */
                MPI_INT, /* data item is an integer */
                status.MPI_SOURCE, /* to who we just received from */
                WORKTAG, /* user chosen message tag */
                MPI_COMM_WORLD); /* default communicator */
        } else {
            finished++;
            MPI_Send(0, 0, MPI_INT, status.MPI_SOURCE, DIETAG, MPI_COMM_WORLD);
        }
    }
    MPI_Request requests[env0.ntasks];
    for (rank = 1; rank < env0.ntasks; ++rank) {
        //printf("Send finish empty msg to finish with %d\n",rank);
        MPI_Isend(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD, &requests[rank]);
    }
}

void CooperativeParallelDGVNS::slave()
{

    MPI_Status status;

    // Structure de voisinage basée sur la notion des clusters
    //ParallelRandomClusterChoice* h = NeighborhoodStructure::NeighborhoodStructureFactory(VariableHeuristic(hname), static_cast<WCSP*>(wcsp), this);
    ParallelRandomClusterChoice* h = new ParallelRandomClusterChoice();
    h->init(wcsp, this);
    if (env0.ntasks > h->getSize()) {
        env0.processes = h->getSize();
    } else {
        env0.processes = env0.ntasks;
    }

    // initializing the timer begin for the hole search of the slave process, in seconds
    double btime = cpuTime();

    wcsp->setUb(MAX_COST);
    while (true) {
        /* Receive a message from the master */
        MPI_Recv(&env0.recvbuff[0], /* message buffer */
            env0.buffsize, /* one data item */
            MPI_INT, /* of type integer */
            0, /* receive from master */
            MPI_ANY_TAG, /* any type of message */
            MPI_COMM_WORLD, /* default communicator */
            &status); /* info about the received message */
        //cout << env0.myrank <<" slave begin" << endl ;
        /* Check the tag of the received message. */
        if (status.MPI_TAG == DIETAG) {
            return;
        }

        /* launch Vns/Lds+cp */
        VnsLdsCP(env0, btime, h);

        /* Send the result back */
        MPI_Request request;
        MPI_Isend(&env0.sendbuff[0], env0.buffsize, MPI_INT, 0, 0,
            MPI_COMM_WORLD, &request);

        //cout << env0.myrank <<" slave end" << endl ;
    }
}

//------------------- Base functions -----------------------//

void CooperativeParallelDGVNS::VnsLdsCP(MPIEnv& env0, double btime, ParallelRandomClusterChoice* h)
{
    uint currentcluster = 0;
    uint numberclu = 0;
    int kinit = ToulBar2::vnsKmin;
    int kmax = ToulBar2::vnsKmax;

    MsgToSol(env0, wcsp->numberOfVariables(), currentcluster, numberclu, kinit,
        kmax, bestUb, BestTimeS, BestTimeMS, bestSolution);
    for (map<int, Value>::iterator it = bestSolution.begin();
         it != bestSolution.end(); ++it)
        lastSolution[(*it).first] = (*it).second;
    if (ToulBar2::verbose >= 1)
        cout << "VNS :: Initial Solution cost " << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->Cost2ADCost(bestUb) << std::setprecision(DECIMAL_POINT) << endl;
    lastUb = bestUb;
    if (ToulBar2::verbose >= 5)
        cout << "VNS :: Initial Solution" << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->Cost2ADCost(bestUb) << std::setprecision(DECIMAL_POINT) << endl;

    //vns/lds+cp
    int k = kinit;
    //cout << "taille maximal du cluster "<< currentcluster<< " "<< numberclu<< " "<< kmax << endl ;
    //cout << env0.myrank <<" slave 1" << endl ;
    //cout << k <<"<="<< kmax <<"&&"<< k <<"<="<< wcsp->numberOfVariables() << "&&"<< ToulBar2::vns_optimum <<"<"  << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->Cost2ADCost(bestUb) << std::setprecision(DECIMAL_POINT) <<"&&"<< (cpuTime()-lbtime)<<endl ;
    for (; k <= kmax && k <= unassignedVars->getSize() && ToulBar2::vnsOptimum < bestUb;) {
        //neighborhood and partial instantiation
        //cout <<"neighborhood"<< " "<<currentcluster<< " "<< numberclu << " " << k << endl ;
        set<int> neighborhood = h->SlaveGetNeighborhood(currentcluster, numberclu, k);

        if (ToulBar2::verbose >= 1) {
            cout << env0.myrank << ": LDS " << ToulBar2::lds << " Neighborhood " << k;
            for (set<int>::iterator it = neighborhood.begin(); it != neighborhood.end(); it++)
                cout << " " << *it;
            cout << endl;
        }
        vector<int> variables;
        variables.reserve(unassignedVars->getSize());
        vector<int> values;
        values.reserve(unassignedVars->getSize());
        for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
            int v = *iter;
            if (neighborhood.find(v) == neighborhood.end()) {
                variables.push_back(v);
                values.push_back(lastSolution[v]);
            }
        }
        //repair
        if (ToulBar2::lds)
            repair_recursiveSolve(ToulBar2::lds, variables, values, bestUb);
        else
            repair_recursiveSolve(variables, values, bestUb);
        //updating
        //cout <<"updating "<< k<< endl ;
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
        }
    }
    //cout << env0.myrank <<" slave 2" << endl ;
    // return the best solution found
    SolToMsg(env0, currentcluster, numberclu, k, kmax, bestUb, BestTimeS,
        BestTimeMS, bestSolution);
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

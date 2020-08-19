/*
 * tb2rpdgvns.cpp
 *
 *      Author: Abdelkader Ouali
 *      PhD Student: LITIO, University of Oran ; GREYC, University of Caen.
 */

#include "tb2rpdgvns.hpp"
#include "core/tb2wcsp.hpp"
#ifdef OPENMPI

// Conversion Tools
// Solution to Message
void ReplicatedParallelDGVNS::SolToMsg(
    MPIEnv& env0, int cluster, int k, int discrepancy, Cost bestUb,
    map<int, Value>& bestSolution)
{
    int i = 0;
    env0.sendbuff[i] = cluster;
    i++;
    env0.sendbuff[i] = k;
    i++;
    env0.sendbuff[i] = discrepancy;
    i++;
    stringstream nfile;
    nfile << bestUb;
    string temp = nfile.str();
    env0.sendbuff[i] = temp.size();
    i++;
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

void ReplicatedParallelDGVNS::SolToMsg2(
    MPIEnv& env0, Cost bestUb, map<int, Value>& bestSolution)
{
    int i = 0;
    stringstream ss_bestUb;
    ss_bestUb << bestUb;
    string temp = ss_bestUb.str();
    env0.sendbuff[i] = temp.size();
    i++;
    for (string::iterator it = temp.begin(); it != temp.end(); ++it) {
        env0.sendbuff[i] = *it - '0';
        i++;
    }
    ss_bestUb.clear();
    for (map<int, Value>::iterator it = bestSolution.begin();
         it != bestSolution.end(); ++it) {
        env0.sendbuff[i] = it->second;
        i++;
    }
}

//Message to solution
void ReplicatedParallelDGVNS::MsgToSol(
    MPIEnv& env0, int nov, int& cluster, int& k, int& discrepancy,
    Cost& bestUb, map<int, Value>& bestSolution)
{
    int j = 0;
    cluster = env0.recvbuff[j];
    j++;
    k = env0.recvbuff[j];
    j++;
    discrepancy = env0.recvbuff[j];
    j++;
    uint size = env0.recvbuff[j];
    j++;
    bestUb = env0.recvbuff[j];
    j++;
    for (uint it = 1; it < size; it++) {
        bestUb = (bestUb * 10) + env0.recvbuff[j];
        j++;
    }
    for (int i = 0; i < nov; i++) {
        bestSolution[i] = env0.recvbuff[j];
        j++;
    }
}

void ReplicatedParallelDGVNS::MsgToSol2(
    MPIEnv& env0, int nov, Cost& bestUb, map<int, Value>& bestSolution)
{
    int j = 0;
    uint size = env0.recvbuff[0];
    j++;
    bestUb = env0.recvbuff[1];
    j++;
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
void timeOut()
{
    MPI_Abort(MPI_COMM_WORLD, 0);
    //    MPI_Finalize();
    //    exit(0);
}

bool ReplicatedParallelDGVNS::solve(bool first)
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
        if (env0.myrank == 0) {
            if (lastUb < MAX_COST)
                wcsp->setSolution(lastUb, &lastSolution);
            endSolve(lastUb < MAX_COST, lastUb, true);
        }
        /* Shut down MPI */
        MPI_Finalize();
        return (lastUb < MAX_COST);
    }

    assert((int)wcsp->numberOfUnassignedVariables() == unassignedVars->getSize());
    ToulBar2::vnsLDSmax = min(ToulBar2::vnsLDSmax, (int)wcsp->getDomainSizeSum() - (int)wcsp->numberOfUnassignedVariables());
    ToulBar2::vnsLDSmin = min(ToulBar2::vnsLDSmin, ToulBar2::vnsLDSmax);
    ToulBar2::vnsKmax = min(ToulBar2::vnsKmax, (int)wcsp->numberOfUnassignedVariables());
    ToulBar2::vnsKmin = min(ToulBar2::vnsKmin, ToulBar2::vnsKmax);
    assert(ToulBar2::vnsLDSmin >= 0);
    assert(ToulBar2::vnsLDSmax >= 0);
    assert(ToulBar2::vnsLDSmin <= ToulBar2::vnsLDSmax);
    assert(ToulBar2::vnsKmin >= 0);
    assert(ToulBar2::vnsKmax >= 0);
    assert(ToulBar2::vnsKmin <= ToulBar2::vnsKmax);

    bool complete = false;

    mysrand(abs(ToulBar2::seed) + env0.myrank);

    MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
    startTime = MPI_Wtime();
    ToulBar2::timeOut = timeOut;

    env0.buffsize = (int)wcsp->numberOfVariables() + 108; // 3 : cluster + k + cost, second time, msecond time, localtime,the rest is the size of solution
    env0.sendbuff = new int[env0.buffsize];
    env0.recvbuff = new int[env0.buffsize];
    if (env0.myrank == 0) {
        if (!ToulBar2::vnsParallelSync) {
            complete = radgvns();
            endSolve(bestUb < MAX_COST, bestUb, complete);
        } else {
            complete = rsdgvns();
            endSolve(bestUb < MAX_COST, bestUb, complete);
        }
    } else {
        ToulBar2::showSolutions = false;
        complete = slave();
    }

    if (ToulBar2::vnsOutput)
        ToulBar2::vnsOutput << "Search end"
                            << " " << cpuTime() << endl;
    if (ToulBar2::vnsOutput)
        ToulBar2::vnsOutput.close(); // close the output file

    /* Shut down MPI */
    MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
    double elapsedTime = MPI_Wtime() - startTime;
    double elapsedCPUTime = cpuTime() - ToulBar2::startCpuTime;
    double totalCPUTime = 0.;
    MPI_Reduce(&elapsedCPUTime, &totalCPUTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (ToulBar2::verbose >= 0 && env0.myrank == 0) { /* use time on master node */
        cout << "Total CPU time = " << totalCPUTime << " seconds" << endl;
        cout << "Solving real-time = " << elapsedTime << " seconds (not including preprocessing time)" << endl;
    }
    MPI_Finalize();

    return (bestUb < MAX_COST);
}

//----------------------- Model Master/Slave -----------------------//
bool ReplicatedParallelDGVNS::radgvns()
{
    if (ToulBar2::verbose >= 1)
        cout << " RADGVNS kinit=" << ToulBar2::vnsKmin << " kmax=" << ToulBar2::vnsKmax
             << " discrepancyinit=" << ToulBar2::vnsLDSmin << " discrepancymax=" << ToulBar2::vnsLDSmax << " neighbor change if improved=" << ToulBar2::vnsNeighborChange
             << " neighbor size synchronization=" << ToulBar2::vnsNeighborSizeSync << " limit number processes=" << ToulBar2::vnsParallelLimit << endl;

    // cluster tree initialized BEFORE generating initial solution
    ParallelRandomClusterChoice* h = new ParallelRandomClusterChoice();
    h->init(wcsp, this);
    if (ToulBar2::verbose >= 0 && ToulBar2::vnsNeighborVarHeur == MASTERCLUSTERRAND && ((ClustersNeighborhoodStructure*)h)->getSize() > 1) {
        ClustersNeighborhoodStructure* ch = (ClustersNeighborhoodStructure*)h;
        if (ToulBar2::verbose >= 1 || ToulBar2::debug)
            ch->printClusters(cout);
        cout << "Problem decomposition in " << ch->getSize() << " clusters with size distribution: min: " << ch->getMinClusterSize() << " median: " << ch->getMedianClusterSize() << " mean: " << ch->getMeanClusterSize() << " max: " << ch->getMaxClusterSize() << endl;
    }

    // MPI data
    MPI_Status status;

    // verify the number of processes and number of clusters
    env0.processes = env0.ntasks - 1;
    //    cout << "number of processes=" << env0.processes << endl;
    // verify the number of processes and number of clusters
    int npr = 0;
    if (ToulBar2::vnsParallelLimit) {
        if (env0.ntasks > h->getSize()) {
            npr = h->getSize();
        } else {
            npr = env0.processes;
        }
    } else {
        npr = env0.processes;
    }

    if (ToulBar2::vnsOutput)
        ToulBar2::vnsOutput << "#param " << ToulBar2::vnsNeighborVarHeur << " " << ToulBar2::vnsKmin << " " << ToulBar2::vnsKmax << " " << ToulBar2::vnsLDSmin << " " << ToulBar2::vnsInitSol << " " << env0.processes << " " << env0.ntasks << " " << h->getSize() << endl;

    // Generation of initial Solution
    file = h->getClustersIndex();
    //    clusterKmax = vector<bool>(file.size(), false);
    //    cout << "reading Time=" <<  cpuTime() << endl;
    //    cout << "initial solution mode : " << ToulBar2::vnsInitSol << endl;
    bool complete = false;
    bestSolution.clear();
    bestUb = generateInitSolution(ToulBar2::vnsInitSol, bestSolution, complete);
    DumpBestSol(ToulBar2::vnsInitSol != LS_INIT_DFBB && ToulBar2::vnsInitSol < LS_INIT_LDS0);

    // Get all clusters from the tree decomposition of constraint graph
    int c = 0;
    if (!complete && bestUb > ToulBar2::vnsOptimum) {
        for (int p = 1; p < npr + 1; ++p) {
            pr pr_p = pr();
            pr_p.cl = c;
            pr_p.k = ToulBar2::vnsKmin;
            pr_p.lds = ToulBar2::vnsLDSmin;
            pr_p.synch = false;
            vecPR.push_back(pr_p);
            SolToMsg(env0, pr_p.cl, pr_p.k, pr_p.lds, bestUb, bestSolution);
            MPI_Send(&env0.sendbuff[0], env0.buffsize, MPI_INT, p, WORKTAG, MPI_COMM_WORLD);
            c = (c + 1) % file.size();
        }
    }
    Cost pbestUb = MAX_COST;
    map<int, Value> pbestSolution;
    while (npr && !complete && bestUb > ToulBar2::vnsOptimum) {
        MPI_Recv(&env0.recvbuff[0], env0.buffsize, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        int pindex = status.MPI_SOURCE - 1;
        complete = (status.MPI_TAG == WORKTAG);
        if ((vecPR[pindex].lds >= ToulBar2::vnsLDSmax || ToulBar2::restart == 1) && vecPR[pindex].k >= ToulBar2::vnsKmax)
            npr = 0;
        MsgToSol2(env0, wcsp->numberOfVariables(), pbestUb, pbestSolution);
        NeighborhoodChange(ToulBar2::vnsNeighborChange, pindex, c, ToulBar2::vnsKmin, ((ClustersNeighborhoodStructure*)h)->getMaxClusterSize() + ((ClustersNeighborhoodStructure*)h)->getSize() - 1, ToulBar2::vnsKmax, ToulBar2::vnsLDSmin, ToulBar2::vnsLDSmax, ToulBar2::vnsNeighborSizeSync, pbestUb, pbestSolution);
        //        if (ToulBar2::restart==1 && find(clusterKmax.begin(), clusterKmax.end(), false) == clusterKmax.end()) npr = 0;
        if (!complete && bestUb > ToulBar2::vnsOptimum) {
            SolToMsg(env0, vecPR[pindex].cl, vecPR[pindex].k, vecPR[pindex].lds, bestUb, bestSolution);
            MPI_Send(&env0.sendbuff[0], env0.buffsize, MPI_INT, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);
        }
    }

    MPI_Request requests[env0.ntasks];
    for (int p = 1; p < env0.ntasks; ++p) {
        //        MPI_Isend(0, 0, MPI_INT, p, DIETAG, MPI_COMM_WORLD);

        MPI_Isend(NULL, 0, MPI_INT, p, DIETAG, MPI_COMM_WORLD, &requests[p]);
    }

    return complete && (bestUb < MAX_COST);
}

bool ReplicatedParallelDGVNS::rsdgvns()
{
    if (ToulBar2::verbose >= 1)
        cout << " RSDGVNS kinit=" << ToulBar2::vnsKmin << " kmax=" << ToulBar2::vnsKmax
             << " discrepancyinit=" << ToulBar2::vnsLDSmin << " discrepancymax=" << ToulBar2::vnsLDSmax << " neighbor change if improved=" << ToulBar2::vnsNeighborChange
             << " neighbor size synchronization=" << ToulBar2::vnsNeighborSizeSync << " limit number processes=" << ToulBar2::vnsParallelLimit << endl;

    // cluster tree initialized BEFORE generating initial solution
    ParallelRandomClusterChoice* h = new ParallelRandomClusterChoice();
    h->init(wcsp, this);
    if (ToulBar2::verbose >= 0 && ToulBar2::vnsNeighborVarHeur == MASTERCLUSTERRAND && ((ClustersNeighborhoodStructure*)h)->getSize() > 1) {
        ClustersNeighborhoodStructure* ch = (ClustersNeighborhoodStructure*)h;
        if (ToulBar2::verbose >= 1 || ToulBar2::debug)
            ch->printClusters(cout);
        cout << "Problem decomposition in " << ch->getSize() << " clusters with size distribution: min: " << ch->getMinClusterSize() << " median: " << ch->getMedianClusterSize() << " mean: " << ch->getMeanClusterSize() << " max: " << ch->getMaxClusterSize() << endl;
    }

    // MPI data
    MPI_Status status;

    // Initialization
    map<int, Value> bestInterSolution;
    Cost bestInterUb = MAX_COST;

    // verify the number of processes and number of clusters
    env0.processes = env0.ntasks - 1;
    cout << "number of processes=" << env0.processes << endl;
    // verify the number of processes and number of clusters
    int npr = 0;
    if (ToulBar2::vnsParallelLimit) {
        if (env0.ntasks > h->getSize()) {
            npr = h->getSize();
        } else {
            npr = env0.processes;
        }
    } else {
        npr = env0.processes;
    }

    if (ToulBar2::vnsOutput)
        ToulBar2::vnsOutput << "#param " << ToulBar2::vnsNeighborVarHeur << " " << ToulBar2::vnsKmin << " " << ToulBar2::vnsKmax << " " << ToulBar2::vnsLDSmin << " " << ToulBar2::vnsInitSol << " " << env0.processes << " " << env0.ntasks << " " << h->getSize() << endl;

    // Generation of initial Solution
    file = h->getClustersIndex();
    //    clusterKmax = vector<bool>(file.size(), false);
    if (ToulBar2::vnsOutput)
        ToulBar2::vnsOutput << "cslsize " << file.size() << " npr= " << npr << endl;

    bool complete = false;
    bestSolution.clear();
    bestUb = generateInitSolution(ToulBar2::vnsInitSol, bestSolution, complete);
    for (map<int, Value>::iterator it = bestSolution.begin();
         it != bestSolution.end(); ++it)
        bestInterSolution[(*it).first] = (*it).second;
    bestInterUb = bestUb;
    DumpBestSol(ToulBar2::vnsInitSol != LS_INIT_DFBB && ToulBar2::vnsInitSol < LS_INIT_LDS0);

    // Get all clusters from the tree decomposition of constraint graph
    int c = 0;
    bool stop = false;
    Long nbRestart = 1;
    Long restart = 1;
    int lds = ToulBar2::vnsLDSmin;
    while (npr && !stop && !complete && bestUb > ToulBar2::vnsOptimum) {
        if (ToulBar2::verbose >= 0 && ToulBar2::restart > 1 && ToulBar2::lds)
            cout << "****** Restart " << nbRestart << " with " << lds << " discrepancies and UB=" << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->Cost2ADCost(bestUb) << std::setprecision(DECIMAL_POINT) << " ****** (" << nbNodes << " nodes)" << endl;
        Long rank = 1;
        int k = ToulBar2::vnsKmin;
        while (npr && !complete && k <= ToulBar2::vnsKmax && bestUb > ToulBar2::vnsOptimum) {

            for (int p = 1; p < npr + 1; ++p) {
                SolToMsg(env0, c, k, lds, bestUb, bestSolution);
                MPI_Send(&env0.sendbuff[0], env0.buffsize, MPI_INT, p, WORKTAG, MPI_COMM_WORLD);
                c = (c + 1) % file.size();
            }

            int finished = 0;
            //S"
            Cost pBestUb;
            map<int, Value> pBestSolution;

            while (finished < npr) {
                MPI_Recv(&env0.recvbuff[0], env0.buffsize, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                MsgToSol2(env0, wcsp->numberOfVariables(), pBestUb, pBestSolution);
                if (pBestUb < bestInterUb) {
                    bestInterUb = pBestUb;
                    for (int v = 0; v < (int)wcsp->numberOfVariables(); v++) {
                        bestInterSolution[v] = pBestSolution[v];
                    }
                }
                finished++;
            }

            if (bestInterUb < bestUb) {
                rank = 1;
                k = ToulBar2::vnsKmin;
                restart = 1;
                lds = ToulBar2::vnsLDSmin;
                bestUb = bestInterUb;
                for (int v = 0; v < (int)wcsp->numberOfVariables(); v++) {
                    bestSolution[v] = bestInterSolution[v];
                }
                DumpBestSol();
            } else {
                rank++;
                if (k < ToulBar2::vnsKmax) {
                    switch (ToulBar2::vnsKinc) {
                    case VNS_ADD1:
                        k++;
                        break;
                    case VNS_MULT2:
                        k *= 2;
                        break;
                    case VNS_LUBY:
                        k = ToulBar2::vnsKmin * (int)luby(rank);
                        break;
                    case VNS_ADD1JUMP:
                        if (ToulBar2::vnsNeighborVarHeur == RANDOMVAR || k < ((ClustersNeighborhoodStructure*)h)->getMaxClusterSize() + ((ClustersNeighborhoodStructure*)h)->getSize() - 1)
                            k++;
                        else
                            k = ToulBar2::vnsKmax;
                        break;
                    default:
                        cerr << "Unknown neighborhood size increment strategy inside VNS (see option -kinc)!" << endl;
                        exit(EXIT_FAILURE);
                    }
                    k = min(k, ToulBar2::vnsKmax);
                } else
                    k++;
            }
        }
        if (!complete && bestUb > ToulBar2::vnsOptimum) {
            nbRestart++;
            restart++;
            if (nbRestart <= ToulBar2::restart && lds < ToulBar2::vnsLDSmax) {
                if (ToulBar2::lds) {
                    switch (ToulBar2::vnsLDSinc) {
                    case VNS_ADD1:
                        lds++;
                        break;
                    case VNS_MULT2:
                        lds *= 2;
                        break;
                    case VNS_LUBY:
                        lds = ToulBar2::vnsLDSmin * (int)luby(restart);
                        break;
                    default:
                        cerr << "Unknown LDS increment strategy inside VNS (see option -linc)!" << endl;
                        exit(EXIT_FAILURE);
                    }
                    lds = min(lds, ToulBar2::vnsLDSmax);
                }
            } else
                stop = true;
        }
    }

    MPI_Request requests[env0.ntasks];
    for (int p = 1; p < env0.ntasks; ++p) {
        //        MPI_Isend(0, 0, MPI_INT, p, DIETAG, MPI_COMM_WORLD);
        MPI_Isend(NULL, 0, MPI_INT, p, DIETAG, MPI_COMM_WORLD, &requests[p]);
    }

    return complete && (bestUb < MAX_COST);
}

bool ReplicatedParallelDGVNS::slave()
{
    // cluster tree initialized BEFORE generating initial solution
    ParallelRandomClusterChoice* h = new ParallelRandomClusterChoice();
    h->init(wcsp, this);

    // MPI data
    MPI_Status status;

    // loop for getting new initial solution to initialize the search process
    bool complete = false;
    while (true) {
        MPI_Recv(&env0.recvbuff[0], env0.buffsize, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == DIETAG)
            return complete; // warning! always wait for DIETAG before stopping
        /* launch Vns/Lds+cp */
        complete = complete || VnsLdsCP(env0, h);
        /* Send the result back */
        MPI_Request request;
        MPI_Isend(&env0.sendbuff[0], env0.buffsize, MPI_INT, 0, complete, MPI_COMM_WORLD, &request);
    }
    return complete;
}

//------------------- Base functions -----------------------//

bool ReplicatedParallelDGVNS::VnsLdsCP(MPIEnv& env0, ParallelRandomClusterChoice* h)
{
    int cluster, k, discrepancy;
    MsgToSol(env0, wcsp->numberOfVariables(), cluster, k, discrepancy, bestUb,
        bestSolution);
    for (map<int, Value>::iterator it = bestSolution.begin();
         it != bestSolution.end(); ++it)
        lastSolution[(*it).first] = (*it).second;
    lastUb = bestUb;
    //vns/lds+cp
    set<int> neighborhood = h->SlaveGetNeighborhood(cluster, k); // based shuffle
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
    if (ToulBar2::verbose >= 1) {
        cout << env0.myrank << ": LDS " << discrepancy << " Neighborhood " << k << ": ";
        for (set<int>::iterator it = neighborhood.begin(); it != neighborhood.end(); it++)
            cout << " " << *it;
        cout << endl;
    }

    //repair
    bool complete = false;
    ToulBar2::vnsKcur = k;
    ToulBar2::vnsLDScur = (ToulBar2::lds) ? discrepancy : -1;
    if (ToulBar2::lds)
        complete = repair_recursiveSolve(discrepancy, variables, values, bestUb);
    else
        complete = repair_recursiveSolve(variables, values, bestUb);
    SolToMsg2(env0, lastUb, lastSolution);
    return complete;
}

void ReplicatedParallelDGVNS::NeighborhoodChange(
    int strategy, int p, int& c, int kinit, int kjump, int kmax, int ldsmin, int ldsmax,
    bool synch, Cost pBestUb, map<int, Value>& pBestSolution)
{
    switch (strategy) {
    case 0:
        ChangeClusterAlways(p, c, kinit, kjump, kmax, ldsmin, ldsmax, synch, pBestUb,
            pBestSolution);
        break;
    case 1:
        ChangeClusterWhenNotImproved(p, c, kinit, kjump, kmax, ldsmin, ldsmax, synch,
            pBestUb, pBestSolution);
        break;
    default:
        ChangeClusterAlways(p, c, kinit, kjump, kmax, ldsmin, ldsmax, synch, pBestUb,
            pBestSolution);
        break;
    }
}

void ReplicatedParallelDGVNS::ChangeClusterAlways(
    int p, int& c, int kinit, int kjump, int kmax, int ldsmin, int ldsmax, bool synch,
    Cost pBestUb, map<int, Value>& pBestSolution)
{
    bool improved = (pBestUb < bestUb);
    if (pBestUb <= bestUb) {
        bestUb = pBestUb;
        for (int v = 0; v < (int)wcsp->numberOfVariables(); v++) {
            bestSolution[v] = pBestSolution[v];
        }
        //		ToulBar2::vnsOutput << "# ------------Change A------------------------------" << endl;
        //		ToulBar2::vnsOutput << "cluster=" << c <<" process="<< p << " k="<< vecPR[p].k << " Time " << BestTime << endl;
        if (improved)
            DumpBestSol();
        for (int i = 0; synch && i < (int)vecPR.size(); i++) {
            vecPR[i].synch = true;
        }
        vecPR[p].synch = false;
        vecPR[p].cl = c;
        vecPR[p].k = kinit;
        vecPR[p].lds = ldsmin;
    } else {
        vecPR[p].cl = c;
        if (synch && vecPR[p].synch == true) {
            vecPR[p].k = kinit;
            vecPR[p].lds = ldsmin;
            vecPR[p].synch = false;
        } else {
            if (vecPR[p].k < kmax) {
                switch (ToulBar2::vnsKinc) {
                case VNS_ADD1:
                case VNS_LUBY:
                    vecPR[p].k++;
                    break;
                case VNS_MULT2:
                    vecPR[p].k *= 2;
                    break;
                case VNS_ADD1JUMP:
                    if (ToulBar2::vnsNeighborVarHeur == RANDOMVAR || vecPR[p].k < kjump)
                        vecPR[p].k++;
                    else
                        vecPR[p].k = kmax;
                    break;
                default:
                    cerr << "Unknown neighborhood size increment strategy inside VNS (see option -kinc)!" << endl;
                    exit(EXIT_FAILURE);
                }
                vecPR[p].k = min(vecPR[p].k, kmax);
            } else if (ToulBar2::restart > 1) { // Warning, unbounded number of restarts...
                vecPR[p].k = kinit;
                if (ToulBar2::lds) {
                    switch (ToulBar2::vnsLDSinc) {
                    case VNS_ADD1:
                    case VNS_LUBY:
                        vecPR[p].lds++;
                        break;
                    case VNS_MULT2:
                        vecPR[p].lds *= 2;
                        break;
                    default:
                        cerr << "Unknown LDS increment strategy inside VNS (see option -linc)!" << endl;
                        exit(EXIT_FAILURE);
                    }
                    vecPR[p].lds = min(vecPR[p].lds, ToulBar2::vnsLDSmax);
                }
                //            } else {
                //                clusterKmax[c] = true;
            }
        }
    }
    c = (c + 1) % file.size();
}

void ReplicatedParallelDGVNS::ChangeClusterWhenNotImproved(
    int p, int& c, int kinit, int kjump, int kmax, int ldsmin, int ldsmax, bool synch,
    Cost pBestUb, map<int, Value>& pBestSolution)
{
    bool improved = (pBestUb < bestUb);
    if (pBestUb <= bestUb) {
        bestUb = pBestUb;
        for (int v = 0; v < (int)wcsp->numberOfVariables(); v++) {
            bestSolution[v] = pBestSolution[v];
        }
        //		ToulBar2::vnsOutput << "# -------NOIMP-----------------------------------" << endl;
        //		ToulBar2::vnsOutput << "cluster="<< c << " process="<< p << " k="<< vecPR[p].k << " Time=" << BestTime <<endl;
        if (improved)
            DumpBestSol();
        for (int i = 0; synch && i < (int)vecPR.size(); i++) {
            vecPR[i].synch = true;
        }
        vecPR[p].synch = false;
        // keep the current cluster
        vecPR[p].k = kinit;
        vecPR[p].lds = ldsmin;
    } else {
        vecPR[p].cl = c; // change cluster when no improvement is considered
        if (synch && vecPR[p].synch == true) {
            vecPR[p].k = kinit;
            vecPR[p].lds = ldsmin;
            vecPR[p].synch = false;
        } else {
            if (vecPR[p].k < kmax) {
                switch (ToulBar2::vnsKinc) {
                case VNS_ADD1:
                case VNS_LUBY:
                    vecPR[p].k++;
                    break;
                case VNS_MULT2:
                    vecPR[p].k *= 2;
                    break;
                case VNS_ADD1JUMP:
                    if (ToulBar2::vnsNeighborVarHeur == RANDOMVAR || vecPR[p].k < kjump)
                        vecPR[p].k++;
                    else
                        vecPR[p].k = kmax;
                    break;
                default:
                    cerr << "Unknown neighborhood size increment strategy inside VNS (see option -kinc)!" << endl;
                    exit(EXIT_FAILURE);
                }
                vecPR[p].k = min(vecPR[p].k, kmax);
            } else if (ToulBar2::restart > 1) { // Warning, unbounded number of restarts...
                vecPR[p].k = kinit;
                if (ToulBar2::lds) {
                    switch (ToulBar2::vnsLDSinc) {
                    case VNS_ADD1:
                    case VNS_LUBY:
                        vecPR[p].lds++;
                        break;
                    case VNS_MULT2:
                        vecPR[p].lds *= 2;
                        break;
                    default:
                        cerr << "Unknown LDS increment strategy inside VNS (see option -linc)!" << endl;
                        exit(EXIT_FAILURE);
                    }
                    vecPR[p].lds = min(vecPR[p].lds, ToulBar2::vnsLDSmax);
                }
                //            } else {
                //                clusterKmax[c] = true;
            }
        }
        c = (c + 1) % file.size(); // change cluster when no improvement is considered
    }
}

void ReplicatedParallelDGVNS::DumpBestSol(bool improved)
{
    // Save
    wcsp->setSolution(bestUb, &bestSolution);
    if (ToulBar2::vnsOutput) {
        //cout << "# ------------------------------------------" << endl;
        ToulBar2::vnsOutput << "InstanceVnsBestTime " << MPI_Wtime() - startTime << endl;
        ToulBar2::vnsOutput << "Cost " << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->Cost2ADCost(bestUb) << std::setprecision(DECIMAL_POINT) << endl;
        ToulBar2::vnsOutput << "Solution ";
        for (map<int, Value>::iterator it = bestSolution.begin();
             it != bestSolution.end(); ++it)
            ToulBar2::vnsOutput << (*it).first << "=" << (*it).second << " ";
        ToulBar2::vnsOutput << endl;
        ToulBar2::vnsOutput << "# ------------------------------------------"
                            << endl;
    } else {
        //        cout << "# ------------------------------------------" << endl;
        //        cout << "InstanceVnsBestTime " << BestTime << endl;
        //        cout << "Cost "  << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->Cost2ADCost(bestUb) << std::setprecision(DECIMAL_POINT) << endl;
        //        cout << "Solution ";
        //        for (map<int, Value>::iterator it = bestSolution.begin();
        //                it != bestSolution.end(); ++it)
        //            cout << (*it).first << "=" << (*it).second << " ";
        //        cout << endl;
        //        cout << "# ------------------------------------------" << endl;
    }
    if (improved && ToulBar2::verbose >= 0) {
        if (!ToulBar2::bayesian)
            cout << "New solution: " << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->Cost2ADCost(bestUb) << std::setprecision(DECIMAL_POINT) << " in " << (MPI_Wtime() - startTime) << " seconds." << endl;
        else
            cout << "New solution: " << bestUb << " energy: " << -(wcsp->Cost2LogProb(bestUb) + ToulBar2::markov_log) << " prob: " << wcsp->Cost2Prob(bestUb) * Exp(ToulBar2::markov_log) << " in " << MPI_Wtime() - startTime << " seconds." << endl;
    }
    if (improved && ToulBar2::showSolutions) {
        wcsp->printSolution();
        cout << endl;
    }
    if (ToulBar2::writeSolution && ToulBar2::solutionFile != NULL) {
        fseek(ToulBar2::solutionFile, ToulBar2::solutionFileRewindPos, SEEK_SET);
        wcsp->printSolution(ToulBar2::solutionFile);
        fprintf(ToulBar2::solutionFile, "\n");
    }
    if (ToulBar2::xmlflag) {
        cout << "o " << bestUb << endl;
    }
    if (ToulBar2::maxsateval) {
        cout << "o " << bestUb << endl;
    }
    if (ToulBar2::uaieval && !ToulBar2::isZ) {
        ((WCSP*)wcsp)->solution_UAI(bestUb);
    }

    if (ToulBar2::newsolution)
        (*ToulBar2::newsolution)(wcsp->getIndex(), wcsp->getSolver());
}
#endif

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

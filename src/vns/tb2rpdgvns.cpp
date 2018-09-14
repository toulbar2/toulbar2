/*
 * tb2rpdgvns.cpp
 *
 *      Author: Abdelkader Ouali
 *      PhD Student: LITIO, University of Oran ; GREYC, University of Caen.
 */

#include "tb2rpdgvns.hpp"
#include "tb2wcsp.hpp"
#ifdef OPENMPI

// Same as Luby in tb2solver.cpp 
// Luby is used to increase the neighborhood's size (k) 
Long lubyIncrK(Long r) {
  int j = cost2log2(r+1);
  if (r+1 == (1L << j)) return (1L << (j-1));
  else return lubyIncrK(r - (1L << j) + 1);
}

// Conversion Tools
// Solution to Message
void ReplicatedParallelDGVNS::SolToMsg(
        MPIEnv &env0, int cluster, int k, int discrepancy, Cost bestUb,
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
        MPIEnv &env0, Cost bestUb, map<int, Value>& bestSolution)
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
        MPIEnv &env0, int nov, int &cluster, int &k, int &discrepancy,
        Cost &bestUb, map<int, Value>& bestSolution)
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
        MPIEnv &env0, int nov, Cost &bestUb, map<int, Value>& bestSolution)
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

bool ReplicatedParallelDGVNS::solveLS()
{ // main algorithm
    bool complete = false;

    mysrand(abs(ToulBar2::seed) + env0.myrank);

    MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
    startTime = MPI_Wtime();
    ToulBar2::timeOut = timeOut;

    env0.buffsize = (int) wcsp->numberOfVariables() + 108; // 3 : cluster + k + cost, second time, msecond time, localtime,the rest is the size of solution
    env0.sendbuff = new int[env0.buffsize];
    env0.recvbuff = new int[env0.buffsize];
    if (env0.myrank == 0) {
        if (!ToulBar2::vnsParallelSync) {
            complete = radgvns();
        } else {
            complete = rsdgvns();
        }
    } else {
        complete = slave();
    }

    /* Shut down MPI */
    if (ToulBar2::vnsOutput) ToulBar2::vnsOutput << "Search end" << " " << cpuTime() << endl;
    if (ToulBar2::vnsOutput) ToulBar2::vnsOutput.close(); // close the output file

    MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
    double elapsedTime = MPI_Wtime() - startTime;
    double elapsedCPUTime = cpuTime() - ToulBar2::startCpuTime;
    double totalCPUTime = 0.;
    MPI_Reduce( &elapsedCPUTime, &totalCPUTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    if (ToulBar2::verbose >= 0 && env0.myrank == 0) { /* use time on master node */
        cout << "Total CPU time = " << totalCPUTime << " seconds" << endl;
        cout << "Solving real-time = " << elapsedTime << " seconds" << endl;
    }
    MPI_Finalize();

    return complete && (bestUb < MAX_COST);
}

//----------------------- Model Master/Slave -----------------------//
bool ReplicatedParallelDGVNS::radgvns()
{
    if (ToulBar2::verbose >= 1) cout << " RADGVNS kinit=" << ToulBar2::vnsKmin << " kmax=" << ToulBar2::vnsKmax
         << " discrepancy=" << ToulBar2::lds << " neighbor change if improved=" << ToulBar2::vnsNeighborChange
         << " neighbor size synchronization=" << ToulBar2::vnsNeighborSizeSync << " limit number processes=" << ToulBar2::vnsParallelLimit << endl;

    ParallelRandomClusterChoice* h = new ParallelRandomClusterChoice();
    h->init(wcsp, this);
    // MPI data
    MPI_Status status;

    // Initialization
    if (ToulBar2::DEE)
        ToulBar2::DEE_ = ToulBar2::DEE; // enforces PSNS after closing the model
//    wcsp->propagate();        // no initial propagation in the master
//    wcsp->setUb(MAX_COST);
    initVarHeuristic();

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

    if (ToulBar2::vnsOutput) ToulBar2::vnsOutput << "#param " << ToulBar2::vnsNeighborVarHeur << " " << ToulBar2::vnsKmin << " " << ToulBar2::vnsKmax << " " << ToulBar2::lds << " " << ToulBar2::vnsInitSol << " " << env0.processes << " " << env0.ntasks << " " << h->getSize() << endl;

    // Generation of initial Solution
    file = h->getClustersIndex();
    clusterKmax = vector<bool>(file.size(), false);
//    cout << "reading Time=" <<  cpuTime() << endl;
//    cout << "initial solution mode : " << ToulBar2::vnsInitSol << endl;
    bool complete = false;
    bestUb = generateInitSolution(ToulBar2::vnsInitSol, bestSolution, complete);
    DumpBestSol();
    // Get all clusters from the tree decomposition of constraint graph

    int c = 0;
    for (int p = 1; p < npr + 1; ++p) {
        pr pr_p = pr();
        pr_p.cl = c;
        pr_p.k = ToulBar2::vnsKmin;
        pr_p.lds = ToulBar2::lds;
        pr_p.synch = false;
		pr_p.rankNeighb = 0;
        vecPR.push_back(pr_p);
        SolToMsg(env0, pr_p.cl, pr_p.k, ToulBar2::lds, bestUb, bestSolution);
        MPI_Send(&env0.sendbuff[0], env0.buffsize, MPI_INT, p, WORKTAG, MPI_COMM_WORLD);
        c = (c + 1) % file.size();
    }
    Cost pbestUb = MAX_COST;
    map<int, Value> pbestSolution;
    while (npr && !complete && bestUb > ToulBar2::vnsOptimum) {
        MPI_Recv(&env0.recvbuff[0], env0.buffsize, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        int pindex = status.MPI_SOURCE - 1;
        complete = (status.MPI_TAG == WORKTAG);
        if (ToulBar2::restart==1 && vecPR[pindex].k == (int) wcsp->numberOfVariables()) npr = 0;
        MsgToSol2(env0, wcsp->numberOfVariables(), pbestUb, pbestSolution);
        NeighborhoodChange(ToulBar2::vnsNeighborChange, pindex, c, ToulBar2::vnsKmin, ToulBar2::vnsKmax, ToulBar2::lds, ToulBar2::vnsNeighborSizeSync, pbestUb, pbestSolution);
        if (ToulBar2::restart==1 && find(clusterKmax.begin(), clusterKmax.end(), false) == clusterKmax.end()) npr = 0;
        if (!complete && bestUb > ToulBar2::vnsOptimum) {
            SolToMsg(env0, vecPR[pindex].cl, vecPR[pindex].k, vecPR[pindex].lds, bestUb, bestSolution);
            MPI_Send(&env0.sendbuff[0], env0.buffsize, MPI_INT, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);
        }
    }
    for (int p = 1; p < env0.ntasks; ++p) {
        MPI_Send(0, 0, MPI_INT, p, DIETAG, MPI_COMM_WORLD);
    }

    if (ToulBar2::verbose >= 0) {
        if (bestUb < MAX_COST) {
            if (ToulBar2::bayesian) cout << (complete?"Optimum: ":"Best upper-bound: ") << bestUb << " log10like: " << (wcsp->Cost2LogProb(bestUb) + ToulBar2::markov_log)/Log(10.) << " prob: " << wcsp->Cost2Prob(bestUb) * Exp(ToulBar2::markov_log) << " in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << ((ToulBar2::DEE)?(" ( "+to_string(wcsp->getNbDEE())+" removals by DEE)"):"") << " and " << cpuTime() - ToulBar2::startCpuTime << " seconds (master process)." << endl;
            else cout << (complete?"Optimum: ":"Best upper-bound: ") << bestUb << " in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << ((ToulBar2::DEE)?(" ( "+to_string(wcsp->getNbDEE())+" removals by DEE)"):"") << " and " << cpuTime() - ToulBar2::startCpuTime << " seconds (master process)." << endl;
        } else cout << "No solution found by " << ((ToulBar2::vnsNeighborVarHeur==RANDOMVAR)?"":"DG") << "VNS in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << ((ToulBar2::DEE)?(" ( "+to_string(wcsp->getNbDEE())+" removals by DEE)"):"") << " and " << cpuTime() - ToulBar2::startCpuTime << " seconds (master process)." << endl;
    }

    return complete && (bestUb < MAX_COST);
}

bool ReplicatedParallelDGVNS::rsdgvns()
{
    if (ToulBar2::verbose >= 1) cout << " RSDGVNS kinit=" << ToulBar2::vnsKmin << " kmax=" << ToulBar2::vnsKmax
         << " discrepancy=" << ToulBar2::lds << " neighbor change if improved=" << ToulBar2::vnsNeighborChange
         << " neighbor size synchronization=" << ToulBar2::vnsNeighborSizeSync << " limit number processes=" << ToulBar2::vnsParallelLimit << endl;

    ParallelRandomClusterChoice* h = new ParallelRandomClusterChoice();
    h->init(wcsp, this);

    // MPI data
    MPI_Status status;

    // Initialization
    map<int, Value> bestInterSolution;
    Cost bestInterUb = MAX_COST;
    if (ToulBar2::DEE)
        ToulBar2::DEE_ = ToulBar2::DEE; // enforces PSNS after closing the model
    //    wcsp->propagate();        // no initial propagation in the master
    //    wcsp->setUb(MAX_COST);
    initVarHeuristic();

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
    int discrepancyInit = ToulBar2::lds;

    if (ToulBar2::vnsOutput) ToulBar2::vnsOutput << "#param " << ToulBar2::vnsNeighborVarHeur << " " << ToulBar2::vnsKmin << " " << ToulBar2::vnsKmax << " " << ToulBar2::lds << " " << ToulBar2::vnsInitSol << " " << env0.processes << " " << env0.ntasks << " " << h->getSize() << endl;

    // Generation of initial Solution
    file = h->getClustersIndex();
    clusterKmax = vector<bool>(file.size(), false);
    if (ToulBar2::vnsOutput) ToulBar2::vnsOutput << "cslsize " << file.size() << " npr= " << npr << endl;

    bool complete = false;
    bestUb = generateInitSolution(ToulBar2::vnsInitSol, bestSolution, complete);
    if (bestUb == wcsp->getLb()) complete = true;
    for (map<int, Value>::iterator it = bestSolution.begin();
            it != bestSolution.end(); ++it)
        bestInterSolution[(*it).first] = (*it).second;
    bestInterUb = bestUb;
    DumpBestSol();
    // Get all clusters from the tree decomposition of constraint graph

    int k = ToulBar2::vnsKmin;
    int c = 0;
    bool getin = true;
    while (npr && k < ToulBar2::vnsKmax && getin) {

        for (int p = 1; p < npr + 1; ++p) {
            SolToMsg(env0, c, k, ToulBar2::lds, bestUb, bestSolution);
            MPI_Send(&env0.sendbuff[0], env0.buffsize, MPI_INT, p, WORKTAG,
                    MPI_COMM_WORLD);
            c = (c + 1) % file.size();
        }

        int finished = 0;
        //S"
        Cost pBestUb;
        map<int, Value> pBestSolution;

        while (finished < npr) {
            MPI_Recv(&env0.recvbuff[0], env0.buffsize, MPI_INT, MPI_ANY_SOURCE,
                    MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            MsgToSol2(env0, wcsp->numberOfVariables(), pBestUb, pBestSolution);
            if (pBestUb < bestInterUb) {
                bestInterUb = pBestUb;
                for (int v = 0; v < (int) wcsp->numberOfVariables(); v++) {
                    bestInterSolution[v] = pBestSolution[v];
                }
            }
            finished++;
        }
        if (bestInterUb < bestUb) {
            k = ToulBar2::vnsKmin;
            ToulBar2::lds = discrepancyInit;
            bestUb = bestInterUb;
            for (int v = 0; v < (int) wcsp->numberOfVariables(); v++) {
                bestSolution[v] = bestInterSolution[v];
            }
            DumpBestSol();
        } else {
            if (k < ToulBar2::vnsKmax && k < (int) wcsp->numberOfVariables())
                k++;
            else if (ToulBar2::restart >= 2) {
                k = ToulBar2::vnsKmin;
                ToulBar2::lds++;
            }
        }
        if (bestUb <= ToulBar2::vnsOptimum) {
            getin = false;
        }
    }
    for (int p = 1; p < env0.ntasks; ++p) {
        MPI_Send(0, 0, MPI_INT, p, DIETAG, MPI_COMM_WORLD);
    }

    if (ToulBar2::verbose >= 0) {
        if (bestUb < MAX_COST) {
            if (ToulBar2::bayesian) cout << (complete?"Optimum: ":"Best upper-bound: ") << bestUb << " log10like: " << (wcsp->Cost2LogProb(bestUb) + ToulBar2::markov_log)/Log(10.) << " prob: " << wcsp->Cost2Prob(bestUb) * Exp(ToulBar2::markov_log) << " in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << ((ToulBar2::DEE)?(" ( "+to_string(wcsp->getNbDEE())+" removals by DEE)"):"") << " and " << cpuTime() - ToulBar2::startCpuTime << " seconds (master process)." << endl;
            else cout << (complete?"Optimum: ":"Best upper-bound: ") << bestUb << " in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << ((ToulBar2::DEE)?(" ( "+to_string(wcsp->getNbDEE())+" removals by DEE)"):"") << " and " << cpuTime() - ToulBar2::startCpuTime << " seconds (master process)." << endl;
        } else cout << "No solution found by " << ((ToulBar2::vnsNeighborVarHeur==RANDOMVAR)?"":"DG") << "VNS in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << ((ToulBar2::DEE)?(" ( "+to_string(wcsp->getNbDEE())+" removals by DEE)"):"") << " and " << cpuTime() - ToulBar2::startCpuTime << " seconds (master process)." << endl;
    }

    return complete && (bestUb < MAX_COST);
}

bool ReplicatedParallelDGVNS::slave()
{
    MPI_Status status;
    ParallelRandomClusterChoice* h = new ParallelRandomClusterChoice();
    h->init(wcsp, this);

    // Initialization (prep data structure)
    if (ToulBar2::DEE)
        ToulBar2::DEE_ = ToulBar2::DEE; // enforces PSNS after closing the model
    wcsp->propagate();
    wcsp->setUb(MAX_COST);
    initVarHeuristic();

    // loop for getting new initial solution to initialize the search process
    bool complete = false;
    while (true) {
        MPI_Recv(&env0.recvbuff[0], env0.buffsize, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == DIETAG) return complete;  // warning! always wait for DIETAG before stopping
        /* launch Vns/Lds+cp */
        complete = complete || VnsLdsCP(env0, h);
        /* Send the result back */
        MPI_Send(&env0.sendbuff[0], env0.buffsize, MPI_INT, 0, complete, MPI_COMM_WORLD);
    }
    return complete;
}

//------------------- Base functions -----------------------//

bool ReplicatedParallelDGVNS::VnsLdsCP(MPIEnv &env0, ParallelRandomClusterChoice* h)
{
    int cluster, k, discrepancy;
    MsgToSol(env0, wcsp->numberOfVariables(), cluster, k, discrepancy, bestUb,
            bestSolution);
    ToulBar2::vnsKcur = k;
    ToulBar2::lds = discrepancy;
    for (map<int, Value>::iterator it = bestSolution.begin();
            it != bestSolution.end(); ++it)
        lastSolution[(*it).first] = (*it).second;
    lastUb = bestUb;
    //vns/lds+cp
    set<int> neighborhood = h->SlaveGetNeighborhood(cluster, k); // based shuffle
    vector<int> variables;
    variables.reserve(wcsp->numberOfVariables());
    vector<int> values;
    values.reserve(wcsp->numberOfVariables());
    for (int v = 0; v < (int) wcsp->numberOfVariables(); v++) {
        if (neighborhood.find(v) == neighborhood.end()) {
            variables.push_back(v);
            values.push_back(lastSolution[v]);
        }
    }

    //repair
    bool complete = false;
    if (ToulBar2::lds) complete = repair_recursiveSolve(discrepancy, variables, values, bestUb);
    else complete = repair_recursiveSolve(variables, values, bestUb);
    SolToMsg2(env0, lastUb, lastSolution);
    return complete;
}

void ReplicatedParallelDGVNS::NeighborhoodChange(
        int strategy, int p, int &c, int kinit, int kmax, int discrepancy,
        bool synch, Cost pBestUb, map<int, Value>& pBestSolution)
{
    switch (strategy) {
    case 0:
        ChangeClusterAlways(p, c, kinit, kmax, discrepancy, synch, pBestUb,
                pBestSolution);
        break;
    case 1:
        ChangeClusterWhenNotImproved(p, c, kinit, kmax, discrepancy, synch,
                pBestUb, pBestSolution);
        break;
    default:
        ChangeClusterAlways(p, c, kinit, kmax, discrepancy, synch, pBestUb,
                pBestSolution);
        break;
    }
}

void ReplicatedParallelDGVNS::ChangeClusterAlways(
        int p, int &c, int kinit, int kmax, int discrepancy, bool synch,
        Cost pBestUb, map<int, Value>& pBestSolution)
{
    bool improved = (pBestUb < bestUb);
    if (pBestUb <= bestUb) {
        bestUb = pBestUb;
        for (int v = 0; v < (int) wcsp->numberOfVariables(); v++) {
            bestSolution[v] = pBestSolution[v];
        }
//		ToulBar2::vnsOutput << "# ------------Change A------------------------------" << endl;
//		ToulBar2::vnsOutput << "cluster=" << c <<" process="<< p << " k="<< vecPR[p].k << " Time " << BestTime << endl;
        if (improved)
            DumpBestSol();
        for (int i = 0; synch && i < (int) vecPR.size(); i++) {
            vecPR[i].synch = true;
        }
        vecPR[p].synch = false;
        vecPR[p].cl = c;
        vecPR[p].k = kinit;
        vecPR[p].lds = discrepancy;
    } else {
        vecPR[p].cl = c;
        if (synch && vecPR[p].synch == true) {
            vecPR[p].k = kinit;
            vecPR[p].lds = discrepancy;
            vecPR[p].synch = false;
        } else {
            if (vecPR[p].k < kmax && vecPR[p].k < (int) wcsp->numberOfVariables()) {
				if (ToulBar2::incr_k_algo == 1) { 
					if (lubyIncrK(vecPR[p].rankNeighb) == 8) {  // 8 because Luby returns 2^k values
					++(vecPR[p].k);
				}
			} else {
				++(vecPR[p].k);
			}
	    }
            else if (ToulBar2::restart > 1) {
                vecPR[p].k = kinit;
                if (vecPR[p].lds >= 4) vecPR[p].lds *= 2;
                else vecPR[p].lds++;
                vecPR[p].lds = max(1LL, min((Long) vecPR[p].lds, (Long) wcsp->numberOfVariables() * (wcsp->getMaxDomainSize()-1)));
            } else {
                clusterKmax[c] = true;
            }
        }
    }
    c = (c + 1) % file.size();
    ++(vecPR[p].rankNeighb);
}

void ReplicatedParallelDGVNS::ChangeClusterWhenNotImproved(
        int p, int &c, int kinit, int kmax, int discrepancy, bool synch,
        Cost pBestUb, map<int, Value>& pBestSolution)
{
    bool improved = (pBestUb < bestUb);
    if (pBestUb <= bestUb) {
        bestUb = pBestUb;
        for (int v = 0; v < (int) wcsp->numberOfVariables(); v++) {
            bestSolution[v] = pBestSolution[v];
        }
//		ToulBar2::vnsOutput << "# -------NOIMP-----------------------------------" << endl;
//		ToulBar2::vnsOutput << "cluster="<< c << " process="<< p << " k="<< vecPR[p].k << " Time=" << BestTime <<endl;
        if (improved)
            DumpBestSol();
        for (int i = 0; synch && i < (int) vecPR.size(); i++) {
            vecPR[i].synch = true;
        }
        vecPR[p].synch = false;
        // keep the current cluster
        vecPR[p].k = kinit;
        vecPR[p].lds = discrepancy;
    } else {
        vecPR[p].cl = c; // change cluster when no improvement is considered
        if (synch && vecPR[p].synch == true) {
            vecPR[p].k = kinit;
            vecPR[p].lds = discrepancy;
            vecPR[p].synch = false;
        } else {
            if (vecPR[p].k < kmax && vecPR[p].k < (int) wcsp->numberOfVariables()) {
                if (ToulBar2::incr_k_algo == 1) { 
					if (lubyIncrK(vecPR[p].rankNeighb) == 8) {  // 8 because Luby returns 2^k values
						++(vecPR[p].k);
					}
				} else {
					++(vecPR[p].k);
				}
			}
            else if (ToulBar2::restart >= 2) {
                vecPR[p].k = kinit;
                vecPR[p].lds++;
            }
        }
        c = (c + 1) % file.size(); // change cluster when no improvement is considered
    }
    ++(vecPR[p].rankNeighb);
}

void ReplicatedParallelDGVNS::DumpBestSol()
{
    // Save
    if (ToulBar2::vnsOutput) {
        //cout << "# ------------------------------------------" << endl;
        ToulBar2::vnsOutput << "InstanceVnsBestTime " << MPI_Wtime() - startTime << endl;
        ToulBar2::vnsOutput << "Cost " << bestUb << endl;
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
//        cout << "Cost " << bestUb << endl;
//        cout << "Solution ";
//        for (map<int, Value>::iterator it = bestSolution.begin();
//                it != bestSolution.end(); ++it)
//            cout << (*it).first << "=" << (*it).second << " ";
//        cout << endl;
//        cout << "# ------------------------------------------" << endl;
    }
    if (ToulBar2::verbose >= 0) cout << "Best solution: " << bestUb << " in " << MPI_Wtime() - startTime << " seconds." << endl;
    if (ToulBar2::uaieval) {
        wcsp->setSolution(&bestSolution);
        ((WCSP*) wcsp)->solution_UAI(bestUb);
    }

}
#endif

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

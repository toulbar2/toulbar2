/*
 * \file tb2cpdgvns.hpp
 * \brief cooperative parallel decomposition guided variable neighborhood search method
 *
 *      Author: Abdelkader Ouali
 *      PhD Student: LITIO, University of Oran ; GREYC, University of Caen.
 */

#ifndef TB2CPDGVNS_HPP_
#define TB2CPDGVNS_HPP_
#ifdef OPENMPI

//#include <omp.h> // include OpenMP Header for parallel execution
#include <mpi.h> // include MPI Header for parallel execution

#include "tb2vns.hpp"

// MPI definitions
const int WORKTAG = 1; // also used to return the search was complete
const int DIETAG = 2;

inline bool MPI_interrupted()
{
    MPI_Status status;
    int flag = 0;
    MPI_Iprobe(MPI_ANY_SOURCE, DIETAG, MPI_COMM_WORLD, &flag, &status);
    return flag;
}

typedef int* unit_of_work_t;
typedef int* unit_result_t;

class MPIEnv {
public:
    int ntasks;
    int myrank; // rang du processus
    int processes; // number of process to be initialize according to the number of cluster
    int buffsize;
    int* sendbuff;
    int* recvbuff;
};

class NeighborhoodStructure;
class ParallelRandomClusterChoice;

class CooperativeParallelDGVNS : public LocalSearch {
protected:
    MPIEnv env0;
    vector<int> file;
    int BestTimeS;
    int BestTimeMS;

public:
    CooperativeParallelDGVNS(Cost initUpperBound, MPIEnv env0Global)
        : LocalSearch(initUpperBound)
    {
        env0 = env0Global;
        BestTimeS = 0;
        BestTimeMS = 0;
    }
    ~CooperativeParallelDGVNS() {}

    bool solve(bool first = true);
    // Model
    void master();
    void slave();
    // Base functions
    void VnsLdsCP(MPIEnv& env0, double btime, ParallelRandomClusterChoice* h);
    // Clustering functions
    uint getCluster();

    // Conversions tools
    void SolToMsg(MPIEnv& env0, uint cluster, uint numberclu, int kinit, int kmax, Cost bestUb, int sec, int msec, map<int, int>& bestSolution);
    void MsgToSol(MPIEnv& env0, int nov, uint& cluster, uint& numberclu, int& k, int& kmax, Cost& bestUb, int& sec, int& msec, map<int, int>& bestSolution);
};

#endif
#endif /* TB2CPDGVNS_HPP_ */

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

/*
 * \file tb2rpdgvns.hpp
 * \brief replicated parallel decomposition guided variable neighborhood search method
 *
 *      Author: Abdelkader Ouali
 *      PhD Student: LITIO, University of Oran ; GREYC, University of Caen.
 */

#ifndef TB2RPDGVNS_HPP_
#define TB2RPDGVNS_HPP_
#ifdef OPENMPI

#include "tb2cpdgvns.hpp"

struct pr { // use it later
    int cl;
    int k;
    int lds;
    bool synch;
};

typedef vector<pr> PR;

class ReplicatedParallelDGVNS : public LocalSearch {
protected:
    MPIEnv env0;
    vector<int> file;
    PR vecPR;
    //    vector<bool> clusterKmax;  // clusterKmax[c] is true if cluster c has its k = kmax
    double startTime;

public:
    ReplicatedParallelDGVNS(Cost initUpperBound, MPIEnv env0Global)
        : LocalSearch(initUpperBound)
        , env0(env0Global)
        , startTime(.0)
    {
    }
    ~ReplicatedParallelDGVNS() {}

    bool solve(bool first = true);
    // Model
    bool radgvns();
    bool rsdgvns();
    bool slave();
    void NeighborhoodChange(int strategy, int p, int& c, int kinit, int kjump, int kmax, int ldsmin, int ldsmax, bool synch, Cost pBestUb, map<int, Value>& pBestSolution);
    void DumpBestSol(bool improved = true);
    bool VnsLdsCP(MPIEnv& env0, ParallelRandomClusterChoice* h);

    // strategies

    void ChangeClusterAlways(int p, int& c, int kinit, int kjump, int kmax, int ldsmin, int ldsmax, bool synch, Cost pBestUb, map<int, Value>& pBestSolution);
    void ChangeClusterWhenNotImproved(int p, int& c, int kinit, int kjump, int kmax, int ldsmin, int ldsmax, bool synch, Cost pBestUb, map<int, Value>& pBestSolution);

    //Conversions tools
    void SolToMsg(MPIEnv& env0, int cluster, int k, int discrepancy, Cost bestUb, map<int, Value>& bestSolution);
    void SolToMsg2(MPIEnv& env0, Cost bestUb, map<int, Value>& bestSolution);
    void MsgToSol(MPIEnv& env0, int nov, int& cluster, int& k, int& discrepancy, Cost& bestUb, map<int, Value>& bestSolution);
    void MsgToSol2(MPIEnv& env0, int nov, Cost& bestUb, map<int, Value>& bestSolution);
};

#endif
#endif /* TB2RPDGVNS_HPP_ */

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

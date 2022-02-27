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

class SolMsg {
private:
    friend class serialization::access;
    template <class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar& cluster;
        ar& k;
        ar& lds;
        ar& bestUb;
        ar& bestSolution;
    }

public:
    int cluster;
    int k;
    int lds;
    Cost bestUb;
    map<int, Value> bestSolution;

    SolMsg() {} // default constructor added to avoid boost/serialization/access.hpp:130:9: error
    SolMsg(int cluster_, int k_, int lds_, Cost bestUb_, map<int, Value>& bestSolution_)
        : cluster(cluster_)
        , k(k_)
        , lds(lds_)
        , bestUb(bestUb_)
        , bestSolution(bestSolution_)
    {
    }
    void get(int& cluster_, int& k_, int& lds_, Cost& bestUb_, map<int, Value>& bestSolution_)
    {
        cluster_ = cluster;
        k_ = k;
        lds_ = lds;
        bestUb_ = bestUb;
        bestSolution_ = bestSolution;
    }
};

class SolMsg2 {
private:
    friend class serialization::access;
    template <class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar& bestUb;
        ar& bestSolution;
    }

public:
    Cost bestUb;
    map<int, Value> bestSolution;

    SolMsg2() {} // default constructor added to avoid boost/serialization/access.hpp:130:9: error
    SolMsg2(Cost bestUb_, map<int, Value>& bestSolution_)
        : bestUb(bestUb_)
        , bestSolution(bestSolution_)
    {
    }
    void get(Cost& bestUb_, map<int, Value>& bestSolution_)
    {
        bestUb_ = bestUb;
        bestSolution_ = bestSolution;
    }
};

class ReplicatedParallelDGVNS : public LocalSearch {
protected:
    vector<int> file;
    PR vecPR;
    //    vector<bool> clusterKmax;  // clusterKmax[c] is true if cluster c has its k = kmax

public:
    ReplicatedParallelDGVNS(Cost initUpperBound)
        : LocalSearch(initUpperBound)
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
    bool VnsLdsCP(SolMsg& solmsg, ParallelRandomClusterChoice* h, SolMsg2& solmsg2);

    // strategies
    void ChangeClusterAlways(int p, int& c, int kinit, int kjump, int kmax, int ldsmin, int ldsmax, bool synch, Cost pBestUb, map<int, Value>& pBestSolution);
    void ChangeClusterWhenNotImproved(int p, int& c, int kinit, int kjump, int kmax, int ldsmin, int ldsmax, bool synch, Cost pBestUb, map<int, Value>& pBestSolution);
};

#endif
#endif /* TB2RPDGVNS_HPP_ */

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

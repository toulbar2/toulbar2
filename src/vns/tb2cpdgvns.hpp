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

#include "tb2vns.hpp"

typedef int* unit_of_work_t;
typedef int* unit_result_t;

class NeighborhoodStructure;
class ParallelRandomClusterChoice;

class SolutionMessage {
private:
    friend class serialization::access;
    template <class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar& cluster;
        ar& numberclu;
        ar& k;
        ar& kmax;
        ar& sec;
        ar& msec;
        ar& bestUb;
        ar& bestSolution;
    }

public:
    uint cluster;
    uint numberclu;
    int k;
    int kmax;
    int sec;
    int msec;
    Cost bestUb;
    map<int, Value> bestSolution;

    SolutionMessage() {} // default constructor added to avoid boost/serialization/access.hpp:130:9: error
    SolutionMessage(uint cluster_, uint numberclu_, int kinit_, int kmax_, int sec_, int msec_, Cost bestUb_, map<int, Value>& bestSolution_)
        : cluster(cluster_)
        , numberclu(numberclu_)
        , k(kinit_)
        , kmax(kmax_)
        , sec(sec_)
        , msec(msec_)
        , bestUb(bestUb_)
        , bestSolution(bestSolution_)
    {
    }
    void get(uint& cluster_, uint& numberclu_, int& k_, int& kmax_, int& sec_, int& msec_, Cost& bestUb_, map<int, Value>& bestSolution_)
    {
        cluster_ = cluster;
        numberclu_ = numberclu;
        k_ = k;
        kmax_ = kmax;
        sec_ = sec;
        msec_ = msec;
        bestUb_ = bestUb;
        bestSolution_ = bestSolution;
    }
};

class CooperativeParallelDGVNS : public LocalSearch {
protected:
    vector<int> file;
    int BestTimeS;
    int BestTimeMS;

public:
    CooperativeParallelDGVNS(Cost initUpperBound)
        : LocalSearch(initUpperBound)
    {
        BestTimeS = 0;
        BestTimeMS = 0;
    }
    ~CooperativeParallelDGVNS() {}

    bool solve(bool first = true);
    // Model
    void master();
    void slave();
    // Base functions
    void VnsLdsCP(SolutionMessage& solmsg, double btime, ParallelRandomClusterChoice* h);
    // Clustering functions
    uint getCluster();
};

#endif
#endif /* TB2CPDGVNS_HPP_ */

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

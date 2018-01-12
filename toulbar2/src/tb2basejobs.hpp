#ifndef TB2BASEJOBS_HPP
#define TB2BASEJOBS_HPP
#ifdef USEMPI
#include <vector>
#include <string>
#include <mpi.h>

using namespace std;

class BaseJobs {
public:
    BaseJobs(string jobsfile);
    ~BaseJobs();
    string get_jobsfile() { return jobsfile; };
    virtual void init_jobs(string jobsfile);
    unsigned nb_jobs() { return jobs.size(); };
    virtual bool next_job(string& wcsp_id);
    virtual bool find_next_job(int myjob);
    void init_node();
    void shutdown_node();
    int mpi_rank();
    void spin_off_solver(int source);
    virtual void distribute_jobs();
    virtual bool request_job();
    string get_current_job();

protected:
    vector<string> jobs;
    string jobsfile;
    int mpi_rank_;
    int current_job;
    bool started;
    MPI_Status status_;
    int tag_;
    int running;
};

#endif
#endif

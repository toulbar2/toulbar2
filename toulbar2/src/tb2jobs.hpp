#ifndef TB2JOBS_HPP
#define TB2JOBS_HPP
#ifdef USEMPI
#include <vector>
#include <string>
#include <mpi.h>
#include "tb2types.hpp"
#include "tb2sequencehandler.hpp"

using namespace std;

class Jobs
{
public:
  Jobs(string jobsfile);
  ~Jobs();
  void init_jobs(string jobsfile);
  unsigned nb_jobs() { return jobs.size();};
  bool next_job(string & wcsp_id);
  void send_results(Cost new_cost);
  bool request_sequence();
  bool find_next_job(int myjob);
  void init_node();
  void shutdown_node();
  int mpi_rank();
  void send_seqid_and_cost(unsigned seqid, Cost new_cost);
  tuple<unsigned, Cost> receive_seqid_and_cost(int source);
  void process_job_results(int slave_job, int slave_rank);
  void spin_off_solver(int source);
  void send_new_sequence(int source);
  void distribute_jobs();
  bool request_job();
  string get_current_job();
  unsigned get_current_sequence(){ return current_sequence;}
protected:
  vector<string> jobs;
  vector< vector <string> > sequences;
  int mpi_rank_;
  int current_job;
  unsigned current_sequence;
  bool started;
  MPI_Status status_;
  int tag_;
  int running;
};

class GetAnotherJob
{

public:
    GetAnotherJob() {
        if (ToulBar2::verbose >= 2)
            cout << "Mpi mode: looking for a new sequence for current job" << endl;
    }
};

#endif
#endif

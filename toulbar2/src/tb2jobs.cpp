#ifdef USEMPI
#include "tb2jobs.hpp"
#include "tb2types.hpp"
#include <fstream>
#include <sys/stat.h>


// Note: current_job starts from 1
// current_sequence starts from 0

Jobs::Jobs(string jobsfile):mpi_rank_(0),
                            current_job(0),
                            started(false)
{
  init_node();
  init_jobs(jobsfile);
  MPI_Comm_size(MPI_COMM_WORLD, &running);
}

Jobs::~Jobs()
{
}


// Loading wcsp files
void Jobs::init_jobs(string jobsfile)
{
  ifstream is(jobsfile);
  string filename;
  while(is)
    {
      is >> filename;
      if (is)
        {
          ifstream testif(filename);
          if (!testif.good())
            {
              cout << "Couldn't find file " << filename << endl;
              exit(1);
            }
          else
            jobs.push_back(filename);
        }
    }
  cout << jobs.size() << " jobs successfully initialized" << endl;
  
}

bool Jobs::next_job(string & wcsp_id)
{
  if (mpi_rank()==0) // MASTER NODE
    {
      distribute_jobs();
      // Only Master decides when to stop.
      return false;
    }
  else  // SLAVE NODE
    {
      bool hired = request_job();
      if (hired)
        wcsp_id = jobs[current_job-1];
      return hired;
    }
}


// This method sends results ... and request sequence
void Jobs::send_results(Cost new_cost)
{
  MPI_Send( &current_job, 1, MPI_INT, 0, 1, MPI_COMM_WORLD );
  send_seqid_and_cost(current_sequence, new_cost);
}

bool Jobs::request_sequence()
{
  MPI_Recv( &current_sequence, 1, MPI_UNSIGNED, 0, MPI_ANY_TAG, MPI_COMM_WORLD, & status_ );
  if (current_sequence>=0 && status_.MPI_TAG==1)
    return true;
  else
    {
      current_job = 0;
      return false;
    }
}

bool Jobs::find_next_job(int myjob)
{
  if (!jobs.size())
    {
      // If no more jobs, check sequences from ongoing jobs
      if (!ToulBar2::sequence_handler->distribute_bkb_and_sequence(current_job, current_sequence))
        {
          current_job=0;
          return false;
        }
      else
        return true;
    }
  else
    {
      current_job = jobs.size();
      ToulBar2::sequence_handler->distribute_sequence(current_job-1, current_sequence); // current_job starts from 1...
      jobs.pop_back();
      return true;
    }

}

void Jobs::init_node()
{
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank_);
}

void Jobs::shutdown_node()
{
  started = false;
  current_job=0;
}

int Jobs::mpi_rank()
{
  return mpi_rank_;
}


void Jobs::send_seqid_and_cost(unsigned seqid, Cost new_cost)
{
  MPI_Send( & seqid, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD );
  MPI_Send( & new_cost, 1, MPI_LONG_LONG, 0, 0, MPI_COMM_WORLD );
}

tuple<unsigned, Cost> Jobs::receive_seqid_and_cost(int source)
{
  Cost new_cost;
  unsigned seqid;
  MPI_Recv( &seqid , 1, MPI_UNSIGNED, source, MPI_ANY_TAG, MPI_COMM_WORLD, & status_ );
  MPI_Recv( &new_cost , 1, MPI_LONG_LONG, source, MPI_ANY_TAG, MPI_COMM_WORLD, & status_ );
  return make_tuple(seqid,new_cost);
}

void Jobs::process_job_results(int slave_job, int slave_rank)
{
  Cost new_cost;
  unsigned seqid;
  tie(seqid, new_cost) = receive_seqid_and_cost(slave_rank);
  ToulBar2::sequence_handler->update_sequences(slave_job-1, seqid, new_cost);  // jobs start at one
} 

void Jobs::spin_off_solver(int source)
{
  MPI_Send( &current_job, 1, MPI_INT, source, 0, MPI_COMM_WORLD );
}

void Jobs::send_new_sequence(int source)
{
  MPI_Send( &current_sequence, 1, MPI_UNSIGNED, source, 1, MPI_COMM_WORLD );
}
// TAGS: from slave: 0 -> request from main 1 -> request_from_solver
// from master: 0 -> spin off 1 -> everything is allright
void Jobs::distribute_jobs()
{
  // This is where Master communicates with slaves.
  // We use MPI Tags to determine the type of communication
  while(true)
    {
      bool request_from_solver = false;
      bool found_next_sequence = false;
      int slave_current_job = 0;
      MPI_Recv( &slave_current_job, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, & status_ );
      if (status_.MPI_TAG==1)
        // The request was sent from the solver: if possible,
        // we will just send a sequence and avoid reloading a wcsp and rerun pre-processing
        request_from_solver=true;
      if (slave_current_job>0)
        process_job_results(slave_current_job,status_.MPI_SOURCE);
      if (request_from_solver)
        {
          found_next_sequence = ToulBar2::sequence_handler->distribute_sequence(slave_current_job-1, current_sequence); // current_job starts from 1...
          if (!found_next_sequence) {
            // We have to look for a new job, no sequences left. The solver has nothing left to do, 
            // the main will handle the request
            spin_off_solver(status_.MPI_SOURCE);
          }
          else
            {
              send_new_sequence(status_.MPI_SOURCE);
            }
        }
      if (!request_from_solver)
        {
          bool some_jobs_left = find_next_job(slave_current_job);
          if (!some_jobs_left)
            {
              // Sending tag 0, no more jobs
              MPI_Send( &current_job, 1, MPI_INT, status_.MPI_SOURCE, 0, MPI_COMM_WORLD );
              running--;
            }
          else
            {
              MPI_Send( &current_job, 1, MPI_INT, status_.MPI_SOURCE, 1, MPI_COMM_WORLD );
              MPI_Send( &current_sequence, 1, MPI_UNSIGNED, status_.MPI_SOURCE, 1, MPI_COMM_WORLD );
            }
        }
      if (running<=1)
        break;
    }
}

bool Jobs::request_job()
{
  int new_job;
  // Send request to master
  MPI_Send(&current_job, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
  // Receive signal from master
  MPI_Recv( & new_job, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, & status_ );
  // If tag is 0, no more jobs
  if (status_.MPI_TAG==0)
    {
      return false;
    }
  else
    {
      MPI_Recv( &current_sequence, 1, MPI_UNSIGNED, 0, MPI_ANY_TAG, MPI_COMM_WORLD, & status_ );
      if (new_job != current_job)
        {
          //switch_job(new_job);
          current_job = new_job;
        }
    }
  return true;
}


string Jobs::get_current_job()
{
  return jobs[current_job-1];
}

#endif

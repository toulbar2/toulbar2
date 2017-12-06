#ifdef USEMPI
#include "tb2basejobs.hpp"
#include <fstream>
#include <sys/stat.h>


// Note: current_job starts from 1
// current_sequence starts from 0

BaseJobs::BaseJobs(string jobsfile):mpi_rank_(0),
                            current_job(0),
                            started(false)
{
  init_node();
  init_jobs(jobsfile);
  MPI_Comm_size(MPI_COMM_WORLD, &running);
}

BaseJobs::~BaseJobs()
{
}


// Loading wcsp files
void BaseJobs::init_jobs(string jobsfile)
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
  is.close();
}

bool BaseJobs::next_job(string & wcsp_id)
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


bool BaseJobs::find_next_job(int myjob)
{
  if (!jobs.size())
    {
      return false;
    }
  else
    {
      current_job = jobs.size();
      jobs.pop_back();
      return true;
    }

}

void BaseJobs::init_node()
{
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank_);
}

void BaseJobs::shutdown_node()
{
  started = false;
  current_job=0;
}

int BaseJobs::mpi_rank()
{
  return mpi_rank_;
}


void BaseJobs::spin_off_solver(int source)
{
  MPI_Send( &current_job, 1, MPI_INT, source, 0, MPI_COMM_WORLD );
}

void BaseJobs::distribute_jobs()
{
  // This is where Master communicates with slaves.
  while(true)
    {
      int slave_current_job = 0;
      MPI_Recv( &slave_current_job, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, & status_ );
      // if (slave_current_job>0)
      //   process_job_results(slave_current_job,status_.MPI_SOURCE);
      bool some_jobs_left = find_next_job(slave_current_job);
      if (!some_jobs_left)
        {
          // Sending tag 0, no more jobs
          MPI_Send( &current_job, 1, MPI_UNSIGNED, status_.MPI_SOURCE, 0, MPI_COMM_WORLD );
          running--;
        }
      else
        {
          MPI_Send( &current_job, 1, MPI_UNSIGNED, status_.MPI_SOURCE, 1, MPI_COMM_WORLD );
        }
      if (running<=1)
        break;
    }
}

bool BaseJobs::request_job()
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
  current_job = new_job;
  return true;
}


string BaseJobs::get_current_job()
{
  return jobs[current_job-1];
}

#endif

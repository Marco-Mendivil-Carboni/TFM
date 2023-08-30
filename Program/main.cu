//Libraries

#include <cstdio> //standard input and output library
#include <cstdlib> //standard general utilities library
#include <glob.h> //pathname pattern matching types

#include "inc/cudautil.cuh"
#include "inc/chrsim.cuh"

//Main

int main(int argc, const char **argv)
{
  //check command-line arguments
  if (argc<2){ std::fprintf(stderr,"No arguments.\n"); return EXIT_FAILURE;}
  if (argc>3){ std::fprintf(stderr,"Extra arguments.\n"); return EXIT_FAILURE;}
  if (sizeof(argv[1])>128)
  {
    std::fprintf(stderr,"Directory path too long.\n"); return EXIT_FAILURE;
  }

  //declare auxiliary variables
  char sim_dir[128]; //simulation directory
  std::snprintf(sim_dir,sizeof(sim_dir),"%s",argv[1]);
  FILE *f_ptr_par; //pointer to parameter file
  FILE *f_ptr_out; //pointer to output file
  FILE *f_ptr_log; //pointer to log file
  char f_str[256]; //file name string

  try
  {
    //open log file
    std::snprintf(f_str,sizeof(f_str),"%s/current-progress.log",sim_dir);
    f_ptr_log = mmcc::fopen(f_str,"wt");

    //read parameters and initialize simulation
    std::snprintf(f_str,sizeof(f_str),"%s/adjustable-parameters.dat",sim_dir);
    f_ptr_par = mmcc::fopen(f_str,"rt");
    mmcc::chrsim sim(f_ptr_par);
    std::fclose(f_ptr_par);
    // setup_PRNG<<<n_blocks,threads_block>>>(time(nullptr),state); //put this inside initial configuration?
    // cuda_check( cudaDeviceSynchronize());

    // begin new simulation or continue a previous one
    int sim_idx = 0; //simulation index
    // int tpf_idx = 0; //trajectory positions file index
    // float t = 0.0; //simulation time
    if (argc==2)
    {
      glob_t prev_sims;
      std::snprintf(f_str,sizeof(f_str),"%s/initial-configuration-*",sim_dir);
      if (glob(f_str,0,nullptr,&prev_sims)==0)
      {
        sim_idx = prev_sims.gl_pathc;
      }
      globfree(&prev_sims);

      // print_time(logfile); std::fprintf(logfile,"New simulation started.\n");
      // std::fprintf(logfile,"sim_idx=%03d tpf_idx=%03d\n",sim_idx,tpf_idx); std::fflush(logfile);

      sim.generate_initial_configuration();

      std::snprintf(f_str,sizeof(f_str),"%s/initial-configuration-%03d.gro",sim_dir,sim_idx);
      f_ptr_out = mmcc::fopen(f_str,"wt");
      sim.write_initial_configuration(f_ptr_out);
      std::fclose(f_ptr_out);
    }
    else
    {
    }

    //close log file
    std::fclose(f_ptr_log);
  }
  catch (const mmcc::error& error)
  {
    //do something
  }

  //exit program successfully
  return EXIT_SUCCESS;
}

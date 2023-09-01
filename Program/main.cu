//Includes

#include "inc/chrsim.cuh"
#include "inc/util.hpp"

#include <iostream> //standard input/output stream objects

#include <glob.h> //pathname pattern matching types

//Functions

int main(const int argc, const char **argv)
{
  //check command-line arguments
  if (argc<2){ std::cout<<"no arguments"<<std::endl; return EXIT_FAILURE;}
  if (argc>3){ std::cout<<"extra arguments"<<std::endl; return EXIT_FAILURE;}

  //declare auxiliary variables
  const std::string sim_dir = argv[1]; //simulation directory
  std::ifstream f_par; //parameter file
  std::ofstream f_out; //output file

  //open log file in simulation directory
  mmcc::logger::set_file(sim_dir+"/history.log");

  try
  {
    //read parameters and initialize simulation
    f_par.open(sim_dir+"/adjustable-parameters.dat");
    // mmcc::chrsim sim(f_ptr_par);
    f_par.close();

    // begin new simulation or continue a previous one
    int sim_idx = 0; //simulation index
    int tpf_idx = 0; //trajectory positions file index
    float t = 0.0; //simulation time

    if (argc==2)
    {
      glob_t prev_sims;
      std::string pattern = sim_dir+"/initial-configuration-*";
      if (glob(pattern.c_str(),0,nullptr,&prev_sims)==0)
      {
        sim_idx = prev_sims.gl_pathc;
      }
      globfree(&prev_sims);
      // log: "new simulation started"
      // log: sim_idx tpf_idx ?

      // sim.generate_initial_configuration();
      //write some kind of int_to_string utility function
      // f_out.open(sim_dir+"/initial-condition-");//finish...
      //sim.write_initial_configuration(f_ptr_out);
      // f_out.close();
    }
    else
    {
    }

    // f_out.open(sim_dir+"/trajectory-");//finish...
    // f_out.close();
  }
  catch (const mmcc::error& error)
  {
    //exit program unsuccessfully
    mmcc::logger::record(error.what());
    return EXIT_FAILURE;
  }

  //exit program successfully
  return EXIT_SUCCESS;
}

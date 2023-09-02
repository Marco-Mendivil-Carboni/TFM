//Includes

#include "chrsim.cuh"
#include "util.hpp"

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
  std::string f_path; //file path string
  int sim_idx = 0; //simulation index
  int tpf_idx = 0; //trajectory positions file index
  float t = 0.0; //simulation time ---------------- move to simulation class member variable ---------------------------

  //open log file inside sim_dir
  f_path = sim_dir+"/complete-history.log";
  mmcc::logger::set_file(f_path);

  try
  {
    //read parameters and initialize simulation
    f_path = sim_dir+"/adjustable-parameters.dat";
    f_par.open(f_path); mmcc::check_file(f_par,f_path);
    mmcc::chrsim sim(f_par);
    f_par.close();

    if (argc==2) //begin new simulation
    {
      glob_t prev_sims;
      std::string pattern = sim_dir+"/initial-configuration-*";
      if (glob(pattern.c_str(),0,nullptr,&prev_sims)==0)
      {
        sim_idx = prev_sims.gl_pathc;
      }
      globfree(&prev_sims); // move all this to some function
      mmcc::logger::record("new simulation started");

      sim.generate_initial_configuration();

      f_path = sim_dir+"/initial-configuration-";
      f_path += mmcc::cnfs(sim_idx,3)+".gro";
      f_out.open(f_path); mmcc::check_file(f_out,f_path);
      sim.write_initial_configuration(f_out);
      f_out.close();
    }
    else //continue previous simulation
    {
      sim_idx = std::stoi(argv[2]);
    }
    // log: sim_idx tpf_idx ?

    //perform simulation
    f_path = sim_dir+"/trajectory-positions-";
    f_path += mmcc::cnfs(sim_idx,3)+"-"+mmcc::cnfs(tpf_idx,3)+".trr";
    f_out.open(f_path,std::ios::binary); mmcc::check_file(f_out,f_path);
    //simulation
    f_out.close();
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

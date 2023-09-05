//Includes

#include "chrsim.cuh"
#include "util.hpp"

//Functions

//main function
int main(const int argc, const char **argv)
{
  //check command-line arguments
  if (argc<2){ std::cout<<"no arguments"<<std::endl; return EXIT_FAILURE;}
  if (argc>3){ std::cout<<"extra arguments"<<std::endl; return EXIT_FAILURE;}

  //declare auxiliary variables
  const std::string sim_dir = argv[1]; //simulation directory
  bool new_sim = (argc==2) ? true : false; //make new simulation
  std::ifstream f_inp; //input file
  std::ofstream f_out; //output file
  std::string f_path; //file path string
  std::string pattern; //file path pattern
  int sim_idx; //simulation index
  int tpf_idx; //trajectory positions file index

  //open log file inside simulation directory
  f_path = sim_dir+"/all.log";
  mmcc::logger::set_file(f_path);

  //main try block
  try
  {
    //read parameters and initialize simulation
    f_path = sim_dir+"/adjustable-parameters.dat";
    f_inp.open(f_path);
    mmcc::check_file(f_inp,f_path);
    mmcc::chrsim sim(f_inp); //simulation
    f_inp.close();

    if (new_sim) //begin new simulation
    {
      //set sim_idx and tpf_idx
      pattern = sim_dir+"/initial-condition-*";
      sim_idx = mmcc::glob_count(pattern);
      tpf_idx = 0;

      //generate and write initial condition
      sim.generate_initial_condition();
      f_path = sim_dir+"/initial-condition-";
      f_path += mmcc::cnfs(sim_idx,3)+".gro";
      f_out.open(f_path);
      mmcc::check_file(f_out,f_path);
      sim.write_initial_condition(f_out);
      f_out.close();
    }
    else //continue previous simulation
    {
      //set sim_idx and tpf_idx
      sim_idx = std::stoi(argv[2]);
      pattern = sim_dir+"/trajectory-";
      pattern += mmcc::cnfs(sim_idx,3)+"*";
      tpf_idx = mmcc::glob_count(pattern);

      //load checkpoint
      f_path = sim_dir+"/checkpoint-";
      f_path += mmcc::cnfs(sim_idx,3)+".bin";
      f_inp.open(f_path,std::ios::binary);
      mmcc::check_file(f_inp,f_path);
      sim.load_checkpoint(f_inp);
      f_inp.close();
    }

    //record indexes
    std::string msg = "indexes:";
    msg += " sim_idx = "+mmcc::cnfs(sim_idx,3);
    msg += " tpf_idx = "+mmcc::cnfs(tpf_idx,3);
    mmcc::logger::record(msg);

    //perform simulation
    f_path = sim_dir+"/trajectory-";
    f_path += mmcc::cnfs(sim_idx,3)+"-";
    f_path += mmcc::cnfs(tpf_idx,3)+".trr";
    f_out.open(f_path,std::ios::binary);
    mmcc::check_file(f_out,f_path);
    for (int i_f = 0; i_f<sim.ap.F; ++i_f)
    {
      float prog_pc = (100.0*i_f)/(sim.ap.F); //progress percentage
      mmcc::logger::show_prog_pc(prog_pc);
      // sim.advance_to_next_frame();
      sim.write_trajectory(f_out,i_f);
    }
    f_out.close();

    //save checkpoint
    f_path = sim_dir+"/checkpoint-";
    f_path += mmcc::cnfs(sim_idx,3)+".bin";
    f_out.open(f_path,std::ios::binary);
    mmcc::check_file(f_out,f_path);
    sim.save_checkpoint(f_out);
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

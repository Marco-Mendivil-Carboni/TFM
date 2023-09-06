//Includes

#include "chrsim.cuh" //chromatin simulation
#include "util.hpp" //utilities

//Functions

//main function
int main(
  const int argc, //argument count
  const char **argv) //argument vector
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
  int t_f_idx; //trajectory file index

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
    mmcc::chrsim sim{f_inp}; //simulation
    f_inp.close();

    if (new_sim) //begin new simulation
    {
      //set sim_idx and t_f_idx
      pattern = sim_dir+"/initial-condition-*";
      sim_idx = mmcc::glob_count(pattern);
      t_f_idx = 0;

      //generate and write initial condition
      sim.generate_initial_condition();
      f_path = sim_dir+"/initial-condition-";
      f_path += mmcc::cnfs(sim_idx,3,'0')+".gro";
      f_out.open(f_path);
      mmcc::check_file(f_out,f_path);
      sim.write_initial_condition(f_out);
      f_out.close();
    }
    else //continue previous simulation
    {
      //set sim_idx and t_f_idx
      sim_idx = std::stoi(argv[2]);
      pattern = sim_dir+"/trajectory-";
      pattern += mmcc::cnfs(sim_idx,3,'0')+"*";
      t_f_idx = mmcc::glob_count(pattern);

      //load checkpoint
      f_path = sim_dir+"/checkpoint-";
      f_path += mmcc::cnfs(sim_idx,3,'0')+".bin";
      f_inp.open(f_path,std::ios::binary);
      mmcc::check_file(f_inp,f_path);
      sim.load_checkpoint(f_inp);
      f_inp.close();
    }

    //record indexes
    std::string msg = "indexes:"; //message
    msg += " sim_idx = "+mmcc::cnfs(sim_idx,3,'0');
    msg += " t_f_idx = "+mmcc::cnfs(t_f_idx,3,'0');
    mmcc::logger::record(msg);

    //run simulation
    f_path = sim_dir+"/trajectory-";
    f_path += mmcc::cnfs(sim_idx,3,'0')+"-";
    f_path += mmcc::cnfs(t_f_idx,3,'0')+".trr";
    f_out.open(f_path,std::ios::binary);
    mmcc::check_file(f_out,f_path);
    sim.run_simulation(f_out);
    f_out.close();

    //save checkpoint
    f_path = sim_dir+"/checkpoint-";
    f_path += mmcc::cnfs(sim_idx,3,'0')+".bin";
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

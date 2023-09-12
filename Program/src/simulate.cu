//Includes

#include "chrsim.cuh" //chromatin simulation

#include <iostream> //standard input/output stream objects

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
  f_path = sim_dir+"/all-messages.log";
  mmc::logger::set_file(f_path);

  //main try block
  try
  {
    //read parameters and initialize simulation
    f_path = sim_dir+"/adjustable-parameters.dat";
    f_inp.open(f_path);
    mmc::check_file(f_inp,f_path);
    mmc::parmap par(f_inp); //parameters
    f_inp.close();
    mmc::chrsim sim(par); //simulation

    if (new_sim) //begin new simulation
    {
      //set sim_idx and t_f_idx
      pattern = sim_dir+"/initial-condition-*";
      sim_idx = mmc::glob_count(pattern);
      t_f_idx = 0;

      //generate and write initial condition
      sim.generate_initial_condition();
      f_path = sim_dir+"/initial-condition-";
      f_path += mmc::cnfs(sim_idx,3,'0')+".gro";
      f_out.open(f_path);
      mmc::check_file(f_out,f_path);
      sim.write_initial_condition(f_out);
      f_out.close();
    }
    else //continue previous simulation
    {
      //set sim_idx and t_f_idx
      sim_idx = std::stoi(argv[2]);
      pattern = sim_dir+"/trajectory-";
      pattern += mmc::cnfs(sim_idx,3,'0')+"*";
      t_f_idx = mmc::glob_count(pattern);

      //load checkpoint
      f_path = sim_dir+"/checkpoint-";
      f_path += mmc::cnfs(sim_idx,3,'0')+".bin";
      f_inp.open(f_path,std::ios::binary);
      mmc::check_file(f_inp,f_path);
      sim.load_checkpoint(f_inp);
      f_inp.close();
    }

    //record indexes
    std::string msg = "indexes:"; //message
    msg += " sim_idx = "+mmc::cnfs(sim_idx,3,'0');
    msg += " t_f_idx = "+mmc::cnfs(t_f_idx,3,'0');
    mmc::logger::record(msg);

    //run simulation
    f_path = sim_dir+"/trajectory-";
    f_path += mmc::cnfs(sim_idx,3,'0')+"-";
    f_path += mmc::cnfs(t_f_idx,3,'0')+".trr";
    f_out.open(f_path,std::ios::binary);
    mmc::check_file(f_out,f_path);
    sim.run_simulation(f_out);
    f_out.close();

    //save checkpoint
    f_path = sim_dir+"/checkpoint-";
    f_path += mmc::cnfs(sim_idx,3,'0')+".bin";
    f_out.open(f_path,std::ios::binary);
    mmc::check_file(f_out,f_path);
    sim.save_checkpoint(f_out);
    f_out.close();
  }
  catch (const mmc::error &err) //caught error
  {
    //exit program unsuccessfully
    mmc::logger::record(err.what());
    return EXIT_FAILURE;
  }

  //exit program successfully
  return EXIT_SUCCESS;
}

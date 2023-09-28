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
  bool new_sim = (argc==2)?true:false; //make new simulation
  std::ifstream inp_f; //input file
  std::ofstream out_f; //output file
  std::string pathstr; //file path string
  std::string pathpat; //file path pathpat
  uint sim_idx; //simulation index
  uint t_f_idx; //trajectory file index

  //open log file inside simulation directory
  pathstr = sim_dir+"/all-messages.log";
  mmc::logger::set_file(pathstr);

  //main try block
  try
  {
    //read parameters and initialize simulation
    pathstr = sim_dir+"/adjustable-parameters.dat";
    inp_f.open(pathstr);
    mmc::check_file(inp_f,pathstr);
    mmc::parmap par(inp_f); //parameters
    inp_f.close();
    mmc::chrsim sim(par); //simulation

    if (new_sim) //begin new simulation
    {
      //set sim_idx and t_f_idx
      pathpat = sim_dir+"/initial-condition-*";
      sim_idx = mmc::glob_count(pathpat);
      t_f_idx = 0;

      //generate and write initial condition
      sim.generate_initial_condition();
      pathstr = sim_dir+"/initial-condition-";
      pathstr += mmc::cnfs(sim_idx,3,'0')+".gro";
      out_f.open(pathstr);
      mmc::check_file(out_f,pathstr);
      sim.write_frame_txt(out_f);
      out_f.close();
    }
    else //continue previous simulation
    {
      //set sim_idx and t_f_idx
      sim_idx = std::stoi(argv[2]);
      pathpat = sim_dir+"/trajectory-";
      pathpat += mmc::cnfs(sim_idx,3,'0')+"*";
      t_f_idx = mmc::glob_count(pathpat);

      //load checkpoint
      pathstr = sim_dir+"/checkpoint-";
      pathstr += mmc::cnfs(sim_idx,3,'0')+".bin";
      inp_f.open(pathstr,std::ios::binary);
      mmc::check_file(inp_f,pathstr);
      sim.load_checkpoint(inp_f);
      inp_f.close();
    }

    //record indexes
    std::string msg = "indexes:"; //message
    msg += " sim_idx = "+mmc::cnfs(sim_idx,3,'0');
    msg += " t_f_idx = "+mmc::cnfs(t_f_idx,3,'0');
    mmc::logger::record(msg);

    //run simulation
    pathstr = sim_dir+"/trajectory-";
    pathstr += mmc::cnfs(sim_idx,3,'0')+"-";
    pathstr += mmc::cnfs(t_f_idx,3,'0')+".trr";
    out_f.open(pathstr,std::ios::binary);
    mmc::check_file(out_f,pathstr);
    sim.run_simulation(out_f);
    out_f.close();

    //save checkpoint
    pathstr = sim_dir+"/checkpoint-";
    pathstr += mmc::cnfs(sim_idx,3,'0')+".bin";
    out_f.open(pathstr,std::ios::binary);
    mmc::check_file(out_f,pathstr);
    sim.save_checkpoint(out_f);
    out_f.close();
  }
  catch (const mmc::error &err) //caught error
  {
    //exit program unsuccessfully
    mmc::logger::record(err.what());
    mmc::logger::record("program exited unsuccessfully\n");
    return EXIT_FAILURE;
  }

  //exit program successfully
  mmc::logger::record("program exited successfully\n");
  return EXIT_SUCCESS;
}

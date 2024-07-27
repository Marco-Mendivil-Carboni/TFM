// Includes

#include "chrsim.cuh" // chromatin simulation

// Functions

// main function
int main(
    const int argc, // argument count
    const char **argv) // argument vector
{
  // check command-line arguments
  if (argc < 2)
  {
    std::cout << "no arguments\n";
    return EXIT_FAILURE;
  }
  if (argc > 3)
  {
    std::cout << "extra arguments\n";
    return EXIT_FAILURE;
  }

  // declare auxiliary variables
  const std::string sim_dir = argv[1]; // simulation directory
  bool new_sim = (argc == 2) ? true : false; // make new simulation
  std::ifstream inp_f; // input file
  std::ofstream out_f; // output file
  std::string pathstr; // file path string
  std::string pathpat; // file path pattern
  uint i_s; // simulation index
  uint i_t_f; // trajectory file index

  // create log file in current working directory
  time_t t_s = time(nullptr); // starting time
  pathstr = std::to_string(t_s) + ".log";
  mmc::logger::set_file(pathstr);

  // main try block
  try
  {
    // read parameters and initialize simulation
    pathstr = sim_dir + "/adjustable-parameters.dat";
    inp_f.open(pathstr);
    mmc::check_file(inp_f, pathstr);
    mmc::parmap par(inp_f); // parameters
    inp_f.close();
    mmc::chrsim sim(par); // simulation

    if (new_sim) // begin new simulation
    {
      // set i_s and i_t_f
      pathpat = sim_dir + "/initial-condition-*";
      i_s = mmc::glob_count(pathpat);
      i_t_f = 0;

      // generate and write initial condition
      sim.generate_initial_condition();
      pathstr = sim_dir + "/initial-condition-";
      pathstr += mmc::cnfs(i_s, 3, '0') + ".gro";
      out_f.open(pathstr);
      mmc::check_file(out_f, pathstr);
      sim.write_frame_txt(out_f);
      out_f.close();

      // write lamina binding sites
      pathstr = sim_dir + "/lamina-binding-sites-";
      pathstr += mmc::cnfs(i_s, 3, '0') + ".gro";
      out_f.open(pathstr);
      mmc::check_file(out_f, pathstr);
      sim.write_lbs_txt(out_f);
      out_f.close();
    }
    else // continue previous simulation
    {
      // set i_s and i_t_f
      i_s = std::stoi(argv[2]);
      pathpat = sim_dir + "/trajectory-";
      pathpat += mmc::cnfs(i_s, 3, '0') + "*";
      i_t_f = mmc::glob_count(pathpat);

      // load checkpoint
      pathstr = sim_dir + "/checkpoint-";
      pathstr += mmc::cnfs(i_s, 3, '0') + ".bin";
      inp_f.open(pathstr, std::ios::binary);
      mmc::check_file(inp_f, pathstr);
      sim.load_checkpoint(inp_f);
      inp_f.close();
    }

    // record i_s and i_t_f
    std::string msg = ""; // message
    msg += "i_s = " + mmc::cnfs(i_s, 3, '0') + " ";
    msg += "i_t_f = " + mmc::cnfs(i_t_f, 3, '0') + " ";
    mmc::logger::record(msg);

    // run simulation
    pathstr = sim_dir + "/trajectory-";
    pathstr += mmc::cnfs(i_s, 3, '0') + "-";
    pathstr += mmc::cnfs(i_t_f, 3, '0') + ".trr";
    out_f.open(pathstr, std::ios::binary);
    mmc::check_file(out_f, pathstr);
    sim.run_simulation(out_f);
    out_f.close();

    // save checkpoint
    pathstr = sim_dir + "/checkpoint-";
    pathstr += mmc::cnfs(i_s, 3, '0') + ".bin";
    out_f.open(pathstr, std::ios::binary);
    mmc::check_file(out_f, pathstr);
    sim.save_checkpoint(out_f);
    out_f.close();
  }
  catch (const mmc::error &err) // caught error
  {
    // exit program unsuccessfully
    mmc::logger::record(err.what());
    return EXIT_FAILURE;
  }

  // remove log file
  mmc::logger::set_file("/dev/null");
  pathstr = std::to_string(t_s) + ".log";
  std::remove(pathstr.c_str());

  // exit program successfully
  return EXIT_SUCCESS;
}

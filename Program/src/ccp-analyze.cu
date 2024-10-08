// Includes

#include "chrana.cuh" // chromatin analysis

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
  if (argc > 2)
  {
    std::cout << "extra arguments\n";
    return EXIT_FAILURE;
  }

  // declare auxiliary variables
  const std::string sim_dir = argv[1]; // simulation directory
  std::ifstream inp_f; // input file
  std::ofstream txt_out_f; // text output file
  std::ofstream bin_out_f; // binary output file
  std::string pathstr; // file path string
  std::string pathpat; // file path pathpat
  uint n_s; // number of simulations
  uint n_t_f; // number of trajectory files
  std::string msg; // message

  // create log file in current working directory
  time_t t_s = time(nullptr); // starting time
  pathstr = std::to_string(t_s) + ".log";
  mmc::logger::set_file(pathstr);

  // main try block
  try
  {
    // read parameters and initialize analysis
    pathstr = sim_dir + "/adjustable-parameters.dat";
    inp_f.open(pathstr);
    mmc::check_file(inp_f, pathstr);
    mmc::parmap par(inp_f); // parameters
    inp_f.close();
    mmc::chrana ana(par); // analysis

    // find the number of simulations
    pathpat = sim_dir + "/initial-condition-*";
    n_s = mmc::glob_count(pathpat);

    // analyse all simulations
    for (uint i_s = 0; i_s < n_s; ++i_s) // simulation index
    {
      // add initial condition to analysis
      pathstr = sim_dir + "/initial-condition-";
      pathstr += mmc::cnfs(i_s, 3, '0') + ".gro";
      inp_f.open(pathstr);
      mmc::check_file(inp_f, pathstr);
      ana.add_initial_condition(inp_f);
      inp_f.close();

      // find the number of trajectory files
      pathpat = sim_dir + "/trajectory-";
      pathpat += mmc::cnfs(i_s, 3, '0') + "*";
      n_t_f = mmc::glob_count(pathpat);

      // add trajectory to analysis
      for (uint i_t_f = 0; i_t_f < n_t_f; ++i_t_f) // trajectory file index
      {
        // add trajectory file to analysis
        pathstr = sim_dir + "/trajectory-";
        pathstr += mmc::cnfs(i_s, 3, '0') + "-";
        pathstr += mmc::cnfs(i_t_f, 3, '0') + ".trr";
        inp_f.open(pathstr, std::ios::binary);
        mmc::check_file(inp_f, pathstr);
        ana.add_trajectory_file(inp_f);
        inp_f.close();

        // record success message
        msg = "added trajectory file ";
        msg += mmc::cnfs(i_s, 3, '0') + "-";
        msg += mmc::cnfs(i_t_f, 3, '0');
        mmc::logger::record(msg);
      }

      // calculate last individual simulation statistics
      ana.calc_last_is_stat();

      // save last individual simulation statistics
      pathstr = sim_dir + "/analysis-";
      pathstr += mmc::cnfs(i_s, 3, '0') + ".dat";
      txt_out_f.open(pathstr);
      mmc::check_file(txt_out_f, pathstr);
      ana.save_last_is_stat(txt_out_f);
      txt_out_f.close();

      // clear individual simulation variables
      ana.clear_is_var();
    }

    // calculate combined simulations final statistics
    ana.calc_cs_final_stat();

    // save combined simulations final statistics
    pathstr = sim_dir + "/analysis-fin.dat";
    txt_out_f.open(pathstr);
    mmc::check_file(txt_out_f, pathstr);
    ana.save_cs_final_stat(txt_out_f);
    txt_out_f.close();

    // save contact map average values to binary file
    pathstr = sim_dir + "/contact-map.bin";
    bin_out_f.open(pathstr, std::ios::binary);
    mmc::check_file(bin_out_f, pathstr);
    ana.save_cm_avg_bin(bin_out_f);
    bin_out_f.close();
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

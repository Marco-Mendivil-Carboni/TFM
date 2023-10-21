//Includes

#include "chrana.cuh" //chromatin analysis

#include <iostream> //standard input/output stream objects

//Functions

//main function
int main(
  const int argc, //argument count
  const char **argv) //argument vector
{
  //check command-line arguments
  if (argc<3){ std::cout<<"missing arguments"; return EXIT_FAILURE;}
  if (argc>3){ std::cout<<"extra arguments"; return EXIT_FAILURE;}

  //declare auxiliary variables
  const std::string sim_dir = argv[1]; //simulation directory
  const uint i_s = std::stoi(argv[2]); //simulation index
  std::ifstream inp_f; //input file
  std::ofstream out_f; //output file
  std::string pathstr; //file path string
  std::string pathpat; //file path pathpat
  uint n_t_f; //number of trajectory files

  //main try block
  try
  {
    //read parameters and initialize analysis
    pathstr = sim_dir+"/adjustable-parameters.dat";
    inp_f.open(pathstr);
    mmc::check_file(inp_f,pathstr);
    mmc::parmap par(inp_f); //parameters
    inp_f.close();
    mmc::chrana ana(par); //analysis

    //add initial condition to analysis
    pathstr = sim_dir+"/initial-condition-";
    pathstr += mmc::cnfs(i_s,3,'0')+".gro";
    inp_f.open(pathstr);
    mmc::check_file(inp_f,pathstr);
    ana.add_initial_condition(inp_f);
    inp_f.close();

    //find the number of trajectory files
    pathpat = sim_dir+"/trajectory-";
    pathpat += mmc::cnfs(i_s,3,'0')+"*";
    n_t_f = mmc::glob_count(pathpat);

    //add trajectory to analysis
    for (uint i_t_f = 0; i_t_f<n_t_f; ++i_t_f) //trajectory file index
    {
      //add trajectory file to analysis
      pathstr = sim_dir+"/trajectory-";
      pathstr += mmc::cnfs(i_s,3,'0')+"-";
      pathstr += mmc::cnfs(i_t_f,3,'0')+".trr";
      inp_f.open(pathstr,std::ios::binary);
      mmc::check_file(inp_f,pathstr);
      ana.add_trajectory_file(inp_f);
      inp_f.close();
    }

    //calculate statistics
    ana.calc_stats();

    //save analysis results
    pathstr = sim_dir+"/analysis-";
    pathstr += mmc::cnfs(i_s,3,'0')+".dat";
    out_f.open(pathstr);
    mmc::check_file(out_f,pathstr);
    ana.save_results(out_f);
    inp_f.close();
  }
  catch (const mmc::error &err) //caught error
  {
    //exit program unsuccessfully
    mmc::logger::record(err.what());
    mmc::logger::record("program exited unsuccessfully");
    return EXIT_FAILURE;
  }

  //exit program successfully
  mmc::logger::record("program exited successfully");
  return EXIT_SUCCESS;
}

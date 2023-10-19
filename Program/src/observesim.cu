//Includes

#include "chrvis.cuh" //chromatin visualization

#include <iostream> //standard input/output stream objects
#include <filesystem> //filesystem manipulation functions

//Functions

//main function
int main(
  const int argc, //argument count
  const char **argv) //argument vector
{
  //check command-line arguments
  if (argc<4){ std::cout<<"missing arguments"; return EXIT_FAILURE;}
  if (argc>4){ std::cout<<"extra arguments"; return EXIT_FAILURE;}

  //declare auxiliary variables
  const std::string sim_dir = argv[1]; //simulation directory
  const uint i_s = std::stoi(argv[2]); //simulation index
  const std::string out_dir = argv[3]; //output directory
  std::ifstream inp_f; //input file
  std::string pathstr; //file path string
  std::string pathpat; //file path pathpat
  uint n_t_f; //number of trajectory files

  //main try block
  try
  {
    //create output directory
    std::filesystem::create_directory(out_dir);

    //read parameters and initialize visualization
    pathstr = sim_dir+"/adjustable-parameters.dat";
    inp_f.open(pathstr);
    mmc::check_file(inp_f,pathstr);
    mmc::parmap par(inp_f); //parameters
    inp_f.close();
    mmc::chrvis vis(par); //visualization

    //read and render initial condition
    pathstr = sim_dir+"/initial-condition-";
    pathstr += mmc::cnfs(i_s,3,'0')+".gro";
    inp_f.open(pathstr);
    mmc::check_file(inp_f,pathstr);
    vis.read_frame_txt(inp_f);
    // vis.render_frame(out_dir);
    inp_f.close();

    //find the number of trajectory files
    pathpat = sim_dir+"/trajectory-";
    pathpat += mmc::cnfs(i_s,3,'0')+"*";
    n_t_f = mmc::glob_count(pathpat);

    //read and render trajectory
    for (uint i_t_f = 0; i_t_f<n_t_f; ++i_t_f) //trajectory file index
    {
      //read and render trajectory file
      pathstr = sim_dir+"/trajectory-";
      pathstr += mmc::cnfs(i_s,3,'0')+"-";
      pathstr += mmc::cnfs(i_t_f,3,'0')+".trr";
      inp_f.open(pathstr,std::ios::binary);
      mmc::check_file(inp_f,pathstr);
      while (!inp_f.eof())
      {
        vis.read_frame_txt(inp_f);
        // vis.render_frame(out_dir);
      }
      inp_f.close();
    }
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

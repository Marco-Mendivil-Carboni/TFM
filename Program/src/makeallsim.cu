//Includes

#include "chrdat.cuh" //chromatin data

#include <iostream> //standard input/output stream objects
#include <filesystem> //facilities for performing operations on file systems

//Functions

//perform missing simulations and analyze them
void perform_missing_sim(
  std::string sim_dir, //simulation directory
  uint n_s, //number of simulations
  uint fps) //files per simulation
{
  //declare auxiliary variables
  std::string pathpat; //file path pattern
  std::string cmd; //command
  bool ns = false; //new simulations

  //create n_s simulations if they don't exist
  pathpat = sim_dir+"/initial-condition-*";
  while (mmc::glob_count(pathpat)<n_s)
  {
    cmd = "./Program/bin/performsim";
    cmd += " "+sim_dir;
    std::system(cmd.c_str());
    ns = true;
  }

  //iterate over all simulations
  for (uint i_s = 0; i_s<n_s; ++i_s) //simulation index
  {
    //generate fps simulation files if they don't exist
    pathpat = sim_dir+"/trajectory-";
    pathpat += mmc::cnfs(i_s,3,'0')+"*";
    while (mmc::glob_count(pathpat)<fps)
    {
      cmd = "./Program/bin/performsim";
      cmd += " "+sim_dir;
      cmd += " "+mmc::cnfs(i_s,3,'0');
      std::system(cmd.c_str());
      ns = true;
    }
  }

  //analyze simulations
  pathpat = sim_dir+"/analysis-*";
  if (ns||mmc::glob_count(pathpat)==0)
  {
    cmd = "./Program/bin/analyzesim";
    cmd += " "+sim_dir;
    std::system(cmd.c_str());
  }
}

//main function
int main(
  const int argc, //argument count
  const char **argv) //argument vector
{
  //check command-line arguments
  if (argc<2){ std::cout<<"no arguments\n"; return EXIT_FAILURE;}
  if (argc>3){ std::cout<<"extra arguments\n"; return EXIT_FAILURE;}

  //declare auxiliary variables
  const std::string srd = "Simulations"; //simulation root directory
  float cvf; //chromatin volume fraction
  float laf; //lbs area fraction
  uint N; //number of particles
  float R; //confinement radius
  uint n_l; //number of lbs
  std::string sim_dir; //simulation directory
  std::ofstream par_f; //parameter file
  std::string pathstr; //file path string
  uint n_s = std::stoi(argv[1]); //number of simulations
  uint fps = (argc==3)?std::stoi(argv[2]):1; //files per simulation

  //main try block
  try
  {
    //iterate over all simulation configurations
    for (N = 4'096; N<65'536; N*=2)
    {
      for (cvf = 0.10; cvf<0.45; cvf+=0.10)
      {
        for (laf = 0.05; laf<0.55; laf+=0.15)
        {
          //create simulation directory if it doesn't exist
          sim_dir = srd;
          sim_dir += "/"+mmc::cnfs(N,5,'0');
          sim_dir += "-"+mmc::cnfs(cvf,5,'0',3);
          sim_dir += "-"+mmc::cnfs(laf,5,'0',3);
          std::filesystem::create_directory(sim_dir);

          //calculate parameters
          R = 0.5+0.5*pow(N/cvf,1.0/3);
          n_l = laf*4.0/pow(mmc::lco/(R-mmc::rco),2.0);

          //write parameter file
          pathstr = sim_dir+"/adjustable-parameters.dat";
          par_f.open(pathstr);
          mmc::check_file(par_f,pathstr);
          par_f<<"number_of_particles "<<mmc::cnfs(N,5,'0')<<"\n";
          par_f<<"confinement_radius "<<mmc::cnfs(R,5,'0',2)<<"\n";
          par_f<<"number_of_lbs "<<mmc::cnfs(n_l,5,'0')<<"\n";
          par_f.close();

          //perform missing simulations and analyze them
          perform_missing_sim(sim_dir,n_s,fps);
        }
      }
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

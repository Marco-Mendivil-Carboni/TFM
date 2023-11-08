//Includes

#include "chrdat.cuh" //chromatin data

#include <iostream> //standard input/output stream objects
#include <filesystem> //facilities for performing operations on file systems

//Functions

//perform and analyze missing simulations
void perform_and_analyze_sim(
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
  uint n_s = std::stoi(argv[1]); //number of simulations
  uint fps = (argc==3)?std::stoi(argv[2]):1; //files per simulation
  std::string sim_dir; //simulation directory
  std::string pathstr; //file path string
  std::ofstream par_f; //parameter file
  std::ifstream inp_f; //input file
  std::ofstream out_f; //output file
  float cvf; //chromatin volume fraction
  float laf; //lbs area fraction
  uint N; //number of particles
  float R; //confinement radius
  uint n_l; //number of lbs

  //open log file inside simulation root directory
  mmc::logger::set_file(srd+"/all-messages.log");

  //record n_s and fps
  std::string msg = ""; //message
  msg += "n_s = "+mmc::cnfs(n_s,3,'0')+" ";
  msg += "fps = "+mmc::cnfs(fps,3,'0')+" ";
  mmc::logger::record(msg);

  //main try block
  try
  {
    //open output file
    pathstr = srd+"/analysis-summary.dat";
    out_f.open(pathstr);
    mmc::check_file(out_f,pathstr);

    //write analysis summary header
    out_f<<"#   N   cvf   laf ";
    out_f<<"#dcm:    avg   sqrt(var)         sem ";
    out_f<<"#rg2:    avg   sqrt(var)         sem ";
    out_f<<"#nop:    avg   sqrt(var)         sem\n";

    //iterate over all simulation configurations
    for (N = 4'096; N<4097; N*=2)
    {
      for (cvf = 0.10; cvf<0.25; cvf+=0.10)
      {
        for (laf = 0.05; laf<0.15; laf+=0.15)
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

          //perform and analyze missing simulations
          perform_and_analyze_sim(sim_dir,n_s,fps);

          //read simulation analysis and write analysis summary
          pathstr = sim_dir+"/analysis-fin.dat";
          inp_f.open(pathstr);
          mmc::check_file(inp_f,pathstr);
          std::string fl; //file line
          out_f<<mmc::cnfs(N,5,'0')<<" ";
          out_f<<mmc::cnfs(cvf,5,'0',3)<<" ";
          out_f<<mmc::cnfs(laf,5,'0',3)<<" ";
          getline(inp_f,fl); getline(inp_f,fl); getline(inp_f,fl);
          getline(inp_f,fl); out_f<<fl<<" "; getline(inp_f,fl);
          getline(inp_f,fl); out_f<<fl<<" "; getline(inp_f,fl);
          getline(inp_f,fl); out_f<<fl<<"\n";
          inp_f.close();

          //record success message
          mmc::logger::record(sim_dir+" done");
        }
      }
    }

    //close output file
    out_f.close();
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

#ifndef MMC_CHROMATIN_H
#define MMC_CHROMATIN_H

//Includes

#include <fstream> //file stream classes

#include <curand_kernel.h> //cuRAND device functions

//Namespace

namespace mmc //Marco Mendívil Carboni
{

//Aliases

using prng = curandStatePhilox4_32_10; //PRNG type

//Classes

class chrsim //chromatin simulation
{
  public:

  //Functions

  //chrsim constructor
  chrsim(std::ifstream &f_par); //parameter file

  //chrsim destructor
  ~chrsim();

  //generate a random initial condition
  void generate_initial_condition();

  //write initial condition to file in gro format
  void write_initial_condition(std::ofstream &f_ic); //IC file

  //save simulation state to binary file
  void save_checkpoint(std::ofstream &f_chkp); //checkpoint file

  //load simulation state from binary file
  void load_checkpoint(std::ifstream &f_chkp); //checkpoint file

  //run simulation and write trajectory file
  void run_simulation(std::ofstream &f_traj); //trajectory file

  private:

  //Parameters and Variables

  struct //adjustable parameters
  {
    uint N; //number of particles
    float T; //temperature
    float R; //confinement radius
    uint f_f; //frames per file
    uint f_s = 1*2048; //steps per frame
  } ap;

  uint thd_blk = 256; //threads per block
  uint n_p_blk; //number of particle blocks

  uint i_f = 0; //frame index
  float t = 0.0; //simulation time
  float sig = 1.0; //particle LJ size

  float4 *r_2; //position array 2
  float4 *r_1; //position array 1

  float4 *f_2; //force array 2
  float4 *f_1; //force array 1

  float sd; //random number standard deviation
  float4 *n_r; //normal random numbers

  prng *state; //device PRNG state array

  //Functions

  //read adjustable parameters from file
  void read_parameters(std::ifstream &f_par); //parameter file

  //make one Runge-Kutta iteration
  void make_RK_iteration();

  //write trajectory frame to binary file in trr format
  void write_trajectory_frame(std::ofstream &f_traj); //trajectory file
};

//Functions

//check for errors in cuda runtime API call
void cuda_check(cudaError_t rtn_val); //cuda runtime API call return value

} //namespace mmc

#endif //MMC_CHROMATIN_H

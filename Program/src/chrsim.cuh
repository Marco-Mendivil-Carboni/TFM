#ifndef MMCC_CHRSIM_H
#define MMCC_CHRSIM_H

//Includes

#include <fstream> //file stream classes

#include <curand_kernel.h> //cuRAND device functions

//Namespace

namespace mmcc //Marco Mendívil Carboni code
{

//Aliases

using prng = curandStatePhilox4_32_10; //PRNG type

//Classes

class chrsim //chromatin simulation
{
  public:

  //Parameters and Variables

  struct //adjustable parameters
  {
    float T; //temperature
    int N; //number of particles
    float R; //radius of sphere
    int F; //frames per file
  } ap;

  size_t thd_blk = 256; //threads per block
  size_t n_p_blk; //number of particle blocks
  size_t n_p_thd; //number of particle threads

  int i_f = 0; //frame index
  float t = 0.0; //simulation time

  float sig = 1.0; //LJ particle size

  float c_rn; //random number constant

  float4 *r_2; //positions 2
  float4 *r_1; //positions 1

  float4 *f_2; //forces 2
  float4 *f_1; //forces 1

  float4 *nrn; //normal random numbers

  prng *state; //device PRNG state

  //Functions

  //chrsim constructor
  chrsim(std::ifstream &f_par); //parameter file

  //chrsim destructor
  ~chrsim();

  //generate a random initial condition
  void generate_initial_condition();

  //write initial condition to file in gro format
  void write_initial_condition(std::ofstream &f_i_c); //initial condition file

  //save simulation state to binary file
  void save_checkpoint(std::ofstream &f_chkp); //checkpoint file

  //load simulation state from binary file
  void load_checkpoint(std::ifstream &f_chkp); //checkpoint file

  //run simulation and write trajectory file
  void run_simulation(std::ofstream &f_traj); //trajectory file

  private:

  //Functions

  //read adjustable parameters from file
  void read_parameters(std::ifstream &f_par); //parameter file

  //take RK step----------------------------------------------------------------tmp
  void take_step();//-----------------------------------------------------------tmp

  //write trajectory frame to binary file in trr format
  void write_trajectory_frame(std::ofstream &f_traj); //trajectory file
};

//Functions

//check for errors in cuda runtime API call
void cuda_check(cudaError_t rtn_val); //cuda runtime API call return value

} //namespace mmcc

#endif //MMCC_CHRSIM_H

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
  chrsim(std::ifstream &f_par);

  //chrsim destructor
  ~chrsim();

  //generate a random initial condition
  void generate_initial_condition();

  //write initial condition to file in gro format
  void write_initial_condition(std::ofstream &f_i_c);

  //save simulation state to binary file
  void save_checkpoint(std::ofstream &f_chkp);

  //load simulation state from binary file
  void load_checkpoint(std::ifstream &f_chkp);

  //write trajectory to binary file in trr format
  void write_trajectory(std::ofstream &f_traj, int i_f);

  private:

  //Functions

  //read adjustable parameters from file
  void read_parameters(std::ifstream &f_par);
};

//Functions

//check for errors in cuda runtime API call
void cuda_check(cudaError_t result);

} //namespace mmcc

#endif //MMCC_CHRSIM_H

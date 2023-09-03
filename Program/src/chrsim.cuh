#ifndef MMCC_CHRSIM_H
#define MMCC_CHRSIM_H

//Includes

#include <fstream> //file stream classes

#include <curand_kernel.h> //cuRAND device functions

//Namespace

namespace mmcc //Marco Mendívil Carboni code
{

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
  using PRNGstate = curandStatePhilox4_32_10; //PRNG state alias

  float4 *r_2; //positions 2
  float4 *r_1; //positions 1

  float4 *f_2; //forces 2
  float4 *f_1; //forces 1

  float4 *nrn; //normal random numbers

  PRNGstate *state; //PRNG state array

  //Functions

  //chrsim constructor
  chrsim(std::ifstream &f_par);

  //chrsim destructor
  ~chrsim();

  //generate a random initial configuration of chromatin
  void generate_initial_configuration();

  //write initial configuration to file in gro format
  void write_initial_configuration(std::ofstream &f_out);

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

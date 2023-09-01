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

  float4 *r_2; //positions 2
  float4 *r_1; //positions 1

  float4 *f_2; //forces 2
  float4 *f_1; //forces 1

  float c_rn; //random number constant
  float4 *nrn; //normal random numbers

  using PRNGstate = curandStatePhilox4_32_10; //PRNG state alias
  PRNGstate *state; //PRNG state array

  float sig; //LJ particle size

  size_t threads_block = 256; //threads per block
  size_t n_blocks; //number of blocks
  size_t n_threads; //nuumber of threads

  //Functions

  //chrsim constructor
  chrsim(FILE *f_ptr_par);

  //chrsim destructor
  ~chrsim();

  //generate a random initial configuration of chromatin
  void generate_initial_configuration();

  //write initial configuration to file in gro format
  void write_initial_configuration(FILE *f_ptr);

  private:

  //Functions

  //read adjustable parameters from file
  void read_parameters(FILE *f_ptr_par);
};

//Functions

//check for errors in cuda runtime API call
void cuda_check(cudaError_t result);

} //namespace mmcc

#endif //MMCC_CHRSIM_H

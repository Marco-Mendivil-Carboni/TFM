#ifndef MMCC_CHRSIM_H
#define MMCC_CHRSIM_H

//Includes

#include "cudautil.cuh"

#include <curand_kernel.h> //cuRAND device functions

//Namespace

namespace mmcc //Marco Mendívil Carboni code
{

//Classes

class chrsim //chromatin simulation
{
  public:

  //Constructor and Destructor

  chrsim(FILE *f_ptr_par);

  ~chrsim();

  //Parameters and Variables

  struct //adjustable parameters
  {
    float T; //temperature
    int N; //number of particles
    float R; //radius of sphere
    int F; //frames per file
  } ap;

  soa3 r_2; //positions 2
  soa3 r_1; //positions 1

  soa3 f_2; //forces 2
  soa3 f_1; //forces 1

  float c_rn; //random number constant
  soa3 nrn; //random numbers

  using PRNGstate = curandStatePhilox4_32_10; //PRNG state alias
  PRNGstate *state; //PRNG state array

  float sig; //LJ particle size

  size_t threads_block = 256; //threads per block
  size_t n_blocks; //number of blocks
  size_t n_threads; //nuumber of threads

  //Host Functions

  void generate_initial_configuration();

  void write_initial_configuration(FILE *f_ptr);

  private:

  //Host Functions

  void read_parameters(FILE *f_ptr_par);

  //Device Functions
};

} //namespace mmcc

#endif //MMCC_CHRSIM_H

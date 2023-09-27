#ifndef MMC_CHRSIM_H
#define MMC_CHRSIM_H

//Includes

#include "chrdat.cuh" //chromatin data
#include "sugrid.cuh" //sorted uniform grid

#include <curand_kernel.h> //cuRAND device functions

//Namespace

namespace mmc //Marco Mendívil Carboni
{

//Aliases

using prng = curandStatePhilox4_32_10; //PRNG type

//Classes

class chrsim : public chrdat //chromatin simulation
{
  public:

  //Functions

  //chromatin simulation constructor
  chrsim(parmap &par); //parameters

  //chromatin simulation destructor
  ~chrsim();

  //generate a random initial condition
  void generate_initial_condition();

  //save simulation state to binary file
  void save_checkpoint(std::ofstream &bin_out_f); //binary output file

  //load simulation state from binary file
  void load_checkpoint(std::ifstream &bin_inp_f); //binary input file

  //run simulation and write trajectory to binary file
  void run_simulation(std::ofstream &bin_out_f); //binary output file

  private:

  //Parameters and Variables

  const uint fpf; //frames per file
  const uint spf; //steps per frame
  const uint tpb; //threads per block

  float4 *er; //extra position array
  float4 *ef; //extra force array

  const float sd; //standard deviation
  float4 *rn; //random number array
  prng *ps; //PRNG state array

  sugrid ljg; //LJ grid
  sugrid *ljp; //LJ grid pointer

  //Functions

  //set random particle types
  void set_particle_types(curandGenerator_t &gen); //host PRNG

  //perform a confined random walk
  void perform_random_walk(curandGenerator_t &gen); //host PRNG

  //make one Runge-Kutta iteration
  void make_RK_iteration();
};

} //namespace mmc

#endif //MMC_CHRSIM_H

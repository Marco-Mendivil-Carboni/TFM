#ifndef MMC_CHRSIM_H
#define MMC_CHRSIM_H

//Includes

#include "chrdat.cuh" //chromatin data

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

  //chrsim constructor
  chrsim(parmap &par); //parameters

  //chrsim destructor
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

  const int framepf; //frames per file
  const int spframe; //steps per frame

  const int thd_blk; //threads per block
  const int n_p_blk; //number of particle blocks

  int *dpt; //device particle type array
  float4 *dr2; //device position array 2
  float4 *dr1; //device position array 1
  float4 *df2; //device force array 2
  float4 *df1; //device force array 1

  float sd; //random number standard deviation
  float4 *drn; //device random number array

  prng *dps; //device PRNG state array

  //Functions

  //make one Runge-Kutta iteration
  void make_RK_iteration();
};

} //namespace mmc

#endif //MMC_CHRSIM_H

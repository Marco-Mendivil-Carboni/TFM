#ifndef MMC_CHRSIM_H
#define MMC_CHRSIM_H

//Includes

#include "chrdat.cuh" //chromatin data

#include <curand_kernel.h> //cuRAND device functions

#include <cub/device/device_radix_sort.cuh> //cub parallel radix sort

//Namespace

namespace mmc //Marco Mendívil Carboni
{

//Aliases

using prng = curandStatePhilox4_32_10; //PRNG type

using sa = cub::DeviceRadixSort; //sorting algorithm

//Structures

struct sugrid //sorted uniform grid
{
  //Variables

  const float csl; //grid cell side length
  const uint cps; //grid cells per side

  uint *uci; //unsorted grid cell index array
  uint *sci; //sorted grid cell index array
  uint *upi; //unsorted particle index array
  uint *spi; //sorted particle index array

  uint *beg; //grid cell beginning array
  uint *end; //grid cell end array

  float4 *sr; //sorted position array

  // cudaTextureObject_t srt; //sorted position texture

  void *eb; //extra buffer
  size_t ebs; //extra buffer size

  //Functions

  //sorted uniform grid constructor
  sugrid(
    const uint N, //number of particles
    const float csl, //grid cell side length
    const uint cps); //grid cells per side

  //sorted uniform grid destructor
  ~sugrid();
};

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

  const int fpf; //frames per file
  const int spf; //steps per frame
  const int tpb; //threads per block

  float4 *er; //extra position array
  float4 *ef; //extra force array

  const float sd; //random number standard deviation
  float4 *rn; //random number array
  prng *ps; //PRNG state array

  sugrid ljg; //LJ grid
  sugrid *ljp; //LJ grid pointer
  const uint ljc; //number of LJ grid cells

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

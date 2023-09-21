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

//Structures

struct sugrid : dev_s //sorted uniform grid
{
  //Variables

  const float csl; //grid cell side length
  const int n_cps; //number of grid cells per side

  int *i_c[2]; //grid cell index arrays
  int *i_p[2]; //particle index arrays

  int *beg; //grid cell beginning array
  int *end; //grid cell end array

  void *eb; //extra buffer
  size_t ebs; //extra buffer size

  //Functions

  //sugrid constructor
  sugrid(
    const int N, //number of particles
    const float csl, //grid cell side length
    const int n_cps); //number of grid cells per side

  //sugrid destructor
  ~sugrid();
};

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

  const int thdpblk; //threads per block
  const int n_blk; //number of blocks

  float4 *er; //extra position array
  float4 *ef; //extra force array

  const float sd; //random number standard deviation
  float4 *rn; //random number array
  prng *ps; //PRNG state array

  sugrid *gp; //grid pointer

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

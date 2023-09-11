#ifndef MMC_CHROMATIN_H
#define MMC_CHROMATIN_H

//Includes

#include "utilities.hpp"

#include <curand_kernel.h> //cuRAND device functions

//Namespace

namespace mmc //Marco Mendívil Carboni
{

//Aliases

using prng = curandStatePhilox4_32_10; //PRNG type

//Classes

class chrdat //chromatin data
{
  public:

  //Functions

  //chrdat constructor
  chrdat(parmap &par); //parameters

  //chrdat destructor
  ~chrdat();

  protected:

  //Parameters and Variables

  const int N; //number of particles
  const float R; //confinement radius
  const float T; //temperature

  float t; //time

  int *pt; //particle type array
  float4 *r; //position array
  float4 *f; //force array

  float sig; //particle LJ size
};

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

  const int f_f; //frames per file
  const int s_f; //steps per frame

  const int thd_blk; //threads per block
  const int n_p_blk; //number of particle blocks

  int i_f = 0; //frame index//tmp------------------------------------------------

  int *dpt; //device particle type array

  float4 *r_2; //device position array 2
  float4 *r_1; //device position array 1

  float4 *f_2; //device force array 2
  float4 *f_1; //device force array 1

  float sd; //random number standard deviation
  float4 *n_r; //device random number array

  prng *state; //device PRNG state array

  //Functions

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

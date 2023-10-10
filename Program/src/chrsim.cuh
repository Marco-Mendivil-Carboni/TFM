#ifndef MMC_CHRSIM_H
#define MMC_CHRSIM_H

//Includes

#include "chrdat.cuh" //chromatin data
#include "sugrid.cuh" //sorted uniform grid

//Namespace

namespace mmc //Marco Mendívil Carboni
{

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
  void *vps; //void PRNG state array

  sugrid pg; //particle grid
  sugrid *pgp; //particle grid pointer

  sugrid lg; //lbs grid
  sugrid *lgp; //lbs grid pointer

  // chrsim *sdp; //simulation device pointer -----------------------------------?

  //Functions

  //set random lbs positions
  void set_lbs_positions();

  //set random particle types
  void set_particle_types();

  //perform a confined random walk
  void perform_random_walk();

  //count particle overlaps
  uint particle_overlaps();
};

} //namespace mmc

#endif //MMC_CHRSIM_H

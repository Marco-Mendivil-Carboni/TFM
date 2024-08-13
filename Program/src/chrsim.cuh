#ifndef MMC_CHRSIM_H
#define MMC_CHRSIM_H

// Includes

#include "chrdat.cuh" // chromatin data
#include "sugrid.cuh" // sorted uniform grid

#include <curand_kernel.h> // cuRAND device functions

// Namespace

namespace mmc // Marco Mendívil Carboni
{

// Constants

static constexpr float dt = 1.0 / 2048; // timestep
static constexpr float mis = 0.838732; // minimum initial separation

// Aliases

using prng = curandStatePhilox4_32_10; // PRNG type

// Enumerations

enum stype // simulation type
{
  DST, // default simulation type
  ICG, // initial condition generation
};

// Classes

class chrsim : public chrdat // chromatin simulation
{
public:
  // Functions

  // chromatin simulation constructor
  chrsim(parmap &par); // parameters

  // chromatin simulation destructor
  ~chrsim();

  // generate a random initial condition
  void generate_initial_condition();

  // save simulation state to binary file
  void save_checkpoint(std::ofstream &bin_out_f); // binary output file

  // load simulation state from binary file
  void load_checkpoint(std::ifstream &bin_inp_f); // binary input file

  // run simulation and write trajectory to binary file
  void run_simulation(std::ofstream &bin_out_f); // binary output file

private:
  // Parameters and Variables

  const uint spf; // steps per frame
  const uint tpb; // threads per block

  const float sd; // rn standard deviation

  vec3f *er; // extra position array
  vec3f *ef; // extra force array

  vec3f *rn; // random number array
  prng *ps; // PRNG state array

  sugrid pg; // particle grid
  sugrid *pgp; // particle grid pointer

  sugrid lg; // lbs grid
  sugrid *lgp; // lbs grid pointer

  // Functions

  // set random lbs positions
  void set_lbs_positions(curandGenerator_t gen); // host PRNG

  // set particle type sequence
  void set_particle_types(curandGenerator_t gen); // host PRNG

  // set random particle positions
  void set_particle_positons(curandGenerator_t gen); // host PRNG

  // perform confined random walk
  void perform_random_walk(
      curandGenerator_t gen, // host PRNG
      uint i_s, // starting index
      uint i_e, // ending index
      vec3f dir); // direction

  // count particle overlaps
  uint particle_overlaps();
};

} // namespace mmc

#endif // MMC_CHRSIM_H

#ifndef MMC_CHRDAT_H
#define MMC_CHRDAT_H

// Includes

#include "util.cuh" // general utilities

// Namespace

namespace mmc // Marco Mendívil Carboni
{

// Constants

static constexpr float k_B = 0.003356; // Boltzmann constant
static constexpr float l_0 = 1.000000; // bond natural length
static constexpr float k_e = 128.0000; // elastic constant
static constexpr float k_b = 6.000000; // bending constant
static constexpr float e_p = 0.500000; // particle energy
static constexpr float rco = 1.154701; // repulsive cutoff
static constexpr float aco = 2.000000; // attractive cutoff
static constexpr float mdl = 22.13471; // mean domain length
static constexpr float e_l = 8.000000; // lbs energy
static constexpr float lco = 0.500000; // lbs cutoff

// Enumerations

enum ptype // particle type
{
  LADh, // lamina associated domain heterochromatin-like (A)
  LNDe, // lamina non-associated domain euchromatin-like (B)
};

// Structures

struct cngeom // cell nucleus geometry
{
  // Variables

  float R_n; // nucleus radius
  float R_o; // opening radius
  float R_b; // bleb radius

  float noc; // nucleus opening cosine
  float boc; // bleb opening cosine

  float nod; // nucleus opening distance
  float bod; // bleb opening distance

  float d_b; // distance to bleb
  float d_m; // maximum distance

  // Functions

  // cell nucleus geometry constructor
  cngeom(parmap &par); // parameters
};

// Classes

class chrdat // chromatin data
{
public:
  // Functions

  // chromatin data constructor
  chrdat(parmap &par); // parameters

  // chromatin data destructor
  ~chrdat();

  // write frame to text file
  void write_frame_txt(std::ofstream &txt_out_f); // text output file

  // read frame from text file
  void read_frame_txt(std::ifstream &txt_inp_f); // text input file

  // write frame to binary file
  void write_frame_bin(std::ofstream &bin_out_f); // binary output file

  // read frame from binary file
  void read_frame_bin(std::ifstream &bin_inp_f); // binary input file

  // write lamina binding sites to text file
  void write_lbs_txt(std::ofstream &txt_out_f); // text output file

  // read lamina binding sites from text file
  void read_lbs_txt(std::ifstream &txt_inp_f); // text input file

protected:
  // Parameters and Variables

  const uint N; // number of particles
  const cngeom ng; // nucleus geometry
  const float T; // temperature
  const uint n_l; // number of lbs

  const uint fpf; // frames per file
  uint i_f; // frame index

  float t; // time

  ptype *pt; // particle type array
  vec3f *r; // position array
  vec3f *f; // force array
  vec3f *lr; // lbs position array

  ptype *hpt; // host particle type array
  vec3f *hr; // host position array
  vec3f *hf; // host force array
  vec3f *hlr; // host lbs position array
};

} // namespace mmc

#endif // MMC_CHRDAT_H

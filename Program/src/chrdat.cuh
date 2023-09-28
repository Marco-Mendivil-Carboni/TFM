#ifndef MMC_CHRDAT_H
#define MMC_CHRDAT_H

//Includes

#include "util.cuh" //general utilities
#include "vecops.cuh" //vector operations

//Namespace

namespace mmc //Marco Mendívil Carboni
{

//Constants

static constexpr float xi = 1.000000; //damping coefficient
static constexpr float k_B = 0.003356; //Boltzmann constant
static constexpr float l_0 = 1.000000; //bond natural length
static constexpr float k_e = 100.0000; //elastic constant
static constexpr float k_b = 2.000000; //bending constant

//Enumerations

enum ptype //particle type
{
  non_LAD, //non lamina associated domain
  LAD //lamina associated domain
};

//Classes

class chrdat //chromatin data
{
  public:

  //Functions

  //chromatin data constructor
  chrdat(parmap &par); //parameters

  //chromatin data destructor
  ~chrdat();

  //write frame to text file
  void write_frame_txt(std::ofstream &txt_out_f); //text output file

  //read frame from text file
  void read_frame_txt(std::ifstream &txt_inp_f); //text input file

  //write frame to binary file
  void write_frame_bin(std::ofstream &bin_out_f); //binary output file

  //read frame from binary file
  void read_frame_bin(std::ifstream &bin_inp_f); //binary input file

  protected:

  //Parameters and Variables

  const uint N; //number of particles
  const float R; //confinement radius
  const float T; //temperature

  const float eps; //LJ particle energy
  float sig; //LJ particle size

  uint i_f; //frame index
  double t; //time

  ptype *pt; //particle type array
  float4 *r; //position array
  float4 *f; //force array

  ptype *hpt; //host particle type array
  float4 *hr; //host position array
};

} //namespace mmc

#endif //MMC_CHRDAT_H

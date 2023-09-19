#ifndef MMC_CHRDAT_H
#define MMC_CHRDAT_H

//Includes

#include "util.cuh" //general utilities
#include "vecops.cuh" //vector operations

//thrust header files
#include <thrust/sort.h>
#include <thrust/device_vector.h>

#include <thrust/device_ptr.h>

//Namespace

namespace mmc //Marco Mendívil Carboni
{

//Constants

static constexpr float xi = 1.000000; //damping coefficient
static constexpr float k_B = 0.001120; //Boltzmann constant
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

  //chrdat constructor
  chrdat(parmap &par); //parameters

  //chrdat destructor
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

  const int N; //number of particles
  const float R; //confinement radius
  const float T; //temperature

  int i_f; //frame index
  double t; //time

  ptype *pt; //particle type array
  float4 *r; //position array
  float4 *f; //force array

  float sig; //LJ particle size
};

} //namespace mmc

#endif //MMC_CHRDAT_H

#ifndef MMC_CHRDAT_H
#define MMC_CHRDAT_H

//Includes

#include "utilities.hpp" //general utilities

//Namespace

namespace mmc //Marco Mendívil Carboni
{

//Constants

static constexpr float xi  = 1.000000; //damping coefficient
static constexpr float l_0 = 1.000000; //bond natural length
static constexpr float k_e = 100.0000; //elastic constant
static constexpr float k_b = 2.000000; //bending constant
static constexpr float r_c = 1.122462; //LJ cutoff radius//rename this to repulsive cutoff radius (rcr)
//add attractive cutoff radius (acr)

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

  int *pt; //particle type array
  float4 *r; //position array
  float4 *f; //force array

  float sig; //particle LJ size
};

} //namespace mmc

#endif //MMC_CHRDAT_H

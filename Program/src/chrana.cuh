#ifndef MMC_CHRANA_H
#define MMC_CHRANA_H

//Includes

#include "chrdat.cuh" //chromatin data

#include <vector> //vector container class

//Namespace

namespace mmc //Marco Mendívil Carboni
{

//Classes

class chrana : public chrdat //chromatin analysis
{
  public:

  //Functions

  //chromatin analysis constructor
  chrana(parmap &par); //parameters

  //chromatin analysis destructor
  ~chrana();

  //add initial condition to analysis
  void add_initial_condition(std::ifstream &txt_inp_f); //text input file

  //add trajectory file to analysis
  void add_trajectory_file(std::ifstream &bin_inp_f); //binary input file

  private:

  //Parameters and Variables

  const uint fpf; //frames per file

  std::vector<vec3f> cmr_s; //center of mass position series

  std::vector<float> rg2_s; //gyration radius squared series

  //Functions
  
  //center of mass position
  vec3f center_of_mass();

  //gyration radius squared
  float gyration_radius_sq(vec3f cmr); //center of mass position
};

} //namespace mmc

#endif //MMC_CHRANA_H

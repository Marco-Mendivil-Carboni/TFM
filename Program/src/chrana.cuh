#ifndef MMC_CHRANA_H
#define MMC_CHRANA_H

//Includes

#include "chrdat.cuh" //chromatin data

#include <vector> //vector container class

//Namespace

namespace mmc //Marco Mendívil Carboni
{

//Structures

struct ranvar //random variable
{
  //Variables

  std::vector<float> val; //values

  double m_1; //1st moment (mean)
  double m_2; //2nd moment
  double var; //variance
  double sem; //standard error of the mean

  uint i_ter; //termalization index
  uint n_b; //number of blocks

  //Functions

  //include value
  void inc_val(float c_val); //current value

  //calculate statistics
  void calc_stats();

  //check if random variable has termalized
  bool has_termalized();
};

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

  //calculate statistics
  void calc_stats();

  //save analysis results
  void save_results(std::ofstream &txt_out_f); //text output file

  private:

  //Parameters and Variables

  const uint fpf; //frames per file

  ranvar dcm; //distance to the center of mass
  ranvar rg2; //gyration radius squared

  //Functions

  //calculate random variables
  void calc_ranvar();
};

} //namespace mmc

#endif //MMC_CHRANA_H

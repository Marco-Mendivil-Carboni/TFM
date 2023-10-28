#ifndef MMC_CHRANA_H
#define MMC_CHRANA_H

//Includes

#include "chrdat.cuh" //chromatin data

#include <vector> //vector container class

//Namespace

namespace mmc //Marco Mendívil Carboni
{

//Constants

static constexpr uint n_b = 128; //number of bins

//Structures

struct gdstat //generic data statistics
{
  double avg; //average
  double var; //variance
  double sem; //standard error of the mean
};

struct idstat : gdstat {}; //independent data statistics

struct cdstat : gdstat //correlated data statistics
{
  uint f_n_b; //final number of blocks
};

struct tdstat : cdstat //time series data statistics
{
  uint i_t; //termalization index
  bool ter; //termalized
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

  //calculate individual simulation statistics
  void calc_ind_sim_stat();

  //save individual simulation analysis results
  void save_ind_sim_results(std::ofstream &txt_out_f); //text output file

  //clear individual simulation analysis data
  void clear_ind_sim_data();

  //calculate final statistics
  void calc_fin_stat();

  //save final analysis results
  void save_fin_results(std::ofstream &txt_out_f); //text output file

  private:

  //Parameters and Variables

  const uint fpf; //frames per file

  std::vector<float> t_v; //time vector

  std::vector<float> dcm_v; //center of mass distance vector
  std::vector<float> rg2_v; //gyration radius squared vector
  std::vector<float> nop_v; //nematic order parameter vector
  std::vector<float> rcd_v[n_b]; //radial chromatin density vector
  std::vector<float> *msd_v; //mean spatial distance vector

  std::vector<tdstat> dcm_s_v; //dcm statistics vector
  std::vector<tdstat> rg2_s_v; //rg2 statistics vector
  std::vector<tdstat> nop_s_v; //nop statistics vector
  std::vector<tdstat> rcd_s_v[n_b]; //rcd statistics vector
  std::vector<tdstat> *msd_s_v; //msd statistics vector

  idstat dcm_f_s; //dcm final statistics
  idstat rg2_f_s; //rg2 final statistics
  idstat nop_f_s; //nop final statistics
  idstat rcd_f_s[n_b]; //rcd final statistics
  idstat *msd_f_s; //msd final statistics

  //Function

  //calculate observables
  void calc_observables();
};

//Functions

//calculate statistics
void calc_stats(
  const std::vector<float> &v, //vector
  idstat &s); //statistics

//calculate statistics
void calc_stats(
  const std::vector<float> &v, //vector
  cdstat &s); //statistics

//calculate statistics
void calc_stats(
  const std::vector<float> &v, //vector
  tdstat &s); //statistics

//calculate statistics
void calc_stats(
  const std::vector<tdstat> &v, //vector
  idstat &s); //statistics

} //namespace mmc

#endif //MMC_CHRANA_H

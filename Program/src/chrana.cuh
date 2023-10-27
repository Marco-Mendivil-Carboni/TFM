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

struct idstat //independent data statistics
{
  double avg; //average
  double var; //variance
  double sem; //standard error of the mean
};

struct cdstat : idstat //correlated data statistics
{
  uint f_n_b; //final number of blocks
};

struct tdstat : cdstat //time series data statistics
{
  uint i_t; //termalization index
  bool ter; //termalized
};

struct avgsem //average and standard error of the mean
{
  double avg; //average
  double sem; //standard error of the mean
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

  //calculate all simulations statistics
  void calc_all_sim_stat();

  //save all simulations analysis results
  void save_all_sim_results(std::ofstream &txt_out_f); //text output file

  private:

  //Parameters and Variables

  const uint fpf; //frames per file

  std::vector<float> t_v; //time vector

  std::vector<float> dcm_v; //center of mass distance vector
  tdstat dcm_s; //center of mass distance statistics
  std::vector<float> rg2_v; //gyration radius squared vector
  tdstat rg2_s; //gyration radius squared statistics
  std::vector<float> nop_v; //nematic order parameter vector
  tdstat nop_s; //nematic order parameter statistics
  std::vector<float> rcd_v[n_b]; //radial chromatin density vector
  tdstat rcd_s[n_b]; //radial chromatin density statistics

  std::vector<avgsem> dcmrv; //center of mass distance results vector
  idstat dcmrs; //center of mass distance results statistics
  std::vector<avgsem> rg2rv; //gyration radius squared results vector
  idstat rg2rs; //gyration radius squared results statistics
  std::vector<avgsem> noprv; //nematic order parameter results vector
  idstat noprs; //nematic order parameter results statistics
  std::vector<avgsem> rcdrv[n_b]; //radial chromatin density results vector
  idstat rcdrs[n_b]; //radial chromatin density results statistics

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
  const std::vector<avgsem> &v, //vector
  idstat &s); //statistics

} //namespace mmc

#endif //MMC_CHRANA_H

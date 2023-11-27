#ifndef MMC_CHRANA_H
#define MMC_CHRANA_H

//Includes

#include "chrdat.cuh" //chromatin data
#include "sostat.cuh" //simulation observable statistics

//Namespace

namespace mmc //Marco Mendívil Carboni
{

//Constants

static constexpr uint n_b = 128; //number of bins

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

  simobs t_o; //time observable

  simobs dcm_o; //center of mass distance observable
  simobs rg2_o; //gyration radius squared observable
  simobs nop_o; //nematic order parameter observable

  simobs rcd_o[3][n_b]; //radial chromatin density observable

  simobs *msd_o; //mean spatial distance observable
  const uint lma; //length of msd arrays
  float *ma; //msd array
  float *hma; //host msd array

  //Function

  //calculate observables
  void calc_observables();
};

} //namespace mmc

#endif //MMC_CHRANA_H

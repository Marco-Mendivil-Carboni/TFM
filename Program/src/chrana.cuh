#ifndef MMC_CHRANA_H
#define MMC_CHRANA_H

// Includes

#include "chrdat.cuh" // chromatin data
#include "sostat.cuh" // simulation observable statistics

// Namespace

namespace mmc // Marco Mendívil Carboni
{

// Constants

static constexpr uint n_b = 64; // number of bins
static constexpr uint px_sz = 4; // pixel size
static constexpr float cf = 1'000.0; // contact factor

static constexpr uint n_chr_h = n_chr / 2; // haploid number of chromosomes

// Classes

class chrana : public chrdat // chromatin analysis
{
public:
  // Functions

  // chromatin analysis constructor
  chrana(parmap &par); // parameters

  // chromatin analysis destructor
  ~chrana();

  // add initial condition to analysis
  void add_initial_condition(std::ifstream &txt_inp_f); // text input file

  // add trajectory file to analysis
  void add_trajectory_file(std::ifstream &bin_inp_f); // binary input file

  // calculate last individual simulation statistics
  void calc_last_is_stat();

  // save last individual simulation statistics
  void save_last_is_stat(std::ofstream &txt_out_f); // text output file

  // clear individual simulation variables
  void clear_is_var();

  // calculate combined simulations final statistics
  void calc_cs_final_stat();

  // save combined simulations final statistics
  void save_cs_final_stat(std::ofstream &txt_out_f); // text output file

  // save contact map average values to binary file
  void save_cm_avg_bin(std::ofstream &bin_out_f); // binary output file

private:
  // Parameters and Variables

  const uint tpb; // threads per block

  simobs t_o; // time observable

  simobs dcm_o; // center of mass distance observable
  simobs rg2_o; // gyration radius squared observable
  simobs nop_o; // nematic order parameter observable
  simobs ncf_o; // nucleus chromatin fraction observable

  simobs rcd_o[3][n_b]; // radial chromatin density observable

  simobs_b *sd_bo[n_chr_h]; // spatial distance basic observable
  float *sd[n_chr_h]; // spatial distance arrays
  float *hsd[n_chr_h]; // host spatial distance arrays

  simobs_b *cp_bo[n_chr_h]; // contact probability basic observable
  float *cp[n_chr_h]; // contact probability arrays
  float *hcp[n_chr_h]; // host contact probability arrays

  uint lsdcp[n_chr_h]; // length of sd arrays and cp arrays

  simobs_b *cm_bo; // contact map basic observable
  float *cm; // contact map array
  float *hcm; // host contact map array

  const uint cms; // contact map side
  const uint lcm; // length of cm array

  // Function

  // calculate observables
  void calc_observables();
};

} // namespace mmc

#endif // MMC_CHRANA_H

#ifndef MMC_SOSTAT_H
#define MMC_SOSTAT_H

// Includes

#include "util.cuh" // general utilities

#include <vector> // vector container class

// Namespace

namespace mmc // Marco Mendívil Carboni
{

// Structures

struct gdstat // generic data statistics
{
  double avg; // average
  double var; // variance
  double sem; // standard error of the mean
};

struct idstat : gdstat
{
}; // independent data statistics

struct cdstat : gdstat // correlated data statistics
{
  uint f_n_b; // final number of blocks
};

struct tdstat : cdstat // time series data statistics
{
  uint i_t; // thermalization index
  bool thm; // thermalized
};

struct simobs // simulation observable
{
  std::vector<float> is_ts; // individual simulation time series
  std::vector<tdstat> is_sv; // individual simulation statistics vector

  idstat cs_fs; // combined simulations final statistics

  // calculate last individual simulation statistics
  void calc_last_is_stat();

  // save last individual simulation statistics
  void save_last_is_stat(std::ofstream &txt_out_f); // text output file

  // save last individual simulation statistics summary
  void save_last_is_stat_s(std::ofstream &txt_out_f); // text output file

  // calculate combined simulations final statistics
  void calc_cs_final_stat();

  // save combined simulations final statistics
  void save_cs_final_stat(std::ofstream &txt_out_f); // text output file
};

// Functions

// calculate statistics
void calc_stats(const std::vector<float> &v, // vector
    idstat &s); // statistics

// calculate statistics
void calc_stats(const std::vector<float> &v, // vector
    cdstat &s); // statistics

// calculate statistics
void calc_stats(const std::vector<float> &v, // vector
    tdstat &s); // statistics

// calculate statistics
void calc_stats(const std::vector<tdstat> &v, // vector
    idstat &s); // statistics

} // namespace mmc

#endif // MMC_SOSTAT_H

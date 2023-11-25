#ifndef MMC_SDSTAT_H
#define MMC_SDSTAT_H

//Includes

#include "util.cuh" //general utilities

#include <vector> //vector container class

//Namespace

namespace mmc //Marco Mendívil Carboni
{

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

struct obsdat //observable data
{
  std::vector<float> cs_ts; //current simulation time series
  std::vector<tdstat> is_sv; //individual simulation statistics vector
  idstat cs_fs; //combined simulations final statistics
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

#endif //MMC_SDSTAT_H

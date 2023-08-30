#ifndef MMCC_CUDAUTIL_H
#define MMCC_CUDAUTIL_H

//Includes

#include <stdexcept> //standard exceptions classes

//Namespace

namespace mmcc //Marco Mendívil Carboni code
{

//Structures

struct soa3 //structure of (3) arrays
{
  float *x; //1st component
  float *y; //2nd component
  float *z; //3rd component
};

//Classes

class error : public std::runtime_error //generic exception type
{
  public:
    error(const std::string &msg);
};

//Functions

void cuda_check(cudaError_t result);

void allocate_soa3(soa3 &soa3_ref, int n_elements);

void free_soa3(soa3 &soa3_ref);

FILE *fopen(const char *filename, const char *mode);

} //namespace mmcc

#endif //MMCC_CUDAUTIL_H

//Includes

#include <cstdio> //standard input and output library
#include <ctime> //time utilities library
#include <stdexcept> //standard exceptions classes

#include "../inc/utilities.cuh"

//Namespace

namespace mmcc //Marco Mendívil Carboni code
{

//Functions

error::error(const std::string &msg) : std::runtime_error{msg} {}

void cuda_check(cudaError_t result)
{
  if (result!=cudaSuccess)
  {
    char msg[512];
    std::snprintf(msg,sizeof(msg),"CUDA: %s",cudaGetErrorString(result));
    throw error(msg);
  }
}

void allocate_soa3(soa3 &soa3_ref, int n_elements)
{
  cuda_check(cudaMallocManaged(&(soa3_ref.x),n_elements*sizeof(float)));
  cuda_check(cudaMallocManaged(&(soa3_ref.y),n_elements*sizeof(float)));
  cuda_check(cudaMallocManaged(&(soa3_ref.z),n_elements*sizeof(float)));
}

void free_soa3(soa3 &soa3_ref)
{
  cudaFree(soa3_ref.x);
  cudaFree(soa3_ref.y);
  cudaFree(soa3_ref.z);
}

FILE *fopen(const char *filename, const char *mode)
{
  FILE *f_ptr = std::fopen(filename,mode);
  if (f_ptr==nullptr)
  {
    char msg[512];
    std::snprintf(msg,sizeof(msg),"unable to open %s",filename);
    throw error(msg);
  }
  return f_ptr;
}

} //namespace mmcc

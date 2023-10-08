//Includes

#include "sugrid.cuh" //chromatin simulation

#include "util.cuh" //general utilities
#include "vecops.cuh" //vector operations

#include <cub/device/device_radix_sort.cuh> //cub parallel radix sort

//Namespace

namespace mmc //Marco Mendívil Carboni
{

//Global Functions

//calculate cell and particle indexes
__global__ void calc_indexes(
  const uint n_o, //number of objects
  const float csl, //cell side length
  const uint cps, //cells per side
  uint *uci, //unsorted cell index array
  uint *upi, //unsorted particle index array
  float4 *r) //position array
{
  //calculate particle index
  int i_p = blockIdx.x*blockDim.x+threadIdx.x; //particle index
  if (i_p>=n_o){ return;}
  upi[i_p] = i_p;

  //calculate auxiliary variables
  float3 r_i = make_float3(r[i_p]); //particle position
  int3 ir = floorf(r_i/csl); //integer coordinates
  int iofst = (cps/2)*(1+cps+cps*cps); //index offset

  //calculate cell index
  uci[i_p] = iofst+ir.x+ir.y*cps+ir.z*cps*cps;
}

//set cells empty
__global__ void set_cells_empty(
  const uint n_c, //number of cells
  uint *beg, //cell beginning array
  uint *end) //cell end array
{
  //calculate limit array index
  int lai = blockIdx.x*blockDim.x+threadIdx.x; //limit array index
  if (lai>=n_c){ return;}

  //set beginning and end of cells
  beg[lai] = 0;
  end[lai] = 0;
}

//find beginning and end of each cell
__global__ void find_cells_limits(
  const uint n_o, //number of objects
  uint *sci, //sorted cell index array
  uint *beg, //cell beginning array
  uint *end) //cell end array
{
  //calculate sorted array index
  int sai = blockIdx.x*blockDim.x+threadIdx.x; //sorted array index
  if (sai>=n_o){ return;}

  //set beginning and end of cells
  int ci_curr = sci[sai]; //current cell index
  if (sai==0)
  {
    beg[ci_curr] = sai;
    return;
  }
  int ci_prev = sci[sai-1]; //previous cell index
  if (ci_prev!=ci_curr)
  {
    beg[ci_curr] = sai;
    end[ci_prev] = sai;
  }
  if (sai==n_o-1)
  {
    end[ci_curr] = sai+1;
  }
}

//Host Functions

//sorted uniform grid constructor
sugrid::sugrid(
    const uint n_o, //number of objects
    const float csl, //cell side length
    const uint cps) //cells per side
  : n_o {n_o}
  , csl {csl}
  , cps {cps}
  , n_c {cps*cps*cps}
{
  //allocate arrays
  cuda_check(cudaMalloc(&uci,n_o*sizeof(uint)));
  cuda_check(cudaMalloc(&sci,n_o*sizeof(uint)));
  cuda_check(cudaMalloc(&upi,n_o*sizeof(uint)));
  cuda_check(cudaMalloc(&spi,n_o*sizeof(uint)));
  cuda_check(cudaMalloc(&beg,n_c*sizeof(uint)));
  cuda_check(cudaMalloc(&end,n_c*sizeof(uint)));

  //allocate extra buffer
  cub::DeviceRadixSort::SortPairs(nullptr,ebs,uci,sci,upi,spi,n_o);
  cuda_check(cudaMalloc(&eb,ebs));
}

//sorted uniform grid delegating constructor
sugrid::sugrid(
    const uint n_o, //number of objects
    const sugrid &g) //grid
  : sugrid(n_o,g.csl,g.cps) {}

//sorted uniform grid destructor
sugrid::~sugrid()
{
  //deallocate arrays
  cuda_check(cudaFree(uci));
  cuda_check(cudaFree(sci));
  cuda_check(cudaFree(upi));
  cuda_check(cudaFree(spi));
  cuda_check(cudaFree(beg));
  cuda_check(cudaFree(end));

  //deallocate extra buffer
  cuda_check(cudaFree(eb));
}

//generate grid arrays
void sugrid::generate_arrays(
  int tpb, //threads per block
  float4 *r) //position array
{
  calc_indexes<<<(n_o+tpb-1)/tpb,tpb>>>(n_o,csl,cps,uci,upi,r);
  cub::DeviceRadixSort::SortPairs(eb,ebs,uci,sci,upi,spi,n_o);
  set_cells_empty<<<(n_c+tpb-1)/tpb,tpb>>>(n_c,beg,end);
  find_cells_limits<<<(n_o+tpb-1)/tpb,tpb>>>(n_o,sci,beg,end);
}

} //namespace mmc

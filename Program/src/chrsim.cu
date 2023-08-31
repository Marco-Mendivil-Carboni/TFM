//Includes

#include <cstdio> //standard input and output library
#include <cmath> //mathematical functions library
#include <ctime> //time utilities library

#include <curand_kernel.h> //cuRAND device functions
#include </usr/local/cuda/samples/common/inc/helper_math.h> //float4 utilities

#include "../inc/util.hpp"
#include "../inc/chrsim.cuh"

//Namespace

namespace mmcc //Marco Mendívil Carboni code
{

//Constants

static constexpr float xi  = 1.000000; //damping coefficient
static constexpr float k_B = 0.001120; //Boltzmann constant
static constexpr float l_0 = 1.000000; //bond natural length
static constexpr float k_e = 100.0000; //elastic constant
static constexpr float k_b = 2.000000; //bending constant
static constexpr float r_c = 1.122462; //LJ cutoff radius
static constexpr float dt  = 1.0/2048; //timestep

static constexpr int n_s = 1*2048; //MD steps between frames

//Host Functions

//check for errors in cuda runtime API call
void cuda_check(cudaError_t result)
{
  if (result!=cudaSuccess)
  {
    char msg[512];
    std::snprintf(msg,sizeof(msg),"CUDA: %s",cudaGetErrorString(result));
    throw error(msg);
  }
}

//chrsim constructor
chrsim::chrsim(FILE *f_ptr_par)
{
  read_parameters(f_ptr_par);

  n_blocks = (ap.N+threads_block-1)/threads_block;
  n_threads = n_blocks*threads_block;

  cuda_check(cudaMallocManaged(&r_2,ap.N*sizeof(float4)));
  cuda_check(cudaMallocManaged(&r_1,ap.N*sizeof(float4)));

  cuda_check(cudaMallocManaged(&f_2,ap.N*sizeof(float4)));
  cuda_check(cudaMallocManaged(&f_1,ap.N*sizeof(float4)));

  cuda_check(cudaMallocManaged(&nrn,n_threads*sizeof(float4)));

  cuda_check(cudaMallocManaged(&state,n_threads*sizeof(PRNGstate)));
}

//chrsim destructor
chrsim::~chrsim()
{
  cudaFree(r_2);
  cudaFree(r_1);

  cudaFree(f_2);
  cudaFree(f_1);

  cudaFree(nrn);

  cudaFree(state);
}

//generate a random initial configuration of chromatin
void chrsim::generate_initial_configuration()
{
  //initialize PRNG
  curandGenerator_t gen;
  curandCreateGeneratorHost(&gen,CURAND_RNG_PSEUDO_DEFAULT);
  curandSetPseudoRandomGeneratorSeed(gen,time(nullptr));

  //declare auxiliary variables
  float beta = 1.0/(k_B*ap.T);
  float rand, theta, phi, bondlen, bondangle;
  float3 randdir, olddir, newdir, perdir;

  //place first particle
  r_2[0] = make_float4(0.0);
  curandGenerateUniform(gen,&rand,1); theta = acos(1.0-2.0*rand);
  curandGenerateUniform(gen,&rand,1); phi = 2.0*M_PI*rand;
  randdir = make_float3(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
  olddir = randdir;

  //set initial sigma value
  sig = 1.0/2;

  //perform random walk
  int n_failures = 0;
  for (int i_p = 1; i_p<ap.N; ++i_p)
  {
    //generate random direction perpendicular to old direction
    curandGenerateUniform(gen,&rand,1); theta = acos(1.0-2.0*rand);
    curandGenerateUniform(gen,&rand,1); phi = 2.0*M_PI*rand;
    randdir = make_float3(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
    perdir = cross(olddir,randdir);
    perdir = normalize(perdir);

    //generate random bond angle and calculate new direction
    curandGenerateUniform(gen,&rand,1);
    bondangle = acos(1.0+log(1.0-(1.0-exp(-2.0*(k_b*beta)))*rand)/(k_b*beta));
    newdir = cos(bondangle)*olddir+sin(bondangle)*perdir;

    //calculate position of next particle
    curandGenerateUniform(gen,&rand,1);
    bondlen = l_0+sqrt(2.0/(k_e*beta))*erfinv(2.0*rand-1.0);
    r_2[i_p] = make_float4(bondlen*newdir)+r_2[i_p-1];

    //check if new position is acceptable
    int accept = 1;
    if (!isfinite(r_2[i_p].x)){ accept = 0;}
    if (!isfinite(r_2[i_p].y)){ accept = 0;}
    if (!isfinite(r_2[i_p].z)){ accept = 0;}
    for (int j_p = 0; j_p<(i_p-1); ++j_p)
    {
      float dist = length(make_float3(r_2[j_p]-r_2[i_p]));
      if (dist<(r_c*sig)){ accept = 0; break;}
    }
    float d_r = length(make_float3(r_2[i_p]));
    if ((ap.R-d_r)<(r_c*sig)){ accept = 0;}

    //continue if it is accepted
    if (accept)
    {
      olddir = newdir;
      n_failures = 0;
    }
    else
    {
      ++n_failures;
      if( n_failures>1024){ i_p = 0;}
      else{ i_p--;}
    }
  }

  //free PRNG state
  curandDestroyGenerator(gen);
}

//write initial configuration to file in gro format
void chrsim::write_initial_configuration(FILE *f_ptr)
{
  std::fprintf(f_ptr,"Chromatin chrsim, t=0.0\n");
  std::fprintf(f_ptr,"%5d\n",ap.N);
  for( int i_p = 0; i_p<ap.N; ++i_p)
  {
    std::fprintf(f_ptr,"%5d%-5s%5s%5d",i_p+1,"X","X",i_p+1);
    std::fprintf(f_ptr,"%8.3f%8.3f%8.3f\n",r_2[i_p].x,r_2[i_p].y,r_2[i_p].z);
  }
  std::fprintf(f_ptr,"%10.5f%10.5f%10.5f\n",0.0,0.0,0.0);
}

//read adjustable parameters from file
void chrsim::read_parameters(FILE *f_ptr_par)
{
  if (std::fscanf(f_ptr_par,"T\t%f\n",&(ap.T))!=1
    ||std::fscanf(f_ptr_par,"N\t%d\n",&(ap.N))!=1
    ||std::fscanf(f_ptr_par,"R\t%f\n",&(ap.R))!=1
    ||std::fscanf(f_ptr_par,"F\t%d\n",&(ap.F))!=1)
  {
    throw error("unable to read parameters");
  }
  if ((ap.T)<__FLT_MIN__){ throw error("T must be positive");}
  if ((ap.N)<__FLT_MIN__){ throw error("N must be positive");}
  if ((ap.R)<__FLT_MIN__){ throw error("R must be positive");}
  if ((ap.F)<__FLT_MIN__){ throw error("F must be positive");}
}

//Device Functions

// __device__ void example_function(float3 &r)
// {
//   r += ...;
// }

//Global Functions

// __global__ void example_kernel(int N, float4 *r)
// {
//   int i_p = ...;
//   if (i_p<N)
//   {
//     float3 r_int = make_float3(r[i_p]);
//     use r_int by reference in device functions onward
//     example_function(r_int);
//   }
// }

} //namespace mmcc

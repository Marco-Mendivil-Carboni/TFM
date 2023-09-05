//Includes

#include "chrsim.cuh"
#include "util.hpp"

#include <ctime> //time utilities library

#include </usr/local/cuda/samples/common/inc/helper_math.h> //float4 utilities

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

static constexpr int n_s = 1*2048; //RK steps between frames

//Device Functions

//Global Functions

//initialize device PRNG state
__global__ void init_PRNG(prng *state, int seed)
{
  int i_p = blockIdx.x*blockDim.x+threadIdx.x; //particle index
  curand_init(seed,i_p,0,&state[i_p]);
}

//begin Runge-Kutta iteration
__global__ void begin_iter(int N, float c_rn, 
  float4 *f_2, float4 *f_1, float4 *nrn, prng *state)
{
  int i_p = blockIdx.x*blockDim.x+threadIdx.x; //particle index
  float4 nrn_thd = nrn[i_p]; //thread normal random numbers
  nrn_thd.x = c_rn*curand_normal(&state[i_p]);
  nrn_thd.y = c_rn*curand_normal(&state[i_p]);
  nrn_thd.z = c_rn*curand_normal(&state[i_p]);
  nrn[i_p] = nrn_thd;
  if (i_p<N)
  {
    f_2[i_p] = make_float4(0.0);
    f_1[i_p] = make_float4(0.0);
  }
}

//execute 1st stage of the Runge-Kutta method
__global__ void exec_RK_1(int N, 
  float4 *r_2, float4 *r_1, float4 *f_2, float4 *nrn)
{
  int i_p = blockIdx.x*blockDim.x+threadIdx.x; //particle index
  if (i_p<N)
  {

  }
}

//execute 2nd stage of the Runge-Kutta method
__global__ void exec_RK_2(int N, 
  float4 *r_2, float4 *r_1, float4 *f_2, float4 *f_1, float4 *nrn)
{
  int i_p = blockIdx.x*blockDim.x+threadIdx.x; //particle index
  if (i_p<N)
  {

  }
}

//Host Functions

//check for errors in cuda runtime API call
void cuda_check(cudaError_t result)
{
  if (result!=cudaSuccess)
  {
    std::string msg = "cuda: ";
    msg += cudaGetErrorString(result);
    throw error(msg);
  }
}

//chrsim constructor
chrsim::chrsim(std::ifstream &f_par)
{
  //initialize parameters and variables
  read_parameters(f_par);
  n_p_blk = (ap.N+thd_blk-1)/thd_blk;
  n_p_thd = n_p_blk*thd_blk;
  c_rn = sqrt(2.0*xi*k_B*ap.T*dt);

  //allocate unified memory
  cuda_check(cudaMallocManaged(&r_2,ap.N*sizeof(float4)));
  cuda_check(cudaMallocManaged(&r_1,ap.N*sizeof(float4)));
  cuda_check(cudaMallocManaged(&f_2,ap.N*sizeof(float4)));
  cuda_check(cudaMallocManaged(&f_1,ap.N*sizeof(float4)));
  cuda_check(cudaMallocManaged(&nrn,n_p_thd*sizeof(float4)));
  cuda_check(cudaMallocManaged(&state,n_p_thd*sizeof(prng)));

  //initialize PRNG
  init_PRNG<<<n_p_blk,thd_blk>>>(state,time(nullptr));
  cuda_check(cudaDeviceSynchronize());
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

//generate a random initial condition
void chrsim::generate_initial_condition()
{
  //initialize PRNG
  curandGenerator_t gen; //host PRNG
  curandCreateGeneratorHost(&gen,CURAND_RNG_PSEUDO_DEFAULT);
  curandSetPseudoRandomGeneratorSeed(gen,time(nullptr));

  //declare auxiliary variables
  float beta = 1.0/(k_B*ap.T); //inverse temperature
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
    bool p_a = true; //position is acceptable
    if (!isfinite(r_2[i_p].x)){ p_a = false;}
    if (!isfinite(r_2[i_p].y)){ p_a = false;}
    if (!isfinite(r_2[i_p].z)){ p_a = false;}
    for (int j_p = 0; j_p<(i_p-1); ++j_p)
    {
      float dist = length(make_float3(r_2[j_p]-r_2[i_p]));
      if (dist<(r_c*sig)){ p_a = false; break;}
    }
    float d_r = length(make_float3(r_2[i_p]));
    if ((ap.R-d_r)<(r_c*sig)){ p_a = false;}

    if (p_a) //continue
    {
      olddir = newdir;
      n_failures = 0;
    }
    else //go back
    {
      ++n_failures;
      if( n_failures>1024){ i_p = 1;}
      else{ i_p--;}
    }
  }

  //free host PRNG
  curandDestroyGenerator(gen);

  //record success message
  logger::record("initial condition generated");
}

//write initial condition to file in gro format
void chrsim::write_initial_condition(std::ofstream &f_i_c)
{
  f_i_c<<"Chromatin simulation, i_f = 0 t = 0.0\n";
  f_i_c<<cnfs(ap.N,5,' ')<<"\n";
  for (int i_p = 0; i_p<ap.N; ++i_p)
  {
    f_i_c<<std::setw(5)<<i_p+1<<std::left<<std::setw(5)<<"X"<<std::right;
    f_i_c<<std::setw(5)<<"X"<<std::setw(5)<<i_p+1;
    f_i_c<<cnfs(r_2[i_p].x,8,' ',3);
    f_i_c<<cnfs(r_2[i_p].y,8,' ',3);
    f_i_c<<cnfs(r_2[i_p].z,8,' ',3);
    f_i_c<<"\n";
  }
  f_i_c<<cnfs(0.0,10,' ',5);
  f_i_c<<cnfs(0.0,10,' ',5);
  f_i_c<<cnfs(0.0,10,' ',5);
  f_i_c<<"\n";
}

//save simulation state to binary file
void chrsim::save_checkpoint(std::ofstream &f_chkp)
{
  f_chkp.write(reinterpret_cast<char *>(&i_f),sizeof(i_f));
  f_chkp.write(reinterpret_cast<char *>(&t),sizeof(t));
  f_chkp.write(reinterpret_cast<char *>(r_2),ap.N*sizeof(float4));
  f_chkp.write(reinterpret_cast<char *>(state),n_p_thd*sizeof(prng));
  logger::record("simulation checkpoint saved");
}

//load simulation state from binary file
void chrsim::load_checkpoint(std::ifstream &f_chkp)
{
  f_chkp.read(reinterpret_cast<char *>(&i_f),sizeof(i_f));
  f_chkp.read(reinterpret_cast<char *>(&t),sizeof(t));
  f_chkp.read(reinterpret_cast<char *>(r_2),ap.N*sizeof(float4));
  f_chkp.read(reinterpret_cast<char *>(state),n_p_thd*sizeof(prng));
  logger::record("simulation checkpoint loaded");
}

//write trajectory to binary file in trr format
void chrsim::write_trajectory(std::ofstream &f_traj)
{
  ++i_f;
  int32_t header[] = {1993, 13, 12, 
    1599622471, 1601073780, 1701603686, 
    0, 0, 0, 0, 0, 0, 0, 3*ap.N*4, 0, 0, ap.N, i_f, 0, 
    *(reinterpret_cast<int32_t *>(&t)), 0}; //trr file header
  //for more information on the contents of the header see chemfiles
  f_traj.write(reinterpret_cast<char *>(header),sizeof(header));
  for (int i_p = 0; i_p<ap.N; ++i_p)
  {
    f_traj.write(reinterpret_cast<char *>(&(r_2[i_p].x)),4);
    f_traj.write(reinterpret_cast<char *>(&(r_2[i_p].y)),4);
    f_traj.write(reinterpret_cast<char *>(&(r_2[i_p].z)),4);
  }
  //this is a minimal trr file writing routine that doesn't rely on \ 
  //the xdr library but only works with vmd in little endian systems
}

//read adjustable parameters from file
void chrsim::read_parameters(std::ifstream &f_par)
{
  std::string key;
  f_par>>key>>(ap.T); if (key!="T"||ap.T<0){ throw error("error reading T");}
  f_par>>key>>(ap.N); if (key!="N"||ap.N<1){ throw error("error reading N");}
  f_par>>key>>(ap.R); if (key!="R"||ap.R<0){ throw error("error reading R");}
  f_par>>key>>(ap.F); if (key!="F"||ap.F<1){ throw error("error reading F");}
  std::string msg = "parameters:";
  msg += " T = "+cnfs(ap.T,6,'0',2);
  msg += " N = "+cnfs(ap.N,5,'0');
  msg += " R = "+cnfs(ap.R,6,'0',2);
  msg += " F = "+cnfs(ap.F,5,'0');
  logger::record(msg);
  float cvf = ap.N*pow(0.5/(ap.R-0.5),3); //chromatin volume fraction
  if (cvf>0.5){ throw error("chromatin volume fraction above 0.5");}
}

// __device__ void example_function(float3 &r)
// {
//   r += ...;
// }

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

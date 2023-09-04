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

static constexpr int n_s = 1*2048; //MD steps between frames

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
  cuda_check(cudaMallocManaged(&state,n_p_thd*sizeof(PRNGstate)));

  //initialize PRNG
  // setup_PRNG<<<n_p_blk,thd_blk>>>(time(nullptr),state);
  // cuda_check(cudaDeviceSynchronize());
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
  f_i_c<<"Chromatin simulation, t=0.0\n";
  f_i_c<<cnfs(ap.N,5,false)<<"\n";
  for (int i_p = 0; i_p<ap.N; ++i_p)
  {
    f_i_c<<std::setw(5)<<i_p+1<<std::left<<std::setw(5)<<"X"<<std::right;
    f_i_c<<std::setw(5)<<"X"<<std::setw(5)<<i_p+1;
    f_i_c<<cnfs(r_2[i_p].x,8,false,3);
    f_i_c<<cnfs(r_2[i_p].y,8,false,3);
    f_i_c<<cnfs(r_2[i_p].z,8,false,3);
    f_i_c<<"\n";
  }
  f_i_c<<cnfs(0.0,10,false,5);
  f_i_c<<cnfs(0.0,10,false,5);
  f_i_c<<cnfs(0.0,10,false,5);
  f_i_c<<"\n";
}

//save simulation state to binary file
void chrsim::save_checkpoint(std::ofstream &f_chkp)
{
  f_chkp.write(reinterpret_cast<char *>(&t),sizeof(t));
  f_chkp.write(reinterpret_cast<char *>(r_2),ap.N*sizeof(float4));
  f_chkp.write(reinterpret_cast<char *>(state),n_p_thd*sizeof(PRNGstate));
}

//load simulation state from binary file
void chrsim::load_checkpoint(std::ifstream &f_chkp)
{
  f_chkp.read(reinterpret_cast<char *>(&t),sizeof(t));
  f_chkp.read(reinterpret_cast<char *>(r_2),ap.N*sizeof(float4));
  f_chkp.read(reinterpret_cast<char *>(state),n_p_thd*sizeof(PRNGstate));
}

//write trajectory to binary file in trr format
void chrsim::write_trajectory(std::ofstream &f_traj, int i_f)
{
  int array1[] = {1993, 13, 12};
  f_traj.write(reinterpret_cast<char *>(array1),sizeof(array1));
  char trr_version[] = "GMX_trn_file";
  f_traj.write((trr_version),sizeof(trr_version)-1);
  int array2[] = {0, 0, 0, 0, 0, 0, 0};
  f_traj.write(reinterpret_cast<char *>(array2),sizeof(array2));
  int r_size = (3*ap.N*sizeof(float));
  f_traj.write(reinterpret_cast<char *>(&r_size),sizeof(r_size));
  int array3[] = {0, 0};
  f_traj.write(reinterpret_cast<char *>(array3),sizeof(array3));
  int natoms = (ap.N);
  f_traj.write(reinterpret_cast<char *>(&natoms),sizeof(natoms));
  int frameidx = (i_f);
  f_traj.write(reinterpret_cast<char *>(&frameidx),sizeof(frameidx));
  int zero = 0;
  f_traj.write(reinterpret_cast<char *>(&zero),sizeof(zero));
  f_traj.write(reinterpret_cast<char *>(&t),sizeof(t));
  f_traj.write(reinterpret_cast<char *>(&zero),sizeof(zero));
  for (int i_p = 0; i_p<ap.N; ++i_p)
  {
    f_traj.write(reinterpret_cast<char *>(&(r_2[i_p].x)),sizeof(float));
    f_traj.write(reinterpret_cast<char *>(&(r_2[i_p].y)),sizeof(float));
    f_traj.write(reinterpret_cast<char *>(&(r_2[i_p].z)),sizeof(float));
  } 
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
  msg += " T = "+cnfs(ap.T,6,false,2);
  msg += " N = "+cnfs(ap.N,5);
  msg += " R = "+cnfs(ap.R,6,false,2);
  msg += " F = "+cnfs(ap.F,5);
  logger::record(msg);
  float cvf = ap.N*pow(0.5/(ap.R-0.5),3); //chromatin volume fraction
  if (cvf>0.5){ throw error("chromatin volume fraction above 0.5");}
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

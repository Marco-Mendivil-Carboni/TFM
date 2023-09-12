//Includes

#include "chrsim.cuh" //chromatin simulation

#include <time.h> //time utilities library
#include </usr/local/cuda/samples/common/inc/helper_math.h> //float4 utilities

//Namespace

namespace mmc //Marco Mendívil Carboni
{

//Constants

static constexpr float k_B = 0.001120; //Boltzmann constant
static constexpr float dt  = 1.0/2048; //timestep

//Device Functions

//calculate bonded forces
inline __device__ void calc_bonded_f(
  const int N, //number of particles
  int i_p, //particle index
  float4 *r, //position array
  float4 *f) //force array
{
  //declare auxiliary variables
  float3 b_vec[4]; //bond vectors
  float b_il[4]; //bond inverse lengths
  float b_cos[3]; //bond angle cosines
  float3 f_b = make_float3(0.0); //bonded forces

  //calculate bond vectors, inverse lengths and angle cosines
  for (int i_b = 0; i_b<4; ++i_b) //bond index
  {
    if ((i_p+i_b)>=2 && (i_p+i_b)<=N) //calculate variables if bond exists
    {
      b_vec[i_b] = make_float3(r[i_p+i_b-1]-r[i_p+i_b-2]);
      b_il[i_b] = rsqrtf(dot(b_vec[i_b],b_vec[i_b]));
    }
    else //set variables to zero if bond doesn't exist
    {
      b_vec[i_b] = make_float3(0.0);
      b_il[i_b] = 0.0;
    }
  }
  for (int i_c = 0; i_c<3; ++i_c) //cosine index
  {
    b_cos[i_c] = dot(b_vec[i_c+1],b_vec[i_c])*b_il[i_c+1]*b_il[i_c];
  }

  //calculate elastic potential force
  f_b += k_e*(+(1.0-l_0*b_il[2])*b_vec[2]-(1.0-l_0*b_il[1])*b_vec[1]);

  //calculate bending potential force
  f_b += k_b*(+b_il[1]*b_il[0]*b_vec[0]-b_cos[0]*b_vec[1]*b_il[1]*b_il[1]);
  f_b += k_b*(+b_il[1]*b_il[2]*b_vec[2]-b_cos[1]*b_vec[1]*b_il[1]*b_il[1]);
  f_b += k_b*(-b_il[2]*b_il[1]*b_vec[1]+b_cos[1]*b_vec[2]*b_il[2]*b_il[2]);
  f_b += k_b*(-b_il[2]*b_il[3]*b_vec[3]+b_cos[2]*b_vec[2]*b_il[2]*b_il[2]);

  //add result to forces
  f[i_p] += make_float4(f_b);
}

//Global Functions

//initialize device PRNG state array
__global__ void init_PRNG(
  const int N, //number of particles
  prng *dps, //device PRNG state array
  int seed) //PRNG seed
{
  int i_p = blockIdx.x*blockDim.x+threadIdx.x; //particle index
  if (i_p>=N){ return;}
  curand_init(seed,i_p,0,&dps[i_p]);
}

//begin Runge-Kutta iteration
__global__ void begin_iter(
  const int N, //number of particles
  float4 *df2, //device force array 2
  float4 *df1, //device force array 1
  float sd, //random number standard deviation
  float4 *drn, //device random number array
  prng *dps) //device PRNG state array
{
  int i_p = blockIdx.x*blockDim.x+threadIdx.x; //particle index
  if (i_p>=N){ return;}
  drn[i_p] = sd*curand_normal4(&dps[i_p]);
  df2[i_p] = make_float4(0.0);
  df1[i_p] = make_float4(0.0);
}

//execute 1st stage of the Runge-Kutta method
__global__ void exec_RK_1(
  const int N, //number of particles
  float4 *dr2, //device position array 2
  float4 *dr1, //device position array 1
  float4 *df2, //device force array 2
  float4 *drn) //device random number array
{
  int i_p = blockIdx.x*blockDim.x+threadIdx.x; //particle index
  if (i_p>=N){ return;}
  calc_bonded_f(N,i_p,dr2,df2);
  dr1[i_p] = dr2[i_p]+df2[i_p]*dt/xi+drn[i_p]/xi;
}

//execute 2nd stage of the Runge-Kutta method
__global__ void exec_RK_2(
  const int N, //number of particles
  float4 *dr2, //device position array 2
  float4 *dr1, //device position array 1
  float4 *df2, //device force array 2
  float4 *df1, //device force array 1
  float4 *drn) //device random number array
{
  int i_p = blockIdx.x*blockDim.x+threadIdx.x; //particle index
  if (i_p>=N){ return;}
  calc_bonded_f(N,i_p,dr1,df1);
  dr2[i_p] = dr2[i_p]+0.5*(df1[i_p]+df2[i_p])*dt/xi+drn[i_p]/xi;
}

//Host Functions

//chrsim constructor
chrsim::chrsim(parmap &par) //parameters
  : chrdat(par)
  , f_f {par.get_val<int>("frames_per_file",100)}
  , s_f {par.get_val<int>("steps_per_frame",1*2048)}
  , thd_blk {par.get_val<int>("threads_per_block",256)}
  , n_p_blk {(N+thd_blk-1)/thd_blk}
  , sd {sqrtf(2.0*xi*k_B*T*dt)}
{
  //check parameters
  if (f_f<1){ throw error("frames_per_file out of range");}
  if (s_f<1){ throw error("steps_per_frame out of range");}
  if (thd_blk<1){ throw error("threads_per_block out of range");}
  std::string msg = "chrsim parameters:"; //message
  msg += " f_f = "+cnfs(f_f,5,'0');
  msg += " s_f = "+cnfs(s_f,5,'0');
  msg += " thd_blk = "+cnfs(thd_blk,5,'0');
  logger::record(msg);

  //allocate unified memory
  cuda_check(cudaMallocManaged(&dpt,N*sizeof(int)));
  cuda_check(cudaMallocManaged(&dr2,N*sizeof(float4)));
  cuda_check(cudaMallocManaged(&dr1,N*sizeof(float4)));
  cuda_check(cudaMallocManaged(&df2,N*sizeof(float4)));
  cuda_check(cudaMallocManaged(&df1,N*sizeof(float4)));
  cuda_check(cudaMallocManaged(&drn,N*sizeof(float4)));
  cuda_check(cudaMallocManaged(&dps,N*sizeof(prng)));

  //initialize PRNG
  init_PRNG<<<n_p_blk,thd_blk>>>(N,dps,time(nullptr));
  cuda_check(cudaDeviceSynchronize());
}

//chrsim destructor
chrsim::~chrsim()
{
  cudaFree(dpt);
  cudaFree(dr2);
  cudaFree(dr1);
  cudaFree(df2);
  cudaFree(df1);
  cudaFree(drn);
  cudaFree(dps);
}

//generate a random initial condition
void chrsim::generate_initial_condition()
{
  //initialize PRNG
  curandGenerator_t gen; //host PRNG
  curandCreateGeneratorHost(&gen,CURAND_RNG_PSEUDO_DEFAULT);
  curandSetPseudoRandomGeneratorSeed(gen,time(nullptr));

  //declare auxiliary variables
  float beta = 1.0/(k_B*T); //inverse temperature
  float rand; //random number in (0,1]
  float theta; //polar angle
  float phi; //azimuthal angle
  float len_b; //bond length
  float angle_b; //bond angle
  float3 olddir; //old direction
  float3 newdir; //new direction
  float3 randdir; //random direction
  float3 perpdir; //perpendicular direction

  //place first particle
  dr2[0] = make_float4(0.0);
  curandGenerateUniform(gen,&rand,1); theta = acos(1.0-2.0*rand);
  curandGenerateUniform(gen,&rand,1); phi = 2.0*M_PI*rand;
  randdir = make_float3(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
  olddir = randdir;

  //reduce sigma for random walk
  sig = 1.0/2;

  //perform random walk
  int att = 0; //number of attempts
  for (int i_p = 1; i_p<N; ++i_p) //particle index
  {
    //generate random direction perpendicular to old direction
    curandGenerateUniform(gen,&rand,1); theta = acos(1.0-2.0*rand);
    curandGenerateUniform(gen,&rand,1); phi = 2.0*M_PI*rand;
    randdir = make_float3(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
    perpdir = cross(olddir,randdir);
    perpdir = normalize(perpdir);

    //generate random bond angle and calculate new direction
    curandGenerateUniform(gen,&rand,1);
    angle_b = acos(1.0+log(1.0-(1.0-exp(-2.0*(k_b*beta)))*rand)/(k_b*beta));
    newdir = cos(angle_b)*olddir+sin(angle_b)*perpdir;

    //calculate position of next particle
    curandGenerateUniform(gen,&rand,1);
    len_b = l_0+sqrt(2.0/(k_e*beta))*erfinv(2.0*rand-1.0);
    dr2[i_p] = make_float4(len_b*newdir)+dr2[i_p-1];

    //check if new position is acceptable
    bool p_a = true; //position is acceptable
    if (!isfinite(dr2[i_p].x)){ p_a = false;}
    if (!isfinite(dr2[i_p].y)){ p_a = false;}
    if (!isfinite(dr2[i_p].z)){ p_a = false;}
    for (int j_p = 0; j_p<(i_p-1); ++j_p) //secondary particle index
    {
      float d_pp; //particle-particle distance
      d_pp = length(make_float3(dr2[j_p]-dr2[i_p]));
      if (d_pp<(r_c*sig)){ p_a = false; break;}
    }
    float d_r; //radial distance to origin
    d_r = length(make_float3(dr2[i_p]));
    if ((R-d_r)<(r_c*sig)){ p_a = false;}

    if (p_a) //continue
    {
      olddir = newdir;
      att = 0;
    }
    else //go back
    {
      ++att;
      if (att>1024){ i_p = 1;}
      else{ i_p--;}
    }
  }

  //expand beads and reset sigma
  while (sig<1.0)
  {
    make_RK_iteration();
    sig += dt/(32*sig*sig);
  }
  cuda_check(cudaDeviceSynchronize());
  sig = 1.0;

  //free host PRNG
  curandDestroyGenerator(gen);

  //record success message
  logger::record("initial condition generated");
}

//write initial condition to file in gro format
void chrsim::write_initial_condition(std::ofstream &f_ic) //IC file
{
  f_ic<<"Chromatin simulation, i_f = 0, t = 0.0\n";
  f_ic<<cnfs(N,5,' ')<<"\n";
  for (int i_p = 0; i_p<N; ++i_p) //particle index
  {
    f_ic<<std::setw(5)<<i_p+1<<std::left<<std::setw(5)<<"X"<<std::right;
    f_ic<<std::setw(5)<<"X"<<std::setw(5)<<i_p+1;
    f_ic<<cnfs(dr2[i_p].x,8,' ',3);
    f_ic<<cnfs(dr2[i_p].y,8,' ',3);
    f_ic<<cnfs(dr2[i_p].z,8,' ',3);
    f_ic<<"\n";
  }
  f_ic<<cnfs(0.0,10,' ',5);
  f_ic<<cnfs(0.0,10,' ',5);
  f_ic<<cnfs(0.0,10,' ',5);
  f_ic<<"\n";
}

//save simulation state to binary file
void chrsim::save_checkpoint(std::ofstream &f_chkp) //checkpoint file
{
  f_chkp.write(reinterpret_cast<char *>(&i_f),sizeof(i_f));
  f_chkp.write(reinterpret_cast<char *>(&t),sizeof(t));
  f_chkp.write(reinterpret_cast<char *>(dr2),N*sizeof(float4));
  f_chkp.write(reinterpret_cast<char *>(dps),N*sizeof(prng));
  logger::record("simulation checkpoint saved");
}

//load simulation state from binary file
void chrsim::load_checkpoint(std::ifstream &f_chkp) //checkpoint file
{
  f_chkp.read(reinterpret_cast<char *>(&i_f),sizeof(i_f));
  f_chkp.read(reinterpret_cast<char *>(&t),sizeof(t));
  f_chkp.read(reinterpret_cast<char *>(dr2),N*sizeof(float4));
  f_chkp.read(reinterpret_cast<char *>(dps),N*sizeof(prng));
  logger::record("simulation checkpoint loaded");
}

//run simulation and write trajectory file
void chrsim::run_simulation(std::ofstream &f_traj) //trajectory file
{
  for (int i_ff = 0; i_ff<f_f; ++i_ff) //file frame index
  {
    float prog_pc = (100.0*i_ff)/(f_f); //progress percentage
    logger::show_prog_pc(prog_pc);
    for (int i_fs = 0; i_fs<s_f; ++i_fs) //frame step index
    {
      make_RK_iteration();
    }
    cuda_check(cudaDeviceSynchronize());
    ++i_f; t += s_f*dt;
    write_trajectory_frame(f_traj);
  }
}

//make one iteration of the Runge-Kutta method
void chrsim::make_RK_iteration()
{
  begin_iter<<<n_p_blk,thd_blk>>>(N,df2,df1,sd,drn,dps);
  exec_RK_1<<<n_p_blk,thd_blk>>>(N,dr2,dr1,df2,drn);
  exec_RK_2<<<n_p_blk,thd_blk>>>(N,dr2,dr1,df2,df1,drn);
}

//write trajectory frame to binary file in trr format
void chrsim::write_trajectory_frame(std::ofstream &f_traj) //trajectory file
{
  //this is a minimal trr file writing routine that doesn't rely on \ 
  //the xdr library but only works with vmd in little endian systems

  int32_t header[18] = {1993, 1, 0, 
    0, 0, 0, 0, 0, 0, 0, 3*N*4, 0, 0, N, i_f, 0, 
    *(reinterpret_cast<int32_t *>(&t)), 0}; //frame header
  //for more information on the contents of the header see chemfiles
  f_traj.write(reinterpret_cast<char *>(header),sizeof(header));
  for (int i_p = 0; i_p<N; ++i_p) //particle index
  {
    f_traj.write(reinterpret_cast<char *>(&(dr2[i_p].x)),4);
    f_traj.write(reinterpret_cast<char *>(&(dr2[i_p].y)),4);
    f_traj.write(reinterpret_cast<char *>(&(dr2[i_p].z)),4);
  }
}

//check for errors in cuda runtime API call
void cuda_check(cudaError_t rtn_val) //cuda runtime API call return value
{
  if (rtn_val!=cudaSuccess)
  {
    std::string msg = "cuda: "; //error message
    msg += cudaGetErrorString(rtn_val);
    throw error(msg);
  }
}

} //namespace mmc

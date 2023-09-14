//Includes

#include "chrsim.cuh" //chromatin simulation

#include <time.h> //time utilities library

//Namespace

namespace mmc //Marco Mendívil Carboni
{

//Constants

static constexpr float dt = 1.0/2048; //timestep

//Device Functions

//calculate bonded forces
inline __device__ void calc_bonded_f(
  const int N, //number of particles
  int i_p, //particle index
  float4 *r, //position array
  float4 *f) //force array
{
  //declare auxiliary variables
  float3 vec[4]; //bond vectors
  float il[4]; //bond inverse lengths
  float cos[3]; //bond angle cosines
  float3 f_b = {0.0,0.0,0.0}; //bonded forces

  //calculate bond vectors, inverse lengths and angle cosines
  for (int i_b = 0; i_b<4; ++i_b) //bond index
  {
    if ((i_p+i_b)>=2 && (i_p+i_b)<=N) //calculate variables if bond exists
    {
      vec[i_b] = make_float3(r[i_p+i_b-1]-r[i_p+i_b-2]);
      il[i_b] = rsqrtf(dot(vec[i_b],vec[i_b]));
    }
    else //set variables to zero if bond doesn't exist
    {
      vec[i_b] = {0.0,0.0,0.0};
      il[i_b] = 0.0;
    }
  }
  for (int i_c = 0; i_c<3; ++i_c) //cosine index
  {
    cos[i_c] = dot(vec[i_c+1],vec[i_c])*il[i_c+1]*il[i_c];
  }

  //calculate elastic potential force
  f_b += k_e*(+(1.0-l_0*il[2])*vec[2]-(1.0-l_0*il[1])*vec[1]);

  //calculate bending potential force
  f_b += k_b*(+il[1]*il[0]*vec[0]-cos[0]*vec[1]*il[1]*il[1]);
  f_b += k_b*(+il[1]*il[2]*vec[2]-cos[1]*vec[1]*il[1]*il[1]);
  f_b += k_b*(-il[2]*il[1]*vec[1]+cos[1]*vec[2]*il[2]*il[2]);
  f_b += k_b*(-il[2]*il[3]*vec[3]+cos[2]*vec[2]*il[2]*il[2]);

  //add result to force array
  f[i_p] += make_float4(f_b);
}

//Global Functions

//initialize device PRNG state array
__global__ void init_PRNG(
  const int N, //number of particles
  prng *dps, //device PRNG state array
  int ps) //PRNG seed
{
  int i_p = blockIdx.x*blockDim.x+threadIdx.x; //particle index
  if (i_p>=N){ return;}
  curand_init(ps,i_p,0,&dps[i_p]);
}

//execute 1st stage of the Runge-Kutta method
__global__ void exec_RK_1(
  const int N, //number of particles
  float4 *dr2, //device position array 2
  float4 *dr1, //device position array 1
  float4 *df2, //device force array 2
  float sd, //random number standard deviation
  float4 *drn, //device random number array
  prng *dps) //device PRNG state array
{
  int i_p = blockIdx.x*blockDim.x+threadIdx.x; //particle index
  if (i_p>=N){ return;}
  drn[i_p] = sd*curand_normal4(&dps[i_p]);
  df2[i_p] = {0.0,0.0,0.0,0.0};
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
  df1[i_p] = {0.0,0.0,0.0,0.0};
  calc_bonded_f(N,i_p,dr1,df1);
  dr2[i_p] = dr2[i_p]+0.5*(df1[i_p]+df2[i_p])*dt/xi+drn[i_p]/xi;
}

//Host Functions

//chrsim constructor
chrsim::chrsim(parmap &par) //parameters
  : chrdat(par)
  , framepf {par.get_val<int>("frames_per_file",100)}
  , spframe {par.get_val<int>("steps_per_frame",1*2048)}
  , thd_blk {par.get_val<int>("threads_per_block",256)}
  , n_p_blk {(N+thd_blk-1)/thd_blk}
  , sd {sqrtf(2.0*xi*k_B*T*dt)}
{
  //check parameters
  if (framepf<1){ throw error("frames_per_file out of range");}
  if (spframe<1){ throw error("steps_per_frame out of range");}
  if (thd_blk<1){ throw error("threads_per_block out of range");}
  std::string msg = "chrsim:"; //message
  msg += " framepf = "+cnfs(framepf,5,'0');
  msg += " spframe = "+cnfs(spframe,5,'0');
  msg += " thd_blk = "+cnfs(thd_blk,5,'0');
  logger::record(msg);

  //allocate unified memory
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
  cuda_check(cudaFree(dr2));
  cuda_check(cudaFree(dr1));
  cuda_check(cudaFree(df2));
  cuda_check(cudaFree(df1));
  cuda_check(cudaFree(drn));
  cuda_check(cudaFree(dps));
}

//generate a random initial condition
void chrsim::generate_initial_condition()
{
  //initialize host PRNG
  curandGenerator_t gen; //host PRNG
  curandCreateGeneratorHost(&gen,CURAND_RNG_PSEUDO_DEFAULT);
  curandSetPseudoRandomGeneratorSeed(gen,time(nullptr));

  //declare auxiliary variables
  float iT = 1.0/(k_B*T); //inverse temperature
  float ran; //random number in (0,1]
  float theta; //polar angle
  float phi; //azimuthal angle
  float len_b; //bond length
  float angle_b; //bond angle
  float3 old_dir; //old direction
  float3 new_dir; //new direction
  float3 ran_dir; //random direction
  float3 per_dir; //perpendicular direction

  //place first particle
  dr2[0] = {0.0,0.0,0.0,0.0};
  curandGenerateUniform(gen,&ran,1); theta = acos(1.0-2.0*ran);
  curandGenerateUniform(gen,&ran,1); phi = 2.0*M_PI*ran;
  ran_dir = {sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)};
  old_dir = ran_dir;

  //reduce sigma for random walk
  sig = 1.0/2;

  //perform random walk
  int att = 0; //number of attempts
  for (int i_p = 1; i_p<N; ++i_p) //particle index
  {
    //generate random direction perpendicular to old direction
    curandGenerateUniform(gen,&ran,1); theta = acos(1.0-2.0*ran);
    curandGenerateUniform(gen,&ran,1); phi = 2.0*M_PI*ran;
    ran_dir = {sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)};
    per_dir = cross(old_dir,ran_dir);
    per_dir = normalize(per_dir);

    //generate random bond angle and calculate new direction
    curandGenerateUniform(gen,&ran,1);
    angle_b = acos(1.0+log(1.0-(1.0-exp(-2.0*(k_b*iT)))*ran)/(k_b*iT));
    new_dir = cos(angle_b)*old_dir+sin(angle_b)*per_dir;

    //calculate position of next particle
    curandGenerateUniform(gen,&ran,1);
    len_b = l_0+sqrt(2.0/(k_e*iT))*erfinv(2.0*ran-1.0);
    dr2[i_p] = make_float4(len_b*new_dir)+dr2[i_p-1];

    //check if new position is acceptable
    bool p_a = true; //position is acceptable
    if (!isfinite(dr2[i_p].x)){ p_a = false;}
    if (!isfinite(dr2[i_p].y)){ p_a = false;}
    if (!isfinite(dr2[i_p].z)){ p_a = false;}
    for (int j_p = 0; j_p<(i_p-1); ++j_p) //secondary particle index
    {
      float d_p; //particle-particle distance
      d_p = length(make_float3(dr2[j_p]-dr2[i_p]));
      if (d_p<(r_c*sig)){ p_a = false; break;}
    }
    float d_r; //radial distance to origin
    d_r = length(make_float3(dr2[i_p]));
    if ((R-d_r)<(r_c*sig)){ p_a = false;}

    if (p_a) //continue
    {
      old_dir = new_dir;
      att = 0;
    }
    else //go back
    {
      ++att;
      if (att>1024){ i_p = 1;}
      else{ i_p--;}
    }
  }

  //expand beads
  while (sig<1.0)
  {
    make_RK_iteration();
    sig += dt/(32*sig*sig);
  }
  cuda_check(cudaMemcpy(r,dr2,N*sizeof(float4),cudaMemcpyDefault));

  //reset sigma
  sig = 1.0;

  //free host PRNG
  curandDestroyGenerator(gen);

  //record success message
  logger::record("initial condition generated");
}

//save simulation state to binary file
void chrsim::save_checkpoint(std::ofstream &bin_out_f) //binary output file
{
  bin_out_f.write(reinterpret_cast<char *>(&i_f),sizeof(i_f));
  bin_out_f.write(reinterpret_cast<char *>(&t),sizeof(t));
  bin_out_f.write(reinterpret_cast<char *>(dr2),N*sizeof(float4));
  bin_out_f.write(reinterpret_cast<char *>(dps),N*sizeof(prng));
  logger::record("simulation checkpoint saved");
}

//load simulation state from binary file
void chrsim::load_checkpoint(std::ifstream &bin_inp_f) //binary input file
{
  bin_inp_f.read(reinterpret_cast<char *>(&i_f),sizeof(i_f));
  bin_inp_f.read(reinterpret_cast<char *>(&t),sizeof(t));
  bin_inp_f.read(reinterpret_cast<char *>(dr2),N*sizeof(float4));
  bin_inp_f.read(reinterpret_cast<char *>(dps),N*sizeof(prng));
  logger::record("simulation checkpoint loaded");
}

//run simulation and trajectory to binary file
void chrsim::run_simulation(std::ofstream &bin_out_f) //binary output file
{
  for (int ffi = 0; ffi<framepf; ++ffi) //file frame index
  {
    float prog_pc = (100.0*ffi)/(framepf); //progress percentage
    logger::show_prog_pc(prog_pc);
    for (int fsi = 0; fsi<spframe; ++fsi) //frame step index
    {
      make_RK_iteration();
    }
    cuda_check(cudaMemcpy(r,dr2,N*sizeof(float4),cudaMemcpyDefault));
    ++i_f; t += spframe*dt;
    write_frame_bin(bin_out_f);
  }
}

//make one iteration of the Runge-Kutta method
void chrsim::make_RK_iteration()
{
  exec_RK_1<<<n_p_blk,thd_blk>>>(N,dr2,dr1,df2,sd,drn,dps);
  exec_RK_2<<<n_p_blk,thd_blk>>>(N,dr2,dr1,df2,df1,drn);
}

} //namespace mmc

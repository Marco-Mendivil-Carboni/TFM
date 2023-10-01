//Includes

#include "chrsim.cuh" //chromatin simulation

#include <time.h> //time utilities library

#include <curand_kernel.h> //cuRAND device functions

//Namespace

namespace mmc //Marco Mendívil Carboni
{

//Constants

static constexpr float dt = 1.0/2048; //timestep
static constexpr float rco = 1.154701; //Wang-Frenkel repulsive cutoff
static constexpr float aco = 2.000000; //Wang-Frenkel attractive cutoff

//Aliases

using prng = curandStatePhilox4_32_10; //PRNG type

//Enumerations

enum srint //short-range interaction
{
  WFI, //Wang-Frenkel interaction
  SRI //Soft-Repulsive interaction
};

//Device Functions

//calculate bonded forces
inline __device__ void calc_bf(
  const uint N, //number of particles
  uint i_p, //particle index
  float4 *r, //position array
  float4 *f) //force array
{
  //declare auxiliary variables
  float3 vec[4]; //bond vectors
  float il[4]; //bond inverse lengths
  float cos[3]; //bond angle cosines
  float3 bf = {0.0,0.0,0.0}; //bonded forces

  //calculate bond vectors, inverse lengths and angle cosines
  for (uint i_b = 0; i_b<4; ++i_b) //bond index
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
  for (uint i_c = 0; i_c<3; ++i_c) //cosine index
  {
    cos[i_c] = dot(vec[i_c+1],vec[i_c])*il[i_c+1]*il[i_c];
  }

  //calculate elastic potential force
  bf += k_e*(+(1.0-l_0*il[2])*vec[2]-(1.0-l_0*il[1])*vec[1]);

  //calculate bending potential force
  bf += k_b*(+il[1]*il[0]*vec[0]-cos[0]*vec[1]*il[1]*il[1]);
  bf += k_b*(+il[1]*il[2]*vec[2]-cos[1]*vec[1]*il[1]*il[1]);
  bf += k_b*(-il[2]*il[1]*vec[1]+cos[1]*vec[2]*il[2]*il[2]);
  bf += k_b*(-il[2]*il[3]*vec[3]+cos[2]*vec[2]*il[2]*il[2]);

  //add result to force array
  f[i_p] += make_float4(bf);
}

//calculate confinement force
inline __device__ void calc_cf(
  const float R, //confinement radius
  uint i_p, //particle index
  float4 *r, //position array
  float4 *f) //force array
{
  //calculate auxiliary variables
  float d_r; //radial distance to origin
  d_r = length(make_float3(r[i_p]));
  float dwp = (R+0.5)-d_r; //wall-particle distance
  if (dwp>rco){ return;}

  //calculate confinement force
  float3 cf = make_float3(-r[i_p]/d_r); //confinement force
  float d2 = dwp*dwp; //dwp squared
  cf *= (18*d2*d2-96*d2+96)/(d2*d2*d2*d2);

  //add result to force array
  f[i_p] += make_float4(cf);
}

//calculate short-range force
template <srint I> inline __device__ void calc_srf(
  const float eps, //particle energy
  float3 vpp, //particle particle vector
  float dpp, //particle particle distance
  float3 &srf); //short-range forces

//calculate Wang-Frenkel force
template <> inline __device__ void calc_srf<WFI>(
  const float eps, //particle energy
  float3 vpp, //particle particle vector
  float dpp, //particle particle distance
  float3 &srf) //short-range forces
{
  float d2 = dpp*dpp; //dpp squared
  srf += eps*(18*d2*d2-96*d2+96)/(d2*d2*d2*d2)*vpp;
}

//calculate Soft-Repulsive force
template <> inline __device__ void calc_srf<SRI>(
  const float eps, //particle energy
  float3 vpp, //particle particle vector
  float dpp, //particle particle distance
  float3 &srf) //short-range forces
{
  srf += (12-6*dpp)*vpp;
}

//calculate short-range forces with cell's particles
template <srint I> inline __device__ void calc_cell_srf(
  const float eps, //particle energy
  uint i_c, //cell index
  uint i_p, //particle index
  float3 r_i, //particle position
  float4 *r, //position array
  sugrid *srg_p, //short-range grid pointer
  float3 &srf) //short-range forces
{
  //declare auxiliary variables
  uint j_p; //secondary particle index
  float3 r_j; //secondary particle position
  uint beg = srg_p->beg[i_c]; //cell beginning
  uint end = srg_p->end[i_c]; //cell end

  //check cell isn't empty
  if (beg==0xffffffff){ return;}

  //range over cell's particles
  for (uint sai = beg; sai<end; ++sai) //sorted array index
  {
    //get secondary particle index
    j_p = srg_p->spi[sai];

    //calculate force only between non-bonded particles
    if (((j_p>i_p)?j_p-i_p:i_p-j_p)>1)
    {
      //calculate particle particle distance
      r_j = make_float3(r[j_p]);
      float3 vpp = r_i-r_j; //particle particle vector
      float dpp = length(vpp); //particle particle distance
      if (dpp>aco){ continue;}

      //calculate short-range force
      calc_srf<I>(eps,vpp,dpp,srf);
    }
  }
}

//calculate all short-range forces
template <srint I> inline __device__ void calc_all_srf(
  const float eps, //particle energy
  uint i_p, //particle index
  float4 *r, //position array
  float4 *f, //force array
  sugrid *srg_p) //short-range grid pointer
{
  //calculate auxiliary variables
  float3 r_i = make_float3(r[i_p]); //particle position
  const float csl = srg_p->csl; //grid cell side length
  const uint cps = srg_p->cps; //grid cells per side
  const uint n_c = srg_p->n_c; //number of grid cells
  int3 ir = floorf(r_i/csl); //integer coordinates
  uint iofst = (cps/2)*(1+cps+cps*cps); //index offset
  float3 srf = {0.0,0.0,0.0}; //short-range forces

  //range over neighbouring cells
  uint nci; //neighbour cell index
  int3 nir; //neighbour integer coordinates
  for (nir.x = ir.x-1; nir.x<=ir.x+1; ++nir.x)
  {
    for (nir.y = ir.y-1; nir.y<=ir.y+1; ++nir.y)
    {
      for (nir.z = ir.z-1; nir.z<=ir.z+1; ++nir.z)
      {
        //calculate neighbour cell index
        nci = iofst+nir.x+nir.y*cps+nir.z*cps*cps;
        if (nci>=n_c){ continue;}

        //calculate short-range forces with cell's particles
        calc_cell_srf<I>(eps,nci,i_p,r_i,r,srg_p,srf);
      }
    }
  }

  //add result to force array
  f[i_p] += make_float4(srf);
}

//Global Functions

//initialize PRNG state array
__global__ void init_ps(
  const uint N, //number of particles
  void *vps, //void PRNG state array
  uint pseed) //PRNG seed
{
  //calculate particle index
  uint i_p = blockIdx.x*blockDim.x+threadIdx.x; //particle index
  if (i_p>=N){ return;}

  //initialize PRNG state
  prng *ps = static_cast<prng *>(vps); //PRNG state array
  curand_init(pseed,i_p,0,&ps[i_p]);
}

//execute 1st stage of the Runge-Kutta method
template <srint I> __global__ void exec_RK_1(
  const uint N, //number of particles
  const float R, //confinement radius
  const float eps, //particle energy
  float4 *r, //position array
  float4 *f, //force array
  float4 *er, //extra position array
  float sd, //standard deviation
  float4 *rn, //random number array
  void *vps, //void PRNG state array
  sugrid *srg_p) //short-range grid pointer
{
  //calculate particle index
  uint i_p = blockIdx.x*blockDim.x+threadIdx.x; //particle index
  if (i_p>=N){ return;}

  //calculate random numbers
  float3 az; //absolute z-score
  prng *ps = static_cast<prng *>(vps); //PRNG state array
  do
  {
    rn[i_p] = sd*curand_normal4(&ps[i_p]);
    az = fabs(make_float3(rn[i_p])/sd);
  }
  while (az.x>5||az.y>5||az.z>5);

  //calculate forces
  f[i_p] = {0.0,0.0,0.0,0.0};
  calc_bf(N,i_p,r,f);
  calc_cf(R,i_p,r,f);
  calc_all_srf<I>(eps,i_p,r,f,srg_p);

  //calculate extra position
  er[i_p] = r[i_p]+f[i_p]*dt+rn[i_p];
}

//execute 2nd stage of the Runge-Kutta method
template <srint I> __global__ void exec_RK_2(
  const uint N, //number of particles
  const float R, //confinement radius
  const float eps, //particle energy
  float4 *r, //position array
  float4 *f, //force array
  float4 *er, //extra position array
  float4 *ef, //extra force array
  float4 *rn, //random number array
  sugrid *srg_p) //short-range grid pointer
{
  //calculate particle index
  uint i_p = blockIdx.x*blockDim.x+threadIdx.x; //particle index
  if (i_p>=N){ return;}

  //calculate forces
  ef[i_p] = {0.0,0.0,0.0,0.0};
  calc_bf(N,i_p,er,ef);
  calc_cf(R,i_p,er,ef);
  calc_all_srf<I>(eps,i_p,er,ef,srg_p);

  //calculate new position
  r[i_p] = r[i_p]+0.5*(ef[i_p]+f[i_p])*dt+rn[i_p];
}

//Host Functions

//chromatin simulation constructor
chrsim::chrsim(parmap &par) //parameters
  : chrdat(par)
  , fpf {par.get_val<uint>("frames_per_file",100)}
  , spf {par.get_val<uint>("steps_per_frame",1*2048)}
  , tpb {par.get_val<uint>("threads_per_block",256)}
  , sd {sqrtf(2.0*k_B*T*dt)}
  , srg(N,aco,2*ceil(R/aco))
{
  //check parameters
  if (!(1<=fpf&&fpf<10'000)){ throw error("frames_per_file out of range");}
  if (!(1<=spf&&spf<10'000)){ throw error("steps_per_frame out of range");}
  if (!(1<=tpb&&tpb<1'025)){ throw error("threads_per_block out of range");}
  std::string msg = "chrsim:"; //message
  msg += " fpf = "+cnfs(fpf,4,'0');
  msg += " spf = "+cnfs(spf,4,'0');
  msg += " tpb = "+cnfs(tpb,4,'0');
  logger::record(msg);

  //allocate device memory
  cuda_check(cudaMalloc(&er,N*sizeof(float4)));
  cuda_check(cudaMalloc(&ef,N*sizeof(float4)));
  cuda_check(cudaMalloc(&rn,N*sizeof(float4)));
  cuda_check(cudaMalloc(&vps,N*sizeof(prng)));
  cuda_check(cudaMalloc(&srg_p,sizeof(sugrid)));

  //copy short-range grid to device
  cuda_check(cudaMemcpy(srg_p,&srg,sizeof(sugrid),cudaMemcpyHostToDevice));

  //initialize PRNG
  init_ps<<<(N+tpb-1)/tpb,tpb>>>(N,vps,time(nullptr));
}

//chromatin simulation destructor
chrsim::~chrsim()
{
  //deallocate device memory
  cuda_check(cudaFree(er));
  cuda_check(cudaFree(ef));
  cuda_check(cudaFree(rn));
  cuda_check(cudaFree(vps));
  cuda_check(cudaFree(srg_p));
}

//generate a random initial condition
void chrsim::generate_initial_condition()
{
  //set random particle types
  set_particle_types();

  //perform a confined random walk
  perform_random_walk();

  //separate beads
  logger::record("bead separation begun");
  // while () //check beads are separated ---------------------------------------
  // {
  //   make_RK_iteration();
  // }
  cuda_check(cudaMemcpy(hr,r,N*sizeof(float4),cudaMemcpyDeviceToHost));
  logger::record("bead separation ended");

  //record success message
  logger::record("initial condition generated");
}

//save simulation state to binary file
void chrsim::save_checkpoint(std::ofstream &bin_out_f) //binary output file
{
  bin_out_f.write(reinterpret_cast<char *>(&i_f),sizeof(i_f));
  bin_out_f.write(reinterpret_cast<char *>(&t),sizeof(t));
  bin_out_f.write(reinterpret_cast<char *>(hr),N*sizeof(float4));
  logger::record("simulation checkpoint saved");
}

//load simulation state from binary file
void chrsim::load_checkpoint(std::ifstream &bin_inp_f) //binary input file
{
  bin_inp_f.read(reinterpret_cast<char *>(&i_f),sizeof(i_f));
  bin_inp_f.read(reinterpret_cast<char *>(&t),sizeof(t));
  bin_inp_f.read(reinterpret_cast<char *>(hr),N*sizeof(float4));
  cuda_check(cudaMemcpy(r,hr,N*sizeof(float4),cudaMemcpyHostToDevice));
  logger::record("simulation checkpoint loaded");
}

//run simulation and write trajectory to binary file
void chrsim::run_simulation(std::ofstream &bin_out_f) //binary output file
{
  logger::record("simulation begun");
  for (uint ffi = 0; ffi<fpf; ++ffi) //file frame index
  {
    logger::show_prog_pc(100.0*ffi/fpf);
    for (uint fsi = 0; fsi<spf; ++fsi) //frame step index
    {
      make_RK_iteration();
    }
    cuda_check(cudaMemcpy(hr,r,N*sizeof(float4),cudaMemcpyDeviceToHost));
    ++i_f; t += spf*dt;
    write_frame_bin(bin_out_f);
  }
  logger::record("simulation ended");
}

//set random particle types
void chrsim::set_particle_types()
{
  //initialize host PRNG
  curandGenerator_t gen; //host PRNG
  curandCreateGeneratorHost(&gen,CURAND_RNG_PSEUDO_DEFAULT);
  curandSetPseudoRandomGeneratorSeed(gen,time(nullptr));

  //set particle types randomly
  float ran; //random number in (0,1]
  for (uint i_p = 0; i_p<N; ++i_p) //particle index
  {
    curandGenerateUniform(gen,&ran,1);
    if (ran<0.5){ hpt[i_p] = LAD;}
    else{ hpt[i_p] = LND;}
  }

  //copy host particle type array to device
  cuda_check(cudaMemcpy(pt,hpt,N*sizeof(ptype),cudaMemcpyHostToDevice));

  //free host PRNG
  curandDestroyGenerator(gen);
}

//perform a confined random walk
void chrsim::perform_random_walk()
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
  hr[0] = {0.0,0.0,0.0,0.0};
  curandGenerateUniform(gen,&ran,1); theta = acos(1.0-2.0*ran);
  curandGenerateUniform(gen,&ran,1); phi = 2.0*M_PI*ran;
  ran_dir = {sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)};
  old_dir = ran_dir;

  //place the rest of particles
  uint att = 0; //number of attempts
  for (uint i_p = 1; i_p<N; ++i_p) //particle index
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
    hr[i_p] = make_float4(len_b*new_dir)+hr[i_p-1];

    //check if new position is acceptable
    bool p_a = true; //position is acceptable
    if (!isfinite(hr[i_p].x)){ p_a = false;}
    if (!isfinite(hr[i_p].y)){ p_a = false;}
    if (!isfinite(hr[i_p].z)){ p_a = false;}
    for (uint j_p = 0; j_p<(i_p-1); ++j_p) //secondary particle index//remove ---
    {
      float dpp; //particle-particle distance
      dpp = length(make_float3(hr[j_p]-hr[i_p]));
      if (dpp<1.0){ p_a = false; break;}
    }
    float d_r; //radial distance to origin
    d_r = length(make_float3(hr[i_p]));
    if (((R+0.5)-d_r)<1.0){ p_a = false;}

    if (p_a) //continue
    {
      old_dir = new_dir;
      att = 0;
    }
    else //go back
    {
      ++att;
      if (att>1024){ i_p = 1;}
      else{ --i_p;}
    }
  }

  //copy host position array to device
  cuda_check(cudaMemcpy(r,hr,N*sizeof(float4),cudaMemcpyHostToDevice));

  //free host PRNG
  curandDestroyGenerator(gen);
}

//make one iteration of the Runge-Kutta method
void chrsim::make_RK_iteration()
{
  srg.generate_arrays(tpb,r);
  exec_RK_1<WFI><<<(N+tpb-1)/tpb,tpb>>>(N,R,eps,r,f,er,sd,rn,vps,srg_p);
  srg.generate_arrays(tpb,er);
  exec_RK_2<WFI><<<(N+tpb-1)/tpb,tpb>>>(N,R,eps,r,f,er,ef,rn,srg_p);
}

} //namespace mmc

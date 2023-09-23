//Includes

#include "chrsim.cuh" //chromatin simulation

#include <time.h> //time utilities library

//Namespace

namespace mmc //Marco Mendívil Carboni
{

//Constants

static constexpr float dt = 1.0/2048; //timestep
static constexpr float rco = 1.122462; //LJ repulsive cutoff
static constexpr float aco = 2.713283; //LJ attractive cutoff

//Device Functions

//calculate bonded forces
inline __device__ void calc_bf(
  const int N, //number of particles
  int i_p, //particle index
  float4 *r, //position array
  float4 *f) //force array
{
  //declare auxiliary variables
  float3 vec[4]; //bond vectors
  float il[4]; //bond inverse lengths
  float cos[3]; //bond angle cosines
  float3 bf = {0.0,0.0,0.0}; //bonded forces

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
  float sig, //LJ particle size
  int i_p, //particle index
  float4 *r, //position array
  float4 *f) //force array
{
  //calculate auxiliary variables
  float d_r; //radial distance to origin
  d_r = length(make_float3(r[i_p]));
  float dwp = (R+sig/2)-d_r; //wall-particle distance
  if (dwp>(rco*sig)){ return;}

  //calculate confinement force
  float3 cf = make_float3(-r[i_p]/d_r); //confinement force
  float s6 = sig*sig*sig*sig*sig*sig; //sig to the sixth power
  float d6 = dwp*dwp*dwp*dwp*dwp*dwp; //dwp to the sixth power
  cf *= 4.0*(12.0*(s6*s6)/(d6*d6*dwp)-6.0*(s6)/(d6*dwp));

  //add result to force array
  f[i_p] += make_float4(cf);
}

// //calculate single Lennard-Jones force
// inline __device__ float calc_single_ljf(
//   float sig, //LJ particle size
//   int i_p, //particle index
//   int j_p, //secondary particle index
//   float4 *r, //position array
//   float3 &ljf) //Lennard-Jones forces
// {
//   //calculate particle particle distance
//   float3 vpp = make_float3(r[i_p]-r[j_p]);//particle particle vector
//   float dpp = length(vpp); //particle particle distance
//   if (dpp>(rco*sig)){ return ((dpp-rco*sig)/(1.25*l_0));}//tmp------------------change to aco

//   //calculate Lennard-Jones force
//   float s6 = sig*sig*sig*sig*sig*sig; //sig to the sixth power
//   float d6 = dpp*dpp*dpp*dpp*dpp*dpp; //dpp to the sixth power
//   ljf += 4.0*(12.0*(s6*s6)/(d6*d6*dpp)-6.0*(s6)/(d6*dpp))*vpp;
//   return 0;
// }

// //calculate all Lennard-Jones forces
// __device__ void calc_all_ljf(
//   const int N, //number of particles
//   float sig, //LJ particle size
//   int i_p, //particle index
//   float4 *r, //position array
//   float4 *f) //force array
// {
//   //declare auxiliary variables
//   int n_skp; //number of skippable particle
//   float3 ljf = {0.0,0.0,0.0}; //Lennard-Jones forces

//   //traverse polymer both ways
//   for( int j_p = i_p-2; j_p>=0; j_p -= 1+n_skp) //secondary particle index
//   {
//     n_skp = calc_single_ljf(sig,i_p,j_p,r,ljf);
//   }
//   for( int j_p = i_p+2; j_p<N; j_p += 1+n_skp) //secondary particle index
//   {
//     n_skp = calc_single_ljf(sig,i_p,j_p,r,ljf);
//   }

//   //add result to force array
//   f[i_p] += make_float4(ljf);
// }

//calculate Lennard-Jones forces with cell's particles
inline __device__ void calc_cell_ljf(
  float sig, //LJ particle size
  int i_p, //particle index
  int i_c, //cell index
  float4 *r, //position array
  sugrid *ljp, //LJ grid pointer
  float3 &ljf) //Lennard-Jones forces
{
  //declare auxiliary variables
  int j_p; //secondary particle index
  float s6 = sig*sig*sig*sig*sig*sig; //sig to the sixth power

  //check cell isn't empty
  if (ljp->beg[i_c]==0xffffffff){ return;}

  //range over cell's particles
  for (int i_a = ljp->beg[i_c]; i_a<ljp->end[i_c]; ++i_a) //array index
  {
    j_p = ljp->spi[i_a];
    if (abs(j_p-i_p)>1)
    {
      //calculate particle particle distance
      float3 vpp = make_float3(r[i_p]-r[j_p]); //particle particle vector
      float dpp = length(vpp); //particle particle distance
      if (dpp>(rco*sig)){ continue;}//tmp---------------------------------------change to aco

      //calculate Lennard-Jones force
      float d6 = dpp*dpp*dpp*dpp*dpp*dpp; //dpp to the sixth power
      ljf += 4.0*(12.0*(s6*s6)/(d6*d6*dpp)-6.0*(s6)/(d6*dpp))*vpp;
    }
  }
}

//calculate all Lennard-Jones forces
inline __device__ void calc_all_ljf(
  float sig, //LJ particle size
  int i_p, //particle index
  float4 *r, //position array
  float4 *f, //force array
  sugrid *ljp) //LJ grid pointer
{
  //calculate auxiliary variables
  const float csl = ljp->csl; //grid cell side length
  const int cps = ljp->cps; //grid cells per side
  int3 ir = floorf(make_float3(r[i_p])/csl); //integer coordinates
  int iofst = (cps/2)*(1+cps+cps*cps); //index offset
  float3 ljf = {0.0,0.0,0.0}; //Lennard-Jones forces

  //range over neighbouring cells
  int3 nir; //neighbour integer coordinates
  for (nir.x = ir.x-1; nir.x<ir.x+2; ++nir.x)
  {
    for (nir.y = ir.y-1; nir.y<ir.y+2; ++nir.y)
    {
      for (nir.z = ir.z-1; nir.z<ir.z+2; ++nir.z)
      {
        int nci = iofst+nir.x+nir.y*cps+nir.z*cps*cps; //neighbour cell index
        if( nci<0 || nci>=cps*cps*cps){ continue;}
        calc_cell_ljf(sig,i_p,nci,r,ljp,ljf);
      }
    }
  }

  //add result to force array
  f[i_p] += make_float4(ljf);
}

//Global Functions

//initialize PRNG state array
__global__ void init_ps(
  const int N, //number of particles
  prng *ps, //PRNG state array
  int pseed) //PRNG seed
{
  //calculate particle index
  int i_p = blockIdx.x*blockDim.x+threadIdx.x; //particle index
  if (i_p>=N){ return;}

  //initialize PRNG state
  curand_init(pseed,i_p,0,&ps[i_p]);
}

//calculate grid cell and particle indexes
__global__ void calc_indexes(
  const int N, //number of particles
  float4 *r, //position array
  sugrid *gp) //grid pointer
{
  //calculate particle index
  int i_p = blockIdx.x*blockDim.x+threadIdx.x; //particle index
  if (i_p>=N){ return;}
  gp->upi[i_p] = i_p;

  //calculate auxiliary variables
  const float csl = gp->csl; //grid cell side length
  const int cps = gp->cps; //grid cells per side
  int3 ir = floorf(make_float3(r[i_p])/csl); //integer coordinates
  int iofst = (cps/2)*(1+cps+cps*cps); //index offset

  //calculate grid cell index
  gp->uci[i_p] = iofst+ir.x+ir.y*cps+ir.z*cps*cps;
}

//set grid cells empty
__global__ void set_cells_empty(
  const uint gc, //number of grid cells
  sugrid *gp) //grid pointer
{
  //calculate array index
  int i_a = blockIdx.x*blockDim.x+threadIdx.x; //array index
  if (i_a>=gc){ return;}

  //beginning and end of grid cells
  gp->beg[i_a] = 0xffffffff;
  gp->end[i_a] = 0;
}

//find beginning and end of each grid cell
__global__ void find_cells_limits(
  const int N, //number of particles
  sugrid *gp) //grid pointer
{
  //calculate array index
  int i_a = blockIdx.x*blockDim.x+threadIdx.x; //array index
  if (i_a>=N){ return;}

  //set beginning and end of cells
  int ci_curr = gp->sci[i_a]; //current cell index
  if (i_a==0)
  {
    gp->beg[ci_curr] = i_a; return;
  }
  int ci_prev = gp->sci[i_a-1]; //previous cell index
  if (ci_prev!=ci_curr)
  {
    gp->beg[ci_curr] = i_a;
    gp->end[ci_prev] = i_a;
  }
  if (i_a==N-1)
  {
    gp->end[ci_curr] = i_a+1;
  }
}

//execute 1st stage of the Runge-Kutta method
__global__ void exec_RK_1(
  const int N, //number of particles
  const float R, //confinement radius
  float4 *r, //position array
  float4 *f, //force array
  float sig, //LJ particle size
  float4 *er, //extra position array
  float sd, //random number standard deviation
  float4 *rn, //random number array
  prng *ps, //PRNG state array
  sugrid *ljp) //LJ grid pointer
{
  //calculate particle index
  int i_p = blockIdx.x*blockDim.x+threadIdx.x; //particle index
  if (i_p>=N){ return;}

  //calculate random numbers
  rn[i_p] = sd*curand_normal4(&ps[i_p]);

  //calculate forces
  f[i_p] = {0.0,0.0,0.0,0.0};
  calc_bf(N,i_p,r,f);
  calc_cf(R,sig,i_p,r,f);
  calc_all_ljf(sig,i_p,r,f,ljp);

  //calculate extra position
  er[i_p] = r[i_p]+f[i_p]*dt/xi+rn[i_p]/xi;
}

//execute 2nd stage of the Runge-Kutta method
__global__ void exec_RK_2(
  const int N, //number of particles
  const float R, //confinement radius
  float4 *r, //position array
  float4 *f, //force array
  float sig, //LJ particle size
  float4 *er, //extra position array
  float4 *ef, //extra force array
  float4 *rn, //random number array
  sugrid *ljp) //LJ grid pointer
{
  //calculate particle index
  int i_p = blockIdx.x*blockDim.x+threadIdx.x; //particle index
  if (i_p>=N){ return;}

  //calculate forces
  ef[i_p] = {0.0,0.0,0.0,0.0};
  calc_bf(N,i_p,er,ef);
  calc_cf(R,sig,i_p,er,ef);
  calc_all_ljf(sig,i_p,er,ef,ljp);

  //calculate new position
  r[i_p] = r[i_p]+0.5*(ef[i_p]+f[i_p])*dt/xi+rn[i_p]/xi;
}

//Host Functions

//sorted uniform grid constructor
sugrid::sugrid(
    const uint N, //number of particles
    const float csl, //grid cell side length
    const uint cps) //grid cells per side
  : csl {csl}
  , cps {cps}
{
  //check parameters
  if (csl<0.0){ throw error("grid_cell_side_length out of range");}
  if (cps<1){ throw error("grid_cells_per_side out of range");}
  std::string msg = "sugrid:"; //message
  msg += " csl = "+cnfs(csl,6,'0',2);
  msg += " cps = "+cnfs(cps,5,'0');
  logger::record(msg);

  //calculate auxiliary variables
  const uint gc = cps*cps*cps; //number of grid cells

  //allocate arrays
  cuda_check(cudaMalloc(&uci,N*sizeof(int)));
  cuda_check(cudaMalloc(&sci,N*sizeof(int)));
  cuda_check(cudaMalloc(&upi,N*sizeof(int)));
  cuda_check(cudaMalloc(&spi,N*sizeof(int)));
  cuda_check(cudaMalloc(&beg,gc*sizeof(int)));
  cuda_check(cudaMalloc(&end,gc*sizeof(int)));
  cuda_check(cudaMalloc(&sr,N*sizeof(float4)));

  //allocate extra buffer
  cub::DeviceRadixSort::SortPairs(nullptr,ebs,uci,sci,upi,spi,N);
  cuda_check(cudaMalloc(&eb,ebs));
}

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
  cuda_check(cudaFree(sr));

  //deallocate extra buffer
  cuda_check(cudaFree(eb));
}

//sort pairs of indexes
void sugrid::sort_indexes(int N) //number of particles
{
  cub::DeviceRadixSort::SortPairs(eb,ebs,uci,sci,upi,spi,N);
}

//chromatin simulation constructor
chrsim::chrsim(parmap &par) //parameters
  : chrdat(par)
  , fpf {par.get_val<int>("frames_per_file",100)}
  , spf {par.get_val<int>("steps_per_frame",1*2048)}
  , tpb {par.get_val<int>("threads_per_block",256)}
  , sd {sqrtf(2.0*xi*k_B*T*dt)}
  , ljg(N,aco*sig,2*ceilf(R/(aco*sig)))
  , ljc {ljg.cps*ljg.cps*ljg.cps}
{
  //check parameters
  if (fpf<1){ throw error("frames_per_file out of range");}
  if (spf<1){ throw error("steps_per_frame out of range");}
  if (tpb<1){ throw error("threads_per_block out of range");}
  std::string msg = "chrsim:"; //message
  msg += " fpf = "+cnfs(fpf,5,'0');
  msg += " spf = "+cnfs(spf,5,'0');
  msg += " tpb = "+cnfs(tpb,5,'0');
  logger::record(msg);

  //allocate arrays
  cuda_check(cudaMalloc(&er,N*sizeof(float4)));
  cuda_check(cudaMalloc(&ef,N*sizeof(float4)));
  cuda_check(cudaMalloc(&rn,N*sizeof(float4)));
  cuda_check(cudaMallocManaged(&ps,N*sizeof(prng)));

  //allocate structures
  cuda_check(cudaMalloc(&ljp,sizeof(sugrid)));

  //copy LJ grid to device
  cuda_check(cudaMemcpy(ljp,&ljg,sizeof(sugrid),cudaMemcpyDefault));

  //initialize PRNG
  init_ps<<<(N+tpb-1)/tpb,tpb>>>(N,ps,time(nullptr));
  cuda_check(cudaDeviceSynchronize());
}

//chromatin simulation destructor
chrsim::~chrsim()
{
  //deallocate arrays
  cuda_check(cudaFree(er));
  cuda_check(cudaFree(ef));
  cuda_check(cudaFree(rn));
  cuda_check(cudaFree(ps));

  //deallocate structures
  cuda_check(cudaFree(ljp));
}

//generate a random initial condition
void chrsim::generate_initial_condition()
{
  //initialize host PRNG
  curandGenerator_t gen; //host PRNG
  curandCreateGeneratorHost(&gen,CURAND_RNG_PSEUDO_DEFAULT);
  curandSetPseudoRandomGeneratorSeed(gen,time(nullptr));

  //set random particle types
  set_particle_types(gen);

  //reduce sigma for random walk
  sig = 1.0/2;

  //perform a confined random walk
  perform_random_walk(gen);

  //expand beads
  while (sig<1.0)
  {
    make_RK_iteration();
    sig += dt/(32*sig*sig);
  }
  cuda_check(cudaDeviceSynchronize());

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
  bin_out_f.write(reinterpret_cast<char *>(r),N*sizeof(float4));
  bin_out_f.write(reinterpret_cast<char *>(ps),N*sizeof(prng));
  logger::record("simulation checkpoint saved");
}

//load simulation state from binary file
void chrsim::load_checkpoint(std::ifstream &bin_inp_f) //binary input file
{
  bin_inp_f.read(reinterpret_cast<char *>(&i_f),sizeof(i_f));
  bin_inp_f.read(reinterpret_cast<char *>(&t),sizeof(t));
  bin_inp_f.read(reinterpret_cast<char *>(r),N*sizeof(float4));
  bin_inp_f.read(reinterpret_cast<char *>(ps),N*sizeof(prng));
  logger::record("simulation checkpoint loaded");
}

//run simulation and trajectory to binary file
void chrsim::run_simulation(std::ofstream &bin_out_f) //binary output file
{
  for (int ffi = 0; ffi<fpf; ++ffi) //file frame index
  {
    float prog_pc = (100.0*ffi)/(fpf); //progress percentage
    logger::show_prog_pc(prog_pc);
    for (int fsi = 0; fsi<spf; ++fsi) //frame step index
    {
      make_RK_iteration();
    }
    cuda_check(cudaDeviceSynchronize());
    ++i_f; t += spf*dt;
    write_frame_bin(bin_out_f);
  }
}

//set random particle types
void chrsim::set_particle_types(curandGenerator_t &gen) //host PRNG
{
  float ran; //random number in (0,1]
  for (int i_p = 0; i_p<N; ++i_p) //particle index
  {
    curandGenerateUniform(gen,&ran,1);
    if (ran<0.5){ pt[i_p] = LAD;}
    else{ pt[i_p] = non_LAD;}
  }
}

//perform a confined random walk
void chrsim::perform_random_walk(curandGenerator_t &gen) //host PRNG
{
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
  r[0] = {0.0,0.0,0.0,0.0};
  curandGenerateUniform(gen,&ran,1); theta = acos(1.0-2.0*ran);
  curandGenerateUniform(gen,&ran,1); phi = 2.0*M_PI*ran;
  ran_dir = {sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)};
  old_dir = ran_dir;

  //place the rest of particles
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
    r[i_p] = make_float4(len_b*new_dir)+r[i_p-1];

    //check if new position is acceptable
    bool p_a = true; //position is acceptable
    if (!isfinite(r[i_p].x)){ p_a = false;}
    if (!isfinite(r[i_p].y)){ p_a = false;}
    if (!isfinite(r[i_p].z)){ p_a = false;}
    for (int j_p = 0; j_p<(i_p-1); ++j_p) //secondary particle index
    {
      float dpp; //particle-particle distance
      dpp = length(make_float3(r[j_p]-r[i_p]));
      if (dpp<sig){ p_a = false; break;}
    }
    float d_r; //radial distance to origin
    d_r = length(make_float3(r[i_p]));
    if (((R+sig/2)-d_r)<sig){ p_a = false;}

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
}

// //generate LJ grid data
// void chrsim::generate_ljg(float4 *r) //position array
// {
//   calc_indexes<<<(N+tpb-1)/tpb,tpb>>>(N,r,ljp);
//   ljg.sort_indexes(N);
//   set_cells_empty<<<(ljc+tpb-1)/tpb,tpb>>>(ljc,ljp);
//   find_cells_limits<<<(N+tpb-1)/tpb,tpb>>>(N,ljp);
// }

//make one iteration of the Runge-Kutta method
void chrsim::make_RK_iteration()
{
  calc_indexes<<<(N+tpb-1)/tpb,tpb>>>(N,r,ljp);
  ljg.sort_indexes(N);
  set_cells_empty<<<(ljc+tpb-1)/tpb,tpb>>>(ljc,ljp);
  find_cells_limits<<<(N+tpb-1)/tpb,tpb>>>(N,ljp);
  exec_RK_1<<<(N+tpb-1)/tpb,tpb>>>(N,R,r,f,sig,er,sd,rn,ps,ljp);
  calc_indexes<<<(N+tpb-1)/tpb,tpb>>>(N,er,ljp);
  ljg.sort_indexes(N);
  set_cells_empty<<<(ljc+tpb-1)/tpb,tpb>>>(ljc,ljp);
  find_cells_limits<<<(N+tpb-1)/tpb,tpb>>>(N,ljp);
  exec_RK_2<<<(N+tpb-1)/tpb,tpb>>>(N,R,r,f,sig,er,ef,rn,ljp);
}

} //namespace mmc

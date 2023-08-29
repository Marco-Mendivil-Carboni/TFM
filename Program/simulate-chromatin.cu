//Libraries

#include <cstdio> //standard input and output library
#include <cstdlib> //standard general utilities library
#include <cmath> //mathematical functions library
#include <ctime> //time utilities library

#include <glob.h> //pathname pattern matching types

#include <curand_kernel.h> //cuRAND device functions

//Constants

constexpr float xi  {1.000000}; //damping coefficient
constexpr float k_B {0.001120}; //Boltzmann constant
constexpr float l_0 {1.000000}; //natural length of bonds
constexpr float k_e {100.0000}; //elastic constant
constexpr float k_b {2.000000}; //bending constant
constexpr float r_c {1.122462}; //VdW cutoff radius
constexpr float dt  {1.0/2048}; //timestep

constexpr int n_s {1*2048}; //MD steps between frames

//Structs

struct sim_par //simulation parameters
{
  float T; //temperature
  int N; //number of particles
  float R; //radius of sphere
  int F; //frames per file
};
typedef struct sim_par sim_par;

using PRNGstate = curandStatePhilox4_32_10;

//Functions

inline void cuda_check( cudaError_t result)
{
  if( result!=cudaSuccess){ std::fprintf(stderr,"CUDA Error: %s\n",cudaGetErrorString(result)); exit(-1);}
}

inline void print_time( FILE* f)
{
  time_t rt = time(nullptr); struct tm *rti = localtime(&rt);
  std::fprintf(f,"%02d:%02d:%02d ",rti->tm_hour,rti->tm_min,rti->tm_sec);
}

void read_parameters( sim_par& sp, FILE* f)
{
  if( std::fscanf(f,"T\t%f\n",&(sp.T))!=1){ std::fprintf(stderr,"Error reading parameters file.\n"); exit(-1);}
  if( std::fscanf(f,"N\t%d\n",&(sp.N))!=1){ std::fprintf(stderr,"Error reading parameters file.\n"); exit(-1);}
  if( std::fscanf(f,"R\t%f\n",&(sp.R))!=1){ std::fprintf(stderr,"Error reading parameters file.\n"); exit(-1);}
  if( std::fscanf(f,"F\t%d\n",&(sp.F))!=1){ std::fprintf(stderr,"Error reading parameters file.\n"); exit(-1);}
  if( (sp.T)<__FLT_MIN__){ std::fprintf(stderr,"T must be positive.\n"); exit(-1);}
  if( (sp.N)<__FLT_MIN__){ std::fprintf(stderr,"N must be positive.\n"); exit(-1);}
  if( (sp.R)<__FLT_MIN__){ std::fprintf(stderr,"R must be positive.\n"); exit(-1);}
  if( (sp.F)<__FLT_MIN__){ std::fprintf(stderr,"F must be positive.\n"); exit(-1);}
}

void generate_initial_configuration( int N, float T, float R, float sig, float *r)
{
  curandGenerator_t gen;
  curandCreateGeneratorHost(&gen,CURAND_RNG_PSEUDO_DEFAULT);
  curandSetPseudoRandomGeneratorSeed(gen,time(nullptr));
  float random, theta, varphi, bondlen, bondangle;
  float dir_old[3], dir_new[3], perpdir[3], perpdirnorm;
  float beta = 1.0/(k_B*T);
  curandGenerateUniform(gen,&random,1); theta = acos(1.0-2.0*random); 
  curandGenerateUniform(gen,&random,1); varphi = 2.0*M_PI*random;
  dir_old[0]=sin(theta)*cos(varphi);
  dir_old[1]=sin(theta)*sin(varphi);
  dir_old[2]=cos(theta);
  r[0]=r[1]=r[2]=0.0;
  int n_failures = 0;
  for( int i_p = 1; i_p<N; ++i_p)
  {
    curandGenerateUniform(gen,&random,1); theta = acos(1.0-2.0*random);
    curandGenerateUniform(gen,&random,1); varphi = 2.0*M_PI*random;
    perpdir[0] = dir_old[1]*cos(theta)-dir_old[2]*sin(theta)*sin(varphi);
    perpdir[1] = dir_old[2]*sin(theta)*cos(varphi)-dir_old[0]*cos(theta);
    perpdir[2] = dir_old[0]*sin(theta)*sin(varphi)-dir_old[1]*sin(theta)*cos(varphi);
    perpdirnorm = sqrt(perpdir[0]*perpdir[0]+perpdir[1]*perpdir[1]+perpdir[2]*perpdir[2]);
    perpdir[0] /= perpdirnorm; perpdir[1] /= perpdirnorm; perpdir[2] /= perpdirnorm;
    curandGenerateUniform(gen,&random,1);
    bondangle = acos(1.0+log(1.0-(1.0-exp(-2.0*(k_b*beta)))*random)/(k_b*beta));
    dir_new[0] = dir_old[0]*cos(bondangle)+perpdir[0]*sin(bondangle);
    dir_new[1] = dir_old[1]*cos(bondangle)+perpdir[1]*sin(bondangle);
    dir_new[2] = dir_old[2]*cos(bondangle)+perpdir[2]*sin(bondangle);
    curandGenerateUniform(gen,&random,1);
    bondlen = l_0+sqrt(2.0/(k_e*beta))*erfinv(2.0*random-1.0);
    r[3*i_p+0] = bondlen*dir_new[0]+r[3*(i_p-1)+0];
    r[3*i_p+1] = bondlen*dir_new[1]+r[3*(i_p-1)+1];
    r[3*i_p+2] = bondlen*dir_new[2]+r[3*(i_p-1)+2];
    int accept = 1;
    if( !isfinite(r[3*i_p+0])){ accept = 0;}
    if( !isfinite(r[3*i_p+1])){ accept = 0;}
    if( !isfinite(r[3*i_p+2])){ accept = 0;}
    for( int j_p = 0; j_p<(i_p-1); ++j_p)
    {
      float dist = 0.0;
      dist += (r[3*j_p+0]-r[3*i_p+0])*(r[3*j_p+0]-r[3*i_p+0]);
      dist += (r[3*j_p+1]-r[3*i_p+1])*(r[3*j_p+1]-r[3*i_p+1]);
      dist += (r[3*j_p+2]-r[3*i_p+2])*(r[3*j_p+2]-r[3*i_p+2]);
      dist = sqrt(dist);
      if( dist<(r_c*sig)){ accept = 0; break;}
    }
    float d_r = 0.0;
    for( int i_c = 0; i_c<3; ++i_c)
    {
      d_r += r[3*i_p+i_c]*r[3*i_p+i_c];
    }
    d_r = sqrt(d_r);
    if( (R-d_r)<(r_c*sig)){ accept = 0;}
    if( accept)
    {
      dir_old[0] = dir_new[0];
      dir_old[1] = dir_new[1];
      dir_old[2] = dir_new[2];
      n_failures = 0;
    }
    else
    {
      ++n_failures;
      if( n_failures>1024){ i_p = 0;}
      else{ i_p--;}
    }
  }
  curandDestroyGenerator(gen);
}

void write_initial_configuration( int N, float *r, FILE *f)
{
  std::fprintf(f,"Chromatin simulation, t=0.0\n");
  std::fprintf(f,"%5d\n",N);
  for( int i_p = 0; i_p<N; ++i_p)
  {
    std::fprintf(f,"%5d%-5s%5s%5d",i_p+1,"X","X",i_p+1);
    std::fprintf(f,"%8.3f%8.3f%8.3f\n",r[3*i_p+0],r[3*i_p+1],r[3*i_p+2]);
  }
  std::fprintf(f,"%10.5f%10.5f%10.5f\n",0.0,0.0,0.0);
}

void write_trajectory_positions( int N, float *r, float t, int i_f, FILE *f)
{
  int magickvalue = 1993;
  fwrite(&magickvalue,sizeof(int),1,f);
  char trrversion[] = "GMX_trn_file";
  int len_s_a = sizeof(trrversion);
  int len_s_b = sizeof(trrversion)-1;
  fwrite(&len_s_a,sizeof(int),1,f);
  fwrite(&len_s_b,sizeof(int),1,f);
  fwrite(trrversion,sizeof(char),sizeof(trrversion)-1,f);
  int zero = 0;
  for( int i = 0; i<7; ++i)
  {
    fwrite(&zero,sizeof(int),1,f);
  }
  int x_size = 3*N*sizeof(float);
  fwrite(&x_size,sizeof(int),1,f);
  int v_size = 0;
  fwrite(&v_size,sizeof(int),1,f);
  int f_size = 0;
  fwrite(&f_size,sizeof(int),1,f);
  int natoms = N;
  fwrite(&natoms,sizeof(int),1,f);
  int step = i_f;
  fwrite(&step,sizeof(int),1,f);
  float time = t;
  fwrite(&zero,sizeof(int),1,f);
  fwrite(&time,sizeof(float),1,f);
  fwrite(&zero,sizeof(int),1,f);
  fwrite(r,sizeof(float),3*N,f);
}

void save_checkpoint( int N, float *r, float *t, int n_threads, PRNGstate *state, int *tpf_idx, FILE *f)
{
  fwrite(tpf_idx,sizeof(int),1,f);
  fwrite(t,sizeof(float),1,f);
  fwrite(r,sizeof(float),3*N,f);
  fwrite(state,sizeof(PRNGstate),n_threads,f);
}

void load_checkpoint( int N, float *r, float *t, int n_threads, PRNGstate *state, int *tpf_idx, FILE *f)
{
  fread(tpf_idx,sizeof(int),1,f);
  fread(t,sizeof(float),1,f);
  fread(r,sizeof(float),3*N,f);
  fread(state,sizeof(PRNGstate),n_threads,f);
  *tpf_idx+=1;
}

//Kernels

__global__
void setup_PRNG( int seed, PRNGstate *state)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  curand_init(seed, i_p, 0, &state[i_p]);
}

__global__
void call_PRNG( float c_rn, float *nrn, PRNGstate *state)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  for( int i_c = 0; i_c<3; ++i_c)
  {
    nrn[3*i_p+i_c] = c_rn*curand_normal(&state[i_p]);
  }
}

__global__
void calc_extern_f( int N, float *f_c, float *f)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  if( i_p<N)
  {
    for( int i_c = 0; i_c<3; ++i_c)
    {
      f[3*i_p+i_c] = f_c[3*i_p+i_c];
    }
  }
}

__global__
void calc_sphere_f( int N, float R, float sig, float *r, float *f)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  float s2 = sig*sig;
  if( i_p<N)
  {
    float k_LJ;
    float d_r = 0.0;
    for( int i_c = 0; i_c<3; ++i_c)
    {
      d_r += r[3*i_p+i_c]*r[3*i_p+i_c];
    }
    d_r = sqrt(d_r);
    float dsp = (R-d_r);
    float d2 = dsp*dsp;
    if( d2<(r_c*r_c*s2))
    {
      k_LJ = 48.0*(s2*s2*s2*s2*s2*s2)/(d2*d2*d2*d2*d2*d2*d2)-24.0*(s2*s2*s2)/(d2*d2*d2*d2);
      for( int i_c = 0; i_c<3; ++i_c)
      {
        f[3*i_p+i_c] += k_LJ*(-dsp/d_r)*r[3*i_p+i_c];
      }
    }
  }
}

__global__ 
void calc_bonds( int N, float *r, float *b, float *invlen)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  if( i_p<N-1)
  {
    invlen[i_p+2] = 0.0;
    for( int i_c = 0; i_c<3; ++i_c)
    {
      b[3*(i_p+2)+i_c] = r[3*(i_p+1)+i_c]-r[3*i_p+i_c];
      invlen[i_p+2] += b[3*(i_p+2)+i_c]*b[3*(i_p+2)+i_c];
    }
    invlen[i_p+2] = 1.0/sqrt(invlen[i_p+2]);
  }
}

__global__
void calc_cosines( int N, float *b, float *invlen, float *cosine)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  if( i_p<N-2)
  {
    cosine[i_p+3] = 0.0;
    for( int i_c = 0; i_c<3; ++i_c)
    {
      cosine[i_p+3] += b[3*(i_p+3)+i_c]*b[3*(i_p+2)+i_c];
    }
    cosine[i_p+3] *= invlen[i_p+3]*invlen[i_p+2];
  }
}

__global__
void calc_intern_f( int N, float *b, float *invlen, float *cosine, float *f)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  if( i_p<N)
  {
    for( int i_c = 0; i_c<3; ++i_c)
    {
      f[3*i_p+i_c] += k_e*(1.0-l_0*invlen[i_p+1])*(-b[3*(i_p+1)+i_c]);

      f[3*i_p+i_c] += k_e*(1.0-l_0*invlen[i_p+2])*(+b[3*(i_p+2)+i_c]);

      f[3*i_p+i_c] += k_b*(+b[3*(i_p+0)+i_c])*invlen[i_p+1]*invlen[i_p+0];

      f[3*i_p+i_c] += k_b*(-cosine[i_p+2]-cosine[i_p+1])*b[3*(i_p+1)+i_c]*invlen[i_p+1]*invlen[i_p+1];

      f[3*i_p+i_c] += k_b*(+b[3*(i_p+2)+i_c]-b[3*(i_p+1)+i_c])*invlen[i_p+2]*invlen[i_p+1];

      f[3*i_p+i_c] += k_b*(+cosine[i_p+3]+cosine[i_p+2])*b[3*(i_p+2)+i_c]*invlen[i_p+2]*invlen[i_p+2];

      f[3*i_p+i_c] += k_b*(-b[3*(i_p+3)+i_c])*invlen[i_p+3]*invlen[i_p+2];
    }
  }
}

__device__
int calc_LJ_f( int N, float sig, float *r, int i_p, int j_p, float *f)
{
  float k_LJ;
  float d2 = 0.0;
  float s2 = sig*sig;
  for( int i_c = 0; i_c<3; ++i_c)
  {
    d2 += (r[3*i_p+i_c]-r[3*j_p+i_c])*(r[3*i_p+i_c]-r[3*j_p+i_c]);
  }
  if( d2<(r_c*r_c*s2))
  {
    k_LJ = 48.0*(s2*s2*s2*s2*s2*s2)/(d2*d2*d2*d2*d2*d2*d2)-24.0*(s2*s2*s2)/(d2*d2*d2*d2);
    for( int i_c = 0; i_c<3; ++i_c)
    {
      f[3*i_p+i_c] += k_LJ*(r[3*i_p+i_c]-r[3*j_p+i_c]);
    }
    return 0;
  }
  else
  {
    return ((sqrt(d2)-r_c*sig)/(1.25*l_0));
  }
}

__global__
void calc_exclvol_f( int N, float sig, float *r, float *f)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  if( i_p<N)
  {
    int skip;
    for( int j_p = i_p-2; j_p>=0; j_p -= 1+skip)
    {
      skip = calc_LJ_f(N,sig,r,i_p,j_p,f);
    }
    for( int j_p = i_p+2; j_p<N; j_p += 1+skip)
    {
      skip = calc_LJ_f(N,sig,r,i_p,j_p,f);
    }
  }
}

__global__
void RK_stage_1( int N, float *r_1, float *r_2, float *f_2, float *nrn)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  if( i_p<N)
  {
    for( int i_c = 0; i_c<3; ++i_c)
    {
      r_1[3*i_p+i_c] = r_2[3*i_p+i_c]+f_2[3*i_p+i_c]*dt/xi+nrn[3*i_p+i_c]/xi;
    }
  }
}

__global__
void RK_stage_2( int N, float *r_2, float *f_1, float *f_2, float *nrn)
{
  int i_p = blockIdx.x * blockDim.x + threadIdx.x;
  if( i_p<N)
  {
    for( int i_c = 0; i_c<3; ++i_c)
    {
      r_2[3*i_p+i_c] = r_2[3*i_p+i_c]+0.5*(f_1[3*i_p+i_c]+f_2[3*i_p+i_c])*dt/xi+nrn[3*i_p+i_c]/xi;
    }
  }
}

int main( int argc, const char** argv)
{
  if( argc<2){ std::fprintf(stderr,"Input missing.\n"); return EXIT_FAILURE;}
  if( argc>3){ std::fprintf(stderr,"Too many arguments.\n"); return EXIT_FAILURE;}
  if( sizeof(argv[1])>128){ std::fprintf(stderr,"Directory name too long.\n"); return EXIT_FAILURE;}

  char sim_dir[128];
  std::snprintf(sim_dir,sizeof(sim_dir),"%s",argv[1]);

  FILE* file_i1;
  FILE* file_o1;
  FILE* logfile;

  char filename[256];

  std::snprintf(filename,sizeof(filename),"%s/current-progress.log",sim_dir);
  logfile = std::fopen(filename,"wt");
  if( logfile==nullptr){ std::fprintf(stderr,"Error opening the current progress file.\n"); return EXIT_FAILURE;}

  //Simulation parameters and variables

  sim_par sp;

  std::snprintf(filename,sizeof(filename),"%s/adjustable-parameters.dat",sim_dir);
  file_i1 = std::fopen(filename,"rt");
  if( file_i1==nullptr){ std::fprintf(stderr,"Error opening the adjustable parameters file.\n"); return EXIT_FAILURE;}
  read_parameters(sp,file_i1);
  std::fclose(file_i1);

  float cvf = sp.N*pow(0.5/(sp.R-0.5),3); //Calculate the chromatin volume fraction
  if( cvf>0.5){ std::fprintf(stderr,"Chromatin volume fraction too high (above 0.5).\n"); return EXIT_FAILURE;}

  print_time(logfile); std::fprintf(logfile,"Adjustable parameters file read.\n");
  std::fprintf(logfile,"T=%05.1f N=%04d R=%05.2f F=%04d\n",sp.T,sp.N,sp.R,sp.F); std::fflush(logfile);

  size_t threads_block = 256;
  size_t n_blocks = (sp.N+threads_block-1)/threads_block;
  size_t n_threads = n_blocks*threads_block;

  float c_rn = sqrt(2.0*xi*k_B*sp.T*dt);

  float *r_2;
  float *r_1;

  float *f_2;
  float *f_1;

  float *nrn;
  PRNGstate *state;

  float *b;
  float *invlen;
  float *cosine;

  float *f_c;

  //Memory allocation

  cuda_check( cudaMallocManaged( &r_2, 3*sp.N*sizeof(float)));
  cuda_check( cudaMallocManaged( &r_1, 3*sp.N*sizeof(float)));

  cuda_check( cudaMallocManaged( &f_2, 3*sp.N*sizeof(float)));
  cuda_check( cudaMallocManaged( &f_1, 3*sp.N*sizeof(float)));

  cuda_check( cudaMallocManaged( &nrn, 3*n_threads*sizeof(float)));
  cuda_check( cudaMallocManaged( &state, n_threads*sizeof(PRNGstate)));

  cuda_check( cudaMallocManaged( &b, 3*(sp.N+3)*sizeof(float)));
  cuda_check( cudaMallocManaged( &invlen, (sp.N+3)*sizeof(float)));
  cuda_check( cudaMallocManaged( &cosine, (sp.N+4)*sizeof(float)));

  cuda_check( cudaMallocManaged( &f_c, 3*sp.N*sizeof(float)));

  //Exceptions for the polymer ends

  b[0] = b[1] = b[2] = 0.0;
  b[3] = b[4] = b[5] = 0.0;
  b[3*(sp.N+2)+0] = b[3*(sp.N+2)+1] = b[3*(sp.N+2)+2] = 0.0;
  b[3*(sp.N+1)+0] = b[3*(sp.N+1)+1] = b[3*(sp.N+1)+2] = 0.0;

  invlen[0] = invlen[1] = 0.0;
  invlen[sp.N+2] = invlen[sp.N+1] = 0.0;

  cosine[0] = cosine[1] = cosine[2] = 0.0;
  cosine[sp.N+3] = cosine[sp.N+2] = cosine[sp.N+1] = 0.0;

  //Constant force

  for( int i_p = 0; i_p<sp.N; ++i_p)
  {
    for( int i_c = 0; i_c<3; ++i_c)
    {
      f_c[3*i_p+i_c] = 0.0;
    }
  }

  //Initialization

  setup_PRNG<<<n_blocks,threads_block>>>(time(nullptr),state);
  cuda_check( cudaDeviceSynchronize());

  int sim_idx = 0;
  int tpf_idx = 0;
  float t = 0.0;

  if( argc==3)
  {
    sim_idx = atoi(argv[2]);

    std::snprintf(filename,sizeof(filename),"%s/simulation-checkpoint-%03d.bin",sim_dir,sim_idx);
    file_i1 = std::fopen(filename,"rb");
    if( file_i1==nullptr){ std::fprintf(stderr,"Error opening the simulation checkpoint file.\n"); return EXIT_FAILURE;}
    load_checkpoint(sp.N,r_2,&t,n_threads,state,&tpf_idx,file_i1);
    std::fclose(file_i1);

    print_time(logfile); std::fprintf(logfile,"Simulation checkpoint file loaded.\n");
    std::fprintf(logfile,"sim_idx=%03d tpf_idx=%03d\n",sim_idx,tpf_idx); std::fflush(logfile);

    std::snprintf(filename,sizeof(filename),"%s/trajectory-positions-%03d-%03d.trr",sim_dir,sim_idx,tpf_idx);
    file_o1 = std::fopen(filename,"wb");
    if( file_o1==nullptr){ std::fprintf(stderr,"Error opening the trajectory positions file.\n"); return EXIT_FAILURE;}
  }
  else
  {
    glob_t prev_sims;
    std::snprintf(filename,sizeof(filename),"%s/initial-configuration-*",sim_dir);
    if( glob(filename,0,nullptr,&prev_sims)==0)
    {
      sim_idx = prev_sims.gl_pathc;
    }
    globfree(&prev_sims);

    print_time(logfile); std::fprintf(logfile,"New simulation started.\n");
    std::fprintf(logfile,"sim_idx=%03d tpf_idx=%03d\n",sim_idx,tpf_idx); std::fflush(logfile);

    float sig = 1.0/2.0;

    generate_initial_configuration(sp.N,sp.T,sp.R,sig,r_2);

    print_time(logfile); std::fprintf(logfile,"Initial configuration generated.\n"); std::fflush(logfile);

    while( sig<1.0)
    {
      call_PRNG<<<n_blocks,threads_block>>>(c_rn,nrn,state);

      calc_extern_f<<<n_blocks,threads_block>>>(sp.N,f_c,f_2);
      calc_sphere_f<<<n_blocks,threads_block>>>(sp.N,sp.R,sig,r_2,f_2);
      calc_bonds<<<n_blocks,threads_block>>>(sp.N,r_2,b,invlen);
      calc_cosines<<<n_blocks,threads_block>>>(sp.N,b,invlen,cosine);
      calc_intern_f<<<n_blocks,threads_block>>>(sp.N,b,invlen,cosine,f_2);
      // calc_exclvol_f<<<n_blocks,threads_block>>>(sp.N,sig,r_2,f_2);

      RK_stage_1<<<n_blocks,threads_block>>>(sp.N,r_1,r_2,f_2,nrn);

      calc_extern_f<<<n_blocks,threads_block>>>(sp.N,f_c,f_1);
      calc_sphere_f<<<n_blocks,threads_block>>>(sp.N,sp.R,sig,r_1,f_1);
      calc_bonds<<<n_blocks,threads_block>>>(sp.N,r_1,b,invlen);
      calc_cosines<<<n_blocks,threads_block>>>(sp.N,b,invlen,cosine);
      calc_intern_f<<<n_blocks,threads_block>>>(sp.N,b,invlen,cosine,f_1);
      // calc_exclvol_f<<<n_blocks,threads_block>>>(sp.N,sig,r_1,f_1);

      RK_stage_2<<<n_blocks,threads_block>>>(sp.N,r_2,f_1,f_2,nrn);

      sig *= 1.0+1.0/8192.0;
    }

    cuda_check( cudaDeviceSynchronize());

    print_time(logfile); std::fprintf(logfile,"Bead expansion finished.\n"); std::fflush(logfile);

    std::snprintf(filename,sizeof(filename),"%s/initial-configuration-%03d.gro",sim_dir,sim_idx);
    file_o1 = std::fopen(filename,"wt");
    if( file_o1==nullptr){ std::fprintf(stderr,"Error opening the initial configuration file.\n"); return EXIT_FAILURE;}
    write_initial_configuration(sp.N,r_2,file_o1);
    std::fclose(file_o1);

    std::snprintf(filename,sizeof(filename),"%s/trajectory-positions-%03d-%03d.trr",sim_dir,sim_idx,tpf_idx);
    file_o1 = std::fopen(filename,"wb");
    if( file_o1==nullptr){ std::fprintf(stderr,"Error opening the trajectory positions file.\n"); return EXIT_FAILURE;}
  }

  //Simulation

  float sig = 1.0;

  for( int i_f = 0; i_f<sp.F; ++i_f)
  {
    std::fprintf(logfile,"Progress:%05.1lf%%",(100.0*i_f)/(1.0*sp.F));
    std::fseek(logfile,-15,SEEK_CUR);

    for( int i_s = 0; i_s<n_s; ++i_s)
    {
      call_PRNG<<<n_blocks,threads_block>>>(c_rn,nrn,state);

      calc_extern_f<<<n_blocks,threads_block>>>(sp.N,f_c,f_2);
      calc_sphere_f<<<n_blocks,threads_block>>>(sp.N,sp.R,sig,r_2,f_2);
      calc_bonds<<<n_blocks,threads_block>>>(sp.N,r_2,b,invlen);
      calc_cosines<<<n_blocks,threads_block>>>(sp.N,b,invlen,cosine);
      calc_intern_f<<<n_blocks,threads_block>>>(sp.N,b,invlen,cosine,f_2);
      // calc_exclvol_f<<<n_blocks,threads_block>>>(sp.N,sig,r_2,f_2);

      RK_stage_1<<<n_blocks,threads_block>>>(sp.N,r_1,r_2,f_2,nrn);

      calc_extern_f<<<n_blocks,threads_block>>>(sp.N,f_c,f_1);
      calc_sphere_f<<<n_blocks,threads_block>>>(sp.N,sp.R,sig,r_1,f_1);
      calc_bonds<<<n_blocks,threads_block>>>(sp.N,r_1,b,invlen);
      calc_cosines<<<n_blocks,threads_block>>>(sp.N,b,invlen,cosine);
      calc_intern_f<<<n_blocks,threads_block>>>(sp.N,b,invlen,cosine,f_1);
      // calc_exclvol_f<<<n_blocks,threads_block>>>(sp.N,sig,r_1,f_1);

      RK_stage_2<<<n_blocks,threads_block>>>(sp.N,r_2,f_1,f_2,nrn);
    }

    cuda_check( cudaDeviceSynchronize());

    t += n_s*dt;

    write_trajectory_positions(sp.N,r_2,t,i_f,file_o1);
  }

  std::fclose(file_o1);

  print_time(logfile); std::fprintf(logfile,"Simulation finished.\n"); std::fflush(logfile);

  std::snprintf(filename,sizeof(filename),"%s/simulation-checkpoint-%03d.bin",sim_dir,sim_idx);
  file_o1 = std::fopen(filename,"wb");
  if( file_o1==nullptr){ std::fprintf(stderr,"Error opening the simulation checkpoint file.\n"); return EXIT_FAILURE;}
  save_checkpoint(sp.N,r_2,&t,n_threads,state,&tpf_idx,file_o1);
  std::fclose(file_o1);

  print_time(logfile); std::fprintf(logfile,"Simulation checkpoint file saved.\n"); std::fflush(logfile);

  std::fclose(logfile);

  //Memory deallocation

  cudaFree(r_2);
  cudaFree(r_1);

  cudaFree(f_2);
  cudaFree(f_1);

  cudaFree(nrn);
  cudaFree(state);

  cudaFree(b);
  cudaFree(invlen);
  cudaFree(cosine);

  cudaFree(f_c);

  return EXIT_SUCCESS;
}

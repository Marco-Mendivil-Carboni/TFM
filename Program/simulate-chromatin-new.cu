//Libraries

#include <cstdio> //standard input and output library
#include <cstdlib> //standard general utilities library
#include <cmath> //mathematical functions library
#include <ctime> //time utilities library
#include <stdexcept> //standard exceptions classes
#include <glob.h> //pathname pattern matching types

#include <curand_kernel.h> //cuRAND device functions

//Namespaces

namespace mmcc //Marco Mendívil Carboni code
{

//Namespaces

namespace util //utilities
{

//Classes

class error : public std::runtime_error
{
  public:
    error(const std::string &msg) : std::runtime_error{msg} {}
};

//Functions

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

} //namespace util

//Non-member Functions

inline void cuda_check( cudaError_t result)
{
  if (result!=cudaSuccess)
  {
    char msg[512];
    std::snprintf(msg,sizeof(msg),"CUDA: %s",cudaGetErrorString(result));
    throw util::error(msg);
  }
}

//Classes

class simulation
{
  public:

  //Constructor and Destructor

  simulation(FILE *f_ptr_par);

  ~simulation();

  //Host Functions

  void generate_initial_configuration();

  void write_initial_configuration(FILE *f_ptr);

  private:

  //Constants

  static constexpr float xi  = 1.000000; //damping coefficient
  static constexpr float k_B = 0.001120; //Boltzmann constant
  static constexpr float l_0 = 1.000000; //bond natural length
  static constexpr float k_e = 100.0000; //elastic constant
  static constexpr float k_b = 2.000000; //bending constant
  static constexpr float r_c = 1.122462; //LJ cutoff radius
  static constexpr float dt  = 1.0/2048; //timestep

  static constexpr int n_s = 1*2048; //MD steps between frames

  //Aliases

  using PRNGstate = curandStatePhilox4_32_10; //PRNG state type

  //Structures and Variables

  struct vecarr { //3D vector array structure
    float *x; //1st component
    float *y; //2nd component
    float *z; //3rd component
  };

  struct { //adjustable parameters
    float T; //temperature
    int N; //number of particles
    float R; //radius of sphere
    int F; //frames per file
  } ap;

  vecarr r_2; //position vector array 2
  vecarr r_1; //position vector array 1

  vecarr f_2; //force vector array 2
  vecarr f_1; //force vector array 1

//This structs are probably not neccessary if I declare the necessary variables local to the calc_bonded_f kernel

  // struct { //bond vector arrays
  //   vecarr p2; //two prior
  //   vecarr p1; //one prior
  //   vecarr n1; //one next
  //   vecarr n2; //two next
  // } b;

  // struct { //bond inverse length arrays
  //   float *p2; //two prior
  //   float *p1; //one prior
  //   float *n1; //one next
  //   float *n2; //two next
  // } invlen;

  // struct { //bond angle cosine arrays
  //   float *p; //prior
  //   float *c; //current
  //   float *n; //next
  // } cosine;

//---------------------------------------------------------------------------------------------------------------

  float c_rn; //random number constant
  vecarr nrn; //random number vector array
  PRNGstate *state; //PRNG state array

  float sig; //LJ particle size

  size_t threads_block; //threads per block
  size_t n_blocks; //number of blocks
  size_t n_threads; //nuumber of threads

  //Host Functions

  void read_parameters(FILE *f_ptr_par);

  void allocate_vecarr(vecarr &vecarr_ref, int n_elements)
  {
    cuda_check(cudaMallocManaged(&(vecarr_ref.x),n_elements*sizeof(float)));
    cuda_check(cudaMallocManaged(&(vecarr_ref.y),n_elements*sizeof(float)));
    cuda_check(cudaMallocManaged(&(vecarr_ref.z),n_elements*sizeof(float)));
  }

  void free_vecarr(vecarr &vecarr_ref)
  {
    cudaFree(vecarr_ref.x);
    cudaFree(vecarr_ref.y);
    cudaFree(vecarr_ref.z);
  }

  //Device Functions

  //Global Functions

};

//Member Functions

//Constructor and Destructor

simulation::simulation(FILE *f_ptr_par)
{
  read_parameters(f_ptr_par);

  threads_block = 256;
  n_blocks = (ap.N+threads_block-1)/threads_block;
  n_threads = n_blocks*threads_block;

  allocate_vecarr(r_2,ap.N);
  allocate_vecarr(r_1,ap.N);

  allocate_vecarr(f_2,ap.N);
  allocate_vecarr(f_1,ap.N);

  // allocate_vecarr(b.p2,ap.N);
  // allocate_vecarr(b.p1,ap.N);
  // allocate_vecarr(b.n1,ap.N);
  // allocate_vecarr(b.n2,ap.N);

  // cuda_check(cudaMallocManaged(&invlen.p2,ap.N*sizeof(float)));
  // cuda_check(cudaMallocManaged(&invlen.p1,ap.N*sizeof(float)));
  // cuda_check(cudaMallocManaged(&invlen.n1,ap.N*sizeof(float)));
  // cuda_check(cudaMallocManaged(&invlen.n2,ap.N*sizeof(float)));

  // cuda_check(cudaMallocManaged(&cosine.p,ap.N*sizeof(float)));
  // cuda_check(cudaMallocManaged(&cosine.c,ap.N*sizeof(float)));
  // cuda_check(cudaMallocManaged(&cosine.n,ap.N*sizeof(float)));

  allocate_vecarr(nrn,n_threads);

  cuda_check(cudaMallocManaged(&state,n_threads*sizeof(PRNGstate)));
}

simulation::~simulation()
{
  free_vecarr(r_2);
  free_vecarr(r_1);

  free_vecarr(f_2);
  free_vecarr(f_1);

  // free_vecarr(b.p2);
  // free_vecarr(b.p1);
  // free_vecarr(b.n1);
  // free_vecarr(b.n2);

  // cudaFree(invlen.p2);
  // cudaFree(invlen.p1);
  // cudaFree(invlen.n1);
  // cudaFree(invlen.n2);

  // cudaFree(cosine.p);
  // cudaFree(cosine.c);
  // cudaFree(cosine.n);

  free_vecarr(nrn);

  cudaFree(state);
}

//Host Functions

void simulation::read_parameters(FILE *f_ptr_par)
{
  if (std::fscanf(f_ptr_par,"T\t%f\n",&(ap.T))!=1
    ||std::fscanf(f_ptr_par,"N\t%d\n",&(ap.N))!=1
    ||std::fscanf(f_ptr_par,"R\t%f\n",&(ap.R))!=1
    ||std::fscanf(f_ptr_par,"F\t%d\n",&(ap.F))!=1)
  {
    throw util::error("unable to read parameters");
  }
  if ((ap.T)<__FLT_MIN__){ throw util::error("T must be positive");}
  if ((ap.N)<__FLT_MIN__){ throw util::error("N must be positive");}
  if ((ap.R)<__FLT_MIN__){ throw util::error("R must be positive");}
  if ((ap.F)<__FLT_MIN__){ throw util::error("F must be positive");}
}

void simulation::generate_initial_configuration()
{
  //initialize PRNG
  curandGenerator_t gen;
  curandCreateGeneratorHost(&gen,CURAND_RNG_PSEUDO_DEFAULT);
  curandSetPseudoRandomGeneratorSeed(gen,time(nullptr));

  //declare auxiliary variables
  float beta = 1.0/(k_B*ap.T);
  float rand, theta, varphi, bondlen, bondangle, perdirlen;
  float3 olddir, newdir, perdir;

  //place first particle
  curandGenerateUniform(gen,&rand,1); theta = acos(1.0-2.0*rand);
  curandGenerateUniform(gen,&rand,1); varphi = 2.0*M_PI*rand;
  olddir.x = sin(theta)*cos(varphi);
  olddir.y = sin(theta)*sin(varphi);
  olddir.z = cos(theta);
  r_2.x[0] = r_2.y[0] = r_2.z[0] = 0.0;

  //perform random walk
  int n_failures = 0;
  for (int i_p = 1; i_p<ap.N; ++i_p)
  {
    //generate random direction perpendicular to old direction
    curandGenerateUniform(gen,&rand,1); theta = acos(1.0-2.0*rand);
    curandGenerateUniform(gen,&rand,1); varphi = 2.0*M_PI*rand;
    perdir.x = olddir.y*cos(theta)-olddir.z*sin(theta)*sin(varphi);
    perdir.y = olddir.z*sin(theta)*cos(varphi)-olddir.x*cos(theta);
    perdir.z = olddir.x*sin(theta)*sin(varphi)-olddir.y*sin(theta)*cos(varphi);
    perdirlen = sqrt(perdir.x*perdir.x+perdir.y*perdir.y+perdir.z*perdir.z);
    perdir.x /= perdirlen; perdir.y /= perdirlen; perdir.z /= perdirlen;

    //generate random bond angle and calculate new direction
    curandGenerateUniform(gen,&rand,1);
    bondangle = acos(1.0+log(1.0-(1.0-exp(-2.0*(k_b*beta)))*rand)/(k_b*beta));
    newdir.x = olddir.x*cos(bondangle)+perdir.x*sin(bondangle);
    newdir.y = olddir.y*cos(bondangle)+perdir.y*sin(bondangle);
    newdir.z = olddir.z*cos(bondangle)+perdir.z*sin(bondangle);

    //calculate position of next particle
    curandGenerateUniform(gen,&rand,1);
    bondlen = l_0+sqrt(2.0/(k_e*beta))*erfinv(2.0*rand-1.0);
    r_2.x[i_p] = bondlen*newdir.x+r_2.x[i_p-1];
    r_2.y[i_p] = bondlen*newdir.y+r_2.y[i_p-1];
    r_2.z[i_p] = bondlen*newdir.z+r_2.z[i_p-1];

    //check if new position is acceptable
    int accept = 1;
    if (!isfinite(r_2.x[i_p])){ accept = 0;}
    if (!isfinite(r_2.y[i_p])){ accept = 0;}
    if (!isfinite(r_2.z[i_p])){ accept = 0;}
    for (int j_p = 0; j_p<(i_p-1); ++j_p)
    {
      float dist = 0.0;
      dist += (r_2.x[j_p]-r_2.x[i_p])*(r_2.x[j_p]-r_2.x[i_p]);
      dist += (r_2.y[j_p]-r_2.y[i_p])*(r_2.y[j_p]-r_2.y[i_p]);
      dist += (r_2.z[j_p]-r_2.z[i_p])*(r_2.z[j_p]-r_2.z[i_p]);
      dist = sqrt(dist);
      if (dist<(r_c*sig)){ accept = 0; break;}
    }
    float d_r = 0.0;
    d_r += r_2.x[i_p]*r_2.x[i_p];
    d_r += r_2.y[i_p]*r_2.y[i_p];
    d_r += r_2.z[i_p]*r_2.z[i_p];
    d_r = sqrt(d_r);
    if ((ap.R-d_r)<(r_c*sig)){ accept = 0;}

    //continue if it is accepted
    if (accept)
    {
      olddir.x = newdir.x;
      olddir.y = newdir.y;
      olddir.z = newdir.z;
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

void simulation::write_initial_configuration(FILE *f_ptr)
{
  std::fprintf(f_ptr,"Chromatin simulation, t=0.0\n");
  std::fprintf(f_ptr,"%5d\n",ap.N);
  for( int i_p = 0; i_p<ap.N; ++i_p)
  {
    std::fprintf(f_ptr,"%5d%-5s%5s%5d",i_p+1,"X","X",i_p+1);
    std::fprintf(f_ptr,"%8.3f%8.3f%8.3f\n",r_2.x[i_p],r_2.y[i_p],r_2.z[i_p]);
  }
  std::fprintf(f_ptr,"%10.5f%10.5f%10.5f\n",0.0,0.0,0.0);
}

} //namespace mmcc

//Main

//Host & Device functions can use pass by reference
//Global functions cannot use pass by reference

int main(int argc, const char **argv)
{
  //check command-line arguments
  if (argc<2){ std::fprintf(stderr,"No arguments.\n"); return EXIT_FAILURE;}
  if (argc>3){ std::fprintf(stderr,"Extra arguments.\n"); return EXIT_FAILURE;}
  if (sizeof(argv[1])>128)
  {
    std::fprintf(stderr,"Directory path too long.\n"); return EXIT_FAILURE;
  }

  //declare auxiliary variables
  char sim_dir[128]; //simulation directory
  std::snprintf(sim_dir,sizeof(sim_dir),"%s",argv[1]);
  FILE *f_ptr_par; //pointer to parameter file
  FILE *f_ptr_out; //pointer to output file
  FILE *f_ptr_log; //pointer to log file
  char f_str[256]; //file name string

  try
  {
    //open log file
    std::snprintf(f_str,sizeof(f_str),"%s/current-progress.log",sim_dir);
    f_ptr_log = mmcc::util::fopen(f_str,"wt");

    //read parameters and initialize simulation
    std::snprintf(f_str,sizeof(f_str),"%s/adjustable-parameters.dat",sim_dir);
    f_ptr_par = mmcc::util::fopen(f_str,"rt");
    mmcc::simulation sim(f_ptr_par);
    std::fclose(f_ptr_par);
    // setup_PRNG<<<n_blocks,threads_block>>>(time(nullptr),state);
    // cuda_check( cudaDeviceSynchronize());

    // begin new simulation or continue a previous one
    int sim_idx = 0; //simulation index
    // int tpf_idx = 0; //trajectory positions file index
    // float t = 0.0; //simulation time
    if (argc==2)
    {
      glob_t prev_sims;
      std::snprintf(f_str,sizeof(f_str),"%s/initial-configuration-*",sim_dir);
      if (glob(f_str,0,nullptr,&prev_sims)==0)
      {
        sim_idx = prev_sims.gl_pathc;
      }
      globfree(&prev_sims);
      // print_time(logfile); std::fprintf(logfile,"New simulation started.\n");
      // std::fprintf(logfile,"sim_idx=%03d tpf_idx=%03d\n",sim_idx,tpf_idx); std::fflush(logfile);
      sim.generate_initial_configuration();

      std::snprintf(f_str,sizeof(f_str),"%s/initial-configuration-%03d.gro",sim_dir,sim_idx);
      f_ptr_out = mmcc::util::fopen(f_str,"wt");
      sim.write_initial_configuration(f_ptr_out);
      std::fclose(f_ptr_out);
    }
    else
    {

    }

    //close log file
    std::fclose(f_ptr_log);
  }
  catch (const mmcc::util::error& error)
  {
    //do something
  }

  //exit program successfully
  return EXIT_SUCCESS;
}

//Includes

#include "chrdat.cuh" //chromatin data

#include </usr/local/cuda/samples/common/inc/helper_math.h> //float4 utilities

//Namespace

namespace mmc //Marco Mendívil Carboni
{

//Functions

//chrdat constructor
chrdat::chrdat(parmap &par) //parameters
  : N {par.get_val<int>("number_of_particles",0)}
  , R {par.get_val<float>("confinement_radius",-1.0)}
  , T {par.get_val<float>("temperature",298.0)}
  , i_f {0}, t {0.0}
  , sig {1.0}
{
  //check parameters
  if (N<1){ throw error("number_of_particles out of range");}
  if (R<0.0){ throw error("confinement_radius out of range");}
  if (T<0.0){ throw error("temperature out of range");}
  float cvf = N*pow(0.5*sig/(R-0.5*sig),3); //chromatin volume fraction
  if (cvf>0.5){ throw error("chromatin volume fraction above 0.5");}
  std::string msg = "parameters:"; //message
  msg += " N = "+cnfs(N,5,'0');
  msg += " R = "+cnfs(R,6,'0',2);
  msg += " T = "+cnfs(T,6,'0',2);
  logger::record(msg);

  //allocate host memory
  pt = new int[N];
  r = new float4[N];
  f = new float4[N];
}

//chrdat destructor
chrdat::~chrdat()
{
  delete[] pt;
  delete[] r;
  delete[] f;
}

} //namespace mmc

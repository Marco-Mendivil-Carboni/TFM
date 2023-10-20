//Includes

#include "chrana.cuh" //chromatin analysis

//Namespace

namespace mmc //Marco Mendívil Carboni
{

//Functions

//chromatin analysis constructor
chrana::chrana(parmap &par) //parameters
  : chrdat(par)
  , fpf {par.get_val<uint>("frames_per_file",100)}
{

}

//chromatin analysis destructor
chrana::~chrana()
{
  for (const float &f : rg2_s){ printf("%f\n",f);}
}

//add initial condition to analysis
void chrana::add_initial_condition(std::ifstream &txt_inp_f) //text input file
{
  //read initial condition frame
  read_frame_txt(txt_inp_f);
  cmr_s.push_back(center_of_mass());
  rg2_s.push_back(gyration_radius_sq(cmr_s.back()));
}

//add trajectory file to analysis
void chrana::add_trajectory_file(std::ifstream &bin_inp_f) //binary input file
{
  //iterate over all frames per file
  for (uint ffi = 0; ffi<fpf; ++ffi) //file frame index
  {
    //read trajectory frame
    read_frame_bin(bin_inp_f);
    cmr_s.push_back(center_of_mass());
    rg2_s.push_back(gyration_radius_sq(cmr_s.back()));
  }
}

//center of mass position
vec3f chrana::center_of_mass()
{
  vec3f cmr = {0.0,0.0,0.0}; //center of mass position
  for (uint i_p = 0; i_p<N; ++i_p)
  {
    cmr += hr[i_p];
  }
  cmr /= N;
  return cmr;
}

//gyration radius squared
float chrana::gyration_radius_sq(vec3f cmr) //center of mass position
{
  float rg2 = 0.0; //gyration radius squared
  for (uint i_p = 0; i_p<N; ++i_p)
  {
    rg2 += dot(hr[i_p]-cmr,hr[i_p]-cmr);
  }
  rg2 /= N;
  return rg2;
}

} //namespace mmc

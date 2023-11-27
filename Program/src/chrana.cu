//Includes

#include "chrana.cuh" //chromatin analysis

//Namespace

namespace mmc //Marco Mendívil Carboni
{

//Global Functions

//calculate the mean spatial distance
__global__ void calc_msd(
  const uint lma, //length of msd arrays
  const uint N, //number of particles
  vec3f *r, //position array
  float *ma) //msd array
{
  //calculate array index
  int i_a = blockIdx.x*blockDim.x+threadIdx.x; //array index
  if (i_a>=lma){ return;}

  //calculate the mean spatial distance
  ma[i_a] = 0.0;
  uint s = i_a+1; //separation
  for (uint i_p = 0; i_p<(N-s); ++i_p) //particle index
  {
    ma[i_a] += length(r[i_p+s]-r[i_p]);
  }
  ma[i_a] /= (N-s);
}

//Host Functions

//chromatin analysis constructor
chrana::chrana(parmap &par) //parameters
  : chrdat(par)
  , lma {(N/8)-1}
{
  //allocate memory
  msd_o = new simobs[lma];

  //allocate device memory
  cuda_check(cudaMalloc(&ma,lma*sizeof(float)));

  //allocate host memory
  cuda_check(cudaMallocHost(&hma,lma*sizeof(float)));
}

//chromatin analysis destructor
chrana::~chrana()
{
  //deallocate memory
  delete[] msd_o;

  //deallocate device memory
  cuda_check(cudaFree(ma));

  //deallocate host memory
  cuda_check(cudaFreeHost(hma));
}

//add initial condition to analysis
void chrana::add_initial_condition(std::ifstream &txt_inp_f) //text input file
{
  //read initial condition frame
  read_frame_txt(txt_inp_f);

  //calculate observables
  calc_observables();
}

//add trajectory file to analysis
void chrana::add_trajectory_file(std::ifstream &bin_inp_f) //binary input file
{
  //iterate over all frames per file
  for (uint ffi = 0; ffi<fpf; ++ffi) //file frame index
  {
    //read trajectory frame
    read_frame_bin(bin_inp_f);

    //calculate observables
    calc_observables();
  }
}

//calculate individual simulation statistics
void chrana::calc_ind_sim_stat()
{
  //calculate dcm statistics
  tdstat dcm_s; //dcm statistics //make a function for this ---------------------
  calc_stats(dcm_o.is_ts,dcm_s);
  dcm_o.is_sv.push_back(dcm_s);

  //calculate rg2 statistics
  tdstat rg2_s; //rg2 statistics
  calc_stats(rg2_o.is_ts,rg2_s);
  rg2_o.is_sv.push_back(rg2_s);

  //calculate nop statistics
  tdstat nop_s; //nop statistics
  calc_stats(nop_o.is_ts,nop_s);
  nop_o.is_sv.push_back(nop_s);

  //calculate rcd statistics
  tdstat rcd_s[3][n_b]; //rcd statistics
  for (uint i_t = 0; i_t<3; ++i_t) //type index
  {
    for (uint i_b = 0; i_b<n_b; ++i_b) //bin index
    {
      calc_stats(rcd_o[i_t][i_b].is_ts,rcd_s[i_t][i_b]);
      rcd_o[i_t][i_b].is_sv.push_back(rcd_s[i_t][i_b]);
    }
  }

  //calculate msd statistics
  tdstat *msd_s = new tdstat[lma]; //msd statistics
  for (uint i_a = 0; i_a<lma; ++i_a) //array index
  {
    calc_stats(msd_o[i_a].is_ts,msd_s[i_a]);
    msd_o[i_a].is_sv.push_back(msd_s[i_a]);
  }
  delete[] msd_s;
}

//save individual simulation analysis results
void chrana::save_ind_sim_results(std::ofstream &txt_out_f) //text output file
{
  //save dcm, rg2 and nop statistics
  txt_out_f<<"#individual simulation analysis\n";
  txt_out_f<<"#        avg   sqrt(var)         sem   f_n_b     i_t ter\n";
  txt_out_f<<"# center of mass distance:\n";
  tdstat dcm_s = dcm_o.is_sv.back(); //dcm statistics //make function for this -------------
  txt_out_f<<cnfs(dcm_s.avg,12,' ',6)<<cnfs(sqrt(dcm_s.var),12,' ',6);
  txt_out_f<<cnfs(dcm_s.sem,12,' ',6)<<cnfs(dcm_s.f_n_b,8,' ');
  txt_out_f<<cnfs(dcm_s.i_t,8,' ')<<(dcm_s.ter?" yes":"  no")<<"\n";
  txt_out_f<<"# gyration radius squared:\n";
  tdstat rg2_s = rg2_o.is_sv.back(); //rg2 statistics
  txt_out_f<<cnfs(rg2_s.avg,12,' ',6)<<cnfs(sqrt(rg2_s.var),12,' ',6);
  txt_out_f<<cnfs(rg2_s.sem,12,' ',6)<<cnfs(rg2_s.f_n_b,8,' ');
  txt_out_f<<cnfs(rg2_s.i_t,8,' ')<<(rg2_s.ter?" yes":"  no")<<"\n";
  txt_out_f<<"# nematic order parameter:\n";
  tdstat nop_s = nop_o.is_sv.back(); //nop statistics
  txt_out_f<<cnfs(nop_s.avg,12,' ',6)<<cnfs(sqrt(nop_s.var),12,' ',6);
  txt_out_f<<cnfs(nop_s.sem,12,' ',6)<<cnfs(nop_s.f_n_b,8,' ');
  txt_out_f<<cnfs(nop_s.i_t,8,' ')<<(nop_s.ter?" yes":"  no")<<"\n";
  txt_out_f<<"\n\n";

  //save rcd statistics
  tdstat rcd_s[3][n_b]; //rcd statistics
  for (uint i_t = 0; i_t<3; ++i_t) //type index
  {
    txt_out_f<<"#        r_b         avg   sqrt(var)         sem\n";
    txt_out_f<<"    0.000000    0.000000    0.000000    0.000000\n";
    for (uint i_b = 0; i_b<n_b; ++i_b) //bin index
    {
      rcd_s[i_t][i_b] = rcd_o[i_t][i_b].is_sv.back();
      txt_out_f<<cnfs(ng.d_m*pow((i_b+1.0)/n_b,1.0/3),12,' ',6);
      txt_out_f<<cnfs(rcd_s[i_t][i_b].avg,12,' ',9); //make function for this ----------
      txt_out_f<<cnfs(sqrt(rcd_s[i_t][i_b].var),12,' ',9);
      txt_out_f<<cnfs(rcd_s[i_t][i_b].sem,12,' ',9);
      txt_out_f<<"\n";
    }
    txt_out_f<<"\n\n";
  }

  //save msd statistics
  tdstat *msd_s = new tdstat[lma]; //msd statistics //this will dissappear after using a function to print
  txt_out_f<<"#    s         avg   sqrt(var)         sem\n";
  for (uint i_a = 0; i_a<lma; ++i_a) //array index
  {
    msd_s[i_a] = msd_o[i_a].is_sv.back();
    txt_out_f<<cnfs((i_a+1),6,' ');
    txt_out_f<<cnfs(msd_s[i_a].avg,12,' ',6);
    txt_out_f<<cnfs(sqrt(msd_s[i_a].var),12,' ',6);
    txt_out_f<<cnfs(msd_s[i_a].sem,12,' ',6);
    txt_out_f<<"\n";
  }
  txt_out_f<<"\n\n";
  delete[] msd_s;

  //save current simulation time series of scalar observables
  uint n_e = t_o.is_ts.size(); //number of elements
  txt_out_f<<"#          t         dcm         rg2         nop\n";
  for (uint i_e = 0; i_e<n_e; ++i_e) //element index
  {
    txt_out_f<<cnfs(t_o.is_ts[i_e],12,' ',2);
    txt_out_f<<cnfs(dcm_o.is_ts[i_e],12,' ',6);
    txt_out_f<<cnfs(rg2_o.is_ts[i_e],12,' ',6);
    txt_out_f<<cnfs(nop_o.is_ts[i_e],12,' ',6);
    txt_out_f<<"\n";
  }

  //check filestream
  if (txt_out_f.fail())
  {
    throw mmc::error("failed to write results to text file");
  }
}

//clear individual simulation analysis data
void chrana::clear_ind_sim_data()
{
  //clear time vector
  t_o.is_ts.clear();

  //clear center of mass distance vector
  dcm_o.is_ts.clear();

  //clear gyration radius squared vector
  rg2_o.is_ts.clear();

  //clear nematic order parameter vector
  nop_o.is_ts.clear();

  //clear radial chromatin density vector
  for (uint i_t = 0; i_t<3; ++i_t) //type index
  {
    for (uint i_b = 0; i_b<n_b; ++i_b) //bin index
    {
      rcd_o[i_t][i_b].is_ts.clear();
    }
  }

  //clear mean spatial distance vector
  for (uint i_a = 0; i_a<lma; ++i_a) //array index
  {
    msd_o[i_a].is_ts.clear();
  }
}

//calculate final statistics
void chrana::calc_fin_stat()
{
  //calculate dcm final statistics
  calc_stats(dcm_o.is_sv,dcm_o.cs_fs);

  //calculate rg2 final statistics
  calc_stats(rg2_o.is_sv,rg2_o.cs_fs);

  //calculate nop final statistics
  calc_stats(nop_o.is_sv,nop_o.cs_fs);

  //calculate rcd final statistics
  for (uint i_t = 0; i_t<3; ++i_t) //type index
  {
    for (uint i_b = 0; i_b<n_b; ++i_b) //bin index
    {
      calc_stats(rcd_o[i_t][i_b].is_sv,rcd_o[i_t][i_b].cs_fs);
    }
  }

  //calculate msd final statistics
  for (uint i_a = 0; i_a<lma; ++i_a) //array index
  {
    calc_stats(msd_o[i_a].is_sv,msd_o[i_a].cs_fs);
  }
}

//save final analysis results
void chrana::save_fin_results(std::ofstream &txt_out_f) //text output file
{
  //save dcm, rg2 and nop final statistics
  txt_out_f<<"#combined simulations final analysis\n";
  txt_out_f<<"#        avg   sqrt(var)         sem\n";
  txt_out_f<<"# center of mass distance:\n";
  txt_out_f<<cnfs(dcm_o.cs_fs.avg,12,' ',6)<<cnfs(sqrt(dcm_o.cs_fs.var),12,' ',6);
  txt_out_f<<cnfs(dcm_o.cs_fs.sem,12,' ',6)<<"\n";
  txt_out_f<<"# gyration radius squared:\n";
  txt_out_f<<cnfs(rg2_o.cs_fs.avg,12,' ',6)<<cnfs(sqrt(rg2_o.cs_fs.var),12,' ',6);
  txt_out_f<<cnfs(rg2_o.cs_fs.sem,12,' ',6)<<"\n";
  txt_out_f<<"# nematic order parameter:\n";
  txt_out_f<<cnfs(nop_o.cs_fs.avg,12,' ',6)<<cnfs(sqrt(nop_o.cs_fs.var),12,' ',6);
  txt_out_f<<cnfs(nop_o.cs_fs.sem,12,' ',6)<<"\n";
  txt_out_f<<"\n\n";

  //save rcd final statistics
  for (uint i_t = 0; i_t<3; ++i_t) //type index
  {
    txt_out_f<<"#        r_b         avg   sqrt(var)         sem\n";
    txt_out_f<<"    0.000000    0.000000    0.000000    0.000000\n";
    for (uint i_b = 0; i_b<n_b; ++i_b) //bin index
    {
      txt_out_f<<cnfs(ng.d_m*pow((i_b+1.0)/n_b,1.0/3),12,' ',6);
      txt_out_f<<cnfs(rcd_o[i_t][i_b].cs_fs.avg,12,' ',9);
      txt_out_f<<cnfs(sqrt(rcd_o[i_t][i_b].cs_fs.var),12,' ',9);
      txt_out_f<<cnfs(rcd_o[i_t][i_b].cs_fs.sem,12,' ',9);
      txt_out_f<<"\n";
    }
    txt_out_f<<"\n\n";
  }

  //save msd final statistics
  txt_out_f<<"#    s         avg   sqrt(var)         sem\n";
  for (uint i_a = 0; i_a<lma; ++i_a) //array index
  {
    txt_out_f<<cnfs((i_a+1),6,' ');
    txt_out_f<<cnfs(msd_o[i_a].cs_fs.avg,12,' ',6);
    txt_out_f<<cnfs(sqrt(msd_o[i_a].cs_fs.var),12,' ',6);
    txt_out_f<<cnfs(msd_o[i_a].cs_fs.sem,12,' ',6);
    txt_out_f<<"\n";
  }
  txt_out_f<<"\n\n";

  //check filestream
  if (txt_out_f.fail())
  {
    throw mmc::error("failed to write results to text file");
  }
}

//calculate observables
void chrana::calc_observables()
{
  //store time
  t_o.is_ts.push_back(t);

  //calculate the center of mass position
  vec3f cmr = {0.0,0.0,0.0}; //center of mass position
  for (uint i_p = 0; i_p<N; ++i_p) //particle index
  {
    cmr += hr[i_p];
  }
  cmr /= N;

  //calculate the center of mass distance
  float dcm = length(cmr); //center of mass distance
  dcm_o.is_ts.push_back(dcm);

  //calculate the gyration radius squared
  float rg2 = 0.0; //gyration radius squared
  for (uint i_p = 0; i_p<N; ++i_p) //particle index
  {
    rg2 += dot(hr[i_p]-cmr,hr[i_p]-cmr);
  }
  rg2 /= N;
  rg2_o.is_ts.push_back(rg2);

  //calculate the nematic order parameter
  vec3f vec; //bond vector
  float cos; //cosine of the angle
  float nop = 0.0; //nematic order parameter
  for (uint i_p = 0; i_p<(N-1); ++i_p) //particle index
  {
    vec = hr[i_p+1]-hr[i_p];
    cos = dot(vec,hr[i_p])/(length(vec)*length(hr[i_p]));
    nop += 0.5*(3.0*cos*cos-1.0);
  }
  nop /= N-1.0;
  nop_o.is_ts.push_back(nop);

  //calculate the radial chromatin density
  float rcd[3][n_b]; //radial chromatin density
  for (uint i_t = 0; i_t<3; ++i_t) //type index
  {
    for (uint i_b = 0; i_b<n_b; ++i_b) //bin index
    {
      rcd[i_t][i_b] = 0.0;
    }
  }
  for (uint i_p = 0; i_p<N; ++i_p) //particle index
  {
    float d_r = length(hr[i_p]); //radial distance
    uint i_b = n_b*d_r*d_r*d_r/(ng.d_m*ng.d_m*ng.d_m); //bin index
    if (hpt[i_p]==LND){ rcd[0][i_b] += 1.0;}
    if (hpt[i_p]==LAD){ rcd[1][i_b] += 1.0;}
    rcd[2][i_b] += 1.0;
  }
  for (uint i_t = 0; i_t<3; ++i_t) //type index
  {
    for (uint i_b = 0; i_b<n_b; ++i_b) //bin index
    {
      rcd[i_t][i_b] /= (4.0/3.0)*M_PI*ng.d_m*ng.d_m*ng.d_m/n_b;
      rcd_o[i_t][i_b].is_ts.push_back(rcd[i_t][i_b]);
    }
  }

  //calculate the mean spatial distance
  uint tpb = 128; //threads per block
  calc_msd<<<(lma+tpb-1)/tpb,tpb>>>(lma,N,r,ma);
  cuda_check(cudaMemcpy(hma,ma,lma*sizeof(float),cudaMemcpyDeviceToHost));
  for (uint i_a = 0; i_a<lma; ++i_a) //array index
  {
    msd_o[i_a].is_ts.push_back(hma[i_a]);
  }
}

} //namespace mmc

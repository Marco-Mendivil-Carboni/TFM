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
  , fpf {par.get_val<uint>("frames_per_file",100)}
  , lma {(N/2)-1}
{
  //allocate memory
  msd_v = new std::vector<float>[lma];
  msd_s_v = new std::vector<tdstat>[lma];
  msd_f_s = new idstat[lma];

  //allocate device memory
  cuda_check(cudaMalloc(&ma,lma*sizeof(float)));

  //allocate host memory
  cuda_check(cudaMallocHost(&hma,lma*sizeof(float)));
}

//chromatin analysis destructor
chrana::~chrana()
{
  //deallocate memory
  delete[] msd_v;
  delete[] msd_s_v;
  delete[] msd_f_s;

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
  tdstat dcm_s; //dcm statistics
  calc_stats(dcm_v,dcm_s);
  dcm_s_v.push_back(dcm_s);

  //calculate rg2 statistics
  tdstat rg2_s; //rg2 statistics
  calc_stats(rg2_v,rg2_s);
  rg2_s_v.push_back(rg2_s);

  //calculate nop statistics
  tdstat nop_s; //nop statistics
  calc_stats(nop_v,nop_s);
  nop_s_v.push_back(nop_s);

  //calculate rcd statistics
  tdstat rcd_s[3][n_b]; //rcd statistics
  for (uint i_t = 0; i_t<3; ++i_t) //type index
  {
    for (uint i_b = 0; i_b<n_b; ++i_b) //bin index
    {
      calc_stats(rcd_v[i_t][i_b],rcd_s[i_t][i_b]);
      rcd_s_v[i_t][i_b].push_back(rcd_s[i_t][i_b]);
    }
  }

  //calculate msd statistics
  tdstat *msd_s = new tdstat[lma]; //msd statistics
  for (uint i_a = 0; i_a<lma; ++i_a) //array index
  {
    calc_stats(msd_v[i_a],msd_s[i_a]);
    msd_s_v[i_a].push_back(msd_s[i_a]);
  }
  delete[] msd_s;
}

//save individual simulation analysis results
void chrana::save_ind_sim_results(std::ofstream &txt_out_f) //text output file
{
  //save dcm, rg2 and nop statistics
  txt_out_f<<"#        avg   sqrt(var)         sem   f_n_b    t_v[i_t] ter\n";
  txt_out_f<<"# center of mass distance:\n";
  tdstat dcm_s = dcm_s_v.back(); //dcm statistics
  txt_out_f<<cnfs(dcm_s.avg,12,' ',6)<<cnfs(sqrt(dcm_s.var),12,' ',6);
  txt_out_f<<cnfs(dcm_s.sem,12,' ',6)<<cnfs(dcm_s.f_n_b,8,' ');
  txt_out_f<<cnfs(t_v[dcm_s.i_t],12,' ',2)<<(dcm_s.ter?" yes":"  no")<<"\n";
  txt_out_f<<"# gyration radius squared:\n";
  tdstat rg2_s = rg2_s_v.back(); //rg2 statistics
  txt_out_f<<cnfs(rg2_s.avg,12,' ',6)<<cnfs(sqrt(rg2_s.var),12,' ',6);
  txt_out_f<<cnfs(rg2_s.sem,12,' ',6)<<cnfs(rg2_s.f_n_b,8,' ');
  txt_out_f<<cnfs(t_v[rg2_s.i_t],12,' ',2)<<(rg2_s.ter?" yes":"  no")<<"\n";
  txt_out_f<<"# nematic order parameter:\n";
  tdstat nop_s = nop_s_v.back(); //nop statistics
  txt_out_f<<cnfs(nop_s.avg,12,' ',6)<<cnfs(sqrt(nop_s.var),12,' ',6);
  txt_out_f<<cnfs(nop_s.sem,12,' ',6)<<cnfs(nop_s.f_n_b,8,' ');
  txt_out_f<<cnfs(t_v[nop_s.i_t],12,' ',2)<<(nop_s.ter?" yes":"  no")<<"\n";
  txt_out_f<<"\n\n";

  //save rcd statistics
  tdstat rcd_s[3][n_b]; //rcd statistics
  for (uint i_t = 0; i_t<3; ++i_t) //type index
  {
    txt_out_f<<"#        r_b         avg   sqrt(var)         sem   f_n_b ter\n";
    txt_out_f<<"    0.000000    0.000000    0.000000    0.000000       - yes\n";
    for (uint i_b = 0; i_b<n_b; ++i_b) //bin index
    {
      rcd_s[i_t][i_b] = rcd_s_v[i_t][i_b].back();
      txt_out_f<<cnfs(R*pow((i_b+1.0)/n_b,1.0/3),12,' ',6);
      txt_out_f<<cnfs(rcd_s[i_t][i_b].avg,12,' ',9);
      txt_out_f<<cnfs(sqrt(rcd_s[i_t][i_b].var),12,' ',9);
      txt_out_f<<cnfs(rcd_s[i_t][i_b].sem,12,' ',9);
      txt_out_f<<cnfs(rcd_s[i_t][i_b].f_n_b,8,' ');
      txt_out_f<<(rcd_s[i_t][i_b].ter?" yes":"  no");
      txt_out_f<<"\n";
    }
    txt_out_f<<"\n\n";
  }

  //save msd statistics
  tdstat *msd_s = new tdstat[lma]; //msd statistics
  txt_out_f<<"#    s         avg   sqrt(var)         sem   f_n_b ter\n";
  for (uint i_a = 0; i_a<lma; ++i_a) //array index
  {
    msd_s[i_a] = msd_s_v[i_a].back();
    txt_out_f<<cnfs((i_a+1),6,' ');
    txt_out_f<<cnfs(msd_s[i_a].avg,12,' ',6);
    txt_out_f<<cnfs(sqrt(msd_s[i_a].var),12,' ',6);
    txt_out_f<<cnfs(msd_s[i_a].sem,12,' ',6);
    txt_out_f<<cnfs(msd_s[i_a].f_n_b,8,' ');
    txt_out_f<<(msd_s[i_a].ter?" yes":"  no");
    txt_out_f<<"\n";
  }
  txt_out_f<<"\n\n";
  delete[] msd_s;

  //save vectors
  uint n_e = t_v.size(); //number of elements
  for (uint i_e = 0; i_e<n_e; ++i_e) //element index
  {
    txt_out_f<<cnfs(t_v[i_e],12,' ',2);
    txt_out_f<<cnfs(dcm_v[i_e],12,' ',6);
    txt_out_f<<cnfs(rg2_v[i_e],12,' ',6);
    txt_out_f<<cnfs(nop_v[i_e],12,' ',6);
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
  t_v.clear();

  //clear center of mass distance vector
  dcm_v.clear();

  //clear gyration radius squared vector
  rg2_v.clear();

  //clear nematic order parameter vector
  nop_v.clear();

  //clear radial chromatin density vector
  for (uint i_t = 0; i_t<3; ++i_t) //type index
  {
    for (uint i_b = 0; i_b<n_b; ++i_b) //bin index
    {
      rcd_v[i_t][i_b].clear();
    }
  }

  //clear mean spatial distance vector
  for (uint i_a = 0; i_a<lma; ++i_a) //array index
  {
    msd_v[i_a].clear();
  }
}

//calculate final statistics
void chrana::calc_fin_stat()
{
  //calculate dcm final statistics
  calc_stats(dcm_s_v,dcm_f_s);

  //calculate rg2 final statistics
  calc_stats(rg2_s_v,rg2_f_s);

  //calculate nop final statistics
  calc_stats(nop_s_v,nop_f_s);

  //calculate rcd final statistics
  for (uint i_t = 0; i_t<3; ++i_t) //type index
  {
    for (uint i_b = 0; i_b<n_b; ++i_b) //bin index
    {
      calc_stats(rcd_s_v[i_t][i_b],rcd_f_s[i_t][i_b]);
    }
  }

  //calculate msd final statistics
  for (uint i_a = 0; i_a<lma; ++i_a) //array index
  {
    calc_stats(msd_s_v[i_a],msd_f_s[i_a]);
  }
}

//save final analysis results
void chrana::save_fin_results(std::ofstream &txt_out_f) //text output file
{
  //save dcm, rg2 and nop final statistics
  txt_out_f<<"#final analysis\n";
  txt_out_f<<"#        avg   sqrt(var)         sem\n";
  txt_out_f<<"# center of mass distance:\n";
  txt_out_f<<cnfs(dcm_f_s.avg,12,' ',6)<<cnfs(sqrt(dcm_f_s.var),12,' ',6);
  txt_out_f<<cnfs(dcm_f_s.sem,12,' ',6)<<"\n";
  txt_out_f<<"# gyration radius squared:\n";
  txt_out_f<<cnfs(rg2_f_s.avg,12,' ',6)<<cnfs(sqrt(rg2_f_s.var),12,' ',6);
  txt_out_f<<cnfs(rg2_f_s.sem,12,' ',6)<<"\n";
  txt_out_f<<"# nematic order parameter:\n";
  txt_out_f<<cnfs(nop_f_s.avg,12,' ',6)<<cnfs(sqrt(nop_f_s.var),12,' ',6);
  txt_out_f<<cnfs(nop_f_s.sem,12,' ',6)<<"\n";
  txt_out_f<<"\n\n";

  //save rcd final statistics
  for (uint i_t = 0; i_t<3; ++i_t) //type index
  {
    txt_out_f<<"#        r_b         avg   sqrt(var)         sem\n";
    txt_out_f<<"    0.000000    0.000000    0.000000    0.000000\n";
    for (uint i_b = 0; i_b<n_b; ++i_b) //bin index
    {
      txt_out_f<<cnfs(R*pow((i_b+1.0)/n_b,1.0/3),12,' ',6);
      txt_out_f<<cnfs(rcd_f_s[i_t][i_b].avg,12,' ',9);
      txt_out_f<<cnfs(sqrt(rcd_f_s[i_t][i_b].var),12,' ',9);
      txt_out_f<<cnfs(rcd_f_s[i_t][i_b].sem,12,' ',9);
      txt_out_f<<"\n";
    }
    txt_out_f<<"\n\n";
  }

  //save msd final statistics
  txt_out_f<<"#    s         avg   sqrt(var)         sem\n";
  for (uint i_a = 0; i_a<lma; ++i_a) //array index
  {
    txt_out_f<<cnfs((i_a+1),6,' ');
    txt_out_f<<cnfs(msd_f_s[i_a].avg,12,' ',6);
    txt_out_f<<cnfs(sqrt(msd_f_s[i_a].var),12,' ',6);
    txt_out_f<<cnfs(msd_f_s[i_a].sem,12,' ',6);
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
  t_v.push_back(t);

  //calculate the center of mass position
  vec3f cmr = {0.0,0.0,0.0}; //center of mass position
  for (uint i_p = 0; i_p<N; ++i_p) //particle index
  {
    cmr += hr[i_p];
  }
  cmr /= N;

  //calculate the center of mass distance
  float dcm = length(cmr); //center of mass distance
  dcm_v.push_back(dcm);

  //calculate the gyration radius squared
  float rg2 = 0.0; //gyration radius squared
  for (uint i_p = 0; i_p<N; ++i_p) //particle index
  {
    rg2 += dot(hr[i_p]-cmr,hr[i_p]-cmr);
  }
  rg2 /= N;
  rg2_v.push_back(rg2);

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
  nop_v.push_back(nop);

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
    uint i_b = n_b*d_r*d_r*d_r/(R*R*R); //bin index
    if (hpt[i_p]==LND){ rcd[0][i_b] += 1.0;}
    if (hpt[i_p]==LAD){ rcd[1][i_b] += 1.0;}
    rcd[2][i_b] += 1.0;
  }
  for (uint i_t = 0; i_t<3; ++i_t) //type index
  {
    for (uint i_b = 0; i_b<n_b; ++i_b) //bin index
    {
      rcd[i_t][i_b] /= N*(4.0/3.0)*M_PI*R*R*R/n_b;
      rcd_v[i_t][i_b].push_back(rcd[i_t][i_b]);
    }
  }

  //calculate the mean spatial distance
  uint tpb = 128; //threads per block
  calc_msd<<<(lma+tpb-1)/tpb,tpb>>>(lma,N,r,ma);
  cuda_check(cudaMemcpy(hma,ma,lma*sizeof(float),cudaMemcpyDeviceToHost));
  for (uint i_a = 0; i_a<lma; ++i_a) //array index
  {
    msd_v[i_a].push_back(hma[i_a]);
  }
}

//calculate statistics
void calc_stats(
  const std::vector<float> &v, //vector
  idstat &s) //statistics
{
  //calculate the first two raw moments
  double m_1 = 0.0; //1st moment
  double m_2 = 0.0; //2nd moment
  uint n_e = v.size(); //number of elements
  for (uint i_e = 0; i_e<n_e; ++i_e) //element index
  {
    m_1 += v[i_e];
    m_2 += v[i_e]*v[i_e];
  }
  m_1 /= n_e;
  m_2 /= n_e;

  //calculate statistics
  s.avg = m_1;
  s.var = (m_2-m_1*m_1)/(1.0-1.0/n_e);
  s.sem = sqrt(s.var/n_e);
}

//calculate statistics
void calc_stats(
  const std::vector<float> &v, //vector
  cdstat &s) //statistics
{
  //declare auxiliary variables
  std::vector<float> av = v; //auxiliary vector
  idstat ias; //independent auxiliary statistics

  //calculate average and variance
  calc_stats(av,ias);
  s.avg = ias.avg;
  s.var = ias.var;

  //calculate the standard error of the mean (by the blocking method)
  uint n_e = av.size(); //number of elements
  double ivm = ias.var/n_e; //independent variance of the mean
  double ulivm = ivm*(1.0+sqrt(2.0/(n_e-1.0))); //ivm upper limit
  while (n_e>=4)
  {
    //block data
    uint i_b; //block index
    for(i_b = 0; (2*i_b+1)<n_e; ++i_b)
    {
      av[i_b] = 0.5*(av[2*i_b]+av[2*i_b+1]);
    }
    n_e = i_b; av.resize(n_e);

    //calculate sem and f_n_b
    calc_stats(av,ias);
    ivm = ias.var/n_e;
    s.sem = sqrt(ivm);
    s.f_n_b = n_e;

    if (ivm>ulivm) //update the ivm upper limit
    {
      ulivm = ivm*(1.0+sqrt(2.0/(n_e-1.0)));
    }
    else //stop as the method has converged
    {
      break;
    }
  }
}

//calculate statistics
void calc_stats(
  const std::vector<float> &v, //vector
  tdstat &s) //statistics
{
  //declare auxiliary variables
  std::vector<float> av; //auxiliary vector
  idstat ias; //independent auxiliary statistics
  cdstat cas; //correlated auxiliary statistics

  //estimate termalization (by the marginal standard error rule)
  double mse; //marginal standard error
  double min_mse = INFINITY; //minimum mse
  for(uint d = 2; d<128; d*=2) //denominator
  {
    //remove termalization vector elements
    uint i_t = v.size()/d; //termalization index
    av = {v.begin()+i_t,v.end()};
    uint n_e = av.size(); //number of elements

    //calculate the marginal standard error
    calc_stats(av,ias);
    mse = ias.var*(n_e-1)/(n_e*n_e);

    //save the optimal termalization index
    if (mse<min_mse)
    {
      min_mse = mse;
      s.i_t = i_t;
    }
  }

  //determine if data has termalized
  if (s.i_t!=v.size()/2){ s.ter = true;} //termalized
  else{ s.ter = false;} //did not termalize

  //calculate the rest of the statistics
  av = {v.begin()+s.i_t,v.end()};
  calc_stats(av,cas);
  s.avg = cas.avg;
  s.var = cas.var;
  s.sem = cas.sem;
  s.f_n_b = cas.f_n_b;
}

//calculate statistics
void calc_stats(
  const std::vector<tdstat> &v, //vector
  idstat &s) //statistics
{
  //calculate the first two weighted raw moments
  double m_1 = 0.0; //1st moment
  double m_2 = 0.0; //2nd moment
  double w_1 = 0.0; //1st weight sum
  double w_2 = 0.0; //2nd weight sum
  uint n_e = v.size(); //number of elements
  for (uint i_e = 0; i_e<n_e; ++i_e) //element index
  {
    double w = 1.0/(v[i_e].sem*v[i_e].sem); //weight
    if (!isnormal(w)) //skip if weight is not normal
    {
      continue;
    }
    m_1 += w*v[i_e].avg;
    m_2 += w*v[i_e].avg*v[i_e].avg;
    w_1 += w;
    w_2 += w*w;
  }
  m_1 /= w_1;
  m_2 /= w_1;

  //calculate weighted statistics
  s.avg = m_1;
  s.var = (m_2-m_1*m_1)/(1.0-w_2/(w_1*w_1));
  s.sem = sqrt(s.var*w_2/(w_1*w_1));
}

} //namespace mmc

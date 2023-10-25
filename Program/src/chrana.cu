//Includes

#include "chrana.cuh" //chromatin analysis

//Namespace

namespace mmc //Marco Mendívil Carboni
{

//Functions

//chromatin analysis constructor
chrana::chrana(parmap &par) //parameters
  : chrdat(par)
  , fpf {par.get_val<uint>("frames_per_file",100)} {}

//chromatin analysis destructor
chrana::~chrana() {}

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

//calculate observables' statistics
void chrana::calc_observables_stat()
{
  //calculate center of mass distance statistics
  calc_stats(dcm_v,dcm_s);

  //calculate gyration radius squared statistics
  calc_stats(rg2_v,rg2_s);
}

//save analysis results
void chrana::save_results(std::ofstream &txt_out_f) //text output file
{
  //save center of mass distance statistics
  txt_out_f<<"# center of mass distance:\n";
  txt_out_f<<"#        avg   sqrt(var)         sem   f_n_b     i_t ter\n";
  txt_out_f<<cnfs(dcm_s.avg,12,' ',6)<<cnfs(sqrt(dcm_s.var),12,' ',6);
  txt_out_f<<cnfs(dcm_s.sem,12,' ',6)<<cnfs(dcm_s.f_n_b,8,' ');
  txt_out_f<<cnfs(dcm_s.i_t,8,' ')<<(dcm_s.ter?" true":" false")<<"\n";
  txt_out_f<<"\n\n";

  //save gyration radius squared statistics
  txt_out_f<<"# gyration radius squared:\n";
  txt_out_f<<"#        avg   sqrt(var)         sem   f_n_b     i_t ter\n";
  txt_out_f<<cnfs(rg2_s.avg,12,' ',6)<<cnfs(sqrt(rg2_s.var),12,' ',6);
  txt_out_f<<cnfs(rg2_s.sem,12,' ',6)<<cnfs(rg2_s.f_n_b,8,' ');
  txt_out_f<<cnfs(rg2_s.i_t,8,' ')<<(rg2_s.ter?" true":" false")<<"\n";
  txt_out_f<<"\n\n";

  //save center of mass distance vector
  for (float dcm : dcm_v) //center of mass distance
  {
    txt_out_f<<cnfs(dcm,12,' ',6)<<"\n";
  }
  txt_out_f<<"\n\n";

  //save gyration radius squared vector
  for (float rg2 : rg2_v) //gyration radius squared
  {
    txt_out_f<<cnfs(rg2,12,' ',6)<<"\n";
  }
  txt_out_f<<"\n\n";
}

//calculate observables
void chrana::calc_observables()
{
  //calculate the center of mass position
  vec3f cmr = {0.0,0.0,0.0}; //center of mass position
  for (uint i_p = 0; i_p<N; ++i_p)
  {
    cmr += hr[i_p];
  }
  cmr /= N;

  //calculate the center of mass distance
  float dcm = length(cmr); //center of mass distance
  dcm_v.push_back(dcm);

  //calculate the gyration radius squared
  float rg2 = 0.0; //gyration radius squared
  for (uint i_p = 0; i_p<N; ++i_p)
  {
    rg2 += dot(hr[i_p]-cmr,hr[i_p]-cmr);
  }
  rg2 /= N;
  rg2_v.push_back(rg2);
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
      s.i_t = i_t; min_mse = mse;
    }
  }
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
  const std::vector<float[2]> &v, //vector
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
    m_1 += v[i_e][0]/v[i_e][1];
    m_2 += v[i_e][0]*v[i_e][0]/v[i_e][1];
    w_1 += 1.0/v[i_e][1];
    w_2 += 1.0/v[i_e][1]*v[i_e][1];
  }
  m_1 /= w_1;
  m_2 /= w_1;

  //calculate weighted statistics
  s.avg = m_1;
  s.var = (m_2-m_1*m_1)/(1.0-w_2/(w_1*w_1));
  s.sem = sqrt(s.var*w_2/(w_1*w_1));
}

} //namespace mmc

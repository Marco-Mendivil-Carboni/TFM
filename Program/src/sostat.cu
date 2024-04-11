//Includes

#include "sostat.cuh" //simulation observable statistics

//Namespace

namespace mmc //Marco Mendívil Carboni
{

//Functions

//calculate last individual simulation statistics
void simobs::calc_last_is_stat()
{
  tdstat tas; //time series auxiliary statistics
  calc_stats(is_ts,tas);
  is_sv.push_back(tas);
}

//save last individual simulation statistics
void simobs::save_last_is_stat(std::ofstream &txt_out_f) //text output file
{
  tdstat tas = is_sv.back(); //time series auxiliary statistics
  txt_out_f<<cnfs(tas.avg,12,' ',6)<<cnfs(sqrt(tas.var),12,' ',6);
  txt_out_f<<cnfs(tas.sem,12,' ',6)<<cnfs(tas.f_n_b,8,' ');
  txt_out_f<<cnfs(tas.i_t,8,' ')<<(tas.thm?" yes":"  no")<<"\n";
}

//save last individual simulation statistics summary
void simobs::save_last_is_stat_s(std::ofstream &txt_out_f) //text output file
{
  tdstat tas = is_sv.back(); //time series auxiliary statistics
  txt_out_f<<cnfs(tas.avg,12,' ',6)<<cnfs(sqrt(tas.var),12,' ',6);
  txt_out_f<<cnfs(tas.sem,12,' ',6)<<"\n";
}

//calculate combined simulations final statistics
void simobs::calc_cs_final_stat()
{
  calc_stats(is_sv,cs_fs);
}

//save combined simulations final statistics
void simobs::save_cs_final_stat(std::ofstream &txt_out_f) //text output file
{
  txt_out_f<<cnfs(cs_fs.avg,12,' ',6)<<cnfs(sqrt(cs_fs.var),12,' ',6);
  txt_out_f<<cnfs(cs_fs.sem,12,' ',6)<<"\n";
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

  //estimate thermalization (by the marginal standard error rule)
  double mse; //marginal standard error
  double min_mse = INFINITY; //minimum mse
  for(uint d = 2; d<128; d*=2) //denominator
  {
    //remove thermalization vector elements
    uint i_t = v.size()/d; //thermalization index
    av = {v.begin()+i_t,v.end()};
    uint n_e = av.size(); //number of elements

    //calculate the marginal standard error
    calc_stats(av,ias);
    mse = ias.var*(n_e-1)/(n_e*n_e);

    //save the optimal thermalization index
    if (mse<min_mse)
    {
      min_mse = mse;
      s.i_t = i_t;
    }
  }

  //determine if data has thermalized
  if (s.i_t!=v.size()/2){ s.thm = true;} //thermalized
  else{ s.thm = false;} //did not thermalize

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
    if (!isfinite(w)){ w = 1.0;}
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

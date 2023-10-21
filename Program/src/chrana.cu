//Includes

#include "chrana.cuh" //chromatin analysis

//Namespace

namespace mmc //Marco Mendívil Carboni
{

//Functions

//include value
void ranvar::inc_val(float c_val) //current value
{
  val.push_back(c_val);
}

//calculate statistics
void ranvar::calc_stats()
{

}

//check if random variable has termalized
bool ranvar::has_termalized()
{
  return false;
}

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

  //calculate random variables
  calc_ranvar();
}

//add trajectory file to analysis
void chrana::add_trajectory_file(std::ifstream &bin_inp_f) //binary input file
{
  //iterate over all frames per file
  for (uint ffi = 0; ffi<fpf; ++ffi) //file frame index
  {
    //read trajectory frame
    read_frame_bin(bin_inp_f);

    //calculate random variables
    calc_ranvar();
  }
}


//calculate statistics
void chrana::calc_stats()
{

}

//save analysis results
void chrana::save_results(std::ofstream &txt_out_f) //text output file
{
  txt_out_f<<"...";
}

//calculate random variables
void chrana::calc_ranvar()
{
  //calculate the center of mass position
  vec3f cmr = {0.0,0.0,0.0}; //center of mass position
  for (uint i_p = 0; i_p<N; ++i_p)
  {
    cmr += hr[i_p];
  }
  cmr /= N;

  //calculate the distance to the center of mass
  float c_dcm = length(cmr); //current distance to the center of mass
  dcm.inc_val(c_dcm);

  //calculate the gyration radius squared
  float c_rg2 = 0.0; //current gyration radius squared
  for (uint i_p = 0; i_p<N; ++i_p)
  {
    c_rg2 += dot(hr[i_p]-cmr,hr[i_p]-cmr);
  }
  c_rg2 /= N;
  rg2.inc_val(c_rg2);
}

// void block_data(std::vector<float> &x)
// {
//   int i;
//   for( i = 0; (2*i+1)<x.size(); ++i)
//   {
//     x[i] = 0.5*(x[2*i]+x[2*i+1]);
//   }
//   x.resize(i);
// }

// void calc_moments( int n_data, double *x, double *m_1, double *m_2)
// {
//   *m_1 = 0.0;
//   *m_2 = 0.0;
//   for( int i = 0; i<n_data; i++)
//   {
//     *m_1 += x[i];
//     *m_2 += x[i]*x[i];
//   }
//   *m_1 /= n_data;
//   *m_2 /= n_data;
// }

// int calc_sig_mean( int n_data, double *x, double *sem)
// {
//   int n_data_b = n_data;
//   double m_1, m_2, var_m_1, uplim_var_m_1;
//   calc_moments(n_data_b,x,&m_1,&m_2);
//   var_m_1 = (m_2-m_1*m_1)/(n_data_b-1.0);
//   uplim_var_m_1 = var_m_1*(1.0+sqrt(2.0/(n_data_b-1.0)));
//   while( n_data_b>3)
//   {
//     block_data(&n_data_b,x);
//     calc_moments(n_data_b,x,&m_1,&m_2);
//     var_m_1 = (m_2-m_1*m_1)/(n_data_b-1.0);
//     if( var_m_1>uplim_var_m_1)
//     {
//       uplim_var_m_1 = var_m_1*(1.0+sqrt(2.0/(n_data_b-1.0)));
//     }
//     else
//     {
//       *sem = sqrt(var_m_1);
//       return n_data_b;
//     }
//   }
//   *sem = sqrt(var_m_1);
//   return n_data_b; 
// }

// int estimate_term( int n_data, double *x, int *opt_n_term)
// {
//   int i_term;
//   int n_term;
//   double m_1, m_2, smer;
//   double smer_min = INFINITY;
//   for( i_term = 0; i_term<50; i_term++)
//   {
//     n_term = i_term*n_data/100;
//     calc_moments(n_data-n_term,&x[n_term],&m_1,&m_2);
//     smer = (m_2-m_1*m_1)/(n_data-n_term);
//     if( smer<smer_min)
//     {
//       *opt_n_term = n_term;
//       smer_min = smer;
//     }
//   }
//   if( *opt_n_term==n_term)
//   {
//     return 0;
//   }
//   else
//   {
//     return 1;
//   }
// }

// void calc_statistics( int n_data, double *x, rand_var_props *x_p)
// {
//   (*x_p).m_1 = 0.0;
//   (*x_p).m_2 = 0.0;
//   for( int i = 0; i<n_data; i++)
//   {
//     (*x_p).m_1 += x[i];
//     (*x_p).m_2 += x[i]*x[i];
//   }
//   (*x_p).m_1 /= n_data;
//   (*x_p).m_2 /= n_data;
//   (*x_p).var = ((*x_p).m_2-(*x_p).m_1*(*x_p).m_1);
//   (*x_p).var *= n_data/(n_data-1.0);
//   (*x_p).sem = sqrt((*x_p).var/n_data);
// }

} //namespace mmc

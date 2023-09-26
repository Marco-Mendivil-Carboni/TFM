//Includes

#include "chrdat.cuh" //chromatin data

//Namespace

namespace mmc //Marco Mendívil Carboni
{

//Functions

//chromatin data constructor
chrdat::chrdat(parmap &par) //parameters
  : N {par.get_val<int>("number_of_particles",0)}
  , R {par.get_val<float>("confinement_radius",-1.0)}
  , T {par.get_val<float>("temperature",298.0)}
  , i_f {0}, t {0.0}
  , sig {1.0}, eps {par.get_val<float>("epsilon",0.5)}
{
  //check parameters
  if (!within(N,1,100'000)){ throw error("number_of_particles out of range");}
  if (!within(R,0.0,100.0)){ throw error("confinement_radius out of range");}
  if (!within(T,0.0,1'000.0)){ throw error("temperature out of range");}
  if (!within(eps,0.125,2.0)){ throw error("epsilon out of range");}
  float cvf = N*pow(0.5*sig/R,3); //chromatin volume fraction
  if (cvf>0.5){ throw error("chromatin volume fraction above 0.5");}
  std::string msg = "chrdat:"; //message
  msg += " N = "+cnfs(N,5,'0');
  msg += " R = "+cnfs(R,5,'0',2);
  msg += " T = "+cnfs(T,5,'0',1);
  msg += " eps = "+cnfs(eps,5,'0',3);
  logger::record(msg);

  //allocate arrays
  cuda_check(cudaMallocManaged(&pt,N*sizeof(ptype)));
  cuda_check(cudaMallocManaged(&r,N*sizeof(float4)));
  cuda_check(cudaMallocManaged(&f,N*sizeof(float4)));
}

//chromatin data destructor
chrdat::~chrdat()
{
  //deallocate arrays
  cuda_check(cudaFree(pt));
  cuda_check(cudaFree(r));
  cuda_check(cudaFree(f));
}

//write frame to text file
void chrdat::write_frame_txt(std::ofstream &txt_out_f) //text output file
{
  txt_out_f<<"Chromatin simulation, i_f = 0, t = 0.0\n";
  txt_out_f<<cnfs(N,5,' ')<<"\n";
  for (int i_p = 0; i_p<N; ++i_p) //particle index
  {
    txt_out_f<<std::setw(5)<<i_p+1<<std::left<<std::setw(5)<<"X";
    txt_out_f<<std::right<<std::setw(5)<<"X"<<std::setw(5)<<i_p+1;
    txt_out_f<<cnfs(r[i_p].x,8,' ',3);
    txt_out_f<<cnfs(r[i_p].y,8,' ',3);
    txt_out_f<<cnfs(r[i_p].z,8,' ',3);
    txt_out_f<<"\n";
  }
  txt_out_f<<cnfs(0.0,10,' ',5);
  txt_out_f<<cnfs(0.0,10,' ',5);
  txt_out_f<<cnfs(0.0,10,' ',5);
  txt_out_f<<"\n";
}

//read frame from text file
void chrdat::read_frame_txt(std::ifstream &txt_inp_f) //text input file
{

}

//write frame to binary file
void chrdat::write_frame_bin(std::ofstream &bin_out_f) //binary output file
{
  //this is a minimal trr file writing routine that doesn't rely on \ 
  //the xdr library but only works with vmd in little endian systems

  //frame header, for more information on its contents see chemfiles
  int32_t header[18] = {1993, 1, 0, 
    0, 0, 0, 0, 0, 0, 0, 3*N*4, 0, 0, N, i_f, 0, 
    *(reinterpret_cast<int32_t *>(&t)), 0};
  bin_out_f.write(reinterpret_cast<char *>(header),sizeof(header));
  for (int i_p = 0; i_p<N; ++i_p) //particle index
  {
    bin_out_f.write(reinterpret_cast<char *>(&(r[i_p].x)),4);
    bin_out_f.write(reinterpret_cast<char *>(&(r[i_p].y)),4);
    bin_out_f.write(reinterpret_cast<char *>(&(r[i_p].z)),4);
  }
}

//read frame from binary file
void chrdat::read_frame_bin(std::ifstream &bin_inp_f) //binary input file
{

}

} //namespace mmc

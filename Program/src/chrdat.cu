//Includes

#include "chrdat.cuh" //chromatin data

//Namespace

namespace mmc //Marco Mendívil Carboni
{

//Functions

//chromatin data constructor
chrdat::chrdat(parmap &par) //parameters
  : N {par.get_val<uint>("number_of_particles",0)}
  , R {par.get_val<float>("confinement_radius",-1.0)}
  , T {par.get_val<float>("temperature",298.0)}
  , eps {par.get_val<float>("particle_energy",1.0)}
  , n_l {par.get_val<uint>("number_of_lbs",0)}
  , i_f {0}, t {0.0}
{
  //check parameters
  if (!(1<=N&&N<100'000)){ throw error("number_of_particles out of range");}
  if (!(0.0<R&&R<100.0)){ throw error("confinement_radius out of range");}
  if (!(0.0<T&&T<1'000.0)){ throw error("temperature out of range");}
  if (!(0.125<eps&&eps<2.0)){ throw error("particle_energy out of range");}
  if (!(n_l<100'000)){ throw error("number_of_lbs out of range");}
  float cvf = N*pow(0.5/(R-0.5),3.0); //chromatin volume fraction
  if (cvf>0.5){ throw error("chromatin volume fraction above 0.5");}
  float laf = n_l*pow(0.5/(R-1.0),2.0); //lbs area fraction
  if (laf>0.5){ throw error("lbs area fraction above 0.5");}
  std::string msg_1 = ""; //1st message
  msg_1 += "N = "+cnfs(N,5,'0')+" ";
  msg_1 += "R = "+cnfs(R,5,'0',2)+" ";
  msg_1 += "T = "+cnfs(T,5,'0',1)+" ";
  logger::record(msg_1);
  std::string msg_2 = ""; //2nd message
  msg_2 += "eps = "+cnfs(eps,5,'0',3)+" ";
  msg_2 += "n_l = "+cnfs(n_l,5,'0')+" ";
  msg_2 += "cvf = "+cnfs(cvf,5,'0',3)+" ";
  logger::record(msg_2);

  //allocate device memory
  cuda_check(cudaMalloc(&pt,N*sizeof(ptype)));
  cuda_check(cudaMalloc(&r,N*sizeof(float4)));
  cuda_check(cudaMalloc(&f,N*sizeof(float4)));
  cuda_check(cudaMalloc(&lr,n_l*sizeof(float4)));

  //allocate host memory
  cuda_check(cudaMallocHost(&hpt,N*sizeof(ptype)));
  cuda_check(cudaMallocHost(&hr,N*sizeof(float4)));
  cuda_check(cudaMallocHost(&hf,N*sizeof(float4)));
  cuda_check(cudaMallocHost(&hlr,n_l*sizeof(float4)));
}

//chromatin data destructor
chrdat::~chrdat()
{
  //deallocate device memory
  cuda_check(cudaFree(pt));
  cuda_check(cudaFree(r));
  cuda_check(cudaFree(f));
  cuda_check(cudaFree(lr));

  //deallocate host memory
  cuda_check(cudaFreeHost(hpt));
  cuda_check(cudaFreeHost(hr));
  cuda_check(cudaFreeHost(hf));
  cuda_check(cudaFreeHost(hlr));
}

//write frame to text file
void chrdat::write_frame_txt(std::ofstream &txt_out_f) //text output file
{
  txt_out_f<<"Chromatin simulation, i_f = 0, t = 0.0\n";
  txt_out_f<<cnfs(N,5,' ')<<"\n";
  std::string pts; //particle type string
  for (uint i_p = 0; i_p<N; ++i_p) //particle index
  {
    pts = (hpt[i_p]==LAD)?"X":"Y";//tmp -----------------------------------------
    txt_out_f<<std::setw(5)<<i_p+1<<std::left<<std::setw(5)<<pts;
    txt_out_f<<std::right<<std::setw(5)<<pts<<std::setw(5)<<i_p+1;
    txt_out_f<<cnfs(hr[i_p].x,8,' ',3);
    txt_out_f<<cnfs(hr[i_p].y,8,' ',3);
    txt_out_f<<cnfs(hr[i_p].z,8,' ',3);
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
  uint32_t header[18] = {1993, 1, 0, 
    0, 0, 0, 0, 0, 0, 0, 3*N*4, 0, 0, N, i_f, 0, 
    0, 0};
  bin_out_f.write(reinterpret_cast<char *>(header),sizeof(header));
  for (uint i_p = 0; i_p<N; ++i_p) //particle index
  {
    bin_out_f.write(reinterpret_cast<char *>(&(hr[i_p].x)),4);
    bin_out_f.write(reinterpret_cast<char *>(&(hr[i_p].y)),4);
    bin_out_f.write(reinterpret_cast<char *>(&(hr[i_p].z)),4);
  }
}

//read frame from binary file
void chrdat::read_frame_bin(std::ifstream &bin_inp_f) //binary input file
{

}

} //namespace mmc

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
  , n_l {par.get_val<uint>("number_of_lbs",0)}
  , i_f {0}, t {0.0}
{
  //check parameters
  if (!(1<=N&&N<100'000)){ throw error("number_of_particles out of range");}
  if (!(0.0<R&&R<100.0)){ throw error("confinement_radius out of range");}
  if (!(0.0<T&&T<1'000.0)){ throw error("temperature out of range");}
  if (!(n_l<100'000)){ throw error("number_of_lbs out of range");}
  float cvf = N*pow(0.5/(R-0.5),3.0); //chromatin volume fraction
  if (cvf>0.5){ throw error("chromatin volume fraction above 0.5");}
  float laf = n_l*pow(lco/(R-rco),2.0); //lbs area fraction
  if (laf>0.6){ throw error("lbs area fraction above 0.6");}
  std::string msg_1 = ""; //1st message
  msg_1 += "N = "+cnfs(N,5,'0')+" ";
  msg_1 += "R = "+cnfs(R,5,'0',2)+" ";
  msg_1 += "T = "+cnfs(T,5,'0',1)+" ";
  logger::record(msg_1);
  std::string msg_2 = ""; //2nd message
  msg_2 += "n_l = "+cnfs(n_l,5,'0')+" ";
  msg_2 += "cvf = "+cnfs(cvf,5,'0',3)+" ";
  msg_2 += "laf = "+cnfs(laf,5,'0',3)+" ";
  logger::record(msg_2);

  //allocate device memory
  cuda_check(cudaMalloc(&pt,N*sizeof(ptype)));
  cuda_check(cudaMalloc(&r,N*sizeof(vec3f)));
  cuda_check(cudaMalloc(&f,N*sizeof(vec3f)));
  cuda_check(cudaMalloc(&lr,n_l*sizeof(vec3f)));

  //allocate host memory
  cuda_check(cudaMallocHost(&hpt,N*sizeof(ptype)));
  cuda_check(cudaMallocHost(&hr,N*sizeof(vec3f)));
  cuda_check(cudaMallocHost(&hf,N*sizeof(vec3f)));
  cuda_check(cudaMallocHost(&hlr,n_l*sizeof(vec3f)));
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
  //write chromatin data
  txt_out_f<<i_f<<" "<<t<<"\n";
  txt_out_f<<cnfs(N,5,' ')<<"\n";
  char ptc; //particle type character
  for (uint i_p = 0; i_p<N; ++i_p) //particle index
  {
    ptc = (hpt[i_p]==LND)?'A':'B';
    txt_out_f<<std::setw(5)<<i_p+1<<std::left<<std::setw(5)<<ptc;
    txt_out_f<<std::right<<std::setw(5)<<ptc<<std::setw(5)<<i_p+1;
    txt_out_f<<cnfs(hr[i_p].x,8,' ',3);
    txt_out_f<<cnfs(hr[i_p].y,8,' ',3);
    txt_out_f<<cnfs(hr[i_p].z,8,' ',3);
    txt_out_f<<"\n";
  }
  txt_out_f<<cnfs(0.0,10,' ',5);
  txt_out_f<<cnfs(0.0,10,' ',5);
  txt_out_f<<cnfs(0.0,10,' ',5);
  txt_out_f<<"\n";

  //check filestream
  if (txt_out_f.fail())
  {
    throw mmc::error("failed to write frame to text file");
  }
}

//read frame from text file
void chrdat::read_frame_txt(std::ifstream &txt_inp_f) //text input file
{
  //read chromatin data
  std::string aux_str; //auxiliary string
  txt_inp_f>>i_f>>t;
  txt_inp_f>>aux_str;
  char ptc; //particle type character
  for (uint i_p = 0; i_p<N; ++i_p) //particle index
  {
    txt_inp_f>>aux_str;
    txt_inp_f>>ptc>>aux_str;
    hpt[i_p] = (ptc=='A')?LND:LAD;
    txt_inp_f>>hr[i_p].x;
    txt_inp_f>>hr[i_p].y;
    txt_inp_f>>hr[i_p].z;
  }
  txt_inp_f>>aux_str;
  txt_inp_f>>aux_str;
  txt_inp_f>>aux_str;

  //copy host arrays to device
  cuda_check(cudaMemcpy(pt,hpt,N*sizeof(ptype),cudaMemcpyHostToDevice));
  cuda_check(cudaMemcpy(r,hr,N*sizeof(vec3f),cudaMemcpyHostToDevice));

  //check filestream
  if (txt_inp_f.fail())
  {
    throw mmc::error("failed to read frame from text file");
  }
}

//write frame to binary file
void chrdat::write_frame_bin(std::ofstream &bin_out_f) //binary output file
{
  //------------------------------note------------------------------
  //this is a minimal trr file writing routine that doesn't rely on \ 
  //the xdr library but only works with vmd in little endian systems

  //frame header, for more information on its contents see chemfiles
  uint32_t header[18] = {1993, 1, 0, 
    0, 0, 0, 0, 0, 0, 0, 3*N*4, 0, 0, N, i_f, 0, 
    reinterpret_cast<uint32_t &>(t), 0};

  //write chromatin data
  bin_out_f.write(reinterpret_cast<char *>(header),sizeof(header));
  for (uint i_p = 0; i_p<N; ++i_p) //particle index
  {
    bin_out_f.write(reinterpret_cast<char *>(&(hr[i_p].x)),4);
    bin_out_f.write(reinterpret_cast<char *>(&(hr[i_p].y)),4);
    bin_out_f.write(reinterpret_cast<char *>(&(hr[i_p].z)),4);
  }

  //check filestream
  if (bin_out_f.fail())
  {
    throw mmc::error("failed to write frame to binary file");
  }
}

//read frame from binary file
void chrdat::read_frame_bin(std::ifstream &bin_inp_f) //binary input file
{
  //read chromatin data
  uint32_t header[18]; //frame header
  bin_inp_f.read(reinterpret_cast<char *>(header),sizeof(header));
  i_f = header[14]; t = reinterpret_cast<float &>(header[16]);
  for (uint i_p = 0; i_p<N; ++i_p) //particle index
  {
    bin_inp_f.read(reinterpret_cast<char *>(&(hr[i_p].x)),4);
    bin_inp_f.read(reinterpret_cast<char *>(&(hr[i_p].y)),4);
    bin_inp_f.read(reinterpret_cast<char *>(&(hr[i_p].z)),4);
  }

  //copy host position array to device
  cuda_check(cudaMemcpy(r,hr,N*sizeof(vec3f),cudaMemcpyHostToDevice));

  //check filestream
  if (bin_inp_f.fail())
  {
    throw mmc::error("failed to read frame from binary file");
  }
}

} //namespace mmc

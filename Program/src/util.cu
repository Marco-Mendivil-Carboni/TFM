//Includes

#include "util.cuh" //general utilities

#include <iostream> //standard input/output stream objects

#include <time.h> //time utilities library
#include <glob.h> //pathname pattern matching types

//Namespace

namespace mmc //Marco Mendívil Carboni
{

//Functions and Operators

//set log file and open it
void logger::set_file(
  const std::string &pathstr, //file path string
  bool ovr) //overwrite file
{
  logger &sin = get_instance(); //singleton instance
  if (sin.log_f.is_open()){ sin.log_f.close();}
  if (pathstr==""){ sin.w_f = false; return;}
  if (ovr){ sin.log_f.open(pathstr); sin.log_f.close();}
  else{ sin.log_f.open(pathstr,std::ios::app); sin.log_f.close();}
  sin.log_f.open(pathstr,std::ios::in|std::ios::ate);
  if (sin.log_f.is_open()){ sin.w_f = true;}
  else{ sin.w_f = false; std::cout<<"unable to open "<<pathstr<<std::endl;}
}

//log message with timestamp
void logger::record(const std::string &msg) //message
{
  logger &sin = get_instance(); //singleton instance
  time_t now = time(nullptr); //current time
  tm *now_info = localtime(&now); //curent time information
  char timestr[22]; //timestamp C-style string
  strftime(timestr,22,"[%d/%m/%y %H:%M:%S] ",now_info);
  if (sin.w_f)
  {
    sin.log_f<<timestr<<msg<<std::endl;
  }
  std::cout<<timestr<<msg<<std::endl;
}

//show progress percentage
void logger::show_prog_pc(float prog_pc) //progress percentage
{
  logger &sin = get_instance(); //singleton instance
  if (sin.w_f)
  {
    sin.log_f<<"progress: "<<cnfs(prog_pc,5,'0',1)<<"%";
    sin.log_f.seekp(-16,std::ios::cur);
  }
  std::cout<<"progress: "<<cnfs(prog_pc,5,'0',1)<<"%";
  std::cout<<"\r"; std::cout.flush();
}

//logger constructor
logger::logger() {}

//logger destructor
logger::~logger()
{
  log_f.close();
}

//return singleton instance
logger &logger::get_instance()
{
  static logger sin; //singleton instance
  return sin;
}

//error constructor
error::error(const std::string &msg) //error message
  : std::runtime_error(msg) {}

//parmap constructor
parmap::parmap(std::ifstream &par_f) //parameter file
{
  std::string key; //parameter key
  std::string val; //parameter value
  while (par_f>>key>>val)
  {
    insert({key,val});
  }
}

//new operator
void *mngdobj::operator new(size_t objsize) //object size
{
  void *obj_p; //object pointer
  cudaMallocManaged(&obj_p,objsize);
  return obj_p;
}

//delete operator
void mngdobj::operator delete(void *obj_p) //object pointer
{
  cudaFree(obj_p);
}

//check for errors in cuda runtime API call
void cuda_check(cudaError_t rtn_val) //cuda runtime API call return value
{
  if (rtn_val!=cudaSuccess)
  {
    std::string msg = "cuda: "; //error message
    msg += cudaGetErrorString(rtn_val);
    throw error(msg);
  }
}

//count files matching pattern
int glob_count(const std::string &pathpat) //file path pattern
{
  glob_t glob_sr; //glob search result
  int rtn_val = glob(pathpat.c_str(),0,nullptr,&glob_sr); //return value
  if (rtn_val!=0)
  {
    globfree(&glob_sr);
    if (rtn_val==GLOB_NOMATCH){ return 0;}
    else{ throw error("unable to find matches of "+pathpat);}
  }
  else
  {
    globfree(&glob_sr);
    return glob_sr.gl_pathc;
  }
}

} //namespace mmc

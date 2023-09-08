//Includes

#include "util.hpp" //utilities

#include <time.h> //time utilities library
#include <glob.h> //pathname pattern matching types

//Namespace

namespace mmcc //Marco Mendívil Carboni code
{

//Functions

//set log file and open it
void logger::set_file(
  const std::string &path, //log file path
  char mode) //log file openmode
{
  logger &sinlog = get_instance(); //singleton instance
  if (sinlog.file.is_open()){ sinlog.file.close();}
  if (path==""){ sinlog.w_f = false; return;}
  switch (mode)
  {
    case 'a': //create file and if it exists do not overwrite it
      sinlog.file.open(path,std::ios::app); sinlog.file.close();
      break;
    case 'w': //create file and if it exists overwrite it
      sinlog.file.open(path); sinlog.file.close();
      break; 
    default: //return if openmode is unknown
      std::cout<<"unknown log file openmode"<<std::endl;
      return;
  }
  sinlog.file.open(path,std::ios::in|std::ios::ate);
  if (sinlog.file.is_open()){ sinlog.w_f = true;}
  else{ sinlog.w_f = false; std::cout<<"unable to open "<<path<<std::endl;}
}

//log message with timestamp
void logger::record(const std::string &msg) //message
{
  logger &sinlog = get_instance(); //singleton instance
  time_t now = time(nullptr); //current time
  tm *now_info = localtime(&now); //curent time information
  char timestamp[22]; //timestamp C-style string
  strftime(timestamp,22,"[%d/%m/%y %H:%M:%S] ",now_info);
  if (sinlog.w_f)
  {
    sinlog.file<<timestamp<<msg<<std::endl;
  }
  std::cout<<timestamp<<msg<<std::endl;
}

//show progress percentage
void logger::show_prog_pc(float prog_pc) //progress percentage
{
  logger &sinlog = get_instance(); //singleton instance
  if (sinlog.w_f)
  {
    sinlog.file<<"progress: "<<cnfs(prog_pc,5,'0',1)<<"%";
    sinlog.file.seekp(-16,std::ios::cur);
  }
  std::cout<<"progress: "<<cnfs(prog_pc,5,'0',1)<<"%";
  std::cout<<"\r"; std::cout.flush();
}

//logger constructor
logger::logger() {}

//logger destructor
logger::~logger()
{
  file.close();
}

//return singleton instance
logger &logger::get_instance()
{
  static logger sinlog; //singleton logger
  return sinlog;
}

//error constructor
error::error(const std::string &msg) //error message
  : std::runtime_error{msg} {}

//count files matching pattern
int glob_count(const std::string &pattern) //file path pattern
{
  glob_t glob_res; //glob result
  int rtn_val = glob(pattern.c_str(),0,nullptr,&glob_res); //return value
  if (rtn_val!=0)
  {
    globfree(&glob_res);
    if (rtn_val==GLOB_NOMATCH){ return 0;}
    else{ throw error("unable to find matches of "+pattern);}
  }
  else
  {
    globfree(&glob_res);
    return glob_res.gl_pathc;
  }
}

} //namespace mmcc

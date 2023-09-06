//Includes

#include "util.hpp"

#include <time.h> //time utilities library
#include <glob.h> //pathname pattern matching types

//Namespace

namespace mmcc //Marco Mendívil Carboni code
{

//Functions

//set log file and open it
void logger::set_file(const std::string &path)
{
  logger &sinlog = get_instance(); //singleton instance
  if (sinlog.file.is_open()){ sinlog.file.close();}
  sinlog.file.open(path,std::ios::app);
  sinlog.file.close();
  sinlog.file.open(path,std::ios::in|std::ios::ate);
  if (!sinlog.file.is_open()){ std::cout<<"unable to open "<<path<<std::endl;}
}

//log message with timestamp
void logger::record(const std::string &msg)
{
  logger &sinlog = get_instance(); //singleton instance
  time_t now = time(nullptr);
  tm *now_info = localtime(&now); char timestamp[22];
  strftime(timestamp,22,"[%d/%m/%y %H:%M:%S] ",now_info);
  sinlog.file<<timestamp<<msg<<std::endl;
  std::cout<<timestamp<<msg<<std::endl;
}

//show progress percentage
void logger::show_prog_pc(float prog_pc)
{
  logger &sinlog = get_instance(); //singleton instance
  sinlog.file<<"progress: "<<cnfs(prog_pc,5,'0',1)<<"%";
  sinlog.file.seekp(-16,std::ios::cur);
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
error::error(const std::string &msg) : std::runtime_error{msg} {}

//count files matching pattern
int glob_count(std::string &pattern)
{
  glob_t glob_result;
  int return_value = glob(pattern.c_str(),0,nullptr,&glob_result);
  if (return_value!=0)
  {
    globfree(&glob_result);
    if (return_value==GLOB_NOMATCH){ return 0;}
    else{ throw error("unable to find matches of "+pattern);}
  }
  else
  {
    globfree(&glob_result);
    return glob_result.gl_pathc;
  }
}

} //namespace mmcc

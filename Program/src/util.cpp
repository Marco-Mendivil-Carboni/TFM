//Includes

#include "../inc/util.hpp"

#include <iostream> //standard input/output stream objects
#include <ctime> //time utilities library

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
  if (!sinlog.file.is_open()){ std::cout<<"unable to open"<<path<<std::endl;}
}

//log message with timestamp
void logger::record(const std::string &msg)
{
  logger &sinlog = get_instance(); //singleton instance
  std::time_t timestamp = time(nullptr);
  sinlog.file<<std::ctime(&timestamp)<<msg<<std::endl;
  std::cout<<std::ctime(&timestamp)<<msg<<std::endl;
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

} //namespace mmcc

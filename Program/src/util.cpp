//Includes

#include <cstdio> //standard input and output library
#include <ctime> //time utilities library
#include <stdexcept> //standard exceptions classes

#include "../inc/util.hpp"

//Namespace

namespace mmcc //Marco Mendívil Carboni code
{

//Functions

//error constructor
error::error(const std::string &msg) : std::runtime_error{msg} {}

//open file and throw error if it fails
FILE *fopen(const char *filename, const char *mode)
{
  FILE *f_ptr = std::fopen(filename,mode);
  if (f_ptr==nullptr)
  {
    char msg[512];
    std::snprintf(msg,sizeof(msg),"unable to open %s",filename);
    throw error(msg);
  }
  return f_ptr;
}

//write message to log file with timestamp
void log_message(FILE *f_ptr, const char *msg)
{
  time_t rt = time(nullptr); struct tm *rti = localtime(&rt);
  std::fprintf(f_ptr,"%02d:%02d:%02d ",rti->tm_hour,rti->tm_min,rti->tm_sec);
  std::fprintf(stdout,"%s.\n",msg);
  std::fprintf(f_ptr,"%s.\n",msg);
  std::fflush(f_ptr);
}

} //namespace mmcc

#ifndef MMCC_UTIL_H
#define MMCC_UTIL_H

//Includes

#include <iostream> //standard input/output stream objects
#include <fstream> //file stream classes
#include <stdexcept> //standard exceptions classes
#include <sstream> //string stream classes
#include <iomanip> //parametric manipulators

//Namespace

namespace mmcc //Marco Mendívil Carboni code
{

//Classes

class logger //basic logger
{
  public:
  
  //Functions

  //set log file and open it
  static void set_file(const std::string &path); //log file path

  //log message with timestamp
  static void record(const std::string &msg); //message

  //show progress percentage
  static void show_prog_pc(float prog_pc); //progress percentage

  private:

  //Variables

  std::ofstream file; //log file

  //Functions

  //logger constructor
  logger();

  //logger destructor
  ~logger();

  //return singleton instance
  static logger &get_instance();
};

class error : public std::runtime_error //generic exception type
{
  public:

    //Functions

    //error constructor
    error(const std::string &msg); //error message
};

//Functions

//count files matching pattern
int glob_count(std::string &pattern); //file path pattern

//Templates

//convert number to formatted string
template <typename T>
std::string cnfs
  (T num, //number
  int len = 0, //length
  char fillc = ' ', //filler character
  int prc = 0) //precision
{
  std::stringstream num_str; //number stringstream
  if (len>0){ num_str<<std::setw(len);}
  if (fillc!=' '){ num_str<<std::setfill(fillc);}
  if (prc>0){ num_str<<std::setprecision(prc)<<std::fixed;}
  num_str<<num; return num_str.str();
}

//check file is open or else throw
template <typename T>
void check_file
  (T &file, //file stream
  std::string &path) //file path string
{
  if (!file.is_open())
  {
    throw error("unable to open "+path);
  }
}

} //namespace mmcc

#endif //MMCC_UTIL_H

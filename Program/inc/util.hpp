#ifndef MMCC_UTIL_H
#define MMCC_UTIL_H

//Includes

#include <fstream> //file stream classes
#include <stdexcept> //standard exceptions classes
#include <string> //string types and functions

//Namespace

namespace mmcc //Marco Mendívil Carboni code
{

//Classes

class logger //basic logger
{
  public:
  
  //Functions

  //set log file and open it
  static void set_file(const std::string &path);

  //log message with timestamp
  static void record(const std::string &msg);

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
    error(const std::string &msg);
};

//Functions

//convert integer to fixed length string
std::string cits(int num, int len);

//Templates

//check file is open or else throw
template <typename T>
void check_file(T &file, std::string &path)
{
  if (!file.is_open())
  {
    throw error("unable to open "+path);
  }
}

} //namespace mmcc

#endif //MMCC_UTIL_H

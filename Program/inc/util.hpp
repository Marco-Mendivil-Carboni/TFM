#ifndef MMCC_UTIL_H
#define MMCC_UTIL_H

//Includes

#include <stdexcept> //standard exceptions classes

//Namespace

namespace mmcc //Marco Mendívil Carboni code
{

//Classes

class error : public std::runtime_error //generic exception type
{
  public:
    //error constructor
    error(const std::string &msg);
};

//Functions

//open file and throw error if it fails
FILE *fopen(const char *filename, const char *mode);

//write message to log file with timestamp
void log_message(FILE *f_ptr, const char *msg);

} //namespace mmcc

#endif //MMCC_UTIL_H

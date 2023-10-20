#ifndef MMC_CHRANA_H
#define MMC_CHRANA_H

//Includes

#include "chrdat.cuh" //chromatin data

//Namespace

namespace mmc //Marco Mendívil Carboni
{

//Classes

class chrana : public chrdat //chromatin analysis
{
  public:

  //Functions

  //chromatin analysis constructor
  chrana(parmap &par); //parameters

  //chromatin analysis destructor
  ~chrana();

  private:

  //Parameters and Variables

  const uint fpf; //frames per file

  //Functions
};

} //namespace mmc

#endif //MMC_CHRANA_H

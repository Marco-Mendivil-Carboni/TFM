//Includes

#include "chrana.cuh" //chromatin analysis

//Namespace

namespace mmc //Marco Mendívil Carboni
{

//Functions

//chromatin analysis constructor
chrana::chrana(parmap &par) //parameters
  : chrdat(par)
  , fpf {par.get_val<uint>("frames_per_file",100)}
{

}

//chromatin analysis destructor
chrana::~chrana()
{

}

} //namespace mmc

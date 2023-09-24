#ifndef MMC_SUGRID_H
#define MMC_SUGRID_H

//Namespace

namespace mmc //Marco Mendívil Carboni
{

//Structures

struct sugrid //sorted uniform grid
{
  //Variables

  const uint n_o; //number of objects
  const float csl; //grid cell side length
  const uint cps; //grid cells per side
  const uint n_c; //number of grid cells

  uint *uci; //unsorted grid cell index array
  uint *sci; //sorted grid cell index array
  uint *upi; //unsorted particle index array
  uint *spi; //sorted particle index array

  uint *beg; //grid cell beginning array
  uint *end; //grid cell end array

  void *eb; //extra buffer
  size_t ebs; //extra buffer size

  //Functions

  //sorted uniform grid constructor
  sugrid(
    const uint N, //number of particles
    const float csl, //grid cell side length
    const uint cps); //grid cells per side

  //sorted uniform grid destructor
  ~sugrid();

  //generate grid arrays
  void generate_arrays(
    int tpb, //threads per block
    float4 *r); //position array
};

} //namespace mmc

#endif //MMC_SUGRID_H

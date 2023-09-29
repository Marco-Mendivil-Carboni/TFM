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
  const float csl; //cell side length
  const uint cps; //cells per side
  const uint n_c; //number of cells

  uint *uci; //unsorted cell index array
  uint *sci; //sorted cell index array
  uint *upi; //unsorted particle index array
  uint *spi; //sorted particle index array

  uint *beg; //cell beginning array
  uint *end; //cell end array

  void *eb; //extra buffer
  size_t ebs; //extra buffer size

  cudaArray_t arr;
  cudaTextureObject_t tex;

  //Functions

  //sorted uniform grid constructor
  sugrid(
    const uint N, //number of particles
    const float csl, //cell side length
    const uint cps); //cells per side

  //sorted uniform grid destructor
  ~sugrid();

  //generate grid arrays
  void generate_arrays(
    int tpb, //threads per block
    float4 *r); //position array
};

} //namespace mmc

#endif //MMC_SUGRID_H

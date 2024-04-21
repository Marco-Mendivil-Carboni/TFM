#ifndef MMC_SUGRID_H
#define MMC_SUGRID_H

// Includes

#include "util.cuh" // general utilities

// Namespace

namespace mmc // Marco Mendívil Carboni
{

// Structures

struct sugrid // sorted uniform grid
{
  // Variables

  const uint n_o; // number of objects
  const float csl; // cell side length
  const uint cps; // cells per side
  const uint n_c; // number of cells

  uint *uci; // unsorted cell index array
  uint *sci; // sorted cell index array
  uint *uoi; // unsorted object index array
  uint *soi; // sorted object index array
  uint *beg; // cell beginning array
  uint *end; // cell end array

  void *eb; // extra buffer
  size_t ebs; // extra buffer size

  // Functions

  // sorted uniform grid constructor
  sugrid(const uint n_o, // number of objects
      const float csl, // cell side length
      const uint cps); // cells per side

  // sorted uniform grid delegating constructor
  sugrid(const uint n_o, // number of objects
      const sugrid &g); // grid

  // sorted uniform grid destructor
  ~sugrid();

  // generate grid arrays
  void generate_arrays(int tpb, // threads per block
      vec3f *r); // position array
};

} // namespace mmc

#endif // MMC_SUGRID_H

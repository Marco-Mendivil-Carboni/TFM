#ifndef MMC_CHRVIS_H
#define MMC_CHRVIS_H

//Includes

#include "chrdat.cuh" //chromatin data

//Namespace

namespace mmc //Marco Mendívil Carboni
{

//Structures

struct ray //light ray
{
  vec3f org; //origin
  vec3f dir; //direction
};

// struct camera
// {
//   vec3f origin;
//   vec3f lower_left_corner;
//   vec3f horizontal;
//   vec3f vertical;

//   __device__ camera(vec3f lookfrom, vec3f lookat, vec3f vup, float vfov, float aspect);

//   __device__ ray get_ray(float u, float v);
// };

//Classes

class chrvis : public chrdat //chromatin visualization
{
  public:

  //Functions

  //chromatin visualization constructor
  chrvis(parmap &par); //parameters

  //chromatin visualization destructor
  ~chrvis();

  //render frame and save image
  void render_frame(const std::string &out_dir); //output directory

  private:

  //Parameters and Variables

  // int nx = 1200;
  // int ny = 600;
  // int ns = 100;
  // int tx = 8;
  // int ty = 8;

  // vec3f *fb;

  // curandState *d_rand_state;

  // hitable **d_list;
  // hitable **d_world;

  // camera **d_camera;

  //Functions
};

} //namespace mmc

#endif //MMC_CHRVIS_H

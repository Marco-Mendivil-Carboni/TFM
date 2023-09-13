#ifndef MMC_UTIL_H
#define MMC_UTIL_H

//Includes

#include <fstream> //file stream classes
#include <iomanip> //input/output parametric manipulators
#include <map> //map container classes

//Namespace

namespace mmc //Marco Mendívil Carboni
{

//Classes

class logger //basic logger
{
  public:
  
  //Functions

  //set log file and open it
  static void set_file(
    const std::string &pathstr, //file path string
    bool ovr = false); //overwrite file

  //log message with timestamp
  static void record(const std::string &msg); //message

  //show progress percentage
  static void show_prog_pc(float prog_pc); //progress percentage

  private:

  //Variables

  bool w_f = false; //write output to file
  std::ofstream log_f; //log file

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

class parmap : public std::map<std::string,std::string> //parameter map
{
  public:

  //Functions

  //parmap constructor
  parmap(std::ifstream &par_f); //parameter file

  //Templates

  //get parameter value
  template <typename T>
  T get_val(
    std::string key, //parameter key
    T def_val) //default value
  {
    T val; //parameter value
    if (find(key)==end()){ val = def_val;}
    else
    {
      std::stringstream{at(key)}>>val;
    }
    return val;
  }
};

//Functions

//count files matching pattern
int glob_count(const std::string &pathpat); //file path pattern

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
  (T &s_f, //stream file 
  const std::string &pathstr) //file path string
{
  if (!s_f.is_open())
  {
    throw error("unable to open "+pathstr);
  }
}

} //namespace mmc

//Inline Functions and Operators for float3 & float4

inline __host__ __device__ float3 make_float3(float4 a)
{
  return make_float3(a.x,a.y,a.z);
}
inline __host__ __device__ float4 make_float4(float3 a)
{
  return make_float4(a.x,a.y,a.z,0.0f);
}
inline __host__ __device__ float3 operator-(float3 &a)
{
  return make_float3(-a.x,-a.y,-a.z);
}
inline __host__ __device__ float4 operator-(float4 &a)
{
  return make_float4(-a.x,-a.y,-a.z,-a.w);
}
inline __host__ __device__ float3 operator+(float3 a, float3 b)
{
  return make_float3(a.x+b.x,a.y+b.y,a.z+b.z);
}
inline __host__ __device__ float4 operator+(float4 a, float4 b)
{
  return make_float4(a.x+b.x,a.y+b.y,a.z+b.z,a.w+b.w);
}
inline __host__ __device__ void operator+=(float3 &a, float3 b)
{
  a.x+=b.x; a.y+=b.y; a.z+=b.z;
}
inline __host__ __device__ void operator+=(float4 &a, float4 b)
{
  a.x+=b.x; a.y+=b.y; a.z+=b.z; a.w+=b.w;
}
inline __host__ __device__ float3 operator-(float3 a, float3 b)
{
  return make_float3(a.x-b.x,a.y-b.y,a.z-b.z);
}
inline __host__ __device__ float4 operator-(float4 a, float4 b)
{
  return make_float4(a.x-b.x,a.y-b.y,a.z-b.z,a.w-b.w);
}
inline __host__ __device__ void operator-=(float3 &a, float3 b)
{
  a.x-=b.x; a.y-=b.y; a.z-=b.z;
}
inline __host__ __device__ void operator-=(float4 &a, float4 b)
{
  a.x-=b.x; a.y-=b.y; a.z-=b.z; a.w-=b.w;
}
inline __host__ __device__ float3 operator*(float3 a, float b)
{
  return make_float3(a.x*b,a.y*b,a.z*b);
}
inline __host__ __device__ float3 operator*(float b, float3 a)
{
  return make_float3(b*a.x,b*a.y,b*a.z);
}
inline __host__ __device__ float4 operator*(float4 a, float b)
{
  return make_float4(a.x*b,a.y*b,a.z*b,a.w*b);
}
inline __host__ __device__ float4 operator*(float b, float4 a)
{
  return make_float4(b*a.x,b*a.y,b*a.z,b*a.w);
}
inline __host__ __device__ float3 operator/(float3 a, float b)
{
  return make_float3(a.x/b,a.y/b,a.z/b);
}
inline __host__ __device__ float4 operator/(float4 a, float b)
{
  return make_float4(a.x/b,a.y/b,a.z/b,a.w/b);
}
inline __host__ __device__ float dot(float3 a, float3 b)
{
  return a.x*b.x+a.y*b.y+a.z*b.z;
}
inline __host__ __device__ float dot(float4 a, float4 b)
{
  return a.x*b.x+a.y*b.y+a.z*b.z+a.w*b.w;
}
inline __host__ __device__ float length(float3 v)
{
  return sqrtf(dot(v,v));
}
inline __host__ __device__ float3 normalize(float3 v)
{
  float inv_len = rsqrtf(dot(v,v));
  return v*inv_len;
}
inline __host__ __device__ float3 cross(float3 a, float3 b)
{
  return make_float3(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);
}

#endif //MMC_UTIL_H

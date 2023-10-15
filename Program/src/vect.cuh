#ifndef VECOPS_H
#define VECOPS_H

//Structures

struct __align__(16) vec3f //vector of 3 floats
{
  float x, y, z;
};

struct __align__(16) vec3i //vector of 3 ints
{
  int x, y, z;
};

//Operators

inline __host__ __device__ float3 operator-(float3 &a)
{
  return {-a.x,-a.y,-a.z};
}
inline __host__ __device__ float3 operator+(float3 a, float3 b)
{
  return {a.x+b.x,a.y+b.y,a.z+b.z};
}
inline __host__ __device__ void operator+=(float3 &a, float3 b)
{
  a.x+=b.x; a.y+=b.y; a.z+=b.z;
}
inline __host__ __device__ float3 operator-(float3 a, float3 b)
{
  return {a.x-b.x,a.y-b.y,a.z-b.z};
}
inline __host__ __device__ void operator-=(float3 &a, float3 b)
{
  a.x-=b.x; a.y-=b.y; a.z-=b.z;
}
inline __host__ __device__ float3 operator*(float3 a, float b)
{
  return {a.x*b,a.y*b,a.z*b};
}
inline __host__ __device__ float3 operator*(float b, float3 a)
{
  return {b*a.x,b*a.y,b*a.z};
}
inline __host__ __device__ void operator*=(float3 &a, float b)
{
  a.x*=b; a.y*=b; a.z*=b;
}
inline __host__ __device__ float3 operator/(float3 a, float b)
{
  return {a.x/b,a.y/b,a.z/b};
}
inline __host__ __device__ void operator/=(float3 &a, float b)
{
  a.x/=b; a.y/=b; a.z/=b;
}

//Functions

//integer floorf
inline __host__ __device__ int ifloorf(float num) //number
{
  return floorf(num);
}

inline __host__ __device__ float dot(float3 a, float3 b)
{
  return a.x*b.x+a.y*b.y+a.z*b.z;
}
inline __host__ __device__ float length(float3 v)
{
  return sqrtf(dot(v,v));
}
inline __host__ __device__ float3 normalize(float3 v)
{
  return v*rsqrtf(dot(v,v));
}
inline __host__ __device__ float3 cross(float3 a, float3 b)
{
  return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};
}
inline __host__ __device__ int3 ifloorf(float3 v)
{
  return {ifloorf(v.x),ifloorf(v.y),ifloorf(v.z)};
}
inline __host__ __device__ float3 fabs(float3 v)
{
  return {fabs(v.x),fabs(v.y),fabs(v.z)};
}

#endif //VECOPS_H

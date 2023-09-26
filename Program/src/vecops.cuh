#ifndef VECOPS_H
#define VECOPS_H

//Inline Functions and Operators

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
inline __host__ __device__ void operator*=(float3 &a, float b)
{
  a.x*=b; a.y*=b; a.z*=b;
}
inline __host__ __device__ float4 operator*(float4 a, float b)
{
  return make_float4(a.x*b,a.y*b,a.z*b,a.w*b);
}
inline __host__ __device__ float4 operator*(float b, float4 a)
{
  return make_float4(b*a.x,b*a.y,b*a.z,b*a.w);
}
inline __host__ __device__ void operator*=(float4 &a, float b)
{
  a.x*=b; a.y*=b; a.z*=b; a.w*=b;
}
inline __host__ __device__ float3 operator/(float3 a, float b)
{
  return make_float3(a.x/b,a.y/b,a.z/b);
}
inline __host__ __device__ void operator/=(float3 &a, float b)
{
  a.x/=b; a.y/=b; a.z/=b;
}
inline __host__ __device__ float4 operator/(float4 a, float b)
{
  return make_float4(a.x/b,a.y/b,a.z/b,a.w/b);
}
inline __host__ __device__ void operator/=(float4 &a, float b)
{
  a.x/=b; a.y/=b; a.z/=b; a.w/=b;
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
inline __host__ __device__ int3 floorf(float3 v)
{
  return make_int3(floorf(v.x),floorf(v.y),floorf(v.z));
}
inline __host__ __device__ int4 floorf(float4 v)
{
  return make_int4(floorf(v.x),floorf(v.y),floorf(v.z),floorf(v.w));
}

#endif //VECOPS_H

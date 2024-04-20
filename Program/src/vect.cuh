#ifndef MMC_VECT_H
#define MMC_VECT_H

// Namespace

namespace mmc // Marco Mendívil Carboni
{

// Structures

struct vec3f // vector of 3 floats
{
  float x, y, z;
};

struct vec3i // vector of 3 ints
{
  int x, y, z;
};

// Operators and Functions

inline __host__ __device__ vec3f operator-(vec3f &a)
{
  return {-a.x, -a.y, -a.z};
}
inline __host__ __device__ vec3f operator+(vec3f a, vec3f b)
{
  return {a.x + b.x, a.y + b.y, a.z + b.z};
}
inline __host__ __device__ void operator+=(vec3f &a, vec3f b)
{
  a.x += b.x;
  a.y += b.y;
  a.z += b.z;
}
inline __host__ __device__ vec3f operator-(vec3f a, vec3f b)
{
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}
inline __host__ __device__ void operator-=(vec3f &a, vec3f b)
{
  a.x -= b.x;
  a.y -= b.y;
  a.z -= b.z;
}
inline __host__ __device__ vec3f operator*(vec3f a, float b)
{
  return {a.x * b, a.y * b, a.z * b};
}
inline __host__ __device__ vec3f operator*(float b, vec3f a)
{
  return {b * a.x, b * a.y, b * a.z};
}
inline __host__ __device__ void operator*=(vec3f &a, float b)
{
  a.x *= b;
  a.y *= b;
  a.z *= b;
}
inline __host__ __device__ vec3f operator/(vec3f a, float b)
{
  return {a.x / b, a.y / b, a.z / b};
}
inline __host__ __device__ void operator/=(vec3f &a, float b)
{
  a.x /= b;
  a.y /= b;
  a.z /= b;
}
inline __host__ __device__ float dot(vec3f a, vec3f b)
{
  return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline __host__ __device__ float length(vec3f v) { return sqrtf(dot(v, v)); }
inline __host__ __device__ vec3f normalized(vec3f v)
{
  return v * rsqrtf(dot(v, v));
}
inline __host__ __device__ vec3f cross(vec3f a, vec3f b)
{
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline __host__ __device__ int ifloorf(float x) { return floorf(x); }
inline __host__ __device__ vec3i ifloorc(vec3f v)
{
  return {ifloorf(v.x), ifloorf(v.y), ifloorf(v.z)};
}
inline __host__ __device__ vec3f fabsc(vec3f v)
{
  return {fabsf(v.x), fabsf(v.y), fabsf(v.z)};
}

} // namespace mmc

#endif // MMC_VECT_H

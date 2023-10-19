//Includes

#include "chrvis.cuh" //chromatin visualization

#include <float.h> //----------------------------------------
#include <curand_kernel.h> //cuRAND device functions

//Namespace

namespace mmc //Marco Mendívil Carboni
{

using prng = curandStatePhilox4_32_10; //PRNG type

class ray
{
  public:

  __device__ ray() {}
  __device__ ray(const vec3f &a, const vec3f &b) { A = a; B = b; }
  __device__ vec3f origin() const       { return A; }
  __device__ vec3f direction() const    { return B; }
  __device__ vec3f point_at_parameter(float t) const { return A + t*B; }

  vec3f A;
  vec3f B;
};

class camera {
  public:

  __device__ camera(vec3f lookfrom, vec3f lookat, vec3f vup, float vfov, float aspect)
  {
    vec3f u, v, w;
    float theta = vfov*M_PI/180;
    float half_height = tan(theta/2);
    float half_width = aspect * half_height;
    origin = lookfrom;
    w = normalized(lookfrom - lookat);
    u = normalized(cross(vup, w));
    v = cross(w, u);
    lower_left_corner = origin - half_width*u -half_height*v - w;
    horizontal = 2*half_width*u;
    vertical = 2*half_height*v;
  }

  __device__ ray get_ray(float u, float v)
  {
    return ray(origin, lower_left_corner + u*horizontal + v*vertical - origin);
  }

  vec3f origin;
  vec3f lower_left_corner;
  vec3f horizontal;
  vec3f vertical;
};

struct hit_record;

class material;

struct hit_record
{
    float t;
    vec3f p;
    vec3f normal;
    material *mat_ptr;
};

class hitable  {
    public:
        __device__ virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const = 0;
};

class hitable_list: public hitable  {
    public:
        __device__ hitable_list() {}
        __device__ hitable_list(hitable **l, int n) {list = l; list_size = n; }
        __device__ virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
        hitable **list;
        int list_size;
};

__device__ bool hitable_list::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
        hit_record temp_rec;
        bool hit_anything = false;
        float closest_so_far = t_max;
        for (int i = 0; i < list_size; i++) {
            if (list[i]->hit(r, t_min, closest_so_far, temp_rec)) {
                hit_anything = true;
                closest_so_far = temp_rec.t;
                rec = temp_rec;
            }
        }
        return hit_anything;
}

__device__ float schlick(float cosine, float ref_idx) {
    float r0 = (1.0f-ref_idx) / (1.0f+ref_idx);
    r0 = r0*r0;
    return r0 + (1.0f-r0)*pow((1.0f - cosine),5.0f);
}

__device__ bool refract(const vec3f& v, const vec3f& n, float ni_over_nt, vec3f& refracted) {
    vec3f uv = normalized(v);
    float dt = dot(uv, n);
    float discriminant = 1.0f - ni_over_nt*ni_over_nt*(1-dt*dt);
    if (discriminant > 0) {
        refracted = ni_over_nt*(uv - n*dt) - n*sqrt(discriminant);
        return true;
    }
    else
        return false;
}

__device__ vec3f random_in_unit_sphere(curandState *local_rand_state)
{
    vec3f p;
    do
    {  
      p = 2.0f*vec3f{curand_uniform(local_rand_state),curand_uniform(local_rand_state),curand_uniform(local_rand_state)} - vec3f{1,1,1};
    } while (dot(p,p) >= 1.0f);
    return p;
}

__device__ vec3f reflect(const vec3f& v, const vec3f& n) {
     return v - 2.0f*dot(v,n)*n;
}

class material  {
    public:
        __device__ virtual bool scatter(const ray& r_in, const hit_record& rec, vec3f& attenuation, ray& scattered, curandState *local_rand_state) const = 0;
};

class lambertian : public material {
    public:
        __device__ lambertian(const vec3f& a) : albedo(a) {}
        __device__ virtual bool scatter(const ray& r_in, const hit_record& rec, vec3f& attenuation, ray& scattered, curandState *local_rand_state) const  {
             vec3f target = rec.p + rec.normal + random_in_unit_sphere(local_rand_state);
             scattered = ray(rec.p, target-rec.p);
             attenuation = albedo;
             return true;
        }

        vec3f albedo;
};

class metal : public material {
    public:
        __device__ metal(const vec3f& a, float f) : albedo(a) { if (f < 1) fuzz = f; else fuzz = 1; }
        __device__ virtual bool scatter(const ray& r_in, const hit_record& rec, vec3f& attenuation, ray& scattered, curandState *local_rand_state) const  {
            vec3f reflected = reflect(normalized(r_in.direction()), rec.normal);
            scattered = ray(rec.p, reflected + fuzz*random_in_unit_sphere(local_rand_state));
            attenuation = albedo;
            return (dot(scattered.direction(), rec.normal) > 0.0f);
        }
        vec3f albedo;
        float fuzz;
};

class dielectric : public material {
public:
    __device__ dielectric(float ri) : ref_idx(ri) {}
    __device__ virtual bool scatter(const ray& r_in,
                         const hit_record& rec,
                         vec3f& attenuation,
                         ray& scattered,
                         curandState *local_rand_state) const  {
        vec3f outward_normal;
        vec3f reflected = reflect(r_in.direction(), rec.normal);
        float ni_over_nt;
        attenuation = vec3f{1.0, 1.0, 1.0};
        vec3f refracted;
        float reflect_prob;
        float cosine;
        if (dot(r_in.direction(), rec.normal) > 0.0f) {
            outward_normal = -1.0*rec.normal;
            ni_over_nt = ref_idx;
            cosine = dot(r_in.direction(), rec.normal) / length(r_in.direction());
            cosine = sqrt(1.0f - ref_idx*ref_idx*(1-cosine*cosine));
        }
        else {
            outward_normal = rec.normal;
            ni_over_nt = 1.0f / ref_idx;
            cosine = -dot(r_in.direction(), rec.normal) / length(r_in.direction());
        }
        if (refract(r_in.direction(), outward_normal, ni_over_nt, refracted))
            reflect_prob = schlick(cosine, ref_idx);
        else
            reflect_prob = 1.0f;
        if (curand_uniform(local_rand_state) < reflect_prob)
            scattered = ray(rec.p, reflected);
        else
            scattered = ray(rec.p, refracted);
        return true;
    }

    float ref_idx;
};

class sphere: public hitable  {
    public:
        __device__ sphere() {}
        __device__ sphere(vec3f cen, float r, material *m) : center(cen), radius(r), mat_ptr(m)  {};
        __device__ virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
        vec3f center;
        float radius;
        material *mat_ptr;
};

__device__ bool sphere::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
    vec3f oc = r.origin() - center;
    float a = dot(r.direction(), r.direction());
    float b = dot(oc, r.direction());
    float c = dot(oc, oc) - radius*radius;
    float discriminant = b*b - a*c;
    if (discriminant > 0) {
        float temp = (-b - sqrt(discriminant))/a;
        if (temp < t_max && temp > t_min) {
            rec.t = temp;
            rec.p = r.point_at_parameter(rec.t);
            rec.normal = (rec.p - center) / radius;
            rec.mat_ptr = mat_ptr;
            return true;
        }
        temp = (-b + sqrt(discriminant)) / a;
        if (temp < t_max && temp > t_min) {
            rec.t = temp;
            rec.p = r.point_at_parameter(rec.t);
            rec.normal = (rec.p - center) / radius;
            rec.mat_ptr = mat_ptr;
            return true;
        }
    }
    return false;
}

__device__ vec3f color(const ray& r, hitable **world, curandState *local_rand_state) {
    ray cur_ray = r;
    vec3f cur_attenuation = vec3f{1.0,1.0,1.0};
    for(int i = 0; i < 50; i++) {
        hit_record rec;
        if ((*world)->hit(cur_ray, 0.001f, FLT_MAX, rec)) {
            ray scattered;
            vec3f attenuation;
            if(rec.mat_ptr->scatter(cur_ray, rec, attenuation, scattered, local_rand_state)) {
                cur_attenuation = multc(cur_attenuation,attenuation);
                cur_ray = scattered;
            }
            else {
                return vec3f{0.0,0.0,0.0};
            }
        }
        else {
            vec3f unit_direction = normalized(cur_ray.direction());
            float t = 0.5f*(unit_direction.y + 1.0f);
            vec3f c = (1.0f-t)*vec3f{1.0, 1.0, 1.0} + t*vec3f{0.5, 0.7, 1.0};
            return multc(cur_attenuation,c);
        }
    }
    return vec3f{0.0,0.0,0.0}; // exceeded recursion
}

__global__ void render_init(int max_x, int max_y, curandState *rand_state) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    if((i >= max_x) || (j >= max_y)) return;
    int pixel_index = j*max_x + i;
    //Each thread gets same seed, a different sequence number, no offset
    curand_init(1984, pixel_index, 0, &rand_state[pixel_index]);
}

__global__ void render(vec3f *fb, int max_x, int max_y, int ns, camera **cam, hitable **world, curandState *rand_state) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    if((i >= max_x) || (j >= max_y)) return;
    int pixel_index = j*max_x + i;
    curandState local_rand_state = rand_state[pixel_index];
    vec3f col{0,0,0};
    for(int s=0; s < ns; s++) {
        float u = float(i + curand_uniform(&local_rand_state)) / float(max_x);
        float v = float(j + curand_uniform(&local_rand_state)) / float(max_y);
        ray r = (*cam)->get_ray(u,v);
        col += color(r, world, &local_rand_state);
    }
    rand_state[pixel_index] = local_rand_state;
    col /= float(ns);
    col.x = sqrt(col.x);
    col.y = sqrt(col.y);
    col.z = sqrt(col.z);
    fb[pixel_index] = col;
}

__global__ void create_world(hitable **d_list, hitable **d_world, camera **d_camera, int nx, int ny) {
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        d_list[0] = new sphere(vec3f{0,0,-1}, 0.5,
                                new lambertian(vec3f{0.1, 0.2, 0.5}));
        d_list[1] = new sphere(vec3f{0,-100.5,-1}, 100,
                                new lambertian(vec3f{0.8, 0.8, 0.0}));
        d_list[2] = new sphere(vec3f{1,0,-1}, 0.5,
                                new metal(vec3f{0.8, 0.6, 0.2}, 0.0));
        d_list[3] = new sphere(vec3f{-1,0,-1}, 0.5,
                                 new dielectric(1.5));
        d_list[4] = new sphere(vec3f{-1,0,-1}, -0.45,
                                 new dielectric(1.5));
        *d_world = new hitable_list(d_list,5);
        *d_camera   = new camera(vec3f{-2,2,1},
                                 vec3f{0,0,-1},
                                 vec3f{0,1,0},
                                 20.0,
                                 float(nx)/float(ny));
    }
}

__global__ void free_world(hitable **d_list, hitable **d_world, camera **d_camera) {
    for(int i=0; i < 5; i++) {
        delete ((sphere *)d_list[i])->mat_ptr;
        delete d_list[i];
    }
    delete *d_world;
    delete *d_camera;
}

chrvis::chrvis(parmap &par) //parameters
  : chrdat(par)
{
    int num_pixels = nx*ny;
    size_t fb_size = num_pixels*sizeof(vec3f);

    cuda_check(cudaMallocManaged((void **)&fb, fb_size));

    cuda_check(cudaMalloc((void **)&d_rand_state, num_pixels*sizeof(curandState)));

    // make our world of hitables & the camera
    cuda_check(cudaMalloc((void **)&d_list, 5*sizeof(hitable *)));
    cuda_check(cudaMalloc((void **)&d_world, sizeof(hitable *)));
    cuda_check(cudaMalloc((void **)&d_camera, sizeof(camera *)));
    create_world<<<1,1>>>(d_list, d_world, d_camera, nx, ny);
    cuda_check(cudaGetLastError());
    cuda_check(cudaDeviceSynchronize());
}

chrvis::~chrvis()
{
      // clean up
    cuda_check(cudaDeviceSynchronize());
    free_world<<<1,1>>>(d_list,d_world,d_camera);
    cuda_check(cudaGetLastError());
    cuda_check(cudaFree(d_camera));
    cuda_check(cudaFree(d_world));
    cuda_check(cudaFree(d_list));
    cuda_check(cudaFree(d_rand_state));
    cuda_check(cudaFree(fb));
    cudaDeviceReset();
}

void chrvis::render_frame()
{
    // Render our buffer
    dim3 blocks(nx/tx+1,ny/ty+1);
    dim3 threads(tx,ty);
    render_init<<<blocks, threads>>>(nx, ny, d_rand_state);
    cuda_check(cudaGetLastError());
    cuda_check(cudaDeviceSynchronize());
    render<<<blocks, threads>>>(fb, nx, ny,  ns, d_camera, d_world, d_rand_state);
    cuda_check(cudaGetLastError());
    cuda_check(cudaDeviceSynchronize());

    // Output FB as Image
    // std::cout << "P3\n" << nx << " " << ny << "\n255\n";
    // for (int j = ny-1; j >= 0; j--) {
    //     for (int i = 0; i < nx; i++) {
    //         size_t pixel_index = j*nx + i;
    //         int ir = int(255.99*fb[pixel_index].r());
    //         int ig = int(255.99*fb[pixel_index].g());
    //         int ib = int(255.99*fb[pixel_index].b());
    //         std::cout << ir << " " << ig << " " << ib << "\n";
    //     }
    // }
}

} //namespace mmc

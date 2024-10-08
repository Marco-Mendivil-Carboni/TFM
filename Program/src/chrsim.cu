// Includes

#include "chrsim.cuh" // chromatin simulation

// Namespace

namespace mmc // Marco Mendívil Carboni
{

// Device Functions

// calculate bonded forces
inline __device__ void calc_bf(
    const uint N, // number of particles
    ptype *pt, // particle type array
    uint i_p, // particle index
    vec3f *r, // position array
    vec3f *f) // force array
{
  // declare auxiliary variables
  vec3f vec[4]; // bond vectors
  float il[4]; // bond inverse lengths
  float cos[3]; // bond angle cosines
  vec3f bf = {0.0, 0.0, 0.0}; // bonded forces

  // determine bond limits
  uint lim_l = 2; // lower bond limit
  uint lim_u = N; // upper bond limit
  if (N == N_def) // take into account chromosome limits
  {
    for (uint i_c = 0; i_c < n_chr; ++i_c) // chromosome index
    {
      if (i_p < chrla[i_c + 1]) // particle belongs to chromosome
      {
        lim_l = chrla[i_c] + 2;
        lim_u = chrla[i_c + 1];
        break;
      }
    }
  }

  // calculate bond vectors, inverse lengths and angle cosines
  for (uint i_b = 0; i_b < 4; ++i_b) // bond index
  {
    if ((i_p + i_b) >= lim_l &&
        (i_p + i_b) <= lim_u) // calculate variables if bond exists
    {
      vec[i_b] = r[i_p + i_b - 1] - r[i_p + i_b - 2];
      il[i_b] = rsqrtf(dot(vec[i_b], vec[i_b]));
    }
    else // set variables to zero if bond doesn't exist
    {
      vec[i_b] = {0.0, 0.0, 0.0};
      il[i_b] = 0.0;
    }
  }
  for (uint i_c = 0; i_c < 3; ++i_c) // cosine index
  {
    cos[i_c] = dot(vec[i_c + 1], vec[i_c]) * il[i_c + 1] * il[i_c];
  }

  // calculate elastic potential force
  bf += k_e * (+(1.0 - l_0 * il[2]) * vec[2] - (1.0 - l_0 * il[1]) * vec[1]);

  // calculate bending potential force only on LADh particles
  if (pt[i_p] == LADh)
  {
    bf += k_b * (+il[1] * il[0] * vec[0] - cos[0] * vec[1] * il[1] * il[1]);
    bf += k_b * (+il[1] * il[2] * vec[2] - cos[1] * vec[1] * il[1] * il[1]);
    bf += k_b * (-il[2] * il[1] * vec[1] + cos[1] * vec[2] * il[2] * il[2]);
    bf += k_b * (-il[2] * il[3] * vec[3] + cos[2] * vec[2] * il[2] * il[2]);
  }

  // add result to force array
  f[i_p] += bf;
}

// calculate wall-particle force
inline __device__ void calc_wpf(
    vec3f vwp, // wall-particle vector
    vec3f &cf) // confinement force
{
  // calculate Wang-Frenkel force
  float dwp = length(vwp); // wall-particle distance
  if (dwp > rco) { return; }
  float d2 = dwp * dwp; // dwp squared
  cf += (18.0 * d2 * d2 - 96.0 * d2 + 96.0) / (d2 * d2 * d2 * d2) * vwp;
}

// calculate confinement force
template <stype T>
inline __device__ void calc_cf(
    const cngeom ng, // nucleus geometry
    uint i_p, // particle index
    vec3f *r, // position array
    vec3f *f) // force array
{
  // declare auxiliary variables
  vec3f r_i = r[i_p]; // particle position
  vec3f r_b = {0.0, 0.0, ng.nod + ng.bod}; // bleb position
  vec3f r_r = r_i - r_b; // particle position relative to bleb
  vec3f vwp; // wall-particle vector
  float d_r; // radial distance
  vec3f cf = {0.0, 0.0, 0.0}; // confinement force

  // find confinement region and calculate vwp
  if (r_i.z < ng.nod) // inside nucleus sphere
  {
    d_r = length(r_i);
    if (+r_i.z / d_r < ng.noc) // outside opening
    {
      vwp = -r_i * (ng.R_n / d_r - 1.0);
    }
    else // inside opening
    {
      d_r = sqrtf(r_i.x * r_i.x + r_i.y * r_i.y);
      vwp.x = -r_i.x * (ng.R_o / d_r - 1.0);
      vwp.y = -r_i.y * (ng.R_o / d_r - 1.0);
      vwp.z = r_i.z - ng.nod;
    }
  }
  else // inside nuclear bleb
  {
    d_r = length(r_r);
    if (-r_r.z / d_r < ng.boc) // outside opening
    {
      vwp = -r_r * (ng.R_b / d_r - 1.0);
    }
    else // inside opening
    {
      d_r = sqrtf(r_r.x * r_r.x + r_r.y * r_r.y);
      vwp.x = -r_r.x * (ng.R_o / d_r - 1.0);
      vwp.y = -r_r.y * (ng.R_o / d_r - 1.0);
      vwp.z = r_r.z + ng.bod;
    }
  }
  if (!isfinite(length(vwp))) { return; }

  // calculate wall-particle force
  calc_wpf(vwp, cf);

  // add result to force array
  f[i_p] += cf;
}

// calculate confinement force
template <>
inline __device__ void calc_cf<ICG>(
    const cngeom ng, // nucleus geometry
    uint i_p, // particle index
    vec3f *r, // position array
    vec3f *f) // force array
{
  // declare auxiliary variables
  vec3f r_i = r[i_p]; // particle position
  vec3f vwp; // wall-particle vector
  float d_r; // radial distance
  vec3f cf = {0.0, 0.0, 0.0}; // confinement force

  // calculate vwp
  d_r = length(r_i);
  vwp = -r_i * (ng.R_n / d_r - 1.0);
  if (!isfinite(length(vwp))) { return; }

  // calculate wall-particle force
  calc_wpf(vwp, cf);

  // add result to force array
  f[i_p] += cf;
}

// calculate particle-particle force
template <stype T>
inline __device__ void calc_ppf(
    vec3f vpp, // particle-particle vector
    bool hhi, // LADh-LADh interaction
    vec3f &srf) // short-range forces
{
  // calculate Wang-Frenkel force
  float dpp = length(vpp); // particle-particle distance
  if (dpp > aco) { return; }
  if (!hhi && dpp > rco) { return; }
  float d2 = dpp * dpp; // dpp squared
  srf += e_p * (18.0 * d2 * d2 - 96.0 * d2 + 96.0) / (d2 * d2 * d2 * d2) * vpp;
}

// calculate particle-particle force
template <>
inline __device__ void calc_ppf<ICG>( // calculate particle-particle force
    vec3f vpp, // particle-particle vector
    bool hhi, // LADh-LADh interaction
    vec3f &srf) // short-range forces
{
  // calculate Soft-Repulsive force
  float dpp = length(vpp); // particle-particle distance
  if (dpp > rco) { return; }
  srf += 128.0 * (3.0 * rco - 3.0 * dpp) * vpp;
}

// calculate lbs-particle force
template <stype T>
inline __device__ void calc_lpf(
    vec3f vlp, // lbs-particle vector
    vec3f &srf) // short-range forces
{
  // calculate Binding force
  float dlp = length(vlp); // lbs-particle distance
  if (dlp > lco) { return; }
  float nd = (dlp / lco); // normalized dlp
  float d2 = nd * nd; // normalized dlp squared
  srf += e_l * (8.0 / (3.0 * lco * lco)) * (d2 * d2 * d2 - 1.0) * vlp;
}

// calculate lbs-particle force
template <>
inline __device__ void calc_lpf<ICG>(
    vec3f vlp, // lbs-particle vector
    vec3f &srf) // short-range forces
{
}

// calculate short-range forces with cell's objects
template <stype T>
inline __device__ void calc_cell_srf(
    ptype *pt, // particle type array
    vec3f *lr, // lbs position array
    sugrid *pgp, // particle grid pointer
    sugrid *lgp, // lbs grid pointer
    uint i_p, // particle index
    vec3f r_i, // particle position
    uint i_c, // cell index
    vec3f *r, // position array
    vec3f &srf) // short-range forces
{
  // declare auxiliary variables
  uint pgbeg = pgp->beg[i_c]; // particle grid cell beginning
  uint pgend = pgp->end[i_c]; // particle grid cell end
  uint lgbeg = lgp->beg[i_c]; // lbs grid cell beginning
  uint lgend = lgp->end[i_c]; // lbs grid cell end
  uint j_p; // secondary particle index
  uint i_l; // lbs index

  // range over cell's particles
  for (uint sai = pgbeg; sai < pgend; ++sai) // sorted array index
  {
    // get secondary particle index
    j_p = pgp->soi[sai];

    // calculate particle-particle force only between non-bonded particles
    if (((j_p > i_p) ? j_p - i_p : i_p - j_p) > 1)
    {
      // calculate particle-particle force
      vec3f vpp = r_i - r[j_p]; // particle-particle vector
      bool hhi = false; // LADh-LADh interaction
      if (pt[i_p] == LADh && pt[j_p] == LADh) { hhi = true; }
      calc_ppf<T>(vpp, hhi, srf);
    }
  }

  // calculate lbs-particle force only on LADh particles
  if (pt[i_p] == LADh)
  {
    // range over cell's lbs
    for (uint sai = lgbeg; sai < lgend; ++sai) // sorted array index
    {
      // get lbs index
      i_l = lgp->soi[sai];

      // calculate lbs-particle force
      vec3f vlp = r_i - lr[i_l]; // lbs-particle vector
      calc_lpf<T>(vlp, srf);
    }
  }
}

// calculate short-range forces
template <stype T>
inline __device__ void calc_srf(
    ptype *pt, // particle type array
    vec3f *lr, // lbs position array
    sugrid *pgp, // particle grid pointer
    sugrid *lgp, // lbs grid pointer
    uint i_p, // particle index
    vec3f *r, // position array
    vec3f *f) // force array
{
  // declare auxiliary variables
  vec3f r_i = r[i_p]; // particle position
  const float csl = pgp->csl; // cell side length
  const uint cps = pgp->cps; // cells per side
  const uint n_c = pgp->n_c; // number of cells
  vec3i ir = ifloorc(r_i / csl); // integer coordinates
  uint iofst = (cps / 2) * (1 + cps + cps * cps); // index offset
  vec3f srf = {0.0, 0.0, 0.0}; // short-range forces

  // range over neighbouring cells
  uint nci; // neighbour cell index
  vec3i nir; // neighbour integer coordinates
  for (nir.x = ir.x - 1; nir.x <= ir.x + 1; ++nir.x)
  {
    for (nir.y = ir.y - 1; nir.y <= ir.y + 1; ++nir.y)
    {
      for (nir.z = ir.z - 1; nir.z <= ir.z + 1; ++nir.z)
      {
        // calculate neighbour cell index
        nci = iofst + nir.x + nir.y * cps + nir.z * cps * cps;
        if (nci >= n_c) { continue; }

        // calculate short-range forces with cell's objects
        calc_cell_srf<T>(pt, lr, pgp, lgp, i_p, r_i, nci, r, srf);
      }
    }
  }

  // add result to force array
  f[i_p] += srf;
}

// Global Functions

// initialize PRNG state array
__global__ void init_ps(
    const uint N, // number of particles
    prng *ps, // PRNG state array
    uint pseed) // PRNG seed
{
  // calculate particle index
  uint i_p = blockIdx.x * blockDim.x + threadIdx.x; // particle index
  if (i_p >= N) { return; }

  // initialize PRNG state
  curand_init(pseed, i_p, 0, &ps[i_p]);
}

// begin Runge-Kutta iteration
__global__ void begin_iter(
    const uint N, // number of particles
    vec3f *f, // force array
    vec3f *ef, // extra force array
    float sd, // standard deviation
    vec3f *rn, // random number array
    prng *ps) // PRNG state array
{
  // calculate particle index
  uint i_p = blockIdx.x * blockDim.x + threadIdx.x; // particle index
  if (i_p >= N) { return; }

  // calculate random numbers
  vec3f az; // absolute z-score
  do
  {
    rn[i_p].x = sd * curand_normal(&ps[i_p]);
    rn[i_p].y = sd * curand_normal(&ps[i_p]);
    rn[i_p].z = sd * curand_normal(&ps[i_p]);
    az = fabsc(rn[i_p] / sd);
  } while (az.x > 5 || az.y > 5 || az.z > 5);

  // initialize forces to zero
  f[i_p] = {0.0, 0.0, 0.0};
  ef[i_p] = {0.0, 0.0, 0.0};
}

// execute 1st stage of the Runge-Kutta method
template <stype T>
__global__ void exec_RK_1(
    const uint N, // number of particles
    const cngeom ng, // nucleus geometry
    ptype *pt, // particle type array
    vec3f *r, // position array
    vec3f *f, // force array
    vec3f *lr, // lbs position array
    vec3f *er, // extra position array
    vec3f *rn, // random number array
    sugrid *pgp, // particle grid pointer
    sugrid *lgp) // lbs grid pointer
{
  // calculate particle index
  uint i_p = blockIdx.x * blockDim.x + threadIdx.x; // particle index
  if (i_p >= N) { return; }

  // calculate forces
  calc_bf(N, pt, i_p, r, f);
  calc_cf<T>(ng, i_p, r, f);
  calc_srf<T>(pt, lr, pgp, lgp, i_p, r, f);

  // calculate extra position
  er[i_p] = r[i_p] + f[i_p] * dt + rn[i_p];
}

// execute 2nd stage of the Runge-Kutta method
template <stype T>
__global__ void exec_RK_2(
    const uint N, // number of particles
    const cngeom ng, // nucleus geometry
    ptype *pt, // particle type array
    vec3f *r, // position array
    vec3f *f, // force array
    vec3f *lr, // lbs position array
    vec3f *er, // extra position array
    vec3f *ef, // extra force array
    vec3f *rn, // random number array
    sugrid *pgp, // particle grid pointer
    sugrid *lgp) // lbs grid pointer
{
  // calculate particle index
  uint i_p = blockIdx.x * blockDim.x + threadIdx.x; // particle index
  if (i_p >= N) { return; }

  // calculate forces
  calc_bf(N, pt, i_p, er, ef);
  calc_cf<T>(ng, i_p, er, ef);
  calc_srf<T>(pt, lr, pgp, lgp, i_p, er, ef);

  // calculate new position
  r[i_p] = r[i_p] + 0.5 * (ef[i_p] + f[i_p]) * dt + rn[i_p];
}

// Host Functions

// chromatin simulation constructor
chrsim::chrsim(parmap &par) // parameters
    : chrdat(par), spf{par.get_val<uint>("steps_per_frame", 1.0 / dt)},
      tpb{par.get_val<uint>("threads_per_block", 128)},
      sd{static_cast<float>(sqrt(2.0 * k_B * T * dt))},
      pg(N, aco, 2 * ceil(ng.d_m / aco)), lg(n_l, pg)
{
  // check parameters
  if (!(1 <= spf && spf < 10'000))
  {
    throw error("steps_per_frame out of range");
  }
  if (!(1 <= tpb && tpb < 1'025))
  {
    throw error("threads_per_block out of range");
  }
  std::string msg = ""; // message
  msg += "spf = " + cnfs(spf, 4, '0') + " ";
  msg += "tpb = " + cnfs(tpb, 4, '0') + " ";
  logger::record(msg);

  // allocate device memory
  cuda_check(cudaMalloc(&er, N * sizeof(vec3f)));
  cuda_check(cudaMalloc(&ef, N * sizeof(vec3f)));
  cuda_check(cudaMalloc(&rn, N * sizeof(vec3f)));
  cuda_check(cudaMalloc(&ps, N * sizeof(prng)));
  cuda_check(cudaMalloc(&pgp, sizeof(sugrid)));
  cuda_check(cudaMalloc(&lgp, sizeof(sugrid)));

  // copy grids to device
  cuda_check(cudaMemcpy(pgp, &pg, sizeof(sugrid), cudaMemcpyHostToDevice));
  cuda_check(cudaMemcpy(lgp, &lg, sizeof(sugrid), cudaMemcpyHostToDevice));

  // initialize PRNG
  init_ps<<<(N + tpb - 1) / tpb, tpb>>>(N, ps, time(nullptr));
}

// chromatin simulation destructor
chrsim::~chrsim()
{
  // deallocate device memory
  cuda_check(cudaFree(er));
  cuda_check(cudaFree(ef));
  cuda_check(cudaFree(rn));
  cuda_check(cudaFree(ps));
  cuda_check(cudaFree(pgp));
  cuda_check(cudaFree(lgp));
}

// generate a random initial condition
void chrsim::generate_initial_condition()
{
  // initialize host PRNG
  curandGenerator_t gen; // host PRNG
  curandCreateGeneratorHost(&gen, CURAND_RNG_PSEUDO_DEFAULT);
  curandSetPseudoRandomGeneratorSeed(gen, time(nullptr));

  // set random lbs positions
  set_lbs_positions(gen);

  // set particle type sequence
  set_particle_types(gen);

  // set random particle positions
  set_particle_positons(gen);

  // free host PRNG
  curandDestroyGenerator(gen);

  // generate lbs grid arrays
  lg.generate_arrays(tpb, lr);

  // separate beads
  uint po = particle_overlaps(); // particle overlaps
  uint ipo = po; // initial particle overlaps
  while (po > 0)
  {
    // show separation progress
    logger::show_prog_pc(100.0 * (ipo - po) / ipo);

    // iterate over all steps per frame
    for (uint fsi = 0; fsi < spf; ++fsi) // frame step index
    {
      // make one Runge-Kutta iteration
      begin_iter<<<(N + tpb - 1) / tpb, tpb>>>(N, f, ef, 0.0, rn, ps);
      pg.generate_arrays(tpb, r);
      exec_RK_1<ICG>
          <<<(N + tpb - 1) / tpb, tpb>>>(N, ng, pt, r, f, lr, er, rn, pgp, lgp);
      pg.generate_arrays(tpb, er);
      exec_RK_2<ICG><<<(N + tpb - 1) / tpb, tpb>>>(
          N, ng, pt, r, f, lr, er, ef, rn, pgp, lgp);
    }

    // copy position array to host
    cuda_check(cudaMemcpy(hr, r, N * sizeof(vec3f), cudaMemcpyDeviceToHost));

    // count particle overlaps
    po = particle_overlaps();
  }

  // record success message
  logger::record("initial condition generated");
}

// save simulation state to binary file
void chrsim::save_checkpoint(std::ofstream &bin_out_f) // binary output file
{
  // write simulation data
  bin_out_f.write(reinterpret_cast<char *>(&i_f), sizeof(i_f));
  bin_out_f.write(reinterpret_cast<char *>(&t), sizeof(t));
  bin_out_f.write(reinterpret_cast<char *>(hpt), N * sizeof(ptype));
  bin_out_f.write(reinterpret_cast<char *>(hr), N * sizeof(vec3f));
  bin_out_f.write(reinterpret_cast<char *>(hlr), n_l * sizeof(vec3f));

  // record success message
  logger::record("simulation checkpoint saved");
}

// load simulation state from binary file
void chrsim::load_checkpoint(std::ifstream &bin_inp_f) // binary input file
{
  // read simulation data
  bin_inp_f.read(reinterpret_cast<char *>(&i_f), sizeof(i_f));
  bin_inp_f.read(reinterpret_cast<char *>(&t), sizeof(t));
  bin_inp_f.read(reinterpret_cast<char *>(hpt), N * sizeof(ptype));
  bin_inp_f.read(reinterpret_cast<char *>(hr), N * sizeof(vec3f));
  bin_inp_f.read(reinterpret_cast<char *>(hlr), n_l * sizeof(vec3f));

  // copy host arrays to device
  cuda_check(cudaMemcpy(pt, hpt, N * sizeof(ptype), cudaMemcpyHostToDevice));
  cuda_check(cudaMemcpy(r, hr, N * sizeof(vec3f), cudaMemcpyHostToDevice));
  cuda_check(cudaMemcpy(lr, hlr, n_l * sizeof(vec3f), cudaMemcpyHostToDevice));

  // record success message
  logger::record("simulation checkpoint loaded");
}

// run simulation and write trajectory to binary file
void chrsim::run_simulation(std::ofstream &bin_out_f) // binary output file
{
  // generate lbs grid arrays
  lg.generate_arrays(tpb, lr);

  // iterate over all frames per file
  for (uint ffi = 0; ffi < fpf; ++ffi) // file frame index
  {
    // show simulation progress
    logger::show_prog_pc(100.0 * ffi / fpf);

    // iterate over all steps per frame
    for (uint fsi = 0; fsi < spf; ++fsi) // frame step index
    {
      // make one Runge-Kutta iteration
      begin_iter<<<(N + tpb - 1) / tpb, tpb>>>(N, f, ef, sd, rn, ps);
      pg.generate_arrays(tpb, r);
      exec_RK_1<DST>
          <<<(N + tpb - 1) / tpb, tpb>>>(N, ng, pt, r, f, lr, er, rn, pgp, lgp);
      pg.generate_arrays(tpb, er);
      exec_RK_2<DST><<<(N + tpb - 1) / tpb, tpb>>>(
          N, ng, pt, r, f, lr, er, ef, rn, pgp, lgp);
    }

    // copy position array to host
    cuda_check(cudaMemcpy(hr, r, N * sizeof(vec3f), cudaMemcpyDeviceToHost));

    // write trajectory frame
    ++i_f;
    t += spf * dt;
    write_frame_bin(bin_out_f);
  }

  // record success message
  logger::record("simulation ended");
}

// set random lbs positions
void chrsim::set_lbs_positions(curandGenerator_t gen) // host PRNG
{
  // declare auxiliary variables
  float ran; // random number in (0,1]
  float theta; // polar angle
  float phi; // azimuthal angle
  vec3f ran_dir; // random direction

  // set lbs positions randomly
  for (uint i_l = 0; i_l < n_l; ++i_l) // lbs index
  {
    // generate random direction
    curandGenerateUniform(gen, &ran, 1);
    theta = acos(1.0 - 2.0 * ran);
    curandGenerateUniform(gen, &ran, 1);
    phi = 2.0 * M_PI * ran;
    ran_dir = {sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};

    // calculate position of next lbs
    hlr[i_l] = (ng.R_n - rco) * ran_dir;

    // check if position is acceptable
    bool p_a = true; // position is acceptable
    for (uint j_l = 0; j_l < i_l; ++j_l) // secondary lbs index
    {
      float dll = length(hlr[j_l] - hlr[i_l]); // lbs-lbs distance
      if (dll < 2.0 * lco) { p_a = false; }
    }
    if (hlr[i_l].z / length(hlr[i_l]) > ng.noc) { p_a = false; }

    if (!p_a) { --i_l; } // repeat
  }

  // copy host lbs position array to device
  cuda_check(cudaMemcpy(lr, hlr, n_l * sizeof(vec3f), cudaMemcpyHostToDevice));
}

// set particle type sequence
void chrsim::set_particle_types(curandGenerator_t gen) // host PRNG
{
  if (N == N_def) // use experimental data to set particles types
  {
    std::ifstream txt_inp_f; // text input file
    std::string pathstr = seqpath_def; // file path string
    txt_inp_f.open(pathstr);
    check_file(txt_inp_f, pathstr);
    char ptc; // particle type character
    for (uint i_p = 0; i_p < N; ++i_p) // particle index
    {
      txt_inp_f >> ptc;
      hpt[i_p] = (ptc == 'A') ? LADh : LNDe;
    }
    if (txt_inp_f.fail()) { throw mmc::error("failed to read sequence"); }
    txt_inp_f.close();
  }
  else // set particle types randomly
  {
    float ran; // random number in (0,1]
    float edl; // exponential domain length
    uint cde = 0; // current domain end
    ptype cpt = LNDe; // current particle type
    for (uint i_p = 0; i_p < N; ++i_p) // particle index
    {
      if (i_p == cde) // change domain type
      {
        do
        {
          curandGenerateUniform(gen, &ran, 1);
          edl = -mdl * log(ran);
        } while ((edl / mdl) > 5 || edl < 1.0);
        cde = i_p + edl;
        cpt = (cpt == LADh) ? LNDe : LADh;
      }
      hpt[i_p] = cpt;
    }
  }

  // copy host particle type array to device
  cuda_check(cudaMemcpy(pt, hpt, N * sizeof(ptype), cudaMemcpyHostToDevice));
}

// set random particle positions
void chrsim::set_particle_positons(curandGenerator_t gen) // host PRNG
{
  if (N == N_def) // perform a confined random walk for each chromosome
  {
    // declare auxiliary variables
    float phi = (1.0 + sqrt(5.0)) / 2.0; // golden ratio
    float mda = M_PI / 2.0 - atan(phi); // maximum direction angle
    vec3f dir[] = // directions
        {{0.0, +1.0, +phi}, {+1.0, +phi, 0.0}, {+phi, 0.0, +1.0},
         {0.0, +1.0, -phi}, {+1.0, -phi, 0.0}, {-phi, 0.0, +1.0},
         {0.0, -1.0, +phi}, {-1.0, +phi, 0.0}, {+phi, 0.0, -1.0},
         {0.0, -1.0, -phi}, {-1.0, -phi, 0.0}, {-phi, 0.0, -1.0}};
    uint n_dir = sizeof(dir) / sizeof(vec3f); // number of directions
    float ran; // random number in (0,1]
    uint i_r; // random index
    vec3f aux_dir; // auxiliary direction

    // shuffle directions randomly
    for (uint i_d = n_dir - 1; i_d > 0; --i_d) // direction index
    {
      curandGenerateUniform(gen, &ran, 1);
      i_r = (i_d + 1) * (1.0 - ran);
      aux_dir = dir[i_d];
      dir[i_d] = dir[i_r];
      dir[i_r] = aux_dir;
    }

    // perform confined random walks for each chromosome with random directions
    for (uint i_c = 0; i_c < n_chr; ++i_c) // chromosome index
    {
      uint i_d = i_c % n_dir; // direction index
      perform_random_walk(gen, hchrla[i_c], hchrla[i_c + 1], dir[i_d], mda);
    }
  }
  else // perform a single confined random walk
  {
    perform_random_walk(gen, 0, N, {0.0, 0.0, 0.0}, M_PI);
  }
}

// perform confined random walk
void chrsim::perform_random_walk(
    curandGenerator_t gen, // host PRNG
    uint i_s, // starting index
    uint i_e, // ending index
    vec3f dir, // direction
    float mda) // maximum direction angle
{
  // declare auxiliary variables
  float len_d = length(dir); // direction length
  float iT = 1.0 / (k_B * T); // inverse temperature
  float ran; // random number in (0,1]
  float theta; // polar angle
  float phi; // azimuthal angle
  float d_r; // radial distance
  float len_b; // bond length
  float angle_b; // bond angle
  vec3f ran_dir; // random direction
  vec3f old_dir; // old direction
  vec3f per_dir; // perpendicular direction
  vec3f new_dir; // new direction
  bool p_a; // position is acceptable

  // place first particle
  do
  {
    curandGenerateUniform(gen, &ran, 1);
    theta = acos(1.0 - 2.0 * ran);
    curandGenerateUniform(gen, &ran, 1);
    phi = 2.0 * M_PI * ran;
    curandGenerateUniform(gen, &ran, 1);
    d_r = (ng.R_n - mis) * pow(ran, 1.0 / 3);
    ran_dir = {sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};
    hr[i_s] = d_r * ran_dir;
    old_dir = ran_dir;
  } while (acos(dot(hr[i_s], dir) / (d_r * len_d)) > mda);

  // place the rest of particles
  uint att = 0; // number of attempts
  for (uint i_p = i_s + 1; i_p < i_e; ++i_p) // particle index
  {
    // generate random direction perpendicular to old direction
    curandGenerateUniform(gen, &ran, 1);
    theta = acos(1.0 - 2.0 * ran);
    curandGenerateUniform(gen, &ran, 1);
    phi = 2.0 * M_PI * ran;
    ran_dir = {sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};
    per_dir = cross(old_dir, ran_dir);
    per_dir = normalized(per_dir);

    // generate random bond angle and calculate new direction
    curandGenerateUniform(gen, &ran, 1);
    angle_b = acos(1.0 + log(1.0 - (1.0 - exp(-1.0 * iT)) * ran) / (0.5 * iT));
    new_dir = cos(angle_b) * old_dir + sin(angle_b) * per_dir;

    // calculate position of next particle
    curandGenerateUniform(gen, &ran, 1);
    len_b = l_0 + sqrt(2.0 / (k_e * iT)) * erfinv(2.0 * ran - 1.0);
    hr[i_p] = len_b * new_dir + hr[i_p - 1];

    // check if new position is acceptable
    p_a = true;
    if (!isfinite(hr[i_p].x)) { p_a = false; }
    if (!isfinite(hr[i_p].y)) { p_a = false; }
    if (!isfinite(hr[i_p].z)) { p_a = false; }
    d_r = length(hr[i_p]);
    if ((ng.R_n - d_r) < mis) { p_a = false; }
    if (acos(dot(hr[i_p], dir) / (d_r * len_d)) > mda) { p_a = false; }

    if (p_a) // continue
    {
      att = 0;
      old_dir = new_dir;
    }
    else // go back
    {
      ++att;
      if (att > 1024) { i_p = i_s + 1 + (i_p - i_s) * 3 / 4; }
      else { --i_p; }
    }
  }

  // copy host position array to device
  cuda_check(cudaMemcpy(r, hr, N * sizeof(vec3f), cudaMemcpyHostToDevice));
}

// count particle overlaps
uint chrsim::particle_overlaps()
{
  // iterate over all pairs of non-bonded particles
  uint po = 0; // particle overlaps
  float dpp; // particle-particle distance
  for (uint i_p = 0; i_p < N; ++i_p) // particle index
  {
    for (uint j_p = 0; (j_p + 1) < i_p; ++j_p) // secondary particle index
    {
      // check if particles overlap
      dpp = length(hr[j_p] - hr[i_p]);
      if (dpp < mis) { ++po; }
    }
  }
  return po;
}

} // namespace mmc

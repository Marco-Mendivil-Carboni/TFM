// Includes

#include "chrana.cuh" // chromatin analysis

// Namespace

namespace mmc // Marco Mendívil Carboni
{

// Global Functions

// calculate the spatial distance and contact probability
__global__ void calc_sd_cp(
    const uint N, // number of particles
    uint lsdcp, // length of sd array and cp array
    uint i_s, // starting index
    uint i_e, // ending index
    vec3f *r, // position array
    float *sd, // spatial distance array
    float *cp) // contact probability array
{
  // calculate array index
  int i_a = blockIdx.x * blockDim.x + threadIdx.x; // array index
  if (i_a >= lsdcp) { return; }

  // calculate the spatial distance and contact probability
  sd[i_a] = 0.0;
  cp[i_a] = 0.0;
  uint s = i_a + 1; // separation
  uint n_chr_c = (N == N_def) ? 2 : 1; // number of chromosome copies
  float dpp; // particle-particle distance
  for (uint i_c = 0; i_c < n_chr_c; ++i_c) // copy index
  {
    for (uint i_p = i_s; i_p < (i_e - s); ++i_p) // particle index
    {
      dpp = length(r[i_p + s] - r[i_p]);
      sd[i_a] += dpp;
      if (dpp < aco) { cp[i_a] += cf; }
    }
    i_s += N / n_chr_c;
    i_e += N / n_chr_c;
  }
  sd[i_a] /= n_chr_c * (i_e - i_s - s);
  cp[i_a] /= n_chr_c * (i_e - i_s - s);
}

// calculate the contact map
__global__ void calc_cm(
    const uint N, // number of particles
    const uint lcm, // length of cm array
    vec3f *r, // position array
    float *cm) // contact map array
{
  // calculate array index
  int i_a = blockIdx.x * blockDim.x + threadIdx.x; // array index
  if (i_a >= lcm) { return; }

  // calculate the contact map
  cm[i_a] = 0.0;
  uint i_x = sqrtf(2.0 * i_a + 0.25) - 0.5; // x index
  uint i_y = i_a - i_x * (i_x + 1.0) / 2.0; // y index
  uint i_s_x = i_x * px_sz; // starting x index
  uint i_e_x = i_s_x + px_sz; // ending x index
  uint i_s_y = i_y * px_sz; // starting y index
  uint i_e_y = i_s_y + px_sz; // ending y index
  uint n_chr_c = (N == N_def) ? 2 : 1; // number of chromosome copies
  float dpp; // particle-particle distance
  for (uint i_c = 0; i_c < n_chr_c; ++i_c) // copy index
  {
    for (uint i_p_x = i_s_x; i_p_x < i_e_x; ++i_p_x) // x particle index
    {
      for (uint i_p_y = i_s_y; i_p_y < i_e_y; ++i_p_y) // y particle index
      {
        dpp = length(r[i_p_x] - r[i_p_y]);
        if (dpp < aco) { cm[i_a] += cf; }
      }
    }
    i_s_x += N / n_chr_c;
    i_e_x += N / n_chr_c;
    i_s_y += N / n_chr_c;
    i_e_y += N / n_chr_c;
  }
  cm[i_a] /= n_chr_c * px_sz * px_sz;
}

// Host Functions

// chromatin analysis constructor
chrana::chrana(parmap &par) // parameters
    : chrdat(par), tpb{par.get_val<uint>("threads_per_block", 128)},
      cms{(N == N_def) ? (N / 2) / px_sz : N / px_sz}, lcm{cms * (cms + 1) / 2}
{
  // allocate memory
  for (uint i_c = 0; i_c < n_chr_h; ++i_c) // chromosome index
  {
    if (N == N_def) { lsdcp[i_c] = (hchrla[i_c + 1] - hchrla[i_c]) / 2; }
    else if (i_c == 0) { lsdcp[i_c] = N / 2; }
    else { lsdcp[i_c] = 0; }
    sd_bo[i_c] = new simobs_b[lsdcp[i_c]];
    cp_bo[i_c] = new simobs_b[lsdcp[i_c]];
  }
  cm_bo = new simobs_b[lcm];

  // allocate device memory
  for (uint i_c = 0; i_c < n_chr_h; ++i_c) // chromosome index
  {
    cuda_check(cudaMalloc(&sd[i_c], lsdcp[i_c] * sizeof(float)));
    cuda_check(cudaMalloc(&cp[i_c], lsdcp[i_c] * sizeof(float)));
  }
  cuda_check(cudaMalloc(&cm, lcm * sizeof(float)));

  // allocate host memory
  for (uint i_c = 0; i_c < n_chr_h; ++i_c) // chromosome index
  {
    cuda_check(cudaMallocHost(&hsd[i_c], lsdcp[i_c] * sizeof(float)));
    cuda_check(cudaMallocHost(&hcp[i_c], lsdcp[i_c] * sizeof(float)));
  }
  cuda_check(cudaMallocHost(&hcm, lcm * sizeof(float)));
}

// chromatin analysis destructor
chrana::~chrana()
{
  // deallocate memory
  for (uint i_c = 0; i_c < n_chr_h; ++i_c) // chromosome index
  {
    delete[] sd_bo[i_c];
    delete[] cp_bo[i_c];
  }
  delete[] cm_bo;

  // deallocate device memory
  for (uint i_c = 0; i_c < n_chr_h; ++i_c) // chromosome index
  {
    cuda_check(cudaFree(sd[i_c]));
    cuda_check(cudaFree(cp[i_c]));
  }
  cuda_check(cudaFree(cm));

  // deallocate host memory
  for (uint i_c = 0; i_c < n_chr_h; ++i_c) // chromosome index
  {
    cuda_check(cudaFreeHost(hsd[i_c]));
    cuda_check(cudaFreeHost(hcp[i_c]));
  }
  cuda_check(cudaFreeHost(hcm));
}

// add initial condition to analysis
void chrana::add_initial_condition(std::ifstream &txt_inp_f) // text input file
{
  // read initial condition frame
  read_frame_txt(txt_inp_f);

  // calculate observables
  calc_observables();
}

// add trajectory file to analysis
void chrana::add_trajectory_file(std::ifstream &bin_inp_f) // binary input file
{
  // iterate over all frames per file
  for (uint ffi = 0; ffi < fpf; ++ffi) // file frame index
  {
    // show analysis progress
    logger::show_prog_pc(100.0 * ffi / fpf);

    // read trajectory frame
    read_frame_bin(bin_inp_f);

    // calculate observables
    calc_observables();
  }
}

// calculate last individual simulation statistics
void chrana::calc_last_is_stat()
{
  dcm_o.calc_last_is_stat();
  rg2_o.calc_last_is_stat();
  nop_o.calc_last_is_stat();
  ncf_o.calc_last_is_stat();
  for (uint i_t = 0; i_t < 3; ++i_t) // type index
  {
    for (uint i_b = 0; i_b < n_b; ++i_b) // bin index
    {
      rcd_o[i_t][i_b].calc_last_is_stat();
    }
  }
  for (uint i_c = 0; i_c < n_chr_h; ++i_c) // chromosome index
  {
    for (uint i_a = 0; i_a < lsdcp[i_c]; ++i_a) // array index
    {
      sd_bo[i_c][i_a].calc_last_is_stat();
      cp_bo[i_c][i_a].calc_last_is_stat();
    }
  }
  for (uint i_a = 0; i_a < lcm; ++i_a) // array index
  {
    cm_bo[i_a].calc_last_is_stat();
  }
}

// save last individual simulation statistics
void chrana::save_last_is_stat(std::ofstream &txt_out_f) // text output file
{
  // save dcm, rg2, nop and ncf last individual simulation statistics
  txt_out_f << "         avg   sqrt(var)         sem   f_n_b     i_t ter\n";
  dcm_o.save_last_is_stat(txt_out_f);
  rg2_o.save_last_is_stat(txt_out_f);
  nop_o.save_last_is_stat(txt_out_f);
  ncf_o.save_last_is_stat(txt_out_f);
  txt_out_f << "\n\n";

  // save rcd last individual simulation statistics
  for (uint i_t = 0; i_t < 3; ++i_t) // type index
  {
    txt_out_f << "         r_b         avg   sqrt(var)         sem\n";
    txt_out_f << "    0.000000    0.000000    0.000000    0.000000\n";
    for (uint i_b = 0; i_b < n_b; ++i_b) // bin index
    {
      txt_out_f << cnfs(ng.R_n * pow((i_b + 1.0) / n_b, 1.0 / 3), 12, ' ', 6);
      rcd_o[i_t][i_b].save_last_is_stat_s(txt_out_f);
    }
    txt_out_f << "\n\n";
  }

  // save sd last individual simulation statistics
  for (uint i_c = 0; i_c < n_chr_h; ++i_c) // chromosome index
  {
    txt_out_f << "     s         avg\n";
    for (uint i_a = 0; i_a < lsdcp[i_c]; ++i_a) // array index
    {
      txt_out_f << cnfs((i_a + 1), 6, ' ');
      sd_bo[i_c][i_a].save_last_is_stat(txt_out_f);
    }
    txt_out_f << "\n\n";
  }

  // save cp last individual simulation statistics
  for (uint i_c = 0; i_c < n_chr_h; ++i_c) // chromosome index
  {
    txt_out_f << "     s         avg\n";
    for (uint i_a = 0; i_a < lsdcp[i_c]; ++i_a) // array index
    {
      txt_out_f << cnfs((i_a + 1), 6, ' ');
      cp_bo[i_c][i_a].save_last_is_stat(txt_out_f);
    }
    txt_out_f << "\n\n";
  }
}

// clear individual simulation variables
void chrana::clear_is_var()
{
  t_o.is_ts.clear();
  dcm_o.is_ts.clear();
  rg2_o.is_ts.clear();
  nop_o.is_ts.clear();
  ncf_o.is_ts.clear();
  for (uint i_t = 0; i_t < 3; ++i_t) // type index
  {
    for (uint i_b = 0; i_b < n_b; ++i_b) // bin index
    {
      rcd_o[i_t][i_b].is_ts.clear();
    }
  }
  for (uint i_c = 0; i_c < n_chr_h; ++i_c) // chromosome index
  {
    for (uint i_a = 0; i_a < lsdcp[i_c]; ++i_a) // array index
    {
      sd_bo[i_c][i_a].is_o_sum = 0.0;
      sd_bo[i_c][i_a].is_n_dp = 0;
      cp_bo[i_c][i_a].is_o_sum = 0.0;
      cp_bo[i_c][i_a].is_n_dp = 0;
    }
  }
  for (uint i_a = 0; i_a < lcm; ++i_a) // array index
  {
    cm_bo[i_a].is_o_sum = 0.0;
    cm_bo[i_a].is_n_dp = 0;
  }
}

// calculate combined simulations final statistics
void chrana::calc_cs_final_stat()
{
  dcm_o.calc_cs_final_stat();
  rg2_o.calc_cs_final_stat();
  nop_o.calc_cs_final_stat();
  ncf_o.calc_cs_final_stat();
  for (uint i_t = 0; i_t < 3; ++i_t) // type index
  {
    for (uint i_b = 0; i_b < n_b; ++i_b) // bin index
    {
      rcd_o[i_t][i_b].calc_cs_final_stat();
    }
  }
  for (uint i_c = 0; i_c < n_chr_h; ++i_c) // chromosome index
  {
    for (uint i_a = 0; i_a < lsdcp[i_c]; ++i_a) // array index
    {
      sd_bo[i_c][i_a].calc_cs_final_stat();
      cp_bo[i_c][i_a].calc_cs_final_stat();
    }
  }
  for (uint i_a = 0; i_a < lcm; ++i_a) // array index
  {
    cm_bo[i_a].calc_cs_final_stat();
  }
}

// save combined simulations final statistics
void chrana::save_cs_final_stat(std::ofstream &txt_out_f) // text output file
{
  // save dcm, rg2, nop and ncf combined simulations final statistics
  txt_out_f << "         avg   sqrt(var)         sem\n";
  dcm_o.save_cs_final_stat(txt_out_f);
  rg2_o.save_cs_final_stat(txt_out_f);
  nop_o.save_cs_final_stat(txt_out_f);
  ncf_o.save_cs_final_stat(txt_out_f);
  txt_out_f << "\n\n";

  // save rcd combined simulations final statistics
  for (uint i_t = 0; i_t < 3; ++i_t) // type index
  {
    txt_out_f << "         r_b         avg   sqrt(var)         sem\n";
    txt_out_f << "    0.000000    0.000000    0.000000    0.000000\n";
    for (uint i_b = 0; i_b < n_b; ++i_b) // bin index
    {
      txt_out_f << cnfs(ng.R_n * pow((i_b + 1.0) / n_b, 1.0 / 3), 12, ' ', 6);
      rcd_o[i_t][i_b].save_cs_final_stat(txt_out_f);
    }
    txt_out_f << "\n\n";
  }

  // save sd combined simulations final statistics
  for (uint i_c = 0; i_c < n_chr_h; ++i_c) // chromosome index
  {
    txt_out_f << "     s         avg   sqrt(var)         sem\n";
    for (uint i_a = 0; i_a < lsdcp[i_c]; ++i_a) // array index
    {
      txt_out_f << cnfs((i_a + 1), 6, ' ');
      sd_bo[i_c][i_a].save_cs_final_stat(txt_out_f);
    }
    txt_out_f << "\n\n";
  }

  // save cp combined simulations final statistics
  for (uint i_c = 0; i_c < n_chr_h; ++i_c) // chromosome index
  {
    txt_out_f << "     s         avg   sqrt(var)         sem\n";
    for (uint i_a = 0; i_a < lsdcp[i_c]; ++i_a) // array index
    {
      txt_out_f << cnfs((i_a + 1), 6, ' ');
      cp_bo[i_c][i_a].save_cs_final_stat(txt_out_f);
    }
    txt_out_f << "\n\n";
  }
}

// save contact map average values to binary file
void chrana::save_cm_avg_bin(std::ofstream &bin_out_f) // binary output file
{
  for (uint i_a = 0; i_a < lcm; ++i_a) // array index
  {
    bin_out_f.write(reinterpret_cast<char *>(&(cm_bo[i_a].cs_fs.avg)), 8);
  }
}

// calculate observables
void chrana::calc_observables()
{
  // store time
  t_o.is_ts.push_back(t);

  // calculate the center of mass position
  vec3f cmr = {0.0, 0.0, 0.0}; // center of mass position
  for (uint i_p = 0; i_p < N; ++i_p) // particle index
  {
    cmr += hr[i_p];
  }
  cmr /= N;

  // calculate the center of mass distance
  float dcm = length(cmr); // center of mass distance
  dcm_o.is_ts.push_back(dcm);

  // calculate the gyration radius squared
  float rg2 = 0.0; // gyration radius squared
  for (uint i_p = 0; i_p < N; ++i_p) // particle index
  {
    rg2 += dot(hr[i_p] - cmr, hr[i_p] - cmr);
  }
  rg2 /= N;
  rg2_o.is_ts.push_back(rg2);

  // calculate the nematic order parameter
  vec3f vec; // bond vector
  float cos; // cosine of the angle
  float nop = 0.0; // nematic order parameter
  for (uint i_p = 0; i_p < (N - 1); ++i_p) // particle index
  {
    vec = hr[i_p + 1] - hr[i_p];
    cos = dot(vec, hr[i_p]) / (length(vec) * length(hr[i_p]));
    nop += 0.5 * (3.0 * cos * cos - 1.0);
  }
  nop /= N - 1.0;
  nop_o.is_ts.push_back(nop);

  // calculate the nucleus chromatin fraction
  float ncf = 0.0; // nucleus wall pressure
  for (uint i_p = 0; i_p < N; ++i_p) // particle index
  {
    if (hr[i_p].z < ng.nod) { ncf += 1.0; }
    else { continue; }
  }
  ncf /= N;
  ncf_o.is_ts.push_back(ncf);

  // calculate the radial chromatin density
  float r_b[n_b]; // bin radius
  float vol_t; // total volume inside
  float rcd[3][n_b]; // radial chromatin density
  for (uint i_b = 0; i_b < n_b; ++i_b) // bin index
  {
    r_b[i_b] = ng.R_n * pow((i_b + 1.0) / n_b, 1.0 / 3);
  }
  for (uint i_t = 0; i_t < 3; ++i_t) // type index
  {
    for (uint i_b = 0; i_b < n_b; ++i_b) // bin index
    {
      rcd[i_t][i_b] = 0.0;
    }
  }
  for (uint i_p = 0; i_p < N; ++i_p) // particle index
  {
    float d_r = length(hr[i_p]); // radial distance
    for (uint i_b = 0; i_b < n_b; ++i_b) // bin index
    {
      if ((d_r - (rco / 2.0)) > r_b[i_b]) { vol_t = 0.0; }
      else if ((d_r + (rco / 2.0)) < r_b[i_b])
      {
        vol_t = (4.0 / 3.0) * M_PI * (rco / 2.0) * (rco / 2.0) * (rco / 2.0);
      }
      else
      {
        float cos_a = // A side cosine
            (d_r * d_r + (rco / 2.0) * (rco / 2.0) - r_b[i_b] * r_b[i_b]) /
            (2.0 * d_r * (rco / 2.0));
        float cos_b = // B side cosine
            (r_b[i_b] * r_b[i_b] + d_r * d_r - (rco / 2.0) * (rco / 2.0)) /
            (2.0 * r_b[i_b] * d_r);
        float vol_a = // A side volume
            (M_PI / 3.0) * (rco / 2.0) * (rco / 2.0) * (rco / 2.0) *
            (2.0 + cos_a) * (1.0 - cos_a) * (1.0 - cos_a);
        float vol_b = // B side volume
            (M_PI / 3.0) * r_b[i_b] * r_b[i_b] * r_b[i_b] * (2.0 + cos_b) *
            (1.0 - cos_b) * (1.0 - cos_b);
        vol_t = vol_a + vol_b;
      }
      if (hpt[i_p] == LADh) { rcd[0][i_b] += vol_t; }
      if (hpt[i_p] == LNDe) { rcd[1][i_b] += vol_t; }
      rcd[2][i_b] += vol_t;
      if (i_b != (n_b - 1))
      {
        if (hpt[i_p] == LADh) { rcd[0][i_b + 1] -= vol_t; }
        if (hpt[i_p] == LNDe) { rcd[1][i_b + 1] -= vol_t; }
        rcd[2][i_b + 1] -= vol_t;
      }
    }
  }
  for (uint i_t = 0; i_t < 3; ++i_t) // type index
  {
    for (uint i_b = 0; i_b < n_b; ++i_b) // bin index
    {
      rcd[i_t][i_b] /= (4.0 / 3.0) * M_PI * ng.R_n * ng.R_n * ng.R_n / n_b;
      rcd_o[i_t][i_b].is_ts.push_back(rcd[i_t][i_b]);
    }
  }

  // calculate the spatial distance and contact probability
  uint i_e; // ending index
  for (uint i_c = 0; i_c < n_chr_h; ++i_c) // chromosome index
  {
    if (N == N_def) { i_e = hchrla[i_c + 1]; }
    else if (i_c == 0) { i_e = N; }
    else { continue; }
    calc_sd_cp<<<(lsdcp[i_c] + tpb - 1) / tpb, tpb>>>(
        N, lsdcp[i_c], hchrla[i_c], i_e, r, sd[i_c], cp[i_c]);
    cuda_check(cudaMemcpy(
        hsd[i_c], sd[i_c], lsdcp[i_c] * sizeof(float), cudaMemcpyDeviceToHost));
    cuda_check(cudaMemcpy(
        hcp[i_c], cp[i_c], lsdcp[i_c] * sizeof(float), cudaMemcpyDeviceToHost));
    for (uint i_a = 0; i_a < lsdcp[i_c]; ++i_a) // array index
    {
      sd_bo[i_c][i_a].is_o_sum += hsd[i_c][i_a];
      ++sd_bo[i_c][i_a].is_n_dp;
      cp_bo[i_c][i_a].is_o_sum += hcp[i_c][i_a];
      ++cp_bo[i_c][i_a].is_n_dp;
    }
  }

  // calculate the contact map
  calc_cm<<<(lcm + tpb - 1) / tpb, tpb>>>(N, lcm, r, cm);
  cuda_check(cudaMemcpy(hcm, cm, lcm * sizeof(float), cudaMemcpyDeviceToHost));
  for (uint i_a = 0; i_a < lcm; ++i_a) // array index
  {
    cm_bo[i_a].is_o_sum += hcm[i_a];
    ++cm_bo[i_a].is_n_dp;
  }
}

} // namespace mmc

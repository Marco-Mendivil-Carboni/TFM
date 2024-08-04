// Includes

#include "chrana.cuh" // chromatin analysis

// Namespace

namespace mmc // Marco Mendívil Carboni
{

// Global Functions

// calculate the mean spatial distance
__global__ void calc_msd(
    const uint N, // number of particles
    uint lma, // length of msd array
    uint i_c, // chromosome index
    vec3f *r, // position array
    float *ma) // msd array
{
  // calculate array index
  int i_a = blockIdx.x * blockDim.x + threadIdx.x; // array index
  if (i_a >= lma) { return; }

  // calculate the mean spatial distance
  uint chr_len; // chromosome length
  if (N == N_def) { chr_len = chrla[i_c + 1] - chrla[i_c]; }
  else { chr_len = N; }
  ma[i_a] = 0.0;
  uint s = i_a + 1; // separation
  uint i_s = chrla[i_c]; // starting index
  for (uint i_p = i_s; i_p < i_s + (chr_len - s); ++i_p) // particle index
  {
    ma[i_a] += length(r[i_p + s] - r[i_p]);
  }
  ma[i_a] /= (chr_len - s);
}

// Host Functions

// chromatin analysis constructor
chrana::chrana(parmap &par) // parameters
    : chrdat(par)
{
  // allocate memory
  if (N == N_def)
  {
    for (uint i_c = 0; i_c < n_chr; ++i_c) // chromosome index
    {
      lma[i_c] = (hchrla[i_c + 1] - hchrla[i_c]) / 2;
      msd_o[i_c] = new simobs_b[lma[i_c]];
    }
  }
  else
  {
    lma[0] = N / 2;
    msd_o[0] = new simobs_b[lma[0]];
    for (uint i_c = 1; i_c < n_chr; ++i_c) // chromosome index
    {
      lma[i_c] = 0;
      msd_o[i_c] = new simobs_b[lma[i_c]];
    }
  }

  // allocate device memory
  for (uint i_c = 0; i_c < n_chr; ++i_c) // chromosome index
  {
    cuda_check(cudaMalloc(&ma[i_c], lma[i_c] * sizeof(float)));
  }

  // allocate host memory
  for (uint i_c = 0; i_c < n_chr; ++i_c) // chromosome index
  {
    cuda_check(cudaMallocHost(&hma[i_c], lma[i_c] * sizeof(float)));
  }
}

// chromatin analysis destructor
chrana::~chrana()
{
  // deallocate memory
  for (uint i_c = 0; i_c < n_chr; ++i_c) // chromosome index
  {
    delete[] msd_o[i_c];
  }

  // deallocate device memory
  for (uint i_c = 0; i_c < n_chr; ++i_c) // chromosome index
  {
    cuda_check(cudaFree(ma[i_c]));
  }

  // deallocate host memory
  for (uint i_c = 0; i_c < n_chr; ++i_c) // chromosome index
  {
    cuda_check(cudaFreeHost(hma[i_c]));
  }
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
  for (uint i_c = 0; i_c < n_chr; ++i_c) // chromosome index
  {
    for (uint i_a = 0; i_a < lma[i_c]; ++i_a) // array index
    {
      msd_o[i_c][i_a].calc_last_is_stat();
    }
  }
}

// save last individual simulation statistics
void chrana::save_last_is_stat(std::ofstream &txt_out_f) // text output file
{
  // save dcm, rg2, nop and ncf last individual simulation statistics
  txt_out_f << "#individual simulation analysis\n";
  txt_out_f << "#        avg   sqrt(var)         sem   f_n_b     i_t ter\n";
  txt_out_f << "# center of mass distance:\n";
  dcm_o.save_last_is_stat(txt_out_f);
  txt_out_f << "# gyration radius squared:\n";
  rg2_o.save_last_is_stat(txt_out_f);
  txt_out_f << "# nematic order parameter:\n";
  nop_o.save_last_is_stat(txt_out_f);
  txt_out_f << "# nucleus chromatin fraction:\n";
  ncf_o.save_last_is_stat(txt_out_f);
  txt_out_f << "\n\n";

  // save rcd last individual simulation statistics
  for (uint i_t = 0; i_t < 3; ++i_t) // type index
  {
    txt_out_f << "#        r_b         avg   sqrt(var)         sem\n";
    txt_out_f << "    0.000000    0.000000    0.000000    0.000000\n";
    for (uint i_b = 0; i_b < n_b; ++i_b) // bin index
    {
      txt_out_f << cnfs(ng.R_n * pow((i_b + 1.0) / n_b, 1.0 / 3), 12, ' ', 6);
      rcd_o[i_t][i_b].save_last_is_stat_s(txt_out_f);
    }
    txt_out_f << "\n\n";
  }

  // save msd last individual simulation statistics
  for (uint i_c = 0; i_c < n_chr; ++i_c) // chromosome index
  {
    txt_out_f << "#    s         avg\n";
    for (uint i_a = 0; i_a < lma[i_c]; ++i_a) // array index
    {
      txt_out_f << cnfs((i_a + 1), 6, ' ');
      msd_o[i_c][i_a].save_last_is_stat(txt_out_f);
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
  for (uint i_c = 0; i_c < n_chr; ++i_c) // chromosome index
  {
    for (uint i_a = 0; i_a < lma[i_c]; ++i_a) // array index
    {
      msd_o[i_c][i_a].is_o_sum = 0.0;
      msd_o[i_c][i_a].is_n_dp = 0;
    }
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
  for (uint i_c = 0; i_c < n_chr; ++i_c) // chromosome index
  {
    for (uint i_a = 0; i_a < lma[i_c]; ++i_a) // array index
    {
      msd_o[i_c][i_a].calc_cs_final_stat();
    }
  }
}

// save combined simulations final statistics
void chrana::save_cs_final_stat(std::ofstream &txt_out_f) // text output file
{
  // save dcm, rg2, nop and ncf combined simulations final statistics
  txt_out_f << "#combined simulations final analysis\n";
  txt_out_f << "#        avg   sqrt(var)         sem\n";
  txt_out_f << "# center of mass distance:\n";
  dcm_o.save_cs_final_stat(txt_out_f);
  txt_out_f << "# gyration radius squared:\n";
  rg2_o.save_cs_final_stat(txt_out_f);
  txt_out_f << "# nematic order parameter:\n";
  nop_o.save_cs_final_stat(txt_out_f);
  txt_out_f << "# nucleus chromatin fraction:\n";
  ncf_o.save_cs_final_stat(txt_out_f);
  txt_out_f << "\n\n";

  // save rcd combined simulations final statistics
  for (uint i_t = 0; i_t < 3; ++i_t) // type index
  {
    txt_out_f << "#        r_b         avg   sqrt(var)         sem\n";
    txt_out_f << "    0.000000    0.000000    0.000000    0.000000\n";
    for (uint i_b = 0; i_b < n_b; ++i_b) // bin index
    {
      txt_out_f << cnfs(ng.R_n * pow((i_b + 1.0) / n_b, 1.0 / 3), 12, ' ', 6);
      rcd_o[i_t][i_b].save_cs_final_stat(txt_out_f);
    }
    txt_out_f << "\n\n";
  }

  // save msd combined simulations final statistics
  for (uint i_c = 0; i_c < n_chr; ++i_c) // chromosome index
  {
    txt_out_f << "#    s         avg   sqrt(var)         sem\n";
    for (uint i_a = 0; i_a < lma[i_c]; ++i_a) // array index
    {
      txt_out_f << cnfs((i_a + 1), 6, ' ');
      msd_o[i_c][i_a].save_cs_final_stat(txt_out_f);
    }
    txt_out_f << "\n\n";
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
      if ((d_r - 0.5) > r_b[i_b]) { vol_t = 0.0; }
      else if ((d_r + 0.5) < r_b[i_b]) { vol_t = (4.0 / 3.0) * M_PI * 0.125; }
      else
      {
        float cos_a = // A side cosine
            (d_r * d_r + 0.25 - r_b[i_b] * r_b[i_b]) / (2.0 * d_r * 0.5);
        float cos_b = // B side cosine
            (r_b[i_b] * r_b[i_b] + d_r * d_r - 0.25) / (2.0 * r_b[i_b] * d_r);
        float vol_a = // A side volume
            (M_PI / 3.0) * 0.125 * (2.0 + cos_a) * (1.0 - cos_a) *
            (1.0 - cos_a);
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

  // calculate the mean spatial distance
  uint tpb = 128; // threads per block
  for (uint i_c = 0; i_c < n_chr; ++i_c) // chromosome index
  {
    calc_msd<<<(lma[i_c] + tpb - 1) / tpb, tpb>>>(N, lma[i_c], i_c, r, ma[i_c]);
    cuda_check(cudaMemcpy(
        hma[i_c], ma[i_c], lma[i_c] * sizeof(float), cudaMemcpyDeviceToHost));
    for (uint i_a = 0; i_a < lma[i_c]; ++i_a) // array index
    {
      msd_o[i_c][i_a].is_o_sum += hma[i_c][i_a];
      ++msd_o[i_c][i_a].is_n_dp;
    }
  }
}

} // namespace mmc

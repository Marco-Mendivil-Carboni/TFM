// Includes

#include "chrsim.cuh" // chromatin simulation

#include <chrono> // time library

// Functions

// main function
int main(
    const int argc, // argument count
    const char **argv) // argument vector
{
  // check command-line arguments
  if (argc < 2)
  {
    std::cout << "no arguments\n";
    return EXIT_FAILURE;
  }
  if (argc > 2)
  {
    std::cout << "extra arguments\n";
    return EXIT_FAILURE;
  }

  // using declarations
  using std::chrono::duration;
  using std::chrono::system_clock;

  // declare auxiliary variables
  std::string pathstr = argv[1]; // file path string
  std::ofstream out_f; // output file
  std::ofstream n_out_f; // null output file
  const uint n_r = 4; // number of repetitions
  const uint n_t = 8; // number of tests
  const float cvf = 0.125; // chromatin volume fraction
  const float laf = 0.25; // lbs area fraction
  uint N; // number of particles
  float R_n; // nucleus radius
  uint n_l; // number of lbs
  const uint fpf = 16; // frames per file
  system_clock::time_point stp; // starting time point
  system_clock::time_point etp; // ending time point
  duration<double> d_e; // execution duration
  double t_e[n_t] = {}; // execution times

  // main try block
  try
  {
    // open output files
    out_f.open(pathstr);
    n_out_f.open("/dev/null");
    mmc::check_file(out_f, pathstr);

    // iterate over all repetitions
    for (uint i_r = 0; i_r < n_r; ++i_r) // repetition index
    {
      // iterate over all tests
      for (uint i_t = 0; i_t < n_t; ++i_t) // test index
      {
        // set parameters
        std::stringstream par_s; // parameter stream
        N = pow(2.0, 8.0 + i_t);
        R_n = (mmc::rco / 2.0) + (mmc::rco / 2.0) * pow(N / cvf, 1.0 / 3);
        n_l = laf * 4.0 / pow(mmc::lco / (R_n - mmc::rco), 2.0);
        par_s << "number_of_particles"
              << " " << mmc::cnfs(N, 5, '0') << std::endl;
        par_s << "nucleus_radius"
              << " " << mmc::cnfs(R_n, 5, '0', 2) << std::endl;
        par_s << "number_of_lbs"
              << " " << mmc::cnfs(n_l, 5, '0') << std::endl;
        par_s << "frames_per_file"
              << " " << mmc::cnfs(fpf, 4, '0') << std::endl;
        mmc::parmap par(par_s); // parameters

        // initialize simulation and generate initial condition
        mmc::chrsim sim(par); // simulation
        sim.generate_initial_condition();

        // measure execution time
        stp = system_clock::now();
        sim.run_simulation(n_out_f);
        etp = system_clock::now();
        d_e = etp - stp;
        t_e[i_t] += 1'000.0 * d_e.count() / (fpf / mmc::dt);
      }
    }

    // write execution times to output file
    for (uint i_t = 0; i_t < n_t; ++i_t) // test index
    {
      N = pow(2.0, 8.0 + i_t);
      out_f << mmc::cnfs(N, 5, '0') << " ";
      out_f << mmc::cnfs(t_e[i_t] / n_r, 5, '0', 3) << std::endl;
    }

    // close output files
    out_f.close();
    n_out_f.close();
  }
  catch (const mmc::error &err) // caught error
  {
    // exit program unsuccessfully
    mmc::logger::record(err.what());
    return EXIT_FAILURE;
  }

  // exit program successfully
  return EXIT_SUCCESS;
}

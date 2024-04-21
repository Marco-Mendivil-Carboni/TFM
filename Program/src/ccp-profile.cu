// Includes

#include "chrsim.cuh" // chromatin simulation

#include <chrono> // time library

// Functions

// main function
int main(const int argc, // argument count
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
  uint N; // number of particles
  float R_n; // nucleus radius
  uint n_l; // number of lbs
  const uint fpf = 10; // frames per file
  system_clock::time_point stp; // starting time point
  system_clock::time_point etp; // ending time point
  duration<double> d_e; // execution duration
  double t_e; // execution time

  // main try block
  try
  {
    // open output files
    out_f.open(pathstr);
    n_out_f.open("/dev/null");
    mmc::check_file(out_f, pathstr);

    // iterate over all tests
    for (uint i_t = 0; i_t < 8; ++i_t) // test index
    {
      // set parameters
      std::stringstream par_s; // parameter stream
      N = pow(2, 8 + i_t);
      R_n = 0.5 + 0.5 * pow(N / 0.2, 1.0 / 3);
      n_l = 0.2 * 4.0 / pow(mmc::lco / (R_n - mmc::rco), 2.0);
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
      t_e = 1'000.0 * d_e.count() / (2048 * fpf);

      // write execution time to output file
      out_f << mmc::cnfs(N, 5, '0') << " " << mmc::cnfs(t_e, 5, '0', 3)
            << std::endl;
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

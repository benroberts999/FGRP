#include "FGRadPot.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

/*
This is a simple example program that uses FGRadPot.hpp.
Use as an example, of just use the output .txt file
*/

//******************************************************************************
int main(int argc, char *argv[]) {

  std::cout << "Example program.\nWill write electric and magnetic parts of "
               "radiative potential to a text file\n";

  if (argc < 2) {
    std::cout << "Run like: \n"
              << "./example Z r_rms l\n"
              << "where:\nZ is nuclear charge number\nr_rms is "
                 "root-mean-square nuclear charge radius, in "
                 "fempto meters.\nl is orbital angular quantum number\n"
              << "e.g.,:\n"
              << "./example 55 4.8041 0\n";
    return 1;
  }

  // Read in atomic charge number:
  const double z = argc > 1 ? std::stod(argv[1]) : 1.0;
  // Rad nuclear rms radius (in fempto meters) - assumed 0.0
  const double rrms_fm = argc > 2 ? std::stod(argv[2]) : 0.0;
  // Read in l (orbital angular quantum number)
  // Note: Only fitting factors A,B depend on l - assumed l=0 by default
  const int l = argc > 3 ? std::stoi(argv[3]) : 0;

  // Check nuclear radius, convert to atomic units
  // CODATA 2018: Bohr Radius a_B = 0.529177210903(80)e-10 m
  constexpr double aB_fm = 0.529177210903e5;
  const auto rN = std::sqrt(5.0 / 3.0) * rrms_fm / aB_fm;
  if (rrms_fm != 0.0 && (rrms_fm < 0.1 || rrms_fm > 10.0)) {
    // probably a units error
    std::cout << "We want rrms in fm - double check!\n";
  }

  std::cout << "\nRunning for Z=" << z << ", l=" << l << "\n";
  std::cout << "r_rms = " << rrms_fm << " fm; => rN = " << rN << " au\n";

  // Set up radial grid: in this example, from r0 to rmax in num_steps,
  // logarithmically
  const double r0 = 1.0e-5;
  const double rmax = 10.0;
  const int num_steps = 250;
  const double delta = std::log(rmax / r0) / double(num_steps - 1);
  std::vector<double> r;
  for (int i = 0; i < num_steps; ++i) {
    r.push_back(r0 * std::exp(i * delta));
  }

  // Lambda function that fills a std::vector with values of a function f,
  // evaluated at positions stored in vector r
  const auto fill = [z, rN](auto function, const std::vector<double> &r) {
    std::vector<double> v;
    v.reserve(r.size());
    for (auto &ri : r) {
      v.push_back(function(z, ri, rN));
    }
    return v;
  };

  // Fill std::vector with values of radial potential along radial grid:
  const auto V_u = fill(FGRP::V_Uehling, r);
  const auto V_h = fill(FGRP::V_SEh, r);
  const auto V_l = fill(FGRP::V_SEl, r);
  const auto H_m = fill(FGRP::H_Magnetic, r);

  // Write out to text file, called radpot.txt
  std::ofstream out_file("radpot.txt");
  out_file << "#r V_el(r) H_mag(r)\n";
  for (auto i = 0ul; i < r.size(); ++i) {
    const auto vel = V_u.at(i) + FGRP::Fit::A(z, l) * V_h.at(i) +
                     FGRP::Fit::B(z, l) * V_l.at(i);
    out_file << r.at(i) << " " << vel << " " << H_m.at(i) << "\n";
  }
}

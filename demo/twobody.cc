// Error evaluation program
#include "../solver4th.h"
#include "cs.h"
#include <sysexits.h>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <string>
#include <thread>
#include <valarray>
#include <vector>

using namespace std;

void test(double rep, double alpha, double t_stop, double out_dt, const string &filename) {
  const int exp_dt_min = -5;
  const int exp_dt_max = -5;
  const double eps = 0.0;
  const double eta_s = 0.01;
  const double eta = 0.03;

  valarray<double> mass {1.0, 1.0e-3};

  // Calculate initial cartesian coordinates
  cs::OrbitalElement oelem0(1.0, 0.1, M_PI, 0.0, 0.0, M_PI);
  vector<cs::OrbitalElement> init_oelems {oelem0};
  vector<valarray<double>> init_pos(mass.size(), valarray<double>(cs::DIM));
  vector<valarray<double>> init_vel(mass.size(), valarray<double>(cs::DIM));
  cs::ConvertOelemToCartesian(mass, init_oelems, &init_pos, &init_vel);

  // Initialize solver
  hermite::Hermite4thSolver h4sol(mass.size(), rep, exp_dt_min, exp_dt_max, eps, alpha, eta_s, eta);
  h4sol.set_mass(mass);
  h4sol.set_pos(init_pos);
  h4sol.set_vel(init_vel);
  h4sol.setup();

  ofstream ofs(filename);
  if (!ofs) {
    cerr << "File Open Error: \"" << filename << "\"" << endl;
    exit(EX_CANTCREAT);
  }

  double t = 0.0;
  double out_t = 0.0;
  auto pos = h4sol.get_pos();
  auto vel = h4sol.get_vel();
  vector<cs::OrbitalElement> oelems(mass.size() -1);
  cs::WriteOelemHeader(oelems.size(), &ofs);

  while (t < t_stop) {
    if (t >= out_t) {
      oelems = cs::ConvertCartesianToOelem(mass, pos, vel);
      cs::WriteOelemData(t, oelems, &ofs);
      out_t += out_dt;
    }
    h4sol.step_forward(&t);
    pos = h4sol.get_pos();
    vel = h4sol.get_vel();
  }

}


int main() {
  vector<thread> threads;
  threads.emplace_back( []() { test(1, 1.0, 2.0*M_PI*10, std::pow(2.0, -5), "./oel_2b_rep1_alpha1_short.csv"); } );
  threads.emplace_back( []() { test(2, 1.0, 2.0*M_PI*10, std::pow(2.0, -5), "./oel_2b_rep2_alpha1_short.csv"); } );
  threads.emplace_back( []() { test(3, 1.0, 2.0*M_PI*10, std::pow(2.0, -5), "./oel_2b_rep3_alpha1_short.csv"); } );
  threads.emplace_back( []() { test(4, 1.0, 2.0*M_PI*10, std::pow(2.0, -5), "./oel_2b_rep4_alpha1_short.csv"); } );
  threads.emplace_back( []() { test(2, 1.0, 2.0*M_PI*1e5, std::pow(2.0, 7), "./oel_2b_rep2_alpha1_long.csv"); } );
  threads.emplace_back( []() { test(3, 1.0, 2.0*M_PI*1e5, std::pow(2.0, 7), "./oel_2b_rep3_alpha1_long.csv"); } );
  threads.emplace_back( []() { test(4, 1.0, 2.0*M_PI*1e5, std::pow(2.0, 7), "./oel_2b_rep4_alpha1_long.csv"); } );
  threads.emplace_back( []() { test(3, 5.0/6.0, 2.0*M_PI*10, std::pow(2.0, -5), "./oel_2b_rep3_alpha083_short.csv"); } );
  threads.emplace_back( []() { test(3, 7.0/6.0, 2.0*M_PI*10, std::pow(2.0, -5), "./oel_2b_rep3_alpha116_short.csv"); } );
  threads.emplace_back( []() { test(3, 7.0/6.0, 2.0*M_PI*1e6, std::pow(2.0, 9), "./oel_2b_rep3_alpha116_longlong.csv"); } );
  for(auto& t : threads) t.join();

  return 0;
}
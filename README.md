# hermite

Hermite integrator for planetary N-body simulation


## Usage


### Test

``` sh
$ cd {PROJECT_ROOT}/build
$ cmake ..
$ make
$ ctest -V
```


### Demo

``` sh
$ cd {PROJECT_ROOT}/build
$ cmake ..
$ make
$ cd demo
$ ./twobody
$ jupyter-lab --notebook-dir .
```

### Example

``` cpp
#include "hermite.h"
#include "cs.h"
#include <sysexits.h>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <valarray>
#include <vector>

// Parameters
std::size_t rep = 3;
int log2_dt_max = -5;
int log2_dt_min = -10;
double eps = 0.0;
double alpha = 1.0;
double eta_s = 0.01;
double eta = 0.03;

// Calculate initial position and velocity from orbital elements
std::valarray<double> mass = {1.0, 1e-3};
cs::OrbitalElement oelem0(1.0, 0.1, M_PI, 0.0, 0.0, M_PI);
vector<cs::OrbitalElement> init_oelems {oelem0};
vector<valarray<double>> pos_init(mass.size(), valarray<double>(cs::DIM));
vector<valarray<double>> vel_init(mass.size(), valarray<double>(cs::DIM));
cs::ConvertOelemToCartesian(mass, init_oelems, &init_pos, &init_vel);

// Set solver
hermite::Hermite4thSolver h4sol(rep, exp_dt_max, exp_dt_min, eps, alpha, eta_s, eta);
h4sol.set_mass(mass);
h4sol.set_pos(init_pos);
h4sol.set_vel(init_vel);
h4sol.setup();

// Set output file
ofstream ofs(filename);
if (!ofs) {
  cerr << "File Open Error: \"" << filename << "\"" << endl;
  exit(EX_CANTCREAT);
}

// Start orbital calculation
double t = 0.0;
double t_stop = 10.0;
auto pos = h4sol.get_pos();
auto vel = h4sol.get_vel();

while (t < t_stop) {
  cs::WritePosVel(t, pos, vel, &ofs);
  h4sol.step_forward(&t);
  pos = h4sol.get_pos();
  vel = h4sol.get_vel();
}
```


## References

- Kokubo, E., Makino, J., 2004, PASJ, 56, 861
- Kokubo, E., Yoshinaga, K., Makino, J., 1998, MNRAS, 297, 1067

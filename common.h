//
// Created by aoshima on 18/10/26.
//

#ifndef HERMITE_COMMON_H
#define HERMITE_COMMON_H
#include <cstddef>


namespace hermite {

  using namespace std;
  using FLOAT = double;
  constexpr size_t DIM = 3;
  constexpr FLOAT G = 1.0;

  bool CHECK_NEAR(FLOAT val1, FLOAT val2, FLOAT abs_error)
  {
    return (val1 - val2) * (val1 - val2) < (abs_error * abs_error);
  }

};  // namespace hermite

#endif //HERMITE_COMMON_H

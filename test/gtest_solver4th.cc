#include "../solver4th.h"
#include "gtest/gtest.h"
#include "cs.h"
#include <vector>
#include <valarray>
#include <cmath>
#include <fstream>


namespace hermite {


class Solver4thTwoBodyTest : public :: testing::Test {
 protected:

  static void SetUpTestCase()
  {
    p_init_oelems = new vector<cs::OrbitalElement>;
    p_last_oelems = new vector<cs::OrbitalElement>;
    valarray<FLOAT> mass {1.0, 1.0e-10};
    // set initial orbital elements
    cs::OrbitalElement init_oelem0(1.0, 0.01, 3.0, 1.0, 0.5, 1.0);
    (*p_init_oelems).push_back(init_oelem0);
    // calculate initial position and velocity
    vector<valarray<FLOAT>> pos_init(mass.size(), valarray<FLOAT>(cs::DIM));
    vector<valarray<FLOAT>> vel_init(mass.size(), valarray<FLOAT>(cs::DIM));
    cs::ConvertOelemToCartesian(mass, *p_init_oelems, &pos_init, &vel_init);
    // reset initial orbital elements
    *p_init_oelems = cs::ConvertCartesianToOelem(mass, pos_init, vel_init);
    // calculate 100 orbital period
    Hermite4thSolver solver(2, 3, -8, -5, 0.0, 1.0, 0.01, 0.03);
    solver.set_mass(mass);
    solver.set_pos(pos_init);
    solver.set_vel(vel_init);
    solver.setup();
    FLOAT t = 0.0;
    FLOAT t_stop = 2.0 * M_PI * 100;
    while (t < t_stop) solver.step_forward(&t);
    // Get last position and velocity
    auto last_pos = solver.get_pos();
    auto last_vel = solver.get_vel();
    // Get last orbital elements
    *p_last_oelems = cs::ConvertCartesianToOelem(mass, last_pos, last_vel);
  }

  static void TearDownTestCase() {
    delete p_init_oelems;
    delete p_last_oelems;
    p_init_oelems = nullptr;
    p_last_oelems = nullptr;
  }

  static vector<cs::OrbitalElement> *p_init_oelems;
  static vector<cs::OrbitalElement> *p_last_oelems;

};

vector<cs::OrbitalElement> *Solver4thTwoBodyTest::p_init_oelems = nullptr;
vector<cs::OrbitalElement> *Solver4thTwoBodyTest::p_last_oelems = nullptr;


TEST_F(Solver4thTwoBodyTest, ConservationOfSemimajorAxis) {
  EXPECT_NEAR((*p_init_oelems)[0].sa, (*p_last_oelems)[0].sa, 1e-9);
}

TEST_F(Solver4thTwoBodyTest, ConservationOfOrbitalEccentricity) {
  EXPECT_NEAR((*p_init_oelems)[0].oe, (*p_last_oelems)[0].oe, 1e-10);
}

TEST_F(Solver4thTwoBodyTest, ConservationOfArgumentOfPeriapsis) {
  EXPECT_NEAR((*p_init_oelems)[0].peri, (*p_last_oelems)[0].peri, 1e-3);
}

TEST_F(Solver4thTwoBodyTest, ConservationOfInclination) {
  EXPECT_NEAR((*p_init_oelems)[0].incl, (*p_last_oelems)[0].incl, 1e-10);
}

TEST_F(Solver4thTwoBodyTest, ConservationOfLongitudeOfTheAscendingNode) {
  EXPECT_NEAR((*p_init_oelems)[0].node, (*p_last_oelems)[0].node, 1e-4);
}

TEST_F(Solver4thTwoBodyTest, ConservationOfMeanAnomaly) {
 EXPECT_NEAR((*p_init_oelems)[0].l, (*p_last_oelems)[0].l, 1e-1);
}


class Solver4thFiveBodyTest : public :: testing::Test {
 protected:

  static void SetUpTestCase()
  {
    p_init_ene = new FLOAT;
    p_last_ene = new FLOAT;
    p_init_h = new valarray<FLOAT>;
    p_last_h = new valarray<FLOAT>;
    valarray<FLOAT> mass {1.0, 1.0e-10, 1.0e-10, 1.0e-10, 1.0e-10};
    // set initial orbital elements
    cs::OrbitalElement init_oelem0(1.0, 0.01, 0.1, 0.0, 0.0, 0.0);
    cs::OrbitalElement init_oelem1(2.0, 0.03, 0.5, 0.5, 0.5, 0.5);
    cs::OrbitalElement init_oelem2(3.0, 0.04, 1.0, 1.0, 1.0, 1.0);
    cs::OrbitalElement init_oelem3(4.0, 0.05, 1.5, 2.0, 1.5, 1.5);
    vector<cs::OrbitalElement> init_oelems {init_oelem0, init_oelem1, init_oelem2, init_oelem3};
    // calculate initial position and velocity
    vector<valarray<FLOAT>> init_pos(mass.size(), valarray<FLOAT>(cs::DIM));
    vector<valarray<FLOAT>> init_vel(mass.size(), valarray<FLOAT>(cs::DIM));
    cs::ConvertOelemToCartesian(mass, init_oelems, &init_pos, &init_vel);
    // calculate initial energy and angular momentum
    *p_init_ene = cs::CALC_ENERGY(mass, init_pos, init_vel);
    *p_init_h = cs::CALC_ANGULAR_MOMENTUM(mass, init_pos, init_vel);
    // calculate 100 orbital period
    Hermite4thSolver solver(5, 3, -8, -5, 0.0, 1.0, 0.01, 0.03);
    solver.set_mass(mass);
    solver.set_pos(init_pos);
    solver.set_vel(init_vel);
    solver.setup();
    FLOAT t = 0.0;
    FLOAT t_stop = 2.0 * M_PI * 100;
    while (t < t_stop) solver.step_forward(&t);
    // Get last position and velocity
    auto last_pos = solver.get_pos();
    auto last_vel = solver.get_vel();
    // calc last energy
    *p_last_ene = cs::CALC_ENERGY(mass, last_pos, last_vel);
    *p_last_h = cs::CALC_ANGULAR_MOMENTUM(mass, last_pos, last_vel);
  }


  static void TearDownTestCase() {
    delete p_init_ene;
    delete p_last_ene;
    delete p_init_h;
    delete p_last_h;
    p_init_ene = nullptr;
    p_last_ene = nullptr;
    p_init_h = nullptr;
    p_last_h = nullptr;
  }

  static FLOAT* p_init_ene;
  static FLOAT* p_last_ene;
  static valarray<FLOAT>* p_init_h;
  static valarray<FLOAT>* p_last_h;

};

FLOAT* Solver4thFiveBodyTest::p_init_ene = nullptr;
FLOAT* Solver4thFiveBodyTest::p_last_ene = nullptr;
valarray<FLOAT>* Solver4thFiveBodyTest::p_init_h = nullptr;
valarray<FLOAT>* Solver4thFiveBodyTest::p_last_h = nullptr;


TEST_F(Solver4thFiveBodyTest, ConservationOfEnergy) {
  EXPECT_NEAR(*p_init_ene, *p_last_ene, 1e-10);
}


TEST_F(Solver4thFiveBodyTest, ConservationOfAngularMomentum) {
  for (size_t k = 0; k < DIM; k++) EXPECT_NEAR((*p_last_h)[k], (*p_init_h)[k], 1e-10);
}


}  // namespace hermite


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
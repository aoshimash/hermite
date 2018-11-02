// This file implements 4th order modified hermite integrator
#ifndef HERMITE_H
#define HERMITE_H
#include "common.h"
#include <valarray>
#include <algorithm>
#include <cstddef>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <numeric>


namespace hermite {

class Hermite4thSolver {

 public:

  Hermite4thSolver(size_t num,
                   size_t rep,
                   int log2_dt_min,
                   int log2_dt_max,
                   FLOAT eps,
                   FLOAT alpha,
                   FLOAT eta_s,
                   FLOAT eta);

  // Setup variables to start orbital calculation
  void setup();

  // Calculate orbits until "t" reaches "t + dt_max"
  void step_forward(FLOAT *t);

  ////////////
  // Setters
  ////////////
  template <class V>
  void set_mass(const V &mass);
  template <class V2>
  void set_pos(const V2 &pos);
  template <class V2>
  void set_vel(const V2 &vel);

  ////////////
  // Getters
  ////////////
  vector<valarray<FLOAT>> get_pos() const { return pos_; }
  vector<valarray<FLOAT>> get_vel() const { return vel_; }
  vector<valarray<FLOAT>> get_acc() const { return acc_; }
  vector<valarray<FLOAT>> get_dacc() const { return dacc_; }

 protected:

  // Take one step only target particles
  void take_one_step();

  // Calculate all particles accelerations and time derived accelerations
  void calc_acc_and_dacc_all(const vector<valarray<FLOAT>> &pos,
                             const vector<valarray<FLOAT>> &vel,
                             vector<valarray<FLOAT>> *p_acc,
                             vector<valarray<FLOAT>> *p_dacc);

  // Calculate second-order differentiation of acceleration of the particle "i"
  void calc_d2acc(size_t i);

  // Calculate third-order differentiation of acceleration of the particle "i"
  void calc_d3acc(size_t i);

  // Calculate "next_t" of the particle "i" for the initial step
  void calc_next_t_for_the_first_step(size_t i);

  // Calculate "next_t" of the particle "i"
  void calc_next_t(size_t i);

  // Calculate predictor of the particle "i"
  void calc_predictor(size_t i);

  // Calculate corrector of the particle "i"
  void calc_corrector(size_t i);

  //////////////
  // Variables
  //////////////
  valarray<FLOAT> mass_ {};               // Mass of particles
  vector<valarray<FLOAT>> pos_ {};        // Position
  vector<valarray<FLOAT>> vel_ {};        // Velocity
  vector<valarray<FLOAT>> acc_ {};        // Acceleration
  vector<valarray<FLOAT>> dacc_ {};       // First-order differentiation of acceleration
  vector<valarray<FLOAT>> d2acc_ {};      // Second-order differentiation of acceleration
  vector<valarray<FLOAT>> d3acc_ {};      // Third-order differentiation of acceleration
  vector<valarray<FLOAT>> pos_next_ {};   // Predictor or corrector position
  vector<valarray<FLOAT>> vel_next_ {};   // Predictor or corrector velocity
  vector<valarray<FLOAT>> acc_next_ {};   // Predictor or corrector acceleration
  vector<valarray<FLOAT>> dacc_next_ {};  // Predictor or corrector first-order differentiation of acceleration
  vector<valarray<FLOAT>> pos_pred_ {};   // Position of predictor
  vector<valarray<FLOAT>> vel_pred_ {};   // Velocity of predictor
  valarray<FLOAT> t_ {};                  // Time of each particles
  valarray<FLOAT> t_next_ {};             // Next time of each particles
  valarray<int> log2_dt_cand_ {};         // Quantized time-step of each particles (log2(cand_dt_))
  valarray<FLOAT> dt_cand_;               // Candidate of the time-step
  valarray<FLOAT> dt_ {};                 // time-step (global-time - t_[i])
  valarray<FLOAT> dt2_ {};                // The second power of time-step
  valarray<FLOAT> dt3_ {};                // The third power of time-step
  valarray<FLOAT> dt4_ {};                // The fourth power of time-step
  valarray<FLOAT> dt5_ {};                // The fifth power of time-step

  //////////////
  // Constants
  //////////////
  const size_t num_;         // Number of particles
  const size_t rep_;         // Times of the EC iteration
  const int log2_dt_min_;    // Minimum time-step which quantized to integer powers of two
  const int log2_dt_max_;    // Maximum time-step which quantized to integer powers of two
  const FLOAT dt_max_;       // Maximum time-step
  const FLOAT eps2_;         // Square of softening parameter
  const FLOAT alpha_;        // Newmark parameter
  const FLOAT eta_s_;        // Accuracy controller required for calculating the candidate of the initial time-step
  const FLOAT eta_;          // Accuracy controller required for calculating the candidate of the new time-step
};


Hermite4thSolver::Hermite4thSolver(size_t num,
                                   size_t rep,
                                   int log2_dt_min,
                                   int log2_dt_max,
                                   FLOAT eps,
                                   FLOAT alpha,
                                   FLOAT eta_s,
                                   FLOAT eta)
        : num_(num), rep_(rep), log2_dt_min_(log2_dt_min), log2_dt_max_(log2_dt_max),
          dt_max_(pow(2.0, log2_dt_max_)), eps2_(eps*eps), alpha_(alpha), eta_s_(eta_s), eta_(eta)
{
  vector<valarray<FLOAT>> init_2darr(num_, valarray<double>(DIM));
  pos_ = init_2darr;
  vel_ = init_2darr;
  acc_ = init_2darr;
  dacc_ = init_2darr;
  d2acc_  = init_2darr;
  d3acc_ = init_2darr;
  pos_next_ = init_2darr;
  vel_next_ = init_2darr;
  acc_next_ = init_2darr;
  dacc_next_ = init_2darr;
  pos_pred_ = init_2darr;
  vel_pred_ = init_2darr;
  mass_.resize(num_);
  t_.resize(num_);
  dt_.resize(num_);
  dt2_.resize(num_);
  dt3_.resize(num_);
  dt4_.resize(num_);
  dt5_.resize(num_);
  t_next_.resize(num_);
  log2_dt_cand_.resize(num_);
  dt_cand_.resize(num_);
}


// Setup variables necessary to start orbital calculation
void Hermite4thSolver::setup()
{
  calc_acc_and_dacc_all(pos_, vel_, &acc_, &dacc_);
  for (size_t i = 0; i < num_; i++) calc_next_t_for_the_first_step(i);
}


template <class V>
void Hermite4thSolver::set_mass(const V &mass)
{
  for (size_t i = 0; i < num_; i++) mass_[i] = mass[i];
}


template <class V2>
void Hermite4thSolver::set_pos(const V2 &pos)
{
  for (size_t i = 0; i < num_; i++) {
    for (size_t k = 0; k < DIM; k++) {
      pos_[i][k] = pos[i][k];
    }
  }
}


template <class V2>
void Hermite4thSolver::set_vel(const V2 &vel)
{
  for (size_t i = 0; i < num_; i++) {
    for (size_t k = 0; k < DIM; k++) {
      vel_[i][k] = vel[i][k];
    }
  }
}


// Take one step only target particles
void Hermite4thSolver::take_one_step()
{
  // Set global time
  FLOAT t_global = t_next_.min();

  // Calculate time-step and predictor
  for (size_t i = 0; i < num_; i++) {
    dt_[i] = t_global - t_[i];
    dt2_[i] = dt_[i] * dt_[i];
    dt3_[i] = dt2_[i] * dt_[i];
    if (CHECK_NEAR(t_next_[i], t_global, 1e-15)) {
      dt4_[i] = dt3_[i] * dt_[i];
      dt5_[i] = dt4_[i] * dt_[i];
    }
    calc_predictor(i);
  }

  // Calculate corrector
  for (size_t n = 0; n < rep_; n++) {
    calc_acc_and_dacc_all(pos_next_, vel_next_, &acc_next_, &dacc_next_);
    for (size_t i = 0; i < num_; i++) {
      calc_d2acc(i);
      calc_d3acc(i);
      if (CHECK_NEAR(t_next_[i], t_global, 1e-15)) calc_corrector(i);
    }
  }

  // Update times, positions, velocities, accelerations and time derived accelerations
  for (std::size_t i = 0; i < num_; i++) {
    if (CHECK_NEAR(t_next_[i], t_global, 1e-15)) {
      t_[i] = t_global;
      calc_next_t(i);
      swap(pos_[i], pos_next_[i]);
      swap(vel_[i], vel_next_[i]);
      swap(acc_[i], acc_next_[i]);
      swap(dacc_[i], dacc_next_[i]);
    }
  }
}


// Calculate orbits until "t" reaches "t + dt_max"
void Hermite4thSolver::step_forward(FLOAT *t)
{
  while(t_next_[0] <= *t + dt_max_) take_one_step();
  *t += dt_max_;
}


// Calculate all particles accelerations and time derived accelerations
void Hermite4thSolver::calc_acc_and_dacc_all(const vector<valarray<FLOAT>> &pos,
                                             const vector<valarray<FLOAT>> &vel,
                                             vector<valarray<FLOAT>> *p_acc,
                                             vector<valarray<FLOAT>> *p_dacc)
{
  FLOAT r2;
  FLOAT dot_prod;
  FLOAT inv_r3;
  FLOAT inv_r5;
  valarray<FLOAT> part1(DIM);
  valarray<FLOAT> part2(DIM);
  valarray<FLOAT> diff_pos(DIM);
  valarray<FLOAT> diff_vel(DIM);
  // Initialize acceleration and differential coefficient of acceleration
  for (size_t i = 0; i < num_; i++) {
    (*p_acc)[i] = 0.0;
    (*p_dacc)[i] = 0.0;
  }
  // Calculate acceleration and differential coefficient of acceleration
  for (size_t i = 0; i < num_ - 1; i++) {
    for (size_t j = i + 1; j < num_; j++) {
      diff_pos = pos[j] - pos[i];
      diff_vel = vel[j] - vel[i];
      r2 = (diff_pos * diff_pos).sum() + eps2_;
      dot_prod = (diff_pos * diff_vel).sum();
      inv_r3 = 1.0 / (sqrt(r2) * r2);
      inv_r5 = inv_r3 / r2;
      part1 = G * inv_r3 * diff_pos;
      part2 = G * (diff_vel * inv_r3 - 3.0 * dot_prod * inv_r5 * diff_pos);
      (*p_acc)[i] += mass_[j] * part1;
      (*p_acc)[j] -= mass_[i] * part1;
      (*p_dacc)[i] += mass_[j] * part2;
      (*p_dacc)[j] -= mass_[i] * part2;
    }
  }
}


// Calculate second-order differentiation of acceleration of the particle "i"
void Hermite4thSolver::calc_d2acc(size_t i)
{
  d2acc_[i] = (-6.0 * (acc_[i] - acc_next_[i]) - dt_[i] * (4.0 * dacc_[i] + 2.0 * dacc_next_[i])) / dt2_[i];
}


// Calculate third-order differentiation of acceleration of the particle "i"
void Hermite4thSolver::calc_d3acc(size_t i)
{
  d3acc_[i] = (12.0 * (acc_[i] - acc_next_[i]) + 6.0 * dt_[i] * (dacc_[i] + dacc_next_[i])) / dt3_[i];
}


// Calculate "next_t" of the particle "i" for the initial step
void Hermite4thSolver::calc_next_t_for_the_first_step(size_t i) {
  // Calculate the candidate of time-step
  dt_cand_[i] = eta_s_ * sqrt((acc_[i] * acc_[i]).sum() / (dacc_[i] * dacc_[i]).sum());
  // Quantize the candidate of time-step to integer powers of two
  log2_dt_cand_[i] = static_cast<int>(floor(log2(dt_cand_[i])));
  log2_dt_cand_[i] = min(log2_dt_cand_[i], log2_dt_max_);
  log2_dt_cand_[i] = max(log2_dt_cand_[i], log2_dt_min_);
  dt_cand_[i] = pow(2.0, log2_dt_cand_[i]);
  // Calculate the time after the particle take next step
  t_next_[i] = t_[i] + dt_cand_[i];
}


// Calculate "next_t" of the particle "i"
void Hermite4thSolver::calc_next_t(size_t i)
{
  // Calculate the candidate of time-step
  valarray<FLOAT> next_d2acc_ = d2acc_[i] + d3acc_[i] * dt_[i];
  FLOAT acc_abs_next = sqrt((acc_next_[i] * acc_next_[i]).sum());
  FLOAT dacc_abs_next = sqrt((dacc_next_[i] * dacc_next_[i]).sum());
  FLOAT d2acc_abs_next = sqrt((next_d2acc_ * next_d2acc_).sum());
  FLOAT d3acc_abs = sqrt((d3acc_[i] * d3acc_[i]).sum());
  dt_cand_[i] = sqrt(eta_ * (acc_abs_next * d2acc_abs_next + dacc_abs_next * dacc_abs_next)
                     / (dacc_abs_next * d3acc_abs + d2acc_abs_next * d2acc_abs_next ));
  // Quantize the candidate of time-step to integer powers of two
  if (dt_cand_[i] < dt_[i]) {
    log2_dt_cand_[i] = max(log2_dt_cand_[i] - 1, log2_dt_min_);
  } else if (dt_cand_[i] > 2.0 * dt_[i] && CHECK_NEAR(fmod(t_[i], 2.0 * dt_[i]), 0.0, 1e-15)) {
    log2_dt_cand_[i] = min(log2_dt_cand_[i] + 1, log2_dt_max_);
  }
  dt_cand_[i] = pow(2.0, log2_dt_cand_[i]);
  // Calculate the time after the particle take next step
  t_next_[i] = t_[i] + dt_cand_[i];
}


// Calculate predictor of the particle "i"
void Hermite4thSolver::calc_predictor(size_t i)
{
  pos_next_[i] = pos_[i] + vel_[i] * dt_[i] + 0.5 * acc_[i] * dt2_[i] + (dacc_[i] / 6.0) * dt3_[i];
  vel_next_[i] = vel_[i] + acc_[i] * dt_[i] + 0.5 * dacc_[i] * dt2_[i];
  pos_pred_[i] = pos_next_[i];
  vel_pred_[i] = vel_next_[i];
}


// Calculate corrector of the particle "i"
void Hermite4thSolver::calc_corrector(size_t i)
{
  pos_next_[i] = pos_pred_[i] + (d2acc_[i] / 24.0) * dt4_[i] + (alpha_ * d3acc_[i] / 120.0) * dt5_[i];
  vel_next_[i] = vel_pred_[i] + (d2acc_[i] / 6.0) * dt3_[i] + (d3acc_[i] / 24.0) * dt4_[i];
}


}  // namespace hermite

#endif
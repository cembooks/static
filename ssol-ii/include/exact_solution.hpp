/******************************************************************************
 * Copyright (C) Siarhei Uzunbajakau, 2023.
 *
 * This program is free software. You can use, modify, and redistribute it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 or (at your option) any later version.
 * This program is distributed without any warranty.
 *
 * Refer to COPYING.LESSER for more details.
 ******************************************************************************/

#ifndef ExactSolutionsSSOLII_H__
#define ExactSolutionsSSOLII_H__

#include "constants.hpp"
#include "settings.hpp"
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>

#include <cmath>

using namespace dealii;

inline Tensor<1, 3>
volume_free_current_density(double x,
                            double y,
                            double z,
                            double K_0,
                            double a,
                            double b)
{
  Tensor<1, 3> J;
  double r;

  r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));

  if ((r > a) && (r < b)) {
    J[0] = -K_0 * y;
    J[1] = K_0 * x;
    J[2] = 0.0;
  } else {
    J[0] = 0.0;
    J[1] = 0.0;
    J[2] = 0.0;
  }

  return J;
}

inline Tensor<1, 3>
magnetic_field(double x,
               double y,
               double z,
               double K_0,
               double mu_0,
               double a,
               double b)
{
  double cos_theta;
  double sin_theta;

  double cos_phi;
  double sin_phi;

  Tensor<1, 3> r_hat;
  Tensor<1, 3> theta_hat;

  Tensor<1, 3> F1;
  Tensor<1, 3> F2;
  Tensor<1, 3> B;

  double r;

  r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));

  cos_theta = z / r;
  sin_theta = sqrt(pow(x, 2) + pow(y, 2)) / r;

  cos_phi = x / sqrt(pow(x, 2) + pow(y, 2));
  sin_phi = y / sqrt(pow(x, 2) + pow(y, 2));

  r_hat[0] = x / r;
  r_hat[1] = y / r;
  r_hat[2] = z / r;

  theta_hat[0] = cos_theta * cos_phi;
  theta_hat[1] = cos_theta * sin_phi;
  theta_hat[2] = -sin_theta;

  F1 = (2.0 / 3.0) * mu_0 * K_0 * (cos_theta * r_hat - sin_theta * theta_hat);
  F2 = (2.0 / 3.0) * mu_0 * K_0 *
       (cos_theta * r_hat + 0.5 * sin_theta * theta_hat) / pow(r, 3);

  if (r < a) {
    B = 0.5 * (pow(b, 2) - pow(a, 2)) * F1;
  } else if (r > b) {
    B = 0.2 * (pow(b, 5) - pow(a, 5)) * F2;
  } else {
    B = 0.5 * (pow(b, 2) - pow(r, 2)) * F1 + 0.2 * (pow(r, 5) - pow(a, 5)) * F2;
  }

  return B;
}

/**
 * \brief Describes the exact solution, \f$\vec{J}_f\f$, of the
 * *Thick spherical coil* [(ssol-ii/)](@ref page_ssol_ii) numerical experiment.
 *****************************************************************************/
class ExactSolutionSSOLII_Jf
  : public Function<3>
  , public SettingsSSOLII
{
public:
  ExactSolutionSSOLII_Jf();

  virtual void vector_value_list(
    const std::vector<Point<3>>& r,
    std::vector<Vector<double>>& values) const override final;
};

/**
 * \brief Describes the exact solution, \f$\vec{B}\f$, of the
 * *Thick spherical coil* [(ssol-ii/)](@ref page_ssol_ii) numerical experiment.
 *****************************************************************************/
class ExactSolutionSSOLII_B
  : public Function<3>
  , public SettingsSSOLII
{
public:
  ExactSolutionSSOLII_B();

  virtual void vector_value_list(
    const std::vector<Point<3>>& r,
    std::vector<Vector<double>>& values) const override final;
};

#endif

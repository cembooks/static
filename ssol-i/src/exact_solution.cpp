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

#include "exact_solution.hpp"
#include <cmath>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"

using namespace dealii;

ExactSolutionSSOLI_B::ExactSolutionSSOLI_B()
  : Function<3>(3)
  , B_0(2.0 * mu_0 * K_0 / 3.0)
{
}

void
ExactSolutionSSOLI_B::vector_value_list(
  const std::vector<Point<3>>& r,
  std::vector<Vector<double>>& values) const
{
  Assert(values.size() == r.size(),
         ExcDimensionMismatch(values.size(), r.size()));

  double cos_theta;
  double sin_theta;

  double cos_phi;
  double sin_phi;

  Tensor<1, 3> r_hat;
  Tensor<1, 3> theta_hat;

  Tensor<1, 3> B;

  auto v = values.begin();
  for (auto p : r) {
    cos_theta = p(2) / p.norm();
    sin_theta = sqrt(pow(p(0), 2) + pow(p(1), 2)) / p.norm();

    cos_phi = p(0) / sqrt(pow(p(0), 2) + pow(p(1), 2));
    sin_phi = p(1) / sqrt(pow(p(0), 2) + pow(p(1), 2));

    r_hat[0] = p(0) / p.norm();
    r_hat[1] = p(1) / p.norm();
    r_hat[2] = p(2) / p.norm();

    theta_hat[0] = cos_theta * cos_phi;
    theta_hat[1] = cos_theta * sin_phi;
    theta_hat[2] = -sin_theta;

    if (p.norm() < a) {
      B = B_0 * a * (cos_theta * r_hat - sin_theta * theta_hat);
    } else {
      B = B_0 * (pow(a, 4) / pow(p.norm(), 3)) *
          (cos_theta * r_hat + 0.5 * sin_theta * theta_hat);
    }

    (*v)[0] = B[0];
    (*v)[1] = B[1];
    (*v)[2] = B[2];

    v++;
  }
}

#pragma GCC diagnostic pop

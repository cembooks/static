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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

using namespace dealii;
using namespace std;

double
ExactSolutionSSOLIAXI_A::value(const Point<2>& p,
                               const unsigned int component) const
{
  double s = p.norm();
  double r = p[0];
  double M = mu_0 * K_0 * a * pow(r, 2) / 3.0;

  if (s < a) {
    return M;
  } else {
    return M * pow(a, 3) / pow(s, 3);
  }
}

Tensor<1, 2>
ExactSolutionSSOLIAXI_A::gradient(const Point<2>& p,
                                  const unsigned int component) const
{

  double s = p.norm();
  double r = p[0];
  double z = p[1];
  double M1 = mu_0 * K_0 * a / 3.0;
  double M2 = M1 * pow(a, 3);

  Point<2> grad_A;

  if (s < a) {
    grad_A(0) = 2.0 * M1 * r;
    grad_A(1) = 0.0;
  } else {
    grad_A(0) = M2 * (2.0 * r / pow(s, 3) - 3.0 * pow(r, 3) / pow(s, 5));
    grad_A(1) = -M2 * (3.0 * z * pow(r, 2) / pow(s, 5));
  }

  return grad_A;
}

void
ExactSolutionSSOLIAXI_B::vector_value_list(
  const std::vector<Point<2>>& r,
  std::vector<Vector<double>>& values) const
{
  double cos_theta;
  double sin_theta;

  double cos_phi;
  double sin_phi;

  Tensor<1, 3> r_hat;
  Tensor<1, 3> theta_hat;

  Tensor<1, 3> B;

  auto v = values.begin();
  for (auto pp : r) {
    Point<3> p(0.0, pp(0), pp(1));

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

    (*v)[0] = pp[0] * B[1];
    (*v)[1] = pp[0] * B[2];

    v++;
  }
}

#pragma GCC diagnostic pop

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
ExactSolutionSSOLIIAXI_A::value(const Point<2>& p,
                                const unsigned int component) const
{
  double r = p.norm();
  double sin_theta = p[0] / r;

  double G1 = mu_0 * K_0 * p[0] / 3.0;
  double G2 = mu_0 * K_0 * sin_theta / (3.0 * pow(r, 2));

  if (r < a)
    return p[0] * 0.5 * (pow(b, 2) - pow(a, 2)) * G1;

  if (r < b)
    return p[0] * (0.5 * (pow(b, 2) - pow(r, 2)) * G1 +
                   0.2 * (pow(r, 5) - pow(a, 5)) * G2);

  return p[0] * 0.2 * (pow(b, 5) - pow(a, 5)) * G2;
}

Tensor<1, 2>
ExactSolutionSSOLIIAXI_A::gradient(const Point<2>& r,
                                   const unsigned int component) const
{
  Point<2> grad_A;

  if (r.norm() <= a) {
    grad_A(0) = 0.0;
    grad_A(1) = 0.0;
  } else {
    grad_A(0) = 0.0;
    grad_A(1) = 0.0;
  }

  return grad_A;
}

void
ExactSolutionSSOLIIAXI_B::vector_value_list(
  const std::vector<Point<2>>& r,
  std::vector<Vector<double>>& values) const
{
  Assert(values.size() == r.size(),
         ExcDimensionMismatch(values.size(), r.size()));

  Tensor<1, 3> B;

  for (unsigned int i = 0; i < values.size(); i++) {
    B = magnetic_field(0.0, r.at(i)[0], r.at(i)[1], K_0, mu_0, a, b);

    values.at(i)[0] = r.at(i)[0] * B[1];
    values.at(i)[1] = r.at(i)[0] * B[2];
  }
}

void
ExactSolutionSSOLIIAXI_H::vector_value_list(
  const std::vector<Point<2>>& r,
  std::vector<Vector<double>>& values) const
{
  Assert(values.size() == r.size(),
         ExcDimensionMismatch(values.size(), r.size()));

  B.vector_value_list(r, values);

  for (unsigned int i = 0; i < r.size(); i++) {
    values.at(i)[0] = values.at(i)[0] / mu_0;
    values.at(i)[1] = values.at(i)[1] / mu_0;
  }
}

#pragma GCC diagnostic pop

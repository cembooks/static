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
  double s = p.norm();
  double r = p[0];

  double rG1 = mu_0 * K_0 * pow(r, 2) / 3.0;
  double rG2 = mu_0 * K_0 * pow(r, 2) / (3.0 * pow(s, 3));

  if (s < a)
    return 0.5 * (pow(b, 2) - pow(a, 2)) * rG1;

  if (s < b)
    return (0.5 * (pow(b, 2) - pow(s, 2)) * rG1 +
            0.2 * (pow(s, 5) - pow(a, 5)) * rG2);

  return 0.2 * (pow(b, 5) - pow(a, 5)) * rG2;
}

Tensor<1, 2>
ExactSolutionSSOLIIAXI_A::gradient(const Point<2>& p,
                                   const unsigned int component) const
{
  Point<2> grad_A, L1, L2;

  double s = p.norm();
  double r = p[0];
  double z = p[1];
  double M = (1.0 / 3.0) * mu_0 * K_0;

  L1(0) = 2.0 * M * r;
  L1(1) = 0.0;

  L2(0) = M * (2.0 * r / pow(s, 3) - 3.0 * pow(r, 3) / pow(s, 5));
  L2(1) = -M * 3.0 * pow(r, 2) * z / pow(s, 5);

  if (s < a) {
    grad_A = 0.5 * (pow(b, 2) - pow(a, 2)) * L1;
  } else if (s < b) {
    grad_A =
      0.5 * (pow(b, 2) - pow(s, 2)) * L1 + 0.2 * (pow(s, 5) - pow(a, 5)) * L2;
  } else {
    grad_A = 0.2 * (pow(b, 5) - pow(a, 5)) * L2;
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
    B = magnetic_field(0.0, r[i][0], r[i][1], K_0, mu_0, a, b);

    values[i][0] = r[i][0] * B[1];
    values[i][1] = r[i][0] * B[2];
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
    values[i][0] = values[i][0] / mu_0;
    values[i][1] = values[i][1] / mu_0;
  }
}

#pragma GCC diagnostic pop

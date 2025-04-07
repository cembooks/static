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

void
ExactSolutionSSOLIIIAXI_B::vector_value_list(
  const std::vector<Point<2>>& r,
  std::vector<Vector<double>>& values) const
{
  Assert(values.size() == r.size(),
         ExcDimensionMismatch(values.size(), r.size()));

  Tensor<1, 3> B;

  for (unsigned int i = 0; i < values.size(); i++) {
    B = magnetic_field_coil(0.0, r[i][0], r[i][1], K_0, mu_0, a2, b2) +
        magnetic_field_core(0.0, r[i][0], r[i][1], H_0, mu_r, mu_0, a1, b1);

    values[i][0] = r[i][0] * B[1];
    values[i][1] = r[i][0] * B[2];
  }
}

void
ExactSolutionSSOLIIIAXI_H::vector_value_list(
  const std::vector<Point<2>>& r,
  std::vector<Vector<double>>& values) const
{
  Assert(values.size() == r.size(),
         ExcDimensionMismatch(values.size(), r.size()));

  B.vector_value_list(r, values);

  double coef;

  for (unsigned int i = 0; i < r.size(); i++) {
    if ((r[i].norm() > a1) && (r[i].norm() < b1)) {
      coef = mu;
    } else {
      coef = mu_0;
    }

    values[i][0] = values[i][0] / coef;
    values[i][1] = values[i][1] / coef;
  }
}

#pragma GCC diagnostic pop

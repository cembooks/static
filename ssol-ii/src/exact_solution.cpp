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

ExactSolutionSSOLII_Jf::ExactSolutionSSOLII_Jf()
  : Function<3>(3)
{
}

void
ExactSolutionSSOLII_Jf::vector_value_list(
  const std::vector<Point<3>>& r,
  std::vector<Vector<double>>& values) const
{
  Assert(values.size() == r.size(),
         ExcDimensionMismatch(values.size(), r.size()));

  Tensor<1, 3> Jf;

  for (unsigned int i = 0; i < values.size(); i++) {
    Jf = volume_free_current_density(
      r.at(i)[0], r.at(i)[1], r.at(i)[2], K_0, a, b);

    values.at(i)[0] = Jf[0];
    values.at(i)[1] = Jf[1];
    values.at(i)[2] = Jf[2];
  }
}

ExactSolutionSSOLII_B::ExactSolutionSSOLII_B()
  : Function<3>(3)
{
}

void
ExactSolutionSSOLII_B::vector_value_list(
  const std::vector<Point<3>>& r,
  std::vector<Vector<double>>& values) const
{
  Assert(values.size() == r.size(),
         ExcDimensionMismatch(values.size(), r.size()));

  Tensor<1, 3> B;

  for (unsigned int i = 0; i < values.size(); i++) {
    B = magnetic_field(r.at(i)[0], r.at(i)[1], r.at(i)[2], K_0, mu_0, a, b);

    values.at(i)[0] = B[0];
    values.at(i)[1] = B[1];
    values.at(i)[2] = B[2];
  }
}

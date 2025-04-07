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

ExactSolutionMMSVTI_Jf::ExactSolutionMMSVTI_Jf()
  : Function<3>(3)
{
}

void
ExactSolutionMMSVTI_Jf::vector_value_list(
  const std::vector<Point<3>>& r,
  std::vector<Vector<double>>& values) const
{
  Assert(values.size() == r.size(),
         ExcDimensionMismatch(values.size(), r.size()));

  Tensor<1, 3> Jf;

  for (unsigned int i = 0; i < values.size(); i++) {
    Jf = volume_free_current_density(r[i][0], r[i][1], mu_0, k);

    values[i][0] = Jf[0];
    values[i][1] = Jf[1];
    values[i][2] = Jf[2];
  }
}

ExactSolutionMMSVTI_B::ExactSolutionMMSVTI_B()
  : Function<3>(3)
{
}

void
ExactSolutionMMSVTI_B::vector_value_list(
  const std::vector<Point<3>>& r,
  std::vector<Vector<double>>& values) const
{
  Assert(values.size() == r.size(),
         ExcDimensionMismatch(values.size(), r.size()));

  Tensor<1, 3> B;

  for (unsigned int i = 0; i < values.size(); i++) {
    B = magnetic_field(r[i][0], r[i][1], k);

    values[i][0] = B[0];
    values[i][1] = B[1];
    values[i][2] = B[2];
  }
}

DirichletBC_MMSVTI_T::DirichletBC_MMSVTI_T()
  : Function<3>(3)
{
}

void
DirichletBC_MMSVTI_T::vector_value_list(
  const std::vector<Point<3>>& r,
  std::vector<Vector<double>>& values) const
{
  Assert(values.size() == r.size(),
         ExcDimensionMismatch(values.size(), r.size()));

  Tensor<1, 3> T;

  for (unsigned int i = 0; i < values.size(); i++) {
    T = current_vector_potential(r[i][0], r[i][1], mu_0, k);

    values[i][0] = T[0];
    values[i][1] = T[1];
    values[i][2] = T[2];
  }
}

DirichletBC_MMSVTI_A::DirichletBC_MMSVTI_A()
  : Function<3>(3)
{
}

void
DirichletBC_MMSVTI_A::vector_value_list(
  const std::vector<Point<3>>& r,
  std::vector<Vector<double>>& values) const
{
  Assert(values.size() == r.size(),
         ExcDimensionMismatch(values.size(), r.size()));

  Tensor<1, 3> A;

  for (unsigned int i = 0; i < values.size(); i++) {
    A = magnetic_vector_potential(r[i][0], r[i][1], k);

    values[i][0] = A[0];
    values[i][1] = A[1];
    values[i][2] = A[2];
  }
}

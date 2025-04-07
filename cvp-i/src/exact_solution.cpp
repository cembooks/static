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

ExactSolutionCVPI_Jf::ExactSolutionCVPI_Jf()
  : Function<3>(3)
{
}

void
ExactSolutionCVPI_Jf::vector_value_list(
  const std::vector<Point<3>>& r,
  std::vector<Vector<double>>& values) const
{
  Assert(values.size() == r.size(),
         ExcDimensionMismatch(values.size(), r.size()));

  auto v = values.begin();
  for (auto p : r) {
    if ((p.norm() > SettingsCVPI::a) && (p.norm() < SettingsCVPI::b)) {
      (*v)[0] = -p[1];
      (*v)[1] = p[0];
      (*v)[2] = 0.0;
    } else {
      (*v)[0] = 0.0;
      (*v)[1] = 0.0;
      (*v)[2] = 0.0;
    }
    v++;
  }
}

DirichletBC_CVPI::DirichletBC_CVPI()
  : Function<3>(3)
{
}

void
DirichletBC_CVPI::vector_value_list(const std::vector<Point<3>>& r,
                                    std::vector<Vector<double>>& values) const
{
  Assert(values.size() == r.size(),
         ExcDimensionMismatch(values.size(), r.size()));

  for (unsigned int i = 0; i < values.size(); i++) {
    values[i][0] = 0.0;
    values[i][1] = 0.0;
    values[i][2] = 0.0;
  }
}

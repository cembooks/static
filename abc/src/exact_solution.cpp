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

template<>
ExactSolutionABC_PHI<2>::ExactSolutionABC_PHI()
  : Function<2>()
{
  Assert(R < x0, ExcInternalError());
  Assert(R > 0, ExcInternalError());
  Assert(x0 > 0, ExcInternalError());
  Assert(a > 0, ExcInternalError());
  Assert(a < 1.5 * (R + x0), ExcInternalError());

  d = sqrt(pow(x0, 2) - pow(R, 2));
  lambda = 1.0 / log(x0 / R + sqrt(pow(x0 / R, 2) - 1));
}

template<>
double
ExactSolutionABC_PHI<2>::value(const Point<2>& r, unsigned int component) const
{
  double xp = pow(r[0] + d, 2) + pow(r[1], 2);
  double xm = pow(r[0] - d, 2) + pow(r[1], 2);

  return (0.5 * lambda * log(xp / xm));
}

template<>
Tensor<1, 2>
ExactSolutionABC_PHI<2>::gradient(const Point<2>& r,
                                  unsigned int component) const
{
  Tensor<1, 2> p;

  double xp = pow(r[0] + d, 2) + pow(r[1], 2);
  double xm = pow(r[0] - d, 2) + pow(r[1], 2);

  p[0] = (r[0] + d) / xp - (r[0] - d) / xm;

  p[1] = r[1] / xp - r[1] / xm;

  return lambda * p;
}

template<>
ExactSolutionABC_PHI<3>::ExactSolutionABC_PHI()
  : Function<3>()
{
  Assert(R < x0, ExcInternalError());
  Assert(R > 0, ExcInternalError());
  Assert(x0 > 0, ExcInternalError());
  Assert(a > 0, ExcInternalError());
  Assert(a < 1.5 * (R + x0), ExcInternalError());
}

template<>
double
ExactSolutionABC_PHI<3>::value(const Point<3>& r, unsigned int componet) const
{
  return a / r.norm();
}

template<>
Tensor<1, 3>
ExactSolutionABC_PHI<3>::gradient(const Point<3>& r,
                                  unsigned int component) const
{
  return -r * a / pow(r.norm(), 3);
}

#pragma GCC diagnostic pop

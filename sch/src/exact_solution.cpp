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
double
ExactSolutionSCH_PHI<2>::value(const Point<2>& r,
                               const unsigned int component) const
{
  if (r.norm() > a) {
    return (log(b) - log(sqrt(r.square()))) / (log(b) - log(a));
  }

  return 1.0;
}

template<>
Tensor<1, 2>
ExactSolutionSCH_PHI<2>::gradient(const Point<2>& r,
                                  const unsigned int component) const
{
  if (r.norm() > a) {
    return -1.0 / (log(b) - log(a)) * r / r.square();
  }

  return Point<2>();
}

template<>
double
ExactSolutionSCH_PHI<3>::value(const Point<3>& r,
                               const unsigned int component) const
{
  if (r.norm() > a) {
    return (a * b / (b - a)) * (1 / sqrt(r.square()) - 1 / b);
  }

  return 1.0;
}

template<>
Tensor<1, 3>
ExactSolutionSCH_PHI<3>::gradient(const Point<3>& r,
                                  const unsigned int component) const
{
  if (r.norm() > a) {
    return -a * b / (b - a) * r / pow(sqrt(r.square()), 3);
  }

  return Point<3>();
}

#pragma GCC diagnostic pop

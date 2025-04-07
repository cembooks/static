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

using namespace dealii;
using namespace std;

template<>
ExactSolutionINTAXI_PHI<true>::ExactSolutionINTAXI_PHI()
{
  alpha = 1.0 / (ep_1 * log(b / d) + ep_2 * log(d / a));
  beta = (ep_1 / ep_2) * log(b / d);
}

template<>
double
ExactSolutionINTAXI_PHI<true>::value(const Point<2>& p,
                                     const unsigned int component) const
{
  Point<2> r;
  r[0] = p[0];
  r[1] = 0.0;

  if (r.norm() < d)
    return (-alpha * ep_2 * (log(r.norm() / d) - beta));
  else
    return (-alpha * ep_1 * log(r.norm() / b));
}

template<>
Tensor<1, 2>
ExactSolutionINTAXI_PHI<true>::gradient(const Point<2>& p,
                                        const unsigned int component) const
{
  Point<2> r;
  r[0] = p[0];
  r[1] = 0.0;

  if (r.norm() < d)
    return -alpha * ep_2 * (1 / r.square()) * r;
  else
    return -alpha * ep_1 * (1 / r.square()) * r;
}

template<>
ExactSolutionINTAXI_PHI<false>::ExactSolutionINTAXI_PHI()
{
  alpha = 1.0 / (ep_1 * (1.0 / b - 1.0 / d) + ep_2 * (1.0 / d - 1.0 / a));
  beta = 1.0 / d + (ep_1 / ep_2) * (1.0 / b - 1.0 / d);
}

template<>
double
ExactSolutionINTAXI_PHI<false>::value(const Point<2>& p,
                                      const unsigned int component) const
{
  Point<2> r;
  r[0] = p[0];
  r[1] = p[1];

  if (r.norm() < d)
    return (-alpha * ep_2 * (1.0 / r.norm() - beta));
  else
    return (-alpha * ep_1 * (1.0 / r.norm() - 1.0 / b));
}

template<>
Tensor<1, 2>
ExactSolutionINTAXI_PHI<false>::gradient(const Point<2>& p,
                                         const unsigned int component) const
{
  Point<2> r;
  r[0] = p[0];
  r[1] = p[1];

  if (r.norm() < d)
    return alpha * ep_2 * r / pow(r.norm(), 3);
  else
    return alpha * ep_1 * r / pow(r.norm(), 3);
}

template<>
void
ExactSolutionINTAXI_E<true>::vector_value_list(
  const std::vector<Point<2>>& r,
  std::vector<Vector<double>>& values) const
{
  Assert(values.size() == r.size(),
         ExcDimensionMismatch(values.size(), r.size()));

  Tensor<1, 2> phi;

  for (unsigned int i = 0; i < r.size(); i++) {
    phi = PHI.gradient(r[i]);

    values[i][0] = -phi[0];
    values[i][1] = -phi[1];
  }
}

template<>
void
ExactSolutionINTAXI_E<false>::vector_value_list(
  const std::vector<Point<2>>& r,
  std::vector<Vector<double>>& values) const
{
  Assert(values.size() == r.size(),
         ExcDimensionMismatch(values.size(), r.size()));

  Tensor<1, 2> phi;

  for (unsigned int i = 0; i < r.size(); i++) {
    phi = PHI.gradient(r[i]);

    values[i][0] = -phi[0];
    values[i][1] = -phi[1];
  }
}

template<>
void
ExactSolutionINTAXI_D<true>::vector_value_list(
  const std::vector<Point<2>>& r,
  std::vector<Vector<double>>& values) const
{
  Assert(values.size() == r.size(),
         ExcDimensionMismatch(values.size(), r.size()));

  double ep;
  Tensor<1, 2> phi;

  for (unsigned int i = 0; i < r.size(); i++) {
    ep = ep_2;

    if (r[i][0] < d)
      ep = ep_1;

    phi = PHI.gradient(r[i]);

    values[i][0] = -ep * phi[0];
    values[i][1] = -ep * phi[1];
  }
}

template<>
void
ExactSolutionINTAXI_D<false>::vector_value_list(
  const std::vector<Point<2>>& r,
  std::vector<Vector<double>>& values) const
{
  Assert(values.size() == r.size(),
         ExcDimensionMismatch(values.size(), r.size()));

  double ep;
  Tensor<1, 2> phi;

  for (unsigned int i = 0; i < r.size(); i++) {
    ep = ep_2;

    if (r[i].norm() < d)
      ep = ep_1;

    phi = PHI.gradient(r[i]);

    values[i][0] = -ep * phi[0];
    values[i][1] = -ep * phi[1];
  }
}

#pragma GCC diagnostic pop

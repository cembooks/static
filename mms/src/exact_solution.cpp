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
ExactSolutionMMS_PHI<2>::value(const Point<2>& r, unsigned int component) const
{
  return (cos(k * r[0]) + cos(k * r[1]));
}

template<>
Tensor<1, 2>
ExactSolutionMMS_PHI<2>::gradient(const Point<2>& r,
                                  unsigned int component) const
{
  Tensor<1, 2> p;

  p[0] = -k * sin(k * r[0]);
  p[1] = -k * sin(k * r[1]);

  return p;
}

template<>
double
ExactSolutionMMS_PHI<3>::value(const Point<3>& r, unsigned int componet) const
{
  return (cos(k * r[0]) + cos(k * r[1]) + cos(k * r[2]));
}

template<>
Tensor<1, 3>
ExactSolutionMMS_PHI<3>::gradient(const Point<3>& r,
                                  unsigned int component) const
{
  Tensor<1, 3> p;

  p[0] = -k * sin(k * r[0]);
  p[1] = -k * sin(k * r[1]);
  p[2] = -k * sin(k * r[2]);

  return p;
}

template<>
void
ExactSolutionMMS_E<2>::vector_value_list(
  const std::vector<Point<2>>& r,
  std::vector<Vector<double>>& values) const
{
  Assert(values.size() == r.size(),
         ExcDimensionMismatch(values.size(), r.size()));

  for (unsigned int i = 0; i < r.size(); i++) {
    values[i][0] = k * sin(k * r[i][0]);
    values[i][1] = k * sin(k * r[i][1]);
  }
}

template<>
void
ExactSolutionMMS_E<3>::vector_value_list(
  const std::vector<Point<3>>& r,
  std::vector<Vector<double>>& values) const
{
  Assert(values.size() == r.size(),
         ExcDimensionMismatch(values.size(), r.size()));

  for (unsigned int i = 0; i < r.size(); i++) {
    values[i][0] = k * sin(k * r[i][0]);
    values[i][1] = k * sin(k * r[i][1]);
    values[i][2] = k * sin(k * r[i][2]);
  }
}

template<>
void
ExactSolutionMMS_D<2>::vector_value_list(
  const std::vector<Point<2>>& r,
  std::vector<Vector<double>>& values) const
{
  Assert(values.size() == r.size(),
         ExcDimensionMismatch(values.size(), r.size()));

  double alpha;

  for (unsigned int i = 0; i < r.size(); i++) {
    alpha = ep_0 * k * (pow(r[i][0], 2.0) * pow(r[i][1], 2.0) + 1.0);

    values[i][0] = alpha * sin(k * r[i][0]);
    values[i][1] = alpha * sin(k * r[i][1]);
  }
}

template<>
void
ExactSolutionMMS_D<3>::vector_value_list(
  const std::vector<Point<3>>& r,
  std::vector<Vector<double>>& values) const
{
  Assert(values.size() == r.size(),
         ExcDimensionMismatch(values.size(), r.size()));

  double alpha;

  for (unsigned int i = 0; i < r.size(); i++) {
    alpha = ep_0 * k *
            (pow(r[i][0], 2.0) * pow(r[i][1], 2.0) * pow(r[i][2], 2.0) + 1.0);

    values[i][0] = alpha * sin(k * r[i][0]);
    values[i][1] = alpha * sin(k * r[i][1]);
    values[i][2] = alpha * sin(k * r[i][2]);
  }
}

#pragma GCC diagnostic pop

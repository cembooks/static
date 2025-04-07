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

template<>
ExactSolutionMMSV_A<3>::ExactSolutionMMSV_A()
  : Function<3>(3)
{
}

template<>
void
ExactSolutionMMSV_A<3>::vector_value_list(
  const std::vector<Point<3>>& r,
  std::vector<Vector<double>>& values) const
{
  Assert(values.size() == r.size(),
         ExcDimensionMismatch(values.size(), r.size()));

  auto v = values.begin();
  for (auto p : r) {
    (*v)(0) = -sin(k * p[1]) / k;
    (*v)(1) = sin(k * p[0]) / k;
    (*v)(2) = 0.0;
    v++;
  }
}

template<>
ExactSolutionMMSV_A<2>::ExactSolutionMMSV_A()
  : Function<2>(2)
{
}

template<>
void
ExactSolutionMMSV_A<2>::vector_value_list(
  const std::vector<Point<2>>& r,
  std::vector<Vector<double>>& values) const
{
  Assert(values.size() == r.size(),
         ExcDimensionMismatch(values.size(), r.size()));

  auto v = values.begin();
  for (auto p : r) {
    (*v)(0) = -sin(k * p[1]) / k;
    (*v)(1) = sin(k * p[0]) / k;

    v++;
  }
}

template<>
ExactSolutionMMSV_B<3>::ExactSolutionMMSV_B()
  : Function<3>(3)
{
}

template<>
void
ExactSolutionMMSV_B<3>::vector_value_list(
  const std::vector<Point<3>>& r,
  std::vector<Vector<double>>& values) const
{
  Assert(values.size() == r.size(),
         ExcDimensionMismatch(values.size(), r.size()));

  auto v = values.begin();
  for (auto p : r) {
    (*v)(0) = 0.0;
    (*v)(1) = 0.0;
    (*v)(2) = cos(k * p[0]) + cos(k * p[1]);
    v++;
  }
}

// Is not used. This prevents messages from linker.
template<>
double
ExactSolutionMMSV_B<3>::value(const Point<3>& r,
                              const unsigned int component) const
{
  return 0.0;
}

template<>
ExactSolutionMMSV_B<2>::ExactSolutionMMSV_B()
{
}

// Is not used. This prevents messages from linker.
template<>
void
ExactSolutionMMSV_B<2>::vector_value_list(
  const std::vector<Point<2>>& r,
  std::vector<Vector<double>>& values) const
{
}

template<>
double
ExactSolutionMMSV_B<2>::value(const Point<2>& r,
                              const unsigned int component) const
{
  return cos(k * r[0]) + cos(k * r[1]);
}

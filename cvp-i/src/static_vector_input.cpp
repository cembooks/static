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

#include <deal.II/base/types.h>
#define BOOST_ALLOW_DEPRECATED_HEADERS

#include "static_vector_input.hpp"
#include <math.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"

template<>
void
StaticVectorSolver::TheCoefficient<3>::value_list(
  const std::vector<Point<3>>& r,
  types::material_id mid,
  unsigned int cuid,
  std::vector<double>& values) const
{
  Assert(r.size() == values.size(),
         ExcDimensionMismatch(r.size(), values.size()));

  for (unsigned int i = 0; i < values.size(); i++)
    values[i] = 1.0;
}

template<>
void
StaticVectorSolver::PdeRhs<3>::value_list(
  const std::vector<Point<3>>& r,
  types::material_id mid,
  unsigned int cuid,
  std::vector<Tensor<1, 3>>& values) const
{
  Assert(r.size() == values.size(),
         ExcDimensionMismatch(r.size(), values.size()));

  auto v = values.begin();
  for (auto p : r) {
    if (mid == mid_2) {
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

template<>
void
StaticVectorSolver::Gamma<3>::value_list(const std::vector<Point<3>>& r,
                                         const std::vector<Tensor<1, 3>>& n,
                                         types::boundary_id bid,
                                         types::material_id mid,
                                         unsigned int cuid,
                                         unsigned int fuid,
                                         std::vector<double>& values) const
{
  Assert(r.size() == values.size(),
         ExcDimensionMismatch(r.size(), values.size()));

  for (unsigned int i = 0; i < values.size(); i++)
    values[i] = 0.0;
}

template<>
void
StaticVectorSolver::RobinRhs<3>::value_list(
  const std::vector<Point<3>>& r,
  const std::vector<Tensor<1, 3>>& n,
  types::boundary_id bid,
  types::material_id mid,
  unsigned int cuid,
  unsigned int fuid,
  std::vector<Tensor<1, 3>>& values) const
{
  Assert(r.size() == values.size(),
         ExcDimensionMismatch(r.size(), values.size()));

  for (unsigned int i = 0; i < values.size(); i++) {
    values[i][0] = 0.0;
    values[i][1] = 0.0;
    values[i][2] = 0.0;
  }
}

template<>
void
StaticVectorSolver::FreeSurfaceCurrent<3>::value_list(
  const std::vector<Point<3>>& r,
  const std::vector<Tensor<1, 3>>& n,
  types::material_id mid,
  unsigned int cuid,
  unsigned int fuid,
  std::vector<Tensor<1, 3>>& values) const
{
  Assert(r.size() == values.size(),
         ExcDimensionMismatch(r.size(), values.size()));

  for (unsigned int i = 0; i < values.size(); i++) {
    values[i][0] = 0.0;
    values[i][1] = 0.0;
    values[i][2] = 0.0;
  }
}

template<>
double
StaticVectorSolver::Weight<3>::value(const Point<3>& r,
                                     const unsigned int component) const
{
  return 1.0;
}

#pragma GCC diagnostic pop

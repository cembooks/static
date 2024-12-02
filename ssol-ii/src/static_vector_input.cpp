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

#include "exact_solution.hpp"
#include "static_vector_input.hpp"
#include <math.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"

//-----------------------------------------------------------------------------
//-------- Stage 0. Calculating current vector potential, T, given Jf ---------
//-----------------------------------------------------------------------------

template<>
void
StaticVectorSolver::TheCoefficient<3, 0>::value_list(
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
StaticVectorSolver::PdeRhs<3, 0>::value_list(
  const std::vector<Point<3>>& r,
  types::material_id mid,
  unsigned int cuid,
  std::vector<Tensor<1, 3>>& values) const
{
  Assert(r.size() == values.size(),
         ExcDimensionMismatch(r.size(), values.size()));

  if (mid == mid_1)
    for (unsigned int i = 0; i < values.size(); i++) {
      values.at(i)[0] = 0.0;
      values.at(i)[1] = 0.0;
      values.at(i)[2] = 0.0;
    }

  if (mid == mid_2)
    for (unsigned int i = 0; i < values.size(); i++) {
      values.at(i)[0] = -K_0 * r.at(i)[1];
      values.at(i)[1] = K_0 * r.at(i)[0];
      values.at(i)[2] = 0.0;
    }
}

template<>
void
StaticVectorSolver::Gamma<3, 0>::value_list(const std::vector<Point<3>>& r,
                                            const std::vector<Tensor<1, 3>>& n,
                                            types::boundary_id bid,
                                            types::material_id mid,
                                            unsigned int cuid,
                                            unsigned int fuid,
                                            std::vector<double>& values) const
{
  Assert(r.size() == values.size(),
         ExcDimensionMismatch(r.size(), values.size()));

  Assert(r.size() == n.size(), ExcDimensionMismatch(r.size(), n.size()));

  for (unsigned int i = 0; i < values.size(); i++)
    values[i] = 1.0 / r.at(i).norm();
}

template<>
void
StaticVectorSolver::RobinRhs<3, 0>::value_list(
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

  Assert(r.size() == n.size(), ExcDimensionMismatch(r.size(), n.size()));

  Tensor<1, 3> Jf;

  for (unsigned int i = 0; i < r.size(); i++) {
    values.at(i)[0] = 0.0;
    values.at(i)[1] = 0.0;
    values.at(i)[2] = 0.0;
  }
}

template<>
void
StaticVectorSolver::FreeSurfaceCurrent<3, 0>::value_list(
  const std::vector<Point<3>>& r,
  const std::vector<Tensor<1, 3>>& n,
  types::material_id mid,
  unsigned int cuid,
  unsigned int fuid,
  std::vector<Tensor<1, 3>>& values) const
{
  Assert(r.size() == values.size(),
         ExcDimensionMismatch(r.size(), values.size()));

  Assert(r.size() == n.size(), ExcDimensionMismatch(r.size(), n.size()));

  for (unsigned int i = 0; i < values.size(); i++) {
    values[i][0] = 0.0;
    values[i][1] = 0.0;
    values[i][2] = 0.0;
  }
}

template<>
double
StaticVectorSolver::Weight<3, 0>::value(const Point<3>& r,
                                        const unsigned int component) const
{
  return 1.0;
}

//-----------------------------------------------------------------------------
// Stage 1. Calculating Jf given T (projection) to check if the Stage 0 was ok.
//-----------------------------------------------------------------------------

template<>
double
StaticVectorSolver::Weight<3, 1>::value(const Point<3>& r,
                                        const unsigned int component) const
{
  return 1.0;
}

//-----------------------------------------------------------------------------
// Stage 2. Calculating vector potential, A, given T --------------------------
//-----------------------------------------------------------------------------

template<>
void
StaticVectorSolver::TheCoefficient<3, 2>::value_list(
  const std::vector<Point<3>>& r,
  types::material_id mid,
  unsigned int cuid,
  std::vector<double>& values) const
{
  Assert(r.size() == values.size(),
         ExcDimensionMismatch(r.size(), values.size()));

  for (unsigned int i = 0; i < values.size(); i++)
    values[i] = mu_0;
}

template<>
void
StaticVectorSolver::PdeRhs<3, 2>::value_list(
  const std::vector<Point<3>>& r,
  types::material_id mid,
  unsigned int cuid,
  std::vector<Tensor<1, 3>>& values) const
{
  Assert(r.size() == values.size(),
         ExcDimensionMismatch(r.size(), values.size()));

  for (unsigned int i = 0; i < values.size(); i++) {
    values.at(i)[0] = 0.0;
    values.at(i)[1] = 0.0;
    values.at(i)[2] = 0.0;
  }
}

template<>
void
StaticVectorSolver::Gamma<3, 2>::value_list(const std::vector<Point<3>>& r,
                                            const std::vector<Tensor<1, 3>>& n,
                                            types::boundary_id bid,
                                            types::material_id mid,
                                            unsigned int cuid,
                                            unsigned int fuid,
                                            std::vector<double>& values) const
{
  Assert(r.size() == values.size(),
         ExcDimensionMismatch(r.size(), values.size()));

  Assert(r.size() == n.size(), ExcDimensionMismatch(r.size(), n.size()));

  for (unsigned int i = 0; i < values.size(); i++)
    values[i] = 1.0 / (mu_0 * r.at(i).norm());
}

template<>
void
StaticVectorSolver::RobinRhs<3, 2>::value_list(
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

  Assert(r.size() == n.size(), ExcDimensionMismatch(r.size(), n.size()));

  for (unsigned int i = 0; i < r.size(); i++) {
    values.at(i) = 0.0;
  }
}

template<>
void
StaticVectorSolver::FreeSurfaceCurrent<3, 2>::value_list(
  const std::vector<Point<3>>& r,
  const std::vector<Tensor<1, 3>>& n,
  types::material_id mid,
  unsigned int cuid,
  unsigned int fuid,
  std::vector<Tensor<1, 3>>& values) const
{
  Assert(r.size() == values.size(),
         ExcDimensionMismatch(r.size(), values.size()));

  Assert(r.size() == n.size(), ExcDimensionMismatch(r.size(), n.size()));

  for (unsigned int i = 0; i < values.size(); i++) {
    values[i][0] = 0.0;
    values[i][1] = 0.0;
    values[i][2] = 0.0;
  }
}

template<>
double
StaticVectorSolver::Weight<3, 2>::value(const Point<3>& r,
                                        const unsigned int component) const
{
  return 1.0;
}
//-----------------------------------------------------------------------------
//-------- Stage 3. Calculating B given A -------------------------------------
//-----------------------------------------------------------------------------

template<>
double
StaticVectorSolver::Weight<3, 3>::value(const Point<3>& r,
                                        const unsigned int component) const
{
  return 1.0;
}

#pragma GCC diagnostic pop

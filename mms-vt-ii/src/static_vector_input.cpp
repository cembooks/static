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
// Stage 2. Calculating vector potential, A, given T --------------------------
//-----------------------------------------------------------------------------

template<>
void
StaticVectorSolver::TheCoefficient<2, 2>::value_list(
  const std::vector<Point<2>>& r,
  types::material_id mid,
  unsigned int cuid,
  std::vector<double>& values) const
{
  Assert(r.size() == values.size(),
         ExcDimensionMismatch(r.size(), values.size()));

  for (unsigned int i = 0; i < values.size(); i++)
    values[i] = permeability(r.at(i)[0], r.at(i)[1], mu_0);
}

template<>
void
StaticVectorSolver::PdeRhs<2, 2>::value_list(
  const std::vector<Point<2>>& r,
  types::material_id mid,
  unsigned int cuid,
  std::vector<Tensor<1, 2>>& values) const
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
StaticVectorSolver::Gamma<2, 2>::value_list(const std::vector<Point<2>>& r,
                                            const std::vector<Tensor<1, 2>>& n,
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
    values[i] = robin_gamma(r.at(i)[0], r.at(i)[1], mu_0);
}

template<>
void
StaticVectorSolver::RobinRhs<2, 2>::value_list(
  const std::vector<Point<2>>& r,
  const std::vector<Tensor<1, 2>>& n,
  types::boundary_id bid,
  types::material_id mid,
  unsigned int cuid,
  unsigned int fuid,
  std::vector<Tensor<1, 2>>& values) const
{
  Assert(r.size() == values.size(),
         ExcDimensionMismatch(r.size(), values.size()));

  Assert(r.size() == n.size(), ExcDimensionMismatch(r.size(), n.size()));

  Tensor<1, 2> A;
  double mu, gamma, T;

  for (unsigned int i = 0; i < r.size(); i++) {
    mu = permeability(r.at(i)[0], r.at(i)[1], mu_0);
    gamma = robin_gamma(r.at(i)[0], r.at(i)[1], mu_0);
    T = current_vector_potential(r.at(i)[0], r.at(i)[1], mu_0, k);
    A = vector_potential(r.at(i)[0], r.at(i)[1], k);

    values.at(i)[0] =
      n.at(i)[1] * T +
      gamma * n.at(i)[1] * (n.at(i)[0] * A[1] - n.at(i)[1] * A[0]);
    values.at(i)[1] =
      -n.at(i)[0] * T -
      gamma * n.at(i)[0] * (n.at(i)[0] * A[1] - n.at(i)[1] * A[0]);
  }
}

template<>
void
StaticVectorSolver::FreeSurfaceCurrent<2, 2>::value_list(
  const std::vector<Point<2>>& r,
  const std::vector<Tensor<1, 2>>& n,
  types::material_id mid,
  unsigned int cuid,
  unsigned int fuid,
  std::vector<Tensor<1, 2>>& values) const
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
StaticVectorSolver::Weight<2, 2>::value(const Point<2>& r,
                                        const unsigned int component) const
{
  return 1.0;
}

//-----------------------------------------------------------------------------
//-------- Stage 3. Calculating B given A -------------------------------------
//-----------------------------------------------------------------------------

template<>
double
StaticVectorSolver::Weight<2, 3>::value(const Point<2>& r,
                                        const unsigned int component) const
{
  return 1.0;
}

#pragma GCC diagnostic pop

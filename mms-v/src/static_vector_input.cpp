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
#include <iostream>
#include <math.h>

using namespace StaticVectorSolver;
using namespace std;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"

template<>
void
TheCoefficient<3>::value_list(const std::vector<Point<3>>& r,
                              types::material_id mid,
                              unsigned int cuid,
                              std::vector<double>& values) const
{
  Assert(r.size() == values.size(),
         ExcDimensionMismatch(r.size(), values.size()));

  auto v = values.begin();
  for (auto p : r) {
    *v = mu_0 * (pow(p[0], 2) + pow(p[1], 2) + 1);
    v++;
  }
}

template<>
void
TheCoefficient<2>::value_list(const std::vector<Point<2>>& r,
                              types::material_id mid,
                              unsigned int cuid,
                              std::vector<double>& values) const
{
  Assert(r.size() == values.size(),
         ExcDimensionMismatch(r.size(), values.size()));

  auto v = values.begin();
  for (auto p : r) {
    *v = mu_0 * (pow(p[0], 2) + pow(p[1], 2) + 1);
    v++;
  }
}

template<>
void
PdeRhs<3>::value_list(const std::vector<Point<3>>& r,
                      types::material_id mid,
                      unsigned int cuid,
                      std::vector<Tensor<1, 3>>& values) const
{
  Assert(r.size() == values.size(),
         ExcDimensionMismatch(r.size(), values.size()));

  auto v = values.begin();
  for (auto p : r) {
    (*v)[0] = 0.0;
    (*v)[1] = 0.0;
    (*v)[2] = (cos(k * p[0]) + cos(k * p[1])) /
              ((mu_0 * (pow(p[0], 2) + pow(p[1], 2) + 1.0)));
    v++;
  }
}

template<>
void
PdeRhs<2>::value_list(const std::vector<Point<2>>& r,
                      types::material_id mid,
                      unsigned int cuid,
                      std::vector<Tensor<1, 2>>& values) const
{
  Assert(r.size() == values.size(),
         ExcDimensionMismatch(r.size(), values.size()));

  auto v = values.begin();
  for (auto p : r) {
    (*v)[0] = (cos(k * p[0]) + cos(k * p[1])) /
              ((mu_0 * (pow(p[0], 2) + pow(p[1], 2) + 1.0)));
    (*v)[1] = 0.0;
    v++;
  }
}

template<>
void
Gamma<3>::value_list(const std::vector<Point<3>>& r,
                     const std::vector<Tensor<1, 3>>& n,
                     types::boundary_id bid,
                     types::material_id mid,
                     unsigned int cuid,
                     unsigned int fuid,
                     std::vector<double>& values) const
{
  Assert(r.size() == values.size(),
         ExcDimensionMismatch(r.size(), values.size()));

  auto v = values.begin();
  for (auto p : r) {
    double mu = mu_0 * (pow(p[0], 2) + pow(p[1], 2) + 1.0);
    *v = (1.0 / mu) * (sqrt(pow(p[0], 2) + pow(p[1], 2)) + 2.0);
    v++;
  }
}

template<>
void
Gamma<2>::value_list(const std::vector<Point<2>>& r,
                     const std::vector<Tensor<1, 2>>& n,
                     types::boundary_id bid,
                     types::material_id mid,
                     unsigned int cuid,
                     unsigned int fuid,
                     std::vector<double>& values) const
{
  Assert(r.size() == values.size(),
         ExcDimensionMismatch(r.size(), values.size()));

  auto v = values.begin();
  for (auto p : r) {
    double mu = mu_0 * (pow(p[0], 2) + pow(p[1], 2) + 1.0);
    *v = (1.0 / mu) * (sqrt(pow(p[0], 2) + pow(p[1], 2)) + 2.0);
    v++;
  }
}

template<>
void
RobinRhs<3>::value_list(const std::vector<Point<3>>& r,
                        const std::vector<Tensor<1, 3>>& n,
                        types::boundary_id bid,
                        types::material_id mid,
                        unsigned int cuid,
                        unsigned int fuid,
                        std::vector<Tensor<1, 3>>& values) const
{
  Assert(r.size() == values.size(),
         ExcDimensionMismatch(r.size(), values.size()));

  double mu;
  double gamma;
  Tensor<1, 3> T;
  Tensor<1, 3> A;
  Tensor<1, 3> Q;

  auto v = values.begin();
  auto nn = n.begin();
  for (auto p : r) {
    mu = mu_0 * (pow(p[0], 2) + pow(p[1], 2) + 1.0);
    gamma = (1 / mu) * (sqrt(pow(p[0], 2) + pow(p[1], 2)) + 2.0);

    T[0] = 0.0;
    T[1] = 0.0;
    T[2] = (cos(k * p[0]) + cos(k * p[1])) / mu;

    A[0] = -sin(k * p[1]) / k;
    A[1] = sin(k * p[0]) / k;
    A[2] = 0.0;

    Q = cross_product_3d(*nn, T) +
        gamma * cross_product_3d(*nn, cross_product_3d(*nn, A));

    (*v)[0] = Q[0];
    (*v)[1] = Q[1];
    (*v)[2] = Q[2];
    v++;
    nn++;
  }
}

template<>
void
RobinRhs<2>::value_list(const std::vector<Point<2>>& r,
                        const std::vector<Tensor<1, 2>>& n,
                        types::boundary_id bid,
                        types::material_id mid,
                        unsigned int cuid,
                        unsigned int fuid,
                        std::vector<Tensor<1, 2>>& values) const
{
  Assert(r.size() == values.size(),
         ExcDimensionMismatch(r.size(), values.size()));

  double mu;
  double gamma;
  double T;
  Tensor<1, 2> A;
  Tensor<1, 2> Q;

  auto v = values.begin();
  auto nn = n.begin();
  for (auto p : r) {
    mu = mu_0 * (pow(p[0], 2) + pow(p[1], 2) + 1.0);
    gamma = (1 / mu) * (sqrt(pow(p[0], 2) + pow(p[1], 2)) + 2.0);
    T = (cos(k * p[0]) + cos(k * p[1])) / mu;

    A[0] = -sin(k * p[1]) / k;
    A[1] = sin(k * p[0]) / k;

    Q[0] =
      (*nn)[1] * T + gamma * (*nn)[1] * ((*nn)[0] * A[1] - (*nn)[1] * A[0]);
    Q[1] =
      -(*nn)[0] * T - gamma * (*nn)[0] * ((*nn)[0] * A[1] - (*nn)[1] * A[0]);

    (*v)[0] = Q[0];
    (*v)[1] = Q[1];

    v++;
    nn++;
  }
}

template<>
void
FreeSurfaceCurrent<3>::value_list(const std::vector<Point<3>>& r,
                                  const std::vector<Tensor<1, 3>>& n,
                                  types::material_id mid,
                                  unsigned int cuid,
                                  unsigned int fuid,
                                  std::vector<Tensor<1, 3>>& values) const
{
  Assert(r.size() == values.size(),
         ExcDimensionMismatch(r.size(), values.size()));

  auto v = values.begin();
  for (auto p : r) {
    (*v)[0] = 0.0;
    (*v)[1] = 0.0;
    (*v)[2] = 0.0;
    v++;
  }
}

template<>
void
FreeSurfaceCurrent<2>::value_list(const std::vector<Point<2>>& r,
                                  const std::vector<Tensor<1, 2>>& n,
                                  types::material_id mid,
                                  unsigned int cuid,
                                  unsigned int fuid,
                                  std::vector<Tensor<1, 2>>& values) const
{
  Assert(r.size() == values.size(),
         ExcDimensionMismatch(r.size(), values.size()));

  auto v = values.begin();
  for (auto p : r) {
    (*v)[0] = 0.0;
    (*v)[1] = 0.0;
    v++;
  }
}

template<>
double
Weight<3>::value(const Point<3>& r, const unsigned int component) const
{
  return 1.0;
}

template<>
double
Weight<2>::value(const Point<2>& r, const unsigned int component) const
{
  return 1.0;
}

#pragma GCC diagnostic pop

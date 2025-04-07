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

#ifndef ExactSolutionsSLDI_H__
#define ExactSolutionsSLDI_H__

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/vector.h>

#include <cmath>

#include "constants.hpp"
#include "settings.hpp"

using namespace dealii;

/**
 * \brief Describes exact solution, \f$\Psi\f$, of the
 * *Magnetostatic shield - 1* [(sld-i/)](@ref page_sld_i)
 * numerical experiment in two and three dimensions.
 *****************************************************************************/
template<int dim>
class ExactSolutionSLDI_PSI
  : public Function<dim>
  , public SettingsSLDI
{
public:
  ExactSolutionSLDI_PSI();

  virtual double value(const Point<dim>& r,
                       const unsigned int component = 0) const override final;

  virtual Tensor<1, dim> gradient(
    const Point<dim>& r,
    const unsigned int component = 0) const override final;

private:
  double OMEGA, alpha_1, beta_1, gamma_1, delta_1;
};

/**
 * \brief Describes exact solution, \f$\vec{H}\f$, of the
 * *Magnetostatic shield - 1* [(sld-i/)](@ref page_sld_i)
 * numerical experiment in two and three dimensions.
 *****************************************************************************/
template<int dim>
class ExactSolutionSLDI_H
  : public Function<dim>
  , public SettingsSLDI
{
public:
  ExactSolutionSLDI_H()
    : Function<dim>(dim)
  {
  }

  virtual void vector_value_list(
    const std::vector<Point<dim>>& r,
    std::vector<Vector<double>>& values) const final;

private:
  ExactSolutionSLDI_PSI<dim> psi;
};

/**
 * \brief Describes exact solution, \f$\vec{B}\f$, of the
 * *Magnetostatic shield - 1* [(sld-i/)](@ref page_sld_i)
 * numerical experiment in two and three dimensions.
 *****************************************************************************/
template<int dim>
class ExactSolutionSLDI_B
  : public Function<dim>
  , public SettingsSLDI
{
public:
  ExactSolutionSLDI_B()
    : Function<dim>(dim)
  {
  }

  virtual void vector_value_list(
    const std::vector<Point<dim>>& r,
    std::vector<Vector<double>>& values) const final;

private:
  ExactSolutionSLDI_PSI<dim> psi;
};

template<int dim>
void
ExactSolutionSLDI_H<dim>::vector_value_list(
  const std::vector<Point<dim>>& r,
  std::vector<Vector<double>>& values) const
{
  Assert(values.size() == r.size(),
         ExcDimensionMismatch(values.size(), r.size()));

  Tensor<1, dim> grad;

  for (unsigned int i = 0; i < r.size(); i++) {
    grad = psi.gradient(r[i]);

    for (unsigned int j = 0; j < dim; j++)
      values[i][j] = -grad[j];
  }
}

template<int dim>
void
ExactSolutionSLDI_B<dim>::vector_value_list(
  const std::vector<Point<dim>>& r,
  std::vector<Vector<double>>& values) const
{
  Assert(values.size() == r.size(),
         ExcDimensionMismatch(values.size(), r.size()));

  Tensor<1, dim> grad;
  double mu;

  for (unsigned int i = 0; i < r.size(); i++) {
    if (r[i].norm() < a) {
      // Inside the shield.
      mu = mu_1;
    } else if (r[i].norm() < b) {
      // In the wall of the shield.
      mu = mu_2;
    } else {
      // Outside the shield.
      mu = mu_3;
    }

    grad = psi.gradient(r[i]);

    for (unsigned int j = 0; j < dim; j++)
      values[i][j] = -mu * grad[j];
  }
}

#endif

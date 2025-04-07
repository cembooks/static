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

#ifndef ExactSolutionsMMSVTII_H__
#define ExactSolutionsMMSVTII_H__

#include "constants.hpp"
#include "settings.hpp"
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>

#include <cmath>

using namespace dealii;

inline double
permeability(double x, double y, double mu_0)
{
  return mu_0 * (pow(x, 2) + pow(y, 2) + 1.0);
}

inline double
robin_gamma(double x, double y, double mu_0)
{
  return (sqrt(pow(x, 2) + pow(y, 2)) + 2.0) / permeability(x, y, mu_0);
}

inline Tensor<1, 2>
volume_free_current_density(double x, double y, double mu_0, double k)
{
  Tensor<1, 2> J;
  const double mu = permeability(x, y, mu_0);

  J[0] = (1.0 / mu) * ((mu_0 / mu) * (-2.0 * y * (cos(k * x) + cos(k * y))) -
                       k * sin(k * y));
  J[1] = -(1.0 / mu) * ((mu_0 / mu) * (-2.0 * x * (cos(k * x) + cos(k * y))) -
                        k * sin(k * x));

  return J;
}

inline Tensor<1, 2>
magnetic_vector_potential(double x, double y, double k)
{
  Tensor<1, 2> A;

  A[0] = -sin(k * y) / k;
  A[1] = sin(k * x) / k;

  return A;
}

inline double
magnetic_field(double x, double y, double k)
{
  return (cos(k * x) + cos(k * y));
}

inline double
current_vector_potential(double x, double y, double mu_0, double k)
{
  double T;
  const double mu = permeability(x, y, mu_0);
  const double B = magnetic_field(x, y, k);

  T = B / mu;

  return T;
}

/**
 * \brief Describes exact solution, \f$ T \f$, of the
 * *Method of manufactured solutions, vector potential*
 * [(mms-vt-ii/)](@ref page_mms_vt_ii) numerical experiment.
 *****************************************************************************/
class ExactSolutionMMSVTII_T
  : public Function<2>
  , public SettingsMMSVTII
{
public:
  ExactSolutionMMSVTII_T(){};

  virtual double value(const Point<2>& r,
                       const unsigned int component = 0) const override final;

  virtual Tensor<1, 2> gradient(
    const Point<2>& r,
    const unsigned int component = 0) const override final;
};

/**
 * \brief Describes exact solution, \f$\vec{J}_{f}\f$, of the
 * *Method of manufactured solutions, vector potential*
 * [(mms-vt-ii/)](@ref page_mms_vt_ii) numerical experiment.
 *****************************************************************************/
class ExactSolutionMMSVTII_Jf
  : public Function<2>
  , public SettingsMMSVTII
{
public:
  ExactSolutionMMSVTII_Jf()
    : Function<2>(2)
  {
  }

  virtual void vector_value_list(
    const std::vector<Point<2>>& r,
    std::vector<Vector<double>>& values) const override final;
};

/**
 * \brief Describes exact solution, \f$ B \f$, of the
 * *Method of manufactured solutions, vector potential*
 * [(mms-vt-ii/)](@ref page_mms_vt_ii) numerical experiment.
 *****************************************************************************/
class ExactSolutionMMSVTII_B
  : public Function<2>
  , public SettingsMMSVTII
{
public:
  ExactSolutionMMSVTII_B() {}

  virtual void value_list(const std::vector<Point<2>>& r,
                          std::vector<double>& values,
                          const unsigned int component) const override final;
};

/**
 * \brief Describes the Dirichlet boundary condition for \f$ T \f$, in the
 * *Method of manufactured solutions, vector potential*
 * [(mms-vt-ii/)](@ref page_mms_vt_ii) numerical experiment.
 *****************************************************************************/
class DirichletBC_MMSVTII_T
  : public Function<2>
  , public SettingsMMSVTII
{
public:
  DirichletBC_MMSVTII_T(){};

  virtual void value_list(const std::vector<Point<2>>& r,
                          std::vector<double>& values,
                          const unsigned int component) const override final;
};

/**
 * \brief Describes the Dirichlet boundary condition for \f$\vec{A}\f$, in the
 * *Method of manufactured solutions, vector potential*
 * [(mms-vt-ii/)](@ref page_mms_vt_ii) numerical experiment.
 *****************************************************************************/
class DirichletBC_MMSVTII_A
  : public Function<2>
  , public SettingsMMSVTII
{
public:
  DirichletBC_MMSVTII_A()
    : Function<2>(2)
  {
  }

  virtual void vector_value_list(
    const std::vector<Point<2>>& r,
    std::vector<Vector<double>>& values) const override final;
};

#endif

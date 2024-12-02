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

#ifndef ExactSolutions_H__
#define ExactSolutions_H__

#include "constants.hpp"
#include "settings.hpp"
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>

using namespace dealii;

/**
 * \brief Describes exact solutions, \f$\vec{A}\f$, of the
 * [Method of manufactured solutions, vector potential (mms-v/)](@ref
 *page_mms_v) numerical experiment.
 *****************************************************************************/
template<int dim>
class ExactSolutionMMSV_A
  : public Function<dim>
  , public SettingsMMSV
{
public:
  ExactSolutionMMSV_A();

  virtual void vector_value_list(
    const std::vector<Point<dim>>& r,
    std::vector<Vector<double>>& values) const override final;
};

/**
 * \brief Describes exact solutions, \f$\vec{B}\f$, of the
 * [Method of manufactured solutions, vector potential (mms-v/)](@ref
 *page_mms_v) numerical experiment.
 *****************************************************************************/
template<int dim>
class ExactSolutionMMSV_B
  : public Function<dim>
  , public SettingsMMSV
{
public:
  ExactSolutionMMSV_B();

  // To be used only if dim = 3. Magnetic field in three dimensions is a vector.
  virtual void vector_value_list(
    const std::vector<Point<dim>>& r,
    std::vector<Vector<double>>& values) const override final;

  // To be used only if dim = 2. Magnetic field in two dimensions is a scalar,
  // i.e., out-of-plane vector.
  virtual double value(const Point<dim>& r,
                       const unsigned int component = 0) const override final;
};

#endif

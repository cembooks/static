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

#ifndef ExactSolutionsCVPI_H__
#define ExactSolutionsCVPI_H__

#include "constants.hpp"
#include "settings.hpp"
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>

#include <cmath>

using namespace dealii;

/**
 * \brief Describes the given volume free-current density, \f$\vec{J}_f\f$,
 * in the *Current vector potential* [(cvp-i/)](@ref page_cvp_i)
 * numerical experiment.
 *****************************************************************************/
class ExactSolutionCVPI_Jf
  : public Function<3>
  , public SettingsCVPI
{
public:
  ExactSolutionCVPI_Jf();

  virtual void vector_value_list(
    const std::vector<Point<3>>& r,
    std::vector<Vector<double>>& values) const override final;
};

/**
 * \brief Describes the Dirichlet boundary condition
 * in the *Current vector potential* [(cvp-i/)](@ref page_cvp_i)
 * numerical experiment.
 *****************************************************************************/
class DirichletBC_CVPI
  : public Function<3>
  , public SettingsCVPI
{
public:
  DirichletBC_CVPI();

  virtual void vector_value_list(
    const std::vector<Point<3>>& r,
    std::vector<Vector<double>>& values) const override final;
};
#endif

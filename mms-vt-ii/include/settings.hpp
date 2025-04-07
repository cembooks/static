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

#ifndef SettingsMMSVTII_H__
#define SettingsMMSVTII_H__

#include "constants.hpp"
#include <deal.II/base/types.h>

using namespace dealii;

/**
 * \brief Global settings for the
 * *Method of manufactured solutions, vector potential*
 * [(mms-vt-ii/)](@ref page_mms_vt_ii) numerical experiment.
 *****************************************************************************/
class SettingsMMSVTII : public Constants::Physics
{
public:
  SettingsMMSVTII(){};

  /**
   * \brief If greater than zero, limits the amount of threads used in the
   * simulations.
   *****************************************************************************/
  const unsigned int nr_threads_max = 0;

  /**
   * \brief The permeability of free space.
   *****************************************************************************/
  const double mu_0 = permeability_fs;

  /**
   * \brief The scaling parameter.
   *****************************************************************************/
  const double k = 1.0 * pi;

  /**
   * \brief The half-side of the square in the middle of the circular mesh.
   *****************************************************************************/
  const double d1 = 0.25;

  /**
   * \brief The radius of the circle (sphere) that encloses the square (cube) in
   * the middle of the mesh.
   *****************************************************************************/
  const double rd1 = sqrt(2) * d1;

  /**
   * \brief The Dirichlet boundary condition will be applied to the boundaries
   * with ID = 1.
   *****************************************************************************/
  const types::boundary_id bid_dirichlet = 1;

  /**
   * \brief The Robin boundary condition will be applied to the boundaries
   * with ID = 2.
   *****************************************************************************/
  const types::boundary_id bid_robin = 2;

  /**
   * \brief Two values in double format are considered to be equal if the
   * absolute value of their difference is less than eps.
   *****************************************************************************/
  const double eps = 1e-12;

  /**
   * \brief If set to true, the program will print time tables on the
   * screen.
   *****************************************************************************/
  const bool print_time_tables = false;

  /**
   * \brief If set to true, the program will project the exact solution.
   *
   * The exact solution will be modeled on the same mesh and by the same finite
   * elements that are used to model the solution. The projected exact solution
   * will be saved in the vtu file next to the solution. This option can be
   * useful when debugging.
   *****************************************************************************/
  const bool project_exact_solution = true;

  /**
   * \brief If set to true, saves the residual at each iteration of the
   * CG solver. The names of the files fit the following wildcard
   * *_cg_convergence.csv
   *
   * The residuals are saved into the subdirectories of ./ Data/ directory.
   *****************************************************************************/
  const bool log_cg_convergence = false;
};

using Settings = SettingsMMSVTII;

#endif

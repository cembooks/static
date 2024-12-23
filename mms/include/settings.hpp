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

#ifndef SettingsMMS_H__
#define SettingsMMS_H__

#include "constants.hpp"
#include <deal.II/base/types.h>

using namespace dealii;

/**
 * \brief Global settings for the Method of manufactured solutions
 * [(mms/)](@ref page_mms) numerical experiment.
 *****************************************************************************/
class SettingsMMS : public Constants::Physics
{
public:
  SettingsMMS(){};

  /**
   * \brief If greater than zero, limits the amount of threads used in the
   * simulations.
   *****************************************************************************/
  const unsigned int nr_threads_max = 8;

  /**
   * \brief The permittivity of free space.
   *****************************************************************************/
  const double ep_0 = permittivity_fs;

  /**
   * \brief The scaling parameter.
   *****************************************************************************/
  const double k = pi;

  /**
   * \brief The Dirichlet boundary condition will be applied to the boundaries
   * with ID = 1. The boundary IDs are set in the geo files that are located in
   * the mms/gmsh directory.
   *****************************************************************************/
  const types::boundary_id bid_dirichlet = 1;

  /**
   * \brief The Robin boundary condition will be applied to the boundaries with
   * even boundary IDs. This constant is not used in the code. The code applies
   * the Robin boundary condition on all boundaries with even IDs. The boundary
   * IDs are set in the geo files that are located in the mms/gmsh directory.
   *****************************************************************************/
  const types::boundary_id bid_robin = 2;

  /**
   * \brief Two values in double format are considered to be equal if the
   * absolute value of their difference is less than <code>eps</code>.
   *****************************************************************************/
  const double eps = 1e-12;

  /**
   * \brief If set to true, the program will print the time tables on the
   * screen.
   *****************************************************************************/
  const bool print_time_tables = false;

  /**
   * \brief If set to true, the program will project the exact solution.
   *
   * The exact solution will be modeled on the same mesh and by the same finite
   * elements that are used to model the solution. The projected exact solution
   * will be saved in the vtk file next to the solution. This option can be
   * useful when debugging.
   *****************************************************************************/
  const bool project_exact_solution = false;

  /**
   * \brief If set to true, saves the residual at each iteration of the
   * CG solver. The names of the files fit the following wildcard
   * *_cg_convergence.csv
   *
   * The residuals are saved into the subdirectories of ./ Data/ directory.
   *****************************************************************************/
  const bool log_cg_convergence = false;
};

using Settings = SettingsMMS;

#endif

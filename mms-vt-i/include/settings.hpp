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

#ifndef SettingsMMSVTI_H__
#define SettingsMMSVTI_H__

#include "constants.hpp"
#include <deal.II/base/types.h>

using namespace dealii;

/**
 * \brief Global settings for the
 * [Method of manufactured solutions, vector potential (mms-vt-i/)](@ref
 *page_mms_vt_i) numerical experiment.
 *****************************************************************************/
class SettingsMMSVTI : public Constants::Physics
{
public:
  SettingsMMSVTI(){};

  /**
   * \brief If greater than zero, limits the amount of threads used in the
   * simulations.
   *****************************************************************************/
  const unsigned int nr_threads_max = 8;

  /**
   * \brief The permeability of free space.
   *****************************************************************************/
  const double mu_0 = permeability_fs;

  /**
   * \brief The scaling parameter.
   *****************************************************************************/
  const double k = 1.0 * pi;

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

using Settings = SettingsMMSVTI;

#endif

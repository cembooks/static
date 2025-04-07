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

#ifndef SettingsABC_H__
#define SettingsABC_H__

#include "constants.hpp"
#include <deal.II/base/types.h>

using namespace dealii;

/**
 * \brief Global settings for the
 * *Asymptotic boundary conditions* [(abc/)](@ref page_abc)
 * numerical experiment.
 *****************************************************************************/
class SettingsABC : public Constants::Physics
{
public:
  SettingsABC(){};

  /**
   * \brief If greater than zero, limits the amount of threads used in the
   * simulations.
   *****************************************************************************/
  const unsigned int nr_threads_max = 0;

  /**
   * \brief The permittivity of free space.
   *****************************************************************************/
  const double ep_0 = permittivity_fs;

  /**
   * \brief The radius of the sphere. Used only in the three- dimensional
   * version of the experiment.
   *****************************************************************************/
  const double a = 0.1;

  /**
   * \brief The offset of the conductors from the origin. Used only in the
   * two- dimensional version of the experiment.
   *****************************************************************************/
  const double x0 = 0.1;

  /**
   * \brief The radius of each conductor. Used in the two- dimensional version
   * of the program.
   *****************************************************************************/
  const double R = 0.05;

  /**
   * \brief The outer radius of the [fixed region](@ref abc_fixed_reg)
   * in the three- dimensional version of the problem.
   *****************************************************************************/
  const double R_mid = 1.5 * (R + x0);

  /**
   * \brief The ID of the circular boundary that represents the left
   * conductor in the two- dimensional version of the experiment.
   *****************************************************************************/
  const types::boundary_id bid_left = 3;

  /**
   * \brief The ID of the circular boundary that represents the right
   * conductor in the two- dimensional version of the experiment.
   *****************************************************************************/
  const types::boundary_id bid_right = 1;

  /**
   * \brief The ID of the inner boundary of the problem domain in the
   * three- dimensional version of the experiment.
   *****************************************************************************/
  const types::boundary_id bid_in = 1;

  /**
   * \brief The ID of the spherical manifold attached to the left boundary
   * (2D only).
   *****************************************************************************/
  const types::manifold_id mfid_left = 3;

  /**
   * \brief The ID of the spherical manifold attached to the right boundary
   * (2D only).
   *****************************************************************************/
  const types::manifold_id mfid_right = 1;

  /**
   * \brief The ID of the spherical manifold attached to the boundary that
   * represents the infinity (both, 2D and 3D).
   *****************************************************************************/
  const types::manifold_id mfid_infty = 2;

  /**
   * \brief The ID of the boundary that represents the infinity. Used in both,
   * two- and three- dimensional versions of the experiment.
   *****************************************************************************/

#if BC_TYPE__ == 0
  const types::boundary_id bid_infty = 0; // Neumann boundary condition on the
                                          // outer boundary.
#endif

#if BC_TYPE__ == 1
  const types::boundary_id bid_infty = 5; // Dirichlet boundary condition on the
                                          // outer boundary.
#endif

#if BC_TYPE__ == 2
  const types::boundary_id bid_infty = 2; // ABC on the outer boundary.
#endif

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

using Settings = SettingsABC;

#endif

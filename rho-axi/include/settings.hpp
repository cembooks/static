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

#ifndef SettingsRHOAXI_H__
#define SettingsRHOAXI_H__

#include "constants.hpp"
#include <deal.II/base/types.h>

using namespace dealii;

/**
 * \brief Global settings for the
 * *Axisymmetric - volume charge* [(rho-axi/)](@ref page_rho_axi)
 * numerical experiment.
 *****************************************************************************/
class SettingsRHOAXI : public Constants::Physics
{
public:
  SettingsRHOAXI(){};

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
   * \brief The volume free-charge density.
   *****************************************************************************/
  const double rho = 1.0 * ep_0;

  /**
   * \brief The width of the rectangle around the origin in the second mesh
   * (half-circle).
   *****************************************************************************/
  const double d = 0.2;

  /**
   * \brief The radius of the circle that encloses the rectangle in the middle
   * of the mesh.
   *****************************************************************************/
  const double rd = sqrt(2) * d;

  /**
   * \brief The radius of the tube in the two- dimensional version of the
   * problem or the radius of the sphere in the three- dimensional problem
   * in which the free- charge is located.
   *****************************************************************************/
  double a = 0.5;

  /**
   * \brief The radius of the outer boundary of the problem domain.
   *****************************************************************************/
  double b = 1.0;

  /**
   * \brief The ID of the material inside the charged region,
   * \f$\rho_f \ne 0.0\f$.
   *****************************************************************************/
  const types::material_id mid_1 = 1;

  /**
   * \brief The ID of the material outside the charged region, \f$\rho_f=0.0\f$.
   *****************************************************************************/
  const types::material_id mid_2 = 2;

  /**
   * \brief The ID of the boundary of the problem domain.
   *****************************************************************************/
  const types::boundary_id bid = 1;

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

using Settings = SettingsRHOAXI;

#endif

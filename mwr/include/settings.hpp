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

#ifndef SettingsMWR_H__
#define SettingsMWR_H__

#include "constants.hpp"
#include <deal.II/base/types.h>

using namespace dealii;

/**
 * \brief Global settings for the
 * *Magnetic wire* [(mwr/)](@ref page_mwr)
 * numerical experiment.
 *****************************************************************************/
class SettingsMWR : public Constants::Physics
{
public:
  SettingsMWR(){};

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
   * \brief The current density.
   *****************************************************************************/
  const double J_f = 1.0;

  /**
   * \brief The relative permeability of the material of the wire.
   *****************************************************************************/
  const double mu_r = 4.0;

  /**
   * \brief The permeability of the material of the wire.
   *****************************************************************************/
  const double mu = mu_r * mu_0;

  /**
   * \brief The half-side of the square in the middle of the mesh.
   *****************************************************************************/
  const double d1 = 0.2;

  /**
   * \brief The radius of the circle that encloses the square in the middle of
   * the mesh.
   *****************************************************************************/
  const double rd1 = sqrt(2) * d1;

  /**
   * \brief The radius of the wire.
   *****************************************************************************/
  const double a = 0.5;

  /**
   * \brief The radius of the outer boundary of the problem domain.
   *****************************************************************************/
  const double b = 1.0;

  /**
   * \brief A constant used to compute the exact solution.
   *****************************************************************************/
  const double A_0 =
    0.25 * mu_0 * J_f * a * a * log(exp(mu_r) * (b * b) / (a * a));

  /**
   * \brief The ID of the material inside the wire.
   *****************************************************************************/
  const types::material_id mid_1 = 1;

  /**
   * \brief The ID of the material outside the wire.
   *****************************************************************************/
  const types::material_id mid_2 = 2;

  /**
   * \brief The ID of the boundary of the problem domain. The boundary ID is set
   * in the geo files that are located in the mwr/gmsh directory.
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
   * The exact solutions will be modeled on the same mesh and by the same finite
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

using Settings = SettingsMWR;

#endif

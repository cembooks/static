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

#ifndef SettingsSSOLI_H__
#define SettingsSSOLI_H__

#include <deal.II/base/types.h>
#include "constants.hpp"

using namespace dealii;

/**
 * \brief Global settings for the
 * [Thin spherical coil (ssol-i/)](@ref page_ssol_i)
 * numerical experiment
 *****************************************************************************/
class SettingsSSOLI : public Constants::Physics
{
public:
	SettingsSSOLI() {};

	/**
	 * \brief If greater than zero, limits the amount of threads used in the
	 * simulations.
	 *************************************************************************/
	const unsigned int nr_threads_max = 8;

/**
 * \brief The permeability of free space.
 *****************************************************************************/
	const double mu_0 = permeability_fs;

/**
 * \brief The radius of the coil.
 *****************************************************************************/
	const double a = 0.5;

/**
 * \brief A constant that defines the magnitude of the surface free-current 
 * density.
 *****************************************************************************/
	const double K_0 = 1.0;

/**
 * \brief Two values in double format are considered to be equal if the
 * absolute value of their difference is less than **eps**.
 *****************************************************************************/
	const double eps = 1e-12;

/**
 * \brief If set to true, the program will print the time tables on the
 * screen.
 *****************************************************************************/
	const bool print_time_tables = false; // true;

/**
 * \brief If set to true, the program will project the exact solution.
 *
 * The exact solutions will be modeled on the same mesh
 * and by the same finite elements that are used to model the solution.
 * The projected exact solution will be saved in the vtk file next to the
 * solution. This option can be useful when debugging.
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

using Settings = SettingsSSOLI;

#endif


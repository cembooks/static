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

#ifndef SettingsPLS_H__

#define SettingsPLS_H__

#include <deal.II/base/types.h>
#include "constants.hpp"

using namespace dealii;

/**
 * \brief Global settings for the
 * [Planes of symmetry](@ref page_pls)
 * numerical experiment
 *****************************************************************************/
class SettingsPLS : public Constants::Physics
{
public:
	SettingsPLS() {};

/**
 * \brief If greater than zero, limits the amount of threads used in the
 * simulations.
 *****************************************************************************/
	const unsigned int nr_threads_max = 8;

/**
 * \brief The permittivity of free space.
 *
 * This variable, ep_0, is used throughout the program,
 * not permittivity_fs. So one can scale the system of linear equation
 * differently by setting ep_0 = 1.0.
 *****************************************************************************/
	const double ep_0 = permittivity_fs;

/**
 * \brief The [scaling parameter](@ref mms_k_parameter).
 *****************************************************************************/
	const double k = pi;

/**
 * \brief The [Dirichlet boundary condition](@ref mms_bcs) will be applied
 * to the boundaries marked by ID = 1.
 *****************************************************************************/
	const types::boundary_id bid_dirichlet = 1;

/**
 * \brief The [Robin boundary condition](@ref mms_bcs) will be applied
 * to the boundaries marked by ID = 2.
 *****************************************************************************/
	const types::boundary_id bid_robin = 2;

/**
 * \brief Two values in double format are considered to be equal if the
 * absolute value of their difference is less than eps.
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
 * The exact solution, i.e., \f$\Phi\f$, will be modeled on the same mesh
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

using Settings = SettingsPLS;

#endif


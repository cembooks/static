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

#ifndef SettingsSCHAXI_H__
#define SettingsSCHAXI_H__

#include <deal.II/base/types.h>
#include "constants.hpp"

using namespace dealii;

/**
 * \brief Global settings for the
 * [Axisymmetric - surface charge (sch-axi/)](@ref page_sch_axi)
 * numerical experiment.
 *****************************************************************************/
class SettingsSCHAXI : public Constants::Physics
{
public:
	SettingsSCHAXI() {};

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
 * \brief The radius of the tube in the two- dimensional version of the
 * problem or the radius of the sphere in the three- dimensional problem
 * on which the free- charge is located.
 *****************************************************************************/
	double a = 0.5;

/**
 * \brief The radius of the outer boundary of the problem domain.
 *****************************************************************************/
	double b = 1.0;

/**
 * \brief The surface free-charge density.
 *****************************************************************************/
#if CYLINDER__ == 1
	const double kappa_f = ep_0 / (a * std::log(b / a));
#endif

#if CYLINDER__ == 0
	const double kappa_f = ep_0 * b / (a*(b-a));
#endif

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

using Settings = SettingsSCHAXI;

#endif


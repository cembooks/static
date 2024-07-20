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

#ifndef SettingsFLCAXI_H__
#define SettingsFLCAXI_H__

#include <deal.II/base/types.h>
#include "constants.hpp"

using namespace dealii;

/**
 * \brief Global settings for the
 * [Axisymmetric - floating conductor](@ref page_flc_axi)
 * numerical experiment
 *****************************************************************************/
class SettingsFLCAXI : public Constants::Physics
{
public:
	SettingsFLCAXI() {};

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
 * \brief The radius of the inner boundary of the problem domain.
 *****************************************************************************/
	double a = 0.4;

/**
 * \brief The radius of the outer boundary of the problem domain.
 *****************************************************************************/
	double b = 1.0;

/**
 * \brief The radius of the inner interface between dissimilar materials.
 *****************************************************************************/
	double d_1 = 0.6;

/**
 * \brief The radius of the outer interface between dissimilar materials.
 *****************************************************************************/
	double d_2 = 0.8;

/**
 * \brief The ID of the inner boundary of the problem domain.
 *****************************************************************************/
	const types::boundary_id bid_in = 1;

/**
 * \brief The ID of the outer boundary of the problem domain.
 *****************************************************************************/
	const types::boundary_id bid_out = 3;

/**
 * \brief The material ID of the inner dielectric tube.
 *****************************************************************************/
	const types::material_id mid_1 = 1;

/**
 * \brief The material ID of the outer dielectric tube.
 *****************************************************************************/
	const types::material_id mid_2 = 2;

/**
 * \brief The material ID of the middle dielectric tube that represents the
 * floating conductor.
 *****************************************************************************/
	const types::material_id mid_3 = 3;

/**
 * \brief Relative permittivity of the inner dielectric.
 *****************************************************************************/
	const double ep_1 = 32.0*ep_0;

/**
 * \brief Relative permittivity of the outer dielectric.
 *****************************************************************************/
	const double ep_2 = 4.0*ep_0;

/**
 * \brief Relative permittivity of the middle dielectric that replaces the
 * floating conductor.
 *****************************************************************************/
	const double ep_3 = 1e9*ep_0;

/**
 * \brief If greater than zero, limits the amount of threads used in the
 * simulations.
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

using Settings = SettingsFLCAXI;

#endif

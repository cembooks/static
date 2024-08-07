/*****************************************************************************
* Copyright (C) Siarhei Uzunbajakau, 2023.
*
* This program is free software. You can use, modify, and redistribute it under
* the terms of the GNU Lesser General Public License as published by the Free
* Software Foundation, either version 3 or (at your option) any later version.
* This program is distributed without any warranty.
*
* Refer to COPYING.LESSER for more details.
******************************************************************************/

#ifndef SettingsSLDI_H__
#define SettingsSLDI_H__

#include <deal.II/base/types.h>
#include "constants.hpp"

using namespace dealii;

/**
 * \brief Global settings for the
 * [Magnetostatic shield - 1 (sld-i/)](@ref page_sld_i)
 * numerical experiment.
 *****************************************************************************/
class SettingsSLDI : public Constants::Physics
{
public:
	SettingsSLDI() {};

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
 * \brief The inner radius of the shield.
 *****************************************************************************/
	const double a = 0.2;

/**
 * \brief The outer radius of the shield.
 *****************************************************************************/
	const double b = 0.4;

/**
 * \brief The half- side length of the square (cube) in which the error norms 
 * are computed.
 *****************************************************************************/
	const double d_2 = 0.8;

/**
 * \brief The half- side length of the square (cube) that represents the outer
 * boundary.
 *****************************************************************************/
	const double d_3 = 2.0;

/**
 * \brief The ID of the only boundary of the problem domain.
 *****************************************************************************/
	const types::boundary_id bid = 1;

/**
 * \brief The ID of the material inside the shield.
 *****************************************************************************/
	const types::material_id mid_1 = 1;

/**
 * \brief The ID of the material of the shield.
 *****************************************************************************/
	const types::material_id mid_2 = 2;

/**
 * \brief The ID of the material outside the shield.
 *****************************************************************************/
	const types::material_id mid_3 = 3;

/**
 * \brief Relative permeability of the material inside the shield.
 *****************************************************************************/
	const double mur_1 = 1.0;

/**
 * \brief Relative permeability of the material of the shield.
 *****************************************************************************/
	const double mur_2 = 4.0;

/**
 * \brief Relative permeability of the material outside the shield.
 *****************************************************************************/
	const double mur_3 = 1.0;

/**
 * \brief Permeability of the material inside the shield.
 * simulations.
 *****************************************************************************/
	const double mu_1 = mur_1*mu_0;

/**
 * \brief Permeability of the material of the shield.
 *****************************************************************************/
	const double mu_2 = mur_2*mu_0;

/**
 * \brief Permeability of the material outside the shield.
 *****************************************************************************/
	const double mu_3 = mur_3*mu_0;

/**
 * \brief The magnitude of the uniform auxiliary field H at the infinity,
 * i.e., in absence of the magnetic shield.
 *****************************************************************************/
	const double H_0 = 1.0;

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
 * The exact solutions will be modeled on the same mesh and by the same finite
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

using Settings = SettingsSLDI;

#endif


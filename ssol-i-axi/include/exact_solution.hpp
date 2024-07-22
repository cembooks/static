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

#ifndef ExactSolutionSSOLIAXI_H__
#define ExactSolutionSSOLIAXI_H__

#include <deal.II/base/function.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/vector.h>

#include <cmath>

#include "constants.hpp"
#include "settings.hpp"

using namespace dealii;

/**
 * \brief Describes exact solution, \f$ A' \f$, of the
 * [Axisymmetric - thin spherical coil (ssol-i-axi)](@ref page_ssol_i_axi)
 * numerical experiment.
 *****************************************************************************/
class ExactSolutionSSOLIAXI_A : public Function<2>, public SettingsSSOLIAXI
{
public:

	ExactSolutionSSOLIAXI_A(){};

	virtual double value(const Point<2> & p,
		const unsigned int component = 0) const override final;

	virtual Tensor<1, 2> gradient(const Point<2> & p,
				const unsigned int component = 0) const override final;
};

/**
 * \brief Describes exact solution, \f$\vec{B}'\f$, of the
 * [Axisymmetric - thin spherical coil (ssol-i-axi)](@ref page_ssol_i_axi)
 * numerical experiment.
 *****************************************************************************/
class ExactSolutionSSOLIAXI_B : public Function<2>, public SettingsSSOLIAXI
{
public:

	ExactSolutionSSOLIAXI_B(): Function<2>(2), B_0(2.0*mu_0*K_0/3.0) {}

	virtual void vector_value_list(const std::vector<Point<2>> & r,
		std::vector<Vector<double>> &values) const final;

private:
	const double B_0;
};

/**
 * \brief Describes exact solution, \f$\vec{H}'\f$, of the
 * [Axisymmetric - thin spherical coil (ssol-i-axi)](@ref page_ssol_i_axi)
 * numerical experiment.
 *****************************************************************************/
class ExactSolutionSSOLIAXI_H: public Function<2>, public SettingsSSOLIAXI
{
public:

	ExactSolutionSSOLIAXI_H(): Function<2>(2) {}

	virtual void vector_value_list(const std::vector<Point<2>> & r,
		std::vector<Vector<double>> &values) const final;

private:

	ExactSolutionSSOLIAXI_B B;
};

#endif


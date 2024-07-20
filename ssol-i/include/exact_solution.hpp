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

#ifndef ExactSolutionsSSOLI_H__
#define ExactSolutionsSSOLI_H__

#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include "constants.hpp"
#include "settings.hpp"

using namespace dealii;

/**
 * \brief Describes exact solution, \f$\vec{B}\f$, of the
 * [Thin spherical coil (ssol-i/)](@ref page_ssol_i)
 * numerical experiment.
 *****************************************************************************/
class ExactSolutionSSOLI_B : public Function<3>, public SettingsSSOLI
{
public:

	ExactSolutionSSOLI_B();

	virtual void vector_value_list(const std::vector<Point<3>> & r,
		std::vector<Vector<double>>	 &values) const override final;

private:
	const double B_0;
};

#endif


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

#ifndef ExactSolutions_H__
#define ExactSolutions_H__

#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include "constants.hpp"
#include "settings.hpp"

#include <cmath>

using namespace dealii;

/**
 * \brief The exact solution, \f$T\f$, in the
 * [Current vector potential (cvp-ii/)](@ref page_cvp_ii)
 * numerical experiment.
 *****************************************************************************/
class ExactSolutionCVPII_T : public Function<2>, public SettingsCVPII
{
public:

	ExactSolutionCVPII_T() {};

	virtual double value(const Point<2> & r,
		const unsigned int component = 0) const override final;

	virtual Tensor<1, 2>
		gradient(const Point<2> & r,
			const unsigned int component = 0) const override final;
};

/**
 * \brief Describes the given volume free-current density, \f$\vec{J}_f\f$,
 * in the
 * [Current vector potential (cvp-ii/)](@ref page_cvp_ii)
 * numerical experiment.
 *****************************************************************************/
class ExactSolutionCVPII_Jf : public Function<2>, public SettingsCVPII
{
public:

	ExactSolutionCVPII_Jf();

	virtual void vector_value_list(const std::vector<Point<2>> & r,
		std::vector<Vector<double>>	 &values) const override final;
};

#endif


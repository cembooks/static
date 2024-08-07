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

#ifndef ExactSolutionsABC_H__
#define ExactSolutionsABC_H__

#include <deal.II/base/point.h>
#include <deal.II/base/function.h>

#include <cmath>

#include "constants.hpp"
#include "settings.hpp"

using namespace dealii;

/**
 * \brief Describes exact solution, \f$\Phi\f$, of the
 * [Asymptotic boundary condition (abc/)](@ref page_abc)
 * numerical experiment.
 *****************************************************************************/
template<int dim>
class ExactSolutionABC_PHI : public Function<dim>, public SettingsABC
{
public:

	ExactSolutionABC_PHI();

	virtual double value(const Point<dim> & r,
		const unsigned int component = 0) const override final;

	virtual Tensor<1, dim>
		gradient(const Point<dim> & r,
			const unsigned int component = 0) const override final;

private:
	// Half- separation between the image lines of charge.
	double d = 0.0;

	// Charge density of the image lines of charge.
	double lambda = 0.0;
};

#endif


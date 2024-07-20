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

#ifndef ExactSolutionCBND_H__
#define ExactSolutionCBND_H__

#include <deal.II/base/point.h>
#include <deal.II/base/function.h>

#include <cmath>

#include "constants.hpp"
#include "settings.hpp"

using namespace dealii;

/**
 * \brief Describes exact solutions, \f$\Phi\f$, of the
 *	[Effect of curved boundaries](@ref page_cbnd)
 *	numerical experiment in
 *	[two](@ref cbnd_exact_solution_2D)
 *	and
 *	[three](@ref cbnd_exact_solution_3D)
 *	dimensions.
 *****************************************************************************/
template<int dim>
class ExactSolutionCBND_PHI : public Function<dim>, public SettingsCBND
{
public:

	ExactSolutionCBND_PHI() {};

	virtual double value(const Point<dim> & r,
		const unsigned int component = 0) const override final;

	virtual Tensor<1, dim>
		gradient(const Point<dim> & r,
				const unsigned int component = 0) const override final;
};

#endif


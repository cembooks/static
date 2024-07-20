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

#ifndef ExactSolutionFLC_H__
#define ExactSolutionFLC_H__

#include <deal.II/base/function.h>
#include "constants.hpp"
#include "settings.hpp"

using namespace dealii;

/**
 * \brief Describes exact solutions, \f$\Phi\f$, of the
 * [Floating conductor](@ref page_flc)
 * numerical experiment in
 * [two](@ref flc_exact_solution_2D)
 * and
 * [three](@ref flc_exact_solution_3D)
 * dimensions.
 *****************************************************************************/
template<int dim>
class ExactSolutionFLC_PHI : public Function<dim>, public SettingsFLC
{
public:

	ExactSolutionFLC_PHI();

	virtual double value(const Point<dim> & r,
		const unsigned int component = 0) const override final;

	virtual Tensor<1, dim>
		gradient(const Point<dim> & r,
				const unsigned int component = 0) const override final;

private:

	double alpha;
	double beta;
	double phi_d;
};

#endif

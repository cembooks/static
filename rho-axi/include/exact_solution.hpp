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

#ifndef ExactSolutionRHOAXI_H__
#define ExactSolutionRHOAXI_H__

#include <deal.II/base/function.h>
#include "constants.hpp"
#include "settings.hpp"

using namespace dealii;

/**
 * \brief Describes exact solutions, \f$\Phi\f$, of the
 * [Axisymmetric - volume charge](@ref page_rho_axi)
 * numerical experiment.
 *****************************************************************************/
template<bool is_cylinder>
class ExactSolutionRHOAXI_PHI : public Function<2>, public SettingsRHOAXI
{
public:

	ExactSolutionRHOAXI_PHI(){};

	virtual double value(const Point<2> & p,
		const unsigned int component = 0) const override final;

	virtual Tensor<1, 2>
		gradient(const Point<2> & p,
			const unsigned int component = 0) const override final;
};

#endif

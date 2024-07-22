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

#ifndef ExactSolutionsPLS_H__
#define ExactSolutionsPLS_H__

#include <deal.II/base/function.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/vector.h>

#include <cmath>

#include "constants.hpp"
#include "settings.hpp"

using namespace dealii;

/**
 * \brief Describes exact solution, \f$\Phi\f$, of the
 * [Planes of symmetry (pls/)](@ref page_pls)
 * numerical experiment in two and three dimensions.
 *****************************************************************************/
template<int dim>
class ExactSolutionPLS_PHI : public Function<dim>, public SettingsPLS
{
public:

	ExactSolutionPLS_PHI() {};

	virtual double value(const Point<dim> & r,
		const unsigned int component = 0) const override final;

	virtual Tensor<1, dim>
		gradient(const Point<dim> & r,
				const unsigned int component = 0) const override final;
};

#endif

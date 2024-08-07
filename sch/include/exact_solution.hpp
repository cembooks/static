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

#ifndef ExactSolutionSCH_H__
#define ExactSolutionSCH_H__

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
 * [Surface charge (sch/)](@ref page_sch)
 * numerical experiment.
 *****************************************************************************/
template<int dim>
class ExactSolutionSCH_PHI : public Function<dim>, public SettingsSCH
{
public:

	ExactSolutionSCH_PHI(){};

	virtual double value(const Point<dim> & r,
		const unsigned int component = 0) const override final;

	virtual Tensor<1, dim>
		gradient(const Point<dim> & r,
			const unsigned int component = 0) const override final;
};

#endif


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

#ifndef ExactSolutionFLCAXI_H__
#define ExactSolutionFLCAXI_H__

#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/vector.h>

#include "constants.hpp"
#include "settings.hpp"

using namespace dealii;

/**
 * \brief Describes the exact solution, \f$\Phi\f$, of the
 * [Axisymmetric - floating conductor (flc-axi/)](@ref page_flc_axi)
 * numerical experiment.
 *****************************************************************************/
template<bool is_cylinder>
class ExactSolutionFLCAXI_PHI : public Function<2>, public SettingsFLCAXI
{
public:

	ExactSolutionFLCAXI_PHI();

	virtual double value(const Point<2> & r,
		const unsigned int component = 0) const override final;

	virtual Tensor<1, 2>
		gradient(const Point<2> & r,
				const unsigned int component = 0) const override final;

private:

	double alpha;
	double beta;
	double phi_d;
};

#endif


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

#ifndef ExactSolutionINTAXI_H__
#define ExactSolutionINTAXI_H__

#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/vector.h>

#include <cmath>

#include "constants.hpp"
#include "settings.hpp"

using namespace dealii;

/**
 * \brief Describes exact solution, \f$\Phi\f$, of the
 * [Axisymmetric - interface between dielectrics (int-axi/)](@ref page_int_axi)
 * numerical experiment in two and three dimensions.
 *****************************************************************************/
template<bool is_cylinder>
class ExactSolutionINTAXI_PHI : public Function<2>, public SettingsINTAXI
{
public:

	ExactSolutionINTAXI_PHI();

	virtual double value(const Point<2> & r,
		const unsigned int component = 0) const override final;

	virtual Tensor<1, 2>
		gradient(const Point<2> & r,
				const unsigned int component = 0) const override final;

private:

	double alpha;
	double beta;

};

/**
 * \brief Describes exact solution, \f$\vec{E}\f$, of the
 * [Axisymmetric - interface between dielectrics (int-axi/)](@ref page_int_axi)
 * numerical experiment in two and three dimensions.
 *****************************************************************************/
template<bool is_cylinder>
class ExactSolutionINTAXI_E: public Function<2>, public SettingsINTAXI
{
public:

	ExactSolutionINTAXI_E(): Function<2>(2) {}

	virtual void vector_value_list(const std::vector<Point<2>> & r,
		std::vector<Vector<double>> &values) const final;

private:

	ExactSolutionINTAXI_PHI<is_cylinder> PHI;
};

/**
 * \brief Describes exact solution, \f$\vec{D}\f$, of the
 * [Axisymmetric - interface between dielectrics (int-axi/)](@ref page_int_axi)
 * numerical experiment in two and three dimensions.
 *****************************************************************************/
template<bool is_cylinder>
class ExactSolutionINTAXI_D : public Function<2>, public SettingsINTAXI
{
public:

	ExactSolutionINTAXI_D(): Function<2>(2) {}

	virtual void vector_value_list(const std::vector<Point<2>> & r,
		std::vector<Vector<double>> &values) const final;

private:

	ExactSolutionINTAXI_PHI<is_cylinder> PHI;
};

#endif


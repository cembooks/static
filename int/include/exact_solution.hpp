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

#ifndef ExactSolutionINT_H__
#define ExactSolutionINT_H__

#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/vector.h>

#include <cmath>

#include "constants.hpp"
#include "settings.hpp"

using namespace dealii;

/**
 * \brief Describes exact solutions, \f$\Phi\f$, of the
 * [Interface between dielectrics](@ref page_int)
 * numerical experiment in two and three dimensions.
 *****************************************************************************/
template<int dim>
class ExactSolutionINT_PHI : public Function<dim>, public SettingsINT
{
public:

	ExactSolutionINT_PHI();

	virtual double value(const Point<dim> & r,
		const unsigned int component = 0) const override final;

	virtual Tensor<1, dim>
		gradient(const Point<dim> & r,
				const unsigned int component = 0) const override final;

private:

	double alpha;
	double beta;
};

/**
 * \brief Describes exact solutions, \f$\vec{E}\f$, of the
 * [Interface between dielectrics](@ref page_int)
 * numerical experiment in two and three	dimensions.
 *****************************************************************************/
template<int dim>
class ExactSolutionINT_E: public Function<dim>, public SettingsINT
{
public:

	ExactSolutionINT_E(): Function<dim>(dim) {}

	virtual void vector_value_list(const std::vector<Point<dim>> & r,
		std::vector<Vector<double>> &values) const final;

private:

	ExactSolutionINT_PHI<dim> PHI;
};

/**
 * \brief Describes exact solutions, \f$\vec{D}\f$, of the
 * [Interface between dielectrics](@ref page_int)
 * numerical experiment in two and three	dimensions.
 *****************************************************************************/
template<int dim>
class ExactSolutionINT_D : public Function<dim>, public SettingsINT
{
public:

	ExactSolutionINT_D(): Function<dim>(dim) {}

	virtual void vector_value_list(const std::vector<Point<dim>> & r,
		std::vector<Vector<double>> &values) const final;

private:

	ExactSolutionINT_PHI<dim> PHI;
};

#endif


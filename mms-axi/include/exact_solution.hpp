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

#ifndef ExactSolutionMMSAXI_H__
#define ExactSolutionMMSAXI_H__

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
 * [Axisymmetric - method of manufactured solutions (mms-axi/)](@ref page_mms_axi)
 * numerical experiment in two and three dimensions.
 *****************************************************************************/
template<int dim>
class ExactSolutionMMSAXI_PHI : public Function<dim>, public SettingsMMSAXI
{
public:

	ExactSolutionMMSAXI_PHI() {};

	virtual double value(const Point<dim> & p,
		const unsigned int component = 0) const override final;

	virtual Tensor<1, dim>
		gradient(const Point<dim> & p,
				const unsigned int component = 0) const override final;
};

/**
 * \brief Describes exact solution, \f$\vec{E}\f$, of the
 * [Axisymmetric - method of manufactured solutions (mms-axi/)](@ref page_mms_axi)
 * numerical experiment in two and three dimensions.
 *****************************************************************************/
template<int dim>
class ExactSolutionMMSAXI_E : public Function<dim>, public SettingsMMSAXI
{
public:

	ExactSolutionMMSAXI_E(): Function<dim>(dim) {}

	virtual void vector_value_list(const std::vector<Point<dim>> & p,
		std::vector<Vector<double>>	 &values) const final;
};

/**
 * \brief Describes exact solution, \f$\vec{D}\f$, of the
 * [Axisymmetric - method of manufactured solutions (mms-axi/)](@ref page_mms_axi)
 * numerical experiment in two and three dimensions.
 *****************************************************************************/
template<int dim>
class ExactSolutionMMSAXI_D : public Function<dim>, public SettingsMMSAXI
{
public:

	ExactSolutionMMSAXI_D(): Function<dim>(dim) {}

	virtual void vector_value_list(const std::vector<Point<dim>> & p,
		std::vector<Vector<double>>	 &values) const final;
};

#endif

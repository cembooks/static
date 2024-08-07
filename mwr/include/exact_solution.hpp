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

#ifndef ExactSolutionMWR_H__
#define ExactSolutionMWR_H__

#include <deal.II/base/function.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/vector.h>

#include <cmath>

#include "constants.hpp"
#include "settings.hpp"

using namespace dealii;

/**
 * \brief Describes exact solution, \f$ A \f$, of the
 * [Magnetic wire (mwr/)](@ref page_mwr)
 * numerical experiment.
 *****************************************************************************/
class ExactSolutionMWR_A : public Function<2>, public SettingsMWR
{
public:

	ExactSolutionMWR_A(){};

	virtual double value(const Point<2> & r,
		const unsigned int component = 0) const override final;

	virtual Tensor<1, 2>
		gradient(const Point<2> & r,
			const unsigned int component = 0) const override final;
};

/**
 * \brief Describes exact solution, \f$\vec{H}\f$, of the
 * [Magnetic wire (mwr/)](@ref page_mwr)
 * numerical experiment.
 *****************************************************************************/
class ExactSolutionMWR_H: public Function<2>, public SettingsMWR
{
public:

	ExactSolutionMWR_H(): Function<2>(2) {}

	virtual void vector_value_list(const std::vector<Point<2>> & r,
		std::vector<Vector<double>> &values) const final;

private:

	ExactSolutionMWR_A A;
};

/**
 * \brief Describes exact solution, \f$\vec{B}\f$, of the
 * [Magnetic wire (mwr/)](@ref page_mwr)
 * numerical experiment.
 *****************************************************************************/
class ExactSolutionMWR_B : public Function<2>, public SettingsMWR
{
public:

	ExactSolutionMWR_B(): Function<2>(2) {}

	virtual void vector_value_list(const std::vector<Point<2>> & r,
		std::vector<Vector<double>> &values) const final;

private:

	ExactSolutionMWR_A A;
};

#endif


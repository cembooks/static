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

#include <deal.II/base/types.h>
#define BOOST_ALLOW_DEPRECATED_HEADERS

#include <math.h>
#include "static_scalar_input.hpp"
#include "exact_solution.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"

//-----------------------------------------------------------------------------
//-------- Stage 0. Calculating current vector potential, T, given Jf ---------
//-----------------------------------------------------------------------------

template<>
void StaticScalarSolver::TheCoefficient<2,0>::value_list(
	const std::vector<Point<2>> &r,
	types::material_id mid,
	unsigned int cuid,
	std::vector<double> & values) const
{
	Assert(r.size() == values.size(),
		ExcDimensionMismatch(r.size(), values.size()));

	for (unsigned int i = 0 ; i < values.size(); i++)
		values[i] = 1.0;
}

template<>
void StaticScalarSolver::TheCoefficient<2,1>::value_list(
	const std::vector<Point<2>> &r,
	types::material_id mid,
	unsigned int cuid,
	std::vector<double> & values) const
{
	Assert(r.size() == values.size(),
		ExcDimensionMismatch(r.size(), values.size()));

	for (unsigned int i = 0 ; i < values.size(); i++)
		values[i] = 1.0;
}

template<>
void StaticScalarSolver::PdeRhs<2,0>::value_list(
	const std::vector<Point<2>> &r,
	types::material_id mid,
	unsigned int cuid,
	std::vector<double> & values) const
{
	Assert(r.size() == values.size(),
		ExcDimensionMismatch(r.size(), values.size()));

	for (unsigned int i = 0 ; i < values.size(); i++)
		values[i] = 0.0;
}

template<>
void StaticScalarSolver::PdeRhsCvp<2,0>::value_list(
	const std::vector<Point<2>> &r,
	types::material_id mid,
	unsigned int cuid,
	std::vector<Tensor<1, 2>> & values) const
{
	Assert(r.size() == values.size(),
		ExcDimensionMismatch(r.size(), values.size()));

	Tensor<1,2> Jf;

	for (unsigned int i = 0; i < values.size(); i++)
	{
		Jf = volume_free_current_density(r.at(i)[0],r.at(i)[1],mu_0,k);

		values.at(i)[0] = Jf[0];
		values.at(i)[1] = Jf[1];
	}
}

template<>
void StaticScalarSolver::Gamma<2,0>::value_list(
	const std::vector<Point<2>> &r,
	const std::vector<Tensor<1, 2>> & n,
	types::boundary_id bid,
	types::material_id mid,
	unsigned int cuid,
	unsigned int fuid,
	std::vector<double> & values) const
{
	Assert(r.size() == values.size(),
		ExcDimensionMismatch(r.size(), values.size()));

	Assert(r.size() == n.size(),
		ExcDimensionMismatch(r.size(), n.size()));

	for (unsigned int i = 0 ; i < values.size(); i++)
		values[i] = 0;
}

template<>
void StaticScalarSolver::RobinRhs<2,0>::value_list(
	const std::vector<Point<2>> &r,
	const std::vector<Tensor<1, 2>> & n,
	types::boundary_id bid,
	types::material_id mid,
	unsigned int cuid,
	unsigned int fuid,
	std::vector<double> & values) const
{

	Assert(r.size() == values.size(),
		ExcDimensionMismatch(r.size(), values.size()));

	Tensor<1,2> Jf;

	for (unsigned int i = 0; i < r.size(); i++)
	{
		Jf = volume_free_current_density(r.at(i)[0],r.at(i)[1],mu_0,k);

		values.at(i) = -n.at(i)[0]*Jf[1] + n.at(i)[1]*Jf[0];
	}
}

template<>
void StaticScalarSolver::FreeSurfaceCharge<2,0>::value_list(
	const std::vector<Point<2>> &r,
	const std::vector<Tensor<1,2>> & n,
	types::material_id mid,
	unsigned int cuid,
	unsigned int fuid,
	std::vector<double> & values) const
{
	Assert(r.size() == values.size(),
		ExcDimensionMismatch(r.size(), values.size()));

	for (unsigned int i = 0; i < values.size(); i++)
		values[i] = 0.0;
}

template<>
double StaticScalarSolver::Weight<2,0>::value(const Point<2> & r,
	const unsigned int component) const
{
	return 1.0;
}

//-----------------------------------------------------------------------------
// Stage 1. Calculating Jf given T (projection) to check if the Stage 0 was ok.
//-----------------------------------------------------------------------------

template<>
double StaticScalarSolver::Weight<2,1>::value(const Point<2> & r,
	const unsigned int component) const
{
	return 1.0;
}

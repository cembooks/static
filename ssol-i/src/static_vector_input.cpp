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

#include <cmath>
#include <deal.II/base/types.h>
#define BOOST_ALLOW_DEPRECATED_HEADERS

#include <math.h>
#include "static_vector_input.hpp"

using namespace StaticVectorSolver;
using namespace std;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"

template<>
void TheCoefficient<3>::value_list(
	const std::vector<Point<3>> &r,
	types::material_id mid,
	unsigned int cuid,
	std::vector<double> & values) const
{
	Assert(r.size() == values.size(),
		ExcDimensionMismatch(r.size(), values.size()));

	auto v = values.begin();
	for (auto p: r)
	{
		*v = mu_0;
		v++;
	}
}

template<>
void TheCoefficient<2>::value_list(
	const std::vector<Point<2>> &r,
	types::material_id mid,
	unsigned int cuid,
	std::vector<double> & values) const
{
	Assert(r.size() == values.size(),
		ExcDimensionMismatch(r.size(), values.size()));

	auto v = values.begin();
	for (auto p: r)
	{
		*v = 0.0;
		v++;
	}
}

template<>
void PdeRhs<3>::value_list(
	const std::vector<Point<3>> &r,
	types::material_id mid,
	unsigned int cuid,
	std::vector<Tensor<1,3>> & values) const
{
	Assert(r.size() == values.size(),
		ExcDimensionMismatch(r.size(), values.size()));

	for (unsigned int i = 0 ; i < r.size(); i++)
	{
		values[i][0] = 0.0;
		values[i][1] = 0.0;
		values[i][2] = 0.0;
	}
}

template<>
void PdeRhs<2>::value_list(
	const std::vector<Point<2>> &r,
	types::material_id mid,
	unsigned int cuid,
	std::vector<Tensor<1,2>> & values) const
{
	Assert(r.size() == values.size(),
		ExcDimensionMismatch(r.size(), values.size()));

	for (unsigned int i = 0 ; i < r.size(); i++)
	{
		values[i][0] = 0.0;
		values[i][1] = 0.0;
	}
}

template<>
void Gamma<3>::value_list(
	const std::vector<Point<3>> &r,
	const std::vector<Tensor<1, 3>> & n,
	types::boundary_id bid,
	types::material_id mid,
	unsigned int cuid,
	unsigned int fuid,
	std::vector<double> & values) const
{
	Assert(r.size() == values.size(),
		ExcDimensionMismatch(r.size(), values.size()));

	for (unsigned int i = 0 ; i < r.size(); i++)
		values[i] = 1.0 / (mu_0 * r.at(i).norm());
}

template<>
void Gamma<2>::value_list(
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

	for (unsigned int i = 0 ; i < r.size(); i++)
		values[i] = 0.0;
}

template<>
void RobinRhs<3>::value_list(
	const std::vector<Point<3>> &r,
	const std::vector<Tensor<1, 3>> & n,
	types::boundary_id  bid,
	types::material_id mid,
	unsigned int cuid,
	unsigned int fuid,
	std::vector<Tensor<1,3>> & values) const
{
	Assert(r.size() == values.size(),
		ExcDimensionMismatch(r.size(), values.size()));

	for (unsigned int i = 0 ; i < r.size(); i++)
	{
		values[i][0] = 0.0;
		values[i][1] = 0.0;
		values[i][2] = 0.0;
	}
}

template<>
void RobinRhs<2>::value_list(
	const std::vector<Point<2>> &r,
	const std::vector<Tensor<1, 2>> & n,
	types::boundary_id bid,
	types::material_id mid,
	unsigned int cuid,
	unsigned int fuid,
	std::vector<Tensor<1,2>> & values) const
{
	Assert(r.size() == values.size(),
		ExcDimensionMismatch(r.size(), values.size()));

	for (unsigned int i = 0 ; i < r.size(); i++)
	{
		values[i][0] = 0.0;
		values[i][1] = 0.0;
	}
}

template<>
void FreeSurfaceCurrent<3>::value_list(
	const std::vector<Point<3>> &r,
	const std::vector<Tensor<1,3>> &n,
	types::material_id mid,
	unsigned int cuid,
	unsigned int fuid,
	std::vector<Tensor<1,3>> & values) const
{
	Assert(r.size() == values.size(),
		ExcDimensionMismatch(r.size(), values.size()));

	if ( (cuid == 1) && (fuid == 1) )
	{
		for (unsigned int i = 0; i < r.size(); i++)
		{
			values.at(i)[0] =-K_0*r.at(i)[1];
			values.at(i)[1] = K_0*r.at(i)[0];
			values.at(i)[2] = 0.0;
		}
	}
	else
	{
		for (unsigned int i = 0; i < values.size(); i++)
		{
			values.at(i)[0] = 0.0;
			values.at(i)[1] = 0.0;
			values.at(i)[2] = 0.0;
		}
	}
}

template<>
void FreeSurfaceCurrent<2>::value_list(
	const std::vector<Point<2>> &r,
	const std::vector<Tensor<1,2>> & n,
	types::material_id mid,
	unsigned int cuid,
	unsigned int fuid,
	std::vector<Tensor<1,2>> & values) const
{
	Assert(r.size() == values.size(),
		ExcDimensionMismatch(r.size(), values.size()));

	for (unsigned int i = 0 ; i < r.size(); i++)
	{
		values[i][0] = 0.0;
		values[i][1] = 0.0;
	}
}

template<>
double Weight<3>::value(const Point<3> & r,
	const unsigned int component) const
{
	return 1.0;
}

template<>
double Weight<2>::value(const Point<2> & r,
	const unsigned int component) const
{
	return 1.0;
}

#pragma GCC diagnostic pop


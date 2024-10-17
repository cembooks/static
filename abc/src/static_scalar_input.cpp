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

#define BOOST_ALLOW_DEPRECATED_HEADERS

#include <math.h>
#include "static_scalar_input.hpp"

using namespace StaticScalarSolver;
using namespace std;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

template<>
void TheCoefficient<2>::value_list(
	const std::vector<Point<2>> &r,
	types::material_id mid,
	unsigned int cuid,
	std::vector<double> & values) const
{
	Assert(r.size() == values.size(),
		ExcDimensionMismatch(r.size(), values.size()));

	for (unsigned int i = 0; i < values.size(); i++)
		values[i] = ep_0;
}

template<>
void TheCoefficient<3>::value_list(
	const std::vector<Point<3>> &r,
	types::material_id mid,
	unsigned int cuid,
	std::vector<double> & values) const
{
	Assert(r.size() == values.size(),
		ExcDimensionMismatch(r.size(), values.size()));

	for (unsigned int i = 0; i < values.size(); i++)
		values[i] = ep_0;
}

template<>
void PdeRhs<2>::value_list(
	const std::vector<Point<2>> &r,
	types::material_id mid,
	unsigned int cuid,
	std::vector<double> & values) const
{
	Assert(r.size() == values.size(),
		ExcDimensionMismatch(r.size(), values.size()));

	for (unsigned int i = 0; i < values.size(); i++)
		values[i] = 0.0;
}

template<>
void PdeRhs<3>::value_list(
	const std::vector<Point<3>> &r,
	types::material_id mid,
	unsigned int cuid,
	std::vector<double> & values) const
{
	Assert(r.size() == values.size(),
		ExcDimensionMismatch(r.size(), values.size()));

	for (unsigned int i = 0; i < values.size(); i++)
		values[i] = 0.0;
}

template<>
void PdeRhsCvp<2>::value_list(
	const std::vector<Point<2>> &r,
	types::material_id mid,
	unsigned int cuid,
	std::vector<Tensor<1, 2>> & values) const
{
	Assert(r.size() == values.size(),
		ExcDimensionMismatch(r.size(), values.size()));

	for (unsigned int i = 0; i < values.size(); i++)
	{
		values.at(i)[0] = 0.0;
		values.at(i)[1] = 0.0;
	}
}

template<>
void PdeRhsCvp<3>::value_list(
	const std::vector<Point<3>> &r,
	types::material_id mid,
	unsigned int cuid,
	std::vector<Tensor<1, 3>> & values) const
{
	Assert(r.size() == values.size(),
		ExcDimensionMismatch(r.size(), values.size()));

	for (unsigned int i = 0; i < values.size(); i++)
	{
		values.at(i)[0] = 0.0;
		values.at(i)[1] = 0.0;
		values.at(i)[2] = 0.0;
	}
}

template<>
void Gamma<2>::value_list(
	const std::vector<Point<2>> &r,
	const std::vector<Tensor<1, 2>> & n,
	types::boundary_id bin,
	types::material_id mid,
	unsigned int cuid,
	unsigned int fuid,
	std::vector<double> & values) const
{
	Assert(r.size() == values.size(),
		ExcDimensionMismatch(r.size(), values.size()));

	if (bid_infty == 2)
	{
		for (unsigned int i = 0; i < values.size(); i++)
			values[i] = ep_0 / r[i].norm();
	}
	else
	{
		for (unsigned int i = 0; i < values.size(); i++)
			values[i] = 0.0;
	}
}

template<>
void  Gamma<3>::value_list(
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

	if (bid_infty == 2)
	{
		for (unsigned int i = 0; i < values.size(); i++)
			values[i] = ep_0 / r[i].norm();
	}
	else
	{
		for (unsigned int i = 0; i < values.size(); i++)
			values[i] = 0.0;
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
	std::vector<double> & values) const
{

	Assert(r.size() == values.size(),
		ExcDimensionMismatch(r.size(), values.size()));

	for (unsigned int i = 0; i < values.size(); i++)
		values[i] = 0.0;
}

template<>
void RobinRhs<3>::value_list(
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

	for (unsigned int i = 0; i < values.size(); i++)
		values[i] = 0.0;
}

template<>
void FreeSurfaceCharge<2>::value_list(
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
void FreeSurfaceCharge<3>::value_list(
	const std::vector<Point<3>> &r,
	const std::vector<Tensor<1,3>> & n,
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
double Weight<2>::value(const Point<2> & r, const unsigned int) const
{
	if ( (r[0] < 2*x0) && (r[0] > -2*x0) &&
		(r[1] < 2*x0) && (r[1] > -2*x0)
		)
		return 1.0;
	else
		return 0.0;
}

template<>
double Weight<3>::value(const Point<3> & r, const unsigned int) const
{
	if ( r.norm() <= R_mid )
		return 1.0;
	else
		return 0.0;
}

#pragma GCC diagnostic pop


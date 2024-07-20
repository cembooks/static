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
		ExcDimensionMismatch(r.size(), values.size()))

	auto v = values.begin();
	for (auto p: r)
	{
		*v = ep_0 * (pow(p[0],2)*pow(p[1],2) + 1);
		v++;
	}
}

template<>
void TheCoefficient<3>::value_list(
	const std::vector<Point<3>> &r,
	types::material_id mid,
	unsigned int cuid,
	std::vector<double> & values) const
{
	Assert(r.size() == values.size(),
		ExcDimensionMismatch(r.size(), values.size()))

	auto v = values.begin();
	for (auto p: r)
	{
		*v = ep_0 * (pow(p[0],2)*pow(p[1],2)*pow(p[2],2) + 1);
		v++;
	}
}

template<>
void PdeRhs<2>::value_list(
	const std::vector<Point<2>> &r,
	types::material_id mid,
	unsigned int cuid,
	std::vector<double> & values) const
{
	Assert(r.size() == values.size(),
		ExcDimensionMismatch(r.size(), values.size()))

	auto v = values.begin();
	for (auto p: r)
	{
		*v = ep_0 * k * (
			2*p[0]*pow(p[1],2)*sin(k*p[0]) +
			2*p[1]*pow(p[0],2)*sin(k*p[1]) +
			k*(pow(p[0],2)*pow(p[1],2) + 1) *
			(cos(k*p[0]) + cos(k*p[1]))
			);
		v++;
	}
}

template<>
void PdeRhs<3>::value_list(
	const std::vector<Point<3>> &r,
	types::material_id mid,
	unsigned int cuid,
	std::vector<double> & values) const
{
	Assert(r.size() == values.size(),
		ExcDimensionMismatch(r.size(), values.size()))

	auto v = values.begin();
	for (auto p: r)
	{
		*v = ep_0 * k * (
			2*p[0]*pow(p[1],2)*pow(p[2],2)*sin(k*p[0]) +
			2*p[1]*pow(p[0],2)*pow(p[2],2)*sin(k*p[1]) +
			2*p[2]*pow(p[1],2)*pow(p[0],2)*sin(k*p[2]) +
			k*(pow(p[0],2)*pow(p[1],2)*pow(p[2],2) + 1) *
			(cos(k*p[0])+cos(k*p[1])+cos(k*p[2]))
			);
		v++;
	}
}

template<>
void PdeRhsCvp<2>::value_list(
	const std::vector<Point<2>> &r,
	types::material_id mid,
	unsigned int cuid,
	std::vector<Tensor<1, 2>> & values) const
{
	Assert(r.size() == values.size(),
		ExcDimensionMismatch(r.size(), values.size()))

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
		ExcDimensionMismatch(r.size(), values.size()))

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
	types::boundary_id bid,
	types::material_id mid,
	unsigned int cuid,
	unsigned int fuid,
	std::vector<double> & values) const
{
	Assert(r.size() == values.size(),
		ExcDimensionMismatch(r.size(), values.size()))

	auto v = values.begin();
	for (auto p: r)
	{
		*v = ( ep_0 * (pow(p[0],2)*pow(p[1],2) + 1) ) *
			(p.norm() + 2);
		v++;
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
		ExcDimensionMismatch(r.size(), values.size()))

	auto v = values.begin();
	for (auto p: r)
	{
		*v = ( ep_0 * (pow(p[0],2)*pow(p[1],2)*pow(p[2],2) + 1) ) *
			(p.norm() + 2);
		v++;
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
		ExcDimensionMismatch(r.size(), values.size()))

	double epsilon;
	double gamma;
	double phi;
	Tensor<1,2> grad_phi;

	auto v = values.begin();
	auto nn = n.begin();
	for (auto p: r)
	{
		epsilon = ep_0 * (pow(p[0],2)*pow(p[1],2) + 1);
		phi = cos(k*p[0]) + cos(k*p[1]);
		grad_phi[0] = -k*sin(k*p[0]);
		grad_phi[1] = -k*sin(k*p[1]);
		gamma = epsilon * (p.norm() + 2);

		*v = epsilon * (*nn * grad_phi) + gamma * phi;
		 v++;
		 nn++;
	}
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
		ExcDimensionMismatch(r.size(), values.size()))

	double epsilon;
	double gamma;
	double phi;
	Tensor<1,3> grad_phi;

	auto v = values.begin();
	auto nn = n.begin();
	for (auto p: r)
	{
		epsilon = ep_0 * (pow(p[0],2)*pow(p[1],2)*pow(p[2],2) + 1);
		phi = (cos(k*p[0]) + cos(k*p[1]) + cos(k*p[2]));
		grad_phi[0] = -k*sin(k*p[0]);
		grad_phi[1] = -k*sin(k*p[1]);
		grad_phi[2] = -k*sin(k*p[2]);
		gamma = epsilon * (p.norm() + 2);
		*v = epsilon * (*nn * grad_phi) + gamma * phi;
		v++;
		nn++;
		}
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
		ExcDimensionMismatch(r.size(), values.size()))

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
		ExcDimensionMismatch(r.size(), values.size()))

	for (unsigned int i = 0; i < values.size(); i++)
		values[i] = 0.0;
}

template<>
double Weight<2>::value(const Point<2> & r,
	const unsigned int component) const
{
	return 1.0;
}

template<>
double Weight<3>::value(const Point<3> & r,
	const unsigned int component) const
{
	return 1.0;
}

#pragma GCC diagnostic pop


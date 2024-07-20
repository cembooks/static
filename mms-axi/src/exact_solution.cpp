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

#include "exact_solution.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

using namespace dealii;
using namespace std;

template<>
double ExactSolutionMMSAXI_PHI<2>::value(const Point<2> & p,
	unsigned int component) const
{
	return (cos(k*p[0]) + cos(k*p[1]));
}

template<>
Tensor<1,2> ExactSolutionMMSAXI_PHI<2>::gradient(const Point<2> & p,
	unsigned int component) const
{
	Tensor<1,2> g;

	g[0] = - k * sin( k*p[0] );
	g[1] = - k * sin( k*p[1] );

	return g;
}

template<>
double ExactSolutionMMSAXI_PHI<3>::value(const Point<3> & p,
	unsigned int componet) const
{
	double r = sqrt( pow(p[0],2) + pow(p[1],2) );

	return (
		cos( k*r )
	+ cos( k*p[2] )
			   );
}

template<>
Tensor<1,3> ExactSolutionMMSAXI_PHI<3>::gradient(const Point<3> & p,
	unsigned int component) const
{
	Tensor<1,3> g;

	double r = sqrt( pow(p[0],2) + pow(p[1],2) );

	if ( r < eps)
	{
		g[0] = 0.0;
		g[1] = 0.0;
	}
	else
	{
		double cos_phi = p[0] / r;
		double sin_phi = p[1] / r;
		g[0] = - k * cos_phi * sin( k*r );
		g[1] = - k * sin_phi * sin( k*r );
	}

	g[2] = - k * sin( k*p[2] );

	return g;
}

template<>
void ExactSolutionMMSAXI_E<2>::vector_value_list(const std::vector<Point<2>> & p,
		std::vector<Vector<double>>	 &values) const
{
	Assert(values.size() == p.size(), ExcDimensionMismatch(values.size(), p.size()));

	auto v = values.begin();
	for (auto q: p)
	{
		(*v)(0) = k*sin(k * q[0]);
		(*v)(1) = k*sin(k * q[1]);
		v++;
	}
}

template<>
void ExactSolutionMMSAXI_E<3>::vector_value_list(const std::vector<Point<3>> & p,
		std::vector<Vector<double>>	 &values) const
{
	Assert(values.size() == p.size(), ExcDimensionMismatch(values.size(), p.size()));

	double r;
	double cos_phi;
	double sin_phi;

	auto v = values.begin();

	for (auto q: p)
	{
		r = sqrt( pow(q[0],2) + pow(q[1],2) );

		if ( r < eps)
		{
			(*v)(0) = 0.0;
			(*v)(1) = 0.0;
		}
		else
		{
			cos_phi = q[0] / r;
			sin_phi = q[1] / r;
			(*v)(0) = k * cos_phi * sin( k*r );
			(*v)(1) = k * sin_phi * sin( k*r );
		}

		(*v)(2)	= k * sin( k * q[2] );

		v++;
	}
}

template<>
void ExactSolutionMMSAXI_D<2>::vector_value_list(const std::vector<Point<2>> & p,
		std::vector<Vector<double>>	 &values) const
{
	Assert(values.size() == p.size(), ExcDimensionMismatch(values.size(), p.size()));

	double epsilon;

	auto v = values.begin();
	for (auto q: p)
	{
		epsilon = ep_0 * ( q[0] * pow(q[1],2) + 1.0 );
		(*v)(0) = epsilon*k*sin(k * q[0]);
		(*v)(1) = epsilon*k*sin(k * q[1]);
		v++;
	}
}

template<>
void ExactSolutionMMSAXI_D<3>::vector_value_list(const std::vector<Point<3>> & p,
		std::vector<Vector<double>>	 &values) const
{
	Assert(values.size() == p.size(), ExcDimensionMismatch(values.size(), p.size()));

	double r;
	double epsilon;
	double cos_phi;
	double sin_phi;

	auto v = values.begin();
	for (auto q: p)
	{
		r = sqrt( pow(q[0],2) + pow(q[1],2) );
		epsilon = ep_0 * ( r * pow(q[2],2) + 1.0 );

		if ( r < eps)
		{
			(*v)(0) = 0.0;
			(*v)(1) = 0.0;
		}
		else
		{
			cos_phi = q[0] / r;
			sin_phi = q[1] / r;
			(*v)(0) = epsilon * k * cos_phi * sin( k*r );
			(*v)(1) = epsilon * k * sin_phi * sin( k*r );
		}

		(*v)(2)	= epsilon * k * sin( k * q[2] );

		v++;
	}
}
#pragma GCC diagnostic pop


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
#include <cmath>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

using namespace dealii;
using namespace std;

template<>
double ExactSolutionPLS_PHI<2>::value(const Point<2> &r,
	unsigned int component) const
{
	return (cos(k*r[0]) + cos(k*r[1]));
}

template<>
Tensor<1,2> ExactSolutionPLS_PHI<2>::gradient(const Point<2> &r,
	unsigned int component) const
{
	Tensor<1,2> p;

	p[0] = -k*sin(k*r[0]);
	p[1] = -k*sin(k*r[1]);

	return p;
}

template<>
double ExactSolutionPLS_PHI<3>::value(const Point<3> &r,
	unsigned int componet) const
{
	return (cos(k*r[0]) + cos(k*r[1]) + cos(k*r[2]));
}

template<>
Tensor<1,3> ExactSolutionPLS_PHI<3>::gradient(const Point<3> &r,
	unsigned int component) const
{
	Tensor<1,3> p;

	p[0]=-k*sin(k*r[0]);
	p[1]=-k*sin(k*r[1]);
	p[2]=-k*sin(k*r[2]);

	return p;
}

#pragma GCC diagnostic pop


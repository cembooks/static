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
double ExactSolutionSCHAXI_PHI<true>::value(const Point<2> & p,
	const unsigned int component) const
{
	Point<2> r;
	r[0] = p[0];
	r[1] = 0.0;

	if ( r.norm() > a )
	{
		return (log(b) - log(sqrt(r.square()))) /
			(log(b) - log(a));
	}

	return 1.0;
}

template<>
Tensor<1, 2> ExactSolutionSCHAXI_PHI<true>::gradient(const Point<2> & p,
	const unsigned int component) const
{
	Point<2> r;
	r[0] = p[0];
	r[1] = 0.0;

	if ( r.norm() > a )
	{
		return - 1.0 / (log(b)-log(a)) * r / r.square();
	}

	return Point<2>();
}

template<>
double ExactSolutionSCHAXI_PHI<false>::value(const Point<2> & p,
	const unsigned int component) const
{
	Point<2> r;
	r[0] = p[0];
	r[1] = p[1];

	if ( r.norm() > a )
	{
		return (a * b/(b - a))*(1 / sqrt(r.square())-1 / b);
	}

	return 1.0;
}

template<>
Tensor<1, 2> ExactSolutionSCHAXI_PHI<false>::gradient(const Point<2> & p,
	const unsigned int component) const
{
	Point<2> r;
	r[0] = p[0];
	r[1] = p[1];

	if ( r.norm() > a )
	{
		return - a * b / ( b - a ) * r / pow(sqrt(r.square()), 3);
	}

	return Point<2>();
}

#pragma GCC diagnostic pop


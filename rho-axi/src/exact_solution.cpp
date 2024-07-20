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
double ExactSolutionRHOAXI_PHI<true>::value(const Point<2> & p,
	const unsigned int component) const
{
	Point<2> r;
	r[0] = p[0];
	r[1] = 0.0;

	if ( r.norm() <= a )
	{
		return rho*pow(a,2)*(1.0+2.0*log(b/a)-pow(r.norm()/a,2))/(4.0*ep_0);
	}else
	{
		return rho*pow(a,2)*log(b/r.norm())/(2.0*ep_0);
	}

	return 0.0;
}

template<>
Tensor<1, 2> ExactSolutionRHOAXI_PHI<true>::gradient(const Point<2> & p,
	const unsigned int component) const
{
	Point<2> r;
	r[0] = p[0];
	r[1] = 0.0;

	if ( r.norm() <= a )
	{
		return -rho*r/(2*ep_0);
	}else
	{
		return -rho*pow(a,2)*r/(2.0*ep_0*r.norm_square());
	}

	return Point<2>();
}

template<>
double ExactSolutionRHOAXI_PHI<false>::value(const Point<2> & p,
	const unsigned int component) const
{
	Point<2> r;
	r[0] = p[0];
	r[1] = p[1];

	if ( r.norm() <= a )
	{
		return rho*(pow(a,2) + 2.0*pow(a,3)*(1.0/a-1.0/b) - r.norm_square())/(6.0*ep_0);
	}else
	{
		return rho*pow(a,3)*(1.0/r.norm() - 1.0/b)/(3.0*ep_0);
	}

	return 0.0;
}

template<>
Tensor<1, 2> ExactSolutionRHOAXI_PHI<false>::gradient(const Point<2> & p,
	const unsigned int component) const
{
	Point<2> r;
	r[0] = p[0];
	r[1] = p[1];

	if ( r.norm() < a )
	{
		return -rho*r/(3.0*ep_0);
	}else
	{
		return -rho*pow(a,3)*r/(3.0*ep_0*pow(r.norm(),3));
	}

	return Point<2>();
}

#pragma GCC diagnostic pop


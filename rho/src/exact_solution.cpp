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
double ExactSolutionRHO_PHI<2>::value(const Point<2> & r,
	const unsigned int component) const
{

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
Tensor<1, 2> ExactSolutionRHO_PHI<2>::gradient(const Point<2> & r,
	const unsigned int component) const
{

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
double ExactSolutionRHO_PHI<3>::value(const Point<3> & r,
	const unsigned int component) const
{
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
Tensor<1, 3> ExactSolutionRHO_PHI<3>::gradient(const Point<3> & r,
	const unsigned int component) const
{
	if ( r.norm() < a )
	{
		return -rho*r/(3.0*ep_0);
	}else
	{
		return -rho*pow(a,3)*r/(3.0*ep_0*pow(r.norm(),3));
	}

	return Point<3>();
}

#pragma GCC diagnostic pop


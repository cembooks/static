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
ExactSolutionFLC_PHI<2>::ExactSolutionFLC_PHI()
{
	alpha = 1.0 / ( ep_1*log(b/d_2) + ep_2*log(d_1/a) );
	beta = (ep_1/ep_2)*log(b/d_2);
	phi_d = alpha*ep_1*log(b/d_2);
}

template<>
double ExactSolutionFLC_PHI<2>::value(const Point<2> & r,
	const unsigned int component) const
{
	if ( r.norm() < d_1 )
		return (- alpha * ep_2 * ( log( r.norm()/d_1 ) - beta ));

	if ( r.norm() > d_2 )
		return (- alpha * ep_1 * (log( r.norm()/b )) );

	return phi_d;
}

template<>
Tensor<1, 2> ExactSolutionFLC_PHI<2>::gradient(const Point<2> & r,
	const unsigned int component) const
{
	if ( r.norm() < d_1 )
		return - alpha * ep_2 * (1 / r.square()) * r;

	if ( r.norm() > d_2 )
		return - alpha * ep_1 * (1 / r.square()) * r;

	return Point<2>();
}

template<>
ExactSolutionFLC_PHI<3>::ExactSolutionFLC_PHI()
{
	alpha = 1.0 / ( ep_1*( 1.0/b - 1.0/d_2 ) + ep_2*( 1.0/d_1 - 1.0/a ) );
	beta = 1.0/d_1 + (ep_1/ep_2)*( 1.0/b - 1.0/d_2 );
	phi_d = alpha*ep_1*(1.0 / b - 1.0 / d_2);
}

template<>
double ExactSolutionFLC_PHI<3>::value(const Point<3> & r,
	const unsigned int component) const
{
	if ( r.norm() < d_1 )
		return (- alpha * ep_2 * ( 1.0/r.norm() - beta ));

	if ( r.norm() > d_2 )
		return (- alpha * ep_1 * ( 1.0/r.norm() - 1.0/b ));

	return phi_d;
}

template<>
Tensor<1, 3> ExactSolutionFLC_PHI<3>::gradient(const Point<3> & r,
	const unsigned int component) const
{
	if ( r.norm() < d_1 )
		return alpha * ep_2 * r / pow(r.norm(),3);

	if ( r.norm() > d_2 )
		return alpha * ep_1 * r / pow(r.norm(),3);

	return Point<3>();
}

#pragma GCC diagnostic pop


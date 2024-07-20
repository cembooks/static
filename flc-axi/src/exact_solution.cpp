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
ExactSolutionFLCAXI_PHI<true>::ExactSolutionFLCAXI_PHI()
{
	alpha = 1.0 / ( ep_1*log(b/d_2) + ep_2*log(d_1/a) );
	beta = (ep_1/ep_2)*log(b/d_2);
	phi_d = alpha*ep_1*log(b/d_2);
}

template<>
double ExactSolutionFLCAXI_PHI<true>::value(const Point<2> & p,
	const unsigned int component) const
{
	Point<2> r;
	r[0] = p[0];
	r[1] = 0.0;

	if ( r.norm() < d_1 )
		return (- alpha * ep_2 * ( log( r.norm()/d_1 ) - beta ));

	if ( r.norm() > d_2 )
		return (- alpha * ep_1 * (log( r.norm()/b )) );

	return phi_d;
}

template<>
Tensor<1, 2> ExactSolutionFLCAXI_PHI<true>::gradient(const Point<2> & p,
	const unsigned int component) const
{
	Point<2> r;
	r[0] = p[0];
	r[1] = 0.0;

	if ( r.norm() < d_1 )
		return - alpha * ep_2 * (1 / r.square()) * r;

	if ( r.norm() > d_2 )
		return - alpha * ep_1 * (1 / r.square()) * r;

	return Point<2>();
}

template<>
ExactSolutionFLCAXI_PHI<false>::ExactSolutionFLCAXI_PHI()
{
	alpha = 1.0 / ( ep_1*( 1.0/b - 1.0/d_2 ) + ep_2*( 1.0/d_1 - 1.0/a ) );
	beta = 1.0/d_1 + (ep_1/ep_2)*( 1.0/b - 1.0/d_2 );
	phi_d = alpha*ep_1*(1.0 / b - 1.0 / d_2);
}

template<>
double ExactSolutionFLCAXI_PHI<false>::value(const Point<2> & p,
	const unsigned int component) const
{
	Point<2> r;
	r[0] = p[0];
	r[1] = p[1];

	if ( r.norm() < d_1 )
		return (- alpha * ep_2 * ( 1.0/r.norm() - beta ));

	if ( r.norm() > d_2 )
		return (- alpha * ep_1 * ( 1.0/r.norm() - 1.0/b ));

	return phi_d;
}

template<>
Tensor<1, 2> ExactSolutionFLCAXI_PHI<false>::gradient(const Point<2> & p,
	const unsigned int component) const
{
	Point<2> r;
	r[0] = p[0];
	r[1] = p[1];

	if ( r.norm() < d_1 )
		return alpha * ep_2 * r / pow(r.norm(),3);

	if ( r.norm() > d_2 )
		return alpha * ep_1 * r / pow(r.norm(),3);

	return Point<2>();
}

#pragma GCC diagnostic pop


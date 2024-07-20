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
double ExactSolutionCBND_PHI<2>::value(const Point<2> & r,
	const unsigned int component) const
{
	return (log(b) - log(r.norm())) /
		(log(b) - log(a));
}

template<>
Tensor<1, 2> ExactSolutionCBND_PHI<2>::gradient(const Point<2> & r,
	const unsigned int component) const
{
	return - 1.0 / (log(b)-log(a)) * r / r.square();
}

template<>
double ExactSolutionCBND_PHI<3>::value(const Point<3> & r,
	const unsigned int component) const
{
	return (a * b/(b - a))*(1 / r.norm()-1 / b);
}

template<>
Tensor<1, 3> ExactSolutionCBND_PHI<3>::gradient(const Point<3> & r,
	const unsigned int component) const
{
	return - a * b / ( b - a ) * r / pow(r.norm(), 3);
}

#pragma GCC diagnostic pop


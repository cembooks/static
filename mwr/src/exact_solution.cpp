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

double ExactSolutionMWR_A::value(const Point<2> & r,
	const unsigned int component) const
{
	if ( r.norm() < a )
	{
		return -0.25*mu*J_f*(r[0]*r[0] + r[1]*r[1]) + A_0;
	}else
	{
		return -0.25*mu_0*J_f*a*a*log(exp(mu_r)*(r[0]*r[0] + r[1]*r[1])/(a*a)) + A_0;
	}
}

Tensor<1, 2> ExactSolutionMWR_A::gradient(const Point<2> & r,
	const unsigned int component) const
{
	Point<2> grad_A;

	if ( r.norm() <= a )
	{
		grad_A(0) = -0.5*mu*J_f*r[0];
		grad_A(1) = -0.5*mu*J_f*r[1];
	}else
	{
		grad_A(0) = -0.5*mu_0*J_f*a*a*r[0]/(r[0]*r[0] + r[1]*r[1]);
		grad_A(1) = -0.5*mu_0*J_f*a*a*r[1]/(r[0]*r[0] + r[1]*r[1]);
	}

	return grad_A;
}

void ExactSolutionMWR_H::vector_value_list(const std::vector<Point<2>> & r,
	std::vector<Vector<double>> &values) const
{
	Assert(values.size() == r.size(), ExcDimensionMismatch(values.size(), r.size()));

	Tensor<1, 2> grad_A;
 	double coef;

	for (unsigned int i = 0; i < r.size(); i++)
	{
		if ( r.at(i).norm() < a )
		{
			coef = mu;
		}else
		{
			coef = mu_0;
		}

		grad_A = A.gradient(r.at(i));
		values.at(i)[0] =  grad_A[1] / coef;
		values.at(i)[1] = -grad_A[0] / coef;
	}
}

void ExactSolutionMWR_B::vector_value_list(const std::vector<Point<2>> & r,
	std::vector<Vector<double>> &values) const
{
	Assert(values.size() == r.size(), ExcDimensionMismatch(values.size(), r.size()));

	Tensor<1, 2> grad_A;

	for (unsigned int i = 0; i < r.size(); i++)
	{
		grad_A = A.gradient(r.at(i));
		values.at(i)[0] =  grad_A[1];
		values.at(i)[1] = -grad_A[0];
	}
}

#pragma GCC diagnostic pop


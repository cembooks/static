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

double ExactSolutionMMSVTII_T::value(const Point<2> & r,
	const unsigned int component) const
{
		return current_vector_potential(r[0], r[1], mu_0, k);
}

Tensor<1, 2> ExactSolutionMMSVTII_T::gradient(const Point<2> & r,
	const unsigned int component) const
{
	double Bm = magnetic_field(r[0], r[1], k);
	double mu = permeability(r[0], r[1], mu_0);

	Point<2> grad_Tm;
	grad_Tm[0] = - 2.0*r[0]*mu_0/pow(mu,2)*Bm - (k / mu)* sin(k*r[0]);
	grad_Tm[1] = - 2.0*r[1]*mu_0/pow(mu,2)*Bm - (k / mu)* sin(k*r[1]);

	return grad_Tm;
}

void ExactSolutionMMSVTII_Jf::vector_value_list(const std::vector<Point<2>> & r,
		std::vector<Vector<double>>	 &values) const
{
	Assert(values.size() == r.size(), ExcDimensionMismatch(values.size(), r.size()));

	Tensor<1,2> Jf;

	for (unsigned int i = 0; i < values.size(); i++)
	{
		Jf = volume_free_current_density(r.at(i)[0], r.at(i)[1], mu_0, k);

		values.at(i)[0] = Jf[0];
		values.at(i)[1] = Jf[1];
	}
}

void ExactSolutionMMSVTII_B::value_list(const std::vector<Point<2>> & r,
	std::vector<double>	 &values,
	const unsigned int component) const
{
	Assert(values.size() == r.size(), ExcDimensionMismatch(values.size(), r.size()));

	for (unsigned int i = 0; i < values.size(); i++)
		values.at(i) = magnetic_field(r.at(i)[0], r.at(i)[1], k);
}

void DirichletBC_T::value_list(const std::vector<Point<2>> & r,
	std::vector<double>	 &values,
	const unsigned int component) const
{
	Assert(values.size() == r.size(), ExcDimensionMismatch(values.size(), r.size()));

	for (unsigned int i = 0; i < values.size(); i++)
		values.at(i)= current_vector_potential(r.at(i)[0], r.at(i)[1], mu_0, k);
}

void DirichletBC_A::vector_value_list(const std::vector<Point<2>> & r,
	std::vector<Vector<double>>	 &values) const
{
	Assert(values.size() == r.size(), ExcDimensionMismatch(values.size(), r.size()));

	Tensor<1,2> A ;

	for (unsigned int i = 0; i < values.size(); i++)
	{
		A = vector_potential(r.at(i)[0], r.at(i)[1], k);

		values.at(i)[0] = A[0];
		values.at(i)[1] = A[1];
	}
}


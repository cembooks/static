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

#ifndef ExactSolutionSSOLIIAXI_H__
#define ExactSolutionSSOLIIAXI_H__

#include <deal.II/base/function.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/vector.h>

#include <cmath>

#include "constants.hpp"
#include "settings.hpp"

using namespace dealii;

inline Tensor<1, 3> magnetic_field(double x, double y, double z,
		double K_0, double mu_0,  double a, double b)
{
	double cos_theta;
	double sin_theta;

	double cos_phi;
	double sin_phi;

	Tensor<1,3> r_hat;
	Tensor<1,3> theta_hat;

	Tensor<1,3> F1;
	Tensor<1,3> F2;
	Tensor<1,3> B;

	double r;

	r = sqrt(pow(x,2) + pow(y,2) + pow(z,2));

	cos_theta = z / r;
	sin_theta = sqrt(pow(x,2)+pow(y,2)) / r;

	cos_phi = x / sqrt(pow(x,2)+pow(y,2));
	sin_phi = y / sqrt(pow(x,2)+pow(y,2));

	r_hat[0] = x / r;
	r_hat[1] = y / r;
	r_hat[2] = z / r;

	theta_hat[0] = cos_theta*cos_phi;
	theta_hat[1] = cos_theta*sin_phi;
	theta_hat[2] = -sin_theta;

	F1 = (2.0 / 3.0) * mu_0 * K_0 * (cos_theta*r_hat - sin_theta*theta_hat);
	F2 = (2.0 / 3.0) * mu_0 * K_0 * (cos_theta*r_hat + 0.5*sin_theta*theta_hat) / pow(r,3);

	if (r < a)
	{
		B = 0.5*(pow(b,2) - pow(a,2))*F1;
	}
	else if (r > b)
	{
		B = 0.2*(pow(b,5) - pow(a,5))*F2;
	}
	else
	{
		B = 0.5*(pow(b,2) - pow(r,2))*F1 + 0.2*(pow(r,5) - pow(a,5))*F2;
	}

	return B;
}

/**
 * \brief Describes exact solution, \f$ A' \f$, of the
 * [Axisymmetric - thick spherical coil (ssol-ii-axi)](@ref page_ssol_ii_axi)
 * numerical experiment.
 *****************************************************************************/
class ExactSolutionSSOLIIAXI_A : public Function<2>, public SettingsSSOLIIAXI
{
public:

	ExactSolutionSSOLIIAXI_A(){};

	virtual double value(const Point<2> & p,
		const unsigned int component = 0) const override final;

	virtual Tensor<1, 2> gradient(const Point<2> & p,
				const unsigned int component = 0) const override final;
};

/**
 * \brief Describes exact solutions, \f$\vec{B}'\f$, of the
 * [Axisymmetric - thick spherical coil (ssol-ii-axi)](@ref page_ssol_ii_axi)
 * numerical experiment.
 *****************************************************************************/
class ExactSolutionSSOLIIAXI_B : public Function<2>, public SettingsSSOLIIAXI
{
public:

	ExactSolutionSSOLIIAXI_B(): Function<2>(2), B_0(2.0*mu_0*K_0/3.0) {}

	virtual void vector_value_list(const std::vector<Point<2>> & r,
		std::vector<Vector<double>> &values) const final;

private:
	const double B_0;
};

/**
 * \brief Describes exact solution, \f$\vec{H}'\f$, of the
 * [Axisymmetric - thick spherical coil (ssol-ii-axi)](@ref page_ssol_ii_axi)
 * numerical experiment.
 *****************************************************************************/
class ExactSolutionSSOLIIAXI_H: public Function<2>, public SettingsSSOLIIAXI
{
public:

	ExactSolutionSSOLIIAXI_H(): Function<2>(2) {}

	virtual void vector_value_list(const std::vector<Point<2>> & r,
		std::vector<Vector<double>> &values) const final;

private:

	ExactSolutionSSOLIIAXI_B B;
};

#endif


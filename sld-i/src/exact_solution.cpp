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
ExactSolutionSLDI_PSI<2>::ExactSolutionSLDI_PSI() : Function<2>(), SettingsSLDI()
{
	double mur = mur_2;
	double a2 = std::pow(a,2.0);
	double b2 = std::pow(b,2.0);

	OMEGA = ((mur-1.0)/(mur+1.0))*(a2/b2);
	gamma_1 = (-2.0*b2*H_0*OMEGA)/((mur+1.0)-(mur-1.0)*OMEGA);
	beta_1 = (mur+1.0)*gamma_1/((mur-1.0)*a2);
	alpha_1 = -b2*H_0+mur*gamma_1-mur*b2*beta_1;
	delta_1 = (mur*a2*beta_1-mur*gamma_1)/a2;
}

template<>
double ExactSolutionSLDI_PSI<2>::value(const Point<2> &r,
	unsigned int component) const
{
	double rn = r.norm();
	double r2 = std::pow(rn,2.0);

	if ( rn < a )
		return delta_1*r[0];

	if ( rn > b )
		return (-H_0+alpha_1/r2)*r[0];

	return (beta_1+gamma_1/r2)*r[0];
}

template<>
Tensor<1,2> ExactSolutionSLDI_PSI<2>::gradient(const Point<2> &r,
	unsigned int component) const
{
	double rn = r.norm();
	double r2 = std::pow(rn,2.0);
	double xx = -2.0*r[0]*r[0]/std::pow(rn,4.0);
	double xy = -2.0*r[0]*r[1]/std::pow(rn,4.0);

	Point<2> grad_phi;

	if ( rn < a )
	{
		grad_phi[0] = delta_1;
		grad_phi[1] = 0.0;
		return grad_phi;
	}

	if ( rn > b )
	{
		grad_phi[0] = -H_0+alpha_1/r2+alpha_1*xx;
		grad_phi[1] = alpha_1*xy;
		return grad_phi;
	}

	grad_phi[0] = beta_1+gamma_1/r2+gamma_1*xx;
	grad_phi[1] = gamma_1*xy;
	return grad_phi;
}

template<>
ExactSolutionSLDI_PSI<3>::ExactSolutionSLDI_PSI() : Function<3>(),
	SettingsSLDI()
{
	double mur = mur_2;
	double a3 = std::pow(a,3.0);
	double b3 = std::pow(b,3.0);

	OMEGA = ((mur-1.0)/(mur+2.0))*(a3/b3);
	gamma_1 = (-3.0*b3*H_0*OMEGA)/((2.0*mur+1.0)-2.0*(mur-1.0)*OMEGA);
	beta_1 = ((2.0*mur+1.0)*gamma_1)/((mur-1.0)*a3);
	alpha_1 = (-b3*H_0+2.0*mur*gamma_1-mur*b3*beta_1)/2.0;
	delta_1 = (mur*a3*beta_1-2.0*mur*gamma_1)/a3;
}

template<>
double ExactSolutionSLDI_PSI<3>::value(const Point<3> &r,
	unsigned int componet) const
{
	double rn = r.norm();
	double r3 = std::pow(rn,3.0);

	if ( rn < a )
		return delta_1*r[2];

	if ( rn > b )
		return (-H_0+alpha_1/r3)*r[2];

	return (beta_1+gamma_1/r3)*r[2];
}

template<>
Tensor<1,3> ExactSolutionSLDI_PSI<3>::gradient(const Point<3> &r,
	unsigned int component) const
{
	double rn = r.norm();
	double r3 = std::pow(rn,3.0);
	double zz = -3.0*r[2]*r[2]/std::pow(rn,5.0);
	double xz = -3.0*r[0]*r[2]/std::pow(rn,5.0);
	double yz = -3.0*r[1]*r[2]/std::pow(rn,5.0);

	Point<3> grad_phi;

	if ( rn < a )
	{
		grad_phi[0] = 0.0;
		grad_phi[1] = 0.0;
		grad_phi[2] = delta_1;
		return grad_phi;
	}

	if ( rn > b )
	{
		grad_phi[0] = alpha_1*xz;
		grad_phi[1] = alpha_1*yz;
		grad_phi[2] = -H_0+alpha_1/r3+alpha_1*zz;
		return grad_phi;
	}

	grad_phi[0] = gamma_1*xz;
	grad_phi[1] = gamma_1*yz;
	grad_phi[2] = beta_1+gamma_1/r3+gamma_1*zz;
	return grad_phi;
}

#pragma GCC diagnostic pop


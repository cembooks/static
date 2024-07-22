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

#ifndef ExactSolutionsSSOLIII_H__
#define ExactSolutionsSSOLIII_H__

#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include "constants.hpp"
#include "settings.hpp"

#include <cmath>

using namespace dealii;

inline Tensor<1, 3> volume_free_current_density(double x, double y, double z,
		double K_0, double a, double b)
{
	Tensor<1, 3> J;
	double r;

	r = sqrt(pow(x,2) + pow(y,2) + pow(z,2));

	if ((r>a)&&(r<b))
	{
		J[0] = - K_0 * y;
		J[1] =   K_0 * x;
		J[2] = 0.0;
	}
	else
	{
		J[0] = 0.0;
		J[1] = 0.0;
		J[2] = 0.0;
	}

	return J;
}

inline Tensor<1, 3> magnetic_field_coil(double x, double y, double z,
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

inline Tensor<1, 3> magnetic_field_core(double x, double y, double z,
		double H_0, double mur, double mu0,  double a, double b)
{
	double a3 = std::pow(a,3);
	double b3 = std::pow(b,3);

	double OMEGA = ((mur-1.0)/(mur+2.0))*(a3/b3);
	double gamma_1 = (-3.0*b3*H_0*OMEGA)/((2.0*mur+1.0)-2.0*(mur-1.0)*OMEGA);
	double beta_1 = ((2.0*mur+1.0)*gamma_1)/((mur-1.0)*a3);
	double alpha_1 = (-b3*H_0+2.0*mur*gamma_1-mur*b3*beta_1)/2.0;
	double delta_1 = (mur*a3*beta_1-2.0*mur*gamma_1)/a3;

	double r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
	double r3 = std::pow(r,3.0);
	double zz = -3.0*z*z/std::pow(r,5.0);
	double xz = -3.0*x*z/std::pow(r,5.0);
	double yz = -3.0*y*z/std::pow(r,5.0);

	double mu;

	Tensor<1, 3> grad_psi;
	Tensor<1, 3> B;
	Tensor<1, 3> Happl;


	if (r < a)
	{
		grad_psi[0] = 0.0;
		grad_psi[1] = 0.0;
		grad_psi[2] = delta_1;
		mu = mu0;
	}
	else if ( r > b )
	{
		grad_psi[0] = alpha_1*xz;
		grad_psi[1] = alpha_1*yz;
		grad_psi[2] = -H_0+alpha_1/r3+alpha_1*zz;
		mu = mu0;
	}
	else
	{
		grad_psi[0] = gamma_1*xz;
		grad_psi[1] = gamma_1*yz;
		grad_psi[2] = beta_1+gamma_1/r3+gamma_1*zz;
		mu = mur*mu0;
	}

	Happl[0] = 0.0;
	Happl[1] = 0.0;
	Happl[2] = H_0;

	B = -mu*grad_psi-mu0*Happl;

	return B;
}

/**
 * \brief Describes exact solution, \f$\vec{J}_f\f$, of the
 * [Thick spherical coil with magnetic core (ssol-iii/)](@ref page_ssol_iii)
 * numerical experiment.
 *****************************************************************************/
class ExactSolutionSSOLIII_Jf : public Function<3>, public SettingsSSOLIII
{
public:

	ExactSolutionSSOLIII_Jf();

	virtual void vector_value_list(const std::vector<Point<3>> & r,
		std::vector<Vector<double>> &values) const override final;
};

/**
 * \brief Describes exact solution, \f$\vec{B}\f$, of the
 * [Thick spherical coil with magnetic core (ssol-iii/)](@ref page_ssol_iii)
 * numerical experiment.
 *****************************************************************************/
class ExactSolutionSSOLIII_B : public Function<3>, public SettingsSSOLIII
{
public:

	ExactSolutionSSOLIII_B();

	virtual void vector_value_list(const std::vector<Point<3>> & r,
		std::vector<Vector<double>> &values) const override final;
};

#endif


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

#ifndef ExactSolutionsMMSVTI_H__
#define ExactSolutionsMMSVTI_H__

#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include "constants.hpp"
#include "settings.hpp"

#include <cmath>

using namespace dealii;

inline double permeability(double x, double y, double mu_0)
{
	return mu_0*(pow(x,2)+pow(y,2)+1.0);
}

inline double robin_gamma(double x, double y, double mu_0)
{
	return (sqrt(pow(x,2)+pow(y,2))+2.0)/permeability(x,y,mu_0);
}

inline Tensor<1, 3> volume_free_current_density(double x, double y, double mu_0, double k)
{
	Tensor<1, 3> J;
	const double mu = permeability(x, y, mu_0);

	J[0] = (1/mu)*(
		(mu_0/mu)*(-2*y*(cos(k*x)+cos(k*y)))
		-k*sin(k*y)
		);
	J[1] =
	-(1/mu)*(
		(mu_0/mu)*(-2*x*(cos(k*x)+cos(k*y)))
		-k*sin(k*x)
		);
	J[2] = 0.0;

	return J;
}

inline Tensor<1, 3> current_vector_potential(double x, double y, double mu_0, double k)
{
	Tensor<1, 3> T;
	const double mu = permeability(x, y, mu_0);

	T[0] = 0.0;
	T[1] = 0.0;
	T[2] = (cos(k*x)+cos(k*y))/mu;

	return T;
}

inline Tensor<1, 3> vector_potential(double x, double y, double k)
{
	Tensor<1, 3> A;

	A[0] =-sin(k*y)/k;
	A[1] = sin(k*x)/k;
	A[2] = 0.0;

	return A;
}

inline Tensor<1, 3> magnetic_field(double x, double y, double k)
{
	Tensor<1, 3> B;

	B[0] = 0.0;
	B[1] = 0.0;
	B[2] = (cos(k*x)+cos(k*y));

	return B;
}

/**
 * \brief Describes exact solution, \f$\vec{J}_{f}\f$, of the
 * [Method of manufactured solutions, vector potential (mms-vt-i/)](@ref page_mms_vt_i)
 * numerical experiment.
 *****************************************************************************/
class ExactSolutionMMSVTI_Jf : public Function<3>, public SettingsMMSVTI
{
public:

	ExactSolutionMMSVTI_Jf();

	virtual void vector_value_list(const std::vector<Point<3>> & r,
		std::vector<Vector<double>>	 &values) const override final;
};

/**
 * \brief Describes exact solution, \f$\vec{B}\f$, of the
 * [Method of manufactured solutions, vector potential (mms-vt-i/)](@ref page_mms_vt_i)
 * numerical experiment.
 *****************************************************************************/
class ExactSolutionMMSVTI_B : public Function<3>, public SettingsMMSVTI
{
public:

	ExactSolutionMMSVTI_B();

	virtual void vector_value_list(const std::vector<Point<3>> & r,
		std::vector<Vector<double>>	 &values) const override final;
};

/**
 * \brief Describes the Dirichlet boundary condition for \f$\vec{T}\f$, in the
 * [Method of manufactured solutions, vector potential (mms-vt-i/)](@ref page_mms_vt_i)
 * numerical experiment.
 *****************************************************************************/
class DirichletBC_T : public Function<3>, public SettingsMMSVTI
{
public:

	DirichletBC_T();

	virtual void vector_value_list(const std::vector<Point<3>> & r,
		std::vector<Vector<double>>	 &values) const override final;
};

/**
 * \brief Describes the Dirichlet boundary condition for \f$\vec{A}\f$, in the
 * [Method of manufactured solutions, vector potential (mms-vt-i/)](@ref page_mms_vt_i)
 * numerical experiment.
 *****************************************************************************/
class DirichletBC_A : public Function<3>, public SettingsMMSVTI
{
public:

	DirichletBC_A();

	virtual void vector_value_list(const std::vector<Point<3>> & r,
		std::vector<Vector<double>>	 &values) const override final;
};

#endif


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

#ifndef SolverMMSVTI_H__
#define SolverMMSVTI_H__

#include <deal.II/base/vectorization.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>
//#include <deal.II/grid/manifold_lib.h>

#include <string>

#include "static_vector_solver_i.hpp"
#include "static_vector_solver_ii.hpp"
#include "exact_solution.hpp"
#include "settings.hpp"

using namespace StaticVectorSolver;

/**
 * \brief Implements the solver for current vector potential, \f$\vec{T}\f$, in the
 * [Method of manufactured solutions, vector potential (mms-vt-i/)](@ref page_mms_vt_i)
 * numerical experiment.
 *****************************************************************************/
class SolverMMSVTI_T : public SettingsMMSVTI, public Solver1<3,0>
{
public:

	SolverMMSVTI_T() = delete;

	SolverMMSVTI_T(
		unsigned int p,
		unsigned int mapping_degree,
		unsigned int r,
		std::string fname):
			Solver1<3,0>(
			p,
			mapping_degree,
			3,
			0.0,
			fname,
			nullptr,
			SettingsMMSVTI::print_time_tables,
			false),
			fname(fname),
			r(r)
	{
		Solver1<3,0>::run();
	}

	~SolverMMSVTI_T() = default;

private:

	const std::string fname;
	const unsigned int r;

//	SphericalManifold<3> sphere;

	DirichletBC_T dirichlet_bc;

	virtual void make_mesh() override final;
	virtual void fill_dirichlet_stack() override final;
	virtual void solve() override final;
};

/**
 * \brief Implements the solver for magnetic vector potential, \f$\vec{A}\f$, in the
 * [Method of manufactured solutions, vector potential (mms-vt-i/)](@ref page_mms_vt_i)
 * numerical experiment.
 *****************************************************************************/
class SolverMMSVTI_A : public SettingsMMSVTI, public Solver2<3,2>
{
public:

	SolverMMSVTI_A() = delete;

	SolverMMSVTI_A(
		unsigned int p,
		unsigned int mapping_degree,
		const Triangulation<3> & triangulation_T,
		const DoFHandler<3> & dof_handler_T,
		const Vector<double> & solution_T,
		unsigned int r,
		std::string fname):
			Solver2<3,2>(
			p,
			mapping_degree,
			triangulation_T,
			dof_handler_T,
			solution_T,
			0.0,
			fname,
			& dirichlet_bc,
			SettingsMMSVTI::print_time_tables,
			SettingsMMSVTI::project_exact_solution),
			fname(fname),
			r(r)
	{
		Solver2<3,2>::run();
	}

	~SolverMMSVTI_A() = default;

private:

	const std::string fname;

	const unsigned int r;

	DirichletBC_A dirichlet_bc;
	// Solver2 parasites on triangulations of other solvers.
	// No make_mesh() this time around.
	virtual void fill_dirichlet_stack() override final;
	virtual void solve() override final;
};

#endif


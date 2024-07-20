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

#ifndef SolverMMSVTII_H__
#define SolverMMSVTII_H__

#include <deal.II/base/vectorization.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>

#include <string>

#include "static_scalar_solver.hpp"
#include "static_vector_solver_ii.hpp"
#include "exact_solution.hpp"
#include "settings.hpp"

using namespace StaticScalarSolver;
using namespace StaticVectorSolver;

/**
 * \brief Implements the solver for current vector potential, \f$ T \f$, in the
 * [Method of manufactured solutions, vector potential (mms-vt-ii/)](@ref page_mms_vt_ii)
 * numerical experiment.
 *****************************************************************************/
class SolverMMSVTII_T : public SettingsMMSVTII, public Solver<2,0>
{
public:

	SolverMMSVTII_T() = delete;

	SolverMMSVTII_T(
		unsigned int p,
		unsigned int r,
		std::string fname):
			Solver<2,0>(
				p,
				1,
				2, // The right-hand side is volume free-current density.
				fname,
				& exact_solution,
				false, // Is axisymmetric.
				true, // Is vector potential.
				SettingsMMSVTII::print_time_tables,
				SettingsMMSVTII::project_exact_solution), // Project exact solution.
			fname(fname),
			r(r)
	{
		Solver<2,0>::run();
	}

	~SolverMMSVTII_T() = default;

private:

	const std::string fname;

	const unsigned int r;

	ExactSolutionMMSVTII_T exact_solution;
	DirichletBC_T dirichlet_bc;

	virtual void make_mesh() override final;
	virtual void fill_dirichlet_stack() override final;
	virtual void solve() override final;
};

/**
 * \brief Implements the solver for magnetic vector potential, \f$\vec{A}\f$, in the
 * [Method of manufactured solutions, vector potential (mms-vt-ii/)](@ref page_mms_vt_ii)
 * numerical experiment.
 *****************************************************************************/
class SolverMMSVTII_A : public SettingsMMSVTII, public Solver2<2,2>
{
public:

	SolverMMSVTII_A() = delete;

	SolverMMSVTII_A(
		unsigned int p,
		const Triangulation<2> & triangulation_T,
		const DoFHandler<2> & dof_handler_T,
		const Vector<double> & solution_T,
		unsigned int r,
		std::string fname):
			Solver2<2,2>(
			p,
			1,
			triangulation_T,
			dof_handler_T,
			solution_T,
			0.0,
			fname,
			nullptr,
			SettingsMMSVTII::print_time_tables,
			SettingsMMSVTII::project_exact_solution),
			fname(fname),
			r(r)
	{
		Solver2<2,2>::run();
	}

	~SolverMMSVTII_A() = default;

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


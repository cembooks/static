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

#ifndef SolverSSOLI_H__
#define SolverSSOLI_H__

#include <deal.II/base/vectorization.h>
#include <deal.II/base/tensor_function.h>

#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>

#include <string>

#include "static_vector_solver_i.hpp"
#include "exact_solution.hpp"
#include "settings.hpp"

using namespace StaticVectorSolver;

/**
 * \brief Solves for \f$\vec{A}\f$ in the
 * [Thin spherical coil (ssol-i/)](@ref page_ssol_i)
 * numerical experiment.
 *****************************************************************************/
class SolverSSOLI : public SettingsSSOLI, public Solver1<3>
{
public:

	SolverSSOLI() = delete;

	SolverSSOLI(
	unsigned int p,
	unsigned int mapping_degree,
	unsigned int r,
	std::string fname):
		Solver1<3>(
		p,
		mapping_degree,
		0,
		0.001/mu_0,
		fname,
		nullptr,
		SettingsSSOLI::print_time_tables,
		false),
		r(r),
		fname(fname)
	{
		Solver1<3>::run();
	}

	~SolverSSOLI() = default;

private:
	const unsigned int r;
	const std::string fname;

	virtual void make_mesh() override final;
	virtual void fill_dirichlet_stack() override final;
	virtual void solve() override final;

	void mark_materials();
};

#endif


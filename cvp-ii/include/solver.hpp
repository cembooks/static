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

#ifndef SolverCVP_H__
#define SolverCVP_H__

#include <deal.II/base/vectorization.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/grid/manifold_lib.h>

#include <string>

#include "static_scalar_solver.hpp"
#include "exact_solution.hpp"
#include "settings.hpp"

using namespace StaticScalarSolver;

/**
 * \brief Implements the solver of the
 * [Current vector potential (cvp-ii/)](@ref page_cvp_ii)
 * numerical experiment.
 *****************************************************************************/
class SolverCVPII : public SettingsCVPII, public Solver<2>
{
public:

	SolverCVPII() = delete;

/**
 * The constructor.
 *
 * @param[in] p - The degree of the interpolating polynomials of the Lagrange
 * finite elements,
 * [FE_Q](https://www.dealii.org/current/doxygen/deal.II/classFE__Q.html).
 * @param[in] mapping_degree - The degree of the interpolating polynomials used
 * for mapping. Setting it to 1 will do in the most of the cases. Note, that it
 * makes sense to attach a meaningful manifold to the triangulation if this
 * parameter is greater than 1.
 * @param[in] r - The parameter that encodes the degree of mesh refinement.
 * Must coincide with one of the values set in cvp-ii/gmsh/build. This parameter
 * is used to compose the name of the mesh file to be uploaded from
 * cvp-ii/gmsh/data/.
 * @param[in] fname - The name of the vtk file without extension to save
 * the data.
 *****************************************************************************/
	SolverCVPII(
	unsigned int p,
	unsigned int mapping_degree,
	unsigned int r,
	std::string fname):
		Solver<2>(
		p,
		mapping_degree,
		2,
		fname,
		& exact_solution,
		false,
		true,
		SettingsCVPII::print_time_tables,
		SettingsCVPII::project_exact_solution),
		fname(fname),
		r(r)
	{
		StaticScalarSolver::Solver<2>::run();
	}

	~SolverCVPII() = default;

private:

	const std::string fname;
	const unsigned int r;

	SphericalManifold<2> sphere;

	ExactSolutionCVPII_T exact_solution;

	const dealii::Functions::ZeroFunction<2> dirichlet_bc;
	virtual void make_mesh() override final;
	virtual void fill_dirichlet_stack() override final;
	virtual void solve() override final;
};

#endif


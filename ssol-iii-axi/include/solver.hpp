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

#ifndef SolverSSOLIIIAXI_H__
#define SolverSSOLIIIAXI_H__

#define BOOST_ALLOW_DEPRECATED_HEADERS

#include <deal.II/base/function.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/numerics/fe_field_function.h>

#include "static_scalar_solver.hpp"
#include "exact_solution.hpp"
#include "settings.hpp"

#define TMR(__name) \
	TimerOutput::Scope timer_section(timer, __name)

using namespace StaticScalarSolver;

/**
 * \brief Implements the
 * [Axisymmetric - thick spherical coil with magnetic core (ssol-iii-axi/)](@ref page_ssol_iii_axi)
 * numerical experiment.
 *****************************************************************************/
class SolverSSOLIIIAXI : public SettingsSSOLIIIAXI, public Solver<2>
{
public:

	SolverSSOLIIIAXI() = delete;

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
 * Must coincide with one of the values set in rho/gmsh/build. This parameter
 * is used to compose the name of the mesh file to be uploaded from
 * rho/gmsh/data/.
 * @param[in] fname - The name of the vtk file without extension to save
 * the data.
 *****************************************************************************/
	SolverSSOLIIIAXI(
	unsigned int p,
	unsigned int mapping_degree,
	unsigned int r,
	std::string fname
	):
		Solver<2>(
			p,
			mapping_degree,
			1,
			fname,
			nullptr,
			true,
			true,
			SettingsSSOLIIIAXI::print_time_tables,
			false),
		r(r),
		fname(fname)
	{
		Solver<2>::run();
	}

	~SolverSSOLIIIAXI() = default;

private:
	const unsigned int r;
	const std::string fname;

	const dealii::Functions::ZeroFunction<2> dirichlet_function;

	SphericalManifold<2> sphere;

	virtual void make_mesh() override final;
	virtual void fill_dirichlet_stack() override final;
	virtual void solve() override final;

	void mark_materials();
};

#endif


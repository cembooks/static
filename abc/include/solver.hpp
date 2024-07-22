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

#ifndef SolverABC_H__
#define SolverABC_H__

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
 * \brief Implements the solver of the
 * [Asymptotic boundary condition (abc/)](@ref page_abc)
 * numerical experiment.
 *****************************************************************************/
template<int dim>
class SolverABC : public SettingsABC, public Solver<dim>
{
public:

	SolverABC() = delete;

/**
 * The constructor.
 *
 * @param[in] m - The [factor](@ref abc_m_factor),\f$m\f$, that scales the
 * radius of the surface that represents the infinity.
 * @param[in] p - The degree of the interpolating polynomials of the Lagrange
 * finite elements,
 * [FE_Q](https://www.dealii.org/current/doxygen/deal.II/classFE__Q.html).
 * @param[in] mapping_degree - The degree of the interpolating polynomials used
 * for mapping. Setting it to 1 will do in the most of the cases. Note, that it
 * makes sense to attach a meaningful manifold to the triangulation if this
 * parameter is greater than 1.
 * @param[in] r - The parameter that encodes the degree of mesh refinement.
 * Must coincide with one of the values set in abc/gmsh/build. This parameter
 * is used to compose the name of the mesh file to be uploaded from
 * abc/gmsh/data/.
 * @param[in] fname - The name of the vtk file without extension to save
 * the data.
 *****************************************************************************/
	SolverABC(
	unsigned int m,
	unsigned int p,
	unsigned int mapping_degree,
	unsigned int r,
	std::string fname
	):
		Solver<dim>(
		p,
		mapping_degree,
		0,
		fname,
		& exact_solution,
		false,
		false,
		print_time_tables,
		project_exact_solution),
		m(m),
		r(r),
		fname(fname),
		dirichlet_function_left(-1.0),
		dirichlet_function_right(1.0),
		dirichlet_function_in(1.0),
		circle_left(Point<2>(-x0, 0.0)),
		circle_right(Point<2>(x0, 0.0)),
		circle_infty(Point<2>(0.0, 0.0))
	{
		R_infty = m * R_mid;
		TimerOutput::OutputFrequency tf =
			(print_time_tables) ? TimerOutput::summary : TimerOutput::never;

		TimerOutput timer(
			std::cout,
			tf,
			TimerOutput::cpu_and_wall_times_grouped);

		{TMR("Solver run"); Solver<dim>::run();}
	}

	~SolverABC() = default;

private:
	const unsigned int m;
	const unsigned int r;
	const std::string fname;

	const	ExactSolutionABC_PHI<dim> exact_solution;
	const Functions::ZeroFunction<dim> dirichlet_function_infty;
	const Functions::ConstantFunction<dim> dirichlet_function_left;
	const Functions::ConstantFunction<dim> dirichlet_function_right;
	const Functions::ConstantFunction<dim> dirichlet_function_in;

	const SphericalManifold<2> circle_left;
	const SphericalManifold<2> circle_right;
	const SphericalManifold<2> circle_infty;

	SphericalManifold<3> sphere;

 // Sets the ID of the outer boundary to 0, 2, 5 depending on the value of
 // BC_TYPE__ set CMakeLists.txt. The ID of the outer boundary in the geo
 // file equals 2.
	void renumber_boundaries();
	virtual void make_mesh() override final;
	virtual void fill_dirichlet_stack() override final;
	virtual void solve() override final;
};

template <int dim>
void SolverABC<dim>::renumber_boundaries()
{
	for (auto cell : Solver<dim>::triangulation.active_cell_iterators())
	{
		for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
		{
			if (cell->face(f)->at_boundary())
			{
				if (cell->face(f)->boundary_id() == 2)
					cell->face(f)->set_all_boundary_ids(bid_infty);
			}
		}
	}
}

template<int dim>
void SolverABC<dim>::solve()
{
	ReductionControl control(Solver<dim>::system_rhs.size(), 0.0, 1e-8, false, false);

	if (log_cg_convergence)
		control.enable_history_data();

	GrowingVectorMemory<Vector<double>> memory;
	SolverCG<Vector<double>> cg(control, memory);

	PreconditionJacobi<SparseMatrix<double>> preconditioner;
	preconditioner.initialize(Solver<dim>::system_matrix, 1.0);

	cg.solve(
		Solver<dim>::system_matrix,
		Solver<dim>::solution,
		Solver<dim>::system_rhs,
		preconditioner);

	Solver<dim>::constraints.distribute(Solver<dim>::solution);

	if (log_cg_convergence)
	{
		const std::vector<double> history_data = control.get_history_data();

		std::ofstream ofs(fname + "_cg_convergence.csv");

		unsigned int i = 1;
		for (auto item : history_data)
		{
			ofs << i << ", " << item  << "\n";
			i++;
		}
		ofs.close();
	}
}

#endif


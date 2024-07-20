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

#ifndef SolverMMSAXI_H__
#define SolverMMSAXI_H__

#define BOOST_ALLOW_DEPRECATED_HEADERS

#include <deal.II/base/function.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/numerics/fe_field_function.h>

#include "static_scalar_solver.hpp"
#include "exact_solution.hpp"
#include "settings.hpp"

#define TMR(__name) \
	TimerOutput::Scope timer_section(timer, __name)

using namespace StaticScalarSolver;

/**
 * \brief Implements the
 * [Axisymmetric - method of manufactured solutions](@ref page_mms_axi)
 * numerical experiment.
 *********************************************************/
template<int dim>
class SolverMMSAXI : public SettingsMMSAXI, public Solver<dim>
{
public:

	SolverMMSAXI() = delete;

/**
 * The constructor.
 *
 * @param[in] p - The degree of the interpolating polynomials of the Lagrange
 * finite elements,
 * [FE_Q](https://www.dealii.org/current/doxygen/deal.II/classFE__Q.html).
 * @param[in] r - The parameter that encodes the degree of mesh refinement.
 * Must coincide with one of the values set in mms-axi/gmsh/build. This parameter
 * is used to compose the name of the mesh file to be uploaded from
 * mms-axi/gmsh/data/.
 * @param[in] fname - The name of the vtk file without extension to save
 * the data.
 ***************************************************************************/
	SolverMMSAXI(
	unsigned int p,
	unsigned int r,
	std::string fname):
		Solver<dim>(
		p,
		1,
		1,
		fname,
		& exact_solution,
#if DIMENSION__ == 2
		true,
#endif
#if DIMENSION__ == 3
		false,
#endif
		false,
		print_time_tables,
		project_exact_solution),
		r(r),
		fname(fname)
	{
#if DIMENSION__ == 2
		fname_mesh_in = "../../gmsh/data/cylinder_2d_r"
			+ std::to_string(r) + ".msh";
		fname_mesh_out = "../../gmsh/data/cylinder_2d_r"
			+ std::to_string(r) + "_p" + std::to_string(p) + "_reordered.msh";
#endif
#if DIMENSION__ == 3
		fname_mesh_in = "../../gmsh/data/cylinder_3d_r"
			+ std::to_string(r) + ".msh";
		fname_mesh_out = "../../gmsh/data/cylinder_3d_r"
			+ std::to_string(r) + "_p" + std::to_string(p) + "_reordered.msh";
#endif

		TimerOutput::OutputFrequency tf =
			(print_time_tables) ? TimerOutput::summary : TimerOutput::never;

		TimerOutput timer(
			std::cout,
			tf,
			TimerOutput::cpu_and_wall_times_grouped);

		{TMR("Solver run"); Solver<dim>::run();}
	}

	~SolverMMSAXI() = default;

private:

	const unsigned int r;
	const std::string fname;
	const ExactSolutionMMSAXI_PHI<dim> exact_solution;

	std::string fname_mesh_in;
	std::string fname_mesh_out;

	virtual void make_mesh() override final;
	virtual void fill_dirichlet_stack() override final;
	virtual void solve() override final;
};

template<int dim>
void SolverMMSAXI<dim>::fill_dirichlet_stack()
{
	Solver<dim>::dirichlet_stack = {{bid_dirichlet, & exact_solution}};
}

template <int dim >
void SolverMMSAXI<dim>::make_mesh()
{
	GridIn<dim> gridin;
	Triangulation<dim> tria_tmp;

	gridin.attach_triangulation(tria_tmp);
	std::ifstream ifs(fname_mesh_in);
	gridin.read_msh(ifs);

	std::tuple< std::vector< Point<dim>>, std::vector< CellData<dim> >, SubCellData> mesh_description;

	mesh_description = GridTools::get_coarse_mesh_description(tria_tmp);

	GridTools::invert_all_negative_measure_cells(
		std::get<0>(mesh_description),
		std::get<1>(mesh_description));

	GridTools::consistently_order_cells(std::get<1>(mesh_description));

	Solver<dim>::triangulation.create_triangulation(
		std::get<0>(mesh_description),
		std::get<1>(mesh_description),
		std::get<2>(mesh_description));

	GridOut gridout;
	GridOutFlags::Msh msh_flags(true, true);
	gridout.set_flags(msh_flags);

	std::ofstream ofs(fname_mesh_out);
	gridout.write_msh(Solver<dim>::triangulation, ofs);
}

template<int dim>
void SolverMMSAXI<dim>::solve()
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


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

#ifndef SolverMMS_H__
#define SolverMMS_H__

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include "static_scalar_solver.hpp"
#include "exact_solution.hpp"
#include "settings.hpp"

using namespace StaticScalarSolver;

/**
 * \brief Implements the solver of the
 * [Method of manufactured solutions (mms/)](@ref page_mms)
 * numerical experiment.
 *****************************************************************************/
template<int dim>
class SolverMMS : public SettingsMMS, public Solver<dim>
{
public:

	SolverMMS() = delete;

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
 * Must coincide with one of the values set in mms/gmsh/build. This parameter
 * is used to compose the name of the mesh file to be uploaded from
 * mms/gmsh/data/.
 * @param[in] fname - The name of the vtk file without extension to save
 * the data.
 *****************************************************************************/
	SolverMMS(
	unsigned int p,
	unsigned int mapping_degree,
	unsigned int r,
	std::string fname):
		Solver<dim>(
		p,
		mapping_degree,
		1, // The PDE right-hand side is free-current density.
		fname,
		& exact_solution,
		false, // Is axisymmetric.
		false, // Is vector potential.
		print_time_tables,
		project_exact_solution),
		fname(fname)
	{
		if (DIMENSION__ == 2)
		{
			if (HYPERCUBE__ == 1)
			{
				fname_mesh_in = "../../gmsh/data/square_r"
					+ std::to_string(r) + ".msh";
				fname_mesh_out = "../../gmsh/data/square_r"
					+ std::to_string(r) + "_p" + std::to_string(p) + "_reordered.msh";
			}else
			{
				fname_mesh_in = "../../gmsh/data/circle_r"
					+ std::to_string(r) + ".msh";
					fname_mesh_out = "../../gmsh/data/circle_r"
					+ std::to_string(r) + "_p" + std::to_string(p) + "_reordered.msh";
			}
		}else
		{
			if (HYPERCUBE__ == 1)
			{
				fname_mesh_in = "../../gmsh/data/cube_r"
					+ std::to_string(r) + ".msh";
					fname_mesh_out = "../../gmsh/data/cube_r"
					+ std::to_string(r) + "_p" + std::to_string(p) + "_reordered.msh";
			}else
			{
				fname_mesh_in = "../../gmsh/data/sphere_r"
					+ std::to_string(r) + ".msh";
					fname_mesh_out = "../../gmsh/data/sphere_r"
					+ std::to_string(r) + "_p" + std::to_string(p) + "_reordered.msh";
			}
		}

		Solver<dim>::run();
	}

	~SolverMMS() = default;

private:

	const std::string fname;

	std::string fname_mesh_in;
	std::string fname_mesh_out;

	const ExactSolutionMMS_PHI<dim> exact_solution;

	virtual void make_mesh() override final;
	virtual void fill_dirichlet_stack() override final;
	virtual void solve() override final;
};

template<int dim>
void SolverMMS<dim>::fill_dirichlet_stack()
{
	Solver<dim>::dirichlet_stack = {{bid_dirichlet, & exact_solution}};
}

template<int dim>
void SolverMMS<dim>::solve()
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

template <int dim >
void SolverMMS<dim>::make_mesh()
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
#endif


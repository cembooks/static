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

#ifndef SolverMMSV_H__
#define SolverMMSV_H__

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <string>

#include "static_vector_solver_i.hpp"
#include "exact_solution.hpp"
#include "settings.hpp"

using namespace StaticVectorSolver;

/**
 * \brief Implements the
 * [Method of manufactured solutions, vector potential (mms-v/)](@ref page_mms_v)
 * numerical experiment.
 *****************************************************************************/
template<int dim>
class SolverMMSV : public SettingsMMSV, public Solver1<dim>
{
public:

	SolverMMSV() = delete;

/**
 * The constructor.
 *
 * @param[in] p - The degree of the Nedelec finite elements.
 * @param[in] r - The parameter that encodes the degree of mesh refinement.
 * Must coincide with one of the values set in mms-v/gmsh/build. This parameter
 * is used to compose the name of the mesh file to be uploaded from
 * mms-v/gmsh/data/.
 * @param[in] fname - The name of the vtk file without extension to save
 * the data.
 *****************************************************************************/
	SolverMMSV(
	unsigned int p,
	unsigned int r,
	std::string fname):
		Solver1<dim>(
		p,
		1,
		3,
		0.0, // 0.01/mu_0,
		fname,
		& exact_solution,
		print_time_tables,
		project_exact_solution),
		r(r),
		fname(fname)
	{
		if (DIMENSION__ == 2)
		{
			if (HYPERCUBE__ == 1)
			{
				fname_mesh = "../../gmsh/data/square_r"
					+ std::to_string(r) + ".msh";
			}else
			{
				fname_mesh = "../../gmsh/data/circle_r"
					+ std::to_string(r) + ".msh";
			}
		}else
		{
			if (HYPERCUBE__ == 1)
			{
				fname_mesh = "../../gmsh/data/cube_r"
					+ std::to_string(r) + ".msh";
			}else
			{
				fname_mesh = "../../gmsh/data/sphere_r"
					+ std::to_string(r) + ".msh";
			}
		}

		Solver1<dim>::run();
	}

	~SolverMMSV() = default;

private:

	const unsigned int r;
	const std::string fname;
	std::string fname_mesh;

	const ExactSolutionMMSV_A<dim> exact_solution;

	virtual void make_mesh() override final;
	virtual void fill_dirichlet_stack() override final;
	virtual void solve() override final;
};

template<int dim>
void SolverMMSV<dim>::fill_dirichlet_stack()
{
	Solver1<dim>::dirichlet_stack = {{bid_dirichlet, & exact_solution}};
}

template<int dim>
void SolverMMSV<dim>::solve()
{
	ReductionControl control(Solver1<dim>::system_rhs.size(), 0.0, 1e-12, false, false);

	if (log_cg_convergence)
		control.enable_history_data();

	GrowingVectorMemory<Vector<double>> memory;
	SolverCG<Vector<double>> cg(control, memory);

//	PreconditionJacobi<SparseMatrix<double>> preconditioner;
//	preconditioner.initialize(Solver1<dim>::system_matrix, 1.0);

	PreconditionSSOR<SparseMatrix<double>> preconditioner;
	preconditioner.initialize(Solver1<dim>::system_matrix, 1.2);

	cg.solve(
		Solver1<dim>::system_matrix,
		Solver1<dim>::solution,
		Solver1<dim>::system_rhs,
		preconditioner);

	Solver1<dim>::constraints.distribute(Solver1<dim>::solution);

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
void SolverMMSV<dim>::make_mesh()
{
	GridIn<dim> gridin;
	gridin.attach_triangulation(Solver1<dim>::triangulation);
	
	std::ifstream ifs(fname_mesh);
	gridin.read_msh(ifs);
}
#endif


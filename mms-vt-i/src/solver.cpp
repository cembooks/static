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

#define BOOST_ALLOW_DEPRECATED_HEADERS

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include "solver.hpp"
#include <fstream>

using namespace StaticVectorSolver;

void SolverMMSVTI_T::make_mesh()
{
	GridIn<3> gridin;

	gridin.attach_triangulation(Solver1<3,0>::triangulation);
#if DOMAIN__ == 0
	std::ifstream ifs("../../gmsh/data/cube_r" + std::to_string(r) + ".msh");
#elif DOMAIN__ == 1
	std::ifstream ifs("../../gmsh/data/sphere_r" + std::to_string(r) + ".msh");
#endif
	gridin.read_msh(ifs);
}

void SolverMMSVTI_T::fill_dirichlet_stack()
{
	Solver1<3,0>::dirichlet_stack =
		{{bid_dirichlet, & dirichlet_bc}};
}

void SolverMMSVTI_T::solve()
{
	ReductionControl control(
		Solver1<3,0>::system_rhs.size(), 0.0, 1e-8, false, false);

	if (log_cg_convergence)
		control.enable_history_data();

	GrowingVectorMemory<Vector<double>> memory;
	SolverCG<Vector<double>> cg(control, memory);

	PreconditionSSOR<SparseMatrix<double>> preconditioner;
	preconditioner.initialize(Solver1<3,0>::system_matrix, 1.2);

	cg.solve(
		Solver1<3,0>::system_matrix,
		Solver1<3,0>::solution,
		Solver1<3,0>::system_rhs,
		preconditioner);

	Solver1<3,0>::constraints.distribute(Solver1<3,0>::solution);

	if (SettingsMMSVTI::log_cg_convergence)
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

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

void SolverMMSVTI_A::fill_dirichlet_stack()
{
	Solver2<3,2>::dirichlet_stack =
		{{bid_dirichlet, & dirichlet_bc}};
}

void SolverMMSVTI_A::solve()
{
	ReductionControl control(
		Solver2<3,2>::system_rhs.size(), 0.0, 1e-8, false, false);

	if (log_cg_convergence)
		control.enable_history_data();

	GrowingVectorMemory<Vector<double>> memory;
	SolverCG<Vector<double>> cg(control, memory);

	PreconditionSSOR<SparseMatrix<double>> preconditioner;
	preconditioner.initialize(Solver2<3,2>::system_matrix, 1.2);

	cg.solve(
		Solver2<3,2>::system_matrix,
		Solver2<3,2>::solution,
		Solver2<3,2>::system_rhs,
		preconditioner);

	Solver2<3,2>::constraints.distribute(Solver2<3,2>::solution);

	if (SettingsMMSVTI::log_cg_convergence)
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


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
	Triangulation<3> tria_tmp;

	gridin.attach_triangulation(tria_tmp);
#if DOMAIN__ == 0
	std::ifstream ifs("../../gmsh/data/cube_r" + std::to_string(r) + ".msh");
#elif DOMAIN__ == 1
	std::ifstream ifs("../../gmsh/data/sphere_r" + std::to_string(r) + ".msh");
#endif
	gridin.read_msh(ifs);

	std::tuple< std::vector< Point<3>>, std::vector< CellData<3> >, SubCellData> mesh_description;

	mesh_description = GridTools::get_coarse_mesh_description(tria_tmp);

	GridTools::invert_all_negative_measure_cells(
		std::get<0>(mesh_description),
		std::get<1>(mesh_description));

	GridTools::consistently_order_cells(std::get<1>(mesh_description));

	Solver1<3,0>::triangulation.create_triangulation(
		std::get<0>(mesh_description),
		std::get<1>(mesh_description),
		std::get<2>(mesh_description));
/*
#if DOMAIN__ == 1
	for (auto cell : Solver1<3,0>::triangulation.active_cell_iterators())
	{
		for (unsigned int f = 0; f < GeometryInfo<3>::faces_per_cell; f++)
		{
			double dif_norm_a = 0.0;

			for (unsigned int v = 0; v < GeometryInfo<3>::vertices_per_face; v++ )
			{
				dif_norm_a  += std::abs(cell->face(f)->vertex(v).norm() - a);
			}

			if ( dif_norm_a < eps )
				cell->face(f)->set_all_manifold_ids(1);
		}
	}

	Solver1<3,0>::triangulation.set_manifold(1,sphere);
#endif
*/
	GridOut gridout;
	GridOutFlags::Msh msh_flags(true, true);
	gridout.set_flags(msh_flags);

#if DOMAIN__ == 0
	std::ofstream ofs("../../gmsh/data/cube_r" + std::to_string(r) + "_reordered.msh");
#elif DOMAIN__ == 1
	std::ofstream ofs("../../gmsh/data/sphere_r" + std::to_string(r) + "_reordered.msh");
#endif
	gridout.write_msh(Solver1<3,0>::triangulation, ofs);
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


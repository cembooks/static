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

#include "solver.hpp"

void SolverSSOLIIIAXI::make_mesh()
{
	GridIn<2> gridin;
	gridin.attach_triangulation(triangulation);
	std::ifstream ifs("../../gmsh/data/circle_r" + std::to_string(r) + ".msh");
	gridin.read_msh(ifs);
	mark_materials();
}

void SolverSSOLIIIAXI::fill_dirichlet_stack()
{
		Solver<2>::dirichlet_stack =
		{{bid_axi, & dirichlet_function}};
}

void SolverSSOLIIIAXI::mark_materials()
{
	Solver<2>::triangulation.reset_all_manifolds();

	for (auto cell : Solver<2>::triangulation.active_cell_iterators())
	{
		cell->set_material_id(mid_1);

		if ((cell->center().norm() < b1) &&
				(cell->center().norm() > a1))
			cell->set_material_id(mid_2);

		if ((cell->center().norm() < b2) &&
				(cell->center().norm() > a2))
			cell->set_material_id(mid_3);

		for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; f++)
			if ( std::abs(cell->face(f)->vertex(1).norm()-
					cell->face(f)->vertex(0).norm()) < eps )
				cell->face(f)->set_all_manifold_ids(1);
	}

	Solver<2>::triangulation.set_manifold(1,sphere);
}

void SolverSSOLIIIAXI::solve()
{
	ReductionControl control(Solver<2>::system_rhs.size(), 0.0, 1e-12, false, false);

	if (log_cg_convergence)
		control.enable_history_data();

	GrowingVectorMemory<Vector<double>> memory;
	SolverCG<Vector<double>> cg(control, memory);

	PreconditionJacobi<SparseMatrix<double>> preconditioner;
	preconditioner.initialize(Solver<2>::system_matrix, 1.0);

	cg.solve(
		Solver<2>::system_matrix,
		Solver<2>::solution,
		Solver<2>::system_rhs,
		preconditioner);

	Solver<2>::constraints.distribute(Solver<2>::solution);

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


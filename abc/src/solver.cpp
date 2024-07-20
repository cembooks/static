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

template <>
void SolverABC<2>::make_mesh()
{
	GridIn<2> gridin;
		gridin.attach_triangulation(Solver<2>::triangulation);

	std::ifstream ifs("../../gmsh/data/ppc_m"
		+ std::to_string(m) + "_r"
		+ std::to_string(r) + ".msh");

	gridin.read_msh(ifs);

	renumber_boundaries();

	Solver<2>::triangulation.reset_all_manifolds();

	for (auto cell : Solver<2>::triangulation.active_cell_iterators())
	{
		for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; ++f)
		{
			if ( cell->face(f)->at_boundary() )
			{
				if ( cell->face(f)->boundary_id() == bid_left)
						cell->face(f)->set_all_manifold_ids(mfid_left);

				if ( cell->face(f)->boundary_id() == bid_right)
						cell->face(f)->set_all_manifold_ids(mfid_right);

				if ( cell->face(f)->boundary_id() == bid_infty)
						cell->face(f)->set_all_manifold_ids(mfid_infty);
			}
		}
	}

	Solver<2>::triangulation.set_manifold(mfid_left, circle_left);
	Solver<2>::triangulation.set_manifold(mfid_right, circle_right);
	Solver<2>::triangulation.set_manifold(mfid_infty, circle_infty);
}

template <>
void SolverABC<3>::make_mesh()
{
	GridIn<3> gridin;
		gridin.attach_triangulation(Solver<3>::triangulation);

	std::ifstream ifs("../../gmsh/data/shell_m"
		+ std::to_string(m) + "_r"
		+ std::to_string(r) + ".msh");

	gridin.read_msh(ifs);

	renumber_boundaries();

	Solver<3>::triangulation.set_all_manifold_ids(mfid_infty);
	Solver<3>::triangulation.set_manifold(mfid_infty, sphere);
}

template<>
void SolverABC<2>::fill_dirichlet_stack()
{
// See also #if BC_TYPE__ preprocessor directives
// in abc/include/settings.hpp
#if BC_TYPE__ == 0
// Neumann boundary condition on the outer boundary.
	dirichlet_stack =
		{{bid_left, & dirichlet_function_left},
		 {bid_right, & dirichlet_function_right}};
#endif
#if BC_TYPE__ == 1
// Dirichlet boundary condition on the outer boundary.
	dirichlet_stack =
		{{bid_left, & dirichlet_function_left},
		 {bid_right, & dirichlet_function_right},
		 {bid_infty, & dirichlet_function_infty}};
#endif
#if BC_TYPE__ == 2
// ABC on the outer boundary.
	dirichlet_stack =
		{{bid_left, & dirichlet_function_left},
		 {bid_right, & dirichlet_function_right}};
#endif
}

template<>
void SolverABC<3>::fill_dirichlet_stack()
{
// See also #if BC_TYPE__ preprocessor directives
// in abc/include/settings.hpp
#if BC_TYPE__ == 0
// Neumann boundary condition on the outer boundary.
	dirichlet_stack =
		{{bid_in, & dirichlet_function_in}};
#endif
#if BC_TYPE__ == 1
// Dirichlet boundary condition on the outer boundary.
	dirichlet_stack =
		{{bid_in, & dirichlet_function_in},
		 {bid_infty, & dirichlet_function_infty}};
#endif
#if BC_TYPE__ == 2
// ABC boundary condition on the outer boundary.
	dirichlet_stack =
		{{bid_in, & dirichlet_function_in}};
#endif
}


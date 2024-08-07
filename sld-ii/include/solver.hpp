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

#ifndef SolverSLDII_H__
#define SolverSLDII_H__

#define BOOST_ALLOW_DEPRECATED_HEADERS

#include <deal.II/base/function.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
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
 * [Magnetostatic shield - 2 (sld-ii/)](@ref page_sld_ii)
 * numerical experiment.
 *****************************************************************************/
template<int dim>
class SolverSLDII : public SettingsSLDII, public Solver<dim>
{
public:

	SolverSLDII() = delete;
/**
 * The constructor.

 * @param[in] p - The degree of the interpolating polynomials of the Lagrange
 * finite elements,
 * [FE_Q](https://www.dealii.org/current/doxygen/deal.II/classFE__Q.html).
 * @param[in] mapping_degree - The degree of the interpolating polynomials used
 * for mapping. Setting it to 1 will do in the most of the cases. Note, that it
 * makes sense to attach a meaningful manifold to the triangulation if this
 * parameter is greater than 1.
 * @param[in] r - The parameter that encodes the degree of mesh refinement.
 * Must coincide with one of the values set in sld-ii/gmsh/build. This parameter
 * is used to compose the name of the mesh file to be uploaded from
 * sld-ii/gmsh/data/.
 * @param[in] fname - The name of the vtk file without extension to save
 * the data.
 *****************************************************************************/
	SolverSLDII(
	unsigned int p,
	unsigned int mapping_degree,
	unsigned int r,
	std::string fname):
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
		r(r),
		fname(fname)
	{
		Solver<dim>::run();
	}

	~SolverSLDII() = default;

private:
	const unsigned int r;
	const std::string fname;

	const	ExactSolutionSLDII_THETA<dim> exact_solution;
	const Functions::ZeroFunction<dim> dirichlet;

	SphericalManifold<dim> sphere;

	virtual void make_mesh() override final;
	virtual void fill_dirichlet_stack() override final;
	virtual void solve() override final;

	void mark_materials();
};

template<int dim>
void SolverSLDII<dim>::fill_dirichlet_stack()
{
//	Solver<dim>::dirichlet_stack = {{bid, & dirichlet}};
	Solver<dim>::dirichlet_stack = {{bid, & exact_solution}};
}

template <int dim >
void SolverSLDII<dim>::make_mesh()
{
	GridIn<dim> gridin;
	Triangulation<dim> tria_tmp;
	gridin.attach_triangulation(tria_tmp);

	std::string fname_mesh_in = (dim == 2) ?
		"../../gmsh/data/square_r"+std::to_string(r)+".msh" :
		"../../gmsh/data/cube_r"+std::to_string(r)+".msh";

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

	mark_materials();

	GridOut gridout;
	GridOutFlags::Msh msh_flags(true, true);
	gridout.set_flags(msh_flags);

	std::string fname_mesh_out = (dim == 2) ?
		"../../gmsh/data/square_reordered_r"+std::to_string(r)+".msh" :
		"../../gmsh/data/cube_reordered_r."+std::to_string(r)+".msh";

	std::ofstream ofs(fname_mesh_out);
	gridout.write_msh(Solver<dim>::triangulation, ofs);
}

template<int dim>
void SolverSLDII<dim>::mark_materials()
{
	Solver<dim>::triangulation.reset_all_manifolds();

	double cell_r;
	for (auto cell : Solver<dim>::triangulation.active_cell_iterators())
	{
		cell_r = cell->center().norm();

		if ( cell_r < a )
		{
			// The cell is inside the shield.
			cell->set_material_id(mid_1);
		}else if ( cell_r < b )
		{
			// The cell is in the wall of the shield.
			cell->set_material_id(mid_2);
		}else
		{
			// The cell is outside the shield.
			cell->set_material_id(mid_3);
		}

		for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
		{
			double dif_norm_a = 0.0;
			double dif_norm_b = 0.0;
			for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face; v++ )
			{
				dif_norm_a += std::abs(cell->face(f)->vertex(v).norm()-a);
				dif_norm_b += std::abs(cell->face(f)->vertex(v).norm()-b);
			}

			if 	( dif_norm_a < eps )
			{
				cell->face(f)->set_all_manifold_ids(0);

				if ( std::abs(cell->center().norm()) < a )
				{
					// The face belongs to the interface Gamma_1
					cell->face(f)->set_user_index(1);
					cell->set_user_index(1);
				}
			}

			if 	( dif_norm_b < eps )
			{
				cell->face(f)->set_all_manifold_ids(0);

				if ( std::abs(cell->center().norm()) < b )
				{
					// The face belongs to the interface Gamma_2
					cell->face(f)->set_user_index(2);
					cell->set_user_index(2);
				}
			}
		}
	}

	Solver<dim>::triangulation.set_manifold(0,sphere);
}

template<int dim>
void SolverSLDII<dim>::solve()
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


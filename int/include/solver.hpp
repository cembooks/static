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

#ifndef SolverINT_H__
#define SolverINT_H__

#define BOOST_ALLOW_DEPRECATED_HEADERS

#include <deal.II/base/function.h>
#include <deal.II/grid/manifold_lib.h>
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
 * [Interface between dielectrics](@ref page_int)
 * numerical experiment.
 *****************************************************************************/
template<int dim>
class SolverINT : public SettingsINT, public Solver<dim>
{
public:

	SolverINT() = delete;

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
 * Must coincide with one of the values set in int/gmsh/build. This parameter
 * is used to compose the name of the mesh file to be uploaded from
 * int/gmsh/data/.
 * @param[in] fname - The name of the vtk file without extension to save
 * the data.
 *****************************************************************************/
	SolverINT(
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
				fname(fname),
				dirichlet_function_in(1.0),
				fe_slice(1)
	{
		if (DIMENSION__ == 2)
		{
			fname_mesh_in = "../../gmsh/data/ring_r"
				+ std::to_string(r) + ".msh";
			fname_mesh_out = "../../gmsh/data/ring_r"
				+ std::to_string(r) + "_p" + std::to_string(p) + "_reordered.msh";
		}else
		{
			fname_mesh_in = "../../gmsh/data/shell_r"
				+ std::to_string(r) + ".msh";
			fname_mesh_out = "../../gmsh/data/shell_r"
				+ std::to_string(r) + "_p" + std::to_string(p) + "_reordered.msh";
		}

		TimerOutput::OutputFrequency tf =
			(print_time_tables) ? TimerOutput::summary : TimerOutput::never;

		TimerOutput timer(
			std::cout,
			tf,
			TimerOutput::cpu_and_wall_times_grouped);

		{TMR("Solver run"); Solver<dim>::run();}
		{TMR("Data slice");	data_slice(fname);}
	}

	~SolverINT() = default;

private:
	const unsigned int r;
	const std::string fname;
	std::string fname_mesh_in;
	std::string fname_mesh_out;

	const	ExactSolutionINT_PHI<dim> exact_solution;
	const dealii::Functions::ZeroFunction<dim> dirichlet_function_out;
	const dealii::Functions::ConstantFunction<dim> dirichlet_function_in;

//The amount of global mesh refinements that need to be done to the
//one-dimensional mesh used for the plot of potential vs. x coordinate.
	const unsigned int nr_slice_global_refs = 10;

// These four data members are needed for making the plot of potential
// vs. x coordinate.
	Triangulation<1,dim> triangulation_slice;
	FE_Q<1,dim> fe_slice;
	DoFHandler<1,dim> dof_handler_slice;
	Vector<double> solution_slice;

	SphericalManifold<dim> sphere;

	virtual void make_mesh() override final;
	virtual void fill_dirichlet_stack() override final;
	virtual void solve() override final;

	void mark_materials();

// This function makes the plot of potential vs. x coordinate.
	void data_slice(std::string fname);
};

template< int dim>
void SolverINT<dim>::fill_dirichlet_stack()
{
	Solver<dim>::dirichlet_stack =
		{{bid_in, & dirichlet_function_in},
		 {bid_out, & dirichlet_function_out}};
}

template <int dim >
void SolverINT<dim>::make_mesh()
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

	Point<dim> origin;

	for (auto cell : Solver<dim>::triangulation.active_cell_iterators())
	{
		if (cell->center().norm() < d)
			cell->set_material_id(mid_1);

		if (cell->center().norm() > d)
			cell->set_material_id(mid_2);

		for (unsigned int f = 0; f <  GeometryInfo<dim>::faces_per_cell; ++f)
		{
			if (cell->face(f)->at_boundary())
			{
				if (cell->center().norm() > d  )
				{
					cell->face(f)->set_all_boundary_ids(bid_out);
				}
				else
				{
					cell->face(f)->set_all_boundary_ids(bid_in);
				}
			}
		}
	}

	GridOut gridout;
	GridOutFlags::Msh msh_flags(true, true);
	gridout.set_flags(msh_flags);

	std::ofstream ofs(fname_mesh_out);
	gridout.write_msh(Solver<dim>::triangulation, ofs);

	Solver<dim>::triangulation.set_all_manifold_ids(0);
	Solver<dim>::triangulation.set_manifold(0,sphere);
}

template <int dim>
void SolverINT<dim>::data_slice(std::string fname)
{
	GridGenerator::hyper_cube(triangulation_slice, a + eps, b - eps);
	triangulation_slice.refine_global(nr_slice_global_refs);

	dof_handler_slice.reinit(triangulation_slice);
	dof_handler_slice.distribute_dofs(fe_slice);
	solution_slice.reinit(dof_handler_slice.n_dofs());

	Functions::FEFieldFunction<dim> potential(
		Solver<dim>::dof_handler,
		Solver<dim>::solution);

	VectorTools::interpolate(dof_handler_slice, potential, solution_slice);

//	DataOut<1,DoFHandler<1,dim>> data_out;
	DataOut<1, dim> data_out;

	data_out.attach_dof_handler(dof_handler_slice);
	data_out.add_data_vector(solution_slice, "solution_slice");
	data_out.build_patches();

	std::ofstream out( fname + "_slice" + ".gpi");

	data_out.write_gnuplot(out);
	out.close();
}

template<int dim>
void SolverINT<dim>::solve()
{
	ReductionControl control(Solver<dim>::system_rhs.size(), 0.0, 1e-12, false, false);

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


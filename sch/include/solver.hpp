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

#ifndef SolverSCH_H__
#define SolverSCH_H__

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
 * \brief Implements the
 * [Surface charge](@ref page_sch)
 * numerical experiment.
 *****************************************************************************/
template<int dim>
class SolverSCH : public SettingsSCH, public Solver<dim>
{
public:

	SolverSCH() = delete;

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
 * Must coincide with one of the values set in sch/gmsh/build. This parameter
 * is used to compose the name of the mesh file to be uploaded from
 * sch/gmsh/data/.
 * @param[in] fname - The name of the vtk file without extension to save
 * the data.
 *****************************************************************************/
	SolverSCH(
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
		r(r),
		fname(fname),
		fe_slice(1)
	{
		TimerOutput::OutputFrequency tf =
			(print_time_tables) ? TimerOutput::summary : TimerOutput::never;

		TimerOutput timer(
			std::cout,
			tf,
			TimerOutput::cpu_and_wall_times_grouped);

		{TMR("Solver run"); Solver<dim>::run();}
		{TMR("Data slice");	data_slice(fname);}
	}

	~SolverSCH() = default;

private:
	const unsigned int r;
	const std::string fname;

	const ExactSolutionSCH_PHI<dim> exact_solution;
	const dealii::Functions::ZeroFunction<dim> dirichlet_function;

	//The amount of global mesh refinements that need to be done to the
	//one-dimensional mesh used for the plot of potential vs. x coordinate.
	const unsigned int nr_slice_global_refs = 10;

	// These four data members are needed for making the plot of potential
	// vs. \f$x\f$ coordinate.
	Triangulation<1,dim> triangulation_slice;
	FE_Q<1,dim> fe_slice;
	DoFHandler<1,dim> dof_handler_slice;
	Vector<double> solution_slice;

	SphericalManifold<dim> sphere;

	virtual void make_mesh() override final;
	virtual void fill_dirichlet_stack() override final;
	virtual void solve() override final;

	void mark_materials();

	// This function makes the plot of potential vs. \f$x\f$ coordinate.
	void data_slice(std::string fname);
};

template< int dim>
void SolverSCH<dim>::fill_dirichlet_stack()
{
		Solver<dim>::dirichlet_stack =
		{{bid, & dirichlet_function}};
}

template <int dim>
void SolverSCH<dim>::mark_materials()
{
	Solver<dim>::triangulation.reset_all_manifolds();

	for (auto cell : Solver<dim>::triangulation.active_cell_iterators())
	{
		if ( std::abs(cell->center().norm()) < a )
		{
			for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; f++)
			{
				double dif_norm = 0.0;
				for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face; v++ )
				{
					dif_norm += std::abs(cell->face(f)->vertex(v).norm()-a);
				}

				if ( dif_norm < eps )
				{
					cell->face(f)->set_user_index(1);
					cell->set_user_index(1);

					cell->face(f)->set_all_manifold_ids(1);
				}
			}
		}
	}

	Solver<dim>::triangulation.set_all_manifold_ids_on_boundary(1);
	Solver<dim>::triangulation.set_manifold(1,sphere);
}

template <int dim>
void SolverSCH<dim>::data_slice(std::string fname)
{
	GridGenerator::hyper_cube(triangulation_slice, 0.0 + eps, b - eps);
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
void SolverSCH<dim>::solve()
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


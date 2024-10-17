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

#ifndef SolverSLDI_H__
#define SolverSLDI_H__

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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
/**
 * \brief The Dirichlet boundary condition of the
 * [Magnetostatic shield - 1 (sld-i/)](@ref page_sld_i)
 * numerical experiment. Applied to the outer (the only) boundary of the problem
 * domain.
 *****************************************************************************/
template<int dim>
class DirichletSLDI : public SettingsSLDI, public Function<dim>
{
public:
	DirichletSLDI(){};

	virtual double value(const Point<dim> & r,
		const unsigned int component = 0) const override final;
};

#pragma GCC diagnostic pop

/**
 * \brief Implements the solver that solves for \f$\Psi\f$ in the
 * [Magnetostatic shield - 1 (sld-i/)](@ref page_sld_i)
 * numerical experiment.
 *****************************************************************************/
template<int dim>
class SolverSLDI : public SettingsSLDI, public Solver<dim>
{
public:

	SolverSLDI() = delete;
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
 * Must coincide with one of the values set in sld-i/gmsh/build. This parameter
 * is used to compose the name of the mesh file to be uploaded from
 * sld-i/gmsh/data/.
 * @param[in] fname - The name of the vtk file without extension to save
 * the data.
 *****************************************************************************/
	SolverSLDI(
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

	~SolverSLDI() = default;

private:
	const unsigned int r;
	const std::string fname;

	const	ExactSolutionSLDI_PSI<dim> exact_solution;
	const DirichletSLDI<dim> dirichlet;

	SphericalManifold<dim> sphere;

	virtual void make_mesh() override final;
	virtual void fill_dirichlet_stack() override final;
	virtual void solve() override final;
};

template<int dim>
void SolverSLDI<dim>::fill_dirichlet_stack()
{
//	Solver<dim>::dirichlet_stack = {{bid, & dirichlet}};
	Solver<dim>::dirichlet_stack = {{bid, & exact_solution}};
}

template <int dim >
void SolverSLDI<dim>::make_mesh()
{
	GridIn<dim> gridin;
	gridin.attach_triangulation(Solver<dim>::triangulation);

	std::string fname_mesh = (dim == 2) ?
		"../../gmsh/data/square_r"+std::to_string(r)+".msh" :
		"../../gmsh/data/cube_r"+std::to_string(r)+".msh";

	std::ifstream ifs(fname_mesh);
	gridin.read_msh(ifs);
	
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
			// The cell is inside the wall of the shield.
			cell->set_material_id(mid_2);
		}else
		{
			// The cell is outside the shield.
			cell->set_material_id(mid_3);
		}

		for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; f++)
		{
			double dif_norm = 0.0;
			for (unsigned int v = 1; v < GeometryInfo<dim>::vertices_per_face; v++)
				dif_norm += std::abs(	cell->face(f)->vertex(0).norm() -	cell->face(f)->vertex(v).norm());

			if ( dif_norm < eps )
				cell->face(f)->set_all_manifold_ids(1);
		}
	}

	Solver<dim>::triangulation.set_manifold(1,sphere);
}

template<int dim>
void SolverSLDI<dim>::solve()
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

